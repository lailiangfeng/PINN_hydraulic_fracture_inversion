"""Staggered FEniCSx solver for the public nondimensional demo case."""

from __future__ import annotations

import csv
from pathlib import Path
from typing import Any

import numpy as np
import ufl
import yaml
from dolfinx import fem, io, mesh
from dolfinx.fem.petsc import LinearProblem
from mpi4py import MPI
from petsc4py import PETSc


def _global_max(comm: MPI.Comm, value: float) -> float:
    return float(comm.allreduce(float(value), op=MPI.MAX))


def _elastic_energy(u, lame_lambda: float, shear_mu: float):
    eps = ufl.sym(ufl.grad(u))
    return 0.5 * lame_lambda * ufl.tr(eps) ** 2 + shear_mu * ufl.inner(eps, eps)


def run_simulation(config: dict[str, Any]) -> Path:
    """Run the demonstration and return the output directory."""
    comm = MPI.COMM_WORLD
    rank = comm.rank

    domain_cfg = config["domain"]
    material = config["material"]
    flow = config["flow"]
    time_cfg = config["time"]
    crack = config["precrack"]

    lx = float(domain_cfg["length_x"])
    ly = float(domain_cfg["length_y"])
    nx = int(domain_cfg["cells_x"])
    ny = int(domain_cfg["cells_y"])

    domain = mesh.create_rectangle(
        comm,
        [np.array([-0.5 * lx, -0.5 * ly]), np.array([0.5 * lx, 0.5 * ly])],
        [nx, ny],
        cell_type=mesh.CellType.triangle,
    )

    gdim = domain.geometry.dim
    displacement_space = fem.functionspace(
        domain, ("Lagrange", 1, (gdim,))
    )
    scalar_space = fem.functionspace(domain, ("Lagrange", 1))

    displacement = fem.Function(displacement_space, name="displacement")
    pressure = fem.Function(scalar_space, name="pressure")
    pressure_old = fem.Function(scalar_space, name="pressure_old")
    damage = fem.Function(scalar_space, name="damage")
    damage_old = fem.Function(scalar_space, name="damage_old")
    history = fem.Function(scalar_space, name="history")
    source = fem.Function(scalar_space, name="source")

    l0 = float(material["regularization_length"])
    crack_half = float(crack["half_length"])
    crack_x = float(crack["center_x"])
    crack_y = float(crack["center_y"])

    def initial_crack(x: np.ndarray) -> np.ndarray:
        along = np.maximum(np.abs(x[0] - crack_x) - crack_half, 0.0)
        normal = np.abs(x[1] - crack_y)
        distance = np.sqrt(along**2 + normal**2)
        return np.exp(-((distance / l0) ** 2))

    damage_old.interpolate(initial_crack)
    damage_old.x.array[:] = np.clip(damage_old.x.array, 0.0, 1.0)
    damage_old.x.scatter_forward()
    damage.x.array[:] = damage_old.x.array
    damage.x.scatter_forward()

    source_amplitude = float(flow["source_amplitude"])
    source_width = float(flow["source_width"])

    def gaussian_source(x: np.ndarray) -> np.ndarray:
        radius2 = (x[0] - crack_x) ** 2 + (x[1] - crack_y) ** 2
        return source_amplitude * np.exp(
            -radius2 / (2.0 * source_width**2)
        )

    source.interpolate(gaussian_source)
    source.x.scatter_forward()

    fdim = domain.topology.dim - 1
    bottom_facets = mesh.locate_entities_boundary(
        domain, fdim, lambda x: np.isclose(x[1], -0.5 * ly)
    )
    displacement_dofs = fem.locate_dofs_topological(
        displacement_space, fdim, bottom_facets
    )
    zero_vector = np.zeros(gdim, dtype=PETSc.ScalarType)
    displacement_bc = fem.dirichletbc(
        zero_vector, displacement_dofs, displacement_space
    )

    outer_facets = mesh.locate_entities_boundary(
        domain,
        fdim,
        lambda x: np.isclose(np.abs(x[0]), 0.5 * lx)
        | np.isclose(np.abs(x[1]), 0.5 * ly),
    )
    pressure_dofs = fem.locate_dofs_topological(
        scalar_space, fdim, outer_facets
    )
    pressure_bc = fem.dirichletbc(
        PETSc.ScalarType(0.0), pressure_dofs, scalar_space
    )

    young = float(material["young_modulus"])
    nu = float(material["poisson_ratio"])
    shear_mu = young / (2.0 * (1.0 + nu))
    lame_lambda = young * nu / ((1.0 + nu) * (1.0 - 2.0 * nu))
    gc = float(material["fracture_toughness"])
    residual_stiffness = float(material["residual_stiffness"])
    biot = float(material["biot_coefficient"])

    storage = float(flow["storage"])
    matrix_perm = float(flow["matrix_permeability"])
    fracture_perm = float(flow["fracture_permeability"])
    viscosity = float(flow["viscosity"])

    dt = float(time_cfg["step_size"])
    number_of_steps = int(time_cfg["number_of_steps"])
    tolerance = float(time_cfg["staggered_tolerance"])
    max_iterations = int(time_cfg["staggered_max_iterations"])
    output_interval = int(time_cfg["output_interval"])

    output_dir = Path(config["output"]["directory"]).resolve()
    if rank == 0:
        output_dir.mkdir(parents=True, exist_ok=True)
        with (output_dir / "resolved_config.yaml").open(
            "w", encoding="utf-8"
        ) as stream:
            yaml.safe_dump(config, stream, sort_keys=False)
        history_file = (output_dir / "history.csv").open(
            "w", newline="", encoding="utf-8"
        )
        history_writer = csv.writer(history_file)
        history_writer.writerow(
            [
                "step",
                "time",
                "staggered_iterations",
                "damage_increment",
                "max_damage",
            ]
        )
    else:
        history_file = None
        history_writer = None
    comm.barrier()

    dx = ufl.dx(domain=domain)
    trial_p = ufl.TrialFunction(scalar_space)
    test_p = ufl.TestFunction(scalar_space)
    trial_u = ufl.TrialFunction(displacement_space)
    test_u = ufl.TestFunction(displacement_space)
    trial_d = ufl.TrialFunction(scalar_space)
    test_d = ufl.TestFunction(scalar_space)

    xdmf_u = io.XDMFFile(
        comm, str(output_dir / "displacement.xdmf"), "w"
    )
    xdmf_p = io.XDMFFile(comm, str(output_dir / "pressure.xdmf"), "w")
    xdmf_d = io.XDMFFile(comm, str(output_dir / "damage.xdmf"), "w")
    for xdmf in (xdmf_u, xdmf_p, xdmf_d):
        xdmf.write_mesh(domain)

    xdmf_u.write_function(displacement, 0.0)
    xdmf_p.write_function(pressure, 0.0)
    xdmf_d.write_function(damage, 0.0)

    linear_solver = {"ksp_type": "preonly", "pc_type": "lu"}

    try:
        for step in range(1, number_of_steps + 1):
            current_time = step * dt

            permeability = (
                matrix_perm
                + (fracture_perm - matrix_perm) * damage_old**2
            )
            pressure_a = (
                storage * trial_p * test_p * dx
                + dt
                * (permeability / viscosity)
                * ufl.dot(ufl.grad(trial_p), ufl.grad(test_p))
                * dx
            )
            pressure_l = (
                storage * pressure_old * test_p * dx
                + dt * source * test_p * dx
            )
            pressure_problem = LinearProblem(
                pressure_a,
                pressure_l,
                bcs=[pressure_bc],
                petsc_options=linear_solver,
            )
            pressure_solution = pressure_problem.solve()
            pressure.x.array[:] = pressure_solution.x.array
            pressure.x.scatter_forward()

            damage_increment = np.inf
            iteration = 0

            while (
                iteration < max_iterations
                and damage_increment > tolerance
            ):
                previous_iteration = damage.x.array.copy()

                eps_trial = ufl.sym(ufl.grad(trial_u))
                sigma_trial = (
                    2.0 * shear_mu * eps_trial
                    + lame_lambda
                    * ufl.tr(eps_trial)
                    * ufl.Identity(gdim)
                )
                degradation = (
                    (1.0 - damage) ** 2 + residual_stiffness
                )
                mechanics_a = (
                    degradation
                    * ufl.inner(
                        sigma_trial,
                        ufl.sym(ufl.grad(test_u)),
                    )
                    * dx
                )
                mechanics_l = (
                    biot * pressure * ufl.div(test_u) * dx
                )

                mechanics_problem = LinearProblem(
                    mechanics_a,
                    mechanics_l,
                    bcs=[displacement_bc],
                    petsc_options=linear_solver,
                )
                displacement_solution = mechanics_problem.solve()
                displacement.x.array[:] = (
                    displacement_solution.x.array
                )
                displacement.x.scatter_forward()

                energy_expression = fem.Expression(
                    _elastic_energy(
                        displacement, lame_lambda, shear_mu
                    ),
                    scalar_space.element.interpolation_points(),
                )
                energy = fem.Function(scalar_space)
                energy.interpolate(energy_expression)
                history.x.array[:] = np.maximum(
                    history.x.array, energy.x.array
                )
                history.x.scatter_forward()

                damage_a = (
                    gc
                    * l0
                    * ufl.dot(
                        ufl.grad(trial_d), ufl.grad(test_d)
                    )
                    * dx
                    + (gc / l0 + 2.0 * history)
                    * trial_d
                    * test_d
                    * dx
                )
                damage_l = 2.0 * history * test_d * dx
                damage_problem = LinearProblem(
                    damage_a,
                    damage_l,
                    petsc_options=linear_solver,
                )
                damage_solution = damage_problem.solve()
                damage.x.array[:] = np.maximum(
                    damage_old.x.array,
                    np.clip(
                        damage_solution.x.array, 0.0, 1.0
                    ),
                )
                damage.x.scatter_forward()

                if damage.x.array.size:
                    local_increment = float(
                        np.max(
                            np.abs(
                                damage.x.array
                                - previous_iteration
                            )
                        )
                    )
                else:
                    local_increment = 0.0
                damage_increment = _global_max(
                    comm, local_increment
                )
                iteration += 1

            pressure_old.x.array[:] = pressure.x.array
            pressure_old.x.scatter_forward()
            damage_old.x.array[:] = damage.x.array
            damage_old.x.scatter_forward()

            local_max_damage = (
                float(np.max(damage.x.array))
                if damage.x.array.size
                else 0.0
            )
            max_damage = _global_max(
                comm, local_max_damage
            )

            if rank == 0:
                history_writer.writerow(
                    [
                        step,
                        current_time,
                        iteration,
                        damage_increment,
                        max_damage,
                    ]
                )
                history_file.flush()
                print(
                    f"step={step:04d} "
                    f"time={current_time:.4f} "
                    f"iterations={iteration:02d} "
                    f"damage_increment={damage_increment:.3e} "
                    f"max_damage={max_damage:.4f}",
                    flush=True,
                )

            if (
                step % output_interval == 0
                or step == number_of_steps
            ):
                xdmf_u.write_function(
                    displacement, current_time
                )
                xdmf_p.write_function(
                    pressure, current_time
                )
                xdmf_d.write_function(
                    damage, current_time
                )
    finally:
        xdmf_u.close()
        xdmf_p.close()
        xdmf_d.close()
        if history_file is not None:
            history_file.close()

    return output_dir
