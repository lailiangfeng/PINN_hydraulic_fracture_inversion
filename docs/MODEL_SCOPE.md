# Model scope

The public example uses nondimensional variables and a deliberately simplified
rectangular domain. Let u, p, and d denote displacement, pressure, and damage,
respectively, with d = 0 intact and d = 1 fully damaged.

The elastic degradation function is

    g(d) = (1 - d)^2 + k_res.

The mechanics problem uses

    div(g(d) sigma_0(u) - alpha p I) = 0,

where sigma_0 is the isotropic linear-elastic stress and alpha is the Biot
coefficient. The lower edge is fixed in the demonstration to remove rigid-body
motion.

Pressure evolves according to

    S (p_n - p_(n-1)) / dt - div(k(d) / mu grad(p_n)) = q,

with

    k(d) = k_m + (k_f - k_m) d^2.

The damage subproblem is an AT2-type linearized relation,

    Gc l0 grad(d) . grad(v)
    + (Gc / l0 + 2 H) d v
    = 2 H v,

where the history field H stores the largest elastic-energy value seen at each
interpolation point. After each damage solve, irreversibility is also enforced
algebraically by d_n >= d_(n-1).

## Deliberate limitations

- The values in configs/demo_case.yaml are nondimensional demo values.
- The mesh is uniform and does not reproduce a production mesh-refinement map.
- The source is a smooth Gaussian rather than a paper-specific injection
  schedule or wellbore boundary model.
- No manuscript datasets, trained PINN checkpoints, scalar-normalization
  values, or figure source data are included.
- The demo is not calibrated for predictive engineering use.

These limitations should remain visible when the repository is shared. Do not
describe this example as a full reproduction package for manuscript results.
