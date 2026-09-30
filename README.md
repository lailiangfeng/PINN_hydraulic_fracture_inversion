# FEniCSx Phase-Field Hydraulic-Fracture Demo

This repository contains a compact, runnable demonstration of a staggered
phase-field/poroelastic workflow in FEniCSx. It is intended for code review,
teaching, and interface testing.

The included case uses **nondimensional demonstration parameters** and a
simplified geometry. It is not the production configuration used to generate
the figures or numerical values of the associated manuscript. Paper-specific
meshes, fracture coordinates, loading histories, high-fidelity datasets,
trained checkpoints, and plotting source files are maintained separately.

## Model included in the demo

The program advances three coupled fields:

1. pressure diffusion with damage-dependent permeability;
2. degraded linear elasticity with Biot pressure coupling;
3. an AT2-type phase-field equation driven by an irreversible history field.

The fields are solved sequentially at each time step. The implementation is a
minimal research demo and is not a validated engineering simulator.

## Repository layout

    hf-phasefield-demo/
    |-- configs/demo_case.yaml
    |-- docs/MODEL_SCOPE.md
    |-- src/hf_demo/
    |   |-- config.py
    |   |-- solver.py
    |   |-- cli.py
    |-- run_demo.py
    |-- environment.yml
    |-- pyproject.toml
    |-- CITATION.cff
    |-- LICENSE

## Installation

FEniCSx is most reliably installed from conda-forge.

    conda env create -f environment.yml
    conda activate hf-phasefield-demo
    python -m pip install -e . --no-deps

The environment is pinned to the FEniCSx 0.8 API used by this demo.

## Run the demonstration

Serial execution:

    python run_demo.py --config configs/demo_case.yaml

MPI execution:

    mpiexec -n 4 python run_demo.py --config configs/demo_case.yaml

Results are written to outputs/demo_case/:

- displacement.xdmf
- pressure.xdmf
- damage.xdmf
- history.csv
- resolved_config.yaml

The XDMF files can be opened in ParaView.

## Create another public example

Copy configs/demo_case.yaml, change the nondimensional values, and pass the new
file through --config. Configuration validation rejects missing keys,
nonphysical Poisson ratios, invalid time steps, and out-of-domain pre-cracks.

## Reproducibility boundary

This public demonstration establishes the solver structure and file formats.
It does not claim to reproduce manuscript tables or figures. Requests for
access to restricted research configurations or datasets should be directed
to the corresponding author and handled under the applicable institutional
and publication policies.

## Citation

If this demonstration is used in academic work, cite the associated article
after its bibliographic record is available. A placeholder record is provided
in CITATION.cff and should be updated before creating a public release.

## License

Copyright (c) 2026 Liangfeng Lai. All rights reserved. See LICENSE.
