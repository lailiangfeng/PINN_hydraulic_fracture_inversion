# FEM Simulation Results (XDMF descriptors)

This directory contains the **XDMF descriptor files** of the four finite-element
reference cases (staggered phase-field / poroelastic FEniCSx solver) used in the
associated manuscript.

## Layout

```
result/
|-- case1-fem/
|   |-- displacement.xdmf
|   |-- pressure.xdmf
|   |-- phase_field.xdmf
|   |-- energy.xdmf
|   `-- history.xdmf
|-- case2-fem/   (same structure)
|-- case3-fem/   (same structure)
`-- case4-fem/   (same structure)
```

Each case directory contains five fields:

| File                  | Field                                             |
|-----------------------|---------------------------------------------------|
| `displacement.xdmf`   | Displacement field (vector, nodal)                |
| `pressure.xdmf`       | Fluid pressure field                              |
| `phase_field.xdmf`    | Phase-field damage variable                       |
| `energy.xdmf`         | Energy quantity                                   |
| `history.xdmf`        | History field (phase-field irreversibility)       |

The meshes are triangular (210k–105k elements depending on the case) with
time-series output; the files can be opened in **ParaView**.

## Important: the XDMF files require the paired HDF5 data

XDMF is a lightweight XML descriptor: it only references the actual field data
stored in the companion `.h5` files (e.g. `displacement.h5`,
`pressure.h5`, ...). The `.xdmf` files in this repository **cannot be opened
without the corresponding `.h5` datasets**, which total ≈ 3.7 GB and exceed
GitHub's per-file size limit.

### Data availability

The full HDF5 datasets (`*.h5`, ≈ 3.7 GB in total) are hosted separately and
will be published on a public research-data repository (Figshare / Zenodo).
Once the data record is published, the DOI / access link will be added here:

> **TODO** — add Figshare/Zenodo record URL or DOI for the `.h5` datasets.

Until the data record is available, requests for access to the HDF5 datasets
should be directed to the corresponding author and are handled under the
applicable institutional and publication policies (see the repository
README / CITATION.cff).
