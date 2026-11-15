# CPC revision reproduction data, version 1.1

This package contains independent Lichtenberg II reproduction runs for the November 2026 revisions of the finite-volume Helmholtz/PML method paper and the `acousticHelmholtzFoam` software paper. The implementation is `TwoPhaseFlow` commit `7f00fc55a67aefbcfa88bdb46a4fc378bac0db0c`, also identified by release tag `acousticHelmholtzFoam-cpc-2026-11-15`. The data tag is `helmholtz-pml-data-v1.1`.

The package supersedes the September v1.0 package for assessing the revised numerical implementation. Version 1.0 remains unchanged in `data/publications/helmholtz-pml/`. Its reserved TUdatalib DOI is not represented as an identifier for this new payload. Upload of v1.1 and its archive/version identifier must be confirmed separately; see `metadata-draft.json`.

## Scope and provenance

The paper's numerical plots and tables report Chuanchao Xu's workstation run. This package supplies independently computed numerical counterparts, the current manuscript sources and exact figure assets for identification, and the scripts and execution logs for the HPC runs. It does not claim to recover the original workstation raw fields or timing samples.

The HPC campaign comprises 576 numerical cases: 24 homogeneous-wave cases, 66 layered serial/MPI comparisons, 16 penetrable-sphere comparisons, 44 sphere/PML cases, five piston meshes, three rigid-sphere mesh levels, seven sphere radii, eleven sphere positions, and 400 levitator height points. Twenty compiled integration tests also passed. `VALIDATION.json` records checks against the printed workstation tables, including the small cross-build changes in the rigid-sphere force. The reported workstation performance measurements are retained separately as transcribed aggregates; they are not HPC timing measurements.

The height points are independent. Completed points from four initial jobs were retained; the remaining points were divided among sixteen jobs on the same 100-point grids. `execution/height-shard-plan.json`, the job table and each series' combination record identify that assembly. The original four jobs were intentionally cancelled for this redistribution, not because their solver runs failed. Binary hashes are checked before combining the points.

## Contents and use

- `secondary/`: CSV, TSV and JSON results; layered pressure/reference profiles; piston axis, arc and compressed cell-centre fields; force records; and the four complete height curves.
- `execution/`: exact campaign scripts, build/test log, job identifiers, accounting and shard records.
- `scripts/original/`: the case postprocessing and setup sources used for these runs. Re-running the numerical solver requires OpenFOAM and the pinned dependencies.
- `scripts/validate.py`: portable checks of hashes, completeness, finiteness, the height grid, force conversion and agreement with manuscript tables.
- `scripts/reproduce.py`: regenerates numerical plots from these HPC data without OpenFOAM; plot styling differs from the manuscript.
- `manuscript/`: identification snapshot of both revised manuscripts, bibliographies, table inputs and referenced figure assets.
- `data-dictionary.json`: paths, columns and row counts for uncompressed CSV/TSV files. Column suffixes identify SI units; dimensionless relative errors are fractions unless explicitly marked percent. Complex pressures use the manuscript's time convention. Blank convergence orders denote the coarsest grid. `seconds` in comparison JSON includes case-run overhead and is not a performance benchmark. Piston `.csv.gz` files have the header `r_m,z_m,V_m3,Pre_sim_pa,Pim_sim_pa,Pre_exact_pa,Pim_exact_pa`.
- `software-provenance.json`, `hpc-status.json`, `checksums.json`: software identity, completion and SHA-256 payload identity.

With Python 3.10 or newer and the packages in `requirements.txt`:

```sh
python3 scripts/validate.py
python3 scripts/reproduce.py --output ../hpc-reproduced
python3 scripts/make_release.py --output ../helmholtz-pml-secondary-data-v1.1.zip
```

Use new output paths outside the package. The ZIP preserves the directory structure, with deterministic ordering/timestamps and stored entries. Raw meshes, all distributed fields, solver binaries and the original workstation timing repetitions are excluded. Logs and source hashes identify the runs; they do not establish experimental validation of the levitator or drop stability.

The checksum manifest covers all payload files except itself and the generated `VALIDATION.json` report; running validation recomputes the scientific checks and that report.

See `PROVENANCE.md` for the manuscript mapping and `RIGHTS.md` for licenses.
