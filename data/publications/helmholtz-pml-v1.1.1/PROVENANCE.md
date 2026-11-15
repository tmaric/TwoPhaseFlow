# Provenance and manuscript mapping

Assembled 15 November 2026. Article input: Overleaf commit `af2711fef2765d93116b37e599e83953e1d179a0`, followed by the recorded finalization edits. The manuscript snapshot identifies the final source bytes. The editorial comparison baseline remains article commit `96437f22a802c953361f31b438e2d6a575fda20f`.

| Revised manuscript evidence | Independent HPC counterpart |
|---|---|
| Homogeneous convergence, method Table 1 and convergence plots | `secondary/homogeneous/studyResults/convergence/metrics_summary.csv` |
| Layered forward/reverse pressure curves | `secondary/layered/pressure-N8000-s500000-{forward,reverse}.csv` |
| Damping and aligned-mesh sensitivity, method Tables 2–3 | `secondary/layered/profile-metrics.csv`; the eleven pressure-profile CSV files preserve the numerical and exact complex values |
| Geometric face comparison, method Table 4 | `secondary/sphere/studyResults/results.json` |
| Combined transmission/PML, method Table 5 | `secondary/sphere-pml/studyResults/main/results.json`; `softPML` and `ranks8` preserve the supporting studies |
| Exact Rayleigh piston, method Table 6 and pressure/error plots | `secondary/piston/meshConvergence/metrics_summary.csv` and `N*/analyticalCompare/` |
| Rigid-sphere mesh table and radius/position plots | `secondary/gorkov/{meshConvergence,gorkovStudy,positionSweep}/` |
| Levitator reflector-force application in both papers | `secondary/height-{empty,1,2,3}/sweepResults/heightSweepResults.csv` |
| Software serial/MPI regression | `secondary/layered/comparison/results.json`, strong-damping forward configuration, ranks 1/4/8; the original workstation's two-rank comparison is not represented as an HPC run |
| Workstation strong-scaling table | `secondary/workstation-performance/strong-scaling-transcribed.csv`, explicitly transcribed from the software manuscript |
| Workstation problem-size timing figure | Exact original PDF and the stated endpoint aggregates are preserved in the manuscript snapshot; original timing samples are not in this package |
| Schematics, geometry and mesh illustrations | Exact source/PDF/PNG assets in `manuscript/`; these are identification assets rather than newly computed secondary data |

The paper's workstation results use OpenFOAM-v2606, PETSc 3.25.0 and MUMPS 5.8.2. Independent HPC results use the same implementation and PETSc source, OpenFOAM-v2606, GCC 11.5.0 and Open MPI 4.1.8, built in an isolated directory. Every campaign records the same solver and acoustic-interface library hashes. No solver sources were changed for the reproduction; the MPI include search path was supplied through the environment.

The compiled integration-test job is 56398204. The initial campaign jobs are 56398263–56398272; the four height jobs in that range were intentionally superseded after saving their completed rows. Height shard jobs are listed in `execution/height-shard-jobs.tsv`. The combination records and final accounting distinguish successful studies from intentionally stopped orchestration jobs.

The small differences between the workstation and HPC rigid-sphere forces are recorded by `scripts/validate.py`; values are compared within 0.01% of force, and no mechanism is assigned to this cross-build difference without a controlled investigation. Relative pressure errors and convergence orders are checked against the printed table precision. The source case and its analytical comparison are preserved separately from the plotted observables. The reflector force is not the force on the liquid drop.

The figure assets in the manuscript snapshot retain the workstation rendering. The portable reproduction produces plots from the HPC data and therefore promises numerical reproduction within the documented checks, not identical pixels. No experimental validation, solver superiority, transmitting-versus-sound-hard performance comparison, moving-interface calculation or drop-stability prediction is inferred from these runs.
