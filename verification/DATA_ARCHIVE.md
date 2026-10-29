# Data archive: periodic, non-retracing deformation tests

This archive accompanies D. Bothe and T. Marić, *Periodic, non-retracing
deformation test cases for interface advection*. It holds the code, the
workflows and every result the manuscript reports. Figure and table captions in
the manuscript name the file below that each was produced from.

## Layout

| Path | Contents |
|---|---|
| `src/VoF/movingFrameFlow/`, `src/VoF/movingFrameFlow3D/` | the laboratory-frame flux assembly (stream function, vector potentials) |
| `apps/benchmark/advectorVoF/` | the advection solver |
| `verification/reversed-vortex/` | Section 4: 2D reversed vortex |
| `verification/deformation-sphere3d/` | Section 5: 3D shear flow |
| `verification/leveque3d/` | Section 6: 3D LeVeque deformation |
| `verification/<study>/results/` | every figure, table and check the manuscript quotes |
| `verification/<study>/runs/<case>/` | per run: the error time series and the case's `system/` |
| `MANIFEST.txt` | the commit the archive was built from, and a checksum of every file |

Case names are `<scheme>_<frame>_N<resolution>`: scheme `plicRDF` or `isoAlpha`;
frame `none` (the original, retracing test), `frameA` (one revolution
superposed), `frameId` (the identity frame) or `frameOnly` (the frame without
the base field). `viz_*` cases write more often and supply the interface figures
only.

## Results per study

| File | Manuscript |
|---|---|
| `results/convergence.csv`, `convergence_table.tex` | convergence tables |
| `results/convergence_shape_error.{pdf,png}` | convergence figures |
| `results/interface*_moving_frame.png`, `*_non_moving_frame.png` | interface snapshots |
| `results/timestep_bound.json` | the bound max‖w‖₁ that fixes the time step |
| `results/courant_check.txt` | per run: Courant number, discrete divergence, volume and boundedness errors reached |
| `reversed-vortex/results/frame_identity_check.txt` | identity frame against the original code path |
| `reversed-vortex/results/frame_only_check.txt` | the frame alone (α ≡ 0) |
| `reversed-vortex/results/interface_evolution.{pdf,png,gif}` | snapshots every 0.1T against the exact interface; animation |
| `<3D study>/results/potential_check.txt` | curl of the vector potential against the base field |
| `<3D study>/results/support_extent.txt` | how far the material reaches: the domain bound |

## Reproducing

Each study regenerates all of its results with one command, from its directory:

```bash
./run.sh            # SLURM if sbatch is available, otherwise local
```

See `verification/<study>/README.md` for the protocol and the requirements
(OpenFOAM v2512, Snakemake 9).
