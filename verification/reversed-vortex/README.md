# Reversed-vortex verification: retracing vs non-retracing

Reproduces the convergence study for the 2D reversed-vortex advection test with
and without a superposed periodic Euclidean frame (Bothe & Maric). Without the
frame the exact motion after `T/2` retraces the motion before it, so numerical
errors can cancel between the two legs; with the frame the final-time flow map is
still the identity but the laboratory-frame path no longer retraces, and each
scheme's own convergence rate is measured.

## One command

```bash
./run.sh
```

Picks the SLURM profile when `sbatch` exists and the local profile otherwise.
On Lichtenberg it builds the library, submits one job per case, and produces

```
results/convergence.csv          all errors and observed orders
results/convergence_table.tex    LaTeX tables, one per reconstruction scheme
results/interface_moving_frame.png
results/interface_non_moving_frame.png
```

`SKIP_BUILD=1 ./run.sh` reuses an existing build. Extra arguments are passed to
snakemake, e.g. `./run.sh slurm -n` for a dry run.

## What it runs

`config.yaml` sets the grid: resolutions x {plicRDF, isoAlpha} x {none, frameA}.
`np: 1` keeps each solve serial, which reproduces the published errors of
`testsuite/advection/test-vortexShearedDisc` exactly; the cases themselves run
concurrently, one sbatch each. Raising `np` switches the solve to
decomposePar + MPI + reconstructPar.

The frame is one revolution about the domain centre. That centre is not
arbitrary: because the base field is `alpha(t) u0` with a pure scalar amplitude,
trajectories never leave their `psi0` streamline, so the material support is
confined to a band of radius 0.4417 about `(0.5, 0.5)` and the unit square
suffices with a margin of 0.0583. At `N=32` that margin is only 1.9 cells and
the framed case loses volume through the boundary -- start the refinement study
at `N=64`.

## Cluster notes

`profiles/slurm/config.yaml` follows the conventions verified in `leia/CLUSTER.md`:
account `special00004`, a mandatory `--mem-per-cpu`, partition auto-routed from
the runtime, and `srun --ntasks={np} --overlap --cpu-bind=none` as the MPI
launcher. The build uses its own `WM_PROJECT_USER_DIR`
(`$HOME/OpenFOAM/tpf-rfv-v2512`) so it never mixes with other builds.

## Data archive

The figures and tables in the manuscript are regenerated from this study's
`results/` directory, which is what the archive cited as `figshare2026`
contains:

```
<study>/results/convergence.csv               errors and observed orders
<study>/results/convergence_table.tex         the LaTeX tables
<study>/results/convergence_shape_error.png   the convergence diagram
<study>/results/interface*_*_frame.png        the interface snapshots
<study>/results/frame_identity_check.txt      code-path equivalence: the
                                              moving-frame path with the frame
                                              set to the identity, against the
                                              original retracing test, at every
                                              output time
```

Every one of those is an output of `./run.sh`; none is produced by hand.
