
HPCToolkit was installed under the user home directory and exposed via
`~/.local/bin`:

```bash
hpcrun --version
hpcstruct --version
hpcprof --version
hpcviewer --version
```

The installed HPCToolkit version used here was `2025.1.2`. It has PAPI support,
but this machine has `kernel.perf_event_paranoid=4`, so hardware-counter events
are not available to normal users. The profile therefore used timer sampling
with `CPUTIME`.

## Run HPCToolkit measurement

The successful HPCToolkit measurement was run with 4 MPI ranks:

```bash
ts=20260504-132230
meas="hpctoolkit-acousticHelmholtzFoam-measurements-4rank-${ts}"
log="log.hpcrun.acousticHelmholtzFoam.4rank.${ts}"

OMP_NUM_THREADS=1 \
OPENBLAS_NUM_THREADS=1 \
MKL_NUM_THREADS=1 \
VECLIB_MAXIMUM_THREADS=1 \
NUMEXPR_NUM_THREADS=1 \
mpirun -np 4 \
    hpcrun -t -o "${meas}" \
    acousticHelmholtzFoam -parallel \
    > "${log}" 2>&1
```

Successful output:

```text
hpctoolkit-acousticHelmholtzFoam-measurements-4rank-20260504-132230
log.hpcrun.acousticHelmholtzFoam.4rank.20260504-132230
```

## Build HPCToolkit structure and database

First, structure information was generated for the solver binary:

```bash
hpcstruct "$(command -v acousticHelmholtzFoam)" \
    -o acousticHelmholtzFoam.hpcstruct \
    > log.hpcstruct.acousticHelmholtzFoam 2>&1
```

Then HPCToolkit structure information was generated for all binaries referenced
by the measurement directory:

```bash
hpcstruct hpctoolkit-acousticHelmholtzFoam-measurements-4rank-20260504-132230 \
    > log.hpcstruct.measurements.20260504-132230 2>&1
```

Finally, the viewer database was generated:

```bash
hpcprof \
    --force \
    --only-exe acousticHelmholtzFoam \
    -n "acousticHelmholtzFoam WedgeLevitatorDropAlpha 4 ranks" \
    -o hpctoolkit-acousticHelmholtzFoam-database-4rank-20260504-132230 \
    hpctoolkit-acousticHelmholtzFoam-measurements-4rank-20260504-132230 \
    > log.hpcprof.acousticHelmholtzFoam.4rank.20260504-132230.with-structs 2>&1
```

Final database:

```text
hpctoolkit-acousticHelmholtzFoam-database-4rank-20260504-132230
```

## Open the result

Open the generated database with:

```bash
hpcviewer hpctoolkit-acousticHelmholtzFoam-database-4rank-20260504-132230
```

Useful views in `hpcviewer`:

- `Calling Context View`: call-path hierarchy annotated with `CPUTIME`.
- `Flat View`: top functions/modules by CPU time.
- Trace view data is present in the database, but this install exposes it
  through the same GUI package rather than a separate `hpctraceviewer` command.

## Current profiling limitations

- Only `CPUTIME (s)` is available in this profile.
- Hardware counters such as cache misses or branch mispredictions are blocked by
  the kernel setting `perf_event_paranoid=4`.
- To use hardware counters, an administrator must lower that kernel setting.
