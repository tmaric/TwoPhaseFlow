## The lower bound in parallel: deformation3D N 64, target Co 0.5, PURE, nAlphaBounds 6 (2026-11-24)

Run: job 57115478, 8 ranks, `ringSumBoundTrace 1`, the Part M binary (libVoF sha256 c1da5c86ddf849a9), through `d2_iso7.sbatch` (output `iso_bpar/`). Solver Courant labelled as in Part M. The same case ran serial in the bounds control (min alpha after the last pass -2.3e-17, 6 or 12 passes) and on 8 ranks in Part M (-4.9e-6).

What the trace prints in parallel: the pass-by-pass minima and maxima are global (`gMin`, `gMax`). Everything after the loop is computed from the cells of the master rank only and printed by the master only: the minimum and maximum after the last pass, the counts of cells outside the band, inside and outside needBounding, and the five worst cells (local cell labels, no processor-patch or neighbour information, no reduction over ranks). In 11 of the 960 steps the master's printed minimum after the last pass is 0 with max - 1 = -1 (its subdomain holds no alpha), and in no step did the master's five worst cells print. The per-cell question of this part (which cells below -1e-9 after the last pass have a face on a processor patch, are within one cell of one, or neither) cannot be answered by the existing trace: it has no global cell list and no processor-patch flag. Answering it needs a change of the trace code (reduction over ranks, a flag from the processor patches); this part is stopped here as instructed.

What the global pass-by-pass numbers give:

| quantity | value |
|---|---|
| steps | 960 |
| steps with the global minimum after the last pass below -1e-9 | 21 |
| step of the lowest final alpha | 768 |
| min alpha before the bounding at that step | -0.001907 |
| after pass 1 to 6 | -7.62e-4, -7.56e-4, -3.63e-5, -2.05e-5, -1.03e-5, -4.93e-6 |
| passes run | 6 (the limit; the serial run stopped at 4 passes at its worst step) |
| lowest final alpha of the run | -4.93e-6 (the solver's own `min(alpha1)`) |

HP1 (at least 80 % of the cells below -1e-9 after the last pass have a face on a processor patch or share a face with such a cell; falsifier: fewer than 40 %): not evaluated; the trace cannot report it. The pass-by-pass numbers show the bounding not converged after 6 passes on 8 ranks at 21 of the 960 steps, where the serial run converges (to round-off) in at most 4 passes.
