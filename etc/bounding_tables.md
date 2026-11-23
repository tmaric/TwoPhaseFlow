## The lower bound of the DAVOF rows: the bounding trace (2026-11-24)

deformation3D N 64, target Co 0.5 (the solver's own Courant maximum is printed in the logs), serial, rows plicRDF and PURE, `ringSumBoundTrace 1`; nAlphaBounds 3, clip false, snapTol 0, boundSecondLayer false (the davof rows' bounding covers surfCells_ plus one layer, the first layer of extendMarkedCells). 'inside' = inside needBounding.

| row | steps | min alpha before the bounding (over the run) | after the last pass (over the run) | max alpha - 1 before / after | steps with min alpha < -1e-6 after the last pass | in those: worst cell outside needBounding | in those: more cells above 1e-6 outside than inside | passes run per step (max) | steps where the last pass left a violation above 1e-12 |
|---|---|---|---|---|---|---|---|---|---|
| plicRDF | 960 | -0.0118 | -0.000205 | 0.0179 / 0.000309 | 67 | 17 | 26 | 3 | 943 |
| PURE | 960 | -0.0287 | -0.00527 | 0.0256 / 2e-05 | 265 | 200 | 257 | 3 | 939 |

HB1: in PURE the cells left below zero after the last pass lie outside needBounding in the majority of the steps with min alpha < -1e-6; in plicRDF inside. The columns 'worst cell outside needBounding' and 'more cells above 1e-6 outside than inside' count those steps (of the steps in the previous column).

### The five worst cells of the run, plicRDF (distinct cells, by violation)

| step | cell | violation | alpha before the bounding | alpha after | in needBounding | interface cell at t^n | outflow faces with dVf != 0 |
|---|---|---|---|---|---|---|---|
| 716 | 177064 | 0.000309 | 1.00048 | 1.00031 | yes | yes | 3 |
| 716 | 60136 | 0.000309 | 1.00048 | 1.00031 | yes | yes | 3 |
| 59 | 79062 | 0.000205 | -0.000815946 | -0.000205227 | yes | yes | 3 |
| 59 | 74966 | 0.000157 | 0.000204379 | -0.000156588 | yes | yes | 2 |
| 59 | 78998 | 0.000157 | 0.000204379 | -0.000156588 | yes | yes | 2 |

Pass by pass at the step of the lowest alpha (plicRDF, step 59): before min -0.00664, max - 1 0.000917; after pass 1 min -0.00175, max - 1 2.22e-16; after pass 2 min -0.000778, max - 1 2.22e-16; after pass 3 min -0.000205, max - 1 2.22e-16; after the last pass min -0.000205, max - 1 2.22e-16; cells outside the band inside needBounding 162, outside 2501.

### The five worst cells of the run, PURE (distinct cells, by violation)

| step | cell | violation | alpha before the bounding | alpha after | in needBounding | interface cell at t^n | outflow faces with dVf != 0 |
|---|---|---|---|---|---|---|---|
| 922 | 203739 | 0.00527 | 0.000498583 | -0.00526909 | yes | yes | 1 |
| 892 | 174637 | 0.00441 | -0.0077165 | -0.00441006 | yes | yes | 2 |
| 923 | 203738 | 0.00387 | -0.00194464 | -0.00386851 | yes | no | 0 |
| 893 | 178860 | 0.00315 | 0.0123352 | -0.00315453 | yes | yes | 2 |
| 924 | 203737 | 0.00302 | -0.00139001 | -0.00301846 | yes | no | 3 |

Pass by pass at the step of the lowest alpha (PURE, step 922): before min -0.0225, max - 1 0.00201; after pass 1 min -0.00675, max - 1 0.000561; after pass 2 min -0.00803, max - 1 9.76e-10; after pass 3 min -0.00527, max - 1 9.76e-10; after the last pass min -0.00527, max - 1 9.76e-10; cells outside the band inside needBounding 1007, outside 65627.

