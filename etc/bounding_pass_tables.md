## The bounds control: the number of bounding passes (2026-11-23)

deformation3D N 64, target Co 0.5, serial, `ringSumBoundTrace 1`, clip false, snapTol 0, boundSecondLayer false; the same binary for all runs (the fork build of f9208c4, sha256 616b3500...; the runs at 3 passes of plicRDF and PURE are job 57098598, the others jobs 57105150 to 57105155 with the copy of that build). Rows: plicRDF; PURE; (c) the PURE loop with normalModel plicRDF. nAlphaBounds in the case's fvSolution is the first column. 'inside' = inside needBounding.

| row | nAlphaBounds (fvSolution) | status | solver max Courant | min alpha before the bounding | min alpha after the last pass | max alpha - 1 before | max alpha - 1 after | steps with min alpha < -1e-6 after the last pass | in those: worst cell outside needBounding | E_shape (end) | mass error (%, max) | wall (s) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| plicRDF | 3 | completed | 0.4950 | -0.01176 | -0.0002052 | 0.0179 | 0.000309 | 67 | 17 | 0.002979 | 1.79e-11 | 321 |
| plicRDF | 6 | completed | 0.4950 | -0.01565 | -5.164e-08 | 0.0109 | 6.66e-16 | 0 | 0 | 0.003019 | 3.66e-11 | 316 |
| plicRDF | 12 | completed | 0.4950 | -0.01565 | -7.033e-13 | 0.0109 | 6.66e-16 | 0 | 0 | 0.003019 | 3.57e-11 | 320 |
| PURE | 3 | completed | 0.4950 | -0.02866 | -0.005269 | 0.0256 | 2e-05 | 265 | 200 | 0.00386 | 1.93e-11 | 287 |
| PURE | 6 | completed | 0.4950 | -0.02973 | -2.255e-17 | 0.0256 | 4.44e-16 | 0 | 0 | 0.00386 | 3.19e-11 | 287 |
| PURE | 12 | completed | 0.4950 | -0.02973 | -2.255e-17 | 0.0256 | 4.44e-16 | 0 | 0 | 0.00386 | 3.19e-11 | 285 |
| (c) PURE loop, plicRDF normal | 3 | completed | 0.4950 | -0.01386 | -0.00027 | 0.0171 | 0.000297 | 56 | 14 | 0.003133 | 1.54e-11 | 335 |

### plicRDF, 3 passes

Pass by pass at the step of the lowest alpha (step 59): before min -0.006639; after pass 1 min -0.001751; after pass 2 min -0.0007782; after pass 3 min -0.0002052; after the last pass min -0.0002052, max - 1 2.22e-16; cells outside the band: inside needBounding 162, outside 2501.

| step | cell | violation | alpha before | alpha after | in needBounding | interface cell at t^n | outflow faces with dVf != 0 |
|---|---|---|---|---|---|---|---|
| 716 | 177064 | 0.000309 | 1.00048 | 1.00031 | yes | yes | 3 |
| 716 | 60136 | 0.000309 | 1.00048 | 1.00031 | yes | yes | 3 |
| 59 | 79062 | 0.000205 | -0.000815946 | -0.000205227 | yes | yes | 3 |
| 59 | 74966 | 0.000157 | 0.000204379 | -0.000156588 | yes | yes | 2 |
| 59 | 78998 | 0.000157 | 0.000204379 | -0.000156588 | yes | yes | 2 |

### plicRDF, 6 passes

Pass by pass at the step of the lowest alpha (step 55): before min -0.008078; after pass 1 min -0.001738; after pass 2 min -0.0004582; after pass 3 min -4.035e-05; after pass 4 min -1.939e-05; after pass 5 min -4.763e-06; after pass 6 min -5.164e-08; after the last pass min -5.164e-08, max - 1 2.22e-16; cells outside the band: inside needBounding 1, outside 0.

| step | cell | violation | alpha before | alpha after | in needBounding | interface cell at t^n | outflow faces with dVf != 0 |
|---|---|---|---|---|---|---|---|

### plicRDF, 12 passes

Pass by pass at the step of the lowest alpha (step 897): before min -0.001492; after pass 1 min -8.479e-08; after pass 2 min -1.904e-11; after pass 3 min -7.033e-13; after the last pass min -7.033e-13, max - 1 2.22e-16; cells outside the band: inside needBounding 0, outside 0.

| step | cell | violation | alpha before | alpha after | in needBounding | interface cell at t^n | outflow faces with dVf != 0 |
|---|---|---|---|---|---|---|---|

### PURE, 3 passes

Pass by pass at the step of the lowest alpha (step 922): before min -0.02248; after pass 1 min -0.006748; after pass 2 min -0.00803; after pass 3 min -0.005269; after the last pass min -0.005269, max - 1 9.76e-10; cells outside the band: inside needBounding 1007, outside 65627.

| step | cell | violation | alpha before | alpha after | in needBounding | interface cell at t^n | outflow faces with dVf != 0 |
|---|---|---|---|---|---|---|---|
| 922 | 203739 | 0.00527 | 0.000498583 | -0.00526909 | yes | yes | 1 |
| 892 | 174637 | 0.00441 | -0.0077165 | -0.00441006 | yes | yes | 2 |
| 923 | 203738 | 0.00387 | -0.00194464 | -0.00386851 | yes | no | 0 |
| 893 | 178860 | 0.00315 | 0.0123352 | -0.00315453 | yes | yes | 2 |
| 924 | 203737 | 0.00302 | -0.00139001 | -0.00301846 | yes | no | 3 |

### PURE, 6 passes

Pass by pass at the step of the lowest alpha (step 860): before min -0.0225; after pass 1 min -0.0001034; after pass 2 min -8.315e-05; after pass 3 min -3.484e-05; after pass 4 min -2.255e-17; after the last pass min -2.255e-17, max - 1 0; cells outside the band: inside needBounding 0, outside 0.

| step | cell | violation | alpha before | alpha after | in needBounding | interface cell at t^n | outflow faces with dVf != 0 |
|---|---|---|---|---|---|---|---|

### PURE, 12 passes

Pass by pass at the step of the lowest alpha (step 860): before min -0.0225; after pass 1 min -0.0001034; after pass 2 min -8.315e-05; after pass 3 min -3.484e-05; after pass 4 min -2.255e-17; after the last pass min -2.255e-17, max - 1 0; cells outside the band: inside needBounding 0, outside 0.

| step | cell | violation | alpha before | alpha after | in needBounding | interface cell at t^n | outflow faces with dVf != 0 |
|---|---|---|---|---|---|---|---|

### (c) PURE loop, plicRDF normal, 3 passes

Pass by pass at the step of the lowest alpha (step 59): before min -0.00644; after pass 1 min -0.002303; after pass 2 min -0.0007549; after pass 3 min -0.00027; after the last pass min -0.00027, max - 1 2.22e-16; cells outside the band: inside needBounding 144, outside 2329.

| step | cell | violation | alpha before | alpha after | in needBounding | interface cell at t^n | outflow faces with dVf != 0 |
|---|---|---|---|---|---|---|---|
| 716 | 177064 | 0.000297 | 1.00047 | 1.0003 | yes | yes | 3 |
| 716 | 60136 | 0.000297 | 1.00047 | 1.0003 | yes | yes | 3 |
| 59 | 79062 | 0.00027 | -0.00100071 | -0.000269999 | yes | yes | 3 |
| 308 | 165159 | 0.000204 | -1.73222e-05 | -0.000204198 | yes | yes | 2 |
| 308 | 84519 | 0.000204 | -1.73222e-05 | -0.000204198 | yes | yes | 2 |

### The pre-registered statements

- PB1 (row (c): min alpha before the bounding and after pass 3 within a factor 1.5 of plicRDF's; falsifier: after pass 3 below -1e-3): before -0.01386 against plicRDF's -0.01176 (ratio 1.179); after -0.00027 against -0.0002052 (ratio 1.316).
- PB2 (plicRDF: |min alpha| after the last pass at most 2e-5 with 6 passes and 1e-7 with 12; falsifier: above 1e-4 with 6): -5.164e-08 with 6 passes, -7.033e-13 with 12.
- PB3 (PURE: |min alpha| after the last pass above 5e-4 with 12 passes; falsifier: below 5e-5): -2.255e-17 with 6 passes, -2.255e-17 with 12.
- PB4 (plicRDF: E_shape changes by less than 1 % between 3, 6 and 12 passes): E_shape 0.0029786, 0.0030189, 0.0030186; changes 3 to 6 1.352 %, 6 to 12 -0.010 %, 3 to 12 1.342 %.
- PB4 (PURE: E_shape changes by less than 1 % between 3, 6 and 12 passes): E_shape 0.0038603, 0.0038602, 0.0038602; changes 3 to 6 -0.001 %, 6 to 12 0.000 %, 3 to 12 -0.001 %.

