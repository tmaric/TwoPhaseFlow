## The area-weighted linear fit of the ring's area vectors, static gate (2026-11-24)

Route of the gate of 2026-11-23 (`d2_gateL.sbatch`): the exact initial state by initDavofFields, one step of deltaT/1000, reconstructionError, serial, one state per job, the rows one after the other, every davof key explicit, one binary for all rows (the build with `ringLinear` of commit 1638911: libVoF sha256 a4ef78aaa10a5b31, advectorVoF 616b35007350e317; also in etc/README-davof.md). Rows: plicRDF, ringSum (PURE), ringLinear with `ringLinearPasses 1` and 2 (`ringLinearCorrection true`, `ringLinearDiag true`). 3D sphere N 32/64/128; 2D circle N 64/128/256/512; standard and shifted centre ((0.013 0.007 0.003) in 3D; (0.013 0.007 0) in 2D, the y component lying along the circle's axis). Normal columns of error.dat with the corrected indices (LNormalDiff1, LNormalDiff2, LNormalDiffInf = columns 13 to 15; 1 - n.n_e L1 = LNormalDiff2^2/2, Linf = LNormalDiffInf^2/2). Cells above 1e-3: interface cells with 1 - n.n_e above 1e-3 (print-only key normalCellDiag). Orders: log2 of the ratio between N and 2N.

### 3D sphere, standard centre

| N | row | status | 1-n.n_e L1 | order | ratio to plicRDF | 1-n.n_e Linf | abs(n-n_e) L1 | order | ratio to plicRDF | abs(n-n_e) Linf | cells with 1-n.n_e > 1e-3 | interface cells | diagnostics (y = 0 by solvability / negative weight / max abs(W (q-qbar).y), per pass) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 32 | plicRDF | completed | 0.0004155 | 1.50 | 1.000 | 0.00258 | 0.02072 | 0.75 | 1.000 | 0.07183 | 72 | 416 |  |
| 64 | plicRDF | completed | 0.0001472 | 2.15 | 1.000 | 0.0009808 | 0.01229 | 1.06 | 1.000 | 0.04429 | 0 | 1760 |  |
| 128 | plicRDF | completed | 3.327e-05 |  | 1.000 | 0.0004302 | 0.005892 |  | 1.000 | 0.02933 | 0 | 6920 |  |
| 32 | ringSum (PURE) | completed | 0.0007891 | 1.38 | 1.899 | 0.004124 | 0.03184 | 0.84 | 1.537 | 0.09082 | 120 | 416 |  |
| 64 | ringSum (PURE) | completed | 0.0003035 | 2.09 | 2.061 | 0.002104 | 0.01777 | 1.04 | 1.446 | 0.06487 | 216 | 1760 |  |
| 128 | ringSum (PURE) | completed | 7.111e-05 |  | 2.137 | 0.0006757 | 0.008614 |  | 1.462 | 0.03676 | 0 | 6920 |  |
| 32 | ringLinear, 1 pass | completed | 3.342e-06 | 2.35 | 0.008 | 1.017e-05 | 0.002321 | 1.31 | 0.112 | 0.004511 | 0 | 416 | pass 1: 0 / 0 / 1.0280049192653888 |
| 64 | ringLinear, 1 pass | completed | 6.545e-07 | 1.66 | 0.004 | 5.418e-06 | 0.0009347 | 0.95 | 0.076 | 0.003292 | 0 | 1760 | pass 1: 0 / 96 / 1.6116058417518038 |
| 128 | ringLinear, 1 pass | completed | 2.071e-07 |  | 0.006 | 4.477e-06 | 0.0004831 |  | 0.082 | 0.002992 | 0 | 6920 | pass 1: 0 / 216 / 1.8585873322495188 |
| 32 | ringLinear, 2 passes | completed | 2.461e-06 | 2.36 | 0.006 | 1.487e-05 | 0.001877 | 1.27 | 0.091 | 0.005454 | 0 | 416 | pass 1: 0 / 0 / 1.0268797222548285; pass 2: 0 / 0 / 1.0437530922046536 |
| 64 | ringLinear, 2 passes | completed | 4.8e-07 | 1.38 | 0.003 | 5.153e-06 | 0.000776 | 0.74 | 0.063 | 0.00321 | 0 | 1760 | pass 1: 0 / 96 / 1.6112029176352618; pass 2: 0 / 96 / 1.629210155394492 |
| 128 | ringLinear, 2 passes | completed | 1.845e-07 |  | 0.006 | 3.349e-06 | 0.0004633 |  | 0.079 | 0.002588 | 0 | 6920 | pass 1: 0 / 216 / 1.8585802850665238; pass 2: 0 / 216 / 1.8527174698162774 |
| 32 | ringLinear, correction off (control) | completed | 0.0007891 |  | 1.899 | 0.004124 | 0.03184 |  | 1.537 | 0.09082 | 120 | 416 | pass 1: 0 / 0 / 0 |

### 3D sphere, shifted centre

| N | row | status | 1-n.n_e L1 | order | ratio to plicRDF | 1-n.n_e Linf | abs(n-n_e) L1 | order | ratio to plicRDF | abs(n-n_e) Linf | cells with 1-n.n_e > 1e-3 | interface cells | diagnostics (y = 0 by solvability / negative weight / max abs(W (q-qbar).y), per pass) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 32 | plicRDF | completed | 0.0006393 | 2.09 | 1.000 | 0.007782 | 0.02597 | 1.07 | 1.000 | 0.1248 | 87 | 435 |  |
| 64 | plicRDF | completed | 0.0001497 | 2.08 | 1.000 | 0.001765 | 0.01234 | 1.07 | 1.000 | 0.05941 | 49 | 1737 |  |
| 128 | plicRDF | completed | 3.534e-05 |  | 1.000 | 0.0004644 | 0.005871 |  | 1.000 | 0.03048 | 0 | 6943 |  |
| 32 | ringSum (PURE) | completed | 0.001143 | 1.99 | 1.788 | 0.01237 | 0.03876 | 1.13 | 1.492 | 0.1573 | 135 | 435 |  |
| 64 | ringSum (PURE) | completed | 0.0002881 | 1.99 | 1.925 | 0.003111 | 0.01769 | 1.06 | 1.433 | 0.07888 | 143 | 1737 |  |
| 128 | ringSum (PURE) | completed | 7.272e-05 |  | 2.058 | 0.000753 | 0.008512 |  | 1.450 | 0.03881 | 0 | 6943 |  |
| 32 | ringLinear, 1 pass | completed | 5.61e-05 | 5.34 | 0.088 | 0.009106 | 0.004065 | 1.95 | 0.157 | 0.135 | 3 | 435 | pass 1: 0 / 15 / 1.7478792998077439 |
| 64 | ringLinear, 1 pass | completed | 1.385e-06 | 2.01 | 0.009 | 0.0003747 | 0.001052 | 1.19 | 0.085 | 0.02738 | 0 | 1737 | pass 1: 0 / 80 / 1.9961440838204447 |
| 128 | ringLinear, 1 pass | completed | 3.447e-07 |  | 0.010 | 0.0004923 | 0.0004623 |  | 0.079 | 0.03138 | 0 | 6943 | pass 1: 0 / 349 / 1.8607886598691337 |
| 32 | ringLinear, 2 passes | completed | 6.191e-05 | 5.85 | 0.097 | 0.01415 | 0.003 | 1.71 | 0.116 | 0.1682 | 2 | 435 | pass 1: 0 / 15 / 1.7481619998907516; pass 2: 0 / 13 / 1.6197966314916661 |
| 64 | ringLinear, 2 passes | completed | 1.073e-06 | 1.58 | 0.007 | 0.0004453 | 0.0009151 | 1.06 | 0.074 | 0.02984 | 0 | 1737 | pass 1: 0 / 80 / 1.99626121061193; pass 2: 0 / 79 / 1.9916950883796816 |
| 128 | ringLinear, 2 passes | completed | 3.595e-07 |  | 0.010 | 0.0006606 | 0.0004379 |  | 0.075 | 0.03635 | 0 | 6943 | pass 1: 0 / 349 / 1.8608561732431306; pass 2: 0 / 350 / 1.8577316673399631 |

### 2D circle, standard centre

| N | row | status | 1-n.n_e L1 | order | ratio to plicRDF | 1-n.n_e Linf | abs(n-n_e) L1 | order | ratio to plicRDF | abs(n-n_e) Linf | cells with 1-n.n_e > 1e-3 | interface cells | diagnostics (y = 0 by solvability / negative weight / max abs(W (q-qbar).y), per pass) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 64 | plicRDF | completed | 4.336e-05 | 1.23 | 1.000 | 0.0002698 | 0.005417 | 0.34 | 1.000 | 0.02323 | 0 | 76 |  |
| 128 | plicRDF | completed | 1.847e-05 | 2.53 | 1.000 | 0.0001491 | 0.00428 | 1.33 | 1.000 | 0.01727 | 0 | 156 |  |
| 256 | plicRDF | completed | 3.206e-06 | 1.79 | 1.000 | 2.163e-05 | 0.001706 | 1.02 | 1.000 | 0.006577 | 0 | 308 |  |
| 512 | plicRDF | completed | 9.245e-07 |  | 1.000 | 1.238e-05 | 0.0008428 |  | 1.000 | 0.004975 | 0 | 612 |  |
| 64 | ringSum (PURE) | completed | 0.0001609 | 1.51 | 3.711 | 0.0009958 | 0.01161 | 0.65 | 2.143 | 0.04463 | 0 | 76 |  |
| 128 | ringSum (PURE) | completed | 5.663e-05 | 2.23 | 3.065 | 0.0003312 | 0.007414 | 1.16 | 1.732 | 0.02574 | 0 | 156 |  |
| 256 | ringSum (PURE) | completed | 1.208e-05 | 2.02 | 3.769 | 7.571e-05 | 0.00331 | 1.04 | 1.940 | 0.01231 | 0 | 308 |  |
| 512 | ringSum (PURE) | completed | 2.972e-06 |  | 3.215 | 1.949e-05 | 0.001614 |  | 1.914 | 0.006244 | 0 | 612 |  |
| 64 | ringLinear, 1 pass | completed | 1.021e-06 | 0.65 | 0.024 | 4.587e-06 | 0.001158 | 0.35 | 0.214 | 0.003029 | 0 | 76 | pass 1: 0 / 0 / 0.56797744479790169 |
| 128 | ringLinear, 1 pass | completed | 6.501e-07 | 2.30 | 0.035 | 3.68e-06 | 0.0009065 | 1.18 | 0.212 | 0.002713 | 0 | 156 | pass 1: 0 / 0 / 0.65559518804004413 |
| 256 | ringLinear, 1 pass | completed | 1.318e-07 | 2.44 | 0.041 | 9.521e-07 | 0.0003995 | 1.28 | 0.234 | 0.00138 | 0 | 308 | pass 1: 0 / 0 / 0.60861089734759388 |
| 512 | ringLinear, 1 pass | completed | 2.435e-08 |  | 0.026 | 2.276e-07 | 0.0001643 |  | 0.195 | 0.0006747 | 0 | 612 | pass 1: 0 / 0 / 0.7092396495108656 |
| 64 | ringLinear, 2 passes | completed | 1.152e-06 | 0.89 | 0.027 | 5.194e-06 | 0.00124 | 0.48 | 0.229 | 0.003223 | 0 | 76 | pass 1: 0 / 0 / 0.56810274752071588; pass 2: 0 / 0 / 0.55789102335870377 |
| 128 | ringLinear, 2 passes | completed | 6.219e-07 | 2.25 | 0.034 | 3.861e-06 | 0.0008886 | 1.15 | 0.208 | 0.002779 | 0 | 156 | pass 1: 0 / 0 / 0.65558840142360908; pass 2: 0 / 0 / 0.65688517728879514 |
| 256 | ringLinear, 2 passes | completed | 1.309e-07 | 2.43 | 0.041 | 9.681e-07 | 0.0004005 | 1.28 | 0.235 | 0.001391 | 0 | 308 | pass 1: 0 / 0 / 0.60860896062993597; pass 2: 0 / 0 / 0.60991606131242626 |
| 512 | ringLinear, 2 passes | completed | 2.429e-08 |  | 0.026 | 2.294e-07 | 0.0001646 |  | 0.195 | 0.0006773 | 0 | 612 | pass 1: 0 / 0 / 0.7092394213129255; pass 2: 0 / 0 / 0.70878581806860086 |

### 2D circle, shifted centre

| N | row | status | 1-n.n_e L1 | order | ratio to plicRDF | 1-n.n_e Linf | abs(n-n_e) L1 | order | ratio to plicRDF | abs(n-n_e) Linf | cells with 1-n.n_e > 1e-3 | interface cells | diagnostics (y = 0 by solvability / negative weight / max abs(W (q-qbar).y), per pass) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 64 | plicRDF | completed | 5.715e-05 | 1.91 | 1.000 | 0.0005224 | 0.006393 | 0.89 | 1.000 | 0.03232 | 0 | 76 |  |
| 128 | plicRDF | completed | 1.525e-05 | 1.82 | 1.000 | 0.0002203 | 0.00345 | 0.91 | 1.000 | 0.02099 | 0 | 154 |  |
| 256 | plicRDF | completed | 4.314e-06 | 1.93 | 1.000 | 5.572e-05 | 0.001832 | 0.95 | 1.000 | 0.01056 | 0 | 308 |  |
| 512 | plicRDF | completed | 1.131e-06 |  | 1.000 | 1.599e-05 | 0.0009478 |  | 1.000 | 0.005656 | 0 | 614 |  |
| 64 | ringSum (PURE) | completed | 0.0001656 | 1.75 | 2.897 | 0.001408 | 0.01244 | 0.89 | 1.945 | 0.05307 | 2 | 76 |  |
| 128 | ringSum (PURE) | completed | 4.94e-05 | 1.90 | 3.238 | 0.0004305 | 0.006688 | 0.99 | 1.938 | 0.02934 | 0 | 154 |  |
| 256 | ringSum (PURE) | completed | 1.327e-05 | 1.99 | 3.076 | 0.0001014 | 0.003357 | 1.04 | 1.833 | 0.01424 | 0 | 308 |  |
| 512 | ringSum (PURE) | completed | 3.33e-06 |  | 2.946 | 2.917e-05 | 0.001637 |  | 1.728 | 0.007638 | 0 | 614 |  |
| 64 | ringLinear, 1 pass | completed | 1.093e-06 | 1.43 | 0.019 | 6.232e-06 | 0.001089 | 0.61 | 0.170 | 0.00353 | 0 | 76 | pass 1: 0 / 0 / 0.70000926989489465 |
| 128 | ringLinear, 1 pass | completed | 4.045e-07 | 1.89 | 0.027 | 2.847e-06 | 0.0007123 | 1.02 | 0.206 | 0.002386 | 0 | 154 | pass 1: 0 / 0 / 0.78097782862092535 |
| 256 | ringLinear, 1 pass | completed | 1.089e-07 | 2.01 | 0.025 | 9.281e-07 | 0.0003517 | 1.02 | 0.192 | 0.001362 | 0 | 308 | pass 1: 0 / 0 / 0.75681769917506214 |
| 512 | ringLinear, 1 pass | completed | 2.71e-08 |  | 0.024 | 2.54e-07 | 0.0001732 |  | 0.183 | 0.0007128 | 0 | 614 | pass 1: 0 / 0 / 0.75394425689050037 |
| 64 | ringLinear, 2 passes | completed | 1.097e-06 | 1.46 | 0.019 | 5.982e-06 | 0.001121 | 0.65 | 0.175 | 0.003459 | 0 | 76 | pass 1: 0 / 0 / 0.70004226455299012; pass 2: 0 / 0 / 0.69334259057101943 |
| 128 | ringLinear, 2 passes | completed | 3.991e-07 | 1.86 | 0.026 | 2.843e-06 | 0.0007121 | 1.01 | 0.206 | 0.002385 | 0 | 154 | pass 1: 0 / 0 / 0.78096414538112457; pass 2: 0 / 0 / 0.78016166127537101 |
| 256 | ringLinear, 2 passes | completed | 1.103e-07 | 2.03 | 0.026 | 9.477e-07 | 0.0003539 | 1.03 | 0.193 | 0.001377 | 0 | 308 | pass 1: 0 / 0 / 0.75681764466802137; pass 2: 0 / 0 / 0.75663671824765033 |
| 512 | ringLinear, 2 passes | completed | 2.709e-08 |  | 0.024 | 2.556e-07 | 0.0001732 |  | 0.183 | 0.000715 | 0 | 614 | pass 1: 0 / 0 / 0.75394425650627395; pass 2: 0 / 0 / 0.75387161015067727 |

### Check: the rows plicRDF and ringSum against the gate table of 2026-11-23 (same states, serial)

The last line of error.dat (LNormalDiff1, 2, Inf) of rows (a) and (b) of the new build is identical to the digit to the gate of 2026-11-23 in 28 of 28 states/rows; differences: none.

### Check: 1 rank against 8 ranks, 3D sphere N 64, standard centre, ringLinear 1 pass (job 57155296)

| quantity | 1 rank | 8 ranks | relative difference |
|---|---|---|---|
| LNormalDiff1 (abs(n-n_e) L1) | 0.00093467016829570763 | 0.0011058832732624038 | 0.183 |
| LNormalDiff2 | 0.0011441428992225739 | 0.0022359622609785311 | 0.954 |
| LNormalDiffInf | 0.0032916701276763409 | 0.031654518447855107 | 8.62 |
| nInterface | 1760 | 1760 | 0 |
| 1-n.n_e L1 (LNormalDiff2^2/2) | 6.5453148692071848e-07 | 2.4997636162601123e-06 | 2.82 |

Agreement to 1e-12 relative: NO (largest relative difference 8.62). The columns LNormalInf_ and LNormal1_ of the 8-rank file are 0 (not reduced over the ranks).

### The pre-registered statements

- PL1 (ringLinear, 1 pass, 1-n.n_e L1 at most 0.5 times plicRDF's on the 3D sphere at N 64 and 128, both positions; falsifier: above 1.05 in any of the four): pos 0 N 64: 0.004; pos 0 N 128: 0.006; pos 1 N 64: 0.009; pos 1 N 128: 0.010. At most 0.5 in 4 of 4; above 1.05 in 0 of 4.
- PL2 (its order in 1-n.n_e L1 from N 64 to 128 at least 3, both positions; falsifier: below 2.5): pos 0: 1.66; pos 1: 2.01.
- PL3 (2D circle, ratio to plicRDF at most 0.5 at N 128 and 256, both positions; falsifier: above 1.05): pos 0 N 128: 0.035; pos 0 N 256: 0.041; pos 1 N 128: 0.027; pos 1 N 256: 0.025. At most 0.5 in 4 of 4; above 1.05 in 0 of 4.
- PL4 (two passes change the 3D L1 by less than 20 % against one pass): pos 0 N 32: -26.4 %; pos 0 N 64: -26.7 %; pos 0 N 128: -10.9 %; pos 1 N 32: +10.4 %; pos 1 N 64: -22.5 %; pos 1 N 128: +4.3 %. Less than 20 % in 3 of 6.

Verdicts, literally against the written statements and falsifiers:
- PL1: held (at most 0.5 in 4 of 4; falsifier, above 1.05 in any of the four: not met).
- PL2: false (order at least 3 at 0 of 2 positions; falsifier, below 2.5: met at 2 of 2 positions).
- PL3: held (at most 0.5 in 4 of 4; falsifier, above 1.05: not met).
- PL4: false (less than 20 % in 3 of 6 3D states).

### The GATE (at most 1.05 times plicRDF's 1-n.n_e L1 on the 3D sphere at N 64 and 128, both positions)

| candidate | pos 0 N 64 | pos 0 N 128 | pos 1 N 64 | pos 1 N 128 | passes |
|---|---|---|---|---|---|
| ringLinear, 1 pass | 0.004 | 0.006 | 0.009 | 0.010 | yes |
| ringLinear, 2 passes | 0.003 | 0.006 | 0.007 | 0.010 | yes |

GATE per point (ringLinear, 1 pass, ratio at most 1.05): pos 0 N 64: 0.004 pass; pos 0 N 128: 0.006 pass; pos 1 N 64: 0.009 pass; pos 1 N 128: 0.010 pass.


## Controls of the 1-rank against 8-rank difference of ringLinear (2026-11-24)

3D sphere N 64, standard centre, one step of deltaT/1000, binary libVoF sha256 a4ef78aaa10a5b31, advectorVoF 616b35007350e317. 8 ranks: the cases of the serial rows decomposed with `numberOfSubdomains 8; method scotch;` (the decomposeParDict of the case with no further coefficients; the same as job 57155296; scotch takes no n per direction), `mpirun -np 8 advectorVoF -parallel`; jobs 57168520 (A, B, E0 on 8 ranks), 57168521 (E0 on 1 rank), 57155296 (E1 on 8 ranks, reused). Only the reduced columns (LNormalDiff1, LNormalDiff2, LNormalDiffInf) and the interface-cell count are compared.

| row | quantity | 1 rank | 8 ranks | relative difference |
|---|---|---|---|---|
| H1: (a) plicRDF | LNormalDiff1 | 0.012289435996391918 | 0.012289435996391344 | 4.67e-14 |
| H1: (a) plicRDF | LNormalDiff2 | 0.017158795506465957 | 0.017158795506465273 | 3.98e-14 |
| H1: (a) plicRDF | LNormalDiffInf | 0.044289082677107842 | 0.044289082677123531 | 3.54e-13 |
| H1: (a) plicRDF | nInterface | 1760 | 1760 | 0 |
| H2: (b) ringSum | LNormalDiff1 | 0.017765275221217895 | 0.017883463150394185 | 0.00665 |
| H2: (b) ringSum | LNormalDiff2 | 0.024636078125038204 | 0.024737736716490093 | 0.00413 |
| H2: (b) ringSum | LNormalDiffInf | 0.064873212274940431 | 0.064873212274940431 | 0 |
| H2: (b) ringSum | nInterface | 1760 | 1760 | 0 |
| H3: ringLinear, correction off | LNormalDiff1 | 0.017765275221218353 | 0.017883463150394362 | 0.00665 |
| H3: ringLinear, correction off | LNormalDiff2 | 0.024636078125039276 | 0.024737736716490724 | 0.00413 |
| H3: ringLinear, correction off | LNormalDiffInf | 0.064873212274940431 | 0.064873212274940431 | 0 |
| H3: ringLinear, correction off | nInterface | 1760 | 1760 | 0 |
| H4: ringLinear, 1 pass, correction on | LNormalDiff1 | 0.00093467016829570763 | 0.0011058832732624038 | 0.183 |
| H4: ringLinear, 1 pass, correction on | LNormalDiff2 | 0.0011441428992225739 | 0.0022359622609785311 | 0.954 |
| H4: ringLinear, 1 pass, correction on | LNormalDiffInf | 0.0032916701276763409 | 0.031654518447855107 | 8.62 |
| H4: ringLinear, 1 pass, correction on | nInterface | 1760 | 1760 | 0 |

Diagnostic lines of ringLinear (ringLinearDiag; last call of the run; cells with y = 0 by solvability / cells with a negative weight / max |W (q - qbar).y|):

| row | 1 rank | 8 ranks |
|---|---|---|
| H3: ringLinear, correction off | 0 / 0 / 0 | 0 / 0 / 0 |
| H4: ringLinear, 1 pass, correction on | 0 / 96 / 1.6116058417518038 | 0 / 96 / 1.6582633429496771 |

Verdicts (agreement to 1e-12 relative in all four quantities):

- H1: (a) plicRDF: agrees (largest relative difference 3.54e-13).
- H2: (b) ringSum: does NOT agree (largest relative difference 0.00665).
- H3: ringLinear, correction off: does NOT agree (largest relative difference 0.00665).
- H4: ringLinear, 1 pass, correction on: does NOT agree (largest relative difference 8.62).


## Where the input of the ring sum differs between 1 rank and 8 ranks (2026-11-24)

State: 3D sphere N 64, standard centre, row ringSum (b), the 1-rank case `gateL/3D_N64_p0/B` and the 8-rank case of job 57168520 (`gateL/ctl8_B`, same decomposition: scotch, 8 subdomains), one step of deltaT/1000. No code change. Jobs 57168520 (reused), 57168976 (8 ranks, davofProcFaceRule former). Comparison script `d2_ctl8-cmp.py` (Python/numpy run on the cluster with /usr/bin/python3; it reads the ascii processor fields and maps them to the serial mesh through cellProcAddressing and faceProcAddressing; no reconstructPar).

### A. Workflow order

initDavofFields runs on the serial case, before the decomposition. `d2_gate.sbatch` (which builds the serial case): `(cd $base && blockMesh > log.blockMesh 2>&1 && initDavofFields > log.initDavofFields 2>&1; echo "init rc=$?"; ...)`; the rows copy `$base/system $base/constant $base/0` and run serial. The 8-rank cases (`d2_gateL-ctl.sbatch`, `d2_gateL-par.sbatch`, `d2_gateL-former.sbatch`) copy the serial row's `system constant 0`, then `foamDictionary -entry numberOfSubdomains -set 8 ...; foamDictionary -entry method -set scotch ...; decomposePar -force; mpirun -np 8 advectorVoF -parallel`. The decomposed initial state is therefore the decomposition of the serial exact state.

### B. Which reconstruction reconstructionError scores

The one made AFTER the advection of the step, by the function object itself. `areaFractionError.C`, in `execute`: `const bool last = (t + 0.5*runTime.deltaTValue() >= runTime.endTime().value()); if (reconstructionErrorEveryStep_ || last) { reconstructionSchemes& surf = advPtr->surf(); surf.reconstruct(true); ... recErrPtr_->calcError(funcDict, surf.centre(), surf.normal(), curv, true);`. This reconstruct uses alpha and alphaf as left by the step's advection. The solver's own reconstruction of the step is the first call of `isoAdvection::advect`: `// reconstruct the interface  surf_->reconstruct();` (isoAdvectionTemplates.C line 552), before the advection, from the initial face fractions; the fields written at the end time (areaNormalDavof, interfaceNormal, interfaceCentre) are those of that start-of-step reconstruction (`advectorVoF.C`: `runTime.write()` follows the advection; the function object's reconstruct is made at the execution after it).

### C. Field comparison, 1 rank against 8 ranks, row ringSum

Times: 0 (the initial state as decomposed) and 5.20833e-06 (the end of the step). Fields written: alpha.water, alphaf.water (times 0 and end); areaNormalDavof.water, interfaceNormal.water, interfaceCentre.water (end only: they are not written at time 0 because the first reconstruction happens in the step). ringSum is not written and cannot be written through the dictionaries; it is recomputed offline below from areaNormal over the serial mesh's point-neighbour ring and is NOT the solver's own field.

Decomposition: 8 subdomains; cells with a processor face: 37954 of 524288; cells in the point-neighbour ring of those but without a processor face: 45677; processor faces (serial internal faces between two subdomains): 24423 of 1552384 internal faces.

### Cell fields (serial against the processor fields mapped through cellProcAddressing; difference above 1e-12)

| field | time | cells differing | maximum difference | of those: with a processor face | only in the point-neighbour ring of one | elsewhere |
|---|---|---|---|---|---|---|
| alpha.water | 0 | 0 | 0 | 0 | 0 | 0 |
| alpha.water | 5.20833e-06 | 0 | 0 | 0 | 0 | 0 |
| areaNormalDavof.water | 0 | not written at time 0 | | | | |
| areaNormalDavof.water | 5.20833e-06 | 0 | 0 | 0 | 0 | 0 |
| interfaceNormal.water | 0 | not written at time 0 | | | | |
| interfaceNormal.water | 5.20833e-06 | 0 | 3.58e-17 | 0 | 0 | 0 |
| interfaceCentre.water | 0 | not written at time 0 | | | | |
| interfaceCentre.water | 5.20833e-06 | 0 | 8.52e-15 | 0 | 0 | 0 |
| ringSum recomputed offline from the written areaNormal over the serial mesh's point-neighbour ring (candidate cells; NOT the solver's field) | 5.20833e-06 | 0 of 1760 | 0 | 0 | 0 | 0 |

### Face fractions alphaf.water (serial against the processor fields mapped through faceProcAddressing)

| time | faces differing from serial (above 1e-12) | of those: processor faces | internal faces of a rank | physical boundary faces | max abs(delta alpha_f) on processor faces | processor faces whose two sides differ (above 1e-12) | max difference between the two sides |
|---|---|---|---|---|---|---|---|
| 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| 5.20833e-06 | 4 | 4 | 0 | 0 | 0.294 | 4 | 0.294 |

Processor faces seen from two subdomains: 24423 of 24423 (each holds a value on both sides).

### Area vectors recomputed offline from the end-time face fractions (the state scored by reconstructionError: the reconstruction made after the advection of the step)

Area vectors m_c differing above 1e-12: 4 cells (maximum difference 7.18e-05); of those with a processor face: 4; only in the point-neighbour ring of one: 0; elsewhere: 0. Cells touching the differing face fractions (the faces' owner and neighbour cells): 8.
Ring sums recomputed from these (candidate cells, alpha inside (1e-8, 1-1e-8)): differing above 1e-12: 40 of 1760 cells (maximum difference 7.18e-05); of those with a processor face: 27; only in the point-neighbour ring of one: 13; elsewhere: 0. Interface-cell normals (the ring sums' directions) differ in the same cells.

The written area vectors and interface fields belong to the reconstruction before the advection (B above), where the face fractions are still the identical initial ones. The last two paragraphs recompute the area vectors from the written end-time face fractions (the state the scored reconstruction is made from), with each rank using its own side of a processor face.

### D. 8 ranks with davofProcFaceRule former, row ringSum

| quantity | 1 rank | 8 ranks, symmetric | 8 ranks, former | rel. diff. former vs 1 rank | rel. diff. former vs symmetric |
|---|---|---|---|---|---|
| LNormalDiff1 | 0.017765275221217895 | 0.017883463150394185 | 0.017743439747296478 | 0.00123 | 0.00783 |
| LNormalDiff2 | 0.024636078125038204 | 0.024737736716490093 | 0.024626696770197458 | 0.000381 | 0.00449 |
| LNormalDiffInf | 0.064873212274940431 | 0.064873212274940431 | 0.064873212274940431 | 0 | 0 |
| nInterface | 1760 | 1760 | 1760 | 0 | 0 |

The same comparison of the fields for the former run (extra):

alphaf.water, former (serial against the processor fields mapped through faceProcAddressing)

| time | faces differing from serial (above 1e-12) | of those: processor faces | internal faces of a rank | physical boundary faces | max abs(delta alpha_f) on processor faces | processor faces whose two sides differ (above 1e-12) | max difference between the two sides |
|---|---|---|---|---|---|---|---|
| 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| 5.20833e-06 | 58 | 58 | 0 | 0 | 0.0157 | 58 | 0.0157 |

Processor faces seen from two subdomains: 24423 of 24423 (each holds a value on both sides).

### The hypothesis HM, literally

(HM: the area vectors m_c differ between 1 and 8 ranks in cells with a processor face because alpha_f on processor faces differs from the value the same face has as an internal face in serial.) See the numbers above: (i) the area vectors recomputed from the end-time face fractions differ in the cells counted in paragraph C (all of them cells with a processor face); (ii) the face fractions alphaf differ from serial only on processor faces (4 faces at the end time, none at time 0; no internal face of a rank and no physical boundary face differs), and the two sides of those 4 faces hold different values; (iii) at time 0 nothing differs, after the step 4 processor faces differ.

