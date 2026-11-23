## Normal controls on the 3D shear (2026-11-24)

Normal norms from `reconstructionError` (`LNormalDiff1`/`LNormalDiffInf` = L1 / Linf of |n - n_e|; 1 - n.n_e: L1 = LNormalDiff2^2/2, Linf = LNormalDiffInf^2/2 for unit normals, same orientation). Orders: log2 of the ratio between N and 2N.

### D1: the normals on the exact initial state (shear3D, OpenFOAM Co label 0.5 case, one step of deltaT/1000; ranks 1/8/64)

| row | N | 1 - n.n_e L1 | order | 1 - n.n_e Linf | abs(n - n_e) L1 | order | abs(n - n_e) Linf | L1(alpha_f) | Linf(alpha_f) |
|---|---|---|---|---|---|---|---|---|---|
| plicRDF | 32 | 0.000415 | 1.50 | 0.00259 | 0.0207 | 0.75 | 0.0719 | n/a | n/a |
| plicRDF | 64 | 0.000147 | 2.06 | 0.000984 | 0.0123 | 1.03 | 0.0444 | n/a | n/a |
| plicRDF | 128 | 3.54e-05 |  | 0.000433 | 0.00602 |  | 0.0294 | n/a | n/a |
| PURE (ringSum) | 32 | 0.000788 | 1.38 | 0.00413 | 0.0318 | 0.84 | 0.0909 | 0.00696 | 0.0361 |
| PURE (ringSum) | 64 | 0.000303 | 2.03 | 0.00211 | 0.0178 | 1.01 | 0.0649 | 0.00344 | 0.0232 |
| PURE (ringSum) | 128 | 7.43e-05 |  | 0.000733 | 0.0088 |  | 0.0383 | 0.00186 | 0.0167 |

PN1 (the ringSum normal's 1 - n.n_e L1 within 15 % of plicRDF's at every N; falsifier beyond 30 %): ratios PURE / plicRDF: N 32: 1.898; N 64: 2.061; N 128: 2.099

### D2: L1 and Linf of the transported alpha_f against the exact sphere's face fractions, PURE (areaFractionError)

| N | OpenFOAM Co label | right after the initialisation L1 | Linf | at t = T L1 | order | Linf |
|---|---|---|---|---|---|---|
| 32 | 0.2 | (the 0.5 case row) |  | 0.259 | 0.61 | 1 |
| 64 | 0.2 | (the 0.5 case row) |  | 0.17 | 0.72 | 1 |
| 128 | 0.2 | (the 0.5 case row) |  | 0.103 |  | 1 |
| 32 | 0.5 | 0.00696 | 0.0361 | 0.244 | 0.60 | 1 |
| 64 | 0.5 | 0.00344 | 0.0232 | 0.161 | 0.77 | 1 |
| 128 | 0.5 | 0.00186 | 0.0167 | 0.0943 |  | 1 |

(The initialisation row is the D1 run of the label-0.5 case; the initial state does not depend on the Courant number.)

### D3: plicRDF's normal computed on PURE's final alpha (restart of the plicRDF case from PURE's fields at t = T, one step of deltaT/1000)

| N | OpenFOAM Co label | 1 - n.n_e L1: plicRDF on PURE's alpha | plicRDF's own run at T | PURE's own run at T | ratio to plicRDF's own | abs(n - n_e) L1: on PURE's alpha | plicRDF own | PURE own |
|---|---|---|---|---|---|---|---|---|
| 32 | 0.2 | 0.0868 | 0.0531 | 0.0819 | 1.633 | 0.381 | 0.293 | 0.369 |
| 64 | 0.2 | 0.0194 | 0.0132 | 0.018 | 1.472 | 0.175 | 0.141 | 0.169 |
| 128 | 0.2 | 0.00374 | 0.00232 | 0.00358 | 1.612 | 0.0699 | 0.053 | 0.0688 |
| 32 | 0.5 | 0.0873 | 0.0514 | 0.081 | 1.697 | 0.382 | 0.291 | 0.367 |
| 64 | 0.5 | 0.0185 | 0.0125 | 0.0174 | 1.482 | 0.172 | 0.138 | 0.166 |
| 128 | 0.5 | 0.00367 | 0.00248 | 0.00343 | 1.483 | 0.0701 | 0.0564 | 0.0674 |

PN3 (plicRDF's normal on PURE's final alpha within 20 % of plicRDF's own run; falsifier beyond 50 %): ratios: N 32 co0.2: 1.633; N 64 co0.2: 1.472; N 128 co0.2: 1.612; N 32 co0.5: 1.697; N 64 co0.5: 1.482; N 128 co0.5: 1.483

### Row (c): the PURE loop with plicRDF's normal (normalModel plicRDF, a control, not a method) against the rows (a), (b), (d): E_shape at the end

| N | OpenFOAM Co label | (a) plicRDF | (b) plicRDF ccs | (c) PURE with plicRDF's normal | (d) PURE | (c)/(a) | (c)/(b) | (d)/(a) | status (c) | 1 - n.n_e L1: (a) / (c) / (d) |
|---|---|---|---|---|---|---|---|---|---|---|
| 32 | 0.2 | 0.004707 | 0.004707 | 0.004707 | 0.006285 | 1.000 | 1.000 | 1.335 | completed | 0.0531 / 0.0531 / 0.0819 |
| 64 | 0.2 | 0.001342 | 0.001342 | 0.001342 | 0.001879 | 1.000 | 1.000 | 1.400 | completed | 0.0132 / 0.0132 / 0.018 |
| 128 | 0.2 | 0.0003307 | 0.0003307 | 0.0003307 | 0.0005035 | 1.000 | 1.000 | 1.523 | completed | 0.00232 / 0.00232 / 0.00358 |
| 32 | 0.5 | 0.004602 | 0.004601 | 0.004601 | 0.006241 | 1.000 | 1.000 | 1.356 | completed | 0.0514 / 0.0512 / 0.081 |
| 64 | 0.5 | 0.001312 | 0.001312 | 0.001312 | 0.001848 | 1.000 | 1.000 | 1.409 | completed | 0.0125 / 0.0126 / 0.0174 |
| 128 | 0.5 | 0.0003288 | 0.0003291 | 0.0003292 | 0.000494 | 1.001 | 1.000 | 1.502 | completed | 0.00248 / 0.00255 / 0.00343 |

PN2 (row (c) within 5 % of plicRDF's and of row (b) at every N and Co; falsifier a remaining gap above 15 %): (c)/(a) and (c)/(b): N 32 co0.2: 1.000 / 1.000; N 64 co0.2: 1.000 / 1.000; N 128 co0.2: 1.000 / 1.000; N 32 co0.5: 1.000 / 1.000; N 64 co0.5: 1.000 / 1.000; N 128 co0.5: 1.001 / 1.000

