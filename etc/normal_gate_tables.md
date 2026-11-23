## The static normal gate (2026-11-23)

Route of D1: the exact initial state by initDavofFields (the shear3D sphere R = 0.15 at (0.5 0.75 0.25), the shear2D circle R = 0.15 at (0.5 0.75) in the x-z plane of the one-cell-thick mesh), one step of deltaT/1000 (deltaT of the study case at target Co 0.5, divided by 1000), reconstructionError at that time. Serial, one state per job, the rows one after the other, every davof key explicit, one binary for all rows. Shifted position: the centre moved by (0.013 0.007 0.003) in 3D and by (0.013 0.007 0) in 2D (the y component lies along the circle's axis and moves nothing: the 2D circle is shifted by 0.013 in x). 1 - n.n_e L1 = LNormalDiff2^2/2 and Linf = LNormalDiffInf^2/2 of reconstructionError (L1 and Linf of 1 - n.n_e); |n - n_e| L1 = LNormalDiff1, Linf = LNormalDiffInf. Orders: log2 of the ratio between N and 2N. Ratio: to row (a) plicRDF of the same state, position and N.

### 3D sphere, standard centre

| N | row | status | 1-n.n_e L1 | order | ratio to (a) | 1-n.n_e Linf | abs(n-n_e) L1 | order | ratio to (a) | abs(n-n_e) Linf | interface cells |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 32 | plicRDF | completed | 0.0004155 | 1.50 | 1.000 | 0.00258 | 0.02072 | 0.75 | 1.000 | 0.07183 | 416 |
| 64 | plicRDF | completed | 0.0001472 | 2.15 | 1.000 | 0.0009808 | 0.01229 | 1.06 | 1.000 | 0.04429 | 1760 |
| 128 | plicRDF | completed | 3.327e-05 |  | 1.000 | 0.0004302 | 0.005892 |  | 1.000 | 0.02933 | 6920 |
| 32 | ringSum (PURE) | completed | 0.0007891 | 1.38 | 1.899 | 0.004124 | 0.03184 | 0.84 | 1.537 | 0.09082 | 416 |
| 64 | ringSum (PURE) | completed | 0.0003035 | 2.09 | 2.061 | 0.002104 | 0.01777 | 1.04 | 1.446 | 0.06487 | 1760 |
| 128 | ringSum (PURE) | completed | 7.111e-05 |  | 2.137 | 0.0006757 | 0.008614 |  | 1.462 | 0.03676 | 6920 |
| 32 | ring quadric (ringFit) | completed | 9.261e-06 | 3.33 | 0.022 | 0.0001157 | 0.002763 | 1.57 | 0.133 | 0.01521 | 416 |
| 64 | ring quadric (ringFit) | completed | 9.213e-07 | -6.01 | 0.006 | 2.015e-05 | 0.0009289 | -1.07 | 0.076 | 0.006349 | 1760 |
| 128 | ring quadric (ringFit) | completed | 5.921e-05 |  | 1.779 | 0.01374 | 0.001955 |  | 0.332 | 0.1658 | 6920 |
| 32 | ringSumRDF, 0 iterations | completed | 0.0007891 |  | 1.899 | 0.004124 | 0.03184 |  | 1.537 | 0.09082 | 416 |
| 32 | ringSumRDF, 1 iteration | completed | 0.0004886 | 1.67 | 1.176 | 0.001697 | 0.02468 | 0.99 | 1.191 | 0.05826 | 416 |
| 64 | ringSumRDF, 1 iteration | completed | 0.000154 | 2.18 | 1.046 | 0.001129 | 0.01247 | 1.09 | 1.014 | 0.04752 | 1760 |
| 128 | ringSumRDF, 1 iteration | completed | 3.39e-05 |  | 1.019 | 0.0004435 | 0.005841 |  | 0.991 | 0.02978 | 6920 |
| 32 | ringSumRDF, 2 iterations | completed | 0.00047 | 1.68 | 1.131 | 0.001611 | 0.02435 | 0.99 | 1.176 | 0.05675 | 416 |
| 64 | ringSumRDF, 2 iterations | completed | 0.0001468 | 2.18 | 0.997 | 0.001066 | 0.01222 | 1.08 | 0.995 | 0.04618 | 1760 |
| 128 | ringSumRDF, 2 iterations | completed | 3.25e-05 |  | 0.977 | 0.0004137 | 0.005779 |  | 0.981 | 0.02876 | 6920 |

### 3D sphere, shifted centre

| N | row | status | 1-n.n_e L1 | order | ratio to (a) | 1-n.n_e Linf | abs(n-n_e) L1 | order | ratio to (a) | abs(n-n_e) Linf | interface cells |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 32 | plicRDF | completed | 0.0006393 | 2.09 | 1.000 | 0.007782 | 0.02597 | 1.07 | 1.000 | 0.1248 | 435 |
| 64 | plicRDF | completed | 0.0001497 | 2.08 | 1.000 | 0.001765 | 0.01234 | 1.07 | 1.000 | 0.05941 | 1737 |
| 128 | plicRDF | completed | 3.534e-05 |  | 1.000 | 0.0004644 | 0.005871 |  | 1.000 | 0.03048 | 6943 |
| 32 | ringSum (PURE) | completed | 0.001143 | 1.99 | 1.788 | 0.01237 | 0.03876 | 1.13 | 1.492 | 0.1573 | 435 |
| 64 | ringSum (PURE) | completed | 0.0002881 | 1.99 | 1.925 | 0.003111 | 0.01769 | 1.06 | 1.433 | 0.07888 | 1737 |
| 128 | ringSum (PURE) | completed | 7.272e-05 |  | 2.058 | 0.000753 | 0.008512 |  | 1.450 | 0.03881 | 6943 |
| 32 | ring quadric (ringFit) | completed | 0.0002332 | -0.15 | 0.365 | 0.0404 | 0.006473 | 1.16 | 0.249 | 0.2843 | 435 |
| 64 | ring quadric (ringFit) | completed | 0.0002595 | 4.34 | 1.734 | 0.2371 | 0.002887 | 2.07 | 0.234 | 0.6886 | 1737 |
| 128 | ring quadric (ringFit) | completed | 1.28e-05 |  | 0.362 | 0.02331 | 0.0006889 |  | 0.117 | 0.2159 | 6943 |
| 32 | ringSumRDF, 1 iteration | completed | 0.0006833 | 2.17 | 1.069 | 0.007005 | 0.02848 | 1.21 | 1.097 | 0.1184 | 435 |
| 64 | ringSumRDF, 1 iteration | completed | 0.0001523 | 2.07 | 1.018 | 0.002035 | 0.01228 | 1.07 | 0.995 | 0.06379 | 1737 |
| 128 | ringSumRDF, 1 iteration | completed | 3.62e-05 |  | 1.024 | 0.0005194 | 0.005863 |  | 0.999 | 0.03223 | 6943 |
| 32 | ringSumRDF, 2 iterations | completed | 0.0006617 | 2.18 | 1.035 | 0.006819 | 0.02821 | 1.22 | 1.086 | 0.1168 | 435 |
| 64 | ringSumRDF, 2 iterations | completed | 0.0001465 | 2.08 | 0.978 | 0.001927 | 0.01212 | 1.07 | 0.982 | 0.06208 | 1737 |
| 128 | ringSumRDF, 2 iterations | completed | 3.474e-05 |  | 0.983 | 0.0004921 | 0.005779 |  | 0.984 | 0.03137 | 6943 |

### 2D circle, standard centre

| N | row | status | 1-n.n_e L1 | order | ratio to (a) | 1-n.n_e Linf | abs(n-n_e) L1 | order | ratio to (a) | abs(n-n_e) Linf | interface cells |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 64 | plicRDF | completed | 4.336e-05 | 1.23 | 1.000 | 0.0002698 | 0.005417 | 0.34 | 1.000 | 0.02323 | 76 |
| 128 | plicRDF | completed | 1.847e-05 | 2.53 | 1.000 | 0.0001491 | 0.00428 | 1.33 | 1.000 | 0.01727 | 156 |
| 256 | plicRDF | completed | 3.206e-06 | 1.79 | 1.000 | 2.163e-05 | 0.001706 | 1.02 | 1.000 | 0.006577 | 308 |
| 512 | plicRDF | completed | 9.245e-07 |  | 1.000 | 1.238e-05 | 0.0008428 |  | 1.000 | 0.004975 | 612 |
| 64 | ringSum (PURE) | completed | 0.0001609 | 1.51 | 3.711 | 0.0009958 | 0.01161 | 0.65 | 2.143 | 0.04463 | 76 |
| 128 | ringSum (PURE) | completed | 5.663e-05 | 2.23 | 3.065 | 0.0003312 | 0.007414 | 1.16 | 1.732 | 0.02574 | 156 |
| 256 | ringSum (PURE) | completed | 1.208e-05 | 2.02 | 3.769 | 7.571e-05 | 0.00331 | 1.04 | 1.940 | 0.01231 | 308 |
| 512 | ringSum (PURE) | completed | 2.972e-06 |  | 3.215 | 1.949e-05 | 0.001614 |  | 1.914 | 0.006244 | 612 |
| 64 | ring quadric (ringFit) | completed | 0.0001609 | 1.51 | 3.711 | 0.0009958 | 0.01161 | 0.64 | 2.143 | 0.04463 | 76 |
| 128 | ring quadric (ringFit) | completed | 5.662e-05 | 2.23 | 3.065 | 0.0003312 | 0.007431 | 1.17 | 1.736 | 0.02574 | 156 |
| 256 | ring quadric (ringFit) | completed | 1.208e-05 | 2.02 | 3.769 | 7.571e-05 | 0.00331 | 1.04 | 1.940 | 0.01231 | 308 |
| 512 | ring quadric (ringFit) | completed | 2.972e-06 |  | 3.215 | 1.949e-05 | 0.001614 |  | 1.914 | 0.006244 | 612 |
| 64 | ringSumRDF, 1 iteration | completed | 4.808e-05 | 1.40 | 1.109 | 0.0003236 | 0.005663 | 0.40 | 1.045 | 0.02544 | 76 |
| 128 | ringSumRDF, 1 iteration | completed | 1.818e-05 | 2.38 | 0.984 | 0.0001357 | 0.004288 | 1.29 | 1.002 | 0.01648 | 156 |
| 256 | ringSumRDF, 1 iteration | completed | 3.49e-06 | 1.93 | 1.088 | 2.66e-05 | 0.001756 | 1.07 | 1.030 | 0.007293 | 308 |
| 512 | ringSumRDF, 1 iteration | completed | 9.148e-07 |  | 0.990 | 1.519e-05 | 0.0008348 |  | 0.990 | 0.005512 | 612 |
| 64 | ringSumRDF, 2 iterations | completed | 4.177e-05 | 1.34 | 0.963 | 0.0002697 | 0.005208 | 0.34 | 0.961 | 0.02322 | 76 |
| 128 | ringSumRDF, 2 iterations | completed | 1.655e-05 | 2.42 | 0.896 | 0.0001261 | 0.00412 | 1.31 | 0.963 | 0.01588 | 156 |
| 256 | ringSumRDF, 2 iterations | completed | 3.103e-06 | 1.92 | 0.968 | 2.377e-05 | 0.001662 | 1.07 | 0.974 | 0.006895 | 308 |
| 512 | ringSumRDF, 2 iterations | completed | 8.228e-07 |  | 0.890 | 1.45e-05 | 0.0007917 |  | 0.939 | 0.005386 | 612 |

### 2D circle, shifted centre

| N | row | status | 1-n.n_e L1 | order | ratio to (a) | 1-n.n_e Linf | abs(n-n_e) L1 | order | ratio to (a) | abs(n-n_e) Linf | interface cells |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 64 | plicRDF | completed | 5.715e-05 | 1.91 | 1.000 | 0.0005224 | 0.006393 | 0.89 | 1.000 | 0.03232 | 76 |
| 128 | plicRDF | completed | 1.525e-05 | 1.82 | 1.000 | 0.0002203 | 0.00345 | 0.91 | 1.000 | 0.02099 | 154 |
| 256 | plicRDF | completed | 4.314e-06 | 1.93 | 1.000 | 5.572e-05 | 0.001832 | 0.95 | 1.000 | 0.01056 | 308 |
| 512 | plicRDF | completed | 1.131e-06 |  | 1.000 | 1.599e-05 | 0.0009478 |  | 1.000 | 0.005656 | 614 |
| 64 | ringSum (PURE) | completed | 0.0001656 | 1.75 | 2.897 | 0.001408 | 0.01244 | 0.89 | 1.945 | 0.05307 | 76 |
| 128 | ringSum (PURE) | completed | 4.94e-05 | 1.90 | 3.238 | 0.0004305 | 0.006688 | 0.99 | 1.938 | 0.02934 | 154 |
| 256 | ringSum (PURE) | completed | 1.327e-05 | 1.99 | 3.076 | 0.0001014 | 0.003357 | 1.04 | 1.833 | 0.01424 | 308 |
| 512 | ringSum (PURE) | completed | 3.33e-06 |  | 2.946 | 2.917e-05 | 0.001637 |  | 1.728 | 0.007638 | 614 |
| 64 | ring quadric (ringFit) | completed | 0.0001656 | 1.75 | 2.897 | 0.001408 | 0.01244 | 0.89 | 1.945 | 0.05307 | 76 |
| 128 | ring quadric (ringFit) | completed | 4.94e-05 | 1.90 | 3.238 | 0.0004305 | 0.006688 | 0.99 | 1.938 | 0.02934 | 154 |
| 256 | ring quadric (ringFit) | completed | 1.327e-05 | 1.99 | 3.076 | 0.0001014 | 0.003357 | 1.04 | 1.833 | 0.01424 | 308 |
| 512 | ring quadric (ringFit) | completed | 3.33e-06 |  | 2.945 | 2.917e-05 | 0.001637 |  | 1.728 | 0.007638 | 614 |
| 64 | ringSumRDF, 1 iteration | completed | 5.808e-05 | 1.85 | 1.016 | 0.0005082 | 0.006543 | 0.89 | 1.023 | 0.03188 | 76 |
| 128 | ringSumRDF, 1 iteration | completed | 1.615e-05 | 1.87 | 1.059 | 0.0002249 | 0.003522 | 0.95 | 1.021 | 0.02121 | 154 |
| 256 | ringSumRDF, 1 iteration | completed | 4.423e-06 | 2.00 | 1.025 | 5.816e-05 | 0.001828 | 0.99 | 0.998 | 0.01079 | 308 |
| 512 | ringSumRDF, 1 iteration | completed | 1.109e-06 |  | 0.981 | 1.414e-05 | 0.0009185 |  | 0.969 | 0.005317 | 614 |
| 64 | ringSumRDF, 2 iterations | completed | 5.169e-05 | 1.84 | 0.904 | 0.0004411 | 0.006186 | 0.89 | 0.968 | 0.0297 | 76 |
| 128 | ringSumRDF, 2 iterations | completed | 1.448e-05 | 1.85 | 0.949 | 0.0002144 | 0.003349 | 0.93 | 0.971 | 0.02071 | 154 |
| 256 | ringSumRDF, 2 iterations | completed | 4.011e-06 | 2.00 | 0.930 | 5.389e-05 | 0.001752 | 0.99 | 0.957 | 0.01038 | 308 |
| 512 | ringSumRDF, 2 iterations | completed | 1.001e-06 |  | 0.885 | 1.305e-05 | 0.0008831 |  | 0.932 | 0.005109 | 614 |

### Check: rdfIterations 0 against ringSum

- 3D N 32 position 0: the last line of error.dat (all columns) is identical for ringSum and rdfIterations 0 (ringSum: 1-n.n_e L1 0.0007890892907721482, rdfIterations 0: 0.0007890892907721482).

### plicRDF's RDF iterations in the step (DebugSwitches plicRDF 1, row (a); the printed value is the loop index at the exit, iterations used = index + 1)

| state | position | N | printed 'iterations =' values (one per reconstruct call) |
|---|---|---|---|
| 3D sphere | 0 | 32 | 1, 1 |
| 3D sphere | 0 | 64 | 1, 1 |
| 3D sphere | 0 | 128 | 2, 1 |
| 3D sphere | 1 | 32 | 1, 1 |
| 3D sphere | 1 | 64 | 1, 1 |
| 3D sphere | 1 | 128 | 2, 1 |
| 2D circle | 0 | 64 | 2, 1 |
| 2D circle | 0 | 128 | 2, 1 |
| 2D circle | 0 | 256 | 2, 1 |
| 2D circle | 0 | 512 | 2, 1 |
| 2D circle | 1 | 64 | 2, 1 |
| 2D circle | 1 | 128 | 2, 1 |
| 2D circle | 1 | 256 | 2, 1 |
| 2D circle | 1 | 512 | 2, 1 |

### The pre-registered statements

- PG1 (ring quadric 1-n.n_e L1 at most 0.5 times plicRDF's on the 3D sphere at N 64 and 128, both positions; falsifier: above 1.0 in any of the four): ratios pos 0 N 64: 0.006; pos 0 N 128: 1.779; pos 1 N 64: 1.734; pos 1 N 128: 0.362. At most 0.5 in 2 of 4; above 1.0 in 2 of 4.
- PG2 (ringSumRDF with one iteration within 10 % of plicRDF's 1-n.n_e L1 at every N, both states, both positions; falsifier: beyond 25 % anywhere): ratios 3D sphere pos 0 N 32: 1.176; 3D sphere pos 0 N 64: 1.046; 3D sphere pos 0 N 128: 1.019; 3D sphere pos 1 N 32: 1.069; 3D sphere pos 1 N 64: 1.018; 3D sphere pos 1 N 128: 1.024; 2D circle pos 0 N 64: 1.109; 2D circle pos 0 N 128: 0.984; 2D circle pos 0 N 256: 1.088; 2D circle pos 0 N 512: 0.990; 2D circle pos 1 N 64: 1.016; 2D circle pos 1 N 128: 1.059; 2D circle pos 1 N 256: 1.025; 2D circle pos 1 N 512: 0.981. Within 10 % in 12 of 14; beyond 25 % in 0.
- PG3 (2D circle, ratio of |n - n_e| L1 of the ring-sum normal to plicRDF's at least 1.2 at N 128 and 256; falsifier: at most 1.0 at N 128 and 256): pos 0 N 128: 1.732; pos 0 N 256: 1.940; pos 1 N 128: 1.938; pos 1 N 256: 1.833. At least 1.2 in 4 of 4; at most 1.0 in 0 of 4.

### The GATE (a candidate passes if its 1-n.n_e L1 is at most 1.05 times plicRDF's on the 3D sphere at N 64 and 128 in both positions)

| candidate | pos 0 N 64 | pos 0 N 128 | pos 1 N 64 | pos 1 N 128 | passes |
|---|---|---|---|---|---|
| ringSum (PURE) | 2.061 | 2.137 | 1.925 | 2.058 | no |
| ring quadric (ringFit) | 0.006 | 1.779 | 1.734 | 0.362 | no |
| ringSumRDF, 1 iteration | 1.046 | 1.019 | 1.018 | 1.024 | yes |
| ringSumRDF, 2 iterations | 0.997 | 0.977 | 0.978 | 0.983 | yes |

