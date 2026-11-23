## The ring quadric's outliers (2026-11-23)

The 3D sphere (shear3D initial state, R = 0.15), N 64, the standard centre and the centre shifted by (0.013 0.007 0.003); the static state of the gate (one step of deltaT/1000), rows plicRDF (a), ringSum (b), ringFit (c) (every davof key explicit); print only: key `normalCellDiag` of reconstructionError (one line per interface cell with 1 - n.n_e of the final plane against the exact normal at its centre) and key `ringFitDiag` of the model ringFit (one line per interface cell with alpha_k, |m_k|/A_k, A_k the sum of the cell's face areas, the ring planes entering the fit (ring interface cells with a plane and a positive weight), those among them with n.e3 > alignTol (the test of the fit's own drop; the fit's count after the drop is internal), the fallback flag). The fit computes no conditioning measure (mongeFit returns only the success of its Cholesky solve); none is added. Build: the branch with the two print-only keys (the Q build, `TwoPhaseFlow-Q`), one binary for all rows.

### shifted centre

Interface cells (plicRDF, ringSum, ringFit): 1737, 1737, 1737. 1 - n.n_e L1: plicRDF 0.0001497, ringSum 0.0002881, ringFit 0.0002595 (ratio to plicRDF 1.734). Cells with 1 - n.n_e above 1e-3: plicRDF 49, ringSum 143, ringFit 20. ringFit's L1 with the ten worst cells removed: 1.575e-05 (0.105 times plicRDF's). ringFit's fallbacks: 3 cells.

| rank | cell | alpha_k | abs(m_k)/A_k | ring planes entering the fit | of those n.e3 > alignTol | fit fell back | 1-n.n_e ringFit | 1-n.n_e ringSum | 1-n.n_e plicRDF |
|---|---|---|---|---|---|---|---|---|---|
| 1 | 40091 | 1 | 4.005e-06 | 14 | 11 | no | 0.2371 | 0.0008714 | 0.0007825 |
| 2 | 47715 | 1 | 5.047e-05 | 14 | 11 | no | 0.04962 | 0.001027 | 0.0006769 |
| 3 | 31905 | 0.99999 | 9.56e-05 | 13 | 9 | no | 0.03497 | 0.002287 | 0.001735 |
| 4 | 76952 | 1 | 7.333e-05 | 14 | 13 | no | 0.03162 | 0.002405 | 0.001565 |
| 5 | 56601 | 1 | 6.724e-06 | 14 | 14 | no | 0.02368 | 0.001673 | 0.0007796 |
| 6 | 55909 | 1 | 3.867e-05 | 14 | 14 | no | 0.01219 | 0.00108 | 0.0006886 |
| 7 | 36130 | 0.99999 | 0.0001793 | 14 | 14 | no | 0.01117 | 0.002309 | 0.0007584 |
| 8 | 44380 | 0.99999 | 8.742e-05 | 15 | 15 | no | 0.00924 | 0.0001051 | 8.529e-05 |
| 9 | 88922 | 1 | 1.682e-05 | 14 | 14 | no | 0.007823 | 0.0005362 | 0.0002911 |
| 10 | 88671 | 0.99998 | 0.000251 | 14 | 14 | no | 0.006138 | 0.001194 | 0.000748 |

### standard centre

Interface cells (plicRDF, ringSum, ringFit): 1760, 1760, 1760. 1 - n.n_e L1: plicRDF 0.0001472, ringSum 0.0003035, ringFit 9.213e-07 (ratio to plicRDF 0.006). Cells with 1 - n.n_e above 1e-3: plicRDF 0, ringSum 216, ringFit 0. ringFit's L1 with the ten worst cells removed: 8.192e-07 (0.006 times plicRDF's). ringFit's fallbacks: 0 cells.

| rank | cell | alpha_k | abs(m_k)/A_k | ring planes entering the fit | of those n.e3 > alignTol | fit fell back | 1-n.n_e ringFit | 1-n.n_e ringSum | 1-n.n_e plicRDF |
|---|---|---|---|---|---|---|---|---|---|
| 1 | 84903 | 0.99984 | 0.001797 | 14 | 14 | no | 2.015e-05 | 0.001687 | 0.0006505 |
| 2 | 97188 | 0.99984 | 0.001797 | 14 | 14 | no | 2.015e-05 | 0.001687 | 0.0006505 |
| 3 | 60711 | 0.99984 | 0.001797 | 14 | 14 | no | 2.015e-05 | 0.001687 | 0.0006505 |
| 4 | 97566 | 0.99984 | 0.001797 | 14 | 14 | no | 2.015e-05 | 0.001687 | 0.0006505 |
| 5 | 85470 | 0.99984 | 0.001797 | 14 | 14 | no | 2.015e-05 | 0.001687 | 0.0006505 |
| 6 | 60900 | 0.99984 | 0.001797 | 14 | 14 | no | 2.015e-05 | 0.001687 | 0.0006505 |
| 7 | 84510 | 0.99984 | 0.001102 | 14 | 14 | no | 1.673e-05 | 0.001796 | 0.0006505 |
| 8 | 60696 | 0.99984 | 0.001102 | 14 | 14 | no | 1.673e-05 | 0.001796 | 0.0006505 |
| 9 | 84888 | 0.99984 | 0.001102 | 14 | 14 | no | 1.673e-05 | 0.001796 | 0.0006505 |
| 10 | 35748 | 0.99984 | 0.001102 | 14 | 14 | no | 1.673e-05 | 0.001796 | 0.0006505 |

### The pre-registered statements

- HQ1 (at least eight of the ten worst cells of the shifted state have alpha_k outside [0.01, 0.99]; falsifier: fewer than five): 10 of 10 (alpha_k of the ten: 1, 1, 1, 1, 1, 1, 1, 1, 1, 1).
- HQ2 (without the ten worst cells ringFit's 1 - n.n_e L1 in the shifted state at most 0.1 times plicRDF's; falsifier: above 0.5): 1.575e-05 against plicRDF's 0.0001497, ratio 0.105 (with all cells: 1.734).

