## The curved standing matrix with the fork of this branch: plicRDF, PURE, ONERING, ENTERING (2026-11-23)

Cells: 3 flows x 3 N x OpenFOAM Co 0.2 and 0.5 (the cases' TARGET_ADVECTIVE_CO, index 0 and 1), hexahedral meshes, the ranks of the standing rows, advectorVoF of this branch (every davof key explicit, davofProcFaceRule symmetric), one repetition per case (wall clock: one run per row, rows run one after the other on one exclusive i01 node), through iso-time.sbatch of the 2026-11-09 matrix (copy `d2_iso-time-ent.sbatch`; the prepared cases and the scripts of leia ac5fd7f5 read from the swept clone, the outputs in `agents/forkdefault/iso/time`). Orders are log2 of the ratio of the end values between N and 2N.
Normal norms at the end of the run: |n - n_e| L1 and Linf (LNormalDiff1, LNormalDiffInf of reconstructionError) and 1 - n.n_e L1 and Linf (= LNormalDiff2^2/2 and LNormalDiffInf^2/2 for unit normals).

| flow | N | ranks | OpenFOAM Co label | solver max Courant | row | status | E_shape | order | L1(x) | order | abs(n-n_e) L1 | order | abs(n-n_e) Linf | 1-n.n_e L1 | 1-n.n_e Linf | mass error (%, max) | min alpha | max alpha - 1 | interface cells (last) | wall (s) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| shear2D | 64 | 1 | 0.2 | 0.200 | plicRDF | completed | 0.0003199 | 3.85 | 39.45 | -0.49 | 0.334 | 0.84 | 1.55 | 0.0972 | 1.2 | 1.31e-12 | -8.78e-18 | 2.22e-16 |  | 28 |
| shear2D | 128 | 1 | 0.2 | 0.200 | plicRDF | completed | 2.215e-05 | 3.16 | 55.36 | 0.61 | 0.186 | 1.79 | 2 | 0.0816 | 1.99 | 3.73e-12 | -2.06e-05 | 4.44e-16 |  | 117 |
| shear2D | 256 | 8 | 0.2 | 0.200 | plicRDF | completed | 2.482e-06 |  | 36.38 |  | 0.054 |  | 0.595 | 0.00658 | 0.177 | 2.51e-12 | -8.8e-05 | 5.78e-05 |  | 119 |
| shear2D | 64 | 1 | 0.2 | 0.200 | PURE | completed | 0.00017 | 3.16 | 39.85 | 1.24 | 0.228 | 1.20 | 0.816 | 0.0401 | 0.333 | 1.04e-12 | -5e-13 | 2.22e-16 | 88 | 12 |
| shear2D | 128 | 1 | 0.2 | 0.200 | PURE | completed | 1.906e-05 | 3.11 | 16.85 | -3.26 | 0.0993 | 1.26 | 0.428 | 0.0109 | 0.0917 | 3.5e-12 | -2.51e-05 | 2.22e-16 | 164 | 64 |
| shear2D | 256 | 8 | 0.2 | 0.200 | PURE | completed | 2.201e-06 |  | 161 |  | 0.0416 |  | 0.281 | 0.00275 | 0.0393 | 2.07e-12 | -9.69e-05 | 1.82e-08 | 319 | 75 |
| shear2D | 64 | 1 | 0.2 | 0.200 | ONERING | completed | 0.0002352 | 3.34 | 85.48 | 0.70 | 0.267 | 1.23 | 0.603 | 0.0484 | 0.182 | 7.66e-13 | -5.38e-14 | 2.22e-16 | 88 | 12 |
| shear2D | 128 | 1 | 0.2 | 0.200 | ONERING | completed | 2.326e-05 | 3.23 | 52.47 | 0.07 | 0.114 | 1.38 | 0.615 | 0.0155 | 0.189 | 3.43e-12 | -2.43e-05 | 0.000125 | 165 | 65 |
| shear2D | 256 | 8 | 0.2 | 0.200 | ONERING | completed | 2.48e-06 |  | 49.81 |  | 0.0439 |  | 0.307 | 0.00326 | 0.0471 | 2.2e-12 | -2.95e-05 | 2.29e-05 | 318 | 76 |
| shear2D | 64 | 1 | 0.2 | 0.200 | ENTERING | completed | 0.0001933 | 3.31 | 43.72 | 1.23 | 0.234 | 1.24 | 0.908 | 0.0404 | 0.412 | 8.66e-13 | -7.24e-14 | 4.44e-16 | 84 | 13 |
| shear2D | 128 | 1 | 0.2 | 0.200 | ENTERING | completed | 1.943e-05 | 3.03 | 18.58 | -0.51 | 0.0988 | 1.14 | 0.428 | 0.0106 | 0.0916 | 3.3e-12 | -2.44e-05 | 9.55e-05 | 165 | 69 |
| shear2D | 256 | 8 | 0.2 | 0.200 | ENTERING | completed | 2.381e-06 |  | 26.4 |  | 0.0447 |  | 0.364 | 0.00349 | 0.0661 | 2.02e-12 | -0.000101 | 3.54e-05 | 316 | 85 |
| shear2D | 64 | 1 | 0.5 | 0.499 | plicRDF | completed | 0.0002046 | 3.17 | 64.31 | -1.70 | 0.325 | 0.75 | 1.98 | 0.128 | 1.96 | 7.55e-13 | -0.000122 | 5.95e-05 |  | 13 |
| shear2D | 128 | 1 | 0.5 | 0.499 | plicRDF | completed | 2.275e-05 | 3.15 | 208.5 | -2.06 | 0.193 | 1.77 | 1.97 | 0.0772 | 1.95 | 3.12e-12 | -0.000448 | 0.000661 |  | 54 |
| shear2D | 256 | 8 | 0.5 | 0.500 | plicRDF | completed | 2.557e-06 |  | 867.7 |  | 0.0566 |  | 0.442 | 0.00476 | 0.0975 | 1.69e-12 | -0.000437 | 0.000859 |  | 54 |
| shear2D | 64 | 1 | 0.5 | 0.499 | PURE | completed | 0.0001309 | 2.72 | 60.34 | 0.47 | 0.18 | 0.82 | 1.29 | 0.0288 | 0.83 | 7.55e-13 | -0.000117 | 4.44e-16 | 87 | 6 |
| shear2D | 128 | 1 | 0.5 | 0.499 | PURE | completed | 1.986e-05 | 2.95 | 43.52 | -0.45 | 0.102 | 1.28 | 1.92 | 0.0213 | 1.85 | 2.65e-12 | -0.000457 | 2.29e-05 | 174 | 28 |
| shear2D | 256 | 8 | 0.5 | 0.500 | PURE | completed | 2.562e-06 |  | 59.48 |  | 0.0418 |  | 0.324 | 0.00277 | 0.0526 | 1.41e-12 | -0.00034 | 0.00033 | 350 | 41 |
| shear2D | 64 | 1 | 0.5 | 0.499 | ONERING | completed | 0.0002225 | 3.03 | 37.59 | 1.21 | 0.244 | 1.22 | 0.56 | 0.0382 | 0.157 | 9.99e-13 | -0.000145 | 2.55e-09 | 94 | 6 |
| shear2D | 128 | 1 | 0.5 | 0.499 | ONERING | completed | 2.718e-05 | 3.26 | 16.27 | -0.61 | 0.105 | 1.22 | 0.596 | 0.013 | 0.178 | 3.35e-12 | -0.000336 | 0.000154 | 174 | 28 |
| shear2D | 256 | 8 | 0.5 | 0.500 | ONERING | completed | 2.841e-06 |  | 24.86 |  | 0.045 |  | 0.404 | 0.00375 | 0.0817 | 1.77e-12 | -0.000259 | 0.000326 | 355 | 33 |
| shear2D | 64 | 1 | 0.5 | 0.499 | ENTERING | completed | 0.0001484 | 2.58 | 46.45 | 0.08 | 0.169 | 0.76 | 0.45 | 0.0207 | 0.101 | 7.55e-13 | -0.000139 | 4.44e-16 | 87 | 6 |
| shear2D | 128 | 1 | 0.5 | 0.499 | ENTERING | completed | 2.485e-05 | 3.24 | 43.8 | 0.04 | 0.0998 | 1.26 | 0.51 | 0.0125 | 0.13 | 3.12e-12 | -0.000487 | 0.00032 | 175 | 30 |
| shear2D | 256 | 8 | 0.5 | 0.500 | ENTERING | completed | 2.635e-06 |  | 42.5 |  | 0.0415 |  | 0.396 | 0.00317 | 0.0786 | 2.14e-12 | -0.000274 | 0.000792 | 351 | 36 |
| shear3D | 32 | 1 | 0.2 | 0.217 | plicRDF | completed | 0.004707 | 1.81 | 0.01765 | 1.87 | 0.293 | 1.05 | 0.736 | 0.0531 | 0.271 | 1.7e-08 | -1.76e-06 | 4.44e-16 |  | 29 |
| shear3D | 64 | 8 | 0.2 | 0.217 | plicRDF | completed | 0.001342 | 2.02 | 0.00482 | 2.05 | 0.141 | 1.42 | 0.361 | 0.0132 | 0.0651 | 3.71e-10 | -3.02e-05 | 4.44e-16 |  | 127 |
| shear3D | 128 | 64 | 0.2 | 0.217 | plicRDF | completed | 0.0003307 |  | 0.001167 |  | 0.053 |  | 0.177 | 0.00232 | 0.0157 | 1.33e-09 | -4.58e-05 | 5.78e-06 |  | 464 |
| shear3D | 32 | 1 | 0.2 | 0.217 | PURE | completed | 0.006285 | 1.74 | 0.0239 | 1.82 | 0.369 | 1.13 | 0.799 | 0.0819 | 0.319 | 5.93e-10 | -4.85e-06 | 2.22e-16 | 480 | 25 |
| shear3D | 64 | 8 | 0.2 | 0.217 | PURE | completed | 0.001879 | 1.90 | 0.006774 | 1.93 | 0.169 | 1.30 | 0.43 | 0.018 | 0.0925 | 5.89e-10 | -2.08e-05 | 4.44e-16 | 1824 | 72 |
| shear3D | 128 | 64 | 0.2 | 0.217 | PURE | completed | 0.0005035 |  | 0.001784 |  | 0.0688 |  | 0.224 | 0.00358 | 0.0251 | 2.24e-10 | -5.92e-05 | 7.09e-06 | 7128 | 276 |
| shear3D | 32 | 1 | 0.2 | 0.217 | ONERING | completed | 0.00629 | 1.74 | 0.0239 | 1.82 | 0.37 | 1.13 | 0.8 | 0.0819 | 0.32 | 2.23e-10 | -5.64e-06 | 4.44e-16 | 482 | 33 |
| shear3D | 64 | 8 | 0.2 | 0.217 | ONERING | completed | 0.001879 | 1.90 | 0.006776 | 1.93 | 0.169 | 1.30 | 0.43 | 0.018 | 0.0926 | 5.97e-10 | -1.99e-05 | 4.44e-16 | 1819 | 73 |
| shear3D | 128 | 64 | 0.2 | 0.217 | ONERING | completed | 0.0005035 |  | 0.001784 |  | 0.0687 |  | 0.224 | 0.00358 | 0.0251 | 9.67e-11 | -5.86e-05 | 5.35e-06 | 7131 | 280 |
| shear3D | 32 | 1 | 0.2 | 0.217 | ENTERING | completed | 0.006286 | 1.74 | 0.0239 | 1.82 | 0.369 | 1.13 | 0.799 | 0.0818 | 0.319 | 3.39e-10 | -4.69e-06 | 4.44e-16 | 480 | 29 |
| shear3D | 64 | 8 | 0.2 | 0.217 | ENTERING | completed | 0.001879 | 1.90 | 0.006778 | 1.92 | 0.169 | 1.30 | 0.43 | 0.018 | 0.0926 | 5.31e-10 | -1.98e-05 | 4.44e-16 | 1822 | 81 |
| shear3D | 128 | 64 | 0.2 | 0.217 | ENTERING | completed | 0.0005035 |  | 0.001785 |  | 0.0688 |  | 0.224 | 0.00358 | 0.0251 | 1.3e-10 | -4.92e-05 | 5.31e-06 | 7131 | 303 |
| shear3D | 32 | 1 | 0.5 | 0.541 | plicRDF | completed | 0.004602 | 1.81 | 0.01742 | 1.86 | 0.291 | 1.08 | 0.701 | 0.0514 | 0.246 | 3.56e-08 | -0.000141 | 1.67e-06 |  | 15 |
| shear3D | 64 | 8 | 0.5 | 0.543 | plicRDF | completed | 0.001312 | 2.00 | 0.004804 | 1.97 | 0.138 | 1.29 | 0.387 | 0.0125 | 0.0749 | 1.21e-07 | -0.000343 | 0.000287 |  | 60 |
| shear3D | 128 | 64 | 0.5 | 0.543 | plicRDF | completed | 0.0003288 |  | 0.001225 |  | 0.0564 |  | 0.285 | 0.00248 | 0.0405 | 1.83e-07 | -0.000949 | 0.000655 |  | 216 |
| shear3D | 32 | 1 | 0.5 | 0.541 | PURE | completed | 0.006241 | 1.76 | 0.02397 | 1.84 | 0.367 | 1.14 | 0.777 | 0.081 | 0.302 | 1.78e-08 | -0.000786 | 4.44e-16 | 490 | 13 |
| shear3D | 64 | 8 | 0.5 | 0.543 | PURE | completed | 0.001848 | 1.90 | 0.006691 | 1.91 | 0.166 | 1.30 | 0.423 | 0.0174 | 0.0894 | 1.95e-06 | -0.00147 | 2.51e-05 | 1879 | 42 |
| shear3D | 128 | 64 | 0.5 | 0.543 | PURE | completed | 0.000494 |  | 0.001776 |  | 0.0674 |  | 0.223 | 0.00343 | 0.0248 | 1.32e-07 | -0.000491 | 0.000122 | 7371 | 126 |
| shear3D | 32 | 1 | 0.5 | 0.541 | ONERING | completed | 0.006249 | 1.76 | 0.02406 | 1.85 | 0.366 | 1.14 | 0.777 | 0.0804 | 0.302 | 2.2e-08 | -0.000808 | 4.44e-16 | 491 | 13 |
| shear3D | 64 | 8 | 0.5 | 0.543 | ONERING | completed | 0.001848 | 1.90 | 0.006695 | 1.92 | 0.166 | 1.30 | 0.427 | 0.0173 | 0.0911 | 2.02e-06 | -0.0015 | 9.96e-05 | 1885 | 35 |
| shear3D | 128 | 64 | 0.5 | 0.543 | ONERING | completed | 0.0004941 |  | 0.001775 |  | 0.0673 |  | 0.223 | 0.00343 | 0.0249 | 1.36e-07 | -0.000488 | 0.000198 | 7354 | 128 |
| shear3D | 32 | 1 | 0.5 | 0.541 | ENTERING | completed | 0.006246 | 1.76 | 0.02407 | 1.85 | 0.368 | 1.15 | 0.776 | 0.0815 | 0.301 | 1.72e-08 | -0.000796 | 4.44e-16 | 488 | 15 |
| shear3D | 64 | 8 | 0.5 | 0.543 | ENTERING | completed | 0.001849 | 1.90 | 0.006678 | 1.91 | 0.166 | 1.31 | 0.427 | 0.0174 | 0.0911 | 1.91e-06 | -0.0014 | 2.52e-05 | 1889 | 38 |
| shear3D | 128 | 64 | 0.5 | 0.543 | ENTERING | completed | 0.0004941 |  | 0.001774 |  | 0.0672 |  | 0.223 | 0.00342 | 0.0249 | 1.6e-07 | -0.000476 | 0.000169 | 7353 | 137 |
| deformation3D | 32 | 1 | 0.2 | 0.197 | plicRDF | completed | 0.00793 | 1.37 | 0.07906 | 0.37 | 0.692 | 0.16 | 2 | 0.337 | 1.99 | 4.87e-12 | -1.15e-08 | 4.44e-16 |  | 63 |
| deformation3D | 64 | 8 | 0.2 | 0.198 | plicRDF | completed | 0.003072 | 2.16 | 0.06121 | 2.95 | 0.618 | 1.66 | 2 | 0.364 | 2 | 2.83e-12 | -5.93e-05 | 4.44e-16 |  | 194 |
| deformation3D | 128 | 64 | 0.2 | 0.198 | plicRDF | completed | 0.0006893 |  | 0.007923 |  | 0.196 |  | 2 | 0.0895 | 2 | 4.57e-12 | -2.44e-05 | 1.63e-06 |  | 861 |
| deformation3D | 32 | 1 | 0.2 | 0.197 | PURE | completed | 0.009977 | 1.34 | 0.1123 | 0.44 | 0.826 | 0.15 | 2 | 0.45 | 2 | 4.64e-12 | -0.000145 | 4.44e-16 | 887 | 49 |
| deformation3D | 64 | 8 | 0.2 | 0.198 | PURE | completed | 0.003936 | 2.37 | 0.08275 | 2.64 | 0.746 | 2.01 | 2 | 0.466 | 2 | 3.04e-12 | -0.000922 | 7.31e-06 | 3284 | 128 |
| deformation3D | 128 | 64 | 0.2 | 0.198 | PURE | completed | 0.000761 |  | 0.01324 |  | 0.185 |  | 2 | 0.0773 | 2 | 4.42e-12 | -0.00085 | 7.89e-07 | 7594 | 462 |
| deformation3D | 32 | 1 | 0.2 | 0.197 | ONERING | completed | 0.009889 | 1.32 | 0.1006 | -0.16 | 0.778 | -0.08 | 2 | 0.402 | 2 | 4.43e-12 | -2.63e-05 | 4.44e-16 | 765 | 47 |
| deformation3D | 64 | 8 | 0.2 | 0.198 | ONERING | completed | 0.003965 | 2.38 | 0.1124 | 2.06 | 0.824 | 1.57 | 2 | 0.54 | 2 | 3.22e-12 | -0.000214 | 4.27e-06 | 3731 | 131 |
| deformation3D | 128 | 64 | 0.2 | 0.198 | ONERING | completed | 0.0007608 |  | 0.02694 |  | 0.277 |  | 2 | 0.153 | 2 | 4.5e-12 | -0.000257 | 7.94e-07 | 8219 | 463 |
| deformation3D | 32 | 1 | 0.2 | 0.197 | ENTERING | completed | 0.009983 | 1.34 | 0.143 | -0.29 | 0.931 | -0.01 | 2 | 0.557 | 2 | 4.87e-12 | -0.000397 | 4.44e-16 | 1132 | 50 |
| deformation3D | 64 | 8 | 0.2 | 0.198 | ENTERING | completed | 0.003937 | 2.37 | 0.1746 | 0.03 | 0.94 | 0.25 | 2 | 0.646 | 2 | 3.22e-12 | -0.00182 | 7.48e-06 | 4950 | 138 |
| deformation3D | 128 | 64 | 0.2 | 0.198 | ENTERING | completed | 0.0007607 |  | 0.1716 |  | 0.791 |  | 2 | 0.56 | 2 | 4.46e-12 | -0.00285 | 7.88e-07 | 16401 | 499 |
| deformation3D | 32 | 1 | 0.5 | 0.492 | plicRDF | completed | 0.007777 | 1.38 | 0.08439 | 0.51 | 0.719 | 0.27 | 2 | 0.365 | 1.99 | 3.82e-12 | -0.00016 | 4.44e-16 |  | 28 |
| deformation3D | 64 | 8 | 0.5 | 0.495 | plicRDF | completed | 0.002979 | 2.11 | 0.05941 | 2.96 | 0.597 | 1.65 | 2 | 0.347 | 2 | 2.07e-12 | -0.000205 | 0.000309 |  | 83 |
| deformation3D | 128 | 64 | 0.5 | 0.496 | plicRDF | completed | 0.0006897 |  | 0.007658 |  | 0.19 |  | 2 | 0.0853 | 1.99 | 4.35e-12 | -0.00023 | 8.03e-05 |  | 366 |
| deformation3D | 32 | 1 | 0.5 | 0.492 | PURE | completed | 0.009566 | 1.31 | 0.1301 | 0.50 | 0.903 | 0.35 | 2 | 0.541 | 2 | 4.33e-12 | -0.000453 | 4.44e-16 | 947 | 21 |
| deformation3D | 64 | 8 | 0.5 | 0.495 | PURE | completed | 0.00386 | 2.35 | 0.09224 | 2.97 | 0.709 | 2.05 | 2 | 0.431 | 2 | 3.73e-12 | -0.00171 | 2e-05 | 3221 | 54 |
| deformation3D | 128 | 64 | 0.5 | 0.496 | PURE | completed | 0.0007554 |  | 0.01173 |  | 0.171 |  | 2 | 0.0641 | 2 | 4.42e-12 | -0.00546 | 0.000287 | 7733 | 196 |
| deformation3D | 32 | 1 | 0.5 | 0.492 | ONERING | completed | 0.009555 | 1.32 | 0.09622 | 0.22 | 0.788 | 0.19 | 2 | 0.417 | 2 | 4.31e-12 | -4.02e-05 | 4.44e-16 | 755 | 21 |
| deformation3D | 64 | 8 | 0.5 | 0.495 | ONERING | completed | 0.003836 | 2.34 | 0.08271 | 2.42 | 0.689 | 1.65 | 2 | 0.413 | 2 | 2.56e-12 | -0.00497 | 2.04e-05 | 3050 | 54 |
| deformation3D | 128 | 64 | 0.5 | 0.496 | ONERING | completed | 0.0007578 |  | 0.01543 |  | 0.219 |  | 2 | 0.103 | 2 | 4.49e-12 | -0.00229 | 0.000299 | 8019 | 199 |
| deformation3D | 32 | 1 | 0.5 | 0.492 | ENTERING | completed | 0.009565 | 1.28 | 0.1617 | 0.57 | 1.05 | 0.34 | 2 | 0.706 | 2 | 4.43e-12 | -0.0024 | 8.1e-05 | 1345 | 23 |
| deformation3D | 64 | 8 | 0.5 | 0.495 | ENTERING | completed | 0.003931 | 2.38 | 0.1091 | 0.11 | 0.831 | 0.41 | 2 | 0.541 | 2 | 2.61e-12 | -0.0103 | 2.02e-05 | 3905 | 58 |
| deformation3D | 128 | 64 | 0.5 | 0.496 | ENTERING | completed | 0.0007564 |  | 0.1012 |  | 0.624 |  | 2 | 0.425 | 2 | 4.62e-12 | -0.0115 | 0.000295 | 12753 | 210 |

### The pre-registered statements

- PE2 (ENTERING completes all 18 cases): 18 of 18 completed.
- PE3 (E_shape of ENTERING <= PURE in at least 14 of 18): 3 of 18.
  ENTERING / plicRDF between 0.9 and 1.5 in 13 of 18; the ratios: shear2D N64 co0.2 0.60; shear2D N64 co0.5 0.73; shear2D N128 co0.2 0.88; shear2D N128 co0.5 1.09; shear2D N256 co0.2 0.96; shear2D N256 co0.5 1.03; shear3D N32 co0.2 1.34; shear3D N32 co0.5 1.36; shear3D N64 co0.2 1.40; shear3D N64 co0.5 1.41; shear3D N128 co0.2 1.52; shear3D N128 co0.5 1.50; deformation3D N32 co0.2 1.26; deformation3D N32 co0.5 1.23; deformation3D N64 co0.2 1.28; deformation3D N64 co0.5 1.32; deformation3D N128 co0.2 1.10; deformation3D N128 co0.5 1.10
- PE4 (wall clock ENTERING <= 1.1 x plicRDF at equal ranks in every case): 18 of 18; the ratios: shear2D N64 co0.2 0.45; shear2D N64 co0.5 0.48; shear2D N128 co0.2 0.59; shear2D N128 co0.5 0.57; shear2D N256 co0.2 0.71; shear2D N256 co0.5 0.67; shear3D N32 co0.2 0.99; shear3D N32 co0.5 0.99; shear3D N64 co0.2 0.64; shear3D N64 co0.5 0.63; shear3D N128 co0.2 0.65; shear3D N128 co0.5 0.63; deformation3D N32 co0.2 0.80; deformation3D N32 co0.5 0.82; deformation3D N64 co0.2 0.71; deformation3D N64 co0.5 0.70; deformation3D N128 co0.2 0.58; deformation3D N128 co0.5 0.57
  ENTERING / PURE in E_shape (3 of 18 cells have a ratio below 1; 10 within 1e-3 of 1): shear2D N64 co0.2 1.1373; shear2D N64 co0.5 1.1337; shear2D N128 co0.2 1.0191; shear2D N128 co0.5 1.2509; shear2D N256 co0.2 1.0817; shear2D N256 co0.5 1.0284; shear3D N32 co0.2 1.0001; shear3D N32 co0.5 1.0008; shear3D N64 co0.2 1.0001; shear3D N64 co0.5 1.0007; shear3D N128 co0.2 0.9999; shear3D N128 co0.5 1.0002; deformation3D N32 co0.2 1.0005; deformation3D N32 co0.5 0.9999; deformation3D N64 co0.2 1.0004; deformation3D N64 co0.5 1.0185; deformation3D N128 co0.2 0.9996; deformation3D N128 co0.5 1.0014
- PE5 (ONERING completes all 18: 18 of 18; between PURE and ENTERING in E_shape in at least 12 of 18): 2 of 18.

### Courant numbers: the solver's own maximum against the target label

| flow | N | label 0.2 | label 0.5 |
|---|---|---|---|
| shear2D | 64 | 0.1996 | 0.4991 |
| shear2D | 128 | 0.1998 | 0.4995 |
| shear2D | 256 | 0.1998 | 0.4996 |
| shear3D | 32 | 0.2165 | 0.5413 |
| shear3D | 64 | 0.2170 | 0.5426 |
| shear3D | 128 | 0.2170 | 0.5426 |
| deformation3D | 32 | 0.1968 | 0.4920 |
| deformation3D | 64 | 0.1980 | 0.4950 |
| deformation3D | 128 | 0.1983 | 0.4956 |

