## Ring-sum inconsistency of the ENTERING variant: tables (2026-11-23)

Case: plane translation, normal arm, U = n = (0.3 0.5 0.8124), OpenFOAM Co 0.5, advectorVoF (serial), the keys written explicitly:
`davofFluxExtension none; davofFaceFractionSource enteringFaces; boundSecondLayer false; sweptBoundedFaceFraction scale; interfaceCellRule alphaAndRingSum; normalModel ringSum; wispTol 0; davofProcFaceRule symmetric;`
(the 2026-11-09 stored case `swept/iso/gate-steps/N20_ENTERING_normal_co0.5` with n_cells, deltaT = 1e-4 x 20/N and endTime = 3 deltaT changed). Fork build of `davof-isoAdvector-ringsum`.

### Step 0: reproduction (job 57087084)

| N | run | check after step 1 | check after step 2 | check after step 3 | interface cells after step 1 / 2 | L1(alpha_f) step 1 / 2 / 3 | Linf(alpha_f) step 1 / 2 / 3 |
|---|---|---|---|---|---|---|---|
| 20 | stored 2026-11-09 | 0.02871 | 0.02653 | 0.0261 | 176 / 143 | 2.84e-08 / 0.000189 / 0.000278 | 1.05e-05 / 0.00613 / 0.0077 |
| 20 | reproduced | 0.02871 | 0.02653 | 0.0261 | 176 / 143 | 2.84e-08 / 0.000189 / 0.000278 | 1.05e-05 / 0.00613 / 0.0077 |
| 40 | reproduced | 0.02577 | 0.02381 | 0.02375 | 768 / 698 | 7.53e-15 / 0.000144 / 0.000208 | 4.87e-12 / 0.00565 / 0.0076 |
| 80 | reproduced | 0.02696 | 0.02491 | 0.02468 | 3223 / 3079 | 2.08e-14 / 0.000162 / 0.000232 | 3.04e-11 / 0.00585 / 0.00735 |

The check at N 20 after step 2 is 0.0287146920884307, the stored value to every digit; 176 interface cells after step 1 (the function object's own count, wispTol 1e-3: 173).
(The step numbering of the 2026-11-09 text, 'at step 2', is the reconstruction after the first step: the first reconstruction after the initial state, with the alpha and the face fractions of t^1.)

### Step 1: the counters at the reconstruction (key ringSumDiag; ENTERING; N 20 and N 40)

| N | reconstruction | interface cells | cells with check above 1e-3 / 1e-6 | faces with alpha_f != 0 not in the wet list | faces in the list with alpha_f == 0 | max abs(alpha_f memory - written) | max abs(alpha_f memory - exact) | cells interface under the active rule, not under normalOnly/wispTol 1e-3 (alpha range) | max abs(alpha_k - alpha of the exact plane) | cells with check above 1e-3 and abs(alpha_k - alpha_plane) above 1e-4 |
|---|---|---|---|---|---|---|---|---|---|---|
| 20 | after step 1 | 176 | 27 / 36 | 0 | 0 | 5.551115123125783e-17 | 1.045810946131986e-05 | 3 1.642216734990143e-06 .. 0.9999964579398068 | 0.001011137809491025 | 18 of 27 |
| 20 | after step 2 | 143 | 49 / 70 | 0 | 0 | 5.551115123125783e-17 | 0.006128055656963616 | 1 1.647159822505273e-06 .. 1.647159822505273e-06 | 0.03543709649790627 | 49 of 49 |
| 40 | after step 1 | 768 | 114 / 170 | 0 | 0 | 5.551115123125783e-17 | 4.874767256524137e-12 | 7 9.494087125631268e-07 .. 0.9999999445211679 | 0.0008596590841809526 | 88 of 114 |
| 40 | after step 2 | 698 | 228 / 354 | 0 | 0 | 5.551115123125783e-17 | 0.005653233422903914 | 6 9.504127009904551e-07 .. 0.9999999459952603 | 0.07311408829376409 | 222 of 228 |

### The ten cells with the largest check, after step 1 (N 20)

| cell | alpha_k | check | abs(m code - m memory all faces)/V^(2/3) | abs(m memory - m written)/V^(2/3) | interface under normalOnly/wispTol 1e-3 | alpha of the exact plane | alpha_k - alpha_plane | abs(1 - n.n_exact) | face sources (alphafSource ids) |
|---|---|---|---|---|---|---|---|---|---|
| 6758 | 0.9897763306 | 0.028715 | 8.336179128497505e-17 | 0 | yes | 0.9887651928 | 1.011e-03 | 7.238654120556021e-14 | 6 6 6 1 1 6 |
| 7137 | 0.9915477751 | 0.025294 | 5.495656578066151e-17 | 0 | yes | 0.990711835 | 8.359e-04 | 1.681987882307112e-13 | 6 6 6 1 1 6 |
| 7192 | 0.9915477751 | 0.025294 | 3.728052640787051e-17 | 0 | yes | 0.990711835 | 8.359e-04 | 2.220446049250313e-16 | 6 6 6 1 1 6 |
| 7571 | 0.9931017794 | 0.02209 | 8.909607841523557e-17 | 0 | yes | 0.9924195332 | 6.822e-04 | 0 | 6 6 6 1 1 6 |
| 7516 | 0.9931017794 | 0.02209 | 2.584939414228214e-17 | 0 | yes | 0.9924195332 | 6.822e-04 | 0 | 6 6 6 1 1 6 |
| 7950 | 0.9944525826 | 0.019103 | 1.550963648536925e-17 | 0 | yes | 0.9939039332 | 5.486e-04 | 0 | 6 6 6 1 1 6 |
| 7895 | 0.9944525826 | 0.019103 | 6.720842476993339e-17 | 0 | yes | 0.9939039332 | 5.486e-04 | 0 | 6 6 6 1 1 6 |
| 6757 | 0.8374391026 | 0.0085274 | 9.248160396692389e-17 | 0 | yes | 0.8384445929 | -1.005e-03 | 0 | 6 21 6 1 1 1 |
| 7136 | 0.8481634867 | 0.0073896 | 1.16981004080458e-16 | 0 | yes | 0.8489994268 | -8.359e-04 | 1.354472090042691e-13 | 6 6 6 1 1 1 |
| 7191 | 0.8481634867 | 0.0073896 | 8.271806125530273e-17 | 0 | yes | 0.8489994268 | -8.359e-04 | 0 | 6 6 6 1 1 1 |

alphafSource ids: 1 downwind face of a surface cell (its plane), 2 upwind face swept by this cell's plane, 3 and 4 boundary faces, 5 second ring, 6 entering sweep, 7 defaults, 8 processor sync, 9 entering sync.

### The dump of the worst cell, after step 1 (N 20, cell 6758)

```
RSDIAG cell 6758 centre (0.003700000000000001 0.0035 0.003300000000000001) alpha 0.9897763306305165 | interface: active rule yes, normalOnly/wispTol 1e-3 yes | |m_k|/V^(2/3) 0.1670421894481113 | check 0.0287146920884307
RSDIAG   m code (2.004512527466092e-09 3.340854212443675e-09 5.428219924378392e-09), m in-memory all faces (2.004512527466091e-09 3.340854212443673e-09 5.428219924378395e-09) (|diff|/V^(2/3) 8.336179128497505e-17), m written (2.004512527466091e-09 3.340854212443673e-09 5.428219924378395e-09) (|diff| 0), m exact faces (2.004512527466164e-09 3.34085421244362e-09 5.428219924378381e-09) (|code - exact| 2.280403019324392e-15, |memory - exact| 2.274370662254519e-15)
RSDIAG   ring sum code (1.160414364052351e-07 1.934026514707992e-07 3.142406281097547e-07), sum of m over the ring computed here (1.160414364052351e-07 1.934026514707992e-07 3.142406281097547e-07), ring cells 27: 6337(|m|/V23 0, alpha -1.568878360558301e-15) 6338(|m|/V23 0.03958990369391873, alpha 0.001296292566028933) 6339(|m|/V23 0.6110978471816431, alpha 0.08507816030630357) 6357(|m|/V23 0.3649273425999374, alpha 0.03627738461113813) 6358(|m|/V23 1.06387469877829, alpha 0.2552386476516439) 6359(|m|/V23 1.200691440199555, alpha 0.6124153192143386) 6377(|m|/V23 1.230912511372589, alpha 0.4901879938143126) 6378(|m|/V23 0.8957051605121893, alpha 0.827526392158569) 6379(|m|/V23 0.188204472609086, alpha 0.986563960479353) 6737(|m|/V23 1.0837752100677, alpha 0.268555422465824) 6738(|m|/V23 1.191326984532489, alpha 0.6272471835541068) 6739(|m|/V23 0.619819041044761, alpha 0.9126117840850511) 
RSDIAG   n used (0.3000005725795037 0.5000016199128711 0.8124026320344333), n exact (0.3000009360043805 0.5000015600073009 0.8124025346998625), |1 - n.n_exact| 7.238654120556021e-14, A (explicit plane) (1.659935147321834e-09 2.766558578869885e-09 4.495104378947705e-09) |A| 5.533099894387209e-09, CCS residual -6.058451752097401e-16, p explicit -0.0001245180340831883
RSDIAG   explicit plane: volume fraction 0.9915337879467099 against alpha_k 0.9897763306305165 (volume defect 0.001757457316193367); exact plane through x0(t) with the same n: p -0.0001208861211652807, volume fraction 0.9887651928210255, area vector (2.00451252746679e-09 3.340854212444836e-09 5.42821992438026e-09) |A| 6.681687577926742e-09, |m_k| 6.681687577924435e-09
RSDIAG   face 19617 memory 1 written 1 exact 1 inList 1 source 6 phi*sgn 1.199999999999996e-08
RSDIAG   face 19618 memory 1 written 1 exact 0.9999999999999999 inList 1 source 6 phi*sgn 1.999999999999987e-08
RSDIAG   face 19619 memory 1 written 1 exact 1 inList 1 source 6 phi*sgn 3.249600000000003e-08
RSDIAG   face 18459 memory 0.8642945018905397 written 0.8642945018905397 exact 0.86429450189054 inList 1 source 1 phi*sgn -3.249599999999985e-08
RSDIAG   face 19559 memory 0.9164786446889081 written 0.9164786446889081 exact 0.9164786446889092 inList 1 source 1 phi*sgn -1.999999999999997e-08
RSDIAG   face 19614 memory 0.9498871868133478 written 0.9498871868133478 exact 0.9498871868133459 inList 1 source 6 phi*sgn -1.200000000000002e-08
RSDIAG cell 7137 centre (0.0035 0.0033 0.0035) alpha 0.9915477751052071 | interface: active rule yes, normalOnly/wispTol 1e-3 yes | |m_k|/V^(2/3) 0.1471416781587094 | check 0.02529380487164045
```

### One key (group) at a time, N 20, Co 0.5, the same diagnostics

| variant (the ENTERING set with ...) | check after step 1 / 2 | cells with check above 1e-3 after step 1 | max abs(alpha_k - alpha_plane) after step 1 / 2 | L1(alpha_f) step 1 / 2 / 3 | job |
|---|---|---|---|---|---|
| nothing changed (ENTERING) | 0.0287 / 0.0265 | 27 | 0.00101 / 0.0354 | 2.84e-08 / 0.000189 / 0.000278 | 57088518 |
| davofFluxExtension secondRing | 0.0287 / 0.0265 | 27 | 0.00101 / 0.0354 | 2.84e-08 / 0.000189 / 0.000278 | 57088519 |
| davofFaceFractionSource secondRing | 0.0287 / 0.0281 | 27 | 0.00101 / 0.0354 | 2.04e-15 / 0.000194 / 0.000284 | 57088709 |
| boundSecondLayer true | 0.0287 / 0.0265 | 27 | 0.00101 / 0.0354 | 2.84e-08 / 0.000189 / 0.000278 | 57088710 |
| sweptBoundedFaceFraction none | 0.0287 / 0.0265 | 27 | 0.00101 / 0.0354 | 2.84e-08 / 0.000189 / 0.000278 | 57088711 |
| interfaceCellRule normalOnly, normalModel cellLocal, wispTol 1e-3 | 0.0323 / 0.0299 | 29 | 0.00121 / 0.0354 | 2.84e-08 / 0.000218 / 0.000703 | 57088712 |
| davofFluxExtension secondRing and davofFaceFractionSource secondRing | 2.47e-12 / 2.51e-12 | 0 | 2.4e-13 / 1.76e-12 | 2.04e-15 / 4.7e-15 / 6.02e-15 | 57088797 |
| as v5 and boundSecondLayer true | 2.47e-12 / 2.51e-12 | 0 | 2.4e-13 / 1.76e-12 | 2.04e-15 / 4.7e-15 / 6.02e-15 | 57088798 |
| ALL FORMER keys (explicit) | 8.82e-13 / 2.96e-12 | 0 | 2.7e-13 / 7.25e-13 | 2.03e-15 / 6.21e-15 / 3.39e-14 | 57088713 |
| ALL PURE keys (explicit, sweepOnly) | 0.622 / 0.622 | 43 | 1 / 1 | 0.007 / 0.00801 / 0.00772 | 57088714 |

