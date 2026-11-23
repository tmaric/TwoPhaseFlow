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


## The trace of the flux (2026-11-23, second task)

Cases as the first task (plane translation, advectorVoF, serial, end time 8e-4, the initial state from initDavofFields for every row, all keys explicit); `ringSumFluxTrace 1`. Face error = dVf - exact swept flux (exact: piecewise quadratic in time, Simpson between the vertex-passage times); V_k = the owner cell's volume; 'plane' = the cell holds a plane at t^n (interface cell with a non-zero normal).

### A: step 1, normal arm, Co 0.5, ENTERING and PURE, N 20 and N 40: the faces with an error above 1e-12 V_k, by class

| row | N | faces | upwind has no plane, downwind none | upwind none, downwind has a plane | upwind has a plane (any) | boundary faces | max abs(err)/V over faces whose upwind cell has a plane | max abs(err)/V | alpha_k deviation max (cell) | alpha^n deviation max |
|---|---|---|---|---|---|---|---|---|---|---|
| entering | 20 | 21 | 19 | 0 | 0 | 2 | 4.82e-15 | 0.001011 | 0.001011 (6758) | 4.829470157119431e-15 |
| entering | 40 | 89 | 84 | 0 | 0 | 5 | 1.49e-13 | 0.0008597 | 0.0008597 (54153) | 1.254552017826427e-14 |
| pure | 20 | 21 | 19 | 0 | 0 | 2 | 4.82e-15 | 0.001011 | 0.001011 (6758) | 4.829470157119431e-15 |
| pure | 40 | 89 | 84 | 0 | 0 | 5 | 1.49e-13 | 0.0008597 | 0.0008597 (54153) | 1.254552017826427e-14 |

Class columns: 'upwind has no plane, downwind none' and 'upwind none, downwind has a plane' count internal faces; 'upwind has a plane' sums the classes up and both; the boundary column all boundary faces (none of them has a plane upwind where non-zero).

#### ENTERING, N 20, step 1: the faces (all 21) and the ten cells with the largest deviation

```
 face 18555 dVf 1.199999999999984e-12 exact 1.199997197471603e-12 err/V 3.503160475670866e-07 up 6397 (1, no plane) down 6398 (no plane) source 0
 face 19614 dVf 1.199999990830972e-12 exact 1.191910888355088e-12 err/V 0.001011137809485464 up 6757 (0.9999999923591413, no plane) down 6758 (no plane) source 0
 face 19615 dVf 1.99999998471828e-12 exact 1.999999984798655e-12 err/V -1.004693248889476e-11 up 6757 (0.9999999923591413, no plane) down 6777 (no plane) source 0
 face 19616 dVf 3.249599975170278e-12 exact 3.249600000000013e-12 err/V -3.103716812260114e-09 up 6757 (0.9999999923591413, no plane) down 7157 (no plane) source 0
 face 19670 dVf 1.199999999999995e-12 exact 1.199999994527511e-12 err/V 6.840604184180416e-10 up 6776 (1, no plane) down 6777 (no plane) source 0
 face 19709 dVf 1.199999999999996e-12 exact 1.198961997277538e-12 err/V 0.0001297503403072261 up 6794 (1, no plane) down 6795 (no plane) source 0
 face 20659 dVf 1.200000000000005e-12 exact 1.199297062493576e-12 err/V 8.786718830358597e-05 up 7118 (1, no plane) down 7119 (no plane) source 0
 face 20712 dVf 1.199999999999996e-12 exact 1.193312479500679e-12 err/V 0.0008359400624146535 up 7136 (1, no plane) down 7137 (no plane) source 0
 face 20821 dVf 1.200000000000001e-12 exact 1.199297062493572e-12 err/V 8.786718830358641e-05 up 7173 (1, no plane) down 7174 (no plane) source 0
 face 20863 dVf 1.200000000000002e-12 exact 1.193312479500685e-12 err/V 0.0008359400624146989 up 7191 (1, no plane) down 7192 (no plane) source 0
 face 21757 dVf 1.2e-12 exact 1.199550935309236e-12 err/V 5.613308634551886e-05 up 7497 (1, no plane) down 7498 (no plane) source 0
 face 21810 dVf 1.200000000000001e-12 exact 1.194542029820824e-12 err/V 0.0006822462723971182 up 7515 (1, no plane) down 7516 (no plane) source 0
 face 21919 dVf 1.199999999999999e-12 exact 1.199550935309234e-12 err/V 5.613308634554416e-05 up 7552 (1, no plane) down 7553 (no plane) source 0
 face 21972 dVf 1.200000000000012e-12 exact 1.194542029820836e-12 err/V 0.0006822462723970884 up 7570 (1, no plane) down 7571 (no plane) source 0
 face 22579 dVf 1.200000000000001e-12 exact 1.199734880929248e-12 err/V 3.31398838441228e-05 up 7876 (1, no plane) down 7877 (no plane) source 0
 face 22614 dVf 1.200000000000004e-12 exact 1.195610804520222e-12 err/V 0.0005486494349727688 up 7894 (1, no plane) down 7895 (no plane) source 0
 face 22686 dVf 1.199999999999996e-12 exact 1.199734880929243e-12 err/V 3.31398838440976e-05 up 7931 (1, no plane) down 7932 (no plane) source 0
 face 22721 dVf 1.199999999999998e-12 exact 1.195610804520216e-12 err/V 0.0005486494349727436 up 7949 (1, no plane) down 7950 (no plane) source 0
 face 22787 dVf 1.200000000000001e-12 exact 1.199734880929248e-12 err/V 3.31398838441228e-05 up 7986 (1, no plane) down 7987 (no plane) source 0
 face 23536(b) dVf 1.199999999999996e-12 exact 1.198961997277538e-12 err/V 0.0001297503403072266 up 6739 (1, no plane) down -1 (no plane) source 0
 face 23591(b) dVf 1.199999999999998e-12 exact 1.195610804520216e-12 err/V 0.000548649434972746 up 7839 (1, no plane) down -1 (no plane) source 0
```

```
 cell 6758 alpha^n 1 alpha^{n+1} 0.9897763306305165 plane has at t^n no D 0.001011137809484919 E_final 0.001011137809485474 E_raw 0.001011137809485474 alpha^n - alpha_plane(t^n) 0 faces with errors: 19614(-0.001011137809485474)
 cell 6757 alpha^n 0.9999999923591413 alpha^{n+1} 0.8374391025542925 plane has at t^n no D -0.001011134695717786 E_final -0.001011134695718917 E_raw -0.001011134695718917 alpha^n - alpha_plane(t^n) -7.771561172376096e-16 faces with errors: 19614(0.001011137809485464) 19615(-1.004693248889476e-11) 19616(-3.103716812260114e-09)
 cell 7191 alpha^n 1 alpha^{n+1} 0.8481634867115421 plane has at t^n no D -0.0008359400624158297 E_final -0.0008359400624140425 E_raw -0.0008359400624140425 alpha^n - alpha_plane(t^n) 0 faces with errors: 20863(0.0008359400624146989)
 cell 7136 alpha^n 1 alpha^{n+1} 0.8481634867115425 plane has at t^n no D -0.0008359400624148305 E_final -0.0008359400624133408 E_raw -0.0008359400624133408 alpha^n - alpha_plane(t^n) 0 faces with errors: 20712(0.0008359400624146535)
 cell 7137 alpha^n 1 alpha^{n+1} 0.9915477751052071 plane has at t^n no D 0.0008359400624143865 E_final 0.0008359400624146983 E_raw 0.0008359400624146983 alpha^n - alpha_plane(t^n) 0 faces with errors: 20712(-0.0008359400624146477)
 cell 7192 alpha^n 1 alpha^{n+1} 0.9915477751052071 plane has at t^n no D 0.0008359400624142754 E_final 0.000835940062414594 E_raw 0.000835940062414594 alpha^n - alpha_plane(t^n) 0 faces with errors: 20863(-0.0008359400624146949)
 cell 7571 alpha^n 1 alpha^{n+1} 0.9931017794370602 plane has at t^n no D 0.0006822462723992118 E_final 0.0006822462723973967 E_raw 0.0006822462723973967 alpha^n - alpha_plane(t^n) 0 faces with errors: 21972(-0.0006822462723970938)
 cell 7516 alpha^n 1 alpha^{n+1} 0.9931017794370595 plane has at t^n no D 0.0006822462723982126 E_final 0.0006822462723969693 E_raw 0.0006822462723969693 alpha^n - alpha_plane(t^n) 0 faces with errors: 21810(-0.0006822462723971208)
 cell 7515 alpha^n 1 alpha^{n+1} 0.8584878448836384 plane has at t^n no D -0.0006822462723961031 E_final -0.0006822462723954268 E_raw -0.0006822462723954268 alpha^n - alpha_plane(t^n) 0 faces with errors: 21810(0.0006822462723971182)
 cell 7570 alpha^n 1 alpha^{n+1} 0.8584878448836365 plane has at t^n no D -0.0006822462723951039 E_final -0.0006822462723969369 E_raw -0.0006822462723969369 alpha^n - alpha_plane(t^n) 0 faces with errors: 21972(0.0006822462723970884)
```

#### ENTERING, N 40, step 1: the faces (all 89) and the ten cells with the largest deviation

```
 face 136875 dVf 1.499999999999998e-13 exact 1.499815227770142e-13 err/V 1.847722298558206e-05 up 46397 (1, no plane) down 46398 (no plane) source 0
 face 141510 dVf 1.499999999999995e-13 exact 1.499916732120401e-13 err/V 8.326787959435818e-06 up 47956 (1, no plane) down 47957 (no plane) source 0
 face 141589 dVf 1.500000000000008e-13 exact 1.496568073199672e-13 err/V 0.0003431926800335455 up 47994 (1, no plane) down 47995 (no plane) source 0
 face 145879 dVf 1.499999999999997e-13 exact 1.49739758495416e-13 err/V 0.0002602415045836272 up 49438 (1, no plane) down 49439 (no plane) source 0
 face 146108 dVf 1.499999999999995e-13 exact 1.499971965459945e-13 err/V 2.803454004957039e-06 up 49515 (1, no plane) down 49516 (no plane) source 0
 face 146221 dVf 1.499999999999995e-13 exact 1.497397584954158e-13 err/V 0.0002602415045836502 up 49553 (1, no plane) down 49554 (no plane) source 0
 face 150477 dVf 1.500000000000008e-13 exact 1.498081346672478e-13 err/V 0.000191865332752996 up 50997 (1, no plane) down 50998 (no plane) source 0
 face 150706 dVf 1.499999999999998e-13 exact 1.499995009294666e-13 err/V 4.990705331943143e-07 up 51074 (1, no plane) down 51075 (no plane) source 0
 face 150819 dVf 1.499999999999992e-13 exact 1.498081346672462e-13 err/V 0.0001918653327530006 up 51112 (1, no plane) down 51113 (no plane) source 0
 face 151019 dVf 1.499999999999995e-13 exact 1.499995009294664e-13 err/V 4.990705331185815e-07 up 51189 (1, no plane) down 51190 (no plane) source 0
 face 154962 dVf 1.500000000000027e-13 exact 1.499999945130476e-13 err/V 5.486955145801013e-09 up 52518 (1, no plane) down 52519 (no plane) source 0
 face 155075 dVf 1.499999999999995e-13 exact 1.498633439860469e-13 err/V 0.0001366560139526238 up 52556 (1, no plane) down 52557 (no plane) source 0
 face 155304 dVf 1.499999999999995e-13 exact 1.499999945130444e-13 err/V 5.486955120557453e-09 up 52633 (1, no plane) down 52634 (no plane) source 0
 face 155417 dVf 1.499999999999998e-13 exact 1.498633439860472e-13 err/V 0.0001366560139525976 up 52671 (1, no plane) down 52672 (no plane) source 0
 face 155646 dVf 1.499999999999989e-13 exact 1.499999945130437e-13 err/V 5.486955145801018e-09 up 52748 (1, no plane) down 52749 (no plane) source 0
 face 155733 dVf 1.499999999999995e-13 exact 1.498633439860469e-13 err/V 0.0001366560139526232 up 52786 (1, no plane) down 52787 (no plane) source 0
 face 159444 dVf 1.500000000000027e-13 exact 1.491403409158978e-13 err/V 0.0008596590841049487 up 54038 (1, no plane) down 54039 (no plane) source 0
 face 159673 dVf 1.500000000000041e-13 exact 1.499067946024105e-13 err/V 9.320539759359618e-05 up 54115 (1, no plane) down 54116 (no plane) source 0
 face 159786 dVf 1.500000000000001e-13 exact 1.491403409158952e-13 err/V 0.0008596590841050009 up 54153 (1, no plane) down 54154 (no plane) source 0
 face 160015 dVf 1.499999999999998e-13 exact 1.499067946024062e-13 err/V 9.320539759362134e-05 up 54230 (1, no plane) down 54231 (no plane) source 0
 face 160128 dVf 1.500000000000021e-13 exact 1.491403409158971e-13 err/V 0.0008596590841049393 up 54268 (1, no plane) down 54269 (no plane) source 0
 face 160357 dVf 1.499999999999995e-13 exact 1.499067946024059e-13 err/V 9.32053975935959e-05 up 54345 (1, no plane) down 54346 (no plane) source 0
 face 160447 dVf 1.500000000000001e-13 exact 1.491403409158952e-13 err/V 0.0008596590841048925 up 54383 (1, no plane) down 54384 (no plane) source 0
 face 164042 dVf 1.500000000000001e-13 exact 1.492970256830326e-13 err/V 0.0007029743169675103 up 55597 (1, no plane) down 55598 (no plane) source 0
 face 164271 dVf 1.500000000000043e-13 exact 1.499398946669173e-13 err/V 6.010533308703882e-05 up 55674 (1, no plane) down 55675 (no plane) source 0
 face 164384 dVf 1.500000000000008e-13 exact 1.492970256830333e-13 err/V 0.0007029743169674893 up 55712 (1, no plane) down 55713 (no plane) source 0
 face 164613 dVf 1.500000000000008e-13 exact 1.499398946669137e-13 err/V 6.01053330870645e-05 up 55789 (1, no plane) down 55790 (no plane) source 0
 face 164726 dVf 1.500000000000008e-13 exact 1.492970256830333e-13 err/V 0.0007029743169674884 up 55827 (1, no plane) down 55828 (no plane) source 0
 face 164955 dVf 1.500000000000001e-13 exact 1.499398946669131e-13 err/V 6.010533308703986e-05 up 55904 (1, no plane) down 55905 (no plane) source 0
 face 165068 dVf 1.500000000000008e-13 exact 1.492970256830332e-13 err/V 0.000702974316967563 up 55942 (1, no plane) down 55943 (no plane) source 0
 face 168527 dVf 1.499999999999964e-13 exact 1.499640523301524e-13 err/V 3.594766984400607e-05 up 57118 (1, no plane) down 57119 (no plane) source 0
 face 168640 dVf 1.499999999999967e-13 exact 1.494334119957666e-13 err/V 0.000566588004230035 up 57156 (1, no plane) down 57157 (no plane) source 0
 face 168869 dVf 1.500000000000015e-13 exact 1.499640523301575e-13 err/V 3.594766984398033e-05 up 57233 (1, no plane) down 57234 (no plane) source 0
 face 168982 dVf 1.499999999999995e-13 exact 1.494334119957695e-13 err/V 0.0005665880042300111 up 57271 (1, no plane) down 57272 (no plane) source 0
 face 169211 dVf 1.499999999999998e-13 exact 1.499640523301558e-13 err/V 3.594766984398037e-05 up 57348 (1, no plane) down 57349 (no plane) source 0
 face 169324 dVf 1.500000000000012e-13 exact 1.494334119957712e-13 err/V 0.0005665880042299765 up 57386 (1, no plane) down 57387 (no plane) source 0
 face 169553 dVf 1.499999999999995e-13 exact 1.499640523301555e-13 err/V 3.59476698439805e-05 up 57463 (1, no plane) down 57464 (no plane) source 0
 face 169666 dVf 1.499999999999989e-13 exact 1.494334119957688e-13 err/V 0.0005665880042300432 up 57501 (1, no plane) down 57502 (no plane) source 0
 face 169877 dVf 1.500000000000008e-13 exact 1.499640523301568e-13 err/V 3.59476698440306e-05 up 57578 (1, no plane) down 57579 (no plane) source 0
 face 173125 dVf 1.500000000000005e-13 exact 1.49980675742725e-13 err/V 1.932425727551565e-05 up 58677 (1, no plane) down 58678 (no plane) source 0
 face 173238 dVf 1.499999999999992e-13 exact 1.495509080046954e-13 err/V 0.0004490919953037188 up 58715 (1, no plane) down 58716 (no plane) source 0
 face 173467 dVf 1.500000000000015e-13 exact 1.499806757427259e-13 err/V 1.932425727551533e-05 up 58792 (1, no plane) down 58793 (no plane) source 0
 face 173580 dVf 1.499999999999965e-13 exact 1.495509080046928e-13 err/V 0.0004490919953036882 up 58830 (1, no plane) down 58831 (no plane) source 0
 face 173809 dVf 1.500000000000023e-13 exact 1.499806757427268e-13 err/V 1.932425727551572e-05 up 58907 (1, no plane) down 58908 (no plane) source 0
 face 173922 dVf 1.499999999999998e-13 exact 1.495509080046961e-13 err/V 0.0004490919953036901 up 58945 (1, no plane) down 58946 (no plane) source 0
 face 174151 dVf 1.499999999999995e-13 exact 1.49980675742724e-13 err/V 1.932425727554103e-05 up 59022 (1, no plane) down 59023 (no plane) source 0
 face 174264 dVf 1.500000000000014e-13 exact 1.495509080046977e-13 err/V 0.0004490919953037313 up 59060 (1, no plane) down 59061 (no plane) source 0
 face 174493 dVf 1.499999999999995e-13 exact 1.49980675742724e-13 err/V 1.932425727551565e-05 up 59137 (1, no plane) down 59138 (no plane) source 0
 face 174591 dVf 1.500000000000005e-13 exact 1.495509080046968e-13 err/V 0.0004490919953036882 up 59175 (1, no plane) down 59176 (no plane) source 0
 face 177723 dVf 1.499999999999967e-13 exact 1.49991173055204e-13 err/V 8.826944792693988e-06 up 60236 (1, no plane) down 60237 (no plane) source 0
 face 177836 dVf 1.500000000000023e-13 exact 1.496509218604028e-13 err/V 0.0003490781395994942 up 60274 (1, no plane) down 60275 (no plane) source 0
 face 178065 dVf 1.499999999999995e-13 exact 1.499911730552068e-13 err/V 8.826944792694002e-06 up 60351 (1, no plane) down 60352 (no plane) source 0
 face 178178 dVf 1.499999999999995e-13 exact 1.496509218603999e-13 err/V 0.0003490781395995463 up 60389 (1, no plane) down 60390 (no plane) source 0
 face 178407 dVf 1.500000000000023e-13 exact 1.499911730552096e-13 err/V 8.826944792693963e-06 up 60466 (1, no plane) down 60467 (no plane) source 0
 face 178520 dVf 1.499999999999995e-13 exact 1.496509218603999e-13 err/V 0.0003490781395995482 up 60504 (1, no plane) down 60505 (no plane) source 0
 face 178749 dVf 1.500000000000008e-13 exact 1.499911730552081e-13 err/V 8.826944792693961e-06 up 60581 (1, no plane) down 60582 (no plane) source 0
 face 178862 dVf 1.499999999999982e-13 exact 1.496509218603987e-13 err/V 0.0003490781395995266 up 60619 (1, no plane) down 60620 (no plane) source 0
 face 179091 dVf 1.499999999999989e-13 exact 1.499911730552062e-13 err/V 8.826944792694077e-06 up 60696 (1, no plane) down 60697 (no plane) source 0
 face 179204 dVf 1.499999999999975e-13 exact 1.49650921860398e-13 err/V 0.0003490781395995268 up 60734 (1, no plane) down 60735 (no plane) source 0
 face 182092 dVf 1.499999999999982e-13 exact 1.497348617134695e-13 err/V 0.0002651382865286355 up 61718 (1, no plane) down 61719 (no plane) source 0
... (29 more in the log)
```

```
 cell 54153 alpha^n 1 alpha^{n+1} 0.84665341946998 plane has at t^n no D -0.0008596590841105645 E_final -0.0008596590841034611 E_raw -0.0008596590841034611 alpha^n - alpha_plane(t^n) 0 faces with errors: 159786(0.0008596590841050009)
 cell 54383 alpha^n 1 alpha^{n+1} 0.8466534194699807 plane has at t^n no D -0.0008596590841088991 E_final -0.0008596590841042613 E_raw -0.0008596590841042613 alpha^n - alpha_plane(t^n) 0 faces with errors: 160447(0.0008596590841048925)
 cell 54154 alpha^n 1 alpha^{n+1} 0.9913079510859701 plane has at t^n no D 0.0008596590841059015 E_final 0.0008596590841063001 E_raw 0.0008596590841063001 alpha^n - alpha_plane(t^n) 0 faces with errors: 159786(-0.0008596590841049874)
 cell 54039 alpha^n 1 alpha^{n+1} 0.9913079510859699 plane has at t^n no D 0.0008596590841045693 E_final 0.0008596590841059042 E_raw 0.0008596590841059042 alpha^n - alpha_plane(t^n) 0 faces with errors: 159444(-0.0008596590841049449)
 cell 54384 alpha^n 1 alpha^{n+1} 0.9913079510859687 plane has at t^n no D 0.000859659084103126 E_final 0.0008596590841050429 E_raw 0.0008596590841050429 alpha^n - alpha_plane(t^n) 0 faces with errors: 160447(-0.0008596590841048915)
 cell 54038 alpha^n 1 alpha^{n+1} 0.8466534194699832 plane has at t^n no D -0.0008596590841013496 E_final -0.0008596590841020709 E_raw -0.0008596590841020709 alpha^n - alpha_plane(t^n) 0 faces with errors: 159444(0.0008596590841049487)
 cell 54269 alpha^n 1 alpha^{n+1} 0.9913079510859701 plane has at t^n no D 0.0008596590841010165 E_final 0.0008596590841064477 E_raw 0.0008596590841064477 alpha^n - alpha_plane(t^n) 0 faces with errors: 160128(-0.0008596590841049331)
 cell 54268 alpha^n 1 alpha^{n+1} 0.8466534194699847 plane has at t^n no D -0.0008596590841004614 E_final -0.0008596590840996634 E_raw -0.0008596590840996634 alpha^n - alpha_plane(t^n) 0 faces with errors: 160128(0.0008596590841049393)
 cell 55943 alpha^n 1 alpha^{n+1} 0.9928921973124363 plane has at t^n no D 0.0007029743169713054 E_final 0.0007029743169677135 E_raw 0.0007029743169677135 alpha^n - alpha_plane(t^n) 0 faces with errors: 165068(-0.000702974316967562)
 cell 55828 alpha^n 1 alpha^{n+1} 0.9928921973124371 plane has at t^n no D 0.0007029743169700842 E_final 0.0007029743169674357 E_raw 0.0007029743169674357 alpha^n - alpha_plane(t^n) 0 faces with errors: 164726(-0.0007029743169674862)
```

#### PURE, N 20, step 1: the faces (all 21) and the ten cells with the largest deviation

```
 face 18555 dVf 1.199999999999984e-12 exact 1.199997197471603e-12 err/V 3.503160475670866e-07 up 6397 (1, no plane) down 6398 (no plane) source 0
 face 19614 dVf 1.199999990830972e-12 exact 1.191910888355088e-12 err/V 0.001011137809485464 up 6757 (0.9999999923591413, no plane) down 6758 (no plane) source 0
 face 19615 dVf 1.99999998471828e-12 exact 1.999999984798655e-12 err/V -1.004693248889476e-11 up 6757 (0.9999999923591413, no plane) down 6777 (no plane) source 0
 face 19616 dVf 3.249599975170278e-12 exact 3.249600000000013e-12 err/V -3.103716812260114e-09 up 6757 (0.9999999923591413, no plane) down 7157 (no plane) source 0
 face 19670 dVf 1.199999999999995e-12 exact 1.199999994527511e-12 err/V 6.840604184180416e-10 up 6776 (1, no plane) down 6777 (no plane) source 0
 face 19709 dVf 1.199999999999996e-12 exact 1.198961997277538e-12 err/V 0.0001297503403072261 up 6794 (1, no plane) down 6795 (no plane) source 0
 face 20659 dVf 1.200000000000005e-12 exact 1.199297062493576e-12 err/V 8.786718830358597e-05 up 7118 (1, no plane) down 7119 (no plane) source 0
 face 20712 dVf 1.199999999999996e-12 exact 1.193312479500679e-12 err/V 0.0008359400624146535 up 7136 (1, no plane) down 7137 (no plane) source 0
 face 20821 dVf 1.200000000000001e-12 exact 1.199297062493572e-12 err/V 8.786718830358641e-05 up 7173 (1, no plane) down 7174 (no plane) source 0
 face 20863 dVf 1.200000000000002e-12 exact 1.193312479500685e-12 err/V 0.0008359400624146989 up 7191 (1, no plane) down 7192 (no plane) source 0
 face 21757 dVf 1.2e-12 exact 1.199550935309236e-12 err/V 5.613308634551886e-05 up 7497 (1, no plane) down 7498 (no plane) source 0
 face 21810 dVf 1.200000000000001e-12 exact 1.194542029820824e-12 err/V 0.0006822462723971182 up 7515 (1, no plane) down 7516 (no plane) source 0
 face 21919 dVf 1.199999999999999e-12 exact 1.199550935309234e-12 err/V 5.613308634554416e-05 up 7552 (1, no plane) down 7553 (no plane) source 0
 face 21972 dVf 1.200000000000012e-12 exact 1.194542029820836e-12 err/V 0.0006822462723970884 up 7570 (1, no plane) down 7571 (no plane) source 0
 face 22579 dVf 1.200000000000001e-12 exact 1.199734880929248e-12 err/V 3.31398838441228e-05 up 7876 (1, no plane) down 7877 (no plane) source 0
 face 22614 dVf 1.200000000000004e-12 exact 1.195610804520222e-12 err/V 0.0005486494349727688 up 7894 (1, no plane) down 7895 (no plane) source 0
 face 22686 dVf 1.199999999999996e-12 exact 1.199734880929243e-12 err/V 3.31398838440976e-05 up 7931 (1, no plane) down 7932 (no plane) source 0
 face 22721 dVf 1.199999999999998e-12 exact 1.195610804520216e-12 err/V 0.0005486494349727436 up 7949 (1, no plane) down 7950 (no plane) source 0
 face 22787 dVf 1.200000000000001e-12 exact 1.199734880929248e-12 err/V 3.31398838441228e-05 up 7986 (1, no plane) down 7987 (no plane) source 0
 face 23536(b) dVf 1.199999999999996e-12 exact 1.198961997277538e-12 err/V 0.0001297503403072266 up 6739 (1, no plane) down -1 (no plane) source 0
 face 23591(b) dVf 1.199999999999998e-12 exact 1.195610804520216e-12 err/V 0.000548649434972746 up 7839 (1, no plane) down -1 (no plane) source 0
```

```
 cell 6758 alpha^n 1 alpha^{n+1} 0.9897763306305165 plane has at t^n no D 0.001011137809484919 E_final 0.001011137809485474 E_raw 0.001011137809485474 alpha^n - alpha_plane(t^n) 0 faces with errors: 19614(-0.001011137809485474)
 cell 6757 alpha^n 0.9999999923591413 alpha^{n+1} 0.8374391025542925 plane has at t^n no D -0.001011134695717786 E_final -0.001011134695718917 E_raw -0.001011134695718917 alpha^n - alpha_plane(t^n) -7.771561172376096e-16 faces with errors: 19614(0.001011137809485464) 19615(-1.004693248889476e-11) 19616(-3.103716812260114e-09)
 cell 7191 alpha^n 1 alpha^{n+1} 0.8481634867115421 plane has at t^n no D -0.0008359400624158297 E_final -0.0008359400624140425 E_raw -0.0008359400624140425 alpha^n - alpha_plane(t^n) 0 faces with errors: 20863(0.0008359400624146989)
 cell 7136 alpha^n 1 alpha^{n+1} 0.8481634867115425 plane has at t^n no D -0.0008359400624148305 E_final -0.0008359400624133408 E_raw -0.0008359400624133408 alpha^n - alpha_plane(t^n) 0 faces with errors: 20712(0.0008359400624146535)
 cell 7137 alpha^n 1 alpha^{n+1} 0.9915477751052071 plane has at t^n no D 0.0008359400624143865 E_final 0.0008359400624146983 E_raw 0.0008359400624146983 alpha^n - alpha_plane(t^n) 0 faces with errors: 20712(-0.0008359400624146477)
 cell 7192 alpha^n 1 alpha^{n+1} 0.9915477751052071 plane has at t^n no D 0.0008359400624142754 E_final 0.000835940062414594 E_raw 0.000835940062414594 alpha^n - alpha_plane(t^n) 0 faces with errors: 20863(-0.0008359400624146949)
 cell 7571 alpha^n 1 alpha^{n+1} 0.9931017794370602 plane has at t^n no D 0.0006822462723992118 E_final 0.0006822462723973967 E_raw 0.0006822462723973967 alpha^n - alpha_plane(t^n) 0 faces with errors: 21972(-0.0006822462723970938)
 cell 7516 alpha^n 1 alpha^{n+1} 0.9931017794370595 plane has at t^n no D 0.0006822462723982126 E_final 0.0006822462723969693 E_raw 0.0006822462723969693 alpha^n - alpha_plane(t^n) 0 faces with errors: 21810(-0.0006822462723971208)
 cell 7515 alpha^n 1 alpha^{n+1} 0.8584878448836384 plane has at t^n no D -0.0006822462723961031 E_final -0.0006822462723954268 E_raw -0.0006822462723954268 alpha^n - alpha_plane(t^n) 0 faces with errors: 21810(0.0006822462723971182)
 cell 7570 alpha^n 1 alpha^{n+1} 0.8584878448836365 plane has at t^n no D -0.0006822462723951039 E_final -0.0006822462723969369 E_raw -0.0006822462723969369 alpha^n - alpha_plane(t^n) 0 faces with errors: 21972(0.0006822462723970884)
```

#### PURE, N 40, step 1: the faces (all 89) and the ten cells with the largest deviation

```
 face 136875 dVf 1.499999999999998e-13 exact 1.499815227770142e-13 err/V 1.847722298558206e-05 up 46397 (1, no plane) down 46398 (no plane) source 0
 face 141510 dVf 1.499999999999995e-13 exact 1.499916732120401e-13 err/V 8.326787959435818e-06 up 47956 (1, no plane) down 47957 (no plane) source 0
 face 141589 dVf 1.500000000000008e-13 exact 1.496568073199672e-13 err/V 0.0003431926800335455 up 47994 (1, no plane) down 47995 (no plane) source 0
 face 145879 dVf 1.499999999999997e-13 exact 1.49739758495416e-13 err/V 0.0002602415045836272 up 49438 (1, no plane) down 49439 (no plane) source 0
 face 146108 dVf 1.499999999999995e-13 exact 1.499971965459945e-13 err/V 2.803454004957039e-06 up 49515 (1, no plane) down 49516 (no plane) source 0
 face 146221 dVf 1.499999999999995e-13 exact 1.497397584954158e-13 err/V 0.0002602415045836502 up 49553 (1, no plane) down 49554 (no plane) source 0
 face 150477 dVf 1.500000000000008e-13 exact 1.498081346672478e-13 err/V 0.000191865332752996 up 50997 (1, no plane) down 50998 (no plane) source 0
 face 150706 dVf 1.499999999999998e-13 exact 1.499995009294666e-13 err/V 4.990705331943143e-07 up 51074 (1, no plane) down 51075 (no plane) source 0
 face 150819 dVf 1.499999999999992e-13 exact 1.498081346672462e-13 err/V 0.0001918653327530006 up 51112 (1, no plane) down 51113 (no plane) source 0
 face 151019 dVf 1.499999999999995e-13 exact 1.499995009294664e-13 err/V 4.990705331185815e-07 up 51189 (1, no plane) down 51190 (no plane) source 0
 face 154962 dVf 1.500000000000027e-13 exact 1.499999945130476e-13 err/V 5.486955145801013e-09 up 52518 (1, no plane) down 52519 (no plane) source 0
 face 155075 dVf 1.499999999999995e-13 exact 1.498633439860469e-13 err/V 0.0001366560139526238 up 52556 (1, no plane) down 52557 (no plane) source 0
 face 155304 dVf 1.499999999999995e-13 exact 1.499999945130444e-13 err/V 5.486955120557453e-09 up 52633 (1, no plane) down 52634 (no plane) source 0
 face 155417 dVf 1.499999999999998e-13 exact 1.498633439860472e-13 err/V 0.0001366560139525976 up 52671 (1, no plane) down 52672 (no plane) source 0
 face 155646 dVf 1.499999999999989e-13 exact 1.499999945130437e-13 err/V 5.486955145801018e-09 up 52748 (1, no plane) down 52749 (no plane) source 0
 face 155733 dVf 1.499999999999995e-13 exact 1.498633439860469e-13 err/V 0.0001366560139526232 up 52786 (1, no plane) down 52787 (no plane) source 0
 face 159444 dVf 1.500000000000027e-13 exact 1.491403409158978e-13 err/V 0.0008596590841049487 up 54038 (1, no plane) down 54039 (no plane) source 0
 face 159673 dVf 1.500000000000041e-13 exact 1.499067946024105e-13 err/V 9.320539759359618e-05 up 54115 (1, no plane) down 54116 (no plane) source 0
 face 159786 dVf 1.500000000000001e-13 exact 1.491403409158952e-13 err/V 0.0008596590841050009 up 54153 (1, no plane) down 54154 (no plane) source 0
 face 160015 dVf 1.499999999999998e-13 exact 1.499067946024062e-13 err/V 9.320539759362134e-05 up 54230 (1, no plane) down 54231 (no plane) source 0
 face 160128 dVf 1.500000000000021e-13 exact 1.491403409158971e-13 err/V 0.0008596590841049393 up 54268 (1, no plane) down 54269 (no plane) source 0
 face 160357 dVf 1.499999999999995e-13 exact 1.499067946024059e-13 err/V 9.32053975935959e-05 up 54345 (1, no plane) down 54346 (no plane) source 0
 face 160447 dVf 1.500000000000001e-13 exact 1.491403409158952e-13 err/V 0.0008596590841048925 up 54383 (1, no plane) down 54384 (no plane) source 0
 face 164042 dVf 1.500000000000001e-13 exact 1.492970256830326e-13 err/V 0.0007029743169675103 up 55597 (1, no plane) down 55598 (no plane) source 0
 face 164271 dVf 1.500000000000043e-13 exact 1.499398946669173e-13 err/V 6.010533308703882e-05 up 55674 (1, no plane) down 55675 (no plane) source 0
 face 164384 dVf 1.500000000000008e-13 exact 1.492970256830333e-13 err/V 0.0007029743169674893 up 55712 (1, no plane) down 55713 (no plane) source 0
 face 164613 dVf 1.500000000000008e-13 exact 1.499398946669137e-13 err/V 6.01053330870645e-05 up 55789 (1, no plane) down 55790 (no plane) source 0
 face 164726 dVf 1.500000000000008e-13 exact 1.492970256830333e-13 err/V 0.0007029743169674884 up 55827 (1, no plane) down 55828 (no plane) source 0
 face 164955 dVf 1.500000000000001e-13 exact 1.499398946669131e-13 err/V 6.010533308703986e-05 up 55904 (1, no plane) down 55905 (no plane) source 0
 face 165068 dVf 1.500000000000008e-13 exact 1.492970256830332e-13 err/V 0.000702974316967563 up 55942 (1, no plane) down 55943 (no plane) source 0
 face 168527 dVf 1.499999999999964e-13 exact 1.499640523301524e-13 err/V 3.594766984400607e-05 up 57118 (1, no plane) down 57119 (no plane) source 0
 face 168640 dVf 1.499999999999967e-13 exact 1.494334119957666e-13 err/V 0.000566588004230035 up 57156 (1, no plane) down 57157 (no plane) source 0
 face 168869 dVf 1.500000000000015e-13 exact 1.499640523301575e-13 err/V 3.594766984398033e-05 up 57233 (1, no plane) down 57234 (no plane) source 0
 face 168982 dVf 1.499999999999995e-13 exact 1.494334119957695e-13 err/V 0.0005665880042300111 up 57271 (1, no plane) down 57272 (no plane) source 0
 face 169211 dVf 1.499999999999998e-13 exact 1.499640523301558e-13 err/V 3.594766984398037e-05 up 57348 (1, no plane) down 57349 (no plane) source 0
 face 169324 dVf 1.500000000000012e-13 exact 1.494334119957712e-13 err/V 0.0005665880042299765 up 57386 (1, no plane) down 57387 (no plane) source 0
 face 169553 dVf 1.499999999999995e-13 exact 1.499640523301555e-13 err/V 3.59476698439805e-05 up 57463 (1, no plane) down 57464 (no plane) source 0
 face 169666 dVf 1.499999999999989e-13 exact 1.494334119957688e-13 err/V 0.0005665880042300432 up 57501 (1, no plane) down 57502 (no plane) source 0
 face 169877 dVf 1.500000000000008e-13 exact 1.499640523301568e-13 err/V 3.59476698440306e-05 up 57578 (1, no plane) down 57579 (no plane) source 0
 face 173125 dVf 1.500000000000005e-13 exact 1.49980675742725e-13 err/V 1.932425727551565e-05 up 58677 (1, no plane) down 58678 (no plane) source 0
 face 173238 dVf 1.499999999999992e-13 exact 1.495509080046954e-13 err/V 0.0004490919953037188 up 58715 (1, no plane) down 58716 (no plane) source 0
 face 173467 dVf 1.500000000000015e-13 exact 1.499806757427259e-13 err/V 1.932425727551533e-05 up 58792 (1, no plane) down 58793 (no plane) source 0
 face 173580 dVf 1.499999999999965e-13 exact 1.495509080046928e-13 err/V 0.0004490919953036882 up 58830 (1, no plane) down 58831 (no plane) source 0
 face 173809 dVf 1.500000000000023e-13 exact 1.499806757427268e-13 err/V 1.932425727551572e-05 up 58907 (1, no plane) down 58908 (no plane) source 0
 face 173922 dVf 1.499999999999998e-13 exact 1.495509080046961e-13 err/V 0.0004490919953036901 up 58945 (1, no plane) down 58946 (no plane) source 0
 face 174151 dVf 1.499999999999995e-13 exact 1.49980675742724e-13 err/V 1.932425727554103e-05 up 59022 (1, no plane) down 59023 (no plane) source 0
 face 174264 dVf 1.500000000000014e-13 exact 1.495509080046977e-13 err/V 0.0004490919953037313 up 59060 (1, no plane) down 59061 (no plane) source 0
 face 174493 dVf 1.499999999999995e-13 exact 1.49980675742724e-13 err/V 1.932425727551565e-05 up 59137 (1, no plane) down 59138 (no plane) source 0
 face 174591 dVf 1.500000000000005e-13 exact 1.495509080046968e-13 err/V 0.0004490919953036882 up 59175 (1, no plane) down 59176 (no plane) source 0
 face 177723 dVf 1.499999999999967e-13 exact 1.49991173055204e-13 err/V 8.826944792693988e-06 up 60236 (1, no plane) down 60237 (no plane) source 0
 face 177836 dVf 1.500000000000023e-13 exact 1.496509218604028e-13 err/V 0.0003490781395994942 up 60274 (1, no plane) down 60275 (no plane) source 0
 face 178065 dVf 1.499999999999995e-13 exact 1.499911730552068e-13 err/V 8.826944792694002e-06 up 60351 (1, no plane) down 60352 (no plane) source 0
 face 178178 dVf 1.499999999999995e-13 exact 1.496509218603999e-13 err/V 0.0003490781395995463 up 60389 (1, no plane) down 60390 (no plane) source 0
 face 178407 dVf 1.500000000000023e-13 exact 1.499911730552096e-13 err/V 8.826944792693963e-06 up 60466 (1, no plane) down 60467 (no plane) source 0
 face 178520 dVf 1.499999999999995e-13 exact 1.496509218603999e-13 err/V 0.0003490781395995482 up 60504 (1, no plane) down 60505 (no plane) source 0
 face 178749 dVf 1.500000000000008e-13 exact 1.499911730552081e-13 err/V 8.826944792693961e-06 up 60581 (1, no plane) down 60582 (no plane) source 0
 face 178862 dVf 1.499999999999982e-13 exact 1.496509218603987e-13 err/V 0.0003490781395995266 up 60619 (1, no plane) down 60620 (no plane) source 0
 face 179091 dVf 1.499999999999989e-13 exact 1.499911730552062e-13 err/V 8.826944792694077e-06 up 60696 (1, no plane) down 60697 (no plane) source 0
 face 179204 dVf 1.499999999999975e-13 exact 1.49650921860398e-13 err/V 0.0003490781395995268 up 60734 (1, no plane) down 60735 (no plane) source 0
 face 182092 dVf 1.499999999999982e-13 exact 1.497348617134695e-13 err/V 0.0002651382865286355 up 61718 (1, no plane) down 61719 (no plane) source 0
... (29 more in the log)
```

```
 cell 54153 alpha^n 1 alpha^{n+1} 0.84665341946998 plane has at t^n no D -0.0008596590841105645 E_final -0.0008596590841034611 E_raw -0.0008596590841034611 alpha^n - alpha_plane(t^n) 0 faces with errors: 159786(0.0008596590841050009)
 cell 54383 alpha^n 1 alpha^{n+1} 0.8466534194699807 plane has at t^n no D -0.0008596590841088991 E_final -0.0008596590841042613 E_raw -0.0008596590841042613 alpha^n - alpha_plane(t^n) 0 faces with errors: 160447(0.0008596590841048925)
 cell 54154 alpha^n 1 alpha^{n+1} 0.9913079510859701 plane has at t^n no D 0.0008596590841059015 E_final 0.0008596590841063001 E_raw 0.0008596590841063001 alpha^n - alpha_plane(t^n) 0 faces with errors: 159786(-0.0008596590841049874)
 cell 54039 alpha^n 1 alpha^{n+1} 0.9913079510859699 plane has at t^n no D 0.0008596590841045693 E_final 0.0008596590841059042 E_raw 0.0008596590841059042 alpha^n - alpha_plane(t^n) 0 faces with errors: 159444(-0.0008596590841049449)
 cell 54384 alpha^n 1 alpha^{n+1} 0.9913079510859687 plane has at t^n no D 0.000859659084103126 E_final 0.0008596590841050429 E_raw 0.0008596590841050429 alpha^n - alpha_plane(t^n) 0 faces with errors: 160447(-0.0008596590841048915)
 cell 54038 alpha^n 1 alpha^{n+1} 0.8466534194699832 plane has at t^n no D -0.0008596590841013496 E_final -0.0008596590841020709 E_raw -0.0008596590841020709 alpha^n - alpha_plane(t^n) 0 faces with errors: 159444(0.0008596590841049487)
 cell 54269 alpha^n 1 alpha^{n+1} 0.9913079510859701 plane has at t^n no D 0.0008596590841010165 E_final 0.0008596590841064477 E_raw 0.0008596590841064477 alpha^n - alpha_plane(t^n) 0 faces with errors: 160128(-0.0008596590841049331)
 cell 54268 alpha^n 1 alpha^{n+1} 0.8466534194699847 plane has at t^n no D -0.0008596590841004614 E_final -0.0008596590840996634 E_raw -0.0008596590840996634 alpha^n - alpha_plane(t^n) 0 faces with errors: 160128(0.0008596590841049393)
 cell 55943 alpha^n 1 alpha^{n+1} 0.9928921973124363 plane has at t^n no D 0.0007029743169713054 E_final 0.0007029743169677135 E_raw 0.0007029743169677135 alpha^n - alpha_plane(t^n) 0 faces with errors: 165068(-0.000702974316967562)
 cell 55828 alpha^n 1 alpha^{n+1} 0.9928921973124371 plane has at t^n no D 0.0007029743169700842 E_final 0.0007029743169674357 E_raw 0.0007029743169674357 alpha^n - alpha_plane(t^n) 0 faces with errors: 164726(-0.0007029743169674862)
```

### B: the per-step counters over the whole run, row x arm x Co x N

Per row: step-1 faces above 1e-12 V (and their classes none/up/down/both), max over faces of abs(dVf - exact)/V over the run, the largest per-step sum abs(err)/V, max abs(alpha_k - alpha_plane) over the run and at the last step, the largest per-step sum abs(D) V (m^3; the volumeFractionError E_shape), the largest explicit-plane check over the run where it exists (davof rows), and for the comparators the reconstruction's normal error at the end (LNormalInf_ / LNormal1_ of the function object).

| row | arm | label abs(U) dt/h (OpenFOAM Co: x1.6124 on the normal arm, so 0.2 = 0.322 and 0.5 = 0.806, outside the range) | N | steps | step 1: faces (none/up/down/both; boundary) | max abs(alpha_k dev) step 1 | max abs(err)/V run | max sum abs(err)/V per step | max abs(alpha_k dev) run | last step | max sum abs(D) V per step (m^3) | max check | normal err at the end (Inf / 1) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| former | normal | 0.2 | 20 | 20 | 0 (0/0/0/0; 0) | 5.66e-15 | 9.66e-15 | 4.57e-13 | 2.04e-14 | 1.69e-14 | 6.05e-24 | 2.19e-12 |  |
| former | normal | 0.2 | 40 | 40 | 0 (0/0/0/0; 0) | 1.22e-14 | 1.41e-13 | 4.32e-12 | 2.09e-13 | 2.04e-13 | 5.37e-24 | 6.23e-12 |  |
| former | normal | 0.2 | 80 | 80 | 0 (0/0/0/0; 0) | 3.05e-14 | 2.48e-10 | 1.23e-08 | 3.61e-10 | 3.46e-10 | 1.14e-21 | 1.19e-08 |  |
| former | normal | 0.5 | 20 | 8 | 0 (0/0/0/0; 0) | 9.2e-15 | 5.35e-09 | 1.4e-08 | 8.62e-09 | 8.62e-09 | 1.39e-19 | 3.48e-07 |  |
| former | normal | 0.5 | 40 | 16 | 0 (0/0/0/0; 0) | 1.91e-14 | 0.321 | 11.5 | 0.861 | 0.578 | 9.45e-12 | 1.71 |  |
| former | normal | 0.5 | 80 | 32 | 0 (0/0/0/0; 0) | 1.77e-13 | 0.406 | 452 | 1 | 1 | 5.15e-11 | 1.73 |  |
| former | oblique | 0.2 | 20 | 20 | 0 (0/0/0/0; 0) | 4.75e-15 | 1.56e-14 | 2.7e-13 | 4.17e-14 | 1.7e-14 | 6.67e-24 | 1.65e-12 |  |
| former | oblique | 0.2 | 40 | 40 | 0 (0/0/0/0; 0) | 1.24e-14 | 1.18e-12 | 6.93e-12 | 3.29e-12 | 2.53e-12 | 2.52e-23 | 5.74e-11 |  |
| former | oblique | 0.2 | 80 | 80 | 0 (0/0/0/0; 0) | 2.9e-14 | 3.23e-10 | 9.11e-09 | 1.39e-09 | 1.39e-09 | 5.02e-21 | 3.97e-08 |  |
| former | oblique | 0.5 | 20 | 8 | 0 (0/0/0/0; 0) | 8.52e-15 | 2.95e-14 | 6.04e-13 | 4.38e-14 | 4.38e-14 | 7.74e-24 | 2.13e-12 |  |
| former | oblique | 0.5 | 40 | 16 | 0 (0/0/0/0; 0) | 1.85e-14 | 3.26e-12 | 1.95e-11 | 3.3e-12 | 3.3e-12 | 3.66e-23 | 1.61e-10 |  |
| former | oblique | 0.5 | 80 | 32 | 0 (0/0/0/0; 0) | 3.53e-14 | 2.54e-09 | 1.69e-08 | 4.98e-09 | 4.98e-09 | 4.84e-21 | 1.15e-07 |  |
| pure | normal | 0.2 | 20 | 20 | 3 (3/0/0/0; 0) | 2.46e-09 | 0.000122 | 0.00317 | 0.000125 | 6.11e-05 | 1.93e-14 | 0.149 |  |
| pure | normal | 0.2 | 40 | 40 | 0 (0/0/0/0; 0) | 1.28e-14 | 0.000132 | 0.00891 | 0.000119 | 7.71e-05 | 8.7e-15 | 0.135 |  |
| pure | normal | 0.2 | 80 | 80 | 0 (0/0/0/0; 0) | 3.21e-14 | 0.000133 | 0.0457 | 0.000127 | 9.28e-05 | 6.19e-15 | 0.141 |  |
| pure | normal | 0.5 | 20 | 8 | 21 (19/0/0/0; 2) | 0.00101 | 0.00351 | 0.147 | 0.0065 | 0.0036 | 7.58e-13 | 1 |  |
| pure | normal | 0.5 | 40 | 16 | 89 (84/0/0/0; 5) | 0.00086 | 0.00382 | 0.58 | 0.00578 | 0.00402 | 4.99e-13 | 1 |  |
| pure | normal | 0.5 | 80 | 32 | 386 (376/0/0/0; 10) | 0.00092 | 0.00418 | 2.79 | 0.00627 | 0.00439 | 3.14e-13 | 1.23 |  |
| pure | oblique | 0.2 | 20 | 20 | 1 (1/0/0/0; 0) | 1.53e-09 | 0.0106 | 0.17 | 0.0472 | 0.0472 | 6.22e-12 | 1 |  |
| pure | oblique | 0.2 | 40 | 40 | 0 (0/0/0/0; 0) | 1.35e-14 | 0.0135 | 0.506 | 0.0578 | 0.0566 | 2.32e-12 | 1 |  |
| pure | oblique | 0.2 | 80 | 80 | 0 (0/0/0/0; 0) | 2.69e-14 | 0.0182 | 1.02 | 0.0793 | 0.0445 | 5.92e-13 | 1 |  |
| pure | oblique | 0.5 | 20 | 8 | 1 (1/0/0/0; 0) | 3.82e-09 | 0.0169 | 0.299 | 0.028 | 0.028 | 4.46e-12 | 1 |  |
| pure | oblique | 0.5 | 40 | 16 | 0 (0/0/0/0; 0) | 1.71e-14 | 0.0256 | 1.02 | 0.0427 | 0.0422 | 1.9e-12 | 1.01 |  |
| pure | oblique | 0.5 | 80 | 32 | 0 (0/0/0/0; 0) | 3.96e-14 | 0.0344 | 2.23 | 0.0573 | 0.0388 | 5.2e-13 | 1.01 |  |
| onering | normal | 0.2 | 20 | 20 | 3 (3/0/0/0; 0) | 2.46e-09 | 0.000122 | 0.00317 | 0.000125 | 6.11e-05 | 1.93e-14 | 0.149 |  |
| onering | normal | 0.2 | 40 | 40 | 0 (0/0/0/0; 0) | 1.28e-14 | 0.000132 | 0.00891 | 0.000119 | 7.71e-05 | 8.7e-15 | 0.135 |  |
| onering | normal | 0.2 | 80 | 80 | 0 (0/0/0/0; 0) | 3.21e-14 | 0.000133 | 0.0457 | 0.000127 | 9.28e-05 | 6.19e-15 | 0.141 |  |
| onering | normal | 0.5 | 20 | 8 | 21 (19/0/0/0; 2) | 0.00101 | 0.00351 | 0.147 | 0.0065 | 0.0036 | 7.58e-13 | 1 |  |
| onering | normal | 0.5 | 40 | 16 | 89 (84/0/0/0; 5) | 0.00086 | 0.00382 | 0.58 | 0.00578 | 0.00402 | 4.99e-13 | 0.625 |  |
| onering | normal | 0.5 | 80 | 32 | 386 (376/0/0/0; 10) | 0.00092 | 0.00418 | 2.79 | 0.00627 | 0.00439 | 3.14e-13 | 1.23 |  |
| onering | oblique | 0.2 | 20 | 20 | 1 (1/0/0/0; 0) | 1.53e-09 | 4.25e-06 | 4.61e-05 | 3.9e-06 | 3.02e-06 | 7.85e-16 | 0.0151 |  |
| onering | oblique | 0.2 | 40 | 40 | 0 (0/0/0/0; 0) | 1.35e-14 | 4.55e-06 | 0.000214 | 5.35e-06 | 4.19e-06 | 6.84e-16 | 0.015 |  |
| onering | oblique | 0.2 | 80 | 80 | 0 (0/0/0/0; 0) | 2.69e-14 | 4.69e-06 | 0.000941 | 5.69e-06 | 5.56e-06 | 4e-16 | 0.015 |  |
| onering | oblique | 0.5 | 20 | 8 | 1 (1/0/0/0; 0) | 3.82e-09 | 9.56e-05 | 0.00177 | 9.46e-05 | 9.46e-05 | 2.22e-14 | 0.0869 |  |
| onering | oblique | 0.5 | 40 | 16 | 0 (0/0/0/0; 0) | 1.71e-14 | 0.000103 | 0.00979 | 0.000122 | 0.000122 | 1.66e-14 | 0.0761 |  |
| onering | oblique | 0.5 | 80 | 32 | 0 (0/0/0/0; 0) | 3.96e-14 | 0.000139 | 0.0576 | 0.000169 | 0.000158 | 1.25e-14 | 0.0807 |  |
| entering | normal | 0.2 | 20 | 20 | 3 (3/0/0/0; 0) | 2.46e-09 | 9.69e-08 | 8.64e-07 | 1.07e-07 | 5.47e-08 | 9.63e-18 | 1.32e-05 |  |
| entering | normal | 0.2 | 40 | 40 | 0 (0/0/0/0; 0) | 1.28e-14 | 1.1e-14 | 2.16e-12 | 2.72e-14 | 2.72e-14 | 3.68e-24 | 7.27e-12 |  |
| entering | normal | 0.2 | 80 | 80 | 0 (0/0/0/0; 0) | 3.21e-14 | 2.16e-14 | 1.2e-11 | 6.5e-14 | 6.42e-14 | 3.43e-24 | 1.71e-11 |  |
| entering | normal | 0.5 | 20 | 8 | 21 (19/0/0/0; 2) | 0.00101 | 0.00101 | 0.0195 | 0.00134 | 0.000521 | 1.56e-13 | 0.0287 |  |
| entering | normal | 0.5 | 40 | 16 | 89 (84/0/0/0; 5) | 0.00086 | 0.00086 | 0.0759 | 0.00122 | 0.000571 | 8.03e-14 | 1 |  |
| entering | normal | 0.5 | 80 | 32 | 386 (376/0/0/0; 10) | 0.00092 | 0.00092 | 0.438 | 0.00141 | 0.00123 | 5.96e-14 | 1 |  |
| entering | oblique | 0.2 | 20 | 20 | 1 (1/0/0/0; 0) | 1.53e-09 | 5.16e-08 | 2.93e-07 | 4.7e-08 | 1.29e-08 | 2.92e-18 | 1.29e-05 |  |
| entering | oblique | 0.2 | 40 | 40 | 0 (0/0/0/0; 0) | 1.35e-14 | 4.14e-08 | 1.44e-06 | 7.76e-08 | 5.65e-08 | 4.34e-18 | 9.84e-06 |  |
| entering | oblique | 0.2 | 80 | 80 | 0 (0/0/0/0; 0) | 2.69e-14 | 5.77e-08 | 4.76e-06 | 7.76e-08 | 7.76e-08 | 1.96e-18 | 1.19e-05 |  |
| entering | oblique | 0.5 | 20 | 8 | 1 (1/0/0/0; 0) | 3.82e-09 | 8.05e-08 | 4.56e-07 | 6.8e-08 | 1.26e-08 | 4.62e-18 | 1.29e-05 |  |
| entering | oblique | 0.5 | 40 | 16 | 0 (0/0/0/0; 0) | 1.71e-14 | 7.61e-14 | 3.44e-12 | 8e-14 | 4.92e-14 | 4.74e-24 | 3.24e-11 |  |
| entering | oblique | 0.5 | 80 | 32 | 0 (0/0/0/0; 0) | 3.96e-14 | 2.8e-13 | 2.32e-11 | 3.17e-13 | 8.42e-14 | 5.06e-24 | 2.08e-10 |  |
| plicRDF | normal | 0.2 | 20 | 20 | 530 (3/57/0/384; 86) | 0.00281 | 0.00564 | 0.0579 | 0.00281 | 1.95e-05 | 2.31e-13 | n/a | 3.47e-10 / 7.6e-11 |
| plicRDF | normal | 0.2 | 40 | 40 | 2043 (0/267/0/1602; 174) | 0.00281 | 0.00634 | 0.106 | 0.00281 | 4.56e-06 | 5.82e-14 | n/a | 3.17e-12 / 1.05e-12 |
| plicRDF | normal | 0.2 | 80 | 80 | 8057 (0/1158/0/6545; 354) | 0.00247 | 0.00504 | 0.197 | 0.00247 | 5.59e-07 | 1.49e-14 | n/a | 5.06e-14 / 5.94e-15 |
| plicRDF | normal | 0.5 | 20 | 8 | 623 (19/114/0/384; 106) | 0.00336 | 0.023 | 0.559 | 0.0382 | 0.0382 | 3.02e-12 | n/a | 0.0314 / 0.00614 |
| plicRDF | normal | 0.5 | 40 | 16 | 2464 (84/561/0/1603; 216) | 0.00429 | 0.0495 | 2.52 | 0.0858 | 0.0816 | 1.96e-12 | n/a | 0.162 / 0.0133 |
| plicRDF | normal | 0.5 | 80 | 32 | 9839 (376/2477/0/6544; 442) | 0.00263 | 0.0724 | 12.4 | 0.114 | 0.0865 | 1.42e-12 | n/a | 0.0945 / 0.00703 |
| plicRDF | oblique | 0.2 | 20 | 20 | 177 (1/8/0/150; 18) | 0.00193 | 0.00228 | 0.0178 | 0.00193 | 0.000378 | 1.47e-13 | n/a | 9.34e-07 / 2.11e-08 |
| plicRDF | oblique | 0.2 | 40 | 40 | 699 (0/25/0/638; 36) | 0.000976 | 0.00144 | 0.0323 | 0.00165 | 0.000499 | 3.56e-14 | n/a | 5.09e-07 / 2.75e-09 |
| plicRDF | oblique | 0.2 | 80 | 80 | 2773 (0/116/0/2584; 73) | 0.00178 | 0.00205 | 0.0639 | 0.00184 | 0.00021 | 8.75e-15 | n/a | 8.92e-08 / 1.18e-10 |
| plicRDF | oblique | 0.5 | 20 | 8 | 171 (1/17/0/136; 17) | 0.00276 | 0.00316 | 0.0237 | 0.00276 | 0.000665 | 1.95e-13 | n/a | 2.65e-06 / 7e-08 |
| plicRDF | oblique | 0.5 | 40 | 16 | 688 (0/63/0/589; 36) | 0.00157 | 0.00218 | 0.0437 | 0.00157 | 0.000717 | 4.7e-14 | n/a | 1.74e-06 / 9.55e-09 |
| plicRDF | oblique | 0.5 | 80 | 32 | 2750 (0/282/0/2396; 72) | 0.0024 | 0.00265 | 0.0866 | 0.0024 | 0.000295 | 1.19e-14 | n/a | 1.28e-07 / 3.22e-10 |
| isoAlpha | normal | 0.2 | 20 | 20 | 543 (3/57/0/384; 99) | 0.0428 | 0.0655 | 2.02 | 0.229 | 0.222 | 3.16e-11 | n/a | 0.19 / 0.0755 |
| isoAlpha | normal | 0.2 | 40 | 40 | 2072 (0/267/0/1606; 199) | 0.0424 | 0.0718 | 4.68 | 0.233 | 0.232 | 1.03e-11 | n/a | 0.213 / 0.0443 |
| isoAlpha | normal | 0.2 | 80 | 80 | 8157 (0/1158/0/6593; 406) | 0.0459 | 0.0717 | 11.3 | 0.232 | 0.231 | 3.06e-12 | n/a | 0.213 / 0.0231 |
| isoAlpha | normal | 0.5 | 20 | 8 | 629 (19/115/0/384; 111) | 0.0917 | 0.13 | 4.72 | 0.19 | 0.187 | 3.08e-11 | n/a | 0.208 / 0.126 |
| isoAlpha | normal | 0.5 | 40 | 16 | 2476 (84/564/0/1606; 222) | 0.0922 | 0.137 | 10.9 | 0.193 | 0.192 | 1.02e-11 | n/a | 0.193 / 0.0981 |
| isoAlpha | normal | 0.5 | 80 | 32 | 9920 (376/2495/0/6593; 456) | 0.0904 | 0.145 | 26.9 | 0.194 | 0.191 | 3.16e-12 | n/a | 0.166 / 0.0683 |
| isoAlpha | oblique | 0.2 | 20 | 20 | 178 (1/8/0/150; 19) | 0.044 | 0.0632 | 0.713 | 0.235 | 0.201 | 2.58e-11 | n/a | 0.218 / 0.0323 |
| isoAlpha | oblique | 0.2 | 40 | 40 | 702 (0/25/0/638; 39) | 0.0412 | 0.0636 | 1.62 | 0.237 | 0.172 | 7.61e-12 | n/a | 0.182 / 0.0167 |
| isoAlpha | oblique | 0.2 | 80 | 80 | 2776 (0/116/0/2584; 76) | 0.0359 | 0.0628 | 3.9 | 0.238 | 0.234 | 2.26e-12 | n/a | 0.201 / 0.00861 |
| isoAlpha | oblique | 0.5 | 20 | 8 | 171 (1/17/0/136; 17) | 0.0541 | 0.134 | 1.47 | 0.222 | 0.2 | 2.16e-11 | n/a | 0.219 / 0.0337 |
| isoAlpha | oblique | 0.5 | 40 | 16 | 691 (0/63/0/589; 39) | 0.0483 | 0.12 | 3.55 | 0.23 | 0.173 | 6.81e-12 | n/a | 0.219 / 0.0172 |
| isoAlpha | oblique | 0.5 | 80 | 32 | 2754 (0/282/0/2396; 76) | 0.0646 | 0.133 | 8.79 | 0.224 | 0.223 | 2.06e-12 | n/a | 0.2 / 0.0086 |


## The plane ladder with the Courant numbers labelled (2026-11-23, third task)

OpenFOAM Co = (|u_x| + |u_y| + |u_z|) dt/h as the solver prints it; |U| dt/h is the label of the earlier tables. Normal arm: OpenFOAM Co = 1.6124 |U| dt/h; oblique arm: equal. The rows with |U| dt/h = 0.2 and 0.5 on the normal arm are OpenFOAM Co 0.322 and 0.806 (0.806 outside the range of the rows run from now on). L1 rel = sum |alpha_k - alpha_plane| V_k / V_domain (V_domain = 0.004^3) at the last step and, in the next column, its maximum over the steps; order = log2 of the ratio of the last-step L1 rel between N and 2N.

| row | arm | OpenFOAM Co | abs(U) dt/h | N | steps | step 1: faces above 1e-12 V | alpha_k dev max step 1 | max abs(err)/V over run | alpha_k dev max over run | L1 rel last | L1 rel max over steps | order (N to 2N) | max check (davof rows) | comparators: abs(n - n_e) L1 / Linf (last step) | 1 - n.n_e: L1 (time avg) / Linf (max over run) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| former | normal | 0.322 | 0.2 | 20 | 20 | 0 | 5.66e-15 | 9.66e-15 | 2.04e-14 | 4.91e-17 | 9.45e-17 | n/a | 2.19e-12 |  |  |
| former | normal | 0.322 | 0.2 | 40 | 40 | 0 | 1.22e-14 | 1.41e-13 | 2.09e-13 | 7.45e-17 | 8.39e-17 | n/a | 6.23e-12 |  |  |
| former | normal | 0.322 | 0.2 | 80 | 80 | 0 | 3.05e-14 | 2.48e-10 | 3.61e-10 | 1.79e-14 | 1.79e-14 |  | 1.19e-08 |  |  |
| former | normal | 0.2 | 0.12404 | 20 | 32 | 0 | 5.55e-15 | 9.59e-15 | 1.44e-14 | 2.8e-17 | 7.97e-17 | n/a | 2.01e-12 |  |  |
| former | normal | 0.2 | 0.12404 | 40 | 64 | 0 | 1.32e-14 | 1.05e-14 | 1.35e-10 | 2.39e-17 | 4.28e-15 | n/a | 3.68e-12 |  |  |
| former | normal | 0.2 | 0.12404 | 80 | 129 | 0 | 2.45e-14 | 1.44e-14 | 1.22e-09 | 2.15e-17 | 2.3e-14 |  | 4e-12 |  |  |
| former | normal | 0.806 (outside the range) | 0.5 | 20 | 8 | 0 | 9.2e-15 | 5.35e-09 | 8.62e-09 | 2.17e-12 | 2.17e-12 | -25.99 | 3.48e-07 |  |  |
| former | normal | 0.806 (outside the range) | 0.5 | 40 | 16 | 0 | 1.91e-14 | 0.321 | 0.861 | 0.000145 | 0.000148 | -1.67 | 1.71 |  |  |
| former | normal | 0.806 (outside the range) | 0.5 | 80 | 32 | 0 | 1.77e-13 | 0.406 | 1 | 0.000462 | 0.000805 |  | 1.73 |  |  |
| former | normal | 0.5 | 0.3101 | 20 | 13 | 0 | 5.43e-15 | 7.61e-09 | 8.28e-09 | 2.48e-12 | 3.07e-12 | -25.65 | 2.95e-07 |  |  |
| former | normal | 0.5 | 0.3101 | 40 | 26 | 0 | 1.8e-14 | 0.214 | 0.955 | 0.000131 | 0.000133 | -1.73 | 1.72 |  |  |
| former | normal | 0.5 | 0.3101 | 80 | 52 | 0 | 3.26e-14 | 0.252 | 1 | 0.000434 | 0.000761 |  | 1.73 |  |  |
| former | oblique | 0.2 | 0.2 | 20 | 20 | 0 | 4.75e-15 | 1.56e-14 | 4.17e-14 | 1.02e-16 | 1.04e-16 | n/a | 1.65e-12 |  |  |
| former | oblique | 0.2 | 0.2 | 40 | 40 | 0 | 1.24e-14 | 1.18e-12 | 3.29e-12 | 3.93e-16 | 3.93e-16 | n/a | 5.74e-11 |  |  |
| former | oblique | 0.2 | 0.2 | 80 | 80 | 0 | 2.9e-14 | 3.23e-10 | 1.39e-09 | 7.84e-14 | 7.84e-14 |  | 3.97e-08 |  |  |
| former | oblique | 0.5 | 0.5 | 20 | 8 | 0 | 8.52e-15 | 2.95e-14 | 4.38e-14 | 1.21e-16 | 1.21e-16 | n/a | 2.13e-12 |  |  |
| former | oblique | 0.5 | 0.5 | 40 | 16 | 0 | 1.85e-14 | 3.26e-12 | 3.3e-12 | 5.71e-16 | 5.71e-16 | n/a | 1.61e-10 |  |  |
| former | oblique | 0.5 | 0.5 | 80 | 32 | 0 | 3.53e-14 | 2.54e-09 | 4.98e-09 | 7.57e-14 | 7.57e-14 |  | 1.15e-07 |  |  |
| pure | normal | 0.322 | 0.2 | 20 | 20 | 3 | 2.46e-09 | 0.000122 | 0.000125 | 4.6e-08 | 3.01e-07 | 1.60 | 0.149 |  |  |
| pure | normal | 0.322 | 0.2 | 40 | 40 | 0 | 1.28e-14 | 0.000132 | 0.000119 | 1.52e-08 | 1.36e-07 | 1.02 | 0.135 |  |  |
| pure | normal | 0.322 | 0.2 | 80 | 80 | 0 | 3.21e-14 | 0.000133 | 0.000127 | 7.53e-09 | 9.67e-08 |  | 0.141 |  |  |
| pure | normal | 0.2 | 0.12404 | 20 | 32 | 3 | 1.53e-09 | 1.12e-05 | 2.69e-05 | 1.36e-08 | 4.89e-08 | 0.95 | 0.0607 |  |  |
| pure | normal | 0.2 | 0.12404 | 40 | 64 | 0 | 1.2e-14 | 1.78e-05 | 2.81e-05 | 7.02e-09 | 3.23e-08 | 1.35 | 0.0602 |  |  |
| pure | normal | 0.2 | 0.12404 | 80 | 129 | 0 | 2.42e-14 | 2e-05 | 3.07e-05 | 2.75e-09 | 2.46e-08 |  | 0.0607 |  |  |
| pure | normal | 0.806 (outside the range) | 0.5 | 20 | 8 | 21 | 0.00101 | 0.00351 | 0.0065 | 2.19e-06 | 1.18e-05 | 1.25 | 1 |  |  |
| pure | normal | 0.806 (outside the range) | 0.5 | 40 | 16 | 89 | 0.00086 | 0.00382 | 0.00578 | 9.25e-07 | 7.8e-06 | 1.16 | 1 |  |  |
| pure | normal | 0.806 (outside the range) | 0.5 | 80 | 32 | 386 | 0.00092 | 0.00418 | 0.00627 | 4.14e-07 | 4.91e-06 |  | 1.23 |  |  |
| pure | normal | 0.5 | 0.3101 | 20 | 13 | 3 | 2.05e-07 | 0.000437 | 0.000724 | 3.37e-07 | 2.02e-06 | 1.28 | 0.333 |  |  |
| pure | normal | 0.5 | 0.3101 | 40 | 26 | 4 | 2.37e-10 | 0.000541 | 0.00091 | 1.39e-07 | 1.14e-06 | 1.42 | 0.335 |  |  |
| pure | normal | 0.5 | 0.3101 | 80 | 52 | 20 | 2.15e-08 | 0.000589 | 0.000914 | 5.2e-08 | 6.97e-07 |  | 1.14 |  |  |
| pure | oblique | 0.2 | 0.2 | 20 | 20 | 1 | 1.53e-09 | 0.0106 | 0.0472 | 9.44e-05 | 9.72e-05 | 1.47 | 1 |  |  |
| pure | oblique | 0.2 | 0.2 | 40 | 40 | 0 | 1.35e-14 | 0.0135 | 0.0578 | 3.42e-05 | 3.62e-05 | 2.24 | 1 |  |  |
| pure | oblique | 0.2 | 0.2 | 80 | 80 | 0 | 2.69e-14 | 0.0182 | 0.0793 | 7.23e-06 | 9.25e-06 |  | 1 |  |  |
| pure | oblique | 0.5 | 0.5 | 20 | 8 | 1 | 3.82e-09 | 0.0169 | 0.028 | 6.97e-05 | 6.97e-05 | 1.26 | 1 |  |  |
| pure | oblique | 0.5 | 0.5 | 40 | 16 | 0 | 1.71e-14 | 0.0256 | 0.0427 | 2.91e-05 | 2.96e-05 | 2.08 | 1.01 |  |  |
| pure | oblique | 0.5 | 0.5 | 80 | 32 | 0 | 3.96e-14 | 0.0344 | 0.0573 | 6.91e-06 | 8.12e-06 |  | 1.01 |  |  |
| onering | normal | 0.322 | 0.2 | 20 | 20 | 3 | 2.46e-09 | 0.000122 | 0.000125 | 4.6e-08 | 3.01e-07 | 1.60 | 0.149 |  |  |
| onering | normal | 0.322 | 0.2 | 40 | 40 | 0 | 1.28e-14 | 0.000132 | 0.000119 | 1.52e-08 | 1.36e-07 | 1.02 | 0.135 |  |  |
| onering | normal | 0.322 | 0.2 | 80 | 80 | 0 | 3.21e-14 | 0.000133 | 0.000127 | 7.53e-09 | 9.67e-08 |  | 0.141 |  |  |
| onering | normal | 0.2 | 0.12404 | 20 | 32 | 3 | 1.53e-09 | 1.12e-05 | 2.69e-05 | 1.36e-08 | 4.89e-08 | 0.95 | 0.0607 |  |  |
| onering | normal | 0.2 | 0.12404 | 40 | 64 | 0 | 1.2e-14 | 1.78e-05 | 2.81e-05 | 7.02e-09 | 3.23e-08 | 1.35 | 0.0602 |  |  |
| onering | normal | 0.2 | 0.12404 | 80 | 129 | 0 | 2.42e-14 | 2e-05 | 3.07e-05 | 2.75e-09 | 2.46e-08 |  | 0.0607 |  |  |
| onering | normal | 0.806 (outside the range) | 0.5 | 20 | 8 | 21 | 0.00101 | 0.00351 | 0.0065 | 2.19e-06 | 1.18e-05 | 1.25 | 1 |  |  |
| onering | normal | 0.806 (outside the range) | 0.5 | 40 | 16 | 89 | 0.00086 | 0.00382 | 0.00578 | 9.25e-07 | 7.8e-06 | 1.16 | 0.625 |  |  |
| onering | normal | 0.806 (outside the range) | 0.5 | 80 | 32 | 386 | 0.00092 | 0.00418 | 0.00627 | 4.14e-07 | 4.91e-06 |  | 1.23 |  |  |
| onering | normal | 0.5 | 0.3101 | 20 | 13 | 3 | 2.05e-07 | 0.000437 | 0.000724 | 3.37e-07 | 2.02e-06 | 1.28 | 0.333 |  |  |
| onering | normal | 0.5 | 0.3101 | 40 | 26 | 4 | 2.37e-10 | 0.000541 | 0.00091 | 1.39e-07 | 1.14e-06 | 1.42 | 0.335 |  |  |
| onering | normal | 0.5 | 0.3101 | 80 | 52 | 20 | 2.15e-08 | 0.000589 | 0.000914 | 5.2e-08 | 6.97e-07 |  | 1.14 |  |  |
| onering | oblique | 0.2 | 0.2 | 20 | 20 | 1 | 1.53e-09 | 4.25e-06 | 3.9e-06 | 1.09e-08 | 1.23e-08 | 0.19 | 0.0151 |  |  |
| onering | oblique | 0.2 | 0.2 | 40 | 40 | 0 | 1.35e-14 | 4.55e-06 | 5.35e-06 | 9.58e-09 | 1.07e-08 | 0.93 | 0.015 |  |  |
| onering | oblique | 0.2 | 0.2 | 80 | 80 | 0 | 2.69e-14 | 4.69e-06 | 5.69e-06 | 5.04e-09 | 6.24e-09 |  | 0.015 |  |  |
| onering | oblique | 0.5 | 0.5 | 20 | 8 | 1 | 3.82e-09 | 9.56e-05 | 9.46e-05 | 3.36e-07 | 3.47e-07 | 0.47 | 0.0869 |  |  |
| onering | oblique | 0.5 | 0.5 | 40 | 16 | 0 | 1.71e-14 | 0.000103 | 0.000122 | 2.42e-07 | 2.6e-07 | 0.49 | 0.0761 |  |  |
| onering | oblique | 0.5 | 0.5 | 80 | 32 | 0 | 3.96e-14 | 0.000139 | 0.000169 | 1.73e-07 | 1.96e-07 |  | 0.0807 |  |  |
| entering | normal | 0.322 | 0.2 | 20 | 20 | 3 | 2.46e-09 | 9.69e-08 | 1.07e-07 | 2.25e-11 | 1.51e-10 | n/a | 1.32e-05 |  |  |
| entering | normal | 0.322 | 0.2 | 40 | 40 | 0 | 1.28e-14 | 1.1e-14 | 2.72e-14 | 8.63e-18 | 5.75e-17 | n/a | 7.27e-12 |  |  |
| entering | normal | 0.322 | 0.2 | 80 | 80 | 0 | 3.21e-14 | 2.16e-14 | 6.5e-14 | 1.14e-17 | 5.36e-17 |  | 1.71e-11 |  |  |
| entering | normal | 0.2 | 0.12404 | 20 | 32 | 3 | 1.53e-09 | 6.39e-08 | 7.12e-08 | 2.53e-11 | 8.46e-11 | 2.62 | 1.46e-05 |  |  |
| entering | normal | 0.2 | 0.12404 | 40 | 64 | 0 | 1.2e-14 | 1.01e-07 | 5.02e-08 | 4.13e-12 | 4.56e-11 | -0.34 | 1.36e-05 |  |  |
| entering | normal | 0.2 | 0.12404 | 80 | 129 | 0 | 2.42e-14 | 9.68e-08 | 7.83e-08 | 5.22e-12 | 2.59e-11 |  | 1.54e-05 |  |  |
| entering | normal | 0.806 (outside the range) | 0.5 | 20 | 8 | 21 | 0.00101 | 0.00101 | 0.00134 | 3.46e-07 | 2.44e-06 | 1.12 | 0.0287 |  |  |
| entering | normal | 0.806 (outside the range) | 0.5 | 40 | 16 | 89 | 0.00086 | 0.00086 | 0.00122 | 1.59e-07 | 1.25e-06 | 0.78 | 1 |  |  |
| entering | normal | 0.806 (outside the range) | 0.5 | 80 | 32 | 386 | 0.00092 | 0.00092 | 0.00141 | 9.23e-08 | 9.31e-07 |  | 1 |  |  |
| entering | normal | 0.5 | 0.3101 | 20 | 13 | 3 | 2.05e-07 | 2.05e-07 | 2.61e-07 | 8.24e-11 | 3.89e-10 | 1.02 | 9.9e-05 |  |  |
| entering | normal | 0.5 | 0.3101 | 40 | 26 | 4 | 2.37e-10 | 1.79e-07 | 2.69e-07 | 4.05e-11 | 1.8e-10 | 1.44 | 9.06e-05 |  |  |
| entering | normal | 0.5 | 0.3101 | 80 | 52 | 20 | 2.15e-08 | 2.13e-07 | 2.75e-07 | 1.49e-11 | 1.52e-10 |  | 0.000102 |  |  |
| entering | oblique | 0.2 | 0.2 | 20 | 20 | 1 | 1.53e-09 | 5.16e-08 | 4.7e-08 | 1.03e-11 | 4.56e-11 | -2.34 | 1.29e-05 |  |  |
| entering | oblique | 0.2 | 0.2 | 40 | 40 | 0 | 1.35e-14 | 4.14e-08 | 7.76e-08 | 5.21e-11 | 6.77e-11 | 0.99 | 9.84e-06 |  |  |
| entering | oblique | 0.2 | 0.2 | 80 | 80 | 0 | 2.69e-14 | 5.77e-08 | 7.76e-08 | 2.62e-11 | 3.06e-11 |  | 1.19e-05 |  |  |
| entering | oblique | 0.5 | 0.5 | 20 | 8 | 1 | 3.82e-09 | 8.05e-08 | 6.8e-08 | 1.04e-11 | 7.22e-11 | n/a | 1.29e-05 |  |  |
| entering | oblique | 0.5 | 0.5 | 40 | 16 | 0 | 1.71e-14 | 7.61e-14 | 8e-14 | 6.78e-17 | 7.41e-17 | n/a | 3.24e-11 |  |  |
| entering | oblique | 0.5 | 0.5 | 80 | 32 | 0 | 3.96e-14 | 2.8e-13 | 3.17e-13 | 6.87e-17 | 7.91e-17 |  | 2.08e-10 |  |  |
| plicRDF | normal | 0.322 | 0.2 | 20 | 20 | 530 | 0.00281 | 0.00564 | 0.00281 | 1.14e-08 | 3.61e-06 | 3.59 | n/a | 8e-12 / 1.23e-05 | 7.6e-11 / 3.47e-10 |
| plicRDF | normal | 0.322 | 0.2 | 40 | 40 | 2043 | 0.00281 | 0.00634 | 0.00281 | 9.47e-10 | 9.09e-07 | 5.10 | n/a | 1e-12 / 1.45e-06 | 1.05e-12 / 3.17e-12 |
| plicRDF | normal | 0.322 | 0.2 | 80 | 80 | 8057 | 0.00247 | 0.00504 | 0.00247 | 2.76e-11 | 2.33e-07 |  | n/a | 1.25e-13 / 1.1e-07 | 5.94e-15 / 5.06e-14 |
| plicRDF | normal | 0.2 | 0.12404 | 20 | 32 | 511 | 0.00233 | 0.00415 | 0.00233 | 8.11e-09 | 2.74e-06 | 3.62 | n/a | 8e-12 / 8.29e-06 | 3.43e-11 / 1.55e-10 |
| plicRDF | normal | 0.2 | 0.12404 | 40 | 64 | 1921 | 0.00253 | 0.00485 | 0.00253 | 6.57e-10 | 6.64e-07 | 5.17 | n/a | 1e-12 / 9.98e-07 | 4.98e-13 / 1.58e-12 |
| plicRDF | normal | 0.2 | 0.12404 | 80 | 129 | 7580 | 0.00166 | 0.00338 | 0.00185 | 1.83e-11 | 1.67e-07 |  | n/a | 1.25e-13 / 7.23e-08 | 2.55e-15 / 2.11e-14 |
| plicRDF | normal | 0.806 (outside the range) | 0.5 | 20 | 8 | 623 | 0.00336 | 0.023 | 0.0382 | 2.75e-05 | 4.71e-05 | 0.36 | n/a | 8e-12 / 0.111 | 0.00614 / 0.0314 |
| plicRDF | normal | 0.806 (outside the range) | 0.5 | 40 | 16 | 2464 | 0.00429 | 0.0495 | 0.0858 | 2.14e-05 | 3.07e-05 | 1.55 | n/a | 1e-12 / 0.163 | 0.0133 / 0.162 |
| plicRDF | normal | 0.806 (outside the range) | 0.5 | 80 | 32 | 9839 | 0.00263 | 0.0724 | 0.114 | 7.3e-06 | 2.22e-05 |  | n/a | 1.25e-13 / 0.119 | 0.00703 / 0.0945 |
| plicRDF | normal | 0.5 | 0.3101 | 20 | 13 | 563 | 0.00226 | 0.00636 | 0.00226 | 1.58e-08 | 4.53e-06 | 3.84 | n/a | 8e-12 / 1.51e-05 | 1.13e-10 / 2.71e-10 |
| plicRDF | normal | 0.5 | 0.3101 | 40 | 26 | 2195 | 0.00287 | 0.0068 | 0.00287 | 1.1e-09 | 1.16e-06 | 2.21 | n/a | 1e-12 / 3.15e-06 | 4.97e-12 / 5.04e-11 |
| plicRDF | normal | 0.5 | 0.3101 | 80 | 52 | 8717 | 0.00356 | 0.00729 | 0.00356 | 2.38e-10 | 2.96e-07 |  | n/a | 1.25e-13 / 5.13e-06 | 1.32e-11 / 7.3e-10 |
| plicRDF | oblique | 0.2 | 0.2 | 20 | 20 | 177 | 0.00193 | 0.00228 | 0.00193 | 6.69e-07 | 2.3e-06 | 3.03 | n/a | 8e-12 / 0.000205 | 2.11e-08 / 9.34e-07 |
| plicRDF | oblique | 0.2 | 0.2 | 40 | 40 | 699 | 0.000976 | 0.00144 | 0.00165 | 8.2e-08 | 5.56e-07 | 4.43 | n/a | 1e-12 / 7.42e-05 | 2.75e-09 / 5.09e-07 |
| plicRDF | oblique | 0.2 | 0.2 | 80 | 80 | 2773 | 0.00178 | 0.00205 | 0.00184 | 3.81e-09 | 1.37e-07 |  | n/a | 1.25e-13 / 1.54e-05 | 1.18e-10 / 8.92e-08 |
| plicRDF | oblique | 0.5 | 0.5 | 20 | 8 | 171 | 0.00276 | 0.00316 | 0.00276 | 1.19e-06 | 3.05e-06 | 2.69 | n/a | 8e-12 / 0.000374 | 7e-08 / 2.65e-06 |
| plicRDF | oblique | 0.5 | 0.5 | 40 | 16 | 688 | 0.00157 | 0.00218 | 0.00157 | 1.84e-07 | 7.34e-07 | 3.83 | n/a | 1e-12 / 0.000138 | 9.55e-09 / 1.74e-06 |
| plicRDF | oblique | 0.5 | 0.5 | 80 | 32 | 2750 | 0.0024 | 0.00265 | 0.0024 | 1.3e-08 | 1.86e-07 |  | n/a | 1.25e-13 / 2.54e-05 | 3.22e-10 / 1.28e-07 |
| isoAlpha | normal | 0.322 | 0.2 | 20 | 20 | 543 | 0.0428 | 0.0655 | 0.229 | 0.000211 | 0.000494 | 1.81 | n/a | 8e-12 / 0.389 | 0.0755 / 0.19 |
| isoAlpha | normal | 0.322 | 0.2 | 40 | 40 | 2072 | 0.0424 | 0.0718 | 0.233 | 6.02e-05 | 0.000161 | 1.88 | n/a | 1e-12 / 0.298 | 0.0443 / 0.213 |
| isoAlpha | normal | 0.322 | 0.2 | 80 | 80 | 8157 | 0.0459 | 0.0717 | 0.232 | 1.63e-05 | 4.78e-05 |  | n/a | 1.25e-13 / 0.215 | 0.0231 / 0.213 |
| isoAlpha | normal | 0.2 | 0.12404 | 20 | 32 | 523 | 0.0283 | 0.0432 | 0.235 | 0.000218 | 0.000501 | 1.81 | n/a | 8e-12 / 0.383 | 0.0732 / 0.186 |
| isoAlpha | normal | 0.2 | 0.12404 | 40 | 64 | 1955 | 0.0265 | 0.0456 | 0.238 | 6.23e-05 | 0.000164 | 1.91 | n/a | 1e-12 / 0.295 | 0.0436 / 0.205 |
| isoAlpha | normal | 0.2 | 0.12404 | 80 | 129 | 7694 | 0.0298 | 0.0456 | 0.238 | 1.66e-05 | 4.85e-05 |  | n/a | 1.25e-13 / 0.215 | 0.0231 / 0.213 |
| isoAlpha | normal | 0.806 (outside the range) | 0.5 | 20 | 8 | 629 | 0.0917 | 0.13 | 0.19 | 0.000206 | 0.000481 | 1.84 | n/a | 8e-12 / 0.503 | 0.126 / 0.208 |
| isoAlpha | normal | 0.806 (outside the range) | 0.5 | 40 | 16 | 2476 | 0.0922 | 0.137 | 0.193 | 5.77e-05 | 0.000159 | 1.90 | n/a | 1e-12 / 0.443 | 0.0981 / 0.193 |
| isoAlpha | normal | 0.806 (outside the range) | 0.5 | 80 | 32 | 9920 | 0.0904 | 0.145 | 0.194 | 1.55e-05 | 4.94e-05 |  | n/a | 1.25e-13 / 0.37 | 0.0683 / 0.166 |
| isoAlpha | normal | 0.5 | 0.3101 | 20 | 13 | 573 | 0.0618 | 0.0894 | 0.216 | 0.000204 | 0.000485 | 1.82 | n/a | 8e-12 / 0.41 | 0.0841 / 0.19 |
| isoAlpha | normal | 0.5 | 0.3101 | 40 | 26 | 2221 | 0.0666 | 0.105 | 0.223 | 5.76e-05 | 0.000158 | 1.88 | n/a | 1e-12 / 0.31 | 0.0482 / 0.216 |
| isoAlpha | normal | 0.5 | 0.3101 | 80 | 52 | 8807 | 0.066 | 0.105 | 0.221 | 1.56e-05 | 4.73e-05 |  | n/a | 1.25e-13 / 0.219 | 0.024 / 0.219 |
| isoAlpha | oblique | 0.2 | 0.2 | 20 | 20 | 178 | 0.044 | 0.0632 | 0.235 | 0.000368 | 0.000402 | 1.81 | n/a | 8e-12 / 0.254 | 0.0323 / 0.218 |
| isoAlpha | oblique | 0.2 | 0.2 | 40 | 40 | 702 | 0.0412 | 0.0636 | 0.237 | 0.000105 | 0.000119 | 1.82 | n/a | 1e-12 / 0.183 | 0.0167 / 0.182 |
| isoAlpha | oblique | 0.2 | 0.2 | 80 | 80 | 2776 | 0.0359 | 0.0628 | 0.238 | 2.97e-05 | 3.53e-05 |  | n/a | 1.25e-13 / 0.131 | 0.00861 / 0.201 |
| isoAlpha | oblique | 0.5 | 0.5 | 20 | 8 | 171 | 0.0541 | 0.134 | 0.222 | 0.000329 | 0.000337 | 1.75 | n/a | 8e-12 / 0.26 | 0.0337 / 0.219 |
| isoAlpha | oblique | 0.5 | 0.5 | 40 | 16 | 691 | 0.0483 | 0.12 | 0.23 | 9.78e-05 | 0.000106 | 1.80 | n/a | 1e-12 / 0.185 | 0.0172 / 0.219 |
| isoAlpha | oblique | 0.5 | 0.5 | 80 | 32 | 2754 | 0.0646 | 0.133 | 0.224 | 2.82e-05 | 3.21e-05 |  | n/a | 1.25e-13 / 0.131 | 0.0086 / 0.2 |

### Evaluation of the pre-registered statements (OpenFOAM Co 0.2 and 0.5, normal arm)

PE1 (ENTERING: 0 step-1 faces above 1e-12 V_k and max alpha_k deviation over the run <= 1e-10 at N >= 40; N 20: <= 1e-6):

| OpenFOAM Co | N | step-1 faces | max alpha_k dev over the run | within PE1 |
|---|---|---|---|---|
| 0.2 | 20 | 3 | 7.12e-08 | yes |
| 0.2 | 40 | 0 | 5.02e-08 | NO |
| 0.2 | 80 | 0 | 7.83e-08 | NO |
| 0.5 | 20 | 3 | 2.61e-07 | yes |
| 0.5 | 40 | 4 | 2.69e-07 | NO |
| 0.5 | 80 | 20 | 2.75e-07 | NO |

PE1b (FORMER below 1e-8 at every N):

| OpenFOAM Co | N | max alpha_k dev over the run | below 1e-8 |
|---|---|---|---|
| 0.2 | 20 | 1.44e-14 | yes |
| 0.2 | 40 | 1.35e-10 | yes |
| 0.2 | 80 | 1.22e-09 | yes |
| 0.5 | 20 | 8.28e-09 | yes |
| 0.5 | 40 | 0.955 | NO |
| 0.5 | 80 | 1 | NO |

PE1c (PURE's L1 deviation below plicRDF's on the normal arm, above it on the oblique arm; L1 rel at the last step / maximum over the steps):

| arm | OpenFOAM Co | N | PURE last / max | plicRDF last / max | PURE below plicRDF (last, max) |
|---|---|---|---|---|---|
| normal | 0.2 | 20 | 1.36e-08 / 4.89e-08 | 8.11e-09 / 2.74e-06 | no, yes |
| normal | 0.2 | 40 | 7.02e-09 / 3.23e-08 | 6.57e-10 / 6.64e-07 | no, yes |
| normal | 0.2 | 80 | 2.75e-09 / 2.46e-08 | 1.83e-11 / 1.67e-07 | no, yes |
| normal | 0.5 | 20 | 3.37e-07 / 2.02e-06 | 1.58e-08 / 4.53e-06 | no, yes |
| normal | 0.5 | 40 | 1.39e-07 / 1.14e-06 | 1.1e-09 / 1.16e-06 | no, yes |
| normal | 0.5 | 80 | 5.2e-08 / 6.97e-07 | 2.38e-10 / 2.96e-07 | no, no |
| normal | 0.322 | 20 | 4.6e-08 / 3.01e-07 | 1.14e-08 / 3.61e-06 | no, yes |
| normal | 0.322 | 40 | 1.52e-08 / 1.36e-07 | 9.47e-10 / 9.09e-07 | no, yes |
| normal | 0.322 | 80 | 7.53e-09 / 9.67e-08 | 2.76e-11 / 2.33e-07 | no, yes |
| normal | 0.806 | 20 | 2.19e-06 / 1.18e-05 | 2.75e-05 / 4.71e-05 | yes, yes |
| normal | 0.806 | 40 | 9.25e-07 / 7.8e-06 | 2.14e-05 / 3.07e-05 | yes, yes |
| normal | 0.806 | 80 | 4.14e-07 / 4.91e-06 | 7.3e-06 / 2.22e-05 | yes, yes |
| oblique | 0.2 | 20 | 9.44e-05 / 9.72e-05 | 6.69e-07 / 2.3e-06 | no, no |
| oblique | 0.2 | 40 | 3.42e-05 / 3.62e-05 | 8.2e-08 / 5.56e-07 | no, no |
| oblique | 0.2 | 80 | 7.23e-06 / 9.25e-06 | 3.81e-09 / 1.37e-07 | no, no |
| oblique | 0.5 | 20 | 6.97e-05 / 6.97e-05 | 1.19e-06 / 3.05e-06 | no, no |
| oblique | 0.5 | 40 | 2.91e-05 / 2.96e-05 | 1.84e-07 / 7.34e-07 | no, no |
| oblique | 0.5 | 80 | 6.91e-06 / 8.12e-06 | 1.3e-08 / 1.86e-07 | no, no |

