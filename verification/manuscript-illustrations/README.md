# Manuscript illustrations

Two stand-alone scripts behind parts of the manuscript that use no simulation.

| Script | Manuscript | Output |
|---|---|---|
| `two_motions_storyboard.py` | Figure 1: base deformation, frame motion alone, both together | `two_motions_storyboard.{pdf,svg,png}`; with `--animate` also `two_motions_animation.gif` |
| `cancellation_1d.py` | Section 2: the closed forms (half-period error, zero terminal error, the frame's terminal error $\beta h^2/(4T)$) and the sine-frame remark | printed residuals, all zero |

```bash
pip install -r ../requirements.txt
python3 two_motions_storyboard.py OUTPUT_DIR [--animate]
python3 cancellation_1d.py
```

Figure 1 is computed from closed-form motions: a time-periodic shear
$k(t) = 1.6\sin^2(\pi t/T)$ of a disc of radius 0.26 centred at (0.55, 0),
carried by one full turn of a rigid frame about the origin. No ODE is integrated.
The script asserts the exact return at $t = T$ and the mirror relations
between $T/4$ and $3T/4$ before it draws.
