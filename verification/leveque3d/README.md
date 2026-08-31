# 3D LeVeque deformation: retracing vs non-retracing

The most stringent of the three studies. One command:

```bash
./run.sh
```

## The field and its potential

The fully three-dimensional deformation field on the unit cube,

    u =  2 sin^2(pi x) sin(2 pi y) sin(2 pi z)
    v = -sin(2 pi x) sin^2(pi y) sin(2 pi z)
    w = -sin(2 pi x) sin(2 pi y) sin^2(pi z)

does not separate into horizontal and axial parts the way the spiralling field
of `../deformation-sphere3d` does, so its vector potential is less obvious. One
exists all the same, obtained by integrating dA3/dx = -v and dA2/dx = w and
fixing the remaining gauge freedom so that curl(A)_x comes out right:

    A = ( 0,
          S(y) [ C(x) s(z) + C(z) ] / 2pi,
         -C(x) s(y) S(z) / 2pi )

with S(a) = sin(2 pi a), C(a) = cos(2 pi a), s(a) = sin^2(pi a). Verified
numerically: max |curl(A) - u| = 1e-10 over 2000 random points in the cube,
i.e. finite-difference truncation.

A is NOT truncated outside the cube. Unlike psi0 it does not vanish on every
face, so a zero extension would put a tangential jump -- a spurious surface
current -- into curl(A); evaluated as written it continues periodically, which
stays smooth and solenoidal.

## Against the reference

The reference case builds phi as a face-centre quadrature and repairs it with
CorrectPhi, a Poisson projection. Unframed, plicRDF, cubic cells:

| N  | quadrature + CorrectPhi | vector potential | difference |
|----|-------------------------|------------------|------------|
| 32 | 7.743e-03               | 7.720e-03        | 0.30%      |
| 64 | 3.079e-03               | 3.001e-03        | 2.5%       |

`max|div(u)| <= 3.4e-13` and volume is conserved to 5e-15, with no linear solve.

## Mesh

The committed `simulationParameter` of `test-Leveque` carries `nz 64` against a
unit-cube domain, i.e. cells with `dz = dx/2`. The testsuite overrides it with
`nz = res`; this workflow does the same (`nz_factor: 1`). Comparing against the
published errors on the anisotropic mesh is what makes the difference look like
17% rather than 0.3%.

## Frame

Rotation about the axis parallel to `e_z` through `(0.5, 0.5)`. The cube is
invariant under the base field (the velocity vanishes on all six faces), and the
material was measured to stay within horizontal radius 0.46 of that axis, so it
remains inside the cube under the rotation. The margin is 0.04 -- 1.3 cells at
N=32 -- so watch `E_mass` at the coarsest level, exactly as in the 2D study.

## Data archive

The figures and tables in the manuscript are regenerated from this study's
`results/` directory, which is what the archive cited as `figshare2026`
contains:

```
<study>/results/convergence.csv               errors and observed orders
<study>/results/convergence_table.tex         the LaTeX tables
<study>/results/convergence_shape_error.png   the convergence diagram
<study>/results/interface*_*_frame.png        the interface snapshots
<study>/results/potential_check.txt           curl(A0) against u0 by
                                              complex-step differentiation
```

Every one of those is an output of `./run.sh`; none is produced by hand.
