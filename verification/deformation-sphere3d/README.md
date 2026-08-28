# 3D spiralling deformation: retracing vs non-retracing

The three-dimensional counterpart of `../reversed-vortex`. One command:

```bash
./run.sh
```

## What is different from 2D

The reference case builds `phi` as a face-centre quadrature `U_f & S_f`, which
is **not** discretely solenoidal, and repairs it with `CorrectPhi` -- a Poisson
projection. That is workable for a fixed field, but under a moving frame the
flux changes shape every step, so the projection would have to be re-solved
every step *and* it perturbs the field away from the analytic one.

The base field admits a closed-form vector potential, so the projection is
unnecessary:

    A0 = psi0(x,y) e_z + (1/2 - 4 rho/3 + rho^2) ( e_z x (p_h - s) )

with psi0 = (1/pi) sin^2(pi x) sin^2(pi y) -- the same stream function as the
2D test -- and rho the distance from the swirl axis s. Then
(curl A0)_z = 1 - 4 rho + 4 rho^2 = (1 - 2 rho)^2, the axial velocity, and the
1/rho of e_theta cancels so A0 is regular on the axis.

Fluxes are the discrete circulation `phi_f = sum 1/2 (A_p1 + A_p2).(x_p2 - x_p1)`
around each face loop, which telescopes to zero over a closed cell.

Measured at N=32 with plicRDF: the potential route gives
`E_shape(T) = 4.477e-03` against the reference `4.483e-03` (0.13% apart), with
`max|div(u)| = 1.3e-13` and volume conserved to 1e-15 -- no Poisson solve.

## Frame choice

The frame rotates about an axis **parallel to e_z** through `(0.5, 0.5)`. The
horizontal motion is two-dimensionally solenoidal with a scalar amplitude, so
`psi0` is constant along trajectories and the material's horizontal support is
a streamline band exactly as in 2D; a z-parallel rotation leaves the axial
drift, which is what the tall `[0,2]` box exists for, untouched. A rotation
about a horizontal axis would sweep that extent through the domain and need a
far larger mesh.

## Dependencies

The snapshot rule needs `scikit-image` (marching cubes) in addition to
numpy/matplotlib/foamlib:

```bash
module load python && python3 -m pip install --user scikit-image
```
