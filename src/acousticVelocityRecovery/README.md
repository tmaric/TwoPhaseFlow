# Acoustic velocity recovery

`acousticHelmholtzFoam` writes the velocity amplitude `Ure`, `Uim` from the
pressure, `V = grad(P)/(i omega rho)`. The cell-centred Gauss gradient divided
by the cell density is wrong beside a gas-liquid interface: next to an
interface on a face (gas cell 2.7 times the exact velocity in the layered 1D
case at NX = 2800) and in mixed cells (0.58 times at NX = 8000).
`reconstructPressureAcceleration.H` reconstructs `(1/rho) grad(p)` from the
face fluxes of the solved pressure operator instead, as `fvc::reconstruct`
does, with density weights in mixed cells. The cell algorithm is the one of
`solver/acousticWaveFoam/reconstructPressureAcceleration.H`.

Select it with `velocityRecovery` in `constant/transportProperties`:

| `velocityRecovery` | Velocity |
|---|---|
| `interfaceCells` (default) | Face-flux recovery in mixed cells and in cells with a face whose volume fraction differs from their own; Gauss gradient elsewhere |
| `allCells` | Face-flux recovery everywhere |
| `gaussGradient` | Gauss gradient everywhere (previous behaviour) |

Without an interface, `interfaceCells` recovers no cell and is identical to
`gaussGradient`. `allCells` halves the velocity error on the warped meshes of
`homogeneousPlaneWave2D`, but next to the curved, non-orthogonal sound-hard
sphere of `WedgeLevitatorDropGorkov` the radiation force deviates from the
Gorkov force by 3.2 % instead of 0.10 % (h_s/a = 1/16); it is therefore not
the default. Mixed cells are cells with `surfCellTol < alpha < 1 - surfCellTol`
of the interface reconstruction. With `flux plicTransmission` the Gauss
gradient is kept, because those interface faces carry their flux in a
separate operator. The pressure solution is not affected by the choice.

Boundary values: on non-coupled patches the cell value is extrapolated and its
normal component is replaced by the face flux, as `fvc::grad` does with the
face gradient.
