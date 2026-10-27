# penetrableSphere case

Controlled comparison of the interface face-area definitions on a curved
transmitting interface. A water sphere (1000 kg/m^3, 1500 m/s) of radius
`a = 1 mm` in air (1.2 kg/m^3, 343 m/s) is centred in the cube
`[-2a, 2a]^3` and enters only through `setAlphaField` volume fractions. An
incident plane wave `exp(i k x)` with `ka = 1` in the gas is scattered by the
sphere. The exact total pressure (regular interior and outgoing exterior
spherical waves, orders 0-24) is prescribed on the cube faces; no PML is
used.

Both modes keep the pressure flux (`flux legacy`) and differ only in the
interface face area:

- `legacy`: `areaFraction legacy`, the original donor-based selection;
- `geometric`: `areaFraction plicAverage`, geometric averaging.

The PML variant of this test is `../penetrableSpherePML`.

## Files and roles

- `sphere_area_comparison.py`: generates, runs and evaluates the cases
  (`blockMesh`, `setAlphaField`, `acousticHelmholtzFoam`). Uses the binaries
  in `PATH`.
- `Allrun`: the full study (`16,24,32,48` cells per direction, serial and
  eight ranks); `Allclean` removes `studyResults/`.

## Usage

Load the OpenFOAM and TwoPhaseFlow environment, then

```bash
./Allrun
```

`SPHERE_SIZES` and `SPHERE_RANKS` (default `1,8`, starting with 1) select
the meshes and rank counts.

## Output

`studyResults/` contains `settings.json` (arguments, Git revision, binary
and script hashes), `results.json`, `area-comparison.csv`, the prepared cases
`prepared-N<N>` and the solved cases `N<N>-<mode>-np<ranks>`.
`pressureRelL2` is the volume-weighted relative complex-pressure L2 error
against the series solution; `liquidRelL2` and `gasRelL2` weight it with
`alpha` and `1 - alpha`, `interfaceRelL2` restricts it to mixed cells.
`parallelDifference` is the volume-weighted relative difference of the
reconstructed parallel pressure from the serial pressure of the same mode.
