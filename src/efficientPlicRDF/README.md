# efficientPlicRDF

Reconstruction scheme `efficientPlicRDF` for OpenFOAM's geometric VoF (`interIsoFoam`,
`libgeometricVoF`). It is OpenFOAM's `plicRDF` with the same algorithm and settings, and
parallel communication that does not scale with the number of processors:

- the cell-point-cell stencil of the interface band is updated without communication
  (`efficientZoneCPCStencil`: neighbours of boundary points are synchronised once per mesh);
- stencil data are exchanged with the stencil neighbour processors only, without a global
  reduction per exchange (`efficientZoneDistribute`);
- the band marking synchronises coupled points through `globalMeshData` instead of a
  reduction on the master (`efficientReconstructedDistanceFunction`).

Results agree with `plicRDF` up to the order of the stencil entries (round-off).

## Build

Compiles with OpenFOAM v2506, v2512 and v2606:

    wmake libso src/efficientPlicRDF

The library links OpenFOAM's `libgeometricVoF`. Use it with OpenFOAM's `interIsoFoam`, not
with the TwoPhaseFlow solvers that load `libVoF`: both libraries define the same
reconstruction classes.

## Usage

`system/controlDict`:

    libs (efficientPlicRDF);

`system/fvSolution`:

    "alpha.water.*"
    {
        advectionScheme      isoAdvection;
        reconstructionScheme efficientPlicRDF;
        ...
    }

Coefficients (`isoFaceTol`, `surfCellTol`, `tol`, `relTol`, `iterations`,
`interpolateNormal`) as for `plicRDF`, optionally in `efficientPlicRDFCoeffs`.
