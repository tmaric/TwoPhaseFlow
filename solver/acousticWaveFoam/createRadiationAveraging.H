const scalar radiationAverageDuration =
    radiationAveragePeriods/driveFrequency;
const label radiationWindowSamples = max
(
    label(1),
    label(round(radiationAverageDuration/runTime.deltaTValue()))
);

Info<< "Radiation averaging window = " << radiationAveragePeriods
    << " periods (" << radiationWindowSamples << " time samples)" << endl;

PtrList<volScalarField> prWindow(radiationWindowSamples);
PtrList<volTensorField> momFluxWindow(radiationWindowSamples);
label radiationWindowCount = 0;
label radiationWindowIndex = 0;

volScalarField prWindowSum
(
    IOobject
    (
        "prWindowSum",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE,
        false
    ),
    mesh,
    dimensionedScalar("prWindowSum", dimPressure, Zero)
);

volTensorField momFluxWindowSum
(
    IOobject
    (
        "momFluxWindowSum",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE,
        false
    ),
    mesh,
    dimensionedTensor
    (
        "momFluxWindowSum",
        dimDensity*dimVelocity*dimVelocity,
        tensor::zero
    )
);
