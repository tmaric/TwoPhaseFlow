Info<< "Reading PML damping field sigma\n" << endl;

volTensorField sigma
(
    IOobject
    (
        "sigma",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    ),
    mesh
);

volVectorField u
(
    IOobject
    (
        "u",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedVector("u", dimPressure/dimLength, vector::zero)
);

volScalarField v
(
    IOobject
    (
        "v",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("v", dimPressure*dimTime, Zero)
);

volScalarField SC0
(
    IOobject
    (
        "SC0",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("SC0", dimTime/(dimLength*dimLength), Zero)
);

volScalarField SC1
(
    IOobject
    (
        "SC1",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("SC1", dimless/(dimLength*dimLength), Zero)
);

volScalarField SC2
(
    IOobject
    (
        "SC2",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar
    (
        "SC2",
        dimless/(dimLength*dimLength*dimTime),
        Zero
    )
);

volTensorField TC0
(
    IOobject
    (
        "TC0",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedTensor("TC0", dimless/dimTime, tensor::zero)
);

volTensorField TC1
(
    IOobject
    (
        "TC1",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedTensor("TC1", dimless/dimTime, tensor::zero)
);

volTensorField TC2
(
    IOobject
    (
        "TC2",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedTensor("TC2", dimless/(dimTime*dimTime), tensor::zero)
);

Info<< "Creating time-domain PML coefficients from sigma\n" << endl;
