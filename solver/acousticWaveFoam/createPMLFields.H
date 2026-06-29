Info<< "Reading PML damping field sigma\n" << endl;

volTensorField sigma
(
    IOobject
    (
        "sigma",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
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

volScalarField C0
(
    IOobject
    (
        "C0",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("C0", dimTime/(dimLength*dimLength), Zero)
);

volScalarField C1
(
    IOobject
    (
        "C1",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("C1", dimless/(dimLength*dimLength), Zero)
);

volScalarField C2
(
    IOobject
    (
        "C2",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar
    (
        "C2",
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

volTensorField T1
(
    IOobject
    (
        "T1",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedTensor("T1", dimless/dimTime, tensor::zero)
);

volTensorField T2
(
    IOobject
    (
        "T2",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedTensor("T2", dimless/(dimTime*dimTime), tensor::zero)
);

Info<< "Creating time-domain PML coefficients from sigma\n" << endl;
