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

const dimensionedScalar sigmaPeak
(
    "sigmaPeak",
    sigma.dimensions(),
    max(max(mag(sigma)).value(), SMALL)
);

const volScalarField pmlImplicitWeight
(
    IOobject
    (
        "pmlImplicitWeight",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE,
        false
    ),
    min(mag(sigma)/sigmaPeak, scalar(1))
);

const volScalarField physicalTimeWeight
(
    IOobject
    (
        "physicalTimeWeight",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE,
        false
    ),
    1.0 - pmlImplicitWeight
);

const surfaceScalarField pmlImplicitWeightf
(
    IOobject
    (
        "pmlImplicitWeightf",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE,
        false
    ),
    fvc::interpolate(pmlImplicitWeight)
);

const surfaceScalarField physicalTimeWeightf
(
    IOobject
    (
        "physicalTimeWeightf",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE,
        false
    ),
    1.0 - pmlImplicitWeightf
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

volScalarField pmlUx
(
    IOobject
    (
        "pmlUx",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE,
        false
    ),
    mesh,
    dimensionedScalar("pmlUx", u.dimensions(), Zero),
    fixedValueFvPatchScalarField::typeName
);

volScalarField pmlUy
(
    IOobject
    (
        "pmlUy",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE,
        false
    ),
    mesh,
    dimensionedScalar("pmlUy", u.dimensions(), Zero),
    fixedValueFvPatchScalarField::typeName
);

volScalarField pmlUz
(
    IOobject
    (
        "pmlUz",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE,
        false
    ),
    mesh,
    dimensionedScalar("pmlUz", u.dimensions(), Zero),
    fixedValueFvPatchScalarField::typeName
);

const surfaceScalarField SfX
(
    IOobject
    (
        "SfX",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE,
        false
    ),
    mesh.Sf().component(vector::X)
);

const surfaceScalarField SfY
(
    IOobject
    (
        "SfY",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE,
        false
    ),
    mesh.Sf().component(vector::Y)
);

const surfaceScalarField SfZ
(
    IOobject
    (
        "SfZ",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE,
        false
    ),
    mesh.Sf().component(vector::Z)
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
    dimensionedScalar("v", dimPressure*dimTime, Zero),
    fixedValueFvPatchScalarField::typeName
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
