    Info<< "Reading transportProperties\n" << endl;

    IOdictionary transportProperties
    (
        IOobject
        (
            "transportProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    dimensionedScalar rhog
    (
        "rhog",
        dimDensity,
        transportProperties
    );

    dimensionedScalar rhol
    (
        "rhol",
        dimDensity,
        transportProperties
    );

    dimensionedScalar kg
    (
        "kg",
        dimCompressibility/dimDensity,
        transportProperties
    );
    dimensionedScalar kl
    (
        "kl",
        dimCompressibility/dimDensity,
        transportProperties
    );

    dimensionedScalar cg
    (
        "cg",
        dimVelocity,
        transportProperties
    );

    dimensionedScalar cl
    (
        "cl",
        dimVelocity,
        transportProperties
    );

    dimensionedScalar f
    (
        "f",
        dimless/dimTime,
        transportProperties
    );
