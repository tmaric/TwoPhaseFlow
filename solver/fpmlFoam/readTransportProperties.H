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

    dimensionedScalar f
    (
        "f",
        dimless/dimTime,
        transportProperties
    );

    Info<< "Creating PML coefficients fields \n" << endl;
    // PML thickness
    dimensionedScalar L
    (
        "L", dimless, transportProperties
    );
    // PML damping profile parameters (COMSOL-like polynomial ramp)
    dimensionedScalar sigmaMax
    (
        "sigmaMax", dimless/dimTime, transportProperties
    );
    dimensionedScalar p
    (
      "p", dimless, transportProperties
    );
    // CFS-PML frequency shift parameter
    dimensionedScalar alphaPml
    (
      "alphaPml", dimless/dimTime, transportProperties
    );
    // global damping shift for high-frequency stabilization
    dimensionedScalar sigmaShift
    (
      "sigmaShift", dimless/dimTime, transportProperties
    );

    // Cylindrical PML parameters
    dimensionedScalar Rmin
    (
        "Rmin", dimLength, transportProperties
    );
    dimensionedScalar Rmax
    (
        "Rmax", dimLength, transportProperties
    );
    dimensionedVector origin
    (
        "origin", dimLength, transportProperties
    );
