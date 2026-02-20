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

    // Optional diagonal shift for Helmholtz stabilization (1/length^2)

    Info<< "Creating PML coefficients fields \n" << endl;
    // PML thickness
    dimensionedScalar L
    (
        "L", dimless, transportProperties
    );
    // PML damping profile parameters (used by cylindrical CFS-PML)
    dimensionedScalar sigmaMax
    (
        "sigmaMax",
        dimless/dimTime,
        transportProperties.lookupOrDefault<scalar>("sigmaMax", 0.0)
    );
    dimensionedScalar p
    (
        "p",
        dimless,
        transportProperties.lookupOrDefault<scalar>("p", 2.0)
    );
    // CFS-PML frequency shift parameter
    dimensionedScalar alphaPml
    (
        "alphaPml",
        dimless/dimTime,
        transportProperties.lookupOrDefault<scalar>("alphaPml", 0.0)
    );

    enum PMLType
    {
        RECTANGLE,
        CYLINDRICAL,
        UNDEFINED
    };

    const word PMLTypeWord = transportProperties.get<word>("PMLType");

    PMLType selectedModel = UNDEFINED;

    if (PMLTypeWord == "rectangle")
    {
        selectedModel = RECTANGLE;
    }
    else if (PMLTypeWord == "cylindrical")
    {
        selectedModel = CYLINDRICAL;
    }
    else
    {
        FatalErrorInFunction
            << "Unknown PMLType '" << PMLTypeWord
            << "' in transportProperties." << nl
            << "Valid types are: rectangle, cylindrical."
            << exit(FatalError);
    }

    // PML parameters (resolved by model)
    dimensionedVector PMLMax;
    dimensionedVector PMLMin;
    dimensionedScalar Rmin;
    dimensionedScalar Rmax;
    dimensionedVector origin;

    switch (selectedModel)
    {
        case RECTANGLE:
        {
            transportProperties.lookup("PMLMax") >> PMLMax;
            transportProperties.lookup("PMLMin") >> PMLMin;
            break;
        }
        case CYLINDRICAL:
        {
            transportProperties.lookup("Rmin") >> Rmin;
            transportProperties.lookup("Rmax") >> Rmax;
            transportProperties.lookup("origin") >> origin;
            break;
        }
        case UNDEFINED:
        {
            FatalErrorInFunction
                << "Internal error: PML type is undefined."
                << exit(FatalError);
            break;
        }
    }
