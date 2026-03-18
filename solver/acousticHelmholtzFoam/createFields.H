    #include "readTransportProperties.H"
    Info<< "Reading field Pre and Pim\n" << endl;

    volScalarField Pre
    (
        IOobject
        (
            "Pre",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    volScalarField Pim
    (
        IOobject
        (
            "Pim",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    volScalarField pa
    (
        IOobject
        (
            "pa",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("pa", dimPressure, 0)
    );


    Info<< "Reading alpha field \n" << endl;
    volScalarField alpha1
    (
        IOobject
        (
            "alpha.water",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar(dimless, Zero)
    );

    surfaceScalarField alphaf
    (
        IOobject
        (
            "alphaf",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        fvc::interpolate(alpha1)
    );

    // create 0 U for instantiate isoAdvector object.
    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedVector("U", dimVelocity, vector::zero)
    );

    // create real and imaginary part of velocity.
    volVectorField Ure
    (
        IOobject
        (
            "Ure",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedVector("Ure", dimVelocity, vector::zero)
    );

    volVectorField Uim
    (
        IOobject
        (
            "Uim",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedVector("Uim", dimVelocity, vector::zero)
    );

    volScalarField rho
    (
        IOobject
        (
            "rho",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        rhol*alpha1 + rhog*(1 - alpha1)
    );

    // radiation force pressure
    volScalarField pr
    (
        IOobject
        (
            "pr",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("pr", dimPressure, 0)
    );

    volTensorField momFlux
    (
        IOobject
        (
            "momFlux",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedTensor("momFlux", dimDensity*dimVelocity*dimVelocity, tensor::zero)
    );

    volScalarField SC0
    (
        IOobject
        (
            "SC0",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    volScalarField SC1
    (
        IOobject
        (
            "SC1",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    volScalarField TC1
    (
        IOobject
        (
            "TC1",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    volScalarField sigma
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

    #include "createPhi.H"

    #include "createFvOptions.H"

    Foam::advection::isoAdvection advector(alpha1, phi, U);

    cutFaceAdvect cutFace(alpha1.mesh(), alpha1);
