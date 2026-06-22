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
        1/(alpha1/rhol + (1 - alpha1)/rhog)
    );

    surfaceScalarField invRhof
    (
        IOobject
        (
            "invRhof",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        alphaf/rhol + (1 - alphaf)/rhog
    );

    volScalarField k2
    (
        IOobject
        (
            "k2",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        sqr(twoPi()*f)*(alpha1/sqr(cl) + (1 - alpha1)/sqr(cg))
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
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("SC0", dimless/sqr(dimLength), Zero)
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
        dimensionedScalar("SC1", dimless/sqr(dimLength), Zero)
    );

    volTensorField T0
    (
        IOobject
        (
            "T0",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedTensor("T0", dimless, tensor::zero)
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
        dimensionedTensor("T1", dimless, tensor::zero)
    );

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

    #include "createPhi.H"

    #include "createFvOptions.H"

    Foam::advection::isoAdvection advector(alpha1, phi, U);

    cutFaceAdvect cutFace(alpha1.mesh(), alpha1);
