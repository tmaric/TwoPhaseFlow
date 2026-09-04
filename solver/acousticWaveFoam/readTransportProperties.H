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

/*    dimensionedScalar mu
    (
//        transportProperties.lookup("nu")
     "mu",
     dimViscosity*dimDensity,
     transportProperties
    );*/

     dimensionedScalar rhog
     (
 //        transportProperties.lookup("nu")
      "rhog",
      dimDensity,
      transportProperties
     );

     dimensionedScalar rhol
     (
 //        transportProperties.lookup("nu")
      "rhol",
      dimDensity,
      transportProperties
     );

     dimensionedScalar kg
     (
     //    thermodynamicProperties.lookup("psi")
      "kg",
      dimCompressibility/dimDensity,
      transportProperties
     );
     dimensionedScalar kl
     (
     //    thermodynamicProperties.lookup("psi")
      "kl",
      dimCompressibility/dimDensity,
      transportProperties
     );
     dimensionedScalar cg
     (
     //    thermodynamicProperties.lookup("psi")
     "cg",
     dimVelocity,
     transportProperties
     );

     dimensionedScalar cl
     (
     //    thermodynamicProperties.lookup("psi")
      "cl",
      dimVelocity,
      transportProperties
     );

const scalar driveFrequency = transportProperties.get<scalar>("f");
const scalar radiationAveragePeriods =
    transportProperties.getOrDefault<scalar>("radiationAveragePeriods", 5);

const dimensionedScalar acousticDampingRate
(
    "acousticDampingRate",
    dimless/dimTime,
    transportProperties.getOrDefault<scalar>("acousticDampingRate", 0)
);

if
(
    driveFrequency <= 0
 || radiationAveragePeriods <= 0
 || acousticDampingRate.value() < 0
)
{
    FatalIOErrorInFunction(transportProperties)
        << "f and radiationAveragePeriods must be positive and "
        << "acousticDampingRate must be nonnegative, but are "
        << driveFrequency << ", " << radiationAveragePeriods << ", and "
        << acousticDampingRate.value()
        << exit(FatalIOError);
}

Info<< "Physical-domain acoustic damping rate = "
    << acousticDampingRate.value() << " 1/s" << endl;

const scalar laplacianCorrection =
    transportProperties.getOrDefault<scalar>("laplacianCorrection", 1.0);

const word timeIntegrationMethod =
    transportProperties.getOrDefault<word>
    (
        "timeIntegrationMethod",
        laplacianCorrection > 0.5 ? "thetaNewmark" : "fullyImplicit"
    );

if
(
    timeIntegrationMethod != "blockCoupled"
 && timeIntegrationMethod != "thetaNewmark"
 && timeIntegrationMethod != "fullyImplicit"
)
{
    FatalIOErrorInFunction(transportProperties)
        << "timeIntegrationMethod must be blockCoupled, thetaNewmark, or "
        << "fullyImplicit, but is " << timeIntegrationMethod
        << exit(FatalIOError);
}

const bool blockCoupledTimeIntegration =
    timeIntegrationMethod == "blockCoupled";

const word timeBlockLinearSolver =
    transportProperties.getOrDefault<word>
    (
        "timeBlockLinearSolver",
        "direct"
    );

if
(
    timeBlockLinearSolver != "direct"
 && timeBlockLinearSolver != "iterative"
)
{
    FatalIOErrorInFunction(transportProperties)
        << "timeBlockLinearSolver must be direct or iterative, but is "
        << timeBlockLinearSolver << exit(FatalIOError);
}

if
(
    mag(laplacianCorrection) > SMALL
 && mag(laplacianCorrection - 1.0) > SMALL
)
{
    FatalIOErrorInFunction(transportProperties)
        << "laplacianCorrection is a legacy switch and must be 0 or 1, but is "
        << laplacianCorrection << exit(FatalIOError);
}

Info<< "Time integration = " << timeIntegrationMethod << endl;

if (blockCoupledTimeIntegration)
{
    Info<< "Time-block linear solver = " << timeBlockLinearSolver << endl;
}

const scalar timeIntegrationTheta =
    transportProperties.getOrDefault<scalar>("timeIntegrationTheta", 0.5);

if (timeIntegrationTheta < 0.5 || timeIntegrationTheta > 1)
{
    FatalIOErrorInFunction(transportProperties)
        << "timeIntegrationTheta must be between 0.5 and 1, but is "
        << timeIntegrationTheta << exit(FatalIOError);
}

Info<< "Time-integration theta = " << timeIntegrationTheta << endl;

scalar MaxCo =
    max(mesh.surfaceInterpolation::deltaCoeffs()*cg).value()
   *runTime.deltaT().value();

Info<< "Max acoustic Courant Number = " << MaxCo << endl;
