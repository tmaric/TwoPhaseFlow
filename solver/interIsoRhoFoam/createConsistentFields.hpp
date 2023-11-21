const dimensionedScalar nu1 = mixture.nuModel1().viscosityProperties().get<dimensionedScalar>("nu");
const dimensionedScalar nu2 = mixture.nuModel2().viscosityProperties().get<dimensionedScalar>("nu");
surfaceScalarField muf(mixture.muf());
surfaceScalarField rhof("rhof", fvc::interpolate(rho));

surfaceScalarField alphaf
(
    IOobject
    (
        "alphaf",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
     ), 
     fvc::interpolate(alpha1)
);

volScalarField rhoFromAlphaf
(
    IOobject
    (
        "rhoFromAlphaf",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    rho
);

//rho.write(); 
//rhoFromAlphaf.write();

// cutFaceAdvect cutFace(alpha1.mesh(), alpha1);
// surfaceIteratorPLIC sIterPLIC(alpha1.mesh(), 1e-8);

//fractional volumetric flux
surfaceScalarField alphaPhi
(
    IOobject
    (
        "alphaPhi",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar(dimVol/dimTime, Zero)
);

//temporary alpha 
volScalarField tAlpha
(
    IOobject
    (
        "tAlpha",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ), 
    alpha1
);

