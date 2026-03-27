    MRF.correctBoundaryVelocity(U);

    volVectorField  nHat = advector->normal()/(mag(advector->normal()) + dimensionedScalar(dimArea, ROOTVSMALL));
    pMomFlux = (momFlux & nHat) & nHat;
    pMomFlux.correctBoundaryConditions();

    fvVectorMatrix UEqn
    (
        fvm::ddt(rho, U) + fvm::div(rhoPhi, U)
      + MRF.DDt(rho, U)
      + turbulence->divDevRhoReff(rho, U)
     ==
        fvOptions(rho, U)
    );

    UEqn.relax();

    fvOptions.constrain(UEqn);

    prFlux = fvc::snGrad(alpha1);
    prDebug = fvc::reconstruct(fvc::interpolate(pr)*fvc::snGrad(alpha1));
    prDebug.correctBoundaryConditions();

    if (pimple.momentumPredictor())
    {
        solve
        (
            UEqn
         ==
            fvc::reconstruct
            (
                (
                  surfForces.surfaceTensionForce()
                  + surfForces.accelerationForce()
                  - fvc::snGrad(p_rgh)
                  + fvc::interpolate(pr)*fvc::snGrad(alpha1)
                  + fvc::interpolate(pMomFlux)*fvc::snGrad(alpha1)
                ) * mesh.magSf()
            )
        );

        fvOptions.correct(U);
    }
