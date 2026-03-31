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
                  - fvc::snGrad(pr)*beta
                ) * mesh.magSf()
            )
        );

        fvOptions.correct(U);
    }
