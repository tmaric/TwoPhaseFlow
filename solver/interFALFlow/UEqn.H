    MRF.correctBoundaryVelocity(U);

    volVectorField  nHat = advector->normal()/(mag(advector->normal()) + dimensionedScalar(dimArea, ROOTVSMALL));
    volScalarField  pMomFlux = (momFlux & nHat) & nHat;

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

    forAll(U, i)
    {
        Fr[i] = pr[i]*advector->normal()[i] / mesh.V()[i];
    }

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
                ) * mesh.magSf()
            )
        );

        fvOptions.correct(U);
    }
