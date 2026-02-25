//#include "acousticMixFoam.H"

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
//      + (momFlux & advector->normal()).ref()/mesh.V()
    );

    UEqn.relax();

    fvOptions.constrain(UEqn);

    forAll(U, i)
    {
        Fr[i] = (pr[i]+pa[i])*advector->normal()[i] / mesh.V()[i];
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
                    //mixture.surfaceTensionForce()
                  surfForces.surfaceTensionForce()
                  //- ghf*fvc::snGrad(rho)
                  + surfForces.accelerationForce()
                  - fvc::snGrad(p_rgh)
                  + fvc::interpolate((pr+pa))*fvc::snGrad(alpha1)
                ) * mesh.magSf()
            )
            //+ Fr // interface normal points into the drop
        );

        fvOptions.correct(U);
    }
