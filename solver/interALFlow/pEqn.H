{
    if (correctPhi)
    {
        rAU.ref() = 1.0/UEqn.A();
    }
    else
    {
        rAU = 1.0/UEqn.A();
    }

    surfaceScalarField rAUf("rAUf", fvc::interpolate(rAU()));
    volVectorField HbyA(constrainHbyA(rAU()*UEqn.H(), U, p_rgh)); // adjust BC for HbyA
    surfaceScalarField phiHbyA
    (
        "phiHbyA",
        fvc::flux(HbyA)
       + MRF.zeroFilter(fvc::interpolate(rho*rAU())*fvc::ddtCorr(U, phi, Uf)) // correct deviation of fvc::flux(U) from phi
       // + MRF.zeroFilter(fvc::interpolate(rAU())*fvc::ddtCorr(rho, U, rhoPhi, Uf))
    );
    MRF.makeRelative(phiHbyA);

    if (p_rgh.needReference())
    {
        fvc::makeRelative(phiHbyA, U);
        adjustPhi(phiHbyA, U, p_rgh);
        fvc::makeAbsolute(phiHbyA, U);
    }

    surfaceScalarField phig
    (
        (
          //  mixture.surfaceTensionForce()
          surfForces.surfaceTensionForce()
          //- ghf*fvc::snGrad(rho)
          + surfForces.accelerationForce()
          + fvc::interpolate((prAver+pa))*fvc::snGrad(alpha1)
          + fvc::interpolate(pMomFlux)*fvc::snGrad(alpha1)
        )*rAUf*mesh.magSf()
        //+ rAUf*fvc::flux(Fr)
    );

    phiHbyA += phig;

    // Update the pressure BCs to ensure flux consistency
    constrainPressure(p_rgh, U, phiHbyA, rAUf, MRF);

    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix p_rghEqn
        (
            fvm::laplacian(rAUf, p_rgh) == fvc::div(phiHbyA)
        );

        p_rghEqn.setReference(pRefCell, getRefCellValue(p_rgh, pRefCell));

        p_rghEqn.solve(mesh.solver(p_rgh.select(pimple.finalInnerIter())));

        if (pimple.finalNonOrthogonalIter())
        {
            phi = phiHbyA - p_rghEqn.flux();

            p_rgh.relax();

            U = HbyA + rAU()*fvc::reconstruct((phig - p_rghEqn.flux())/rAUf);
            U.correctBoundaryConditions();
            fvOptions.correct(U);
        }
/*        if(runTime.writeTime())
        {
            p_rghEqn.flux().ref().write();
        }*/
    }

    #include "continuityErrs.H"

    // Correct Uf if the mesh is moving
    fvc::correctUf(Uf, U, phi);

    // Make the fluxes relative to the mesh motion
    fvc::makeRelative(phi, U);

    p == p_rgh + rho*gh;

    if (p_rgh.needReference())
    {
        p += dimensionedScalar
        (
            "p",
            p.dimensions(),
            pRefValue - getRefCellValue(p, pRefCell)
        );
        p_rgh = p - rho*gh;
    }

/*    if(runTime.writeTime())
    {
        rAU().write();
        rAUf.write();
        UEqn.H().ref().write();
        HbyA.write();
        phiHbyA.write();
    }*/

    if (!correctPhi)
    {
        rAU.clear();
    }
}
