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
    /*surfaceScalarField rAUf(
                            IOobject
                            (
                            "rAUf",
                            runTime.timeName(),
                            mesh,
                            IOobject::NO_READ,
                            IOobject::NO_WRITE
                            ),
                            mesh,
                            dimensionedScalar("rAUf", fvc::interpolate(rAU()).ref().dimensions(), 1.0));*/

    volVectorField HbyA(constrainHbyA(rAU()*UEqn.H(), U, p_rgh));
/*    Info << "##test HbyA, max " << max(HbyA) << ", min " << min(HbyA) << endl;
    Info << "##test rAU, max " << max(rAU()) << ", min " << min(rAU()) << endl;
    Info << "##test UEqn diag(), max " << max(UEqn.diag()) << ", min " << min(UEqn.diag()) << endl;
    Info << "##test UEqn diag(), max " << max(UEqn.D()) << ", min " << min(UEqn.D()) << endl;*/

    surfaceScalarField phiHbyA
    (
        "phiHbyA",
        fvc::flux(HbyA)
      + MRF.zeroFilter(fvc::interpolate(rho*rAU())*fvc::ddtCorr(U, phi, Uf))
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
        )*rAUf*mesh.magSf()
    );
    //Info << "##test rAUf, max " << max(rAUf) << ", min " << min(rAUf) << endl;
    //Info << "##test magSf, max " << max(mesh.magSf()) << ", min " << min(mesh.magSf()) << endl;
  //  Info << "##test st force, max " << max(surfForces.surfaceTensionForce()) << ", min " << min(surfForces.surfaceTensionForce()) << endl;
    // Info << "##test phiHbyA, max " << max(phiHbyA) << ", min " << min(phiHbyA) << endl;
    phiHbyA += phig;

//    volScalarField divPhig("divPhig", fvc::div(phig));
//    divPhig.write();

//    Info << "##test div(phig), max " << max(fvc::div(phig)).value() << ", min " << min(fvc::div(phig)).value() << endl;
    // Update the pressure BCs to ensure flux consistency
    constrainPressure(p_rgh, U, phiHbyA, rAUf, MRF);
    //Info << "##test div(phiHbyA), max " << max(fvc::div(phiHbyA)) << ", min " << min(fvc::div(phiHbyA)) << endl;
    if (args.get<word>("nonOrthogControl") == "PIMPLE")
    {
        while (pimple.correctNonOrthogonal())
        {
            fvScalarMatrix p_rghEqn
            (
                fvm::laplacian(rAUf, p_rgh) == fvc::div(phiHbyA)
            );

            p_rghEqn.setReference(pRefCell, getRefCellValue(p_rgh, pRefCell));

            p_rghEqn.solve(mesh.solver(p_rgh.select(pimple.finalInnerIter())));
            // Info << "##test (phig - p_rghEqn.flux()), max " << max(phig - p_rghEqn.flux()).value() << ", min " << min(phig - p_rghEqn.flux()).value() << endl;

            if (pimple.finalNonOrthogonalIter())
            {
                phi = phiHbyA - p_rghEqn.flux();
                //Info << "##test div(phi), max " << max(fvc::div(phi)).value() << ", min " << min(fvc::div(phi)).value() << endl;
                p_rgh.relax();

                // Info << "##test (phig - p_rghEqn.flux())/rAUf, max " << max((phig - p_rghEqn.flux())/rAUf).value() << ", min " << min((phig - p_rghEqn.flux())/rAUf).value() << endl;
                U = HbyA + rAU()*fvc::reconstruct((phig - p_rghEqn.flux())/rAUf);
                U.correctBoundaryConditions();
                fvOptions.correct(U);
            }
        }
    }

    if (args.get<word>("nonOrthogControl") == "ERCNOC")
    {
        label Nnon = 0;
        bool correctNonOrthogonality = true;

        while (correctNonOrthogonality)
        {
            ++Nnon;
            fvScalarMatrix p_rghEqn
            (
                fvm::laplacian(rAUf, p_rgh) == fvc::div(phiHbyA)
            );

            p_rghEqn.setReference(pRefCell, getRefCellValue(p_rgh, pRefCell));

            auto p_rgh_performance = p_rghEqn.solve(mesh.solver(p_rgh.select(pimple.finalInnerIter())));
            correctNonOrthogonality = (p_rgh_performance.nIterations() > 0); 
            // p_rghEqn.solve(mesh.solver(p_rgh.select(pimple.finalInnerIter())));

            if (!correctNonOrthogonality)
            {
                phi = phiHbyA - p_rghEqn.flux();
                p_rgh.relax();

                U = HbyA + rAU()*fvc::reconstruct((phig - p_rghEqn.flux())/rAUf);
                U.correctBoundaryConditions();
                fvOptions.correct(U);
            }
        }

        Info << "Nnon " << Nnon << endl;
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

    if (!correctPhi)
    {
        rAU.clear();
    }
}
