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

    volVectorField HbyA(constrainHbyA(rAU()*UEqn.H(), U, p_rgh));

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

    phiHbyA += phig;

    if (args.get<word>("nonOrthogControl") == "PIMPLE")
    {
         constrainPressure(p_rgh, U, phiHbyA, rAUf, MRF);

        while (pimple.correctNonOrthogonal())
        {
            // constrainPressure(p_rgh, U, phiHbyA, rAUf, MRF);
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
        }
    }

    if (args.get<word>("nonOrthogControl") == "ERCNOC")
    {
        constrainPressure(p_rgh, U, phiHbyA, rAUf, MRF);
        label Nnon = 0;
        bool correctNonOrthogonality = true;

        while ( correctNonOrthogonality && pimple.correctNonOrthogonal() )
        {
            ++Nnon;
            // constrainPressure(p_rgh, U, phiHbyA, rAUf, MRF);
            
            fvScalarMatrix p_rghEqn
            (
                fvm::laplacian(rAUf, p_rgh) == fvc::div(phiHbyA)
            );
            p_rghEqn.setReference(pRefCell, getRefCellValue(p_rgh, pRefCell));
            p_rghEqn.solve(mesh.solver(p_rgh.select(pimple.finalInnerIter())));
            // const solverPerformance& p_rgh_performance = p_rghEqn.solve(mesh.solver(p_rgh.select(pimple.finalInnerIter())));
            const List<solverPerformance> sp = mesh.data().solverPerformanceDict().lookup("p_rgh");
            correctNonOrthogonality = (sp.last().nIterations() > 0);//(p_rgh_performance.nIterations() > 0); 

            if ((!correctNonOrthogonality) || pimple.finalNonOrthogonalIter())
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
