/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2016 DHI
    Copyright (C) 2017 OpenCFD Ltd.
    Copyright (C) 2018 Johan Roenby
    Copyright (C) 2019-2020 DLR
    Copyright (C) 2020 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    interIsoFoam

Group
    grpMultiphaseSolvers

Description
    Solver derived from interFoam for two incompressible, isothermal immiscible
    fluids using the isoAdvector phase-fraction based interface capturing
    approach, with optional mesh motion and mesh topology changes including
    adaptive re-meshing.

    Reference:
    \verbatim
        Roenby, J., Bredmose, H. and Jasak, H. (2016).
        A computational method for sharp interface advection
        Royal Society Open Science, 3
        doi 10.1098/rsos.160405
    \endverbatim

    isoAdvector code supplied by Johan Roenby, STROMNING (2018)

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "isoAdvection.H"
#include "cutFaceAdvect.H"
#include "surfaceIteratorPLIC.H"
#include "EulerDdtScheme.H"
#include "localEulerDdtScheme.H"
#include "CrankNicolsonDdtScheme.H"
#include "subCycle.H"
#include "immiscibleIncompressibleTwoPhaseMixture.H"
#include "incompressibleInterPhaseTransportModel.H"
#include "pimpleControl.H"
#include "fvOptions.H"
#include "CorrectPhi.H"
#include "fvcSmooth.H"
#include "dynamicRefineFvMesh.H"
#include "reconstructionSchemes.H"
#include "upwind.H"
#include "fvMesh.H"
#include "processorPolyPatch.H"
#include "processorBC.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Solver for two incompressible, isothermal immiscible fluids"
        " using isoAdvector phase-fraction based interface capturing.\n"
        "With optional mesh motion and mesh topology changes including"
        " adaptive re-meshing.\n"
        "The solver is derived from interFoam"
    );

    argList::addOption
    (
        "tScheme",
        "time scheme name",
        "Euler, CrankNicolson"
    );

    // TODO (TT): Does not work with CMake yet.
    #include "postProcess.H"

    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createDynamicFvMesh.H"
    #include "initContinuityErrs.H"
    #include "createDyMControls.H"
    #include "createFields.H"
    #include "createConsistentFields.hpp"
    #include "initCorrectPhi.H"
    #include "createUfIfPresent.H"

    #include "CourantNo.H"
    #include "setInitialDeltaT.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readDyMControls.H"
        #include "CourantNo.H"
        #include "alphaCourantNo.H"
        #include "setDeltaT.H"

        ++runTime;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        // TODO: double-check, discussion points
        // - check what oldTime() does here really
        rhoPhi.oldTime() == rhoPhi; 
        //tAlpha.oldTime() == alpha1;
        // #include "computeRhof.H"

//        #include "alphaEqn.H"  ## set back to orig alpha update approach

        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
            if (pimple.firstIter() || moveMeshOuterCorrectors)
            {
                if (isA<dynamicRefineFvMesh>(mesh))
                {
                    advector.surf().reconstruct();
                }

                mesh.update();

                if (mesh.changing())
                {
                    gh = (g & mesh.C()) - ghRef;
                    ghf = (g & mesh.Cf()) - ghRef;

                    if (isA<dynamicRefineFvMesh>(mesh))
                    {
                        advector.surf().mapAlphaField();
                        alpha2 = 1.0 - alpha1;
                        alpha2.correctBoundaryConditions();
                        rho == alpha1*rho1 + alpha2*rho2;
                        rho.correctBoundaryConditions();
                        rho.oldTime() = rho;
                        alpha2.oldTime() = alpha2;
                    }

                    MRF.update();

                    if (correctPhi)
                    {
                        // Calculate absolute flux
                        // from the mapped surface velocity
                        phi = mesh.Sf() & Uf();

                        #include "correctPhi.H"

                        // Make the flux relative to the mesh motion
                        fvc::makeRelative(phi, U);

                        mixture.correct();
                    }

                    if (checkMeshCourantNo)
                    {
                        #include "meshCourantNo.H"
                    }
                }
            }


            #include "alphaControls.H"
            #include "alphaEqnSubCycle.H" //"alphaEqnSubCycle.H"
            // Jun ,testen bitte

            // Equation 21 in Overleaf for Euler integration, for n
            // What do we want: 
            // symbolic: rhoPhi^{n+1} = (rho1 - rho2)*alphaf^{n+1}*phi^m + rho2*phi^m;
            //
            // code rhoPhi = (rho1 - rho2)*alphaf*phi + rho2*phi;
            // 
            // what we have:
            //
//             rhoPhi = (rho1 - rho2)*advector.alphaPhi() + rho2*phi;
            #include "computeRhof.H"

// ######   compare the difference between alphaPhi from isoAdvector and alphaf*phi
/*     	    const auto& alphafPhi = alphaf*phi;
            DynamicList<scalar> diff;
            forAll(phi, fi)
            {
                if (mag(phi[fi]) > SMALL)
                {
                    diff.append((advector.alphaPhi()[fi]-alphafPhi.ref()[fi])/phi[fi]);
                }
                
            }
            Info << "### The difference between alphaPhi/phi and alphaf: " << gMax(mag(diff)) << endl;
            diff.clear();*/
/*            if (pimple.nCorrPIMPLE() > 1)
            {
                if (!pimple.firstIter())
                {
                    rhoPhi = rhof * 0.5*(phi.prevIter() + phi);
                }
            }
*/
//            rhoPhi == rhof * 0.5*(phi.prevIter() + phi); //phi; 

            // TODO: Mixture update? 

            // in Momentum equation, convective term:
            // rho_f^{n+1} F_f^m v_f^{n+1} - implict Euler 

            // Solve for \rho_c^{n+1} using the mass conservation equation.
            if (args.get<word>("tScheme") == "Euler")
            {
               // Option 1
               fvScalarMatrix rhoEqn
               (
                     fvm::ddt(rho) + fvc::div(rhoPhi)
               );
               rhoEqn.solve();
               
               // Option 2: besser, weil wir keinen Eintrag im system/fvSolution brauchen
               // TODO: rho.oldTime() == rho; 
               // rho == rho.oldTime() - runTime.time().deltaT()*fvc::surfaceIntegrate(rhoPhi); 
               // auskommentieren wenn oben rho = rho.oldTime()... steht.
               //rho.correctBoundaryConditions();
            }

	    Info<< "### density from rhoEqn = "
                << "  Min(" << rho.name() << ") = " << min(rho).value()
                << "  Max(" << rho.name() << ") = " << max(rho).value()
                << endl;

/*            rhoFromAlphaf == rho;
            tAlpha == (rho-rho2)/(rho1-rho2); //temporary alpha field used to bound rhoFromAlphaf
//	    Info << "### Conservation before bounding: " << fvc::domainIntegrate(1-tAlpha) << endl;
            #include "alphaSuSp.H"
//            Info << "### Conservation after bounding: " << fvc::domainIntegrate(1-tAlpha) << endl;
            tAlpha.correctBoundaryConditions();
            rho == tAlpha*(rho1-rho2) + rho2;

            Info<< "### bounded density from rhoEqn = "
                << "  Min(" << rho.name() << ") = " << min(rho).value()
                << "  Max(" << rho.name() << ") = " << max(rho).value()
                << endl;

            Info<< "### difference beween bounded density and orig densityfrom rhoEqn = "
                << "  Min(rho_diff) = " << min(rhoFromAlphaf-rho).value()
                << "  Max(rho_diff) = " << max(rhoFromAlphaf-rho).value()
                << endl;*/ 
            // If rhoPhi is computed and upated in alphaEqnSubcycle.H
            // there is no need to calculate it explicitly, so 
            // TODO: Remove this 
            //alphaface == mag(cutfaceInfo.subFaceArea())*areaDim/mesh.magSf();
            //rhof == alphaface*rho1+(1-alphaface)*rho2;
            //muf == alphaface*rho1*nu1+(1-alphaface)*rho2*nu2;
            //rhoPhi == rhof * phi; //test new rhoPhi


            if (args.get<word>("tScheme") == "CrankNicolson")
            {
               // Option 1
               // CrankNicolson für RhoEqn
               //fvScalarMatrix rhoEqn
               //(
                     //fvm::ddt(rho) + fvm::div(rhoPhi)
               //);
               //rhoEqn.solve();
               // in system/fvSchemes für dddt(rho) Crank Nicolson auswählen.
               
               // Option 2
               // CrankNicolson für RhoEqn
               rho = rho.oldTime() - 
                   0.5*runTime.time().deltaT()*fvc::surfaceIntegrate(rhoPhi.oldTime() + rhoPhi);
            }

            // mixture.correct() corrects sigma*K - we need this for surface tension.
            mixture.correct();
            // but mixture.correct() also overwrites muf, so we need our muf 
//            muf == alphaf*rho1*nu1 + (1 - alphaf)*rho2*nu2;


            if (pimple.frozenFlow())
            {
                continue;
            }

            // Will UEqn use rhoPhi.oldTime() with CrankNicolson?
            #include "UEqn.hpp"

            // --- Pressure corrector loop
            while (pimple.correct())
            {
                #include "pEqn.H"
            }

            if (pimple.turbCorr())
            {
                turbulence->correct();
            }
        }

//        rhoFromAlphaf == rho;
        // Reset the density using cell-centered volume fractions.
        rho == alpha1*rho1 + (1.0 - alpha1)*rho2;
      //  Info << "max difference between rho und rhoFromAlphaf: " << max(mag(rho - rhoFromAlphaf)) <<endl;
     //   Info << ((rho1 - rho2)*alpha1+rho2 - rhoFromAlphaf).value() <<endl;
        runTime.write();

        runTime.printExecutionTime(Info);
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
