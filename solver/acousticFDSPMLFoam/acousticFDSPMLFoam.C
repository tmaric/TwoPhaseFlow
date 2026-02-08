/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
     \\/     M anipulation  |
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
   acousticTDSPMLFoam

Group
    grpAcousticSolvers

Description
    Frequency domain acoustic solver solving the acoustic pressure wave equation and calculating acoustic radiation pressure, 
    considering inhomogenous scattering and non-reflective PML.

Author
    Jun Liu, MMA, TU Darmstadt, 04. Nov. 2025
    Email: liu@mma.tu-darmstadt.de
SourceFiles
    acousticFDSPMLFoam.C

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "fvOptions.H"
#include "simpleControl.H"
#include "emptyPolyPatch.H"
#include "wedgePolyPatch.H"
#include "isoAdvection.H"
#include "cutFaceAdvect.H"
#include "surfaceIteratorPLIC.H"
#include "reconstructionSchemes.H"
#include "upwind.H"
#include "processorPolyPatch.H"
#include "processorBC.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "postProcess.H"
    #include "addCheckCaseOptions.H"
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
//    #include "readTransportProperties.H"
    #include "createFields.H"
    #include "computePMLCoefs.H"
    #include "computeAlphaf.H"
    simpleControl simple(mesh);

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    Info<< "\nStarting time loop\n" << endl;
    float pi = Foam::constant::mathematical::pi; 

    while (simple.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        while (simple.correctNonOrthogonal())
        {
            fvScalarMatrix PimEqn
            (
                sqr(2*pi*f/cg)*fvm::Sp( 1 + ((kl - kg)/kg) * alpha1, Pim) 
                + fvm::laplacian(1 - (rhol - rhog)/rhol*alphaf, Pim)
            );

            solve (
                PimEqn ==
                       (SC0 * Pim) 
                    +  (SC1 * Pre)
                    -  fvc::div(
                       (TC0 & fvc::grad(Pim))
                    -  (TC1 & fvc::grad(Pre)) 
                    )


                );

            

            fvScalarMatrix PreEqn
            (
                sqr(2*pi*f/cg)*fvm::Sp( 1 + ((kl - kg)/kg) * alpha1, Pre) 
                    + fvm::laplacian(1 - (rhol - rhog)/rhol*alphaf, Pre)
            );

            solve( 
                PreEqn ==
                   (SC0 * Pre) 
                -  (SC1 * Pim)
                -  fvc::div(
                   (TC0 & fvc::grad(Pre))
                +  (TC1 & fvc::grad(Pim)) 
                )
            );
            
            Pim.relax();
            Pre.relax();

        }

        Ure == 1/(2*pi*f*rho) * fvc::grad(Pim);
        Uim == -1/(2*pi*f*rho) * fvc::grad(Pre);
        pa == Foam::sqrt(Pim*Pim + Pre*Pre);
        pr == 0.25*(kl*alpha1 + kg*(1-alpha1))*(Pre*Pre + Pim*Pim) - 0.25*rho*((Ure&Ure) + (Uim&Uim));
        momFlux == 0.5*rho*(Ure*Ure + Uim*Uim);
        

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
