/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2020 OpenCFD Ltd.
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
   acousticFoam

Group
    grpAcousticSolvers

Description
    Acoustic solver solving the acoustic pressure wave equation.

    \f[
        \ddt2{pa} - c^2 \laplacian{pa} = 0
    \f]

    where
    \vartable
        c       | Sound speed
        pa      | Acoustic pressure
    \endvartable

SourceFiles
    acousticFoam.C

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "fvOptions.H"
#include "pimpleControl.H"
#include "isoAdvection.H"
#include "cutFaceAdvect.H"
#include "surfaceIteratorPLIC.H"
#include "reconstructionSchemes.H"
#include "upwind.H"
#include "processorPolyPatch.H"
#include "processorBC.H"
#include "emptyPolyPatch.H"
#include "wedgePolyPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Acoustic solver solving the acoustic pressure wave equation."
    );

    #include "postProcess.H"

    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createControl.H"
//    #include "readTransportProperties.H"
    #include "createFields.H"
    #include "createPMLFields.H"


    Info<< "\nStarting time loop\n" << endl;
    #include "computeAlphaf.H"
    #include "computePMLCoefs.H"
    float iter=0.0;
    
    while (runTime.run())
    {
        ++runTime;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        while (pimple.correct())
        {
            #include "pEqn.H"
        }

        iter=iter+1.0;
//        AcousticForcePotential == (iter - 1)/iter * AcousticForcePotential + 1/iter * (2*Foam::constant::mathematical::pi*pow(R,3)*( 1/(3.0*rho*sqr(cg))*p*p - rho/2*(U&U)));
        pr == (iter - 1)/iter * pr + 1/iter * ( 1.0/(2.0*rho*sqr(cg))*p*p - rho/2.0*(U&U));
        momFlux == (iter - 1)/iter * momFlux + 1/iter * ( rho*(U*U));
        runTime.write();

        runTime.printExecutionTime(Info);
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
