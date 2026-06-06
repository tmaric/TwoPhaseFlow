/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           |
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
    acousticHelmholtzFoam

Author
    Chuanchao Xu, MMA, TU Darmstadt
    Email: xu@mma.tu-darmstadt.de

Description
    MPI-capable block-coupled frequency-domain acoustic solver.
    Uses the same pressure block structure as acousticHelmholtzSerialFoam:

        [ A  -(B1 + B2) ] [Pim] = [bPim]
        [ (B1 + B2)  A ] [Pre]   [bPre]

    but assembles A and coupling contributions on decomposed subdomains,
    explicitly including processor-interface couplings in the PETSc matrix.
    PETSc then solves the global distributed linear system (default:
    preonly+lu+mumps).

    Practical difference to acousticHelmholtzSerialFoam:
    - acousticHelmholtzSerialFoam: serial OpenFOAM assembly only (reference).
    - acousticHelmholtzFoam: distributed OpenFOAM assembly + distributed PETSc solve.
\*---------------------------------------------------------------------------*/

#include <petscksp.h>
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
#include "processorLduInterface.H"
#include "processorBC.H"

static inline scalar twoPi() { return constant::mathematical::twoPi; }

#include "diagnosticsHelpers.H"
#include "petscBlockAssembly.H"
#include "petscBlockSolve.H"

int main(int argc, char *argv[])
{
    #include "postProcess.H"
    #include "addCheckCaseOptions.H"
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createFields.H"
    #include "computeAlphaf.H"

    simpleControl simple(mesh);

    PetscInitialize(&argc, &argv, nullptr, nullptr);

    // Global indexing for block system
    globalIndex globalCells(mesh.nCells());
    const PetscInt N      = (PetscInt)globalCells.size();
    const PetscInt nLocal = (PetscInt)mesh.nCells();

    Mat M;
    Vec x, b;
    KSP ksp;
    bool dumpedMatrixStats = false;
    bool dumpedOpStats = false;

    initializePetscSystem(nLocal, N, M, x, b, ksp);

    Info<< "\nStarting time loop\n" << endl;

    rho = 1/(alpha1/rhol + (1 - alpha1)/rhog);
    invRhof = alphaf/rhol + (1 - alphaf)/rhog;

    volScalarField k2
    (
        IOobject
        (
            "k2",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        sqr(twoPi()*f)*(alpha1/sqr(cl) + (1 - alpha1)/sqr(cg))
    );

    while (simple.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        while (simple.correctNonOrthogonal())
        {
            MatZeroEntries(M);
            VecSet(b, 0.0);

            fvScalarMatrix AopPre
            (
              rho*fvm::laplacian(invRhof, Pre)
              + fvm::laplacian(T0, Pre)
              + fvm::Sp(k2 - SC0, Pre)
            );

            fvScalarMatrix AopPim
            (
              rho*fvm::laplacian(invRhof, Pim)
              + fvm::laplacian(T0, Pim)
              + fvm::Sp(k2 - SC0, Pim)
            );

            // Keep B1 (laplacian) and B2 (Sp) as separate operators.
            // Their off-block diagonal handling differs in assembly.
            fvScalarMatrix couplingLaplPre(fvm::laplacian(T1, Pre));  // B1
            fvScalarMatrix couplingMassPre(fvm::Sp(SC1, Pre));         // B2

            fvScalarMatrix couplingLaplPim(fvm::laplacian(T1, Pim));  // B1
            fvScalarMatrix couplingMassPim(fvm::Sp(SC1, Pim));         // B2

            if (!dumpedOpStats)
            {
                reportFvMatrixBreakdown("AopPre beforeBoundaryManipulate", AopPre);
                reportFvMatrixBreakdown("AopPim beforeBoundaryManipulate", AopPim);
                reportFvMatrixBreakdown("couplingLaplPre beforeBoundaryManipulate", couplingLaplPre);
                reportFvMatrixBreakdown("couplingMassPre beforeBoundaryManipulate", couplingMassPre);
            }

            assembleBlockSystem
            (
                M, globalCells, N,
                AopPim, AopPre,
                couplingLaplPre, couplingMassPre,
                couplingLaplPim, couplingMassPim
            );

            if (!dumpedOpStats)
            {
                reportFvMatrixBreakdown("AopPre afterBoundaryManipulate", AopPre);
                reportFvMatrixBreakdown("AopPim afterBoundaryManipulate", AopPim);
                reportFvMatrixBreakdown("couplingLaplPre afterBoundaryManipulate", couplingLaplPre);
                reportFvMatrixBreakdown("couplingMassPre afterBoundaryManipulate", couplingMassPre);
                dumpedOpStats = true;
            }

            scalarField bPim;
            scalarField bPre;
            buildRhs
            (
                AopPim,
                AopPre,
                couplingLaplPre,
                couplingMassPre,
                couplingLaplPim,
                couplingMassPim,
                bPim,
                bPre
            );

            setBlockRhs(b, globalCells, N, bPim, bPre);

            MatAssemblyBegin(M, MAT_FINAL_ASSEMBLY);
            MatAssemblyEnd(M, MAT_FINAL_ASSEMBLY);
            VecAssemblyBegin(b);
            VecAssemblyEnd(b);

            if (!dumpedMatrixStats)
            {
                dumpedMatrixStats = true;
                reportMatrixStats(M, b, N);
                dumpProcessorInterfaceRows(M, AopPim, globalCells, N);
            }

            KSPSolve(ksp, b, x);
            reportKspStats(ksp, b);
            scatterBlockSolution(x, globalCells, N, Pim, Pre);

            Pre.correctBoundaryConditions();
            Pim.correctBoundaryConditions();
        }

        Ure == 1/(2*constant::mathematical::pi*f*rho) * fvc::grad(Pim);
        Uim == -1/(2*constant::mathematical::pi*f*rho) * fvc::grad(Pre);
        pa == Foam::sqrt(Pim*Pim + Pre*Pre);
        pr == 0.25*(kl*alpha1 + kg*(1-alpha1))*(Pre*Pre + Pim*Pim)
            - 0.25*rho*((Ure&Ure) + (Uim&Uim));
        momFlux == 0.5*rho*(Ure*Ure + Uim*Uim);

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    KSPDestroy(&ksp);
    VecDestroy(&x);
    VecDestroy(&b);
    MatDestroy(&M);
    PetscFinalize();

    Info<< "End\n" << endl;
    return 0;
}

// ************************************************************************* //
