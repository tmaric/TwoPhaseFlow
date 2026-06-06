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
    acousticHelmholtzSerialFoam

Author
    Chuanchao Xu, MMA, TU Darmstadt
    Email: xu@mma.tu-darmstadt.de

Description
    Block-coupled frequency-domain acoustic solver for a single, non-
    decomposed mesh. The real/imaginary pressure pair is solved as one
    2x2 block Helmholtz system:

        [ A  -(B1 + B2) ] [Pim] = [bPim]
        [ (B1 + B2)  A ] [Pre]   [bPre]

    where A is the main Helmholtz operator and B1/B2 are PML-induced
    coupling operators. Matrix assembly is local/serial in OpenFOAM and
    the PETSc solve uses a direct backend (default: MUMPS via LU).

    This solver intentionally aborts in decomposed MPI runs. Use it as the
    serial reference for equation form and damping behavior.
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

static inline scalar twoPi() { return constant::mathematical::twoPi; }

// Keep PETSc/helper code
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

    if (Pstream::parRun())
    {
        FatalErrorInFunction
            << "acousticHelmholtzSerialFoam uses serial assembly. "
            << "Run on one MPI rank and use OMP_NUM_THREADS for multicore."
            << exit(FatalError);
    }

    // Global indexing for block system
    globalIndex globalCells(mesh.nCells());
    const PetscInt N      = (PetscInt)globalCells.size();
    const PetscInt nLocal = (PetscInt)mesh.nCells();

    Mat M;
    Vec x, b;
    KSP ksp;
    bool dumpedMatrixStats = false;

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

            // Coupling operators built from the opposite field:
            // Pim-row couples to Pre via (laplacian(T1,Pre) + Sp(SC1,Pre))
            fvScalarMatrix couplingLaplPre(fvm::laplacian(T1, Pre));
            fvScalarMatrix couplingMassPre(fvm::Sp(SC1, Pre));

            // Pre-row couples to Pim via (laplacian(T1,Pim) + Sp(SC1,Pim))
            fvScalarMatrix couplingLaplPim(fvm::laplacian(T1, Pim));
            fvScalarMatrix couplingMassPim(fvm::Sp(SC1, Pim));

            assembleBlockSystem
            (
                M, globalCells, N,
                AopPim, AopPre,
                couplingLaplPre, couplingMassPre,
                couplingLaplPim, couplingMassPim
            );

            // RHS from source terms (includes BC contributions)
            scalarField bPim(sourceWithBoundary(AopPim));
            scalarField bCouplingLaplPre(sourceWithBoundary(couplingLaplPre));
            scalarField bCouplingMassPre(sourceWithBoundary(couplingMassPre));
            bPim -= bCouplingLaplPre;
            bPim -= bCouplingMassPre;

            scalarField bPre(sourceWithBoundary(AopPre));
            scalarField bCouplingLaplPim(sourceWithBoundary(couplingLaplPim));
            scalarField bCouplingMassPim(sourceWithBoundary(couplingMassPim));
            bPre += bCouplingLaplPim;
            bPre += bCouplingMassPim;

            setBlockRhs(b, globalCells, N, bPim, bPre);

            MatAssemblyBegin(M, MAT_FINAL_ASSEMBLY);
            MatAssemblyEnd(M, MAT_FINAL_ASSEMBLY);
            VecAssemblyBegin(b);
            VecAssemblyEnd(b);

            if (!dumpedMatrixStats)
            {
                dumpedMatrixStats = true;
                reportMatrixStats(M, b);
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
