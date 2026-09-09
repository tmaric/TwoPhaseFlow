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
#include "acousticInterface.H"
#include <cmath>

static inline scalar twoPi() { return constant::mathematical::twoPi; }

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
    globalIndex globalCells(mesh.nCells());
    Foam::acoustic::acousticInterface interfaceModel(mesh,globalCells);
    interfaceModel.computeAreas(alpha1,alphaf,phi,cutFace);
    #include "computePMLCoefs.H"

    simpleControl simple(mesh);

    PetscCallAbort(PETSC_COMM_WORLD, PetscInitialize(&argc, &argv, nullptr, nullptr));

    // Global indexing for block system
    interfaceModel.build(sigma,rhol.value(),rhog.value());
    const PetscInt N      = (PetscInt)globalCells.size();
    const PetscInt nLocal = (PetscInt)mesh.nCells();

    Mat M;
    Vec x, b;
    KSP ksp;

    initializePetscSystem(mesh, globalCells, interfaceModel.operators(), nLocal, N, M, x, b, ksp);

    Info<< "\nStarting time loop\n" << endl;

    rho = alpha1*rhol + (1 - alpha1)*rhog;
    compressibility = alpha1*kl + (1 - alpha1)*kg;
    invRhof = 1/(alphaf*rhol + (1 - alphaf)*rhog);
    k2 = sqr(twoPi()*f)*rho*compressibility;
    // Zeroing this coefficient removes both the old implicit face term and
    // its explicit nonorthogonal correction on replaced faces.
    interfaceModel.mask(invRhof);

    while (simple.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        while (simple.correctNonOrthogonal())
        {
            PetscCallAbort(PETSC_COMM_WORLD, MatZeroEntries(M));
            PetscCallAbort(PETSC_COMM_WORLD, VecSet(b, 0.0));

            fvScalarMatrix AopPre
            (
              rho*fvm::laplacian(invRhof, Pre)
              + fvm::laplacian(T0, Pre)
              + fvm::Sp(k2 - C0, Pre)
            );

            fvScalarMatrix AopPim
            (
              rho*fvm::laplacian(invRhof, Pim)
              + fvm::laplacian(T0, Pim)
              + fvm::Sp(k2 - C0, Pim)
            );

            // Keep B1 (laplacian) and B2 (Sp) as separate operators.
            // Their off-block diagonal handling differs in assembly.
            fvScalarMatrix couplingLaplPre(fvm::laplacian(T1, Pre));  // B1
            fvScalarMatrix couplingMassPre(fvm::Sp(C1, Pre));          // B2

            fvScalarMatrix couplingLaplPim(fvm::laplacian(T1, Pim));  // B1
            fvScalarMatrix couplingMassPim(fvm::Sp(C1, Pim));          // B2

            assembleBlockSystem
            (
                M, globalCells,
                AopPim, AopPre,
                couplingLaplPre, couplingMassPre,
                couplingLaplPim, couplingMassPim
            );

            insertTransmissionOperators(M,globalCells,interfaceModel.operators(),rho);

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

            setBlockRhs(b, bPim, bPre);

            PetscCallAbort(PETSC_COMM_WORLD, MatAssemblyBegin(M, MAT_FINAL_ASSEMBLY));
            PetscCallAbort(PETSC_COMM_WORLD, MatAssemblyEnd(M, MAT_FINAL_ASSEMBLY));
            PetscCallAbort(PETSC_COMM_WORLD, VecAssemblyBegin(b));
            PetscCallAbort(PETSC_COMM_WORLD, VecAssemblyEnd(b));

            const scalarField previousPre(Pre.primitiveField());
            const scalarField previousPim(Pim.primitiveField());
            PetscCallAbort(PETSC_COMM_WORLD, KSPSolve(ksp, b, x));
            verifyPetscSolution(M,x,b,ksp);
            scatterBlockSolution(x, globalCells, Pim, Pre);

            Pre.correctBoundaryConditions();
            Pim.correctBoundaryConditions();
            const scalar change = Foam::sqrt
            (
                gSum(sqr(Pre.primitiveField()-previousPre)+sqr(Pim.primitiveField()-previousPim))
               /max(gSum(sqr(Pre.primitiveField())+sqr(Pim.primitiveField())),scalar(VSMALL))
            );
            Info<< "Acoustic nonorthogonal change=" << change << nl;
        }

        Ure == 1/(2*constant::mathematical::pi*f*rho) * fvc::grad(Pim);
        Uim == -1/(2*constant::mathematical::pi*f*rho) * fvc::grad(Pre);
        pa == Foam::sqrt(Pim*Pim + Pre*Pre);
        pr == 0.25*compressibility*(Pre*Pre + Pim*Pim)
            - 0.25*rho*((Ure&Ure) + (Uim&Uim));
        momFlux == 0.5*rho*(Ure*Ure + Uim*Uim);

        runTime.write();
        if (runTime.writeTime()) interfaceModel.write(alphaf,Pre,Pim);

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    PetscCallAbort(PETSC_COMM_WORLD, KSPDestroy(&ksp));
    PetscCallAbort(PETSC_COMM_WORLD, VecDestroy(&x));
    PetscCallAbort(PETSC_COMM_WORLD, VecDestroy(&b));
    PetscCallAbort(PETSC_COMM_WORLD, MatDestroy(&M));
    PetscCallAbort(PETSC_COMM_WORLD, PetscFinalize());

    Info<< "End\n" << endl;
    return 0;
}

// ************************************************************************* //
