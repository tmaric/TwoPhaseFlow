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
    freITBCFoam

Group
    AcousticSolvers

Description
    Iterative block-coupled frequency-domain acoustic solver for a single,
    non-decomposed mesh. The real/imaginary pressure pair is solved as one
    2x2 block Helmholtz system:

        [ A  -(B1 + B2) ] [Pim] = [bPim]
        [ (B1 + B2)  A ] [Pre]   [bPre]

    where A is the main Helmholtz operator and B1/B2 are PML-induced
    coupling operators.

    PETSc is used with two assembled operators:
    - M: physical solve operator (original block system)
    - P: shifted preconditioner operator, controlled by
         shiftedLaplacianBeta in transportProperties

    The default preconditioner is a PETSc PCSHELL layered sweep on P with
    optional post-corrections (diagonal and residual/defect-based updates).
    Sweep ordering is selectable through -itbc_order (y or radial).

    This solver is intentionally serial at OpenFOAM level: decomposed MPI
    runs are rejected to preserve one-rank assembly behavior.
\*---------------------------------------------------------------------------*/

#include <petscksp.h>
#include <algorithm>
#include <cstring>
#include <utility>
#include <vector>
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
    #include "computePMLCoefs.H"
    #include "computeAlphaf.H"

    simpleControl simple(mesh);

    PetscInitialize(&argc, &argv, nullptr, nullptr);

    if (Pstream::parRun())
    {
        FatalErrorInFunction
            << "freITBCFoam uses serial assembly. "
            << "Run on one MPI rank and use OMP_NUM_THREADS for multicore."
            << exit(FatalError);
    }

    // Global indexing for block system
    globalIndex globalCells(mesh.nCells());
    const PetscInt N      = (PetscInt)globalCells.size();
    const PetscInt nLocal = (PetscInt)mesh.nCells();

    // Build a geometric sweep order, interleaved by block.
    // -itbc_order y       : sort by y coordinate
    // -itbc_order radial  : sort by distance to origin
    char sweepMode[32] = "radial";
    PetscBool hasSweepMode = PETSC_FALSE;
    PetscOptionsGetString
    (
        nullptr,
        nullptr,
        "-itbc_order",
        sweepMode,
        sizeof(sweepMode),
        &hasSweepMode
    );
    const bool radialSweep = (std::strcmp(sweepMode, "radial") == 0);

    std::vector<std::pair<scalar, PetscInt>> cellOrder;
    cellOrder.reserve(mesh.nCells());
    const vectorField& Cc = mesh.C();
    forAll(Cc, celli)
    {
        const PetscInt gi = globalCells.toGlobal(celli);
        const scalar metric = radialSweep
            ? mag(Cc[celli] - origin.value())
            : Cc[celli].component(vector::Y);
        cellOrder.emplace_back(metric, gi);
    }
    std::sort
    (
        cellOrder.begin(),
        cellOrder.end(),
        [](const std::pair<scalar, PetscInt>& a, const std::pair<scalar, PetscInt>& b)
        {
            return a.first < b.first;
        }
    );

    List<PetscInt> sweepOrder(2*mesh.nCells());
    label sweepi = 0;
    for (const auto& yi : cellOrder)
    {
        sweepOrder[sweepi++] = yi.second;     // Pim block
        sweepOrder[sweepi++] = yi.second + N; // Pre block
    }

    // Build coarse geometric layers for block sweeping preconditioner.
    scalar metricMin = GREAT;
    scalar metricMax = -GREAT;
    forAll(Cc, celli)
    {
        const scalar metric = radialSweep
            ? mag(Cc[celli] - origin.value())
            : Cc[celli].component(vector::Y);
        metricMin = std::min(metricMin, metric);
        metricMax = std::max(metricMax, metric);
    }
    const scalar metricSpan = std::max(metricMax - metricMin, SMALL);
    PetscInt nGeomLayersOpt = 3;
    PetscOptionsGetInt(nullptr, nullptr, "-itbc_layers", &nGeomLayersOpt, nullptr);
    const label nGeomLayers = std::min<label>(mesh.nCells(), std::max<label>(1, (label)nGeomLayersOpt));

    List<DynamicList<PetscInt>> layerBuckets(nGeomLayers);
    forAll(Cc, celli)
    {
        const PetscInt gi = globalCells.toGlobal(celli);
        const scalar metric = radialSweep
            ? mag(Cc[celli] - origin.value())
            : Cc[celli].component(vector::Y);
        label li = floor(((metric - metricMin)/metricSpan)*nGeomLayers);
        li = std::min<label>(std::max<label>(li, 0), nGeomLayers - 1);
        layerBuckets[li].append(gi);
        layerBuckets[li].append(gi + N);
    }

    label nNonEmpty = 0;
    forAll(layerBuckets, li)
    {
        if (!layerBuckets[li].empty()) ++nNonEmpty;
    }

    List<List<PetscInt>> sweepLayers(nNonEmpty);
    label outLi = 0;
    forAll(layerBuckets, li)
    {
        if (!layerBuckets[li].empty())
        {
            sweepLayers[outLi] = layerBuckets[li];
            ++outLi;
        }
    }

    Mat M, P;
    Vec x, b;
    KSP ksp;
    bool dumpedMatrixStats = false;

    initializePetscSystem(nLocal, N, sweepOrder, sweepLayers, M, P, x, b, ksp);

    Info<< "\nStarting time loop\n" << endl;

    const dimensionedScalar k2
    (
        "k2",
        dimless/dimLength/dimLength,
        sqr(twoPi()*f/cg).value()
    );

    while (simple.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        while (simple.correctNonOrthogonal())
        {
            MatZeroEntries(M);
            MatZeroEntries(P);
            VecSet(b, 0.0);

            const volScalarField n2
            (
                "n2",
                1 + ((kl - kg)/kg)*alpha1
            );

            fvScalarMatrix AopPre
            (
                fvm::laplacian(1 - (rhol - rhog)/rhol*alphaf, Pre)
              + fvm::Sp(k2*n2 - SC0, Pre)
            );

            fvScalarMatrix AopPim
            (
                fvm::laplacian(1 - (rhol - rhog)/rhol*alphaf, Pim)
              + fvm::Sp(k2*n2 - SC0, Pim)
            );

            // Coupling operators built from the opposite field:
            // Pim-row couples to Pre via (laplacian(TC1,Pre) + Sp(SC1,Pre))
            fvScalarMatrix couplingLaplPre(fvm::laplacian(TC1, Pre));
            fvScalarMatrix couplingMassPre(fvm::Sp(SC1, Pre));

            // Pre-row couples to Pim via (laplacian(TC1,Pim) + Sp(SC1,Pim))
            fvScalarMatrix couplingLaplPim(fvm::laplacian(TC1, Pim));
            fvScalarMatrix couplingMassPim(fvm::Sp(SC1, Pim));

            // Preconditioner uses C + beta*k^2*n^2 in off-diagonal coupling.
            fvScalarMatrix couplingMassShiftedPre
            (
                fvm::Sp(SC1 + shiftedLaplacianBeta*k2*n2, Pre)
            );
            fvScalarMatrix couplingMassShiftedPim
            (
                fvm::Sp(SC1 + shiftedLaplacianBeta*k2*n2, Pim)
            );

            assembleBlockSystem
            (
                M, globalCells, N,
                AopPim, AopPre,
                couplingLaplPre, couplingMassPre,
                couplingLaplPim, couplingMassPim
            );
            assembleBlockSystem
            (
                P, globalCells, N,
                AopPim, AopPre,
                couplingLaplPre, couplingMassShiftedPre,
                couplingLaplPim, couplingMassShiftedPim
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
            MatAssemblyBegin(P, MAT_FINAL_ASSEMBLY);
            MatAssemblyEnd(P, MAT_FINAL_ASSEMBLY);
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
    MatDestroy(&P);
    PetscFinalize();

    Info<< "End\n" << endl;
    return 0;
}

// ************************************************************************* //
