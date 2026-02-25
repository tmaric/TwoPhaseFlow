{
    #include "computeAlphaf.H"
    #include "computePMLCoefs.H"

    const dimensionedScalar k2
    (
        "k2",
        dimless/dimLength/dimLength,
        sqr(twoPi()*f/cg).value()
    );

    globalIndex globalCells(mesh.nCells());
    const PetscInt N = (PetscInt)globalCells.size();

    fvScalarMatrix AopPre
    (
        fvm::laplacian(1 - (rhol - rhog)/rhol*alphaf, Pre)
      + fvm::Sp(k2*(1 + ((kl - kg)/kg)*alpha1) - SC0, Pre)
    );

    fvScalarMatrix AopPim
    (
        fvm::laplacian(1 - (rhol - rhog)/rhol*alphaf, Pim)
      + fvm::Sp(k2*(1 + ((kl - kg)/kg)*alpha1) - SC0, Pim)
    );

    fvScalarMatrix couplingLaplPre(fvm::laplacian(TC1, Pre));
    fvScalarMatrix couplingMassPre(fvm::Sp(SC1, Pre));
    fvScalarMatrix couplingLaplPim(fvm::laplacian(TC1, Pim));
    fvScalarMatrix couplingMassPim(fvm::Sp(SC1, Pim));

    labelList localRows;
    labelList localCols;
    scalarField localVals;
    appendBlockTriplets
    (
        globalCells,
        N,
        AopPim,
        AopPre,
        couplingLaplPre,
        couplingMassPre,
        couplingLaplPim,
        couplingMassPim,
        localRows,
        localCols,
        localVals
    );

    // Keep processor-coupled contributions implicit in matrix entries.
    // Must be after boundaryManipulate (called in appendBlockTriplets) to
    // mirror interFALFlow algebra.
    scalarField bPim(sourceWithBoundary(AopPim, false));
    scalarField bCouplingLaplPre(sourceWithBoundary(couplingLaplPre, false));
    scalarField bCouplingMassPre(sourceWithBoundary(couplingMassPre, false));
    bPim -= bCouplingLaplPre;
    bPim -= bCouplingMassPre;

    scalarField bPre(sourceWithBoundary(AopPre, false));
    scalarField bCouplingLaplPim(sourceWithBoundary(couplingLaplPim, false));
    scalarField bCouplingMassPim(sourceWithBoundary(couplingMassPim, false));
    bPre += bCouplingLaplPim;
    bPre += bCouplingMassPim;

    labelList localGids(mesh.nCells());
    forAll(localGids, i)
    {
        localGids[i] = globalCells.toGlobal(i);
    }

    scalarField localPim(mesh.nCells(), Zero);
    scalarField localPre(mesh.nCells(), Zero);

    if (Pstream::parRun())
    {
        // Dedicated tags avoid collisions with other OpenFOAM exchanges.
        const int gatherTag = 71001;
        const int scatterTag = 71002;

        if (Pstream::master())
        {
            const label nProcs = Pstream::nProcs();
            List<labelList> allRows(nProcs), allCols(nProcs), allGids(nProcs);
            List<scalarField> allVals(nProcs), allBPim(nProcs), allBPre(nProcs);
            allRows[0].transfer(localRows);
            allCols[0].transfer(localCols);
            allVals[0].transfer(localVals);
            allGids[0].transfer(localGids);
            allBPim[0].transfer(bPim);
            allBPre[0].transfer(bPre);

            for (label proc = 1; proc < nProcs; ++proc)
            {
                IPstream fromProc(Pstream::commsTypes::blocking, proc, 0, gatherTag);
                fromProc
                    >> allRows[proc]
                    >> allCols[proc]
                    >> allVals[proc]
                    >> allGids[proc]
                    >> allBPim[proc]
                    >> allBPre[proc];
            }

            Mat M;
            Vec x, b;
            KSP ksp;
            MatCreate(PETSC_COMM_SELF, &M);
            MatSetSizes(M, 2*N, 2*N, 2*N, 2*N);
            MatSetType(M, MATSEQAIJ);
            MatSetUp(M);

            VecCreateSeq(PETSC_COMM_SELF, 2*N, &x);
            VecDuplicate(x, &b);

            for (label proc = 0; proc < nProcs; ++proc)
            {
                const labelList& r = allRows[proc];
                const labelList& c = allCols[proc];
                const scalarField& v = allVals[proc];
                forAll(r, i)
                {
                    MatSetValue(M, r[i], c[i], v[i], INSERT_VALUES);
                }
            }

            for (label proc = 0; proc < nProcs; ++proc)
            {
                const labelList& gids = allGids[proc];
                const scalarField& bp = allBPim[proc];
                const scalarField& br = allBPre[proc];
                forAll(gids, i)
                {
                    VecSetValue(b, gids[i], bp[i], INSERT_VALUES);
                    VecSetValue(b, gids[i] + N, br[i], INSERT_VALUES);
                }
            }

            MatAssemblyBegin(M, MAT_FINAL_ASSEMBLY);
            MatAssemblyEnd(M, MAT_FINAL_ASSEMBLY);
            VecAssemblyBegin(b);
            VecAssemblyEnd(b);

            KSPCreate(PETSC_COMM_SELF, &ksp);
            KSPSetOperators(ksp, M, M);
            PetscOptionsSetValue(nullptr, "-ksp_type", "preonly");
            PetscOptionsSetValue(nullptr, "-pc_type", "lu");
            PetscOptionsSetValue(nullptr, "-pc_factor_mat_solver_type", "mumps");
            KSPSetFromOptions(ksp);

            KSPSolve(ksp, b, x);

            scalarField globalPim(N, Zero), globalPre(N, Zero);
            const PetscScalar* xArr = nullptr;
            VecGetArrayRead(x, &xArr);
            for (PetscInt gi = 0; gi < N; ++gi)
            {
                globalPim[gi] = PetscRealPart(xArr[gi]);
                globalPre[gi] = PetscRealPart(xArr[gi + N]);
            }
            VecRestoreArrayRead(x, &xArr);

            for (label proc = 0; proc < nProcs; ++proc)
            {
                const labelList& gids = allGids[proc];
                scalarField sendPim(gids.size(), Zero);
                scalarField sendPre(gids.size(), Zero);
                forAll(gids, i)
                {
                    sendPim[i] = globalPim[gids[i]];
                    sendPre[i] = globalPre[gids[i]];
                }

                if (proc == 0)
                {
                    localPim.transfer(sendPim);
                    localPre.transfer(sendPre);
                }
                else
                {
                    OPstream toProc(Pstream::commsTypes::blocking, proc, 0, scatterTag);
                    toProc << sendPim << sendPre;
                }
            }

            KSPDestroy(&ksp);
            VecDestroy(&x);
            VecDestroy(&b);
            MatDestroy(&M);
        }
        else
        {
            {
                OPstream toMaster
                (
                    Pstream::commsTypes::blocking,
                    Pstream::masterNo(),
                    0,
                    gatherTag
                );
                toMaster
                    << localRows
                    << localCols
                    << localVals
                    << localGids
                    << bPim
                    << bPre;
            }

            IPstream fromMaster
            (
                Pstream::commsTypes::blocking,
                Pstream::masterNo(),
                0,
                scatterTag
            );
            fromMaster >> localPim >> localPre;
        }
    }
    else
    {
        Mat M;
        Vec x, b;
        KSP ksp;

        MatCreate(PETSC_COMM_SELF, &M);
        MatSetSizes(M, 2*N, 2*N, 2*N, 2*N);
        MatSetType(M, MATSEQAIJ);
        MatSetUp(M);

        VecCreateSeq(PETSC_COMM_SELF, 2*N, &x);
        VecDuplicate(x, &b);

        forAll(localRows, i)
        {
            MatSetValue(M, localRows[i], localCols[i], localVals[i], INSERT_VALUES);
        }
        forAll(localGids, i)
        {
            VecSetValue(b, localGids[i], bPim[i], INSERT_VALUES);
            VecSetValue(b, localGids[i] + N, bPre[i], INSERT_VALUES);
        }

        MatAssemblyBegin(M, MAT_FINAL_ASSEMBLY);
        MatAssemblyEnd(M, MAT_FINAL_ASSEMBLY);
        VecAssemblyBegin(b);
        VecAssemblyEnd(b);

        KSPCreate(PETSC_COMM_SELF, &ksp);
        KSPSetOperators(ksp, M, M);
        PetscOptionsSetValue(nullptr, "-ksp_type", "preonly");
        PetscOptionsSetValue(nullptr, "-pc_type", "lu");
        PetscOptionsSetValue(nullptr, "-pc_factor_mat_solver_type", "mumps");
        KSPSetFromOptions(ksp);

        KSPSolve(ksp, b, x);

        const PetscScalar* xArr = nullptr;
        VecGetArrayRead(x, &xArr);
        forAll(localGids, i)
        {
            const label gi = localGids[i];
            localPim[i] = PetscRealPart(xArr[gi]);
            localPre[i] = PetscRealPart(xArr[gi + N]);
        }
        VecRestoreArrayRead(x, &xArr);

        KSPDestroy(&ksp);
        VecDestroy(&x);
        VecDestroy(&b);
        MatDestroy(&M);
    }

    forAll(Pim, i)
    {
        Pim[i] = localPim[i];
        Pre[i] = localPre[i];
    }

    Pre.correctBoundaryConditions();
    Pim.correctBoundaryConditions();

    const tmp<volVectorField> Ure = 1/(twoPi()*f*rho) * fvc::grad(Pim);
    const tmp<volVectorField> Uim = -1/(twoPi()*f*rho) * fvc::grad(Pre);

    pa == Foam::sqrt(Pim*Pim + Pre*Pre);
    pr == 0.25*(kl*alpha1 + kg*(1-alpha1))*(Pre*Pre + Pim*Pim)
        - 0.25*rho*((Ure()&Ure()) + (Uim()&Uim()));
    momFlux == 0.5*rho*(Ure()*Ure() + Uim()*Uim());
}
