// PETSc lifecycle/solve helpers for freITBCFoam.

#include <vector>

struct SweepPCContext
{
    Mat M;
    Mat P;
    PetscInt nGlobal;
    std::vector<PetscInt> layerOf;
    std::vector<std::vector<PetscInt>> layers;
    std::vector<IS> layerIs;
    std::vector<Mat> layerMat;
    std::vector<KSP> layerKsp;
    std::vector<Vec> layerB;
    std::vector<Vec> layerZ;
    PetscInt nSweeps;
    PetscBool symmetric;
    PetscReal omega;
    PetscBool postDiagJacobi;
    PetscReal postDiagOmega;
    PetscBool postGlobalIlu;
    PetscReal postGlobalOmega;
    PetscInt postGlobalSteps;
    PetscBool postDefectCorrection;
    PetscReal defectOmega;
    PetscInt defectSteps;
    KSP globalKsp;
    Vec globalR;
    Vec globalZ;
    Vec diag;
    Vec res;
    Vec py;
};

static PetscErrorCode destroySweepPC(PC pc)
{
    SweepPCContext* ctx = nullptr;
    PCShellGetContext(pc, (void**)&ctx);
    if (!ctx)
    {
        return PETSC_SUCCESS;
    }

    for (std::size_t li = 0; li < ctx->layerKsp.size(); ++li)
    {
        if (ctx->layerKsp[li]) KSPDestroy(&ctx->layerKsp[li]);
        if (ctx->layerMat[li]) MatDestroy(&ctx->layerMat[li]);
        if (ctx->layerIs[li]) ISDestroy(&ctx->layerIs[li]);
        if (ctx->layerB[li]) VecDestroy(&ctx->layerB[li]);
        if (ctx->layerZ[li]) VecDestroy(&ctx->layerZ[li]);
    }
    if (ctx->diag) VecDestroy(&ctx->diag);
    if (ctx->res) VecDestroy(&ctx->res);
    if (ctx->py) VecDestroy(&ctx->py);
    if (ctx->globalR) VecDestroy(&ctx->globalR);
    if (ctx->globalZ) VecDestroy(&ctx->globalZ);
    if (ctx->globalKsp) KSPDestroy(&ctx->globalKsp);

    delete ctx;
    PCShellSetContext(pc, nullptr);
    return PETSC_SUCCESS;
}

static void solveLayerSweepPass
(
    SweepPCContext* ctx,
    const PetscScalar* xArr,
    PetscScalar* yArr,
    const bool forward
)
{
    const PetscInt nLayers = (PetscInt)ctx->layers.size();

    for (PetscInt passi = 0; passi < nLayers; ++passi)
    {
        const PetscInt li = forward ? passi : (nLayers - 1 - passi);

        MatCreateSubMatrix
        (
            ctx->P,
            ctx->layerIs[li],
            ctx->layerIs[li],
            ctx->layerMat[li] ? MAT_REUSE_MATRIX : MAT_INITIAL_MATRIX,
            &ctx->layerMat[li]
        );

        if (!ctx->layerKsp[li])
        {
            KSPCreate(PETSC_COMM_SELF, &ctx->layerKsp[li]);
            KSPSetType(ctx->layerKsp[li], KSPPREONLY);
            PC subPc = nullptr;
            KSPGetPC(ctx->layerKsp[li], &subPc);
            PCSetType(subPc, PCILU);
            PCFactorSetLevels(subPc, 2);
            KSPSetOptionsPrefix(ctx->layerKsp[li], "itbc_sub_");
            KSPSetFromOptions(ctx->layerKsp[li]);
        }

        KSPSetOperators(ctx->layerKsp[li], ctx->layerMat[li], ctx->layerMat[li]);

        PetscScalar* bArr = nullptr;
        VecGetArray(ctx->layerB[li], &bArr);

        const std::vector<PetscInt>& rows = ctx->layers[li];
        for (PetscInt locali = 0; locali < (PetscInt)rows.size(); ++locali)
        {
            const PetscInt row = rows[locali];

            const PetscInt* cols = nullptr;
            const PetscScalar* vals = nullptr;
            PetscInt nCols = 0;
            MatGetRow(ctx->P, row, &nCols, &cols, &vals);

            PetscScalar rhs = xArr[row];
            for (PetscInt k = 0; k < nCols; ++k)
            {
                const PetscInt col = cols[k];
                const PetscScalar a = vals[k];
                const PetscInt lCol = ctx->layerOf[col];

                if (lCol == li)
                {
                    continue;
                }

                if (forward)
                {
                    if (lCol < li) rhs -= a*yArr[col];
                }
                else
                {
                    if (lCol > li) rhs -= a*yArr[col];
                }
            }

            bArr[locali] = rhs;
            MatRestoreRow(ctx->P, row, &nCols, &cols, &vals);
        }

        VecRestoreArray(ctx->layerB[li], &bArr);

        KSPSolve(ctx->layerKsp[li], ctx->layerB[li], ctx->layerZ[li]);

        const PetscScalar* zArr = nullptr;
        VecGetArrayRead(ctx->layerZ[li], &zArr);
        for (PetscInt locali = 0; locali < (PetscInt)rows.size(); ++locali)
        {
            const PetscInt row = rows[locali];
            yArr[row] = (1.0 - ctx->omega)*yArr[row] + ctx->omega*zArr[locali];
        }
        VecRestoreArrayRead(ctx->layerZ[li], &zArr);
    }
}

static PetscErrorCode applySweepPC(PC pc, Vec x, Vec y)
{
    SweepPCContext* ctx = nullptr;
    PCShellGetContext(pc, (void**)&ctx);

    const PetscScalar* xArr = nullptr;
    PetscScalar* yArr = nullptr;
    VecGetArrayRead(x, &xArr);
    VecGetArray(y, &yArr);
    VecSet(y, 0.0);

    for (PetscInt sweep = 0; sweep < ctx->nSweeps; ++sweep)
    {
        solveLayerSweepPass(ctx, xArr, yArr, true);

        if (ctx->symmetric)
        {
            solveLayerSweepPass(ctx, xArr, yArr, false);
        }
    }

    if (ctx->postDiagJacobi)
    {
        if (!ctx->diag)
        {
            VecDuplicate(x, &ctx->diag);
            VecDuplicate(x, &ctx->res);
            VecDuplicate(x, &ctx->py);
        }

        MatGetDiagonal(ctx->P, ctx->diag);
        MatMult(ctx->P, y, ctx->py);
        VecWAXPY(ctx->res, -1.0, ctx->py, x); // res = x - P*y

        const PetscScalar* dArr = nullptr;
        const PetscScalar* rArr = nullptr;
        VecGetArrayRead(ctx->diag, &dArr);
        VecGetArrayRead(ctx->res, &rArr);
        for (PetscInt i = 0; i < ctx->nGlobal; ++i)
        {
            if (PetscAbsScalar(dArr[i]) > 1e-30)
            {
                yArr[i] += ctx->postDiagOmega * (rArr[i]/dArr[i]);
            }
        }
        VecRestoreArrayRead(ctx->res, &rArr);
        VecRestoreArrayRead(ctx->diag, &dArr);
    }

    VecRestoreArray(y, &yArr);
    VecRestoreArrayRead(x, &xArr);

    if (ctx->postGlobalIlu || ctx->postDefectCorrection)
    {
        if (!ctx->globalR)
        {
            VecDuplicate(x, &ctx->globalR);
            VecDuplicate(x, &ctx->globalZ);
            if (!ctx->py) VecDuplicate(x, &ctx->py);
        }
        if (!ctx->globalKsp)
        {
            KSPCreate(PETSC_COMM_WORLD, &ctx->globalKsp);
            KSPSetType(ctx->globalKsp, KSPPREONLY);
            PC gpc = nullptr;
            KSPGetPC(ctx->globalKsp, &gpc);
            PCSetType(gpc, PCILU);
            PCFactorSetLevels(gpc, 1);
            KSPSetOptionsPrefix(ctx->globalKsp, "itbc_global_");
            KSPSetFromOptions(ctx->globalKsp);
        }
        KSPSetOperators(ctx->globalKsp, ctx->P, ctx->P);

        // Optional classical correction based on preconditioner residual.
        for (PetscInt step = 0; step < ctx->postGlobalSteps; ++step)
        {
            MatMult(ctx->P, y, ctx->py);
            VecWAXPY(ctx->globalR, -1.0, ctx->py, x); // r = x - P*y
            KSPSolve(ctx->globalKsp, ctx->globalR, ctx->globalZ);
            VecAXPY(y, ctx->postGlobalOmega, ctx->globalZ);
        }

        // Defect-correction towards M^{-1} using P^{-1} as inner approximate solve.
        for (PetscInt step = 0; step < ctx->defectSteps && ctx->postDefectCorrection; ++step)
        {
            MatMult(ctx->M, y, ctx->py);
            VecWAXPY(ctx->globalR, -1.0, ctx->py, x); // r = x - M*y
            KSPSolve(ctx->globalKsp, ctx->globalR, ctx->globalZ);
            VecAXPY(y, ctx->defectOmega, ctx->globalZ);
        }
    }

    return PETSC_SUCCESS;
}

static void initializePetscSystem
(
    const PetscInt nLocal,
    const PetscInt N,
    const List<PetscInt>& sweepOrder,
    const List<List<PetscInt>>& sweepLayers,
    Mat& M,
    Mat& P,
    Vec& x,
    Vec& b,
    KSP& ksp
)
{
    MatCreate(PETSC_COMM_WORLD, &M);
    MatSetSizes(M, 2*nLocal, 2*nLocal, 2*N, 2*N);
    MatSetType(M, (Pstream::nProcs() == 1) ? MATSEQAIJ : MATMPIAIJ);
    MatSetUp(M);

    MatCreate(PETSC_COMM_WORLD, &P);
    MatSetSizes(P, 2*nLocal, 2*nLocal, 2*N, 2*N);
    MatSetType(P, (Pstream::nProcs() == 1) ? MATSEQAIJ : MATMPIAIJ);
    MatSetUp(P);

    VecCreate(PETSC_COMM_WORLD, &x);
    VecSetSizes(x, 2*nLocal, 2*N);
    VecSetFromOptions(x);
    VecDuplicate(x, &b);

    KSPCreate(PETSC_COMM_WORLD, &ksp);
    KSPSetOperators(ksp, M, P);

    PetscBool hasOpt = PETSC_FALSE;

    PetscOptionsHasName(nullptr, nullptr, "-ksp_type", &hasOpt);
    if (!hasOpt) PetscOptionsSetValue(nullptr, "-ksp_type", "fgmres");

    PetscOptionsHasName(nullptr, nullptr, "-ksp_gmres_restart", &hasOpt);
    if (!hasOpt) PetscOptionsSetValue(nullptr, "-ksp_gmres_restart", "200");

    PetscOptionsHasName(nullptr, nullptr, "-ksp_rtol", &hasOpt);
    if (!hasOpt) PetscOptionsSetValue(nullptr, "-ksp_rtol", "1e-8");

    PetscOptionsHasName(nullptr, nullptr, "-ksp_atol", &hasOpt);
    if (!hasOpt) PetscOptionsSetValue(nullptr, "-ksp_atol", "1e-14");

    PetscOptionsHasName(nullptr, nullptr, "-ksp_max_it", &hasOpt);
    if (!hasOpt) PetscOptionsSetValue(nullptr, "-ksp_max_it", "2000");

    PetscOptionsHasName(nullptr, nullptr, "-pc_type", &hasOpt);
    if (!hasOpt) PetscOptionsSetValue(nullptr, "-pc_type", "shell");

    KSPSetFromOptions(ksp);

    PC pc = nullptr;
    KSPGetPC(ksp, &pc);
    PCType pcType = nullptr;
    PCGetType(pc, &pcType);
    PetscBool isShell = PETSC_FALSE;
    if (pcType)
    {
        PetscStrcmp(pcType, PCSHELL, &isShell);
    }

    if (isShell)
    {
        SweepPCContext* ctx = new SweepPCContext();
        ctx->M = M;
        ctx->P = P;
        ctx->nGlobal = 2*N;
        ctx->layerOf.assign(2*N, -1);
        ctx->nSweeps = 1;
        ctx->symmetric = PETSC_FALSE;
        ctx->omega = 1.0;
        ctx->postDiagJacobi = PETSC_TRUE;
        ctx->postDiagOmega = 0.8;
        ctx->postGlobalIlu = PETSC_TRUE;
        ctx->postGlobalOmega = 0.3;
        ctx->postGlobalSteps = 1;
        ctx->postDefectCorrection = PETSC_FALSE;
        ctx->defectOmega = 0.7;
        ctx->defectSteps = 1;
        ctx->globalKsp = nullptr;
        ctx->globalR = nullptr;
        ctx->globalZ = nullptr;
        ctx->diag = nullptr;
        ctx->res = nullptr;
        ctx->py = nullptr;

        PetscOptionsGetInt(nullptr, nullptr, "-itbc_sweeps", &ctx->nSweeps, nullptr);
        PetscOptionsGetBool(nullptr, nullptr, "-itbc_symmetric", &ctx->symmetric, nullptr);
        PetscOptionsGetReal(nullptr, nullptr, "-itbc_omega", &ctx->omega, nullptr);
        PetscOptionsGetBool(nullptr, nullptr, "-itbc_post_diag_jacobi", &ctx->postDiagJacobi, nullptr);
        PetscOptionsGetReal(nullptr, nullptr, "-itbc_post_diag_omega", &ctx->postDiagOmega, nullptr);
        PetscOptionsGetBool(nullptr, nullptr, "-itbc_post_global_ilu", &ctx->postGlobalIlu, nullptr);
        PetscOptionsGetReal(nullptr, nullptr, "-itbc_post_global_omega", &ctx->postGlobalOmega, nullptr);
        PetscOptionsGetInt(nullptr, nullptr, "-itbc_post_global_steps", &ctx->postGlobalSteps, nullptr);
        PetscOptionsGetBool(nullptr, nullptr, "-itbc_post_defect_correction", &ctx->postDefectCorrection, nullptr);
        PetscOptionsGetReal(nullptr, nullptr, "-itbc_defect_omega", &ctx->defectOmega, nullptr);
        PetscOptionsGetInt(nullptr, nullptr, "-itbc_defect_steps", &ctx->defectSteps, nullptr);

        const label nLayers = sweepLayers.size();
        ctx->layers.resize(nLayers);
        ctx->layerIs.resize(nLayers, nullptr);
        ctx->layerMat.resize(nLayers, nullptr);
        ctx->layerKsp.resize(nLayers, nullptr);
        ctx->layerB.resize(nLayers, nullptr);
        ctx->layerZ.resize(nLayers, nullptr);

        for (label li = 0; li < nLayers; ++li)
        {
            const List<PetscInt>& layerRows = sweepLayers[li];
            ctx->layers[li].assign(layerRows.begin(), layerRows.end());

            for (const PetscInt row : ctx->layers[li])
            {
                ctx->layerOf[row] = li;
            }

            ISCreateGeneral
            (
                PETSC_COMM_SELF,
                (PetscInt)ctx->layers[li].size(),
                ctx->layers[li].data(),
                PETSC_COPY_VALUES,
                &ctx->layerIs[li]
            );

            VecCreateSeq(PETSC_COMM_SELF, (PetscInt)ctx->layers[li].size(), &ctx->layerB[li]);
            VecDuplicate(ctx->layerB[li], &ctx->layerZ[li]);
        }

        PCShellSetContext(pc, (void*)ctx);
        PCShellSetApply(pc, applySweepPC);
        PCShellSetDestroy(pc, destroySweepPC);
        PCShellSetName(pc, "yLayerSweepBlock");

        if (ctx->layers.empty())
        {
            FatalErrorInFunction
                << "Sweeping preconditioner has zero layers."
                << exit(FatalError);
        }
    }
}

static void scatterBlockSolution
(
    Vec& x,
    const globalIndex& globalCells,
    const PetscInt N,
    volScalarField& Pim,
    volScalarField& Pre
)
{
    const PetscScalar* xArr = nullptr;
    VecGetArrayRead(x, &xArr);
    forAll(Pim, i)
    {
        const PetscInt gi = globalCells.toGlobal(i);
        Pim[i] = (scalar)PetscRealPart(xArr[gi + 0]);
        Pre[i] = (scalar)PetscRealPart(xArr[gi + N]);
    }
    VecRestoreArrayRead(x, &xArr);
}
