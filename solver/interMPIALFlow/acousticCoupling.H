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

    // Coupling operators built from opposite field.
    fvScalarMatrix couplingLaplPre(fvm::laplacian(TC1, Pre));
    fvScalarMatrix couplingMassPre(fvm::Sp(SC1, Pre));
    fvScalarMatrix couplingPre(couplingLaplPre + couplingMassPre);

    fvScalarMatrix couplingLaplPim(fvm::laplacian(TC1, Pim));
    fvScalarMatrix couplingMassPim(fvm::Sp(SC1, Pim));
    fvScalarMatrix couplingPim(couplingLaplPim + couplingMassPim);

    const PetscInt nLocal = (PetscInt)mesh.nCells();
    Mat M;
    Vec x, b;
    KSP ksp;

    MatCreate(PETSC_COMM_WORLD, &M);
    MatSetSizes(M, 2*nLocal, 2*nLocal, 2*N, 2*N);
    MatSetType(M, MATAIJ);
    MatSetUp(M);

    VecCreate(PETSC_COMM_WORLD, &x);
    VecSetSizes(x, 2*nLocal, 2*N);
    VecSetFromOptions(x);
    VecDuplicate(x, &b);

    assembleBlockSystem
    (
        M,
        globalCells,
        N,
        AopPim,
        AopPre,
        couplingPre
    );

    // Keep processor-coupled interfaces implicit in matrix so decomposed
    // algebra follows the freBCMFoam path.
    scalarField bPim(sourceWithBoundary(AopPim, false));
    scalarField bCouplingPre(sourceWithBoundary(couplingPre, false));
    bPim -= bCouplingPre;

    scalarField bPre(sourceWithBoundary(AopPre, false));
    scalarField bCouplingPim(sourceWithBoundary(couplingPim, false));
    bPre += bCouplingPim;

    setBlockRhs(b, globalCells, N, bPim, bPre);

    MatAssemblyBegin(M, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(M, MAT_FINAL_ASSEMBLY);
    VecAssemblyBegin(b);
    VecAssemblyEnd(b);

    KSPCreate(PETSC_COMM_WORLD, &ksp);
    KSPSetOperators(ksp, M, M);
    PetscOptionsSetValue(nullptr, "-ksp_type", "preonly");
    PetscOptionsSetValue(nullptr, "-pc_type", "lu");
    PetscOptionsSetValue(nullptr, "-pc_factor_mat_solver_type", "mumps");
    KSPSetFromOptions(ksp);

    KSPSolve(ksp, b, x);

    scatterBlockSolution(x, globalCells, N, Pim, Pre);

    KSPDestroy(&ksp);
    VecDestroy(&x);
    VecDestroy(&b);
    MatDestroy(&M);

    Pre.correctBoundaryConditions();
    Pim.correctBoundaryConditions();

    const tmp<volVectorField> Ure = 1/(twoPi()*f*rho) * fvc::grad(Pim);
    const tmp<volVectorField> Uim = -1/(twoPi()*f*rho) * fvc::grad(Pre);

    pa == Foam::sqrt(Pim*Pim + Pre*Pre);
    pr == 0.25*(kl*alpha1 + kg*(1-alpha1))*(Pre*Pre + Pim*Pim)
        - 0.25*rho*((Ure()&Ure()) + (Uim()&Uim()));
    momFlux == 0.5*rho*(Ure()*Ure() + Uim()*Uim());
}
