{
    // Acoustic coupling sequence for each outer flow iteration:
    // 1) Recompute alpha-face transport and PML coefficients from alpha1.
    // 2) Assemble distributed 2x2 block Helmholtz system:
    //      [ A  -(B1 + B2) ] [Pim] = [bPim]
    //      [ (B1 + B2)  A ] [Pre]   [bPre]
    // 3) Solve with PETSc in distributed mode and scatter solution fields.
    // 4) Recompute pa/pr/momFlux to feed the flow-side forcing terms.
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

    fvScalarMatrix mainOpPre
    (
        fvm::laplacian(1 - (rhol - rhog)/rhol*alphaf, Pre)
      + fvm::Sp(k2*(1 + ((kl - kg)/kg)*alpha1) - SC0, Pre)
    );

    fvScalarMatrix mainOpPim
    (
        fvm::laplacian(1 - (rhol - rhog)/rhol*alphaf, Pim)
      + fvm::Sp(k2*(1 + ((kl - kg)/kg)*alpha1) - SC0, Pim)
    );

    // Cross-coupling operators are built from the opposite field.
    fvScalarMatrix couplingLaplPre(fvm::laplacian(TC1, Pre));
    fvScalarMatrix couplingMassPre(fvm::Sp(SC1, Pre));
    fvScalarMatrix couplingOpPre(couplingLaplPre + couplingMassPre);

    fvScalarMatrix couplingLaplPim(fvm::laplacian(TC1, Pim));
    fvScalarMatrix couplingMassPim(fvm::Sp(SC1, Pim));
    fvScalarMatrix couplingOpPim(couplingLaplPim + couplingMassPim);

    const PetscInt nLocal = (PetscInt)mesh.nCells();
    Mat M;
    Vec x, b;
    KSP ksp;
    initializeAcousticPetscSystem(nLocal, N, M, x, b, ksp);

    assembleBlockSystem
    (
        M,
        globalCells,
        N,
        mainOpPim,
        mainOpPre,
        couplingOpPre
    );

    scalarField bPim;
    scalarField bPre;
    buildAcousticRhs
    (
        mainOpPim,
        mainOpPre,
        couplingOpPre,
        couplingOpPim,
        bPim,
        bPre
    );

    setBlockRhs(b, globalCells, N, bPim, bPre);

    MatAssemblyBegin(M, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(M, MAT_FINAL_ASSEMBLY);
    VecAssemblyBegin(b);
    VecAssemblyEnd(b);

    KSPSolve(ksp, b, x);

    scatterBlockSolution(x, globalCells, N, Pim, Pre);

    KSPDestroy(&ksp);
    VecDestroy(&x);
    VecDestroy(&b);
    MatDestroy(&M);

    Pre.correctBoundaryConditions();
    Pim.correctBoundaryConditions();

    updateDerivedAcousticFields
    (
        pa,
        pr,
        momFlux,
        Pim,
        Pre,
        alpha1,
        rho,
        f,
        kl,
        kg
    );
}
