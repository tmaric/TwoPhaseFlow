{
    // Acoustic coupling sequence for each outer flow iteration:
    // 1) Recompute alpha-face transport. PML coefficient fields are prepared
    //    beforehand by setPMLFields and read as normal case data.
    // 2) Assemble distributed 2x2 block Helmholtz system:
    //      [ A  -(B1 + B2) ] [Pim] = [bPim]
    //      [ (B1 + B2)  A ] [Pre]   [bPre]
    // 3) Solve with PETSc in distributed mode and scatter solution fields.
    // 4) Recompute pr/momFlux to feed the flow-side forcing terms.
    #include "computeAlphaf.H"

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

    fvScalarMatrix B1Pre(fvm::laplacian(TC1, Pre));
    fvScalarMatrix B2Pre(fvm::Sp(SC1, Pre));
    fvScalarMatrix B1Pim(fvm::laplacian(TC1, Pim));
    fvScalarMatrix B2Pim(fvm::Sp(SC1, Pim));

    const PetscInt nLocal = (PetscInt)mesh.nCells();
    Mat M;
    Vec x, b;
    KSP ksp;
    initializePetscSystem(nLocal, N, M, x, b, ksp);

    assembleBlockSystem
    (
        M,
        globalCells,
        N,
        AopPim,
        AopPre,
        B1Pre,
        B2Pre,
        B1Pim,
        B2Pim
    );

    scalarField bPim;
    scalarField bPre;
    buildRhs
    (
        AopPim,
        AopPre,
        B1Pre,
        B2Pre,
        B1Pim,
        B2Pim,
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
        pr, momFlux, Pim, Pre, alpha1, rho, f, kl, kg
    );

    updateLessmannAcousticForce
    (
        mesh,
        alpha1,
        pr,
        momFlux,
        lessmannShiftCoeff,
        lessmannInterfaceAlphaMin,
        lessmannInterfaceAlphaMax,
        lessmannGasAlphaMax,
        lessmannMinGradAlpha,
        prShift,
        acousticTraction,
        acousticForce
    );
}
