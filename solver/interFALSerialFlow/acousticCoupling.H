{
    // Acoustic coupling sequence for each outer flow iteration:
    // 1) Recompute alpha-face transport and PML coefficients.
    // 2) Assemble real-valued 2x2 block system for complex pressure:
    //      [ A  -(B1 + B2) ] [Pim] = [bPim]
    //      [ (B1 + B2)  A ] [Pre]   [bPre]
    // 3) Solve once with PETSc and update Pre/Pim.
    // 4) Reconstruct pr/momFlux used by interfacial flow forcing.
    #include "computeAlphaf.H"

    const dimensionedScalar k2
    (
        "k2",
        dimless/dimLength/dimLength,
        sqr(twoPi()*f/cg).value()
    );

    globalIndex globalCells(mesh.nCells());
    const PetscInt N = (PetscInt)globalCells.size();
    const PetscInt nLocal = (PetscInt)mesh.nCells();

    Mat M;
    Vec x, b;
    KSP ksp;
    initializePetscSystem(nLocal, N, M, x, b, ksp);

    MatZeroEntries(M);
    VecSet(b, 0.0);

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

    assembleBlockSystem
    (
        M, globalCells, N,
        AopPim, AopPre,
        couplingLaplPre, couplingMassPre,
        couplingLaplPim, couplingMassPim
    );

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

    KSPSolve(ksp, b, x);
    scatterBlockSolution(x, globalCells, N, Pim, Pre);

    Pre.correctBoundaryConditions();
    Pim.correctBoundaryConditions();

    const tmp<volVectorField> Ure = 1/(twoPi()*f*rho) * fvc::grad(Pim);
    const tmp<volVectorField> Uim = -1/(twoPi()*f*rho) * fvc::grad(Pre);

    pr == 0.25*(kl*alpha1 + kg*(1-alpha1))*(Pre*Pre + Pim*Pim)
        - 0.25*rho*((Ure()&Ure()) + (Uim()&Uim()));
    momFlux == 0.5*rho*(Ure()*Ure() + Uim()*Uim());

    KSPDestroy(&ksp);
    VecDestroy(&x);
    VecDestroy(&b);
    MatDestroy(&M);
}
