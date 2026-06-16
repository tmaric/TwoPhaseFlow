{
    // Acoustic coupling sequence for each outer flow iteration:
    // 1) Recompute alpha-face transport and formulation-specific PML
    //    coefficients from the sigma field prepared by setPMLFields.
    // 2) Assemble distributed 2x2 block Helmholtz system:
    //      [ A  -(B1 + B2) ] [Pim] = [bPim]
    //      [ (B1 + B2)  A ] [Pre]   [bPre]
    // 3) Solve with PETSc in distributed mode and scatter solution fields.
    // 4) Recompute pr/momFlux to feed the flow-side forcing terms.
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
    const PetscInt nLocal = (PetscInt)mesh.nCells();

    const label nAcousticNonOrthogonalCorrectors =
        mesh.solutionDict().subDict("PIMPLE").lookupOrDefault<label>
        (
            "nAcousticNonOrthogonalCorrectors",
            0
        );

    for
    (
        label acousticNonOrthCorr = 0;
        acousticNonOrthCorr <= nAcousticNonOrthogonalCorrectors;
        ++acousticNonOrthCorr
    )
    {
        fvScalarMatrix AopPre
        (
            fvm::laplacian(1 - (rhol - rhog)/rhol*alphaf, Pre)
          + fvm::laplacian(T0, Pre)
          + fvm::Sp(k2*(1 + ((kl - kg)/kg)*alpha1) - SC0, Pre)
        );

        fvScalarMatrix AopPim
        (
            fvm::laplacian(1 - (rhol - rhog)/rhol*alphaf, Pim)
          + fvm::laplacian(T0, Pim)
          + fvm::Sp(k2*(1 + ((kl - kg)/kg)*alpha1) - SC0, Pim)
        );

        fvScalarMatrix B1Pre(fvm::laplacian(T1, Pre));
        fvScalarMatrix B2Pre(fvm::Sp(SC1, Pre));
        fvScalarMatrix B1Pim(fvm::laplacian(T1, Pim));
        fvScalarMatrix B2Pim(fvm::Sp(SC1, Pim));

        const PetscInt rowNnzUpper = estimateBlockPreallocationUpper(AopPim);

        ensurePetscSystem(nLocal, N, rowNnzUpper, acousticSystem);
        resetPetscSystem(acousticSystem);

        assembleBlockSystem
        (
            acousticSystem.M,
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

        setBlockRhs(acousticSystem.b, globalCells, N, bPim, bPre);

        MatAssemblyBegin(acousticSystem.M, MAT_FINAL_ASSEMBLY);
        MatAssemblyEnd(acousticSystem.M, MAT_FINAL_ASSEMBLY);
        VecAssemblyBegin(acousticSystem.b);
        VecAssemblyEnd(acousticSystem.b);

        KSPSolve(acousticSystem.ksp, acousticSystem.b, acousticSystem.x);

        scatterBlockSolution(acousticSystem, globalCells, Pim, Pre);

        Pre.correctBoundaryConditions();
        Pim.correctBoundaryConditions();
    }

    updateDerivedAcousticFields
    (
        pr, Pa, Ua, momFlux, Pim, Pre, alpha1, rho, f, kl, kg
    );
}
