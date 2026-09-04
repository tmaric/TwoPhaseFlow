{
    const scalar deltaT = runTime.deltaTValue();
    const scalar theta = timeIntegrationTheta;
    const scalar oldWeight = 1.0 - theta;
    const scalarField& cellVolumes = mesh.V();

    if
    (
        timeBlockMatrixAssembled
     && mag(deltaT - timeBlockMatrixDeltaT)
      > SMALL*max(deltaT, timeBlockMatrixDeltaT)
    )
    {
        timeBlockMatrixAssembled = false;
    }

    p.correctBoundaryConditions();
    v.correctBoundaryConditions();
    u.correctBoundaryConditions();

    pmlUx == u.component(vector::X);
    pmlUy == u.component(vector::Y);
    pmlUz == u.component(vector::Z);

    fvScalarMatrix pressureSpatial
    (
        -rho*fvm::laplacian(invRhof, p, "laplacian(invRhof,p)")
      + fvm::Sp(C1, p)
    );

    fvScalarMatrix divUx(fvm::div(SfX, pmlUx, "div(u)"));
    fvScalarMatrix divUy(fvm::div(SfY, pmlUy, "div(u)"));
    fvScalarMatrix divUz(fvm::div(SfZ, pmlUz, "div(u)"));

    fvScalarMatrix gradPx(fvm::div(SfX, p, "div(u)"));
    fvScalarMatrix gradPy(fvm::div(SfY, p, "div(u)"));
    fvScalarMatrix gradPz(fvm::div(SfZ, p, "div(u)"));

    fvScalarMatrix gradVx(fvm::div(SfX, v, "div(u)"));
    fvScalarMatrix gradVy(fvm::div(SfY, v, "div(u)"));
    fvScalarMatrix gradVz(fvm::div(SfZ, v, "div(u)"));

    pressureSpatial.boundaryManipulate(p.boundaryFieldRef());
    divUx.boundaryManipulate(pmlUx.boundaryFieldRef());
    divUy.boundaryManipulate(pmlUy.boundaryFieldRef());
    divUz.boundaryManipulate(pmlUz.boundaryFieldRef());
    gradPx.boundaryManipulate(p.boundaryFieldRef());
    gradPy.boundaryManipulate(p.boundaryFieldRef());
    gradPz.boundaryManipulate(p.boundaryFieldRef());
    gradVx.boundaryManipulate(v.boundaryFieldRef());
    gradVy.boundaryManipulate(v.boundaryFieldRef());
    gradVz.boundaryManipulate(v.boundaryFieldRef());

    scalarField qDivScaleX(mesh.nCells());
    scalarField qDivScaleY(mesh.nCells());
    scalarField qDivScaleZ(mesh.nCells());
    scalarField uxGradPScale(mesh.nCells());
    scalarField uyGradPScale(mesh.nCells());
    scalarField uzGradPScale(mesh.nCells());
    scalarField uxGradVScale(mesh.nCells());
    scalarField uyGradVScale(mesh.nCells());
    scalarField uzGradVScale(mesh.nCells());

    forAll(qDivScaleX, celli)
    {
        qDivScaleX[celli] = -theta;
        qDivScaleY[celli] = -theta;
        qDivScaleZ[celli] = -theta;

        uxGradPScale[celli] = theta*T1[celli].xx();
        uyGradPScale[celli] = theta*T1[celli].yy();
        uzGradPScale[celli] = theta*T1[celli].zz();

        uxGradVScale[celli] = -theta*T2[celli].xx();
        uyGradVScale[celli] = -theta*T2[celli].yy();
        uzGradVScale[celli] = -theta*T2[celli].zz();
    }

    if (!timeBlockMatrixAssembled)
    {
        scalarField pP(mesh.nCells());
        scalarField pQ(mesh.nCells());
        scalarField qQ(mesh.nCells());
        scalarField qV(mesh.nCells());
        scalarField vP(mesh.nCells());
        scalarField vV(mesh.nCells());
        scalarField uxUx(mesh.nCells());
        scalarField uyUy(mesh.nCells());
        scalarField uzUz(mesh.nCells());

        forAll(pP, celli)
        {
            const scalar volume = cellVolumes[celli];
            pP[celli] = volume/deltaT;
            pQ[celli] = -theta*volume;

            qQ[celli] = volume
               *(rho[celli]*compressibility[celli]/deltaT
               + theta*C0[celli]);
            qV[celli] = theta*volume*C2[celli];

            vP[celli] = -theta*volume;
            vV[celli] = volume/deltaT;

            uxUx[celli] = volume
               *(1.0/deltaT + theta*TC0[celli].xx());
            uyUy[celli] = volume
               *(1.0/deltaT + theta*TC0[celli].yy());
            uzUz[celli] = volume
               *(1.0/deltaT + theta*TC0[celli].zz());
        }

        MatZeroEntries(timeBlockMatrix);

        addTimeBlockCellDiagonal
        (
            timeBlockMatrix,
            globalCells,
            mesh,
            pComponent,
            pComponent,
            pP
        );
        addTimeBlockCellDiagonal
        (
            timeBlockMatrix,
            globalCells,
            mesh,
            pComponent,
            qComponent,
            pQ
        );

        addTimeBlockCellDiagonal
        (
            timeBlockMatrix,
            globalCells,
            mesh,
            qComponent,
            qComponent,
            qQ
        );
        insertTimeBlockOperator
        (
            timeBlockMatrix,
            pressureSpatial,
            globalCells,
            qComponent,
            pComponent,
            theta
        );
        addTimeBlockCellDiagonal
        (
            timeBlockMatrix,
            globalCells,
            mesh,
            qComponent,
            vComponent,
            qV
        );
        insertTimeBlockOperator
        (
            timeBlockMatrix,
            divUx,
            globalCells,
            qComponent,
            uxComponent,
            1.0,
            &qDivScaleX
        );
        insertTimeBlockOperator
        (
            timeBlockMatrix,
            divUy,
            globalCells,
            qComponent,
            uyComponent,
            1.0,
            &qDivScaleY
        );
        insertTimeBlockOperator
        (
            timeBlockMatrix,
            divUz,
            globalCells,
            qComponent,
            uzComponent,
            1.0,
            &qDivScaleZ
        );

        addTimeBlockCellDiagonal
        (
            timeBlockMatrix,
            globalCells,
            mesh,
            vComponent,
            pComponent,
            vP
        );
        addTimeBlockCellDiagonal
        (
            timeBlockMatrix,
            globalCells,
            mesh,
            vComponent,
            vComponent,
            vV
        );

        insertTimeBlockOperator
        (
            timeBlockMatrix,
            gradPx,
            globalCells,
            uxComponent,
            pComponent,
            1.0,
            &uxGradPScale
        );
        insertTimeBlockOperator
        (
            timeBlockMatrix,
            gradVx,
            globalCells,
            uxComponent,
            vComponent,
            1.0,
            &uxGradVScale
        );
        addTimeBlockCellDiagonal
        (
            timeBlockMatrix,
            globalCells,
            mesh,
            uxComponent,
            uxComponent,
            uxUx
        );

        insertTimeBlockOperator
        (
            timeBlockMatrix,
            gradPy,
            globalCells,
            uyComponent,
            pComponent,
            1.0,
            &uyGradPScale
        );
        insertTimeBlockOperator
        (
            timeBlockMatrix,
            gradVy,
            globalCells,
            uyComponent,
            vComponent,
            1.0,
            &uyGradVScale
        );
        addTimeBlockCellDiagonal
        (
            timeBlockMatrix,
            globalCells,
            mesh,
            uyComponent,
            uyComponent,
            uyUy
        );

        insertTimeBlockOperator
        (
            timeBlockMatrix,
            gradPz,
            globalCells,
            uzComponent,
            pComponent,
            1.0,
            &uzGradPScale
        );
        insertTimeBlockOperator
        (
            timeBlockMatrix,
            gradVz,
            globalCells,
            uzComponent,
            vComponent,
            1.0,
            &uzGradVScale
        );
        addTimeBlockCellDiagonal
        (
            timeBlockMatrix,
            globalCells,
            mesh,
            uzComponent,
            uzComponent,
            uzUz
        );

        MatAssemblyBegin(timeBlockMatrix, MAT_FINAL_ASSEMBLY);
        MatAssemblyEnd(timeBlockMatrix, MAT_FINAL_ASSEMBLY);
        timeBlockMatrixAssembled = true;
        timeBlockMatrixDeltaT = deltaT;

        Info<< "Assembled monolithic six-field PML matrix" << endl;
    }

    const tmp<volScalarField> oldPressureResidual
    (
        previousPressureLaplacian()
      + C0*pDot.oldTime()
      + C1*p.oldTime()
      + C2*v.oldTime()
      - fvc::div(u.oldTime(), "div(u)")
    );

    const tmp<volVectorField> oldPmlResidual
    (
        (TC0 & u.oldTime())
      + (T1 & previousPressureGradient())
      - (T2 & fvc::grad(v.oldTime()))
    );

    scalarField pRhs(mesh.nCells());
    scalarField qRhs(mesh.nCells());
    scalarField vRhs(mesh.nCells());
    scalarField uxRhs(mesh.nCells());
    scalarField uyRhs(mesh.nCells());
    scalarField uzRhs(mesh.nCells());

    forAll(pRhs, celli)
    {
        const scalar volume = cellVolumes[celli];
        const scalar oldP = p.oldTime()[celli];
        const scalar oldQ = pDot.oldTime()[celli];
        const scalar oldV = v.oldTime()[celli];
        const vector& oldU = u.oldTime()[celli];
        const vector& oldPml = oldPmlResidual()[celli];

        pRhs[celli] = volume
           *(oldP/deltaT + oldWeight*oldQ);
        qRhs[celli] = volume
           *(rho[celli]*compressibility[celli]*oldQ/deltaT
           - oldWeight*oldPressureResidual()[celli]);
        vRhs[celli] = volume
           *(oldV/deltaT + oldWeight*oldP);

        uxRhs[celli] = volume
           *(oldU.x()/deltaT - oldWeight*oldPml.x());
        uyRhs[celli] = volume
           *(oldU.y()/deltaT - oldWeight*oldPml.y());
        uzRhs[celli] = volume
           *(oldU.z()/deltaT - oldWeight*oldPml.z());
    }

    addTimeBlockOperatorSource(qRhs, pressureSpatial, theta);
    addTimeBlockOperatorSource(qRhs, divUx, 1.0, &qDivScaleX);
    addTimeBlockOperatorSource(qRhs, divUy, 1.0, &qDivScaleY);
    addTimeBlockOperatorSource(qRhs, divUz, 1.0, &qDivScaleZ);

    addTimeBlockOperatorSource(uxRhs, gradPx, 1.0, &uxGradPScale);
    addTimeBlockOperatorSource(uxRhs, gradVx, 1.0, &uxGradVScale);
    addTimeBlockOperatorSource(uyRhs, gradPy, 1.0, &uyGradPScale);
    addTimeBlockOperatorSource(uyRhs, gradVy, 1.0, &uyGradVScale);
    addTimeBlockOperatorSource(uzRhs, gradPz, 1.0, &uzGradPScale);
    addTimeBlockOperatorSource(uzRhs, gradVz, 1.0, &uzGradVScale);

    setTimeBlockRhs
    (
        timeBlockRhs,
        pRhs,
        qRhs,
        vRhs,
        uxRhs,
        uyRhs,
        uzRhs
    );
    VecAssemblyBegin(timeBlockRhs);
    VecAssemblyEnd(timeBlockRhs);

    KSPSolve(timeBlockKsp, timeBlockRhs, timeBlockSolution);

    KSPConvergedReason reason;
    KSPGetConvergedReason(timeBlockKsp, &reason);
    if (reason < 0)
    {
        FatalErrorInFunction
            << "PETSc time-block solve diverged with reason " << reason
            << exit(FatalError);
    }

    if (timeBlockLinearSolver == "iterative")
    {
        PetscInt iterations = 0;
        PetscReal residualNorm = 0;
        KSPGetIterationNumber(timeBlockKsp, &iterations);
        KSPGetResidualNorm(timeBlockKsp, &residualNorm);

        Info<< "PETSc time-block solve: iterations = " << iterations
            << ", residual = " << residualNorm << endl;
    }

    scatterTimeBlockSolution
    (
        timeBlockSolution,
        globalCells,
        p,
        pDot,
        v,
        pmlUx,
        pmlUy,
        pmlUz,
        u
    );
}
