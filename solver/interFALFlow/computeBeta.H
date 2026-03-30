{
    const scalar alphaMin = 1e-6;
    const scalar alphaMax = 1.0 - alphaMin;

    beta = dimensionedScalar("beta", dimless, 0.0);

    volScalarField isInterface
    (
        IOobject
        (
            "isInterface",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("isInterface", dimless, 0.0)
    );

    forAll(alpha1, celli)
    {
        if (alpha1[celli] > alphaMin && alpha1[celli] < alphaMax)
        {
            isInterface[celli] = 1.0;
        }
    }

    isInterface.correctBoundaryConditions();

    const labelListList& cellCells = mesh.cellCells();

    forAll(alpha1, celli)
    {
        if (isInterface[celli] < 0.5)
        {
            continue;
        }

        beta[celli] = 1.0;

        forAll(cellCells[celli], nI)
        {
            beta[cellCells[celli][nI]] = 1.0;
        }
    }

    const fvBoundaryMesh& patches = mesh.boundary();

    forAll(patches, patchi)
    {
        const fvPatch& patch = patches[patchi];

        if (!patch.coupled() || !isA<processorPolyPatch>(patch.patch()))
        {
            continue;
        }

        const tmp<Field<scalar>> tNbrInterface =
            isInterface.boundaryField()[patchi].patchNeighbourField();
        const Field<scalar>& nbrInterface = tNbrInterface();
        const labelUList& faceCells = patch.faceCells();

        forAll(faceCells, facei)
        {
            if (nbrInterface[facei] > 0.5)
            {
                beta[faceCells[facei]] = 1.0;
            }
        }
    }

    beta.correctBoundaryConditions();
}
