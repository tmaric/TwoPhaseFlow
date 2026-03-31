{
    const scalar alphaMin = 1e-6;
    const scalar alphaMax = 1.0 - alphaMin;

    // Face mask used for interface-neighbour diagnostics/gating.
    beta = dimensionedScalar("beta", dimless, 0.0);

    auto isInterfaceCell = [&](const label celli) -> bool
    {
        const scalar alphaCell = alpha1[celli];
        return alphaCell > alphaMin && alphaCell < alphaMax;
    };

    scalarField& betaInternal = beta.primitiveFieldRef();
    const labelUList& owner = mesh.owner();
    const labelUList& neighbour = mesh.neighbour();

    // Internal faces: mark the face when either side is an interface cell.
    forAll(neighbour, facei)
    {
        if (isInterfaceCell(owner[facei]) || isInterfaceCell(neighbour[facei]))
        {
            betaInternal[facei] = 1.0;
        }
    }

    const fvBoundaryMesh& patches = mesh.boundary();

    forAll(patches, patchi)
    {
        const fvPatch& patch = patches[patchi];
        scalarField& betaPatch = beta.boundaryFieldRef()[patchi];
        const labelUList& faceCells = patch.faceCells();

        if (!patch.coupled() || !isA<processorPolyPatch>(patch.patch()))
        {
            // Non-processor boundary faces only see the owner-side cell.
            forAll(faceCells, facei)
            {
                if (isInterfaceCell(faceCells[facei]))
                {
                    betaPatch[facei] = 1.0;
                }
            }

            continue;
        }

        const tmp<Field<scalar>> tNbrAlpha =
            alpha1.boundaryField()[patchi].patchNeighbourField();
        const Field<scalar>& nbrAlpha = tNbrAlpha();

        // Processor faces also inspect the remote neighbour alpha value.
        forAll(faceCells, facei)
        {
            const scalar alphaNbr = nbrAlpha[facei];

            if
            (
                isInterfaceCell(faceCells[facei])
             || (alphaNbr > alphaMin && alphaNbr < alphaMax)
            )
            {
                betaPatch[facei] = 1.0;
            }
        }
    }

    beta.correctBoundaryConditions();
}
