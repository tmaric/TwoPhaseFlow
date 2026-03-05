{
const auto& faceOwner = mesh.owner();
const cellList& cellFaces = mesh.cells();

alphaf == upwind<scalar>(phi.mesh(), phi).interpolate(alpha1);

 reconstructionSchemes& surf = mesh.lookupObjectRef<reconstructionSchemes>("reconstructionScheme");
 surf.reconstruct(true);

const scalarField& phiIn = phi.primitiveField();
const scalarField& magSfIn = mesh.magSf().primitiveField();
DynamicList<label> bsFaces(label(0.2*mesh.nBoundaryFaces()));
DynamicList<vector> bsx0(bsFaces.size());
DynamicList<vector> bsn0(bsFaces.size());

const DynamicField<label>& interfaceLabels = surf.interfaceLabels();

forAll(interfaceLabels, i)
{
    const label celli = interfaceLabels[i];
    if (mag(surf.normal()[celli]) != 0)
    {
        const point x0 = surf.centre()[celli];
        vector n0 = surf.normal()[celli];
        n0 /= (mag(n0));
        const cell& celliFaces = cellFaces[celli];
        forAll(celliFaces, fi)
        {
            const label facei = celliFaces[fi];

            if (mesh.isInternalFace(facei))
            {
                bool isDownwindFace = false;

                if (celli == faceOwner[facei])
                {
                    if (phiIn[facei] >= 0)
                    { isDownwindFace = true; }
                }
                else
                {
                    if (phiIn[facei] < 0)
                    { isDownwindFace = true; }
                }

                if (isDownwindFace)
                {
                    cutFace.calcSubFace(facei, n0, x0);
                    alphaf[facei] = mag(cutFace.subFaceArea()) / magSfIn[facei];
                }
            }
            else
            {
                bsFaces.append(facei);
                bsx0.append(x0);
                bsn0.append(n0);
            }
        }
    }
}

const polyBoundaryMesh& boundaryMesh = mesh.boundaryMesh();
const surfaceScalarField::Boundary& phib = phi.boundaryField();
const surfaceScalarField::Boundary& magSfb = mesh.magSf().boundaryField();
surfaceScalarField::Boundary& alphafb = alphaf.boundaryFieldRef();
const label nInternalFaces = mesh.nInternalFaces();

forAll(bsFaces, i)
{
    const label facei = bsFaces[i];
    const label patchi = boundaryMesh.patchID()[facei - nInternalFaces];
    const polyPatch& pp = boundaryMesh[patchi];
    const label start = boundaryMesh[patchi].start();
    const label patchFacei = facei - start;

    if (phib[patchi].size())
    {
        const scalar phiP = phib[patchi][patchFacei];

        if(phiP >= 0)
        {
            const scalar magSf = magSfb[patchi][patchFacei];

            cutFace.calcSubFace(facei, bsn0[i], bsx0[i]);
            alphafb[patchi][patchFacei] = mag(cutFace.subFaceArea()) / magSf;

            const cyclicPolyPatch* cpp = isA<cyclicPolyPatch>(pp);
            if (cpp)
            {
                label neiPatchID = cpp->neighbPatchID();
                alphafb[neiPatchID][patchFacei] = alphafb[patchi][patchFacei];
            }
        }
    }
}
bsFaces.clear();
bsx0.clear();
bsn0.clear();
}

Info<< "### display alphaf: "
                << "  Min(" << alphaf.name() << ") = " << min(alphaf).value()
                << "  Max(" << alphaf.name() << ") = " << max(alphaf).value()
                << endl;
