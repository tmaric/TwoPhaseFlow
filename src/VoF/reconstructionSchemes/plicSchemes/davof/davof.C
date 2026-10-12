/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2026 Tomislav Maric, TU Darmstadt
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "davof.H"
#include "davofInit.H"
#include "addToRunTimeSelectionTable.H"
#include "profiling.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
inline namespace twoPhaseFlow
{
namespace reconstruction
{
    defineTypeNameAndDebug(davof, 0);
    addToRunTimeSelectionTable(reconstructionSchemes, davof, components);
}
} // End inline namespace twoPhaseFlow
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::autoPtr<Foam::surfaceScalarField>
Foam::reconstruction::davof::readAlphafIfAbsent(const volScalarField& alpha1)
{
    const fvMesh& mesh = alpha1.mesh();
    const word name(IOobject::groupName("alphaf", alpha1.group()));

    if (mesh.foundObject<surfaceScalarField>(name))
    {
        return nullptr;
    }

    IOobject io
    (
        name,
        mesh.time().timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    );

    if (!io.typeHeaderOk<surfaceScalarField>(true))
    {
        FatalErrorInFunction
            << "The face-fraction field " << name
            << " is neither registered nor in " << mesh.time().timeName()
            << "/ of the case." << nl
            << "The davof reconstruction needs the face fractions of the"
            << " DAVOF state: run initDavofFields (system/initAlphaFieldDict"
            << " with davofMethod) to write " << name << "."
            << exit(FatalError);
    }

    Info<< "davof: reading " << name << " from " << mesh.time().timeName()
        << endl;
    return autoPtr<surfaceScalarField>::New(io, mesh);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::reconstruction::davof::davof
(
    volScalarField& alpha1,
    const surfaceScalarField& phi,
    const volVectorField& U,
    const dictionary& dict
)
:
    reconstructionSchemes
    (
        typeName,
        alpha1,
        phi,
        U,
        dict
    ),
    mesh_(alpha1.mesh()),
    wispTol_(modelDict().getOrDefault<scalar>("wispTol", 1e-6)),
    surfCellTol_(modelDict().getOrDefault<scalar>("surfCellTol", 1e-8)),
    isoFaceTol_(modelDict().getOrDefault<scalar>("isoFaceTol", 1e-8)),
    sIterPLIC_(mesh_, surfCellTol_),
    cutCell_(mesh_),
    alphafPtr_(readAlphafIfAbsent(alpha1)),
    alphaf_
    (
        mesh_.lookupObjectRef<surfaceScalarField>
        (
            IOobject::groupName("alphaf", alpha1.group())
        )
    ),
    areaNormal_
    (
        IOobject
        (
            IOobject::groupName("areaNormalDavof", alpha1.group()),
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedVector(dimArea, Zero)
    ),
    pPlane_
    (
        IOobject
        (
            IOobject::groupName("pPlane", alpha1.group()),
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimLength, Zero)
    ),
    planeVolDefect_
    (
        IOobject
        (
            IOobject::groupName("planeVolDefect", alpha1.group()),
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimless, Zero)
    ),
    planeAreaDefect_
    (
        IOobject
        (
            IOobject::groupName("planeAreaDefect", alpha1.group()),
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimless, Zero)
    ),
    maxPlaneVolDefect_(0),
    maxPlaneAreaDefect_(0)
{
    reconstruct(true);
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::reconstruction::davof::reconstruct(bool forceUpdate)
{
    addProfilingInFunction(geometricVoF);
    const bool uptodate = alreadyReconstructed(forceUpdate);

    if (uptodate && !forceUpdate)
    {
        return;
    }

    if (mesh_.topoChanging())
    {
        if (interfaceCell_.size() != mesh_.nCells())
        {
            interfaceCell_.resize(mesh_.nCells());
        }
    }
    interfaceCell_ = false;
    interfaceLabels_.clear();

    const vectorField& Sf = mesh_.faceAreas();
    const vectorField& Cf = mesh_.faceCentres();
    const vectorField& C = mesh_.cellCentres();
    const scalarField& V = mesh_.V();
    const labelUList& own = mesh_.faceOwner();
    const labelUList& nei = mesh_.faceNeighbour();
    const label nInt = mesh_.nInternalFaces();
    const cellList& cells = mesh_.cells();

    // m_c = sum_f alpha_f S_f^out and the per-polyMesh-face values.
    vectorField m;
    scalarField alphafAll;
    davofInit::areaNormal(mesh_, alphaf_, m, alphafAll);

    // Faces that carry a value: the internal faces and those of patches
    // whose fvsPatchField has the patch size (empty patches have none).
    boolList hasValue(mesh_.nFaces(), false);
    for (label facei = 0; facei < nInt; ++facei)
    {
        hasValue[facei] = true;
    }
    const surfaceScalarField::Boundary& bf = alphaf_.boundaryField();
    forAll(bf, patchi)
    {
        const polyPatch& pp = mesh_.boundaryMesh()[patchi];
        if (bf[patchi].size() == pp.size())
        {
            forAll(pp, i)
            {
                hasValue[pp.start() + i] = true;
            }
        }
    }

    // The face moment sum_f alpha_f (x_f - x_c).S_f^out of the explicit
    // position.
    scalarField faceMoment(mesh_.nCells(), Zero);
    forAll(Sf, facei)
    {
        if (!hasValue[facei]) continue;
        const scalar af = alphafAll[facei];
        if (af == 0) continue;
        faceMoment[own[facei]] += af*((Cf[facei] - C[own[facei]]) & Sf[facei]);
        if (facei < nInt)
        {
            faceMoment[nei[facei]] -=
                af*((Cf[facei] - C[nei[facei]]) & Sf[facei]);
        }
    }

    // A cell whose face fractions are all equal holds no interface (all
    // wet, all dry): its sum_f S_f vanishes only to round-off. Zero it.
    forAll(cells, celli)
    {
        const labelList& cf = cells[celli];
        scalar lo = GREAT, hi = -GREAT;
        forAll(cf, i)
        {
            if (!hasValue[cf[i]]) continue;
            lo = min(lo, alphafAll[cf[i]]);
            hi = max(hi, alphafAll[cf[i]]);
        }
        if (lo == hi)
        {
            m[celli] = Zero;
        }
        else if (mag(m[celli]) > davofInit::wispThreshold(mesh_, celli, wispTol_))
        {
            interfaceCell_[celli] = true;
            interfaceLabels_.append(celli);
        }
    }

    areaNormal_.primitiveFieldRef() = m;
    normal_.primitiveFieldRef() = Zero;
    centre_.primitiveFieldRef() = Zero;
    pPlane_.primitiveFieldRef() = Zero;
    planeVolDefect_.primitiveFieldRef() = Zero;
    planeAreaDefect_.primitiveFieldRef() = Zero;
    maxPlaneVolDefect_ = 0;
    maxPlaneAreaDefect_ = 0;

    forAll(interfaceLabels_, i)
    {
        const label celli = interfaceLabels_[i];
        const scalar mMag = mag(m[celli]);
        const vector n = m[celli]/mMag;

        // The plane positioned by alpha_c along the DAVOF normal.
        sIterPLIC_.vofCutCell
        (
            celli,
            alpha1_[celli],
            isoFaceTol_,
            100,
            m[celli]
        );

        if (sIterPLIC_.cellStatus() == 0)
        {
            normal_[celli] = sIterPLIC_.surfaceArea();
            centre_[celli] = sIterPLIC_.surfaceCentre();
            if (mag(normal_[celli]) == 0)
            {
                normal_[celli] = Zero;
                centre_[celli] = Zero;
            }
        }
        else
        {
            normal_[celli] = Zero;
            centre_[celli] = Zero;
        }

        // The explicit position and the realizability defects of its plane.
        const scalar p = (faceMoment[celli] - 3.0*alpha1_[celli]*V[celli])/mMag;
        pPlane_[celli] = p;

        const label status = cutCell_.calcSubCell(celli, p, n);
        const scalar alphaPlane =
            (status == 0) ? cutCell_.VolumeOfFluid() : ((status < 0) ? 1 : 0);
        const vector areaPlane = (status == 0) ? cutCell_.faceArea() : vector(Zero);

        planeVolDefect_[celli] = mag(alphaPlane - alpha1_[celli]);
        planeAreaDefect_[celli] =
            mag(areaPlane - m[celli])/Foam::pow(V[celli], 2.0/3.0);
        maxPlaneVolDefect_ = max(maxPlaneVolDefect_, planeVolDefect_[celli]);
        maxPlaneAreaDefect_ = max(maxPlaneAreaDefect_, planeAreaDefect_[celli]);
    }

    maxPlaneVolDefect_ = returnReduce(maxPlaneVolDefect_, maxOp<scalar>());
    maxPlaneAreaDefect_ = returnReduce(maxPlaneAreaDefect_, maxOp<scalar>());
    const label nInterface = returnReduce(interfaceLabels_.size(), sumOp<label>());

    normal_.correctBoundaryConditions();
    centre_.correctBoundaryConditions();
    areaNormal_.correctBoundaryConditions();
    pPlane_.correctBoundaryConditions();
    planeVolDefect_.correctBoundaryConditions();
    planeAreaDefect_.correctBoundaryConditions();

    Info<< "davof: " << nInterface << " interface cells (wispTol " << wispTol_
        << "), explicit plane defects max |alpha_plane - alpha| = "
        << maxPlaneVolDefect_ << ", max |A n - m|/V^(2/3) = "
        << maxPlaneAreaDefect_ << endl;
}


// ************************************************************************* //
