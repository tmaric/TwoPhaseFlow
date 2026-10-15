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
#include "zoneDistribute.H"
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


bool Foam::reconstruction::davof::positionPlane
(
    const label celli,
    const vector& n
)
{
    // The RETURNED status: vofCutCell returns -1 (alpha below the
    // tolerance) or 1 (above 1 - tolerance) before it cuts anything, and
    // cellStatus() then still holds the previous cell's cut (surfaceArea
    // and surfaceCentre with it). An interface cell is selected here by its
    // face fractions, so its alpha can be at 0 or 1; reading cellStatus()
    // gave such a cell the plane of the previously positioned cell (a
    // phantom surface cell whose plane sits elsewhere, which isoAdvection
    // then fluxed; found on the sphere 2026-10-13).
    const label status = sIterPLIC_.vofCutCell
    (
        celli,
        alpha1_[celli],
        isoFaceTol_,
        100,
        n
    );

    if (status == 0)
    {
        normal_[celli] = sIterPLIC_.surfaceArea();
        centre_[celli] = sIterPLIC_.surfaceCentre();
        if (mag(normal_[celli]) != 0)
        {
            return true;
        }
    }

    normal_[celli] = Zero;
    centre_[celli] = Zero;
    return false;
}


void Foam::reconstruction::davof::updateHasValue()
{
    if (hasValue_.size() == mesh_.nFaces() && !mesh_.topoChanging())
    {
        return;
    }

    // Faces that carry a value: the internal faces and those of patches
    // whose fvsPatchField has the patch size (empty patches have none).
    hasValue_.resize(mesh_.nFaces());
    hasValue_ = false;
    for (label facei = 0; facei < mesh_.nInternalFaces(); ++facei)
    {
        hasValue_[facei] = true;
    }
    const surfaceScalarField::Boundary& bf = alphaf_.boundaryField();
    forAll(bf, patchi)
    {
        const polyPatch& pp = mesh_.boundaryMesh()[patchi];
        if (bf[patchi].size() == pp.size())
        {
            forAll(pp, i)
            {
                hasValue_[pp.start() + i] = true;
            }
        }
    }
}


void Foam::reconstruction::davof::formAreaNormal()
{
    addProfilingInFunction(geometricVoF);
    const vectorField& Sf = mesh_.faceAreas();
    const labelUList& own = mesh_.faceOwner();
    const labelUList& nei = mesh_.faceNeighbour();
    const label nInt = mesh_.nInternalFaces();
    const label nCells = mesh_.nCells();

    vectorField& m = areaNormal_.primitiveFieldRef();

    // The per-cell sums of a pass over every face in face order, the
    // dry faces (alpha_f = 0) skipped (davofInit::areaNormal: 0 S_f adds
    // +-0, which changes no sum, and the sums start at +0). The wet list
    // holds exactly the faces with alpha_f != 0, ascending (the internal
    // faces, then the patch faces in patch order, as the full pass
    // visits them), so every cell receives the same terms in the same
    // order. A cell no wet face touches keeps the +0 it holds: the cells
    // of the previous sum are zeroed, every other cell was never written
    // since the last full zeroing.
    const bool useList =
    (
        wetFacesPtr_
     && wetFacesNFaces_ == mesh_.nFaces()
     && !mesh_.topoChanging()
     && touchedCell_.size() == nCells
     && m.size() == nCells
    );

    if (useList)
    {
        for (const label celli : touchedCell_)
        {
            m[celli] = Zero;
        }
        touchedCell_.reset();

        const scalarField& afIn = alphaf_.primitiveField();
        for (const label facei : *wetFacesPtr_)
        {
            const scalar af = (facei < nInt) ? afIn[facei] : faceFraction(facei);
            if (af == 0) continue;
            const vector contrib = af*Sf[facei];
            m[own[facei]] += contrib;
            touchedCell_.set(own[facei]);
            if (facei < nInt)
            {
                m[nei[facei]] -= contrib;
                touchedCell_.set(nei[facei]);
            }
        }
        return;
    }

    m = Zero;
    touchedCell_.resize(nCells);
    touchedCell_.reset();

    const scalarField& afIn = alphaf_.primitiveField();
    for (label facei = 0; facei < nInt; ++facei)
    {
        const scalar af = afIn[facei];
        if (af == 0) continue;
        const vector contrib = af*Sf[facei];
        m[own[facei]] += contrib;
        m[nei[facei]] -= contrib;
        touchedCell_.set(own[facei]);
        touchedCell_.set(nei[facei]);
    }
    const surfaceScalarField::Boundary& bf = alphaf_.boundaryField();
    forAll(bf, patchi)
    {
        const polyPatch& pp = mesh_.boundaryMesh()[patchi];
        const fvsPatchScalarField& pf = bf[patchi];
        forAll(pf, i)
        {
            const label facei = pp.start() + i;
            if (pf[i] == 0) continue;
            m[own[facei]] += pf[i]*Sf[facei];
            touchedCell_.set(own[facei]);
        }
    }
}


Foam::scalar Foam::reconstruction::davof::stateTime()
{
    const Time& runTime = mesh_.time();

    if (runTime.timeIndex() != lastTimeIndex_)
    {
        lastTimeIndex_ = runTime.timeIndex();
        nReconstructInIndex_ = 0;
    }
    ++nReconstructInIndex_;

    // The constructor (start time index) and every reconstruction after
    // the advection step's own (the first one of a time index) see the
    // state at the current time; the advection step's sees the state the
    // step starts from.
    if
    (
        runTime.timeIndex() == runTime.startTimeIndex()
     || nReconstructInIndex_ > 1
    )
    {
        return runTime.value();
    }

    return runTime.value() - runTime.deltaTValue();
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
    wispTol_(modelDict().getOrDefault<scalar>("wispTol", 1e-3)),
    surfCellTol_(modelDict().getOrDefault<scalar>("surfCellTol", 1e-8)),
    isoFaceTol_(modelDict().getOrDefault<scalar>("isoFaceTol", 1e-8)),
    lastTimeIndex_(-1),
    nReconstructInIndex_(0),
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
    hasValue_(),
    wetFacesPtr_(nullptr),
    wetFacesNFaces_(-1),
    touchedCell_(),
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
    sliverFlag_
    (
        IOobject
        (
            IOobject::groupName("sliverFlag", alpha1.group()),
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimless, Zero)
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
    maxPlaneAreaDefect_(0),
    normalModel_(nullptr),
    usedNormal_(0),
    tensorRequested_(false),
    tensorAlignTol_(modelDict().getOrDefault<scalar>("tensorAlignTol", 0.2)),
    tensorMinCells_(modelDict().getOrDefault<label>("tensorMinCells", 6)),
    writeCurvatureTensor_
    (
        modelDict().getOrDefault<bool>("writeCurvatureTensor", false)
    ),
    polyMoment_(0),
    curvatureTensorPtr_(nullptr),
    curvatureOffsetPtr_(nullptr),
    mRingPtr_(nullptr),
    mRingFilled_(),
    tensorCells_(),
    nTensorCells_(0),
    nTensorFallback_(0)
{
    normalModel_ = davofNormalModel::New(modelDict(), *this);

    Info<< "davof: normalModel " << normalModel_->type() << ", applyTo "
        << normalModel_->applyTo() << ", sliverTol "
        << normalModel_->sliverTol() << ", maxAngle "
        << normalModel_->maxAngle() << endl;

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

    // The interface cells of the previous reconstruction: the planes and
    // the per-cell diagnostics are set in interface cells only and are
    // zero everywhere else (zero at construction), so they are reset
    // there instead of over the whole mesh.
    const labelList prevLabels(interfaceLabels_);
    interfaceCell_ = false;
    interfaceLabels_.clear();

    const scalar tState = stateTime();

    const vectorField& Sf = mesh_.faceAreas();
    const vectorField& Cf = mesh_.faceCentres();
    const vectorField& C = mesh_.cellCentres();
    const scalarField& V = mesh_.V();
    const labelUList& own = mesh_.faceOwner();
    const cellList& cells = mesh_.cells();

    updateHasValue();

    // m_c = sum_f alpha_f S_f^out, formed in areaNormal_ itself over the
    // wet faces (the list isoAdvection keeps) or over every face
    formAreaNormal();
    vectorField& m = areaNormal_.primitiveFieldRef();

    for (const label celli : prevLabels)
    {
        normal_[celli] = Zero;
        centre_[celli] = Zero;
        sliverFlag_[celli] = 0;
        pPlane_[celli] = 0;
        planeVolDefect_[celli] = 0;
        planeAreaDefect_[celli] = 0;
    }
    maxPlaneVolDefect_ = 0;
    maxPlaneAreaDefect_ = 0;

    // A cell whose face fractions are all equal holds no interface (all
    // wet, all dry): its sum_f S_f vanishes only to round-off. Zero it.
    // A cell no wet face touched has m exactly +0 (the sums start at +0
    // and skip the dry faces): all dry, or an exact cancellation, in
    // either case not an interface cell and m already what the rule
    // below would leave (the former full sum of 0 S_f terms gave +0 as
    // well, and the all-dry rule assigned +0), so it is skipped without
    // reading its faces: the cells the sum touched are visited, in cell
    // order (the bitSet iterates ascending), the others hold +0.
    for (const label celli : touchedCell_)
    {
        if (m[celli] == vector::zero) continue;

        const labelList& cf = cells[celli];
        scalar lo = GREAT, hi = -GREAT;
        forAll(cf, i)
        {
            if (!hasValue_[cf[i]]) continue;
            const scalar af = faceFraction(cf[i]);
            lo = min(lo, af);
            hi = max(hi, af);
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

    // The faces of one cell in ascending order: primitiveMesh::cells()
    // lists the owned faces first, then the faces the cell neighbours,
    // while the face moment used to be accumulated by one pass over all
    // faces in face order; the same order keeps its partial sums.
    DynamicList<label> sortedFaces(16);

    // The sliver flag, the explicit position and the realizability defects
    // of the explicit plane (always from m_c; cutCell_ cuts the explicit
    // plane, not the positioned one). The models that read the planes
    // positioned by alpha along m_c (exact, ringFit) get them here, in a
    // first pass; for the other models every plane is positioned once,
    // below, along the normal the model selects. The final plane is the
    // same either way: vofCutCell starts from the cell's vertex values
    // and carries nothing over from an earlier cut, so positioning along
    // the selected normal once gives the plane the former two passes
    // ended with (the second pass re-positioned along that normal, the
    // first pass's plane along m_c was only overwritten).
    const bool twoPass = normalModel_->needsPositionedPlanes();

    // The curvature tensor's per-cell records (requested by the face
    // update's model only; otherwise untouched)
    usedNormal_.resize(interfaceLabels_.size());
    if (tensorRequested_)
    {
        polyMoment_.resize(interfaceLabels_.size());
    }

    forAll(interfaceLabels_, i)
    {
        const label celli = interfaceLabels_[i];
        const scalar mMag = mag(m[celli]);
        const vector n = m[celli]/mMag;

        if (normalModel_->isSliver(celli, mMag))
        {
            sliverFlag_[celli] = 1;
        }

        if (twoPass)
        {
            const bool ok = positionPlane(celli, m[celli]);
            if (tensorRequested_)
            {
                polyMoment_[i] = ok ? cutPolygonMoments() : symmTensor(Zero);
            }
        }

        // The face moment sum_f alpha_f (x_f - x_c).S_f^out of the explicit
        // position, over this cell's faces only (it is read for interface
        // cells only), in face order with the per-face terms of the former
        // pass over all faces
        sortedFaces.clear();
        sortedFaces.append(cells[celli]);
        Foam::sort(sortedFaces);
        scalar faceMoment = 0;
        for (const label facei : sortedFaces)
        {
            if (!hasValue_[facei]) continue;
            const scalar af = faceFraction(facei);
            if (af == 0) continue;
            if (own[facei] == celli)
            {
                faceMoment += af*((Cf[facei] - C[celli]) & Sf[facei]);
            }
            else
            {
                faceMoment -= af*((Cf[facei] - C[celli]) & Sf[facei]);
            }
        }

        const scalar p = (faceMoment - 3.0*alpha1_[celli]*V[celli])/mMag;
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

    if (twoPass)
    {
        // The boundary values of the planes for the ring models' exchange
        normal_.correctBoundaryConditions();
        centre_.correctBoundaryConditions();
    }

    // The model's normals; the plane of every interface cell positioned by
    // alpha along the selected normal (nUse[i] is m_c itself where the
    // model does not replace it), or, with two passes, re-positioned where
    // the normal changed
    vectorField nUse;
    boolList changed;
    normalModel_->select(tState, m, nUse, changed);

    forAll(interfaceLabels_, i)
    {
        if (!twoPass || changed[i])
        {
            const bool ok = positionPlane(interfaceLabels_[i], nUse[i]);
            if (tensorRequested_)
            {
                polyMoment_[i] = ok ? cutPolygonMoments() : symmTensor(Zero);
            }
        }
        const scalar nMag = mag(nUse[i]);
        usedNormal_[i] = (nMag > 0) ? nUse[i]/nMag : vector(Zero);
    }

    maxPlaneVolDefect_ = returnReduce(maxPlaneVolDefect_, maxOp<scalar>());
    maxPlaneAreaDefect_ = returnReduce(maxPlaneAreaDefect_, maxOp<scalar>());
    const label nInterface = returnReduce(interfaceLabels_.size(), sumOp<label>());

    normal_.correctBoundaryConditions();
    centre_.correctBoundaryConditions();
    areaNormal_.correctBoundaryConditions();
    sliverFlag_.correctBoundaryConditions();
    pPlane_.correctBoundaryConditions();
    planeVolDefect_.correctBoundaryConditions();
    planeAreaDefect_.correctBoundaryConditions();

    Info<< "davof: " << nInterface << " interface cells (wispTol " << wispTol_
        << "), ";
    normalModel_->report(Info);
    if (normalModel_->type() == "exact")
    {
        Info<< " at t = " << tState;
    }
    Info<< ", explicit plane defects max |alpha_plane - alpha| = "
        << maxPlaneVolDefect_ << ", max |A n - m|/V^(2/3) = "
        << maxPlaneAreaDefect_ << endl;
}


Foam::symmTensor Foam::reconstruction::davof::cutPolygonMoments()
{
    // The second moments (1/A) int (x - c)(x - c)^T dA of the cut polygon
    // about its centroid c, by the fan of triangles (c, p_i, p_{i+1}): for
    // a triangle with the vertices 0, v1, v2 the integral of x x^T is
    // (A_t/12)(2 v1 v1^T + 2 v2 v2^T + v1 v2^T + v2 v1^T)
    const DynamicList<point>& pts = sIterPLIC_.facePoints();
    const point& c = sIterPLIC_.surfaceCentre();
    const label n = pts.size();

    symmTensor M(Zero);
    scalar A = 0;
    for (label i = 0; i < n; ++i)
    {
        const vector v1 = pts[i] - c;
        const vector v2 = pts[(i + 1) % n] - c;
        const scalar At = 0.5*mag(v1 ^ v2);
        M += (At/6.0)*(sqr(v1) + sqr(v2) + symm(v1*v2));
        A += At;
    }
    if (A > 0)
    {
        return M/A;
    }
    return symmTensor(Zero);
}


void Foam::reconstruction::davof::updateCurvatureTensor()
{
    addProfilingInFunction(geometricVoF);

    if (!tensorRequested_ || polyMoment_.size() != interfaceLabels_.size())
    {
        FatalErrorInFunction
            << "The curvature tensor was not requested before the"
            << " reconstruction (requestCurvatureTensor(true))."
            << exit(FatalError);
    }

    if (!curvatureTensorPtr_)
    {
        const IOobject::writeOption wo =
            writeCurvatureTensor_ ? IOobject::AUTO_WRITE : IOobject::NO_WRITE;
        curvatureTensorPtr_.reset
        (
            new volSymmTensorField
            (
                IOobject
                (
                    IOobject::groupName("curvatureTensorDavof", alpha1_.group()),
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    wo
                ),
                mesh_,
                dimensionedSymmTensor(inv(dimLength), Zero)
            )
        );
        curvatureOffsetPtr_.reset
        (
            new volScalarField
            (
                IOobject
                (
                    IOobject::groupName("curvatureOffsetDavof", alpha1_.group()),
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    wo
                ),
                mesh_,
                dimensionedScalar(dimLength, Zero)
            )
        );
        mRingPtr_.reset
        (
            new volVectorField
            (
                IOobject
                (
                    IOobject::groupName("davofRingTensor", alpha1_.group()),
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    IOobject::NO_REGISTER
                ),
                mesh_,
                dimensionedVector(dimArea, Zero)
            )
        );
    }

    volSymmTensorField& Bf = *curvatureTensorPtr_;
    volScalarField& z0f = *curvatureOffsetPtr_;
    symmTensorField& B = Bf.primitiveFieldRef();
    scalarField& z0 = z0f.primitiveFieldRef();

    // Zero at the cells of the previous update (zero elsewhere)
    for (const label celli : tensorCells_)
    {
        B[celli] = Zero;
        z0[celli] = 0;
    }
    tensorCells_ = interfaceLabels_;

    // The ring exchange (as the ring models): m_c, the plane centroids and
    // normals of the interface cells across the processors
    davofNormalModel::fillInterfaceNormals(*this, *mRingPtr_, mRingFilled_);
    zoneDistribute& exchange = zoneDistribute::New(mesh_);
    const boolList& zone = interfaceCell_;
    exchange.setUpCommforZone(zone, true);

    const Map<vector> mapM(exchange.getDatafromOtherProc(zone, *mRingPtr_));
    const Map<vector> mapX(exchange.getDatafromOtherProc(zone, centre_));
    const Map<vector> mapN(exchange.getDatafromOtherProc(zone, normal_));
    const labelListList& stencil = exchange.getStencil();

    const vectorField& m = areaNormal_.primitiveField();
    const vectorField& C = mesh_.cellCentres();
    const scalarField& V = mesh_.V();

    DynamicList<vector> xs(32);
    DynamicList<vector> ns(32);
    DynamicList<scalar> ws(32);

    nTensorCells_ = interfaceLabels_.size();
    nTensorFallback_ = 0;

    forAll(interfaceLabels_, i)
    {
        const label celli = interfaceLabels_[i];

        // The ring's positioned planes, weighted by |m_k|
        xs.clear(); ns.clear(); ws.clear();
        scalar wSum = 0;
        for (const label gblIdx : stencil[celli])
        {
            const vector mk = exchange.getValue(*mRingPtr_, mapM, gblIdx);
            const scalar mkMag = mag(mk);
            if (!(mkMag > 0)) continue;       // not an interface cell
            const vector nk = exchange.getValue(normal_, mapN, gblIdx);
            const scalar nkMag = mag(nk);
            if (!(nkMag > 0)) continue;       // no positioned plane
            const vector xk = exchange.getValue(centre_, mapX, gblIdx);
            xs.append(xk);
            ns.append(nk/nkMag);
            ws.append(mkMag);
            wSum += mkMag;
        }

        // The frame of the used normal (into the liquid) at the cell's
        // plane centroid
        vector e3 = usedNormal_[i];
        if (!(mag(e3) > 0))
        {
            const scalar mMag = mag(m[celli]);
            e3 = (mMag > 0) ? m[celli]/mMag : vector(Zero);
        }
        const point o = (mag(normal_[celli]) > 0) ? centre_[celli] : C[celli];
        const scalar h = Foam::cbrt(V[celli]);

        vector e1, e2;
        scalar q[6];
        if
        (
            !(mag(e3) > 0)
         || xs.size() < tensorMinCells_
         || !(wSum > 0)
         || !davofNormalModel::mongeFit
            (
                e3, o, h, xs, ns, ws, wSum, tensorAlignTol_, tensorMinCells_,
                q, e1, e2
            )
        )
        {
            ++nTensorFallback_;
            continue;
        }

        // The Hessian of z = h (q0 + q1 x/h + ... + q3 x^2/h^2 + q4 x y/h^2
        // + q5 y^2/h^2) in the frame (e1, e2): [[2 q3, q4], [q4, 2 q5]]/h,
        // mapped to the global frame; the vertex offset from the volume
        // match over the cut polygon
        const scalar B11 = 2.0*q[3]/h;
        const scalar B12 = q[4]/h;
        const scalar B22 = 2.0*q[5]/h;
        const symmTensor Bc =
            B11*sqr(e1) + (2.0*B12)*symm(e1*e2) + B22*sqr(e2);
        B[celli] = Bc;
        z0[celli] = -0.5*(Bc && polyMoment_[i]);
    }

    nTensorCells_ = returnReduce(nTensorCells_, sumOp<label>());
    nTensorFallback_ = returnReduce(nTensorFallback_, sumOp<label>());

    Bf.correctBoundaryConditions();
    z0f.correctBoundaryConditions();

    Info<< "davof: curvature tensor on " << nTensorCells_
        << " cells, " << nTensorFallback_ << " tensor fallbacks" << endl;
}


// ************************************************************************* //
