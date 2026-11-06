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
#include "emptyPolyPatch.H"
#include "davofInit.H"
#include "davofCcsPositioning.H"
#include "syncTools.H"
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
    if (positioning_ == "ccs")
    {
        // Consecutive Cubic Spline positioning on the same cut, from the
        // explicit position pPlane_[celli] along n (2026-11-04)
        const scalar nMag = mag(n);
        if (nMag == 0)
        {
            normal_[celli] = Zero;
            centre_[celli] = Zero;
            return false;
        }
        const scalar ac = alpha1_[celli];
        Foam::davof::ccsRecord rec;
        const label status = Foam::davof::ccsPositionPlane
        (
            mesh_, cutCellPos_, celli, ac, n/nMag, pPlane_[celli],
            positioningTol_, positioningMaxIter_, rec
        );
        const bool atEnd = (ac < isoFaceTol_) || (1 - ac < isoFaceTol_);
        ++nPosCells_;
        nPosEvals_ += rec.nEvals;
        maxPosEvals_ = max(maxPosEvals_, rec.nEvals);
        ++posEvalHist_[min(rec.nEvals, posEvalHist_.size() - 1)];
        if (rec.exitCode == 2) ++nPosFloor_;
        if (rec.exitCode == 4) ++nPosMaxIter_;
        if (atEnd) ++nPosEnd_;
        const scalar res = mag(rec.rFinal);
        positioningResidual_[celli] = rec.rFinal;   // signed, of V
        maxPosRes_ = max(maxPosRes_, res);
        sumPosRes_ += res;

        if (status == 0)
        {
            normal_[celli] = cutCellPos_.faceArea();
            centre_[celli] = cutCellPos_.faceCentre();
            if (mag(normal_[celli]) != 0)
            {
                return true;
            }
        }
        if (atEnd) ++nPosEndNoPlane_;
        normal_[celli] = Zero;
        centre_[celli] = Zero;
        return false;
    }

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
    emptyFace_.resize(mesh_.nFaces());
    emptyFace_ = false;
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
        else if (isA<emptyPolyPatch>(pp))
        {
            forAll(pp, i)
            {
                emptyFace_[pp.start() + i] = true;
            }
        }
    }
}


void Foam::reconstruction::davof::setFaceFraction(const label facei, const scalar v)
{
    const label nInt = mesh_.nInternalFaces();
    if (facei < nInt)
    {
        alphaf_.primitiveFieldRef()[facei] = v;
        return;
    }
    const polyBoundaryMesh& pbm = mesh_.boundaryMesh();
    const label patchi = pbm.whichPatch(facei);
    if (alphaf_.boundaryField()[patchi].size() != pbm[patchi].size()) return;
    alphaf_.boundaryFieldRef()[patchi][facei - pbm[patchi].start()] = v;
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
    positioning_(modelDict().getOrDefault<word>("positioning", "ccs")),
    interfaceCellRule_(modelDict().getOrDefault<word>("interfaceCellRule", "normalOnly")),
    positioningTol_(modelDict().getOrDefault<scalar>("positioningTol", 1e-14)),
    positioningMaxIter_(modelDict().getOrDefault<label>("positioningMaxIter", 30)),
    cutCellPos_(alpha1.mesh()),
    positioningResidual_
    (
        IOobject
        (
            IOobject::groupName("positioningResidual", alpha1.group()),
            alpha1.mesh().time().timeName(),
            alpha1.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        alpha1.mesh(),
        dimensionedScalar(dimless, Zero)
    ),
    nPosCells_(0),
    nPosEvals_(0),
    maxPosEvals_(0),
    nPosFloor_(0),
    nPosMaxIter_(0),
    nPosEnd_(0),
    nPosEndNoPlane_(0),
    maxPosRes_(0),
    sumPosRes_(0),
    posEvalHist_(8, Zero),
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
    emptyFace_(),
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
    usedNormal_(0)
{
    if (positioning_ != "vofCutCell" && positioning_ != "ccs")
    {
        FatalIOErrorInFunction(modelDict())
            << "positioning '" << positioning_
            << "' unknown. Valid: vofCutCell, ccs." << exit(FatalIOError);
    }
    // decision 3 (2026-11-07): the alpha condition of the interface-cell rule
    if
    (
        interfaceCellRule_ != "normalOnly" && interfaceCellRule_ != "normalAndAlpha"
     && interfaceCellRule_ != "normalAndAlphaReset"
    )
    {
        FatalIOErrorInFunction(modelDict())
            << "interfaceCellRule " << interfaceCellRule_
            << ": normalOnly | normalAndAlpha | normalAndAlphaReset" << exit(FatalIOError);
    }
    Info<< "davof: interfaceCellRule " << interfaceCellRule_ << endl;
    if (positioning_ == "ccs")
    {
        Info<< "davof: positioning ccs (Consecutive Cubic Spline on the"
            << " cutCellPLIC cut from the explicit position, positioningTol "
            << positioningTol_ << " of V, positioningMaxIter "
            << positioningMaxIter_ << ")" << endl;
    }

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

    // decision 3 (2026-11-07): the alpha condition of the interface-cell rule
    const scalar tolA = surfCellTol_;
    auto alphaInside = [&](const label celli) -> bool
    {
        const scalar a = alpha1_[celli];
        return (a > tolA && a < 1 - tolA);
    };
    nRuleOutside_ = nRuleOutsidePlane_ = nRuleResetCells_ = nRuleResetFaces_ = 0;
    ruleResetMaxDaf_ = 0;
    if (interfaceCellRule_ == "normalAndAlphaReset")
    {
        // the faces of a dry or full cell with |m_c| above the threshold are
        // set to its state; a face shared with a neighbour changes the
        // neighbour's m_c as well; no mass moves (alpha is untouched)
        const cellList& cellsR = mesh_.cells();
        labelList resetFlag(mesh_.nFaces(), 0);   // 1: set to 0, 2: set to 1
        bool any = false;
        for (const label celli : touchedCell_)
        {
            if (m[celli] == vector::zero) continue;
            if (alphaInside(celli)) continue;
            if (!(mag(m[celli]) > davofInit::wispThreshold(mesh_, celli, wispTol_))) continue;
            const scalar state = (alpha1_[celli] <= tolA) ? 0 : 1;
            const labelList& cf = cellsR[celli];
            bool changed = false;
            forAll(cf, i)
            {
                const label f = cf[i];
                if (!hasValue_[f]) continue;
                const scalar af = faceFraction(f);
                if (af == state) continue;
                ruleResetMaxDaf_ = max(ruleResetMaxDaf_, mag(af - state));
                setFaceFraction(f, state);
                resetFlag[f] = (state == 0) ? 1 : 2;
                ++nRuleResetFaces_;
                changed = true;
            }
            if (changed)
            {
                ++nRuleResetCells_;
                any = true;
            }
        }
        if (returnReduce(any, orOp<bool>()))
        {
            if (Pstream::parRun())
            {
                syncTools::syncFaceList(mesh_, resetFlag, maxEqOp<label>());
                const polyBoundaryMesh& pbm = mesh_.boundaryMesh();
                forAll(pbm, patchi)
                {
                    if (!pbm[patchi].coupled()) continue;
                    forAll(pbm[patchi], i)
                    {
                        const label f = pbm[patchi].start() + i;
                        if (resetFlag[f] == 1) setFaceFraction(f, 0);
                        else if (resetFlag[f] == 2) setFaceFraction(f, 1);
                    }
                }
            }
            // the wet-face list no longer describes alphaf: the full pass
            wetFacesNFaces_ = -1;
            formAreaNormal();
        }
    }

    for (const label celli : prevLabels)
    {
        normal_[celli] = Zero;
        centre_[celli] = Zero;
        sliverFlag_[celli] = 0;
        pPlane_[celli] = 0;
        planeVolDefect_[celli] = 0;
        planeAreaDefect_[celli] = 0;
        positioningResidual_[celli] = 0;
    }
    maxPlaneVolDefect_ = 0;
    maxPlaneAreaDefect_ = 0;
    nPosCells_ = nPosEvals_ = maxPosEvals_ = nPosFloor_ = nPosMaxIter_ = 0;
    nPosEnd_ = nPosEndNoPlane_ = 0;
    maxPosRes_ = sumPosRes_ = 0;
    posEvalHist_ = Zero;

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
            const bool inside = alphaInside(celli);
            if (!inside) ++nRuleOutside_;
            if (interfaceCellRule_ == "normalOnly" || inside)
            {
                if (!inside) ++nRuleOutsidePlane_;
                interfaceCell_[celli] = true;
                interfaceLabels_.append(celli);
            }
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

    usedNormal_.resize(interfaceLabels_.size());

    forAll(interfaceLabels_, i)
    {
        const label celli = interfaceLabels_[i];
        const scalar mMag = mag(m[celli]);
        const vector n = m[celli]/mMag;

        if (normalModel_->isSliver(celli, mMag))
        {
            sliverFlag_[celli] = 1;
        }

        // The face moment sum_f alpha_f (x_f - x_c).S_f^out of the explicit
        // position, over this cell's faces only (it is read for interface
        // cells only), in face order with the per-face terms of the former
        // pass over all faces. Computed before the first-pass positioning
        // (2026-11-04: positioning ccs starts from pPlane_); the two do not
        // depend on each other, so the order leaves vofCutCell's result as
        // it was.
        sortedFaces.clear();
        sortedFaces.append(cells[celli]);
        Foam::sort(sortedFaces);
        scalar faceMoment = 0;
        for (const label facei : sortedFaces)
        {
            // A face of an empty patch carries no value (the one-cell-thick
            // mesh of the 2D tests): its cross-section is the cell's, so
            // its fraction is alpha_c and its term closes the Gauss sum
            // (2026-11-07; a closed cell has no such face and keeps its
            // sum bit for bit)
            scalar af;
            if (!hasValue_[facei])
            {
                if (!emptyFace_[facei]) continue;
                af = alpha1_[celli];
            }
            else
            {
                af = faceFraction(facei);
            }
            if (af == 0) continue;
            const scalar term = af*((Cf[facei] - C[celli]) & Sf[facei]);
            if (own[facei] == celli)
            {
                faceMoment += term;
            }
            else
            {
                faceMoment -= term;
            }
        }

        const scalar p = (faceMoment - 3.0*alpha1_[celli]*V[celli])/mMag;
        pPlane_[celli] = p;

        if (twoPass)
        {
            positionPlane(celli, m[celli]);
        }

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
            positionPlane(interfaceLabels_[i], nUse[i]);
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
        << maxPlaneAreaDefect_;
    {
        label nO = nRuleOutside_, nOP = nRuleOutsidePlane_;
        label nRc = nRuleResetCells_, nRf = nRuleResetFaces_;
        scalar dMax = ruleResetMaxDaf_;
        if (Pstream::parRun())
        {
            reduce(nO, sumOp<label>()); reduce(nOP, sumOp<label>());
            reduce(nRc, sumOp<label>()); reduce(nRf, sumOp<label>());
            reduce(dMax, maxOp<scalar>());
        }
        Info<< "; interface-cell rule " << interfaceCellRule_
            << ": cells with |m_c| above the threshold and alpha outside (tol, 1 - tol) "
            << nO << ", with a plane " << nOP;
        if (interfaceCellRule_ == "normalAndAlphaReset")
        {
            Info<< ", reset cells " << nRc << ", faces " << nRf
                << ", max |dalphaf| " << dMax;
        }
    }
    if (positioning_ == "ccs")
    {
        label nC = nPosCells_, nE = nPosEvals_, mE = maxPosEvals_;
        label nF = nPosFloor_, nM = nPosMaxIter_, nEnd = nPosEnd_, nEndNo = nPosEndNoPlane_;
        scalar rMax = maxPosRes_, rSum = sumPosRes_;
        labelList hist(posEvalHist_);
        if (Pstream::parRun())
        {
            reduce(nC, sumOp<label>()); reduce(nE, sumOp<label>());
            reduce(mE, maxOp<label>()); reduce(nF, sumOp<label>());
            reduce(nM, sumOp<label>()); reduce(nEnd, sumOp<label>());
            reduce(nEndNo, sumOp<label>());
            reduce(rMax, maxOp<scalar>()); reduce(rSum, sumOp<scalar>());
            forAll(hist, k) reduce(hist[k], sumOp<label>());
        }
        Info<< "; positioning ccs: " << nC << " positionings, cuts per cell "
            << scalar(nE)/max(nC, label(1)) << " (max " << mE
            << ", histogram 1.." << hist.size() - 1 << "+ " << hist
            << "), residual max " << rMax << " mean " << rSum/max(nC, label(1))
            << ", floor exits " << nF << ", maxIter exits " << nM
            << ", end-of-interval cells " << nEnd << " (" << nEndNo
            << " without a plane)";
    }
    Info<< endl;
}


// ************************************************************************* //
