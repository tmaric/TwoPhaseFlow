/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2017 DHI
    Modified code Copyright (C) 2016-2017 OpenCFD Ltd.
    Modified code Copyright (C) 2019 Johan Roenby
    Modified code Copyright (C) 2019-2020 DLR
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

#include "isoAdvection.H"
#include "volFields.H"
#include "interpolationCellPoint.H"
#include "volPointInterpolation.H"
#include "fvcSurfaceIntegrate.H"
#include "fvcGrad.H"
#include "upwind.H"
#include "cellSet.H"
#include "meshTools.H"
#include "OBJstream.H"
#include "syncTools.H"
#include "profiling.H"
#include "davof.H"
#include <cstring>

#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace advection
{
    defineTypeNameAndDebug(isoAdvection, 0);
    addToRunTimeSelectionTable(advectionSchemes,isoAdvection, components);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::advection::isoAdvection::isoAdvection
(
    volScalarField& alpha1,
    const surfaceScalarField& phi,
    const volVectorField& U
)
:
    advectionSchemes
    (
        typeName,
        alpha1,
        phi,
        U
    ),
    // General data
    mesh_(alpha1.mesh()),
    alpha1In_(alpha1.ref()),
    dVf_
    (
        IOobject
        (
            "dVf_",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimVol, Zero)
    ),
    advectionTime_(0),
    timeIndex_(-1),

    // Tolerances and solution controls
    nAlphaBounds_(modelDict().lookupOrDefault<label>("nAlphaBounds", 3)),
    isoFaceTol_(modelDict().lookupOrDefault<scalar>("isoFaceTol", 1e-8)),
    surfCellTol_(modelDict().lookupOrDefault<scalar>("surfCellTol", 1e-8)),
    writeIsoFacesToFile_(modelDict().lookupOrDefault("writeIsoFaces", false)),

    // Cell cutting data
    surfCells_(label(0.2*mesh_.nCells())),
    advectFace_(alpha1.mesh(), alpha1),
    bsFaces_(label(0.2*mesh_.nBoundaryFaces())),
    bsx0_(bsFaces_.size()),
    bsn0_(bsFaces_.size()),
    bsUn0_(bsFaces_.size()),
    bsU0_(0),

    // Parallel run data
    procPatchLabels_(mesh_.boundary().size()),
    surfaceCellFacesOnProcPatches_(0),

    // DAVOF face fractions
    storeAreaFractions_
    (
        modelDict().getOrDefault<bool>
        (
            "storeAreaFractions",
            surf_->type() == "davof"
        )
    ),
    alphafPtr_(nullptr),
    alphafOwnPtr_(nullptr),
    alphafSetFace_(),
    alphafSetFaces_(0),
    alphafPrevSetFaces_(0),
    alphafFallbackFaces_(0),
    defaultsValid_(false),
    alphaDefault_(0),
    phiPos_(),
    phiEventNo_(-1),
    wetFace_(),
    wetFaces_(0),
    davofPtr_(nullptr),
    nAlphafFromPlane_(0),
    nAlphafSecondRing_(0),
    nAlphafFallback_(0),
    nAlphafClipped_(0),
    nFluxSecondRing_(0),
    boundSecondLayer_(modelDict().getOrDefault<bool>("boundSecondLayer", true)),
    nBoundedBeyondRing_(0),
    contentBound_(modelDict().getOrDefault<bool>("contentBound", true)),
    planeFluxUpCell_(0),
    dVfUpwindPlaneFace_(0),
    nContentBoundFaces_(0),
    nContentBoundCells_(0),
    debugCell_(modelDict().getOrDefault<label>("debugCell", -1)),
    surfX0_(0),
    surfN0_(0),
    surfUn0_(0),
    surfU0_(0),
    cellToSurf_(0),
    faceClippingPtr_(nullptr),
    clipCorrects_(false),
    planeFluxFaces_(0),
    secondRingFaces_(0),
    secondRingBest_(0),
    secondRingDist_(0),
    clippingSurface_(modelDict().getOrDefault<word>("clippingSurface", "none")),
    faceState_(modelDict().getOrDefault<word>("faceState", "scalar")),
    wetVecPtr_(nullptr),
    tensorCentroid_(modelDict().getOrDefault<word>("tensorCentroid", "cut")),
    nCentroidFail_(0),
    curvedTolRel_(modelDict().getOrDefault<scalar>("curvedTolRel", 1e-11)),
    curvedMaxDepth_(modelDict().getOrDefault<label>("curvedMaxDepth", 12)),
    fitPriorWeight_(modelDict().getOrDefault<scalar>("fitPriorWeight", 1e-2)),
    fitTol_(modelDict().getOrDefault<scalar>("fitTol", 1e-10)),
    fitMaxIter_(modelDict().getOrDefault<label>("fitMaxIter", 50)),
    fitAnchor_(modelDict().getOrDefault<word>("fitAnchor", "faces")),
    fitDof_(modelDict().getOrDefault<word>("fitDof", "six")),
    anchorTolRel_(modelDict().getOrDefault<scalar>("anchorTolRel", 1e-11)),
    defectEvery_(modelDict().getOrDefault<label>("defectEvery", 10)),
    curvedDebug_(modelDict().getOrDefault<bool>("curvedDebug", false)),
    anchorDiagnostics_(modelDict().getOrDefault<bool>("anchorDiagnostics", false)),
    anchorMethod_(modelDict().getOrDefault<word>("anchorMethod", "ccs")),
    fluxIntegration_(modelDict().getOrDefault<word>("fluxIntegration", "prism")),
    debugR0Max_(0),
    debugCellR0_(-1),
    debugSurface_(),
    debugAlpha_(0),
    nAnchorFail_(0),
    nAnchorLeaves_(0),
    nAnchorEvals_(0),
    nAnchorEdge_(0),
    nCompare_(0),
    nCompareDiff_(0),
    nCompareDiffEdge_(0),
    nCompareDiffResolved_(0),
    anchorShiftMax_(0),
    anchorShiftSum_(0),
    defectLevels_(modelDict().getOrDefault<label>("defectLevels", 3)),
    writeCurvedFields_(modelDict().getOrDefault<bool>("writeCurvedFields", false)),
    cellSurfaces_(0),
    hasSurface_(),
    volumeRulePtr_(nullptr),
    nRuleFallback_(0),
    nCurvedFaces_(0),
    nCurvedSecondRing_(0),
    nCurvedBoundary_(0),
    nCurvedSkipped_(0),
    nCurvedFallback_(0),
    nCurvedLeaves_(0),
    nFitCells_(0),
    nFitConverged_(0),
    nFitStalled_(0),
    nFitSingular_(0),
    nFitNoBracket_(0),
    nFitNoFrame_(0),
    fitResMax_(0),
    fitResSum_(0),
    defectMax_(0),
    defectSum_(0),
    nDefectCells_(0),
    fitResidualPtr_(nullptr),
    paraboloidDefectPtr_(nullptr)
{
   cutFaceAdvect::debug = debug;

    initAreaFractions();

    if (storeAreaFractions_ && alphafPtr_ && isA<reconstruction::davof>(*surf_))
    {
        davofPtr_ = &refCast<reconstruction::davof>(*surf_);

        // The face-clipping model (default plane: nothing changes)
        faceClippingPtr_ = davofFaceClipping::New(modelDict(), *davofPtr_);
        clipCorrects_ = faceClippingPtr_->corrects();
        if (clipCorrects_)
        {
            advectFace_.keepChord(true);
        }
        if (faceClippingPtr_->needsCurvatureTensor())
        {
            // The scheme keeps what the tensor needs from now on; the
            // tensor itself is fitted after the reconstruction at the
            // start of every step (advect)
            davofPtr_->requestCurvatureTensor(true);
        }
        Info<< "isoAdvection: faceClipping " << faceClippingPtr_->type()
            << endl;

        // The curved clipping surface: the dual advection with one surface
        // per cell replaces the plane sweep (isoAdvectionCurved.C)
        if (clippingSurface_ != "none")
        {
            if (clippingSurface_ != "exactSphere" && clippingSurface_ != "paraboloid")
            {
                FatalIOErrorInFunction(modelDict())
                    << "clippingSurface '" << clippingSurface_
                    << "' unknown. Valid: none, exactSphere, paraboloid."
                    << exit(FatalIOError);
            }
            if (clipCorrects_)
            {
                FatalIOErrorInFunction(modelDict())
                    << "clippingSurface " << clippingSurface_
                    << " replaces the plane sweep; faceClipping must be plane"
                    << exit(FatalIOError);
            }
            if (clippingSurface_ == "paraboloid")
            {
                davofPtr_->requestCurvatureTensor(true);
            }
            Info<< "isoAdvection: clippingSurface " << clippingSurface_
                << " (the volume flux and the face update from one surface"
                << " per cell; curvedTolRel " << curvedTolRel_
                << ", curvedMaxDepth " << curvedMaxDepth_;
            if (clippingSurface_ == "paraboloid")
            {
                if (fitAnchor_ != "faces" && fitAnchor_ != "volume")
                {
                    FatalIOErrorInFunction(modelDict())
                        << "fitAnchor '" << fitAnchor_
                        << "' unknown. Valid: faces, volume." << exit(FatalIOError);
                }
                if (fitDof_ != "six" && fitDof_ != "one" && fitDof_ != "none")
                {
                    FatalIOErrorInFunction(modelDict())
                        << "fitDof '" << fitDof_
                        << "' unknown. Valid: six, one, none." << exit(FatalIOError);
                }
                if (fitDof_ == "none" && fitAnchor_ != "volume")
                {
                    FatalIOErrorInFunction(modelDict())
                        << "fitDof none needs fitAnchor volume (the offset"
                        << " comes from nowhere else)" << exit(FatalIOError);
                }
                if
                (
                    anchorMethod_ != "ccs" && anchorMethod_ != "newton"
                 && anchorMethod_ != "illinois" && anchorMethod_ != "compare"
                )
                {
                    FatalIOErrorInFunction(modelDict())
                        << "anchorMethod '" << anchorMethod_
                        << "' unknown. Valid: ccs, newton, illinois, compare."
                        << exit(FatalIOError);
                }
                Info<< ", fitPriorWeight " << fitPriorWeight_ << ", fitTol "
                    << fitTol_ << ", fitMaxIter " << fitMaxIter_
                    << ", fitAnchor " << fitAnchor_ << ", fitDof " << fitDof_
                    << ", anchorTolRel " << anchorTolRel_ << ", anchorMethod "
                    << anchorMethod_ << ", defectEvery "
                    << defectEvery_ << ", anchorDiagnostics "
                    << (anchorDiagnostics_ ? "on" : "off");
            }
            if
            (
                fluxIntegration_ != "prism" && fluxIntegration_ != "trapezoidal"
             && fluxIntegration_ != "simpson" && fluxIntegration_ != "simpsonSplit"
            )
            {
                FatalIOErrorInFunction(modelDict())
                    << "fluxIntegration '" << fluxIntegration_
                    << "' unknown. Valid: prism, trapezoidal, simpson, simpsonSplit."
                    << exit(FatalIOError);
            }
            Info<< ", fluxIntegration " << fluxIntegration_
                << ", defectLevels " << defectLevels_ << ")" << endl;

            // The rule for the volumes and face fractions under a surface
            volumeRulePtr_ = davofVolumeRule::New(mesh_, modelDict());
        }

        // The face state on warped faces (2026-10-25): the wetted area
        // vector next to the scalar fraction, produced by the curved face
        // update only
        if (faceState_ != "scalar" && faceState_ != "vector")
        {
            FatalIOErrorInFunction(modelDict())
                << "faceState '" << faceState_
                << "' unknown. Valid: scalar, vector." << exit(FatalIOError);
        }
        if (tensorCentroid_ != "cut" && tensorCentroid_ != "anchored")
        {
            FatalIOErrorInFunction(modelDict())
                << "tensorCentroid '" << tensorCentroid_
                << "' unknown. Valid: cut, anchored." << exit(FatalIOError);
        }
        if (faceState_ == "vector")
        {
            if (clippingSurface_ == "none")
            {
                FatalIOErrorInFunction(modelDict())
                    << "faceState vector needs a clippingSurface (the plane"
                    << " sweep stores the scalar fraction only)"
                    << exit(FatalIOError);
            }
            wetVecPtr_ = davofPtr_->wettedAreaPtr();
            if (!wetVecPtr_)
            {
                FatalIOErrorInFunction(modelDict())
                    << "faceState vector: the davof scheme holds no wetted"
                    << " area vector field" << exit(FatalIOError);
            }
            Info<< "isoAdvection: faceState vector (the wetted area vector "
                << wetVecPtr_->name() << " transported with alphaf)" << endl;
        }
        if (tensorCentroid_ == "anchored" && clippingSurface_ != "paraboloid")
        {
            FatalIOErrorInFunction(modelDict())
                << "tensorCentroid anchored needs clippingSurface paraboloid"
                << " (the ring tensor is fitted there only)" << exit(FatalIOError);
        }
        if (tensorCentroid_ == "anchored")
        {
            Info<< "isoAdvection: tensorCentroid anchored (the ring tensor's"
                << " centroids re-anchored by the volume rule)" << endl;
        }
    }
    else if (modelDict().getOrDefault<word>("clippingSurface", "none") != "none")
    {
        FatalIOErrorInFunction(modelDict())
            << "clippingSurface needs the davof reconstruction scheme with"
            << " stored face fractions" << exit(FatalIOError);
    }
    else if (modelDict().getOrDefault<word>("faceClipping", "plane") != "plane")
    {
        FatalIOErrorInFunction(modelDict())
            << "faceClipping other than plane needs the davof reconstruction"
            << " scheme with stored face fractions" << exit(FatalIOError);
    }

    // Prepare lists used in parallel runs
    if (Pstream::parRun())
    {
        // Force calculation of required demand driven data (else parallel
        // communication may crash)
        mesh_.cellCentres();
        mesh_.cellVolumes();
        mesh_.faceCentres();
        mesh_.faceAreas();
        mesh_.magSf();
        mesh_.boundaryMesh().patchID();
        mesh_.cellPoints();
        mesh_.cellCells();
        mesh_.cells();

        // Get boundary mesh and resize the list for parallel comms
        setProcessorPatches();
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::advection::isoAdvection::setProcessorPatches()
{
    const polyBoundaryMesh& patches = mesh_.boundaryMesh();
    surfaceCellFacesOnProcPatches_.clear();
    surfaceCellFacesOnProcPatches_.resize(patches.size());

    // Append all processor patch labels to the list
    procPatchLabels_.clear();
    forAll(patches, patchi)
    {
        if
        (
            isA<processorPolyPatch>(patches[patchi])
            && patches[patchi].size() > 0
        )
        {
            procPatchLabels_.append(patchi);
        }
    }
}

void Foam::advection::isoAdvection::extendMarkedCells
(
    bitSet& markedCell
) const
{
    // Mark faces using any marked cell
    bitSet markedFace(mesh_.nFaces());

    for (const label celli : markedCell)
    {
        markedFace.set(mesh_.cells()[celli]);  // set multiple faces
    }

    syncTools::syncFaceList(mesh_, markedFace, orEqOp<unsigned int>());

    // Update cells using any markedFace: the marked faces are visited
    // (the bitSet iterates its set bits, skipping the zero words)
    // instead of testing every face of the mesh; the resulting SET of
    // cells is the same, and its consumers (boundFlux, the second-layer
    // count) iterate the bitSet in cell order whatever the order the
    // bits were set in
    const label nInt = mesh_.nInternalFaces();
    const labelList& own = mesh_.faceOwner();
    const labelList& nei = mesh_.faceNeighbour();
    for (const label facei : markedFace)
    {
        markedCell.set(own[facei]);
        if (facei < nInt)
        {
            markedCell.set(nei[facei]);
        }
    }
}

void Foam::advection::isoAdvection::timeIntegratedFlux()
{
    // Get time step
    addProfilingInFunction(geometricVoF);
    const scalar dt = mesh_.time().deltaTValue();

    // Create object for interpolating velocity to isoface centres
    interpolationCellPoint<vector> UInterp(U_);

    // For each downwind face of each surface cell we "isoadvect" to find dVf
    label nSurfaceCells = 0;

    // Clear out the data for re-use and reset list containing information
    // whether cells could possibly need bounding
    clearIsoFaceData();

    // Get necessary references
    const scalarField& phiIn = phi_.primitiveField();
    const scalarField& magSfIn = mesh_.magSf().primitiveField();
    scalarField& dVfIn = dVf_.primitiveFieldRef();

    // Get necessary mesh data
    const cellList& cellFaces = mesh_.cells();
    const labelList& own = mesh_.faceOwner();
    const labelList& nei = mesh_.faceNeighbour();
    const volVectorField& normal = surf_->normal();
    const boolList& interfaceCell = surf_->interfaceCell();

    // DAVOF: the faces the planes sweep are set below; every other face
    // keeps the value of the previous step until defaultAreaFractions()
    // evaluates it after the bounding (a sweep that falls back writes
    // the default itself; nothing reads an unset face in between)
    const bool storeAf = (storeAreaFractions_ && alphafPtr_);
    if (storeAf)
    {
        alphafSetFace_.reset();
        alphafSetFace_.resize(mesh_.nFaces());
        alphafSetFaces_.clear();
        alphafFallbackFaces_.clear();
        if
        (
            mesh_.topoChanging()
         || alphaDefault_.size() != mesh_.nCells()
         || phiPos_.size() != mesh_.nInternalFaces()
         || wetFace_.size() != mesh_.nFaces()
        )
        {
            defaultsValid_ = false;
            wetFace_.resize(mesh_.nFaces());
        }
        nAlphafFromPlane_ = 0;
        nAlphafSecondRing_ = 0;
        nAlphafFallback_ = 0;
        nAlphafClipped_ = 0;
        nFluxSecondRing_ = 0;
        nContentBoundFaces_ = 0;
        nContentBoundCells_ = 0;
        surfX0_.clear();
        surfN0_.clear();
        surfUn0_.clear();
        surfU0_.clear();
        if (clipCorrects_)
        {
            faceClippingPtr_->resetCounters();
        }

        // The per-face and per-cell records: allocated and filled on the
        // first call or a changed mesh; otherwise they were reset at the
        // entries the previous step set (the end of this function), the
        // same -1 / 0 everywhere as a fill
        if (planeFluxUpCell_.size() != mesh_.nFaces())
        {
            planeFluxUpCell_.resize(mesh_.nFaces());
            planeFluxUpCell_ = -1;
            dVfUpwindPlaneFace_.resize(mesh_.nFaces());
            dVfUpwindPlaneFace_ = 0;
        }
        planeFluxFaces_.clear();
        if (cellToSurf_.size() != mesh_.nCells())
        {
            cellToSurf_.resize(mesh_.nCells());
            cellToSurf_ = -1;
        }
    }


    // Storage for isoFace points. Only used if writeIsoFacesToFile_
    DynamicList<List<point>> isoFacePts;
    const DynamicField<label>& interfaceLabels = surf_->interfaceLabels();

    // Loop through cells
    forAll(interfaceLabels, i)
    {
        const label celli = interfaceLabels[i];
        if (mag(normal[celli]) != 0)
        {

            // This is a surface cell, increment counter, append and mark cell
            nSurfaceCells++;
            surfCells_.append(celli);

            DebugInfo
                << "\n------------ Cell " << celli << " with alpha1 = "
                << alpha1In_[celli] << " and 1-alpha1 = "
                << 1.0 - alpha1In_[celli] << " ------------"
                << endl;

            // Cell is cut
            const point x0 = surf_->centre()[celli];
            vector n0 = -normal[celli];
            n0 /= (mag(n0));

            // Get the speed of the isoface by interpolating velocity and
            // dotting it with isoface unit normal (the velocity itself is
            // kept for the face-clipping models: the same dot product)
            const vector U0 = UInterp.interpolate(x0, celli);
            const scalar Un0 = U0 & n0;

            DebugInfo
                << "calcIsoFace gives initial surface: \nx0 = " << x0
                << ", \nn0 = " << n0 << ", \nUn0 = "
                << Un0 << endl;

            if (storeAf)
            {
                cellToSurf_[celli] = surfX0_.size();
                surfX0_.append(x0);
                surfN0_.append(n0);
                surfUn0_.append(Un0);
                surfU0_.append(U0);
            }

            // Estimate time integrated flux through each downwind face
            // Note: looping over all cell faces - in reduced-D, some of
            //       these faces will be on empty patches
            const cell& celliFaces = cellFaces[celli];
            forAll(celliFaces, fi)
            {
                const label facei = celliFaces[fi];

                if (mesh_.isInternalFace(facei))
                {
                    bool isDownwindFace = false;

                    if (celli == own[facei])
                    {
                        if (phiIn[facei] >= 0)
                        {
                            isDownwindFace = true;
                        }
                    }
                    else
                    {
                        if (phiIn[facei] < 0)
                        {
                            isDownwindFace = true;
                        }
                    }

                    if (isDownwindFace)
                    {
                        dVfIn[facei] = advectFace_.timeIntegratedFaceFlux
                        (
                            facei,
                            x0,
                            n0,
                            Un0,
                            dt,
                            phiIn[facei],
                            magSfIn[facei]
                        );

                        if (storeAf)
                        {
                            // The downwind face: this cell's plane sets it
                            storeAreaFraction
                            (
                                facei,
                                advectFace_.subFaceAreaEnd(),
                                false,
                                celli, x0, n0, U0, dt
                            );
                        }
                    }
                    else if (storeAf)
                    {
                        // The upwind face: the upwind cell's plane sets it
                        // when that cell is a surface cell (its downwind
                        // face); otherwise this cell's plane sweeps it and
                        // only the end area is kept (the flux stays upwind)
                        const label upCell =
                            (celli == own[facei]) ? nei[facei] : own[facei];
                        if (!(interfaceCell[upCell] && mag(normal[upCell]) != 0))
                        {
                            const scalar dVfPlane =
                                advectFace_.timeIntegratedFaceFlux
                                (
                                    facei,
                                    x0,
                                    n0,
                                    Un0,
                                    dt,
                                    phiIn[facei],
                                    magSfIn[facei]
                                );
                            storeAreaFraction
                            (
                                facei,
                                advectFace_.subFaceAreaEnd(),
                                false,
                                celli, x0, n0, U0, dt
                            );
                            // The flux is the plane's where the plane can
                            // speak for the upwind cell (planeFlux)
                            dVfIn[facei] =
                                planeFlux(facei, upCell, dVfPlane, dVfIn[facei]);
                        }
                    }

                }
                else
                {
                    bsFaces_.append(facei);
                    bsx0_.append(x0);
                    bsn0_.append(n0);
                    bsUn0_.append(Un0);
                    if (storeAf)
                    {
                        bsU0_.append(U0);
                    }

                    // Note: we must not check if the face is on the
                    // processor patch here.
                }
            }
        }
    }

    // Get references to boundary fields
    const polyBoundaryMesh& boundaryMesh = mesh_.boundaryMesh();
    const surfaceScalarField::Boundary& phib = phi_.boundaryField();
    const surfaceScalarField::Boundary& magSfb = mesh_.magSf().boundaryField();
    surfaceScalarField::Boundary& dVfb = dVf_.boundaryFieldRef();
    const label nInternalFaces = mesh_.nInternalFaces();

    // Loop through boundary surface faces
    forAll(bsFaces_, i)
    {
        // Get boundary face index (in the global list)
        const label facei = bsFaces_[i];
        const label patchi = boundaryMesh.patchID()[facei - nInternalFaces];
        const label start = boundaryMesh[patchi].start();

        if (phib[patchi].size())
        {
            const label patchFacei = facei - start;
            const scalar phiP = phib[patchi][patchFacei];

            if (phiP >= 0)
            {
                const scalar magSf = magSfb[patchi][patchFacei];

                dVfb[patchi][patchFacei] = advectFace_.timeIntegratedFaceFlux
                (
                    facei,
                    bsx0_[i],
                    bsn0_[i],
                    bsUn0_[i],
                    dt,
                    phiP,
                    magSf
                );

                // Check if the face is on processor patch and append it to
                // the list if necessary
                checkIfOnProcPatch(facei);

                if (storeAf)
                {
                    storeAreaFraction
                    (
                        facei,
                        advectFace_.subFaceAreaEnd(),
                        false,
                        own[facei], bsx0_[i], bsn0_[i], bsU0_[i], dt
                    );
                }
            }
            else if (storeAf)
            {
                // The inflow boundary face of a surface cell: the owner's
                // plane sweeps it, only the end area is kept (the flux is
                // the boundary condition's); on a processor patch the
                // neighbouring surface cell's value replaces it in the sync
                const scalar magSf = magSfb[patchi][patchFacei];
                advectFace_.timeIntegratedFaceFlux
                (
                    facei,
                    bsx0_[i],
                    bsn0_[i],
                    bsUn0_[i],
                    dt,
                    phiP,
                    magSf
                );
                storeAreaFraction
                (
                    facei,
                    advectFace_.subFaceAreaEnd(),
                    false,
                    own[facei], bsx0_[i], bsn0_[i], bsU0_[i], dt
                );
            }
        }
    }

    // DAVOF: the faces of the neighbours of the surface cells that no plane
    // has set yet, swept by the neighbouring plane
    if (storeAf)
    {
        sweepSecondRing(dt);

        // The content bound of the cells the planes set fluxes for
        if (contentBound_)
        {
            boundCellContent(dt);
        }
    }

    // Synchronize processor patches
    syncProcPatches(dVf_, phi_, false, storeAf ? alphafPtr_ : nullptr);

    if (storeAf)
    {
        // The records of this step reset at the entries it set (they are
        // read in this function only)
        for (const label facei : planeFluxFaces_)
        {
            planeFluxUpCell_[facei] = -1;
            dVfUpwindPlaneFace_[facei] = 0;
        }
        for (const label celli : surfCells_)
        {
            cellToSurf_[celli] = -1;
        }
    }

    DebugInfo << "Number of isoAdvector surface cells = "
        << returnReduce(nSurfaceCells, sumOp<label>()) << endl;

    if (storeAf)
    {
        Info<< "isoAdvection: alphaf set by a plane on "
            << returnReduce(nAlphafFromPlane_, sumOp<label>()) << " faces ("
            << returnReduce(nAlphafSecondRing_, sumOp<label>())
            << " in the second ring; the plane's flux on "
            << returnReduce(nFluxSecondRing_, sumOp<label>())
            << " faces beyond the downwind faces), "
            << returnReduce(nAlphafFallback_, sumOp<label>())
            << " fallbacks, "
            << returnReduce(nAlphafClipped_, sumOp<label>())
            << " clips";
        if (contentBound_)
        {
            Info<< "; content bound on "
                << returnReduce(nContentBoundFaces_, sumOp<label>())
                << " faces of "
                << returnReduce(nContentBoundCells_, sumOp<label>())
                << " cells";
        }
        if (clipCorrects_)
        {
            Info<< "; ";
            faceClippingPtr_->report(Info);
        }
        Info<< endl;
    }
}


void Foam::advection::isoAdvection::initAreaFractions()
{
    const word name(IOobject::groupName("alphaf", alpha1_.group()));

    if (mesh_.foundObject<surfaceScalarField>(name))
    {
        alphafPtr_ = &mesh_.lookupObjectRef<surfaceScalarField>(name);
    }
    else if (storeAreaFractions_)
    {
        Info<< "isoAdvection: creating the face fractions " << name
            << " (upwind) for the DAVOF face update" << endl;
        alphafOwnPtr_.reset
        (
            new surfaceScalarField
            (
                IOobject
                (
                    name,
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                upwind<scalar>(mesh_, phi_).interpolate(alpha1_)
            )
        );
        alphafPtr_ = alphafOwnPtr_.get();
    }
}


void Foam::advection::isoAdvection::setDefaultAreaFraction
(
    const label facei,
    const scalarField& alpha
)
{
    const labelList& own = mesh_.faceOwner();

    if (mesh_.isInternalFace(facei))
    {
        const scalar phif = phi_.primitiveField()[facei];
        const scalar af =
            (phif >= 0) ? alpha[own[facei]] : alpha[mesh_.faceNeighbour()[facei]];
        alphafPtr_->primitiveFieldRef()[facei] = af;
        wetFace_.set(facei, (af != 0) ? 1u : 0u);
        if (wetVecPtr_)
        {
            // the planar default of the vector state: alpha_up S_f (the
            // fan triangles' area vectors sum to S_f)
            wetVecPtr_->primitiveFieldRef()[facei] = af*mesh_.faceAreas()[facei];
        }
        return;
    }

    const polyBoundaryMesh& pbm = mesh_.boundaryMesh();
    const label patchi = pbm.patchID()[facei - mesh_.nInternalFaces()];
    fvsPatchScalarField& pf = alphafPtr_->boundaryFieldRef()[patchi];
    if (pf.size())
    {
        const scalar af = alpha[own[facei]];
        pf[facei - pbm[patchi].start()] = af;
        wetFace_.set(facei, (af != 0) ? 1u : 0u);
        if (wetVecPtr_)
        {
            wetVecPtr_->boundaryFieldRef()[patchi][facei - pbm[patchi].start()] =
                af*mesh_.faceAreas()[facei];
        }
    }
}


void Foam::advection::isoAdvection::defaultAreaFractions()
{
    addProfilingInFunction(geometricVoF);
    const scalarField& alphaIn = alpha1_.primitiveField();
    const scalarField& phiIn = phi_.primitiveField();
    const label nInt = mesh_.nInternalFaces();
    const label nFaces = mesh_.nFaces();

    // The target of every face no plane set in this step is alpha of its
    // upwind cell (boundary faces: the owner), with alpha AFTER the
    // bounding. The invariant kept between the steps: such a face holds
    // alphaDefault_ of that cell, alphaDefault_ being a copy of alpha at
    // the last evaluation, as long as nothing has written the face since
    // (a plane, the sync or a fallback writes a face in the step that
    // marks or lists it). The stored value can then differ from the
    // target only on (a) a face whose upwind direction changed, (b) a
    // face set in the previous step, (c) a fallback of this step (it
    // holds the entry alpha of its upwind cell), (d) a face of a cell
    // whose alpha differs bitwise from alphaDefault_. Every other unset
    // face already holds its target; evaluating a face always writes the
    // target, so visiting a face twice is harmless.

    if (!defaultsValid_)
    {
        for (label facei = 0; facei < nFaces; ++facei)
        {
            if (!alphafSetFace_.test(facei))
            {
                setDefaultAreaFraction(facei, alphaIn);
            }
        }
        phiPos_.resize(nInt);
        phiPos_.reset();
        for (label facei = 0; facei < nInt; ++facei)
        {
            if (phiIn[facei] >= 0)
            {
                phiPos_.set(facei);
            }
        }
        phiEventNo_ = phi_.eventNo();
        alphaDefault_ = alphaIn;

        // The wet faces from the final values (the set faces included)
        wetFace_.resize(nFaces);
        wetFace_.reset();
        const scalarField& afIn = alphafPtr_->primitiveField();
        for (label facei = 0; facei < nInt; ++facei)
        {
            if (afIn[facei] != 0) wetFace_.set(facei);
        }
        const surfaceScalarField::Boundary& bf = alphafPtr_->boundaryField();
        forAll(bf, patchi)
        {
            const label start = mesh_.boundaryMesh()[patchi].start();
            const fvsPatchScalarField& pf = bf[patchi];
            forAll(pf, i)
            {
                if (pf[i] != 0) wetFace_.set(start + i);
            }
        }
        defaultsValid_ = true;
        return;
    }

    // (a) phi written since the last evaluation (its event number moved:
    // every non-const access of a GeometricField records one): re-scan
    // the upwind direction, the faces whose direction changed
    if (phi_.eventNo() != phiEventNo_)
    {
        for (label facei = 0; facei < nInt; ++facei)
        {
            const bool pos = (phiIn[facei] >= 0);
            if (pos == phiPos_.test(facei)) continue;
            if (pos) phiPos_.set(facei); else phiPos_.unset(facei);
            if (!alphafSetFace_.test(facei))
            {
                setDefaultAreaFraction(facei, alphaIn);
            }
        }
        phiEventNo_ = phi_.eventNo();
    }

    // (b) the faces set in the previous step and (c) the fallbacks
    for (const label facei : alphafPrevSetFaces_)
    {
        if (!alphafSetFace_.test(facei))
        {
            setDefaultAreaFraction(facei, alphaIn);
        }
    }
    for (const label facei : alphafFallbackFaces_)
    {
        if (!alphafSetFace_.test(facei))
        {
            setDefaultAreaFraction(facei, alphaIn);
        }
    }

    // (d) the faces of the cells whose alpha changed, compared bitwise
    // (-0 and +0 differ in the written field; NaN never occurs here)
    const cellList& cells = mesh_.cells();
    forAll(alphaIn, celli)
    {
        if (std::memcmp(&alphaIn[celli], &alphaDefault_[celli], sizeof(scalar)) == 0)
        {
            continue;
        }
        alphaDefault_[celli] = alphaIn[celli];
        for (const label facei : cells[celli])
        {
            if (!alphafSetFace_.test(facei))
            {
                setDefaultAreaFraction(facei, alphaIn);
            }
        }
    }
}


void Foam::advection::isoAdvection::publishWetFaces()
{
    wetFaces_ = wetFace_.toc();
    if (davofPtr_)
    {
        davofPtr_->setWetFaces(wetFaces_);
    }
}


void Foam::advection::isoAdvection::storeAreaFraction
(
    const label facei,
    const scalar areaEndPlane,
    const bool secondRing,
    const label planeCell,
    const point& x0,
    const vector& n0,
    const vector& U0,
    const scalar dt
)
{
    if (areaEndPlane < 0)
    {
        // The default the former pass before the sweep had written: the
        // upwind cell's alpha at the start of the step (alpha1_ is
        // updated after the sweep); the face stays unset and is
        // re-evaluated with the new alpha after the bounding, or set by
        // the second ring
        ++nAlphafFallback_;
        setDefaultAreaFraction(facei, alpha1_.primitiveField());
        alphafFallbackFaces_.append(facei);
        return;
    }

    // The face-clipping model: the plane-clipped area corrected (plane:
    // never called, the branch is the only addition to the sweep). n0
    // points OUT of the liquid (isoAdvector's convention, n0 =
    // -normal/|normal|), the model takes the normal INTO the liquid; the
    // plane's centroid at the end of the step is x0 + U0 dt
    scalar areaEnd = areaEndPlane;
    if (clipCorrects_)
    {
        areaEnd = faceClippingPtr_->correct
        (
            facei,
            planeCell,
            -n0,
            x0 + U0*dt,
            advectFace_.chordEnd(),
            areaEndPlane
        );
    }

    scalar af = areaEnd/mag(mesh_.faceAreas()[facei]);
    if (af < -1e-12 || af > 1 + 1e-12)
    {
        ++nAlphafClipped_;
    }
    af = min(max(af, scalar(0)), scalar(1));

    setFaceValue(*alphafPtr_, facei, af);
    wetFace_.set(facei, (af != 0) ? 1u : 0u);
    if (!alphafSetFace_.test(facei))
    {
        alphafSetFaces_.append(facei);
    }
    alphafSetFace_.set(facei);
    if (secondRing)
    {
        ++nAlphafSecondRing_;
    }
    else
    {
        ++nAlphafFromPlane_;
    }
}


void Foam::advection::isoAdvection::sweepSecondRing(const scalar dt)
{
    // A face cut by the moving interface at t = dt whose two cells were not
    // surface cells at t = 0 is swept by no plane of its own: it lies in the
    // slab the planes cross. Its cell shares a face with a surface cell
    // (the slab is thinner than a cell at Courant numbers below one); that
    // cell's plane, extended, gives the end-of-step fraction. Among several
    // candidates the plane closest to the face centre at t = dt is taken
    // (one and the same plane for a planar interface).
    addProfilingInFunction(geometricVoF);
    const cellList& cellFaces = mesh_.cells();
    const labelList& own = mesh_.faceOwner();
    const labelList& nei = mesh_.faceNeighbour();
    const vectorField& Cf = mesh_.faceCentres();
    const vectorField& Sf = mesh_.faceAreas();
    const scalarField& phiIn = phi_.primitiveField();

    // The per-face records are sized once and hold -1 / GREAT outside
    // the candidates of the previous sweep, which were reset below
    if (secondRingBest_.size() != mesh_.nFaces())
    {
        secondRingBest_.resize(mesh_.nFaces());
        secondRingBest_ = -1;
        secondRingDist_.resize(mesh_.nFaces());
        secondRingDist_ = GREAT;
    }
    labelList& bestSurf = secondRingBest_;
    scalarField& bestDist = secondRingDist_;
    DynamicLabelList& candidates = secondRingFaces_;
    candidates.clear();

    // A face of a second-ring cell is reached once per surface cell the
    // cell shares a face with (two cells share one face at most), and
    // the same face is reached again by the SAME plane through another
    // second-ring cell that holds it: that visit recomputes the distance
    // it already holds, which cannot beat itself (strict comparison), so
    // it is skipped. The choice among the planes is unchanged: the same
    // surface-cell order, the same strict rule (the first of equal
    // distances wins).
    forAll(surfCells_, si)
    {
        const label celli = surfCells_[si];
        const cell& cf = cellFaces[celli];
        forAll(cf, fi)
        {
            const label g = cf[fi];
            if (!mesh_.isInternalFace(g)) continue;
            const label cellD = (own[g] == celli) ? nei[g] : own[g];
            if (cellToSurf_[cellD] >= 0) continue;   // has its own plane

            const cell& Df = cellFaces[cellD];
            forAll(Df, dfi)
            {
                const label f = Df[dfi];
                if (alphafSetFace_.test(f)) continue;
                if (bestSurf[f] == si) continue;     // visited by this plane
                const scalar d =
                    mag(((Cf[f] - surfX0_[si]) & surfN0_[si]) - surfUn0_[si]*dt);
                if (d < bestDist[f])
                {
                    if (bestSurf[f] < 0)
                    {
                        candidates.append(f);
                    }
                    bestDist[f] = d;
                    bestSurf[f] = si;
                }
            }
        }
    }

    // The candidate faces in face order, as the former pass over all
    // faces visited them (the faces are independent of each other here;
    // the order fixes the sequence of the proc-patch face list only)
    Foam::sort(candidates);

    scalarField& dVfIn = dVf_.primitiveFieldRef();
    for (const label f : candidates)
    {
        const label si = bestSurf[f];
        if (si < 0) continue;
        const scalar magSf = mag(Sf[f]);
        if (magSf <= VSMALL) continue;
        const scalar phif =
            mesh_.isInternalFace(f) ? phiIn[f] : faceValue(phi_, f);
        const scalar dVf = advectFace_.timeIntegratedFaceFlux
        (
            f,
            surfX0_[si],
            surfN0_[si],
            surfUn0_[si],
            dt,
            phif,
            magSf
        );

        // The face fraction: the plane's end value (the swept part or the
        // side value; the upwind default re-evaluated with the new alpha
        // would hand the cell fraction of a cell the plane has entered to
        // all its faces)
        storeAreaFraction
        (
            f,
            advectFace_.subFaceAreaEnd(),
            true,
            surfCells_[si], surfX0_[si], surfN0_[si], surfU0_[si], dt
        );

        // The flux is the plane's where the plane can speak for the
        // upwind cell (planeFlux); inflow boundary faces keep the boundary
        // condition's flux
        const bool internal = mesh_.isInternalFace(f);
        const label upCell = (internal && phif < 0) ? nei[f] : own[f];
        if (internal)
        {
            dVfIn[f] = planeFlux(f, upCell, dVf, dVfIn[f]);
        }
        else if (phif >= 0)
        {
            setFaceValue(dVf_, f, planeFlux(f, upCell, dVf, faceValue(dVf_, f)));
            checkIfOnProcPatch(f);
        }
    }

    // The records of this sweep reset at its candidates
    for (const label f : candidates)
    {
        bestSurf[f] = -1;
        bestDist[f] = GREAT;
    }
}


void Foam::advection::isoAdvection::boundCellContent(const scalar dt)
{
    addProfilingInFunction(geometricVoF);
    const cellList& cellFaces = mesh_.cells();
    const labelList& own = mesh_.faceOwner();
    const scalarField& V = mesh_.V();
    const polyBoundaryMesh& pbm = mesh_.boundaryMesh();
    const surfaceScalarField::Boundary& phib = phi_.boundaryField();

    // The cells with plane-set fluxes (planeFlux: their outflow faces);
    // planeFluxUpCell_ is >= 0 exactly on the faces planeFlux recorded,
    // and the bitSet is visited in cell order whatever the insertion order
    bitSet marked(mesh_.nCells());
    for (const label facei : planeFluxFaces_)
    {
        marked.set(planeFluxUpCell_[facei]);
    }

    for (const label celli : marked)
    {
        const scalar Vc = V[celli];
        const scalar aTol = 1e-12*Vc;
        const scalar alphaW = alpha1In_[celli];

        // The liquid and gas volumes leaving the cell with the current
        // fluxes (1) and with the plane-set faces at upwind (0), and the
        // volumes entering it
        scalar outL0 = 0, outL1 = 0, inL = 0;
        scalar outG0 = 0, outG1 = 0, inG = 0;

        for (const label facei : cellFaces[celli])
        {
            scalar phif;
            if (mesh_.isInternalFace(facei))
            {
                phif = phi_.primitiveField()[facei];
            }
            else
            {
                const label patchi = pbm.patchID()[facei - mesh_.nInternalFaces()];
                if (!phib[patchi].size()) continue;   // empty patch
                phif = phib[patchi][facei - pbm[patchi].start()];
            }
            const bool outflow = (own[facei] == celli) ? (phif >= 0) : (phif < 0);
            const scalar magPhiDt = mag(phif)*dt;
            const scalar l1 = mag(faceValue(dVf_, facei));

            if (outflow)
            {
                const scalar l0 =
                    (planeFluxUpCell_[facei] == celli)
                  ? mag(dVfUpwindPlaneFace_[facei])
                  : l1;
                outL0 += l0; outL1 += l1;
                outG0 += magPhiDt - l0; outG1 += magPhiDt - l1;
            }
            else
            {
                inL += l1;
                inG += magPhiDt - l1;
            }
        }

        // g(lambda) = out(lambda) - in - content <= 0, linear in lambda
        const scalar gL0 = outL0 - inL - alphaW*Vc;
        const scalar gL1 = outL1 - inL - alphaW*Vc;
        const scalar gG0 = outG0 - inG - (1 - alphaW)*Vc;
        const scalar gG1 = outG1 - inG - (1 - alphaW)*Vc;

        scalar lambda = 1;
        if (gL1 > aTol && gL1 > gL0)
        {
            lambda = min(lambda, max(-gL0/(gL1 - gL0), scalar(0)));
        }
        if (gG1 > aTol && gG1 > gG0)
        {
            lambda = min(lambda, max(-gG0/(gG1 - gG0), scalar(0)));
        }

        if (debug && celli == debugCell_)
        {
            Pout<< "boundCellContent: cell " << celli << " alpha " << alphaW
                << " V " << Vc << " outL0 " << outL0 << " outL1 " << outL1
                << " inL " << inL << " outG0 " << outG0 << " outG1 " << outG1
                << " inG " << inG << " lambda " << lambda << nl;
            for (const label facei : cellFaces[celli])
            {
                Pout<< "    face " << facei << " planeSetFor "
                    << planeFluxUpCell_[facei] << " dVf "
                    << faceValue(dVf_, facei) << " upwind "
                    << dVfUpwindPlaneFace_[facei] << " phi dt "
                    << faceValue(phi_, facei)*dt << nl;
            }
        }

        if (lambda >= 1)
        {
            continue;
        }

        ++nContentBoundCells_;
        for (const label facei : cellFaces[celli])
        {
            if (planeFluxUpCell_[facei] != celli) continue;
            const scalar up = dVfUpwindPlaneFace_[facei];
            const scalar cur = faceValue(dVf_, facei);
            setFaceValue(dVf_, facei, up + lambda*(cur - up));
            ++nContentBoundFaces_;
        }
    }
}


bool Foam::advection::isoAdvection::planeSpeaks(const label upCell) const
{
    // The upwind cell is not a surface cell. A pure cell (alpha 0 or 1)
    // is described by the plane only where the plane crosses the face
    // during the step (the swept part; the side value of an extended
    // plane may contradict the cell's content on a curved interface); a
    // partial cell below the interface-cell threshold (a sliver) has no
    // plane of its own and takes the neighbouring plane's value, the
    // swept part or the side value (the upwind value phi alpha dt is
    // neither).
    const scalar a = alpha1In_[upCell];
    const bool pure = (a <= 0 || a >= 1);
    return !(pure && !advectFace_.faceSwept());
}


Foam::scalar Foam::advection::isoAdvection::planeFlux
(
    const label facei,
    const label upCell,
    const scalar dVfPlane,
    const scalar dVfUpwind
)
{
    if (!planeSpeaks(upCell))
    {
        return dVfUpwind;
    }
    ++nFluxSecondRing_;
    planeFluxUpCell_[facei] = upCell;
    dVfUpwindPlaneFace_[facei] = dVfUpwind;
    planeFluxFaces_.append(facei);
    return dVfPlane;
}


void Foam::advection::isoAdvection::setDownwindFaces
(
    const label celli,
    DynamicLabelList& downwindFaces
) const
{
    DebugInFunction << endl;

    // Get necessary mesh data and cell information
    const labelList& own = mesh_.faceOwner();
    const cellList& cells = mesh_.cells();
    const cell& c = cells[celli];

    downwindFaces.clear();

    // Check all faces of the cell
    forAll(c, fi)
    {
        // Get face and corresponding flux
        const label facei = c[fi];
        const scalar phi = faceValue(phi_, facei);

        if (own[facei] == celli)
        {
            if (phi >= 0)
            {
                downwindFaces.append(facei);
            }
        }
        else if (phi < 0)
        {
            downwindFaces.append(facei);
        }
    }

    downwindFaces.shrink();
}


Foam::scalar Foam::advection::isoAdvection::netFlux
(
    const surfaceScalarField& dVf,
    const label celli
) const
{
    scalar dV = 0;

    // Get face indices
    const cell& c = mesh_.cells()[celli];

    // Get mesh data
    const labelList& own = mesh_.faceOwner();

    forAll(c, fi)
    {
        const label facei = c[fi];
        const scalar dVff = faceValue(dVf, facei);

        if (own[facei] == celli)
        {
            dV += dVff;
        }
        else
        {
            dV -= dVff;
        }
    }

    return dV;
}


Foam::DynamicList<Foam::label>  Foam::advection::isoAdvection::syncProcPatches
(
    surfaceScalarField& dVf,
    const surfaceScalarField& phi,
    bool returnSyncedFaces,
    surfaceScalarField* alphaf,
    surfaceVectorField* wetVec
)
{
    DynamicLabelList syncedFaces(0);
    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

    if (Pstream::parRun())
    {
        PstreamBuffers pBufs(Pstream::commsTypes::nonBlocking);

        // Send
        forAll(procPatchLabels_, i)
        {
            const label patchi = procPatchLabels_[i];

            const processorPolyPatch& procPatch =
                refCast<const processorPolyPatch>(patches[patchi]);

            UOPstream toNbr(procPatch.neighbProcNo(), pBufs);
            const scalarField& pFlux = dVf.boundaryField()[patchi];

            const List<label>& surfCellFacesOnProcPatch =
                surfaceCellFacesOnProcPatches_[patchi];

            const UIndirectList<scalar> dVfPatch
            (
                pFlux,
                surfCellFacesOnProcPatch
            );

            toNbr << surfCellFacesOnProcPatch << dVfPatch;

            if (alphaf)
            {
                const UIndirectList<scalar> afPatch
                (
                    alphaf->boundaryField()[patchi],
                    surfCellFacesOnProcPatch
                );
                toNbr << afPatch;
            }
            if (wetVec)
            {
                const UIndirectList<vector> wvPatch
                (
                    wetVec->boundaryField()[patchi],
                    surfCellFacesOnProcPatch
                );
                toNbr << wvPatch;
            }
        }

        pBufs.finishedSends();


        // Receive and combine
        forAll(procPatchLabels_, patchLabeli)
        {
            const label patchi = procPatchLabels_[patchLabeli];

            const processorPolyPatch& procPatch =
                refCast<const processorPolyPatch>(patches[patchi]);

            UIPstream fromNeighb(procPatch.neighbProcNo(), pBufs);
            List<label> faceIDs;
            List<scalar> nbrdVfs;

            fromNeighb >> faceIDs >> nbrdVfs;
            List<scalar> nbrAf;
            if (alphaf)
            {
                fromNeighb >> nbrAf;
            }
            List<vector> nbrWv;
            if (wetVec)
            {
                fromNeighb >> nbrWv;
            }
            if (returnSyncedFaces)
            {
                List<label> syncedFaceI(faceIDs);
                for (label& faceI : syncedFaceI)
                {
                    faceI += procPatch.start();
                }
                syncedFaces.append(syncedFaceI);
            }

            if (debug)
            {
                Pout<< "Received at time = " << mesh_.time().value()
                    << ": surfCellFacesOnProcPatch = " << faceIDs << nl
                    << "Received at time = " << mesh_.time().value()
                    << ": dVfPatch = " << nbrdVfs << endl;
            }

            // Combine fluxes
            scalarField& localFlux = dVf.boundaryFieldRef()[patchi];

            forAll(faceIDs, i)
            {
                const label facei = faceIDs[i];
                localFlux[facei] = - nbrdVfs[i];
                if (alphaf)
                {
                    // The upwind side's face fraction, no sign change
                    alphaf->boundaryFieldRef()[patchi][facei] = nbrAf[i];
                    if (wetVec)
                    {
                        // the wetted area vector is oriented with each
                        // side's own outward patch face: the sign flips,
                        // as the flux's does
                        wetVec->boundaryFieldRef()[patchi][facei] = -nbrWv[i];
                    }
                    const label meshFacei = procPatch.start() + facei;
                    wetFace_.set(meshFacei, (nbrAf[i] != 0) ? 1u : 0u);
                    if (!alphafSetFace_.test(meshFacei))
                    {
                        alphafSetFaces_.append(meshFacei);
                    }
                    alphafSetFace_.set(meshFacei);
                }
                if (debug && mag(localFlux[facei] + nbrdVfs[i]) > ROOTVSMALL)
                {
                    Pout<< "localFlux[facei] = " << localFlux[facei]
                        << " and nbrdVfs[i] = " << nbrdVfs[i]
                        << " for facei = " << facei << endl;
                }
            }
        }

        if (debug)
        {
            // Write out results for checking
            forAll(procPatchLabels_, patchLabeli)
            {
                const label patchi = procPatchLabels_[patchLabeli];
                const scalarField& localFlux = dVf.boundaryField()[patchi];
                Pout<< "time = " << mesh_.time().value() << ": localFlux = "
                    << localFlux << endl;
            }
        }

        // Reinitialising list used for minimal parallel communication
        forAll(surfaceCellFacesOnProcPatches_, patchi)
        {
            surfaceCellFacesOnProcPatches_[patchi].clear();
        }
    }

    return syncedFaces;
}


void Foam::advection::isoAdvection::checkIfOnProcPatch(const label facei)
{
    if (!mesh_.isInternalFace(facei))
    {
        const polyBoundaryMesh& pbm = mesh_.boundaryMesh();
        const label patchi = pbm.patchID()[facei - mesh_.nInternalFaces()];

        if (isA<processorPolyPatch>(pbm[patchi]) && pbm[patchi].size())
        {
            const label patchFacei = pbm[patchi].whichFace(facei);
            surfaceCellFacesOnProcPatches_[patchi].append(patchFacei);
        }
    }
}




void Foam::advection::isoAdvection::applyBruteForceBounding()
{
    addProfilingInFunction(geometricVoF);
    bool alpha1Changed = false;

    scalar snapAlphaTol = modelDict().lookupOrDefault<scalar>("snapTol", 0);
    if (snapAlphaTol > 0)
    {
        alpha1_ =
            alpha1_
           *pos0(alpha1_ - snapAlphaTol)
           *neg0(alpha1_ - (1.0 - snapAlphaTol))
          + pos0(alpha1_ - (1.0 - snapAlphaTol));

        alpha1Changed = true;
    }

    if (modelDict().lookupOrDefault("clip", true))
    {
        alpha1_ = min(scalar(1), max(scalar(0), alpha1_));
        alpha1Changed = true;
    }

    if (alpha1Changed)
    {
        alpha1_.correctBoundaryConditions();
    }
}

// ************************************************************************* //
