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
    nAlphafFromPlane_(0),
    nAlphafSecondRing_(0),
    nAlphafFallback_(0),
    nAlphafClipped_(0),
    nFluxSecondRing_(0),
    surfX0_(0),
    surfN0_(0),
    surfUn0_(0),
    cellToSurf_(0)
{
   cutFaceAdvect::debug = debug;

    initAreaFractions();

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

    // Update cells using any markedFace
    for (label facei = 0; facei < mesh_.nInternalFaces(); ++facei)
    {
        if (markedFace.test(facei))
        {
            markedCell.set(mesh_.faceOwner()[facei]);
            markedCell.set(mesh_.faceNeighbour()[facei]);
        }
    }
    for (label facei = mesh_.nInternalFaces(); facei < mesh_.nFaces(); ++facei)
    {
        if (markedFace.test(facei))
        {
            markedCell.set(mesh_.faceOwner()[facei]);
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

    // DAVOF: the default face fractions (the upwind cell's alpha) before
    // the planes overwrite the faces they sweep
    const bool storeAf = (storeAreaFractions_ && alphafPtr_);
    if (storeAf)
    {
        defaultAreaFractions(false);
        alphafSetFace_.reset();
        alphafSetFace_.resize(mesh_.nFaces());
        nAlphafFromPlane_ = 0;
        nAlphafSecondRing_ = 0;
        nAlphafFallback_ = 0;
        nAlphafClipped_ = 0;
        nFluxSecondRing_ = 0;
        surfX0_.clear();
        surfN0_.clear();
        surfUn0_.clear();
        cellToSurf_.resize(mesh_.nCells());
        cellToSurf_ = -1;
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
            // dotting it with isoface unit normal
            const scalar Un0 = UInterp.interpolate(x0, celli) & n0;

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
                                false
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
                                false
                            );
                            // The flux is the plane's where the plane can
                            // speak for the upwind cell (planeFlux)
                            dVfIn[facei] =
                                planeFlux(upCell, dVfPlane, dVfIn[facei]);
                        }
                    }

                }
                else
                {
                    bsFaces_.append(facei);
                    bsx0_.append(x0);
                    bsn0_.append(n0);
                    bsUn0_.append(Un0);

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
                        false
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
                storeAreaFraction(facei, advectFace_.subFaceAreaEnd(), false);
            }
        }
    }

    // DAVOF: the faces of the neighbours of the surface cells that no plane
    // has set yet, swept by the neighbouring plane
    if (storeAf)
    {
        sweepSecondRing(dt);
    }

    // Synchronize processor patches
    syncProcPatches(dVf_, phi_, false, storeAf ? alphafPtr_ : nullptr);

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
            << " clips" << endl;
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


void Foam::advection::isoAdvection::defaultAreaFractions
(
    const bool untouchedOnly
)
{
    surfaceScalarField& alphaf = *alphafPtr_;
    const scalarField& alphaIn = alpha1_.primitiveField();
    const scalarField& phiIn = phi_.primitiveField();
    const labelList& own = mesh_.faceOwner();
    const labelList& nei = mesh_.faceNeighbour();

    scalarField& afIn = alphaf.primitiveFieldRef();
    forAll(afIn, facei)
    {
        if (untouchedOnly && alphafSetFace_.test(facei))
        {
            continue;
        }
        afIn[facei] = (phiIn[facei] >= 0) ? alphaIn[own[facei]] : alphaIn[nei[facei]];
    }

    surfaceScalarField::Boundary& bf = alphaf.boundaryFieldRef();
    forAll(bf, patchi)
    {
        const label start = mesh_.boundaryMesh()[patchi].start();
        fvsPatchScalarField& pf = bf[patchi];
        forAll(pf, i)
        {
            const label facei = start + i;
            if (untouchedOnly && alphafSetFace_.test(facei))
            {
                continue;
            }
            pf[i] = alphaIn[own[facei]];
        }
    }
}


void Foam::advection::isoAdvection::storeAreaFraction
(
    const label facei,
    const scalar areaEnd,
    const bool secondRing
)
{
    if (areaEnd < 0)
    {
        ++nAlphafFallback_;
        return;
    }

    scalar af = areaEnd/mag(mesh_.faceAreas()[facei]);
    if (af < -1e-12 || af > 1 + 1e-12)
    {
        ++nAlphafClipped_;
    }
    af = min(max(af, scalar(0)), scalar(1));

    setFaceValue(*alphafPtr_, facei, af);
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
    const cellList& cellFaces = mesh_.cells();
    const labelList& own = mesh_.faceOwner();
    const labelList& nei = mesh_.faceNeighbour();
    const vectorField& Cf = mesh_.faceCentres();
    const vectorField& Sf = mesh_.faceAreas();
    const scalarField& phiIn = phi_.primitiveField();

    labelList bestSurf(mesh_.nFaces(), -1);
    scalarField bestDist(mesh_.nFaces(), GREAT);

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
                const scalar d =
                    mag(((Cf[f] - surfX0_[si]) & surfN0_[si]) - surfUn0_[si]*dt);
                if (d < bestDist[f])
                {
                    bestDist[f] = d;
                    bestSurf[f] = si;
                }
            }
        }
    }

    scalarField& dVfIn = dVf_.primitiveFieldRef();
    forAll(bestSurf, f)
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
        storeAreaFraction(f, advectFace_.subFaceAreaEnd(), true);

        // The flux is the plane's where the plane can speak for the
        // upwind cell (planeFlux); inflow boundary faces keep the boundary
        // condition's flux
        const bool internal = mesh_.isInternalFace(f);
        const label upCell = (internal && phif < 0) ? nei[f] : own[f];
        if (internal)
        {
            dVfIn[f] = planeFlux(upCell, dVf, dVfIn[f]);
        }
        else if (phif >= 0)
        {
            setFaceValue(dVf_, f, planeFlux(upCell, dVf, faceValue(dVf_, f)));
            checkIfOnProcPatch(f);
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
    surfaceScalarField* alphaf
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
                    alphafSetFace_.set(procPatch.start() + facei);
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
