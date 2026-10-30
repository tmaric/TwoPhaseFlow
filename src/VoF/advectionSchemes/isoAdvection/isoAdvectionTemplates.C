/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2017 DHI
    Modified code Copyright (C) 2016-2017 OpenCFD Ltd.
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
#include "fvcSurfaceIntegrate.H"
#include "upwind.H"
#include "davof.H"
#include <chrono>

// ************************************************************************* //

template<typename Type>
Type Foam::advection::isoAdvection::faceValue
(
    const GeometricField<Type, fvsPatchField, surfaceMesh>& f,
    const label facei
) const
{
    if (mesh_.isInternalFace(facei))
    {
        return f.primitiveField()[facei];
    }
    else
    {
        const polyBoundaryMesh& pbm = mesh_.boundaryMesh();

        // Boundary face. Find out which face of which patch
        const label patchi = pbm.patchID()[facei - mesh_.nInternalFaces()];

        if (patchi < 0 || patchi >= pbm.size())
        {
            FatalErrorInFunction
                << "Cannot find patch for face " << facei
                << abort(FatalError);
        }

        // Handle empty patches
        const polyPatch& pp = pbm[patchi];
        if (isA<emptyPolyPatch>(pp) || pp.empty())
        {
            return pTraits<Type>::zero;
        }

        const label patchFacei = pp.whichFace(facei);
        return f.boundaryField()[patchi][patchFacei];
    }
}


template<typename Type>
void Foam::advection::isoAdvection::setFaceValue
(
    GeometricField<Type, fvsPatchField, surfaceMesh>& f,
    const label facei,
    const Type& value
) const
{
    if (mesh_.isInternalFace(facei))
    {
        f.primitiveFieldRef()[facei] = value;
    }
    else
    {
        const polyBoundaryMesh& pbm = mesh_.boundaryMesh();

        // Boundary face. Find out which face of which patch
        const label patchi = pbm.patchID()[facei - mesh_.nInternalFaces()];

        if (patchi < 0 || patchi >= pbm.size())
        {
            FatalErrorInFunction
                << "Cannot find patch for face " << facei
                << abort(FatalError);
        }

        // Handle empty patches
        const polyPatch& pp = pbm[patchi];
        if (isA<emptyPolyPatch>(pp) || pp.empty())
        {
            return;
        }

        const label patchFacei = pp.whichFace(facei);

        f.boundaryFieldRef()[patchi][patchFacei] = value;
    }
}


template<class SpType, class SuType>
void Foam::advection::isoAdvection::limitFluxes
(
    const SpType& Sp,
    const SuType& Su
)
{
    addProfilingInFunction(geometricVoF);
    DebugInFunction << endl;

    const scalar aTol = 1.0e-12;          // Note: tolerances
    scalar maxAlphaMinus1 = gMax(alpha1_) - 1;      // max(alphaNew - 1);
    scalar minAlpha = gMin(alpha1_);           // min(alphaNew);
    const label nOvershoots = 20;         // sum(pos0(alphaNew - 1 - aTol));

    const labelList& owner = mesh_.faceOwner();
    const labelList& neighbour = mesh_.faceNeighbour();

    DebugInfo
        << "isoAdvection: Before conservative bounding: min(alpha) = "
        << minAlpha << ", max(alpha) = 1 + " << maxAlphaMinus1 << endl;

    surfaceScalarField dVfcorrectionValues("dVfcorrectionValues", dVf_*0.0);

    bitSet needBounding(mesh_.nCells(),false);
    needBounding.set(surfCells_);

    extendMarkedCells(needBounding);

    // DAVOF: the second-ring sweep (timeIntegratedFlux) sets plane fluxes
    // on the faces of the cells that share a face with a surface cell,
    // i.e. on faces whose other cell lies one layer further out, and the
    // extended plane of a curved interface can flux a full face out of
    // such a cell although it holds no liquid. Extend the bounding by one
    // more layer so that boundFlux corrects those cells too, and count
    // the cells beyond the first layer it corrects.
    const bool extendBounding =
        (storeAreaFractions_ && alphafPtr_ && boundSecondLayer_);
    bitSet firstLayer;
    nBoundedBeyondRing_ = 0;
    if (extendBounding)
    {
        firstLayer = needBounding;
        extendMarkedCells(needBounding);
    }

    // Loop number of bounding steps
    for (label n = 0; n < nAlphaBounds_; n++)
    {
        if (maxAlphaMinus1 > aTol || minAlpha < -aTol) // Note: tolerances
        {
            DebugInfo << "boundAlpha... " << endl;

            if (extendBounding)
            {
                const scalar bTol = 10*SMALL;   // boundFlux's tolerance
                for (const label celli : needBounding)
                {
                    if
                    (
                        !firstLayer.test(celli)
                     && (alpha1_[celli] < -bTol || alpha1_[celli] > 1 + bTol)
                    )
                    {
                        ++nBoundedBeyondRing_;
                    }
                }
            }

            DynamicList<label> correctedFaces(3*nOvershoots);
            dVfcorrectionValues = dimensionedScalar("0",dimVolume,0.0);
            boundFlux(needBounding, dVfcorrectionValues, correctedFaces,Sp,Su);

            correctedFaces.append
            (
                syncProcPatches(dVfcorrectionValues, phi_,true)
            );

            labelHashSet alreadyUpdated;
            forAll(correctedFaces, fi)
            {
                label facei = correctedFaces[fi];
                if (alreadyUpdated.insert(facei))
                {
                    checkIfOnProcPatch(facei);
                    const label own = owner[facei];

                    alpha1_[own] +=
                        -faceValue(dVfcorrectionValues, facei)/mesh_.V()[own];
                    if (mesh_.isInternalFace(facei))
                    {
                        const label nei = neighbour[facei];
                        alpha1_[nei] -=
                            -faceValue(dVfcorrectionValues, facei)/mesh_.V()[nei];
                    }

                    // Change to treat boundaries consistently
                    scalar corrVf =
                        faceValue(dVf_, facei)
                      + faceValue(dVfcorrectionValues, facei);

                    setFaceValue(dVf_, facei, corrVf);
                }
            }
            syncProcPatches(dVf_, phi_);
        }
        else
        {
            break;
        }

        maxAlphaMinus1 = gMax(alpha1_) - 1;     // max(alphaNew - 1);
        minAlpha = gMin(alpha1_);               // min(alphaNew);

        if (debug)
        {
            // Check if still unbounded
            //scalarField alphaNew(alpha1In_ - fvc::surfaceIntegrate(dVf_)());
            label maxAlphaMinus1 = max(alpha1_.primitiveField() - 1);
            scalar minAlpha = min(alpha1_.primitiveField());
            label nUndershoots = sum(neg0(alpha1_.primitiveField() + aTol));
            label nOvershoots = sum(pos0(alpha1_.primitiveField() - 1 - aTol));

            Info<< "After bounding number " << n + 1 << " of time "
                << mesh_.time().value() << ":" << nl
                << "nOvershoots = " << nOvershoots << " with max(alpha1_-1) = "
                << maxAlphaMinus1 << " and nUndershoots = " << nUndershoots
                << " with min(alpha1_) = " << minAlpha << endl;
        }
    }

    alpha1_.correctBoundaryConditions();

    if (extendBounding)
    {
        Info<< "isoAdvection: bounding extended by a second layer: "
            << returnReduce(nBoundedBeyondRing_, sumOp<label>())
            << " cells corrected beyond the first layer" << endl;
    }
}


template<class SpType, class SuType>
void Foam::advection::isoAdvection::boundFlux
(
    const bitSet& nextToInterface,
    surfaceScalarField& dVfcorrectionValues,
    DynamicList<label>& correctedFaces,
    const SpType& Sp,
    const SuType& Su
)
{
    addProfilingInFunction(geometricVoF);
    DebugInFunction << endl;
    scalar rDeltaT = 1/mesh_.time().deltaTValue();

    correctedFaces.clear();
    scalar aTol = 10*SMALL; // Note: tolerances

    const scalarField& meshV = mesh_.cellVolumes();
    const scalar dt = mesh_.time().deltaTValue();

    DynamicList<label> downwindFaces(10);
    DynamicList<label> facesToPassFluidThrough(downwindFaces.size());
    DynamicList<scalar> dVfmax(downwindFaces.size());
    DynamicList<scalar> phi(downwindFaces.size());
    const volScalarField& alphaOld = alpha1_.oldTime();

    // Loop through alpha cell centred field
    for(label celli: nextToInterface)
    {
        if (alpha1_[celli] < -aTol || alpha1_[celli] > 1 + aTol)
        {
            const scalar Vi = meshV[celli];
            scalar alphaOvershoot = pos0(alpha1_[celli]-1)*(alpha1_[celli]-1)
                + neg0(alpha1_[celli])*alpha1_[celli];
            scalar fluidToPassOn = alphaOvershoot*Vi;
            label nFacesToPassFluidThrough = 1;

            bool firstLoop = true;

            // First try to pass surplus fluid on to neighbour cells that are
            // not filled and to which dVf < phi*dt
            for (label iter=0; iter<10; iter++)
            {
                if(mag(alphaOvershoot) < aTol || nFacesToPassFluidThrough == 0)
                {
                    break;
                }

                DebugInfo
                    << "\n\nBounding cell " << celli
                    << " with alpha overshooting " << alphaOvershoot
                    << endl;

                facesToPassFluidThrough.clear();
                dVfmax.clear();
                phi.clear();

                // Find potential neighbour cells to pass surplus phase to
                setDownwindFaces(celli, downwindFaces);

                scalar dVftot = 0;
                nFacesToPassFluidThrough = 0;

                forAll(downwindFaces, fi)
                {
                    const label facei = downwindFaces[fi];
                    const scalar phif = faceValue(phi_, facei);

                    const scalar dVff =
                        faceValue(dVf_, facei)
                      + faceValue(dVfcorrectionValues, facei);

                    const scalar maxExtraFaceFluidTrans =
                        mag(pos0(fluidToPassOn)*phif*dt - dVff);

                    // dVf has same sign as phi and so if phi>0 we have
                    // mag(phi_[facei]*dt) - mag(dVf[facei]) = phi_[facei]*dt
                    // - dVf[facei]
                    // If phi < 0 we have mag(phi_[facei]*dt) -
                    // mag(dVf[facei]) = -phi_[facei]*dt - (-dVf[facei]) > 0
                    // since mag(dVf) < phi*dt
                    DebugInfo
                        << "downwindFace " << facei
                        << " has maxExtraFaceFluidTrans = "
                        << maxExtraFaceFluidTrans << endl;

                    if (maxExtraFaceFluidTrans/Vi > aTol)
                    {
//                    if (maxExtraFaceFluidTrans/Vi > aTol &&
//                    mag(dVfIn[facei])/Vi > aTol) //Last condition may be
//                    important because without this we will flux through uncut
//                    downwind faces
                        facesToPassFluidThrough.append(facei);
                        phi.append(phif);
                        dVfmax.append(maxExtraFaceFluidTrans);
                        dVftot += mag(phif*dt);
                    }
                }

                DebugInfo
                    << "\nfacesToPassFluidThrough: "
                    << facesToPassFluidThrough << ", dVftot = "
                    << dVftot << " m3 corresponding to dalpha = "
                    << dVftot/Vi << endl;

                forAll(facesToPassFluidThrough, fi)
                {
                    const label facei = facesToPassFluidThrough[fi];
                    scalar fluidToPassThroughFace =
                        mag(fluidToPassOn)*mag(phi[fi]*dt)/dVftot;

                    nFacesToPassFluidThrough +=
                        pos0(dVfmax[fi] - fluidToPassThroughFace);

                    fluidToPassThroughFace =
                        min(fluidToPassThroughFace, dVfmax[fi]);

                    scalar dVff = faceValue(dVfcorrectionValues, facei);

                    dVff +=
                        sign(phi[fi])*fluidToPassThroughFace*sign(fluidToPassOn);

                    setFaceValue(dVfcorrectionValues, facei, dVff);

                    if (firstLoop)
                    {
                        checkIfOnProcPatch(facei);
                        correctedFaces.append(facei);
                    }
                }

                firstLoop = false;

                scalar alpha1New =
                (
                    alphaOld[celli]*rDeltaT  + Su[celli]
                  - netFlux(dVf_, celli)/Vi*rDeltaT
                  - netFlux(dVfcorrectionValues, celli)/Vi*rDeltaT
                )
                /
                (rDeltaT - Sp[celli]);

                alphaOvershoot =
                    pos0(alpha1New-1)*(alpha1New-1)
                  + neg0(alpha1New)*alpha1New;

                fluidToPassOn = alphaOvershoot*Vi;

                DebugInfo
                    << "\nNew alpha for cell " << celli << ": "
                    << alpha1New << endl;
            }
        }
    }

    DebugInfo << "correctedFaces = " << correctedFaces << endl;
}


template<class SpType, class SuType>
void Foam::advection::isoAdvection::advect(const SpType& Sp, const SuType& Su)
{
    addProfilingInFunction(geometricVoF);
    DebugInFunction << endl;

    if (mesh_.topoChanging())
    {
        setProcessorPatches();
    }

    scalar advectionStartTime = mesh_.time().elapsedCpuTime();
    const auto advectWall0 = std::chrono::steady_clock::now();

    scalar rDeltaT = 1/mesh_.time().deltaTValue();

    // reconstruct the interface
    surf_->reconstruct();

    // DAVOF: the curvature tensor of the start-of-step planes for the
    // face-clipping model that corrects the plane clipping (planeSagitta);
    // it serves the whole step
    if (davofPtr_ && faceClippingPtr_ && faceClippingPtr_->needsCurvatureTensor())
    {
        davofPtr_->updateCurvatureTensor();
    }

    if(timeIndex_ < mesh_.time().timeIndex())
    {
        timeIndex_= mesh_.time().timeIndex();
        surf_->normal().oldTime() = surf_->normal();
        surf_->centre().oldTime() = surf_->centre();
    }

    // Initialising dVf with upwind values
    // i.e. phi[facei]*alpha1[upwindCell[facei]]*dt
    dVf_ = upwind<scalar>(mesh_, phi_).flux(alpha1_)*mesh_.time().deltaT();

    // Do the isoAdvection on surface cells; DAVOF with a curved clipping
    // surface: the dual advection with one surface per cell instead
    if (davofPtr_ && clippingSurface_ != "none")
    {
        curvedFlux(mesh_.time().deltaTValue());
        reportCurved();
    }
    else
    {
        timeIntegratedFlux();
    }

    // Adjust alpha for mesh motion
    if (mesh_.moving())
    {
        alpha1In_ *= (mesh_.Vsc0()/mesh_.Vsc());
    }

    // Advect the free surface
    alpha1_.primitiveFieldRef() =
    (
        alpha1_.oldTime().primitiveField()*rDeltaT
      + Su.field()
      - fvc::surfaceIntegrate(dVf_)().primitiveField()*rDeltaT
    )/(rDeltaT - Sp.field());

    alpha1_.correctBoundaryConditions();

    // DAVOF curved loop: the mass sum(alpha V) through the stages of the
    // step, relative to its value at t^n (2026-10-29; a print only)
    const bool massReport = (davofPtr_ && clippingSurface_ != "none");
    scalar massOld = 0, massAdv = 0;
    if (massReport)
    {
        massOld = gSum(alpha1_.oldTime().primitiveField()*mesh_.V());
        massAdv = gSum(alpha1In_*mesh_.V());
    }
    if (sweptClassify_)
    {
        dVfBeforeBound_ = dVf_.primitiveField();
    }

    // Adjust dVf for unbounded cells
    limitFluxes
    (
        Sp,
        Su
    );

    scalar maxAlphaMinus1 = gMax(alpha1In_) - 1;
    scalar minAlpha = gMin(alpha1In_);

    Info<< "isoAdvection: After conservative bounding: min(alpha) = "
        << minAlpha << ", max(alpha) = 1 + " << maxAlphaMinus1 << endl;

    // The first cell that leaves [0, 1] by more than 1e-10 after the
    // conservative bounding (sweptClassify): the cell, its alpha before
    // and after the step, its faces with phi, dVf before and after the
    // bounding, the pre-image class and the other cell's state (serial)
    if (sweptClassify_ && !sweptEventReported_ && !Pstream::parRun())
    {
        label worst = -1;
        scalar excess = 1e-10;
        forAll(alpha1In_, celli)
        {
            const scalar e = max(-alpha1In_[celli], alpha1In_[celli] - 1);
            if (e > excess) { excess = e; worst = celli; }
        }
        if (worst >= 0)
        {
            sweptEventReported_ = true;
            const scalarField& a0 = alpha1_.oldTime().primitiveField();
            const labelList& own = mesh_.faceOwner();
            const labelList& nei = mesh_.faceNeighbour();
            static const char* clsName[6] =
                {"donorOnly", "interfaceNbr", "fullNbr", "emptyNbr", "outside", "notSwept"};
            Info<< "isoAdvection: FIRST UNBOUNDED CELL after the bounding: step "
                << mesh_.time().timeIndex() << " cell " << worst
                << " centre " << mesh_.cellCentres()[worst]
                << " alpha old " << a0[worst] << " new " << alpha1In_[worst]
                << " surface " << (hasSurface_.test(worst) ? "yes" : "no")
                << " V " << mesh_.V()[worst] << nl;
            for (const label facei : mesh_.cells()[worst])
            {
                const bool internal = mesh_.isInternalFace(facei);
                const label other = internal ? (own[facei] == worst ? nei[facei] : own[facei]) : -1;
                const scalar sgn = (own[facei] == worst) ? 1 : -1;   // out of the cell
                const label cls = (sweptClass_.size() > facei) ? sweptClass_[facei] : -1;
                Info<< "    face " << facei << (internal ? "" : " (boundary)")
                    << " phi out " << sgn*faceValue(phi_, facei)
                    << " dVf out before " << sgn*(internal ? dVfBeforeBound_[facei] : faceValue(dVf_, facei))
                    << " after " << sgn*faceValue(dVf_, facei)
                    << " class " << clsName[(cls >= 0 && cls < 5) ? cls : 5];
                if (other >= 0)
                {
                    Info<< " other cell " << other << " alpha old " << a0[other]
                        << " new " << alpha1In_[other]
                        << " surface " << (hasSurface_.test(other) ? "yes" : "no");
                }
                Info<< nl;
            }
            Info<< endl;
        }
    }

    // Apply non-conservative bounding mechanisms (clipping and snapping)
    // Note: We should be able to write out alpha before this is done!
    applyBruteForceBounding();

    if (massReport)
    {
        const scalar massBound = gSum(alpha1In_*mesh_.V());
        const scalar ref = max(mag(massOld), VSMALL);
        // the liquid that left through the non-coupled boundary patches in
        // the step (the advection term's only source besides the explicit
        // source Su, since the interior fluxes cancel in surfaceIntegrate)
        scalar bnd = 0;
        forAll(dVf_.boundaryField(), patchi)
        {
            if (!dVf_.boundaryField()[patchi].coupled())
            {
                bnd += sum(dVf_.boundaryField()[patchi]);
            }
        }
        reduce(bnd, sumOp<scalar>());
        Info<< "isoAdvection: mass change relative to sum(alpha V) at t^n: advection "
            << (massAdv - massOld)/ref << ", bounding " << (massBound - massAdv)/ref
            << ", clipping " << (gSum(alpha1In_*mesh_.V()) - massBound)/ref
            << ", boundary outflow " << bnd/ref << endl;
    }

    // DAVOF: the faces no plane has set carry the upwind cell's new alpha
    if (storeAreaFractions_ && alphafPtr_)
    {
        defaultAreaFractions();

        // The faces set in this step are the previous ones of the next
        alphafPrevSetFaces_.swap(alphafSetFaces_);
        alphafSetFaces_.clear();

        // The faces with a non-zero fraction, ascending, for the davof
        // scheme's area sums
        publishWetFaces();
    }

    // The swept fluxes: the remote cells of the parallel path and the wall
    // times of the step's advection (2026-10-30; the maxima over the ranks,
    // the counts summed)
    if
    (
        davofPtr_ && clippingSurface_ != "none"
     && (fluxIntegration_ == "sweptFace" || fluxIntegration_ == "sweptFaceCells")
    )
    {
        scalar tAdv =
            std::chrono::duration<scalar>(std::chrono::steady_clock::now() - advectWall0).count();
        scalar tCurved = sweptCurvedTime_, tExch = sweptExchangeTime_;
        label nRc = nSweptRemoteCells_, nRf = nSweptRemoteFaces_, nEu = nSweptEulerTrace_;
        reduce(tAdv, maxOp<scalar>());
        reduce(tCurved, maxOp<scalar>());
        reduce(tExch, maxOp<scalar>());
        reduce(nRc, sumOp<label>());
        reduce(nRf, sumOp<label>());
        reduce(nEu, sumOp<label>());
        Info<< "isoAdvection: swept parallel: remote cells " << nRc
            << ", faces with remote cells " << nRf
            << ", Euler traces " << nEu
            << "; wall time max over " << Pstream::nProcs() << " ranks: advect "
            << tAdv << " s, curvedFlux " << tCurved << " s, exchange " << tExch
            << " s" << endl;
    }

    advectionTime_ += (mesh_.time().elapsedCpuTime() - advectionStartTime);
    DebugInfo
        << "isoAdvection: time consumption = "
        << label(100*advectionTime_/(mesh_.time().elapsedCpuTime() + SMALL))
        << "%" << endl;

    alphaPhi_ = dVf_/mesh_.time().deltaT();
}


// ************************************************************************* //
