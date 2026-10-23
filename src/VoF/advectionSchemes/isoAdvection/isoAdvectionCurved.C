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

Description
    The DAVOF dual advection with ONE curved clipping surface per cell
    (2026-10-23; dictionary key clippingSurface: none | exactSphere |
    paraboloid): the volume flux through a face is the liquid volume of the
    prism swept between the face's pre-image f - U dt and the face under the
    surface of the upwind cell, and the new face fraction is the liquid
    fraction of that pre-image under the same surface, so that alpha_k and
    alpha_f are the volume and the boundary moments of one surface at every
    step (the requirement of the sagitta studies, 2026-10-16). Replaces the
    plane sweep of timeIntegratedFlux() when active; serial runs only.

    exactSphere: the sphere of the exactNormal dictionary (origin + U0 t at
    the state's time) is every cell's surface: the machinery gate, exact to
    the quadrature tolerance under uniform translation by construction.

    paraboloid: the six-coefficient paraboloid of each interface cell in the
    frame of the normal the scheme used (offset, two slopes, the curvature
    tensor), all six fitted to the cell's own face fractions by the
    Gauss-Newton of leia's regeneration (2026-10-12) with the ring curvature
    tensor as the prior, from the one-scalar anchoring as the start; cells
    whose fit does not converge keep the best iterate (counted).

SourceFiles
    isoAdvectionCurved.C

\*---------------------------------------------------------------------------*/

#include "isoAdvection.H"
#include "davofCurvedGeometry.H"
#include "interpolationCellPoint.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::advection::isoAdvection::buildCellSurfaces(const scalar tState)
{
    addProfilingInFunction(geometricVoF);
    const label nCells = mesh_.nCells();
    if (cellSurfaces_.size() != nCells)
    {
        cellSurfaces_.setSize(nCells);
    }
    hasSurface_.reset();
    hasSurface_.resize(nCells);

    nFitCells_ = nFitConverged_ = nFitStalled_ = nFitSingular_ = nFitNoBracket_ = 0;
    nFitNoFrame_ = 0;
    debugR0Max_ = 0;
    debugCellR0_ = -1;
    fitResMax_ = 0;
    fitResSum_ = 0;
    nAnchorFail_ = 0;
    nAnchorLeaves_ = 0;
    anchorShiftMax_ = 0;
    anchorShiftSum_ = 0;

    const DynamicField<label>& interfaceLabels = surf_->interfaceLabels();

    if (clippingSurface_ == "exactSphere")
    {
        const dictionary& ed = modelDict().subDict("exactNormal");
        const point origin0 = ed.get<point>("origin");
        const vector U0 = ed.getOrDefault<vector>("U0", vector::zero);
        const scalar R = ed.get<scalar>("radius");
        davof::cellSurface S;
        S.sphere = true;
        S.o = origin0 + U0*tState;
        S.R = R;
        forAll(interfaceLabels, i)
        {
            const label celli = interfaceLabels[i];
            cellSurfaces_[celli] = S;
            hasSurface_.set(celli);
        }
        return;
    }

    // paraboloid: the ring curvature tensor of the start-of-step state
    davofPtr_->updateCurvatureTensor();
    const symmTensorField& B = davofPtr_->curvatureTensor().primitiveField();
    const vectorField& usedN = davofPtr_->usedNormals();
    const vectorField& normal = surf_->normal().primitiveField();
    const vectorField& centre = surf_->centre().primitiveField();
    const vectorField& m = davofPtr_->areaNormal().primitiveField();
    const vectorField& C = mesh_.cellCentres();
    const scalarField& V = mesh_.V();
    const vectorField& Sf = mesh_.faceAreas();
    const cellList& cells = mesh_.cells();

    if (writeCurvedFields_ && !fitResidualPtr_)
    {
        fitResidualPtr_.reset
        (
            new volScalarField
            (
                IOobject
                (
                    IOobject::groupName("fitResidualDavof", alpha1_.group()),
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh_,
                dimensionedScalar(dimless, Zero)
            )
        );
    }
    if (fitResidualPtr_)
    {
        fitResidualPtr_->primitiveFieldRef() = 0;
    }

    forAll(interfaceLabels, i)
    {
        const label celli = interfaceLabels[i];

        // The frame: e3 points out of the liquid (the used normal points in)
        vector nIn = usedN[i];
        if (!(mag(nIn) > 0))
        {
            const scalar mMag = mag(m[celli]);
            if (!(mMag > 0))
            {
                ++nFitNoFrame_;
                continue;
            }
            nIn = m[celli]/mMag;
        }
        davof::cellSurface S;
        S.e3 = -nIn;
        davof::frameOf(S.e3, S.e1, S.e2);
        S.o = (mag(normal[celli]) > 0) ? centre[celli] : C[celli];
        S.c = 0;
        const symmTensor& Bc = B[celli];
        S.b11 = -(S.e1 & (Bc & S.e1));
        S.b12 = -(S.e1 & (Bc & S.e2));
        S.b22 = -(S.e2 & (Bc & S.e2));
        S.paraboloid = true;
        const scalar h = Foam::cbrt(V[celli]);

        // The targets: the cell's own face fractions, weighted by the areas
        const labelList& cf = cells[celli];
        const label F = cf.size();
        scalarList tF(F), wF(F);
        scalar sumS = 0, target = 0;
        forAll(cf, k) sumS += mag(Sf[cf[k]]);
        forAll(cf, k)
        {
            const label f = cf[k];
            tF[k] = faceValue(*alphafPtr_, f);
            wF[k] = mag(Sf[f])/max(sumS, VSMALL);
            target += mag(Sf[f])*tF[k];
        }

        // The fractions of the cell's faces under a surface, in face order
        auto fractionsOf = [&](const davof::cellSurface& T, scalarList& aF)
        {
            aF.setSize(F);
            label nfb = 0;
            forAll(cf, k)
            {
                aF[k] = davof::faceCutDisplaced
                (
                    mesh_, cf[k], vector::zero,
                    [&](const point& x) { return T.value(x); }, nfb
                );
            }
        };

        // The start: the offset that reproduces the total wetted area with
        // the ring tensor and zero slopes (leia's paraboloidFaces anchoring);
        // skipped when the offset comes from the volume alone (fitDof none
        // with fitAnchor volume)
        if (!(fitDof_ == "none" && fitAnchor_ == "volume"))
        {
            auto g = [&](const scalar cc) -> scalar
            {
                davof::cellSurface T = S;
                T.c = cc;
                scalarList aF;
                fractionsOf(T, aF);
                scalar a = 0;
                forAll(cf, k) a += aF[k]*mag(Sf[cf[k]]);
                return a;
            };
            scalar cc = 0;
            bool bracketed = false;
            davof::solveOffset(g, target, -h, h, 1e-13*h, cc, bracketed);
            if (bracketed)
            {
                S.c = cc;
            }
            else
            {
                ++nFitNoBracket_;
            }
        }

        // All six coefficients (fitDof six); fitDof one keeps the one-scalar
        // anchoring, fitDof none the frame and the ring tensor only
        ++nFitCells_;
        if (fitDof_ == "six")
        {
            const davof::paraboloidFitResult res = davof::fitParaboloidFull
            (
                S, h, tF, wF, fractionsOf, fitPriorWeight_, fitTol_, fitMaxIter_
            );
            if (res.converged) ++nFitConverged_;
            if (res.stalled) ++nFitStalled_;
            if (res.singular) ++nFitSingular_;
            fitResMax_ = max(fitResMax_, res.rmax);
            fitResSum_ += res.rmax;
            if (fitResidualPtr_)
            {
                fitResidualPtr_->primitiveFieldRef()[celli] = res.rmax;
            }
        }
        else
        {
            // the face residual of the surface as it stands (the one-scalar
            // anchoring or the plane of the frame with the ring tensor)
            scalarList aF;
            fractionsOf(S, aF);
            scalar rmax = 0;
            forAll(aF, k) rmax = max(rmax, mag(aF[k] - tF[k]));
            ++nFitConverged_;
            fitResMax_ = max(fitResMax_, rmax);
            fitResSum_ += rmax;
            if (fitResidualPtr_)
            {
                fitResidualPtr_->primitiveFieldRef()[celli] = rmax;
            }
        }

        // The volume anchor (fitAnchor volume): the offset c shifted so that
        // the volume under the surface equals alpha_k V, the slopes and the
        // tensor from the face fit; the flux of the step is the volume under
        // the same surface by the same quadrature, so alpha_k stays the
        // volume under the surface and feeds back into it. Newton with the
        // wetted area as the derivative of the volume in the offset.
        if (fitAnchor_ == "volume")
        {
            const scalar Vc = V[celli];
            const scalar target = alpha1In_[celli]*Vc;
            const scalar mMag = max(mag(m[celli]), 1e-3*h*h);
            label nLeaves = 0;
            bool okAnchor = false;
            // the cell volume under the surface at the anchor's own
            // tolerance (anchorTolRel of V; the flux's tolerance is tighter
            // and the anchor's error enters the surface as an offset of
            // anchorTolRel V / A, far below the face-fraction errors)
            auto volumeOf = [&]() -> scalar
            {
                return davof::cellLiquidVolumeQuadric
                (
                    mesh_, celli, [&](const point& x) { return S.value(x); },
                    anchorTolRel_, curvedMaxDepth_, nLeaves
                );
            };
            // V(c) is non-decreasing in c (c larger: more liquid, dV/dc = the
            // wetted area, zero where the surface has left the cell). A
            // Newton guess with the wetted area, then a bracket of the root
            // (the step doubled from h/4 until the residual changes sign;
            // the volume is confined to [0, V], so a bracket exists within
            // the cell's diagonal) and the Illinois rule inside it, to the
            // tolerance or to a bracket of 1e-12 h; the secant alone failed
            // in 3 to 50 cells per step where V(c) is flat or the start far.
            // The volume residual alone is no criterion for the OFFSET: a
            // residual of tol V is an offset error of tol V / A_patch, large
            // where the wetted patch is small (a corner cut), and the plane
            // test showed exactly those faces off by 2e-4. The bracket of
            // the offset is therefore closed to 1e-11 h; the residual
            // criterion applies only where it is round-off.
            const scalar tol = max(anchorTolRel_*Vc, 1e-14*Vc);
            const scalar cTol = 1e-11*h;
            scalar c0 = S.c;
            scalar r0 = volumeOf() - target;
            // the diagnostic of the quadrature (key curvedDebug): on an exact
            // state the surface before the anchor encloses alpha V, so r0 is
            // the quadrature's own error; the worst cell is printed with its
            // pyramids against the exact linear cut of the Detrixhe-Aslam tets
            if (curvedDebug_ && mag(r0)/Vc > debugR0Max_)
            {
                debugR0Max_ = mag(r0)/Vc;
                debugCellR0_ = celli;
                debugSurface_ = S;
                debugAlpha_ = alpha1In_[celli];
            }
            if (mag(r0) <= 1e-14*Vc)
            {
                okAnchor = true;
            }
            else
            {
                // the bracket [lo, hi] with r(lo) < 0 < r(hi)
                scalar lo, hi, rlo, rhi;
                scalar step = min(mag(r0)/mMag, 0.25*h);
                step = max(step, 1e-3*h);
                if (r0 < 0)
                {
                    lo = c0; rlo = r0; hi = c0; rhi = r0;
                    for (label k = 0; k < 12 && rhi < 0; ++k)
                    {
                        hi = lo + step; S.c = hi; rhi = volumeOf() - target;
                        if (rhi < 0) { lo = hi; rlo = rhi; }
                        step *= 2;
                    }
                }
                else
                {
                    hi = c0; rhi = r0; lo = c0; rlo = r0;
                    for (label k = 0; k < 12 && rlo > 0; ++k)
                    {
                        lo = hi - step; S.c = lo; rlo = volumeOf() - target;
                        if (rlo > 0) { hi = lo; rhi = rlo; }
                        step *= 2;
                    }
                }
                if (rlo <= 0 && rhi >= 0)
                {
                    // Illinois (modified regula falsi), closed on the offset
                    label side = 0;
                    for (label it = 0; it < 80; ++it)
                    {
                        if (mag(rlo) <= 1e-14*Vc) { S.c = lo; okAnchor = true; break; }
                        if (mag(rhi) <= 1e-14*Vc) { S.c = hi; okAnchor = true; break; }
                        if (hi - lo <= cTol) { S.c = (rlo*hi - rhi*lo)/(rlo - rhi); if (!(S.c > lo && S.c < hi)) S.c = 0.5*(lo + hi); okAnchor = true; break; }
                        scalar c1 = (rlo*hi - rhi*lo)/(rlo - rhi);
                        // keep the iterate inside the bracket
                        if (!(c1 > lo && c1 < hi)) c1 = 0.5*(lo + hi);
                        S.c = c1;
                        const scalar r1 = volumeOf() - target;
                        if (mag(r1) <= 1e-14*Vc) { okAnchor = true; break; }
                        if (mag(r1) <= tol && hi - lo <= 1e-6*h) { okAnchor = true; break; }
                        if (r1 < 0)
                        {
                            lo = c1; rlo = r1;
                            if (side == -1) rhi *= 0.5;
                            side = -1;
                        }
                        else
                        {
                            hi = c1; rhi = r1;
                            if (side == 1) rlo *= 0.5;
                            side = 1;
                        }
                    }
                }
                if (!okAnchor)
                {
                    // no bracket (a cell the surface cannot fill or empty
                    // within the search) or no convergence: keep the better
                    // bracket end
                    S.c = (mag(rlo) < mag(rhi)) ? lo : hi;
                }
            }
            if (!okAnchor) ++nAnchorFail_;
            nAnchorLeaves_ += nLeaves;
            // the face fractions the anchored surface cuts against the targets
            scalarList aF;
            fractionsOf(S, aF);
            scalar dmax = 0;
            forAll(aF, k) dmax = max(dmax, mag(aF[k] - tF[k]));
            anchorShiftMax_ = max(anchorShiftMax_, dmax);
            anchorShiftSum_ += dmax;
        }

        cellSurfaces_[celli] = S;
        hasSurface_.set(celli);
    }
}


void Foam::advection::isoAdvection::curvedFlux(const scalar dt)
{
    addProfilingInFunction(geometricVoF);

    if (Pstream::parRun())
    {
        FatalErrorInFunction
            << "clippingSurface " << clippingSurface_
            << " runs in serial only (no exchange of the surfaces across"
            << " processor boundaries)" << exit(FatalError);
    }

    // The records of the face fractions, as timeIntegratedFlux() keeps them
    clearIsoFaceData();
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
    nCurvedFaces_ = nCurvedSecondRing_ = nCurvedBoundary_ = nCurvedSkipped_ = 0;
    nCurvedFallback_ = 0;
    nCurvedLeaves_ = 0;
    if (mesh_.time().timeIndex() % max(defectEvery_, label(1)) == 0 || mesh_.time().timeIndex() <= 1)
    {
        defectMax_ = 0;
        defectSum_ = 0;
        nDefectCells_ = 0;
    }

    // The surfaces of the state at t^n: the time has been advanced to
    // t^{n+1} when the advection runs (the scheme's stateTime() counts its
    // calls and cannot be asked again here)
    buildCellSurfaces(mesh_.time().value() - dt);

    interpolationCellPoint<vector> UInterp(U_);

    const labelList& own = mesh_.faceOwner();
    const labelList& nei = mesh_.faceNeighbour();
    const label nInt = mesh_.nInternalFaces();
    const labelListList& cellCells = mesh_.cellCells();
    const vectorField& Cf = mesh_.faceCentres();
    const vectorField& Sf = mesh_.faceAreas();
    const scalarField& phiIn = phi_.primitiveField();
    scalarField& dVfIn = dVf_.primitiveFieldRef();
    const DynamicField<label>& interfaceLabels = surf_->interfaceLabels();
    const polyBoundaryMesh& pbm = mesh_.boundaryMesh();

    // The candidate cells: the interface cells and their face neighbours
    // (the surface moves less than one cell per step); the interface cells
    // are the surface cells of the bounding (limitFluxes marks the cells
    // next to them from surfCells_)
    bitSet cand(mesh_.nCells());
    forAll(interfaceLabels, i)
    {
        const label celli = interfaceLabels[i];
        surfCells_.append(celli);
        cand.set(celli);
        for (const label nb : cellCells[celli])
        {
            cand.set(nb);
        }
    }

    // The surface of a face: the upwind cell's, else the other cell's, else
    // the closest surface cell among the face neighbours of both (second
    // ring); -1 when none
    auto surfaceCellOf = [&](const label up, const label other, const point& xf, bool& secondRing) -> label
    {
        secondRing = false;
        if (hasSurface_.test(up)) return up;
        if (other >= 0 && hasSurface_.test(other)) return other;
        label best = -1;
        scalar dMin = GREAT;
        auto scan = [&](const label c)
        {
            if (c < 0) return;
            for (const label nb : cellCells[c])
            {
                if (!hasSurface_.test(nb)) continue;
                const point& o = cellSurfaces_[nb].sphere ? mesh_.cellCentres()[nb] : cellSurfaces_[nb].o;
                const scalar d = magSqr(xf - o);
                if (d < dMin) { dMin = d; best = nb; }
            }
        };
        scan(up);
        scan(other);
        secondRing = (best >= 0);
        return best;
    };

    auto processFace = [&](const label facei, const label up, const label other, const scalar phiF, const bool boundary)
    {
        bool secondRing = false;
        const label sc = surfaceCellOf(up, other, Cf[facei], secondRing);
        if (sc < 0)
        {
            ++nCurvedSkipped_;
            return;
        }
        const davof::cellSurface& S = cellSurfaces_[sc];
        auto phi = [&](const point& x) { return S.value(x); };

        // The displacement of the step at the face
        const vector d = UInterp.interpolate(Cf[facei], up)*dt;

        // The volume flux: the liquid volume of the prism between f - d and f
        if (!(boundary && phiF < 0))
        {
            scalar vLiq = 0;
            if (mag(d & Sf[facei]) > VSMALL)
            {
                vLiq = davof::faceSweptLiquidVolume
                (
                    mesh_, facei, d, phi, curvedTolRel_, curvedMaxDepth_, nCurvedLeaves_
                );
            }
            const scalar dVfNew = (phiF >= 0 ? vLiq : -vLiq);
            if (boundary)
            {
                setFaceValue(dVf_, facei, dVfNew);
            }
            else
            {
                dVfIn[facei] = dVfNew;
            }
        }

        // The face fraction: the pre-image's liquid fraction under the surface
        label nfb = 0;
        scalar af = davof::faceCutDisplaced(mesh_, facei, -d, phi, nfb);
        nCurvedFallback_ += nfb;
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
        ++nCurvedFaces_;
        if (secondRing) ++nCurvedSecondRing_;
        if (boundary) ++nCurvedBoundary_;
    };

    for (label facei = 0; facei < nInt; ++facei)
    {
        const label o = own[facei];
        const label n = nei[facei];
        if (!cand.test(o) && !cand.test(n)) continue;
        const scalar phiF = phiIn[facei];
        const label up = (phiF >= 0) ? o : n;
        const label other = (phiF >= 0) ? n : o;
        processFace(facei, up, other, phiF, false);
    }
    forAll(pbm, patchi)
    {
        const polyPatch& pp = pbm[patchi];
        if (alphafPtr_->boundaryField()[patchi].size() != pp.size()) continue;
        if (pp.coupled()) continue;
        forAll(pp, i)
        {
            const label facei = pp.start() + i;
            const label o = own[facei];
            if (!cand.test(o)) continue;
            processFace(facei, o, -1, faceValue(phi_, facei), true);
        }
    }

    // The paraboloid defect: alpha_k against the volume under the cell's
    // surface, at the start of the step (diagnostic, every defectEvery-th
    // step and at the first; the per-step line prints the last evaluation)
    if
    (
        defectLevels_ >= 0
     && (mesh_.time().timeIndex() % max(defectEvery_, label(1)) == 0 || mesh_.time().timeIndex() <= 1)
    )
    {
        if (writeCurvedFields_ && !paraboloidDefectPtr_)
        {
            paraboloidDefectPtr_.reset
            (
                new volScalarField
                (
                    IOobject
                    (
                        IOobject::groupName("paraboloidDefectDavof", alpha1_.group()),
                        mesh_.time().timeName(),
                        mesh_,
                        IOobject::NO_READ,
                        IOobject::AUTO_WRITE
                    ),
                    mesh_,
                    dimensionedScalar(dimless, Zero)
                )
            );
        }
        if (paraboloidDefectPtr_)
        {
            paraboloidDefectPtr_->primitiveFieldRef() = 0;
        }
        const scalarField& V = mesh_.V();
        forAll(interfaceLabels, i)
        {
            const label celli = interfaceLabels[i];
            if (!hasSurface_.test(celli)) continue;
            const davof::cellSurface& S = cellSurfaces_[celli];
            // the same quadrature as the flux (defectLevels kept as the switch)
            label nLeavesDefect = 0;
            const scalar vLiq = davof::cellLiquidVolumeQuadric
            (
                mesh_, celli, [&](const point& x) { return S.value(x); },
                curvedTolRel_, curvedMaxDepth_, nLeavesDefect
            );
            const scalar defect = alpha1In_[celli] - vLiq/V[celli];
            defectMax_ = max(defectMax_, mag(defect));
            defectSum_ += mag(defect);
            ++nDefectCells_;
            if (paraboloidDefectPtr_)
            {
                paraboloidDefectPtr_->primitiveFieldRef()[celli] = defect;
            }
        }
    }
}


void Foam::advection::isoAdvection::reportCurved() const
{
    Info<< "isoAdvection: curved clipping (" << clippingSurface_
        << "): alphaf and dVf set on " << nCurvedFaces_ << " faces ("
        << nCurvedSecondRing_ << " second ring, " << nCurvedBoundary_
        << " boundary), " << nCurvedSkipped_ << " skipped, "
        << nAlphafClipped_ << " clips, " << nCurvedFallback_
        << " linear cuts, " << nCurvedLeaves_ << " quadrature leaves";
    if (clippingSurface_ == "paraboloid")
    {
        Info<< "; fit: cells " << nFitCells_ << ", converged " << nFitConverged_
            << ", stalled " << nFitStalled_ << ", singular " << nFitSingular_
            << ", no bracket " << nFitNoBracket_ << ", no frame " << nFitNoFrame_
            << ", residual max " << fitResMax_
            << " mean " << fitResSum_/max(nFitCells_, label(1));
        if (fitAnchor_ == "volume")
        {
            Info<< "; volume anchor: failed " << nAnchorFail_
                << ", face shift max " << anchorShiftMax_
                << " mean " << anchorShiftSum_/max(nFitCells_, label(1))
                << ", leaves " << nAnchorLeaves_;
        }
    }
    if (nDefectCells_ > 0)
    {
        Info<< "; surface defect max " << defectMax_
            << " L1 " << defectSum_/nDefectCells_
            << " over " << nDefectCells_ << " cells";
    }
    Info<< endl;

    if (curvedDebug_ && debugCellR0_ >= 0)
    {
        // the worst cell's pyramids: the split quadrature against the linear
        // cut of the Detrixhe-Aslam tets (exact for a plane)
        const label celli = debugCellR0_;
        const davof::cellSurface& S = debugSurface_;
        auto phi = [&](const point& x) { return S.value(x); };
        const pointField& P = mesh_.points();
        const faceList& faces = mesh_.faces();
        const vectorField& Cf = mesh_.faceCentres();
        const point& xc = mesh_.cellCentres()[celli];
        const labelList& cFaces = mesh_.cells()[celli];
        const scalar Vc = mesh_.V()[celli];
        Info<< "curved debug: cell " << celli << " alpha " << debugAlpha_
            << " |r0|/V " << debugR0Max_ << " centre " << xc
            << " surface o " << S.o << " e3 " << S.e3 << " c " << S.c
            << " b " << S.b11 << " " << S.b12 << " " << S.b22 << " phi(xc) " << phi(xc) << nl;
        scalar vq = 0, vl = 0;
        forAll(cFaces, cf)
        {
            const label faceI = cFaces[cf];
            const face& f = faces[faceI];
            const point& xf = Cf[faceI];
            forAll(f, ip)
            {
                const label ip1 = f.nextLabel(ip);
                const point& p0 = P[f[ip]];
                const point& p1 = P[ip1];
                label nl2 = 0;
                const scalar vQuad = davof::pyramidLiquidVolume
                (
                    xc, xf, p0, p1, phi, anchorTolRel_*Vc, curvedMaxDepth_, nl2
                );
                const scalar vol = mag(((xf - xc) ^ (p0 - xc)) & (p1 - xc))/6.0;
                const scalar vLin = davof::tetNegativeFraction(phi(xc), phi(xf), phi(p0), phi(p1))*vol;
                vq += vQuad; vl += vLin;
                if (mag(vQuad - vLin) > 1e-12*Vc)
                {
                    Info<< "  face " << faceI << " tri " << ip << " quad " << vQuad/Vc
                        << " linear " << vLin/Vc << " diff " << (vQuad - vLin)/Vc
                        << " phi(xf,p0,p1) " << phi(xf) << " " << phi(p0) << " " << phi(p1)
                        << " leaves " << nl2 << nl;
                }
            }
        }
        Info<< "  total quad " << vq/Vc << " linear " << vl/Vc << " alpha " << debugAlpha_ << endl;
    }
}


// ************************************************************************* //
