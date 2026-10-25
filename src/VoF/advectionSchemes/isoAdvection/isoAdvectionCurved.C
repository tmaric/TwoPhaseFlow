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
    plane sweep of timeIntegratedFlux() when active. In parallel (2026-10-24)
    the surfaces of the cells across the processor patches are exchanged per
    patch face; a processor face is computed on one side (the upwind side
    when its cell has a surface, else the downwind side when that cell has
    one, else the upwind side from its second ring) and copied to the other
    by syncProcPatches as the plane sweep's fluxes are; the second-ring
    search of a cell next to a processor patch sees its remote face
    neighbours; the counters of the Info line are reduced.

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
#include "processorPolyPatch.H"
#include "PstreamBuffers.H"

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace
{

// The surface of a cell as 21 scalars (the two flags, the origin, the
// frame, the seven coefficients), for the exchange across processor patches
inline void packSurface(const Foam::davof::cellSurface& S, Foam::UList<Foam::scalar>& v, const Foam::label off)
{
    v[off] = S.paraboloid ? 1 : 0;
    v[off + 1] = S.sphere ? 1 : 0;
    for (Foam::direction k = 0; k < 3; ++k)
    {
        v[off + 2 + k] = S.o[k];
        v[off + 5 + k] = S.e1[k];
        v[off + 8 + k] = S.e2[k];
        v[off + 11 + k] = S.e3[k];
    }
    v[off + 14] = S.b11; v[off + 15] = S.b12; v[off + 16] = S.b22;
    v[off + 17] = S.c; v[off + 18] = S.a1; v[off + 19] = S.a2; v[off + 20] = S.R;
}

inline Foam::davof::cellSurface unpackSurface(const Foam::UList<Foam::scalar>& v, const Foam::label off)
{
    Foam::davof::cellSurface S;
    S.paraboloid = (v[off] > 0.5);
    S.sphere = (v[off + 1] > 0.5);
    for (Foam::direction k = 0; k < 3; ++k)
    {
        S.o[k] = v[off + 2 + k];
        S.e1[k] = v[off + 5 + k];
        S.e2[k] = v[off + 8 + k];
        S.e3[k] = v[off + 11 + k];
    }
    S.b11 = v[off + 14]; S.b12 = v[off + 15]; S.b22 = v[off + 16];
    S.c = v[off + 17]; S.a1 = v[off + 18]; S.a2 = v[off + 19]; S.R = v[off + 20];
    return S;
}

constexpr Foam::label nPack = 21;

} // End anonymous namespace


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
        // (the volume rule's face cut, consistent with its volumes)
        auto fractionsOf = [&](const davof::cellSurface& T, scalarList& aF)
        {
            aF.setSize(F);
            label nfb = 0;
            forAll(cf, k)
            {
                aF[k] = volumeRulePtr_->faceLiquidFraction(cf[k], vector::zero, T, nfb);
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
            ++nFitConverged_;
            // the face residual of the surface as it stands (the one-scalar
            // anchoring or the plane of the frame with the ring tensor): a
            // diagnostic of one face-fraction evaluation per face, computed
            // only when it is printed (anchorDiagnostics) or written
            if (anchorDiagnostics_ || fitResidualPtr_)
            {
                scalarList aF;
                fractionsOf(S, aF);
                scalar rmax = 0;
                forAll(aF, k) rmax = max(rmax, mag(aF[k] - tF[k]));
                fitResMax_ = max(fitResMax_, rmax);
                fitResSum_ += rmax;
                if (fitResidualPtr_)
                {
                    fitResidualPtr_->primitiveFieldRef()[celli] = rmax;
                }
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
            bool okAnchor = false;
            // the cell volume under the surface by the volume rule at the
            // anchor's tolerance (anchorTolRel of V), the same rule as the
            // flux, so that the anchored volume is the flux's volume
            auto volumeOf = [&]() -> scalar
            {
                return volumeRulePtr_->cellLiquidVolume(celli, S, anchorTolRel_);
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
            // the face fractions the anchored surface cuts against the
            // targets (the face shift): a diagnostic, computed when printed
            if (anchorDiagnostics_)
            {
                scalarList aF;
                fractionsOf(S, aF);
                scalar dmax = 0;
                forAll(aF, k) dmax = max(dmax, mag(aF[k] - tF[k]));
                anchorShiftMax_ = max(anchorShiftMax_, dmax);
                anchorShiftSum_ += dmax;
            }
        }

        cellSurfaces_[celli] = S;
        hasSurface_.set(celli);
    }
    // the rule's cost and fallbacks of the fit and the anchor
    nAnchorLeaves_ = volumeRulePtr_->nLeaves();
    nRuleFallback_ += volumeRulePtr_->nFallback();
    volumeRulePtr_->resetCounters();
}


void Foam::advection::isoAdvection::curvedFlux(const scalar dt)
{
    addProfilingInFunction(geometricVoF);

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
    nRuleFallback_ = 0;
    volumeRulePtr_->resetCounters();
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
    const label nPatches = pbm.size();

    // The surfaces of the cells across the processor patches (the face
    // neighbours of the local cells there), one per patch face, exchanged
    // with the neighbouring processor: a processor face takes the surface
    // of its upwind cell as an internal face does, and the second-ring
    // search of a cell next to the boundary sees its remote face neighbours
    List<List<davof::cellSurface>> nbrSurface(nPatches);
    List<boolList> nbrHas(nPatches);
    if (Pstream::parRun())
    {
        PstreamBuffers pBufs(Pstream::commsTypes::nonBlocking);
        forAll(pbm, patchi)
        {
            const polyPatch& pp = pbm[patchi];
            if (!isA<processorPolyPatch>(pp)) continue;
            const processorPolyPatch& procPatch = refCast<const processorPolyPatch>(pp);
            const labelUList& fc = pp.faceCells();
            scalarList packed(nPack*fc.size(), Zero);
            boolList has(fc.size(), false);
            forAll(fc, i)
            {
                if (hasSurface_.test(fc[i]))
                {
                    has[i] = true;
                    packSurface(cellSurfaces_[fc[i]], packed, nPack*i);
                }
            }
            UOPstream toNbr(procPatch.neighbProcNo(), pBufs);
            toNbr << has << packed;
        }
        pBufs.finishedSends();
        forAll(pbm, patchi)
        {
            const polyPatch& pp = pbm[patchi];
            if (!isA<processorPolyPatch>(pp)) continue;
            const processorPolyPatch& procPatch = refCast<const processorPolyPatch>(pp);
            UIPstream fromNbr(procPatch.neighbProcNo(), pBufs);
            scalarList packed;
            fromNbr >> nbrHas[patchi] >> packed;
            nbrSurface[patchi].setSize(pp.size());
            forAll(pp, i)
            {
                if (nbrHas[patchi][i])
                {
                    nbrSurface[patchi][i] = unpackSurface(packed, nPack*i);
                }
            }
        }
    }

    // The candidate cells: the interface cells and their face neighbours
    // (the surface moves less than one cell per step), the face neighbours
    // across the processor patches included; the interface cells are the
    // surface cells of the bounding (limitFluxes marks the cells next to
    // them from surfCells_)
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
    forAll(pbm, patchi)
    {
        if (nbrHas[patchi].empty()) continue;
        const labelUList& fc = pbm[patchi].faceCells();
        forAll(fc, i)
        {
            if (nbrHas[patchi][i]) cand.set(fc[i]);
        }
    }

    // The second exchange: per processor face the candidate flag of the
    // cell across and the closest surface among that cell's face
    // neighbours with its distance to the face centre, so that a processor
    // face is processed when either cell is a candidate and its second ring
    // is the serial one (the face neighbours of both cells)
    List<boolList> nbrCand(nPatches), nbrRingHas(nPatches);
    List<scalarList> nbrRingDist(nPatches);
    List<List<davof::cellSurface>> nbrRing(nPatches);
    if (Pstream::parRun())
    {
        PstreamBuffers pBufs(Pstream::commsTypes::nonBlocking);
        forAll(pbm, patchi)
        {
            const polyPatch& pp = pbm[patchi];
            if (!isA<processorPolyPatch>(pp)) continue;
            const processorPolyPatch& procPatch = refCast<const processorPolyPatch>(pp);
            const labelUList& fc = pp.faceCells();
            boolList isCand(fc.size(), false), ringHas(fc.size(), false);
            scalarList ringDist(fc.size(), GREAT);
            scalarList packed(nPack*fc.size(), Zero);
            forAll(fc, i)
            {
                const label c = fc[i];
                isCand[i] = cand.test(c);
                const point& xf = Cf[pp.start() + i];
                label best = -1;
                scalar dMin = GREAT;
                for (const label nb : cellCells[c])
                {
                    if (!hasSurface_.test(nb)) continue;
                    const point& o = cellSurfaces_[nb].sphere ? mesh_.cellCentres()[nb] : cellSurfaces_[nb].o;
                    const scalar dd = magSqr(xf - o);
                    // the same tie-break as the face's search below
                    bool take = (best < 0 || dd < dMin*(1 - 1e-9));
                    if (!take && dd <= dMin*(1 + 1e-9))
                    {
                        const point& ob = cellSurfaces_[best].o;
                        const point& on = cellSurfaces_[nb].o;
                        for (direction k = 0; k < 3; ++k)
                        {
                            if (on[k] < ob[k]) { take = true; break; }
                            if (on[k] > ob[k]) break;
                        }
                    }
                    if (take) { dMin = min(dMin, dd); best = nb; }
                }
                if (best >= 0)
                {
                    ringHas[i] = true;
                    ringDist[i] = dMin;
                    packSurface(cellSurfaces_[best], packed, nPack*i);
                }
            }
            UOPstream toNbr(procPatch.neighbProcNo(), pBufs);
            toNbr << isCand << ringHas << ringDist << packed;
        }
        pBufs.finishedSends();
        forAll(pbm, patchi)
        {
            const polyPatch& pp = pbm[patchi];
            if (!isA<processorPolyPatch>(pp)) continue;
            const processorPolyPatch& procPatch = refCast<const processorPolyPatch>(pp);
            UIPstream fromNbr(procPatch.neighbProcNo(), pBufs);
            scalarList packed;
            fromNbr >> nbrCand[patchi] >> nbrRingHas[patchi] >> nbrRingDist[patchi] >> packed;
            nbrRing[patchi].setSize(pp.size());
            forAll(pp, i)
            {
                if (nbrRingHas[patchi][i])
                {
                    nbrRing[patchi][i] = unpackSurface(packed, nPack*i);
                }
            }
        }
    }

    // A zero flux (a face parallel to the velocity) leaves "upwind"
    // undefined; the cell with the lexicographically smaller centre takes
    // the face, a convention both sides of a processor patch can evaluate
    // (the serial owner cell, the lower cell index, is not the patch's
    // owning side after decomposition: the four-processor check differed by
    // 1e-3 of L1(alphaf) through the other surface on those faces)
    auto lexLessPoint = [](const point& a, const point& b) -> bool
    {
        for (direction k = 0; k < 3; ++k)
        {
            if (a[k] < b[k]) return true;
            if (a[k] > b[k]) return false;
        }
        return false;
    };
    const vectorField& C = mesh_.cellCentres();

    // The surface of a face: the upwind cell's, else the other cell's, else
    // the closest surface among the face neighbours of both (second ring,
    // the remote face neighbours across the processor patches included);
    // nullptr when none. For a processor face the remote cell's surface
    // (remote, remoteHas) and the closest surface of its face neighbours
    // with its distance (remoteRing, remoteRingHas, remoteRingDist) come
    // from the exchanges.
    auto surfaceOf = [&]
    (
        const label up, const label other,
        const davof::cellSurface* remote, const bool remoteHas,
        const davof::cellSurface* remoteRing, const bool remoteRingHas, const scalar remoteRingDist,
        const point& xf, bool& secondRing
    ) -> const davof::cellSurface*
    {
        secondRing = false;
        if (up >= 0 && hasSurface_.test(up)) return &cellSurfaces_[up];
        if (other >= 0 && hasSurface_.test(other)) return &cellSurfaces_[other];
        if (remoteHas) return remote;
        const davof::cellSurface* best = nullptr;
        scalar dMin = GREAT;
        // the closest surface; ties (two surfaces at the same distance,
        // common on a Cartesian mesh with a symmetric interface) are broken
        // by the surface origin lexicographically, so that the choice does
        // not depend on the order of consideration, which differs between
        // the serial and the parallel scan and amplified round-off to 1e-3
        // of L1(alphaf) in the four-processor check
        auto lexLess = [](const point& a, const point& b) -> bool
        {
            for (direction k = 0; k < 3; ++k)
            {
                if (a[k] < b[k]) return true;
                if (a[k] > b[k]) return false;
            }
            return false;
        };
        auto considerD = [&](const davof::cellSurface& T, const scalar d)
        {
            if (!best || d < dMin*(1 - 1e-9))
            {
                dMin = d; best = &T;
            }
            else if (d <= dMin*(1 + 1e-9) && lexLess(T.o, best->o))
            {
                dMin = min(dMin, d); best = &T;
            }
        };
        auto consider = [&](const davof::cellSurface& T, const point& o)
        {
            considerD(T, magSqr(xf - o));
        };
        if (remoteRingHas) considerD(*remoteRing, remoteRingDist);
        auto scan = [&](const label c)
        {
            if (c < 0) return;
            for (const label nb : cellCells[c])
            {
                if (!hasSurface_.test(nb)) continue;
                consider(cellSurfaces_[nb], cellSurfaces_[nb].sphere ? mesh_.cellCentres()[nb] : cellSurfaces_[nb].o);
            }
            // the remote face neighbours of c across the processor patches
            if (Pstream::parRun())
            {
                for (const label cf : mesh_.cells()[c])
                {
                    if (cf < nInt) continue;
                    const label patchi = pbm.whichPatch(cf);
                    if (patchi < 0 || nbrHas[patchi].empty()) continue;
                    const label pf = cf - pbm[patchi].start();
                    if (!nbrHas[patchi][pf]) continue;
                    const davof::cellSurface& T = nbrSurface[patchi][pf];
                    consider(T, T.o);
                }
            }
        };
        scan(up);
        scan(other);
        secondRing = (best != nullptr);
        return best;
    };

    // facei: the face; local: the local cell whose velocity interpolates the
    // displacement (the upwind cell, or the local cell of a processor face);
    // up/other: the local upwind and downwind cells (-1 when remote or
    // none); remote/remoteHas: the surface across a processor face
    auto processFace = [&]
    (
        const label facei, const label local, const label up, const label other,
        const davof::cellSurface* remote, const bool remoteHas,
        const davof::cellSurface* remoteRing, const bool remoteRingHas, const scalar remoteRingDist,
        const scalar phiF, const bool boundary, const bool procFace
    )
    {
        bool secondRing = false;
        const davof::cellSurface* Sp = surfaceOf
        (
            up, other, remote, remoteHas, remoteRing, remoteRingHas, remoteRingDist,
            Cf[facei], secondRing
        );
        if (!Sp)
        {
            ++nCurvedSkipped_;
            return;
        }
        const davof::cellSurface& S = *Sp;

        // The displacement of the step at the face
        const vector d = UInterp.interpolate(Cf[facei], local)*dt;

        // The volume flux: the liquid volume of the prism between f - d and f
        if (!(boundary && phiF < 0))
        {
            // a face parallel to the displacement sweeps no volume: the
            // test is relative (the round-off product of 1e-28 of an
            // axis-aligned face and displacement built a flat prism whose
            // face sum was 5e-14, five percent of the cell, 2026-10-24; the
            // zero components of the face area vectors carry 1e-14 of
            // round-off, so the threshold is 1e-10: a sweep below 1e-10 of
            // a normal flux is taken as none)
            scalar vLiq = 0;
            if (mag(d & Sf[facei]) > 1e-10*mag(d)*mag(Sf[facei]))
            {
                vLiq = volumeRulePtr_->prismLiquidVolume(facei, d, S, curvedTolRel_);
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
        scalar af = volumeRulePtr_->faceLiquidFraction(facei, -d, S, nfb);
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
        if (procFace) checkIfOnProcPatch(facei);
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
        const bool ownerUp = (phiF > 0) || (phiF == 0 && lexLessPoint(C[o], C[n]));
        const label up = ownerUp ? o : n;
        const label other = ownerUp ? n : o;
        processFace(facei, up, up, other, nullptr, false, nullptr, false, GREAT, phiF, false, false);
    }
    forAll(pbm, patchi)
    {
        const polyPatch& pp = pbm[patchi];
        if (alphafPtr_->boundaryField()[patchi].size() != pp.size()) continue;
        if (pp.coupled() && !isA<processorPolyPatch>(pp)) continue;
        const bool procFace = isA<processorPolyPatch>(pp);
        // the cell centres across a processor patch (the zero-flux convention)
        tmp<vectorField> tCnbr;
        if (procFace)
        {
            tCnbr = mesh_.C().boundaryField()[patchi].patchNeighbourField();
        }
        forAll(pp, i)
        {
            const label facei = pp.start() + i;
            const label o = own[facei];
            const scalar phiF = faceValue(phi_, facei);
            if (!procFace)
            {
                if (!cand.test(o)) continue;
                processFace(facei, o, o, -1, nullptr, false, nullptr, false, GREAT, phiF, true, false);
                continue;
            }
            // A processor face is computed on ONE side and copied to the
            // other by syncProcPatches, as the plane sweep does: on the
            // upwind side when the upwind cell has a surface, else on the
            // downwind side when that cell has one, else on the upwind side
            // from its second ring; both sides know both flags
            const bool localHas = hasSurface_.test(o);
            const bool remoteHas = !nbrHas[patchi].empty() && nbrHas[patchi][i];
            // a zero flux (a face parallel to the velocity) makes both sides
            // "upwind" by phiF >= 0; the cell with the lexicographically
            // smaller centre takes it, as in the serial loop above, and the
            // two sides never both compute
            const bool localUpwind =
                (phiF > 0)
             || (phiF == 0 && lexLessPoint(C[o], tCnbr()[i]));
            // a face with no surface on either side is processed when either
            // cell is a candidate (the serial rule), by the upwind side
            const bool localCand = cand.test(o);
            const bool remoteCand = !nbrCand[patchi].empty() && nbrCand[patchi][i];
            bool computeHere;
            if (localUpwind) computeHere = localHas || (!remoteHas && (localCand || remoteCand));
            else computeHere = localHas && !remoteHas;
            if (!computeHere) continue;
            const davof::cellSurface* remote = remoteHas ? &nbrSurface[patchi][i] : nullptr;
            const bool remoteRingHas = !nbrRingHas[patchi].empty() && nbrRingHas[patchi][i];
            const davof::cellSurface* remoteRing = remoteRingHas ? &nbrRing[patchi][i] : nullptr;
            const scalar remoteRingDist = remoteRingHas ? nbrRingDist[patchi][i] : GREAT;
            processFace
            (
                facei, o, o, -1, remote, remoteHas, remoteRing, remoteRingHas, remoteRingDist,
                phiF, true, true
            );
        }
    }
    if (Pstream::parRun())
    {
        // the computed processor faces to the other side (the flux with the
        // opposite sign, the face fraction as it is)
        syncProcPatches(dVf_, phi_, false, alphafPtr_);
    }
    nCurvedLeaves_ = volumeRulePtr_->nLeaves();
    nRuleFallback_ += volumeRulePtr_->nFallback();
    volumeRulePtr_->resetCounters();

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
            // the same rule as the flux (defectLevels kept as the switch)
            const scalar vLiq = volumeRulePtr_->cellLiquidVolume(celli, S, curvedTolRel_);
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
    // the counters over the processors (sums; the maxima as maxima)
    label nFaces = nCurvedFaces_, nSecond = nCurvedSecondRing_, nBound = nCurvedBoundary_;
    label nSkip = nCurvedSkipped_, nClip = nAlphafClipped_, nLin = nCurvedFallback_;
    label nLeaves = nCurvedLeaves_, nFitC = nFitCells_, nFitConv = nFitConverged_;
    label nFitSt = nFitStalled_, nFitSing = nFitSingular_, nFitNoBr = nFitNoBracket_;
    label nFitNoFr = nFitNoFrame_, nAnchF = nAnchorFail_, nAnchL = nAnchorLeaves_;
    label nDef = nDefectCells_, nRuleFb = nRuleFallback_ + volumeRulePtr_->nFallback();
    scalar fitRMax = fitResMax_, fitRSum = fitResSum_, ashMax = anchorShiftMax_;
    scalar ashSum = anchorShiftSum_, defMax = defectMax_, defSum = defectSum_;
    if (Pstream::parRun())
    {
        reduce(nFaces, sumOp<label>()); reduce(nSecond, sumOp<label>());
        reduce(nBound, sumOp<label>()); reduce(nSkip, sumOp<label>());
        reduce(nClip, sumOp<label>()); reduce(nLin, sumOp<label>());
        reduce(nLeaves, sumOp<label>()); reduce(nFitC, sumOp<label>());
        reduce(nFitConv, sumOp<label>()); reduce(nFitSt, sumOp<label>());
        reduce(nFitSing, sumOp<label>()); reduce(nFitNoBr, sumOp<label>());
        reduce(nFitNoFr, sumOp<label>()); reduce(nAnchF, sumOp<label>());
        reduce(nAnchL, sumOp<label>()); reduce(nDef, sumOp<label>());
        reduce(nRuleFb, sumOp<label>());
        reduce(fitRMax, maxOp<scalar>()); reduce(fitRSum, sumOp<scalar>());
        reduce(ashMax, maxOp<scalar>()); reduce(ashSum, sumOp<scalar>());
        reduce(defMax, maxOp<scalar>()); reduce(defSum, sumOp<scalar>());
    }
    Info<< "isoAdvection: curved clipping (" << clippingSurface_
        << "): alphaf and dVf set on " << nFaces << " faces ("
        << nSecond << " second ring, " << nBound
        << " boundary), " << nSkip << " skipped, "
        << nClip << " clips, " << nLin
        << " linear cuts, " << nLeaves << " quadrature leaves";
    if (clippingSurface_ == "paraboloid")
    {
        Info<< "; fit: cells " << nFitC << ", converged " << nFitConv
            << ", stalled " << nFitSt << ", singular " << nFitSing
            << ", no bracket " << nFitNoBr << ", no frame " << nFitNoFr;
        // the face residual is the fit's own with fitDof six, else a
        // diagnostic computed only when anchorDiagnostics (or the field) is on
        if (fitDof_ == "six" || anchorDiagnostics_ || fitResidualPtr_)
        {
            Info<< ", residual max " << fitRMax
                << " mean " << fitRSum/max(nFitC, label(1));
        }
        if (fitAnchor_ == "volume")
        {
            Info<< "; volume anchor: failed " << nAnchF;
            if (anchorDiagnostics_)
            {
                Info<< ", face shift max " << ashMax
                    << " mean " << ashSum/max(nFitC, label(1));
            }
            Info<< ", leaves " << nAnchL;
        }
    }
    if (nDef > 0)
    {
        Info<< "; surface defect max " << defMax
            << " L1 " << defSum/nDef
            << " over " << nDef << " cells";
    }
    Info<< "; volume rule " << volumeRulePtr_->type()
        << " fallbacks " << nRuleFb;
    Info<< endl;
    volumeRulePtr_->report();

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
        // the selected rule's value against the split quadrature and the linear cut
        const scalar vRule = volumeRulePtr_->cellLiquidVolume(celli, S, anchorTolRel_);
        Info<< "  total quad " << vq/Vc << " linear " << vl/Vc << " alpha " << debugAlpha_
            << " rule(" << volumeRulePtr_->type() << ") " << vRule/Vc << endl;

        // The reference of the prism quadrature on a sphere: the split
        // quadrature against the exact sphere-tetrahedron overlap on every
        // prism of the interface cells' faces, with the sphere of the
        // exactNormal dictionary at the state's time and the uniform
        // displacement U0 dt; the faces whose normal is perpendicular to U0
        // (the sphere's equator) carry the silhouette kink the split rule
        // leaves to refinement (2026-10-24: the Kromer-Bothe rule differed
        // from the split quadrature by 5e-6 of the prism volume there)
        if (modelDict().found("exactNormal") && !Pstream::parRun())
        {
            const dictionary& ed = modelDict().subDict("exactNormal");
            if (ed.found("radius"))
            {
                const scalar dt = mesh_.time().deltaTValue();
                davof::cellSurface Ssph;
                Ssph.sphere = true;
                Ssph.o = ed.get<point>("origin") + ed.getOrDefault<vector>("U0", vector::zero)*(mesh_.time().value() - dt);
                Ssph.R = ed.get<scalar>("radius");
                const vector d = ed.getOrDefault<vector>("U0", vector::zero)*dt;
                dictionary dSplit(modelDict()); dSplit.set("volumeRule", "splitQuadrature");
                dictionary dExact(modelDict()); dExact.set("volumeRule", "exactSphere");
                autoPtr<davofVolumeRule> split(davofVolumeRule::New(mesh_, dSplit));
                autoPtr<davofVolumeRule> exact(davofVolumeRule::New(mesh_, dExact));
                const vectorField& Sf = mesh_.faceAreas();
                const DynamicField<label>& labels = surf_->interfaceLabels();
                bitSet done(mesh_.nFaces());
                scalar maxRel = 0, sumRel = 0, maxRelEq = 0;
                label n = 0, nEq = 0, worst = -1;
                forAll(labels, i)
                {
                    for (const label f : mesh_.cells()[labels[i]])
                    {
                        if (done.test(f)) continue;
                        done.set(f);
                        const scalar pv = mag(d & Sf[f]);
                        if (pv <= 1e-10*mag(d)*mag(Sf[f])) continue;
                        const scalar vs = split->prismLiquidVolume(f, d, Ssph, curvedTolRel_);
                        const scalar ve = exact->prismLiquidVolume(f, d, Ssph, curvedTolRel_);
                        const scalar rel = mag(vs - ve)/pv;
                        sumRel += rel; ++n;
                        if (rel > maxRel) { maxRel = rel; worst = f; }
                        // the equatorial faces: the sphere's radius through
                        // the face centre perpendicular to d within 0.1 h
                        const point& xf = mesh_.faceCentres()[f];
                        if (mag((xf - Ssph.o) & d) < 0.1*Foam::cbrt(mesh_.V()[labels[i]])*mag(d))
                        {
                            maxRelEq = max(maxRelEq, rel); ++nEq;
                        }
                    }
                }
                Info<< "curved debug sphere prisms: |split - exact|/prismVol over " << n
                    << " prisms max " << maxRel << " (face " << worst << ") mean " << sumRel/max(n, label(1))
                    << "; over " << nEq << " equatorial prisms max " << maxRelEq << endl;
            }
        }
    }
}


// ************************************************************************* //
