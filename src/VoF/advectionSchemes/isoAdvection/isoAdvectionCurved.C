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
    whose fit does not converge keep the best iterate (counted). The volume
    anchor (fitAnchor volume) positions the surface by its offset so that
    the volume under it is alpha V, by the Consecutive Cubic Spline
    positioning of T. Maric, Iterative volume-of-fluid interface positioning
    in general polyhedrons with Consecutive Cubic Spline interpolation,
    J. Comput. Phys.: X 11 (2021) 100093, doi:10.1016/j.jcpx.2021.100093,
    arXiv:2012.02227, applied to the paraboloid with the exact volumes and
    slopes of the volume rule (anchorMethod ccs, the default; newton and
    illinois selectable), and takes the edge of the flat root interval for
    a cell whose target is 0 or V.

SourceFiles
    isoAdvectionCurved.C

\*---------------------------------------------------------------------------*/

#include "isoAdvection.H"
#include "davofCurvedGeometry.H"
#include "davofKromerBothe.H"
#include "interpolationCellPoint.H"
#include "volFields.H"
#include "processorPolyPatch.H"
#include "PstreamBuffers.H"
#include <cmath>

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

// The root in (0, 1) of the cubic Hermite interpolant p(s) with p(0) = f0 < 0,
// p(1) = f1 > 0 and the derivatives p'(0) = d0, p'(1) = d1 (in s): Newton on
// the polynomial from the secant point, bisection where it leaves the
// polynomial's own bracket. The volume anchor's Consecutive Cubic Spline step
// (T. Maric, Iterative volume-of-fluid interface positioning in general
// polyhedrons with Consecutive Cubic Spline interpolation, J. Comput. Phys.:
// X 11 (2021) 100093, doi:10.1016/j.jcpx.2021.100093, arXiv:2012.02227).
inline bool hermiteRoot01
(
    const Foam::scalar f0, const Foam::scalar d0,
    const Foam::scalar f1, const Foam::scalar d1,
    Foam::scalar& s
)
{
    using Foam::scalar; using Foam::label; using Foam::mag;
    auto p = [&](const scalar t, scalar& v, scalar& dv)
    {
        const scalar t2 = t*t, t3 = t2*t;
        v = (2*t3 - 3*t2 + 1)*f0 + (t3 - 2*t2 + t)*d0 + (-2*t3 + 3*t2)*f1 + (t3 - t2)*d1;
        dv = (6*t2 - 6*t)*f0 + (3*t2 - 4*t + 1)*d0 + (-6*t2 + 6*t)*f1 + (3*t2 - 2*t)*d1;
    };
    scalar a = 0, b = 1;
    s = f0/(f0 - f1);
    for (label it = 0; it < 60; ++it)
    {
        scalar v, dv;
        p(s, v, dv);
        if (v < 0) a = s; else b = s;
        if (mag(v) <= 1e-15*(mag(f0) + mag(f1)) || b - a <= 1e-15) break;
        scalar sn = (dv > 0) ? s - v/dv : 0.5*(a + b);
        if (!(sn > a && sn < b)) sn = 0.5*(a + b);
        s = sn;
    }
    return (s > 0 && s < 1);
}

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
    nAnchorEvals_ = 0;
    nAnchorEdge_ = 0;
    nCompare_ = 0;
    nCompareDiff_ = 0;
    nCompareDiffEdge_ = 0;
    nCompareDiffResolved_ = 0;
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

    // paraboloid: the ring curvature tensor of the start-of-step state,
    // fitted through the cut centroids or (tensorCentroid anchored) through
    // the centroids re-anchored by the volume rule
    nCentroidFail_ = 0;
    if (tensorCentroid_ == "anchored")
    {
        anchorTensorCentroids();
    }
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
        // the pre-anchor surface and its face residual, for the event print
        // after the anchor
        scalar rmaxPre = 0;
        scalarList aFPre;
        const davof::cellSurface Spre = S;
        DynamicList<Pair<scalar>> anchorIterates(16);
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
            // only when it is printed (anchorDiagnostics) or written; the
            // anchor event print below needs it (anchorDiagnostics on)
            if (anchorDiagnostics_ || fitResidualPtr_)
            {
                scalarList aF;
                fractionsOf(S, aF);
                scalar rmax = 0;
                forAll(aF, k) rmax = max(rmax, mag(aF[k] - tF[k]));
                rmaxPre = rmax;
                aFPre = aF;
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
        // volume under the surface and feeds back into it.
        if (fitAnchor_ == "volume")
        {
            const scalar Vc = V[celli];
            const scalar target = alpha1In_[celli]*Vc;
            const scalar mMag = max(mag(m[celli]), 1e-3*h*h);
            // the cell volume under the surface by the volume rule at the
            // anchor's tolerance (anchorTolRel of V), the same rule as the
            // flux, so that the anchored volume is the flux's volume; with
            // it the slope dV/dc, the area of the wetted patch projected
            // along e3 (the Kromer-Bothe rule walks it as oint t1 dt2 on
            // the arcs it integrates anyway), negative where the rule has
            // none (a fallback, another rule); every evaluation counted
            scalar slope = -1;
            auto volumeOfT = [&](const davof::cellSurface& T) -> scalar
            {
                ++nAnchorEvals_;
                const scalar v = volumeRulePtr_->cellLiquidVolumeSlope(celli, T, anchorTolRel_, slope);
                anchorIterates.append(Pair<scalar>(T.c, v));
                return v;
            };
            // V(c) is non-decreasing in c (c larger: more liquid, dV/dc =
            // the projected wetted patch area, zero where the surface has
            // left the cell). anchorMethod ccs (the default): the
            // Consecutive Cubic Spline positioning of T. Maric, J. Comput.
            // Phys.: X 11 (2021) 100093, doi:10.1016/j.jcpx.2021.100093,
            // arXiv:2012.02227, on the paraboloid: the root of the cubic
            // Hermite interpolant through the bracket ends with the exact
            // volumes and slopes, re-bracketed, bisection-safeguarded (the
            // branch below). anchorMethod illinois (the method until
            // 2026-10-25): a bracket of the root (the step doubled from h/4
            // until the residual changes sign; the volume is confined to
            // [0, V], so a bracket exists within the cell's diagonal) and
            // the Illinois rule inside it (M. Dowell, P. Jarratt, BIT 11
            // (1971) 168-174, doi:10.1007/BF01934364), to the tolerance or
            // to a bracket of 1e-11 h; the secant alone failed in 3 to 50
            // cells per step where V(c) is flat or the start far
            // (2026-10-23).
            // anchorMethod newton: Newton steps on that slope, safeguarded
            // by the bracket the iterates build (an iterate outside it is
            // replaced by the bisection point, as the rtsafe of W. H. Press,
            // S. A. Teukolsky, W. T. Vetterling, B. P. Flannery, Numerical
            // Recipes, 3rd ed., Cambridge University Press, 2007, Sect.
            // 9.4, ISBN 978-0-521-88068-8), the step corrected by the
            // quadratic model of V(c) whose second derivative is the secant
            // of the two latest slopes (Halley's step with that estimate:
            // T. R. Scavo, J. B. Thoo, On the geometry of Halley's method,
            // Amer. Math. Monthly 102 (1995) 417-426, doi:10.2307/2975033);
            // where the slope is not usable (a fallback, a vanishing patch)
            // and no bracket exists yet, the bracket search and the Illinois
            // rule take over from the last iterate. On the sphere's ring row
            // (N = 32, Co 0.2) newton takes 3.7 evaluations per interface
            // cell and step against 6.5 for illinois (2026-10-25).
            // The volume residual alone is no criterion for the OFFSET: a
            // residual of tol V is an offset error of tol V / A_patch, large
            // where the wetted patch is small (a corner cut), and the plane
            // test showed exactly those faces off by 2e-4. The bracket of
            // the offset is therefore closed to 1e-11 h, or Newton's own
            // step r/A_patch is below it; the residual criterion applies
            // only where it is round-off.
            // The root is NOT unique where the target is 0 or V (a full or
            // empty cell flagged as an interface cell by |m_c| > wispTol
            // V^(2/3)): V(c) is flat beyond the offset at which the surface
            // leaves the cell's last corner, and any c there is a root.
            // illinois stops at the first bracket end its doubling reaches
            // (0.75 h or 1.75 h from the start on the sphere), newton near
            // the edge (0.72 to 0.80 h), up to 0.95 h apart; the cell's own
            // volume is the same, but its surface also cuts the prisms of
            // its faces and the second-ring faces, so the runs differ (the
            // sphere's ring row ended 0.5 % lower at N = 32 and 2.5 % lower
            // at N = 64 with newton, 2026-10-25). The offset of such a cell
            // is therefore DEFINED as the edge of the interval, the surface
            // tangent to the cell (kbOffsetRange, exact), in ccs and newton;
            // illinois keeps its bracket end as the ladder's reference path.
            const scalar tol = max(anchorTolRel_*Vc, 1e-14*Vc);
            const scalar cTol = 1e-11*h;
            const scalar cStart = S.c;
            const scalar rStart = volumeOfT(S) - target;
            const scalar AStart = slope;
            // the diagnostic of the quadrature (key curvedDebug): on an exact
            // state the surface before the anchor encloses alpha V, so r0 is
            // the quadrature's own error; the worst cell is printed with its
            // pyramids against the exact linear cut of the Detrixhe-Aslam tets
            if (curvedDebug_ && mag(rStart)/Vc > debugR0Max_)
            {
                debugR0Max_ = mag(rStart)/Vc;
                debugCellR0_ = celli;
                debugSurface_ = S;
                debugAlpha_ = alpha1In_[celli];
            }
            // the root by the given method on the surface T (its offset set),
            // from the start above; exitCode names the criterion that ended
            // the search (0 none: the better bracket end is kept)
            auto anchorSolve = [&]
            (
                const word& method, davof::cellSurface& T, label& exitCode
            ) -> bool
            {
                bool okAnchor = false;
                exitCode = 0;
                T.c = cStart;
                scalar c0 = cStart, r0 = rStart;
                const scalar A0 = AStart;
                auto volumeOf = [&]() -> scalar { return volumeOfT(T); };
                // A cell whose target is 0 or V (within the anchor's
                // residual tolerance 1e-14 V: a full or empty cell that
                // counts as an interface cell) has an interval of roots, on
                // which V(c) is flat. Its offset is DEFINED (Tomislav,
                // 2026-10-25) as the edge of that interval, the surface
                // tangent to the cell: the first offset at which the cell is
                // entirely liquid (cMax of kbOffsetRange) or entirely empty
                // (cMin), moved 1e-12 h onto the flat side so that no vertex
                // lies exactly on the surface. Exact, without an evaluation.
                // The method illinois keeps the bracket end it reaches there
                // (the runs before this definition).
                if
                (
                    method != "illinois"
                 && (target <= 1e-14*Vc || target >= (1 - 1e-14)*Vc)
                )
                {
                    DynamicList<davof::kbFace> cellFaces(24);
                    davof::kbCellFaces(mesh_, celli, cellFaces);
                    scalar cMin, cMax;
                    davof::kbOffsetRange(cellFaces, T, cMin, cMax);
                    const bool full = target >= (1 - 1e-14)*Vc;
                    T.c = full ? cMax + 1e-12*h : cMin - 1e-12*h;
                    exitCode = full ? 11 : 10;
                    return true;
                }
                if (mag(r0) <= 1e-14*Vc)
                {
                    okAnchor = true;
                    exitCode = 1;
                }
                else
                {
                    // the bracket [lo, hi] with r(lo) < 0 < r(hi), by either method
                    scalar lo = c0, hi = c0, rlo = r0, rhi = r0;
                    bool haveLo = false, haveHi = false;
                    if (method == "ccs")
                    {
                        // Consecutive Cubic Spline positioning (T. Maric,
                        // J. Comput. Phys.: X 11 (2021) 100093,
                        // doi:10.1016/j.jcpx.2021.100093, arXiv:2012.02227)
                        // on the paraboloid: once the root is bracketed, the
                        // iterate is the root of the cubic Hermite
                        // interpolant through the two bracket ends with the
                        // exact values V and the exact slopes dV/dc (the
                        // projected wetted patch area), re-bracketed at
                        // every step; an iterate outside the bracket, or a
                        // bracket not halved in two steps, is replaced by
                        // the bisection point. Between the breakpoints of
                        // V(c) (the vertex offsets, where d2V/dc2 jumps) V
                        // is smooth, a cubic for a plane, so the cubic model
                        // is exact there up to the quadrature; a bracket
                        // across a breakpoint shrinks onto one piece within
                        // a step or two. Until the bracket exists the step is
                        // Newton's from the latest point (at most h/2); a
                        // point without a usable slope falls to the bracket
                        // search and the Illinois rule below.
                        scalar c = c0, r = r0, A = A0;
                        scalar Alo = -1, Ahi = -1;
                        scalar wPrev = GREAT, wPrev2 = GREAT;
                        for (label it = 0; it < 60; ++it)
                        {
                            if (r < 0) { lo = c; rlo = r; Alo = A; haveLo = true; }
                            else { hi = c; rhi = r; Ahi = A; haveHi = true; }
                            const bool bracketed = haveLo && haveHi;
                            if (bracketed && hi - lo <= cTol)
                            {
                                T.c = (rlo*hi - rhi*lo)/(rlo - rhi);
                                if (!(T.c > lo && T.c < hi)) T.c = 0.5*(lo + hi);
                                okAnchor = true;
                                exitCode = 3;
                                break;
                            }
                            if (bracketed && mag(r) <= tol && hi - lo <= 1e-6*h)
                            {
                                okAnchor = true;
                                exitCode = 4;
                                break;
                            }
                            const bool slopeOK = std::isfinite(A) && A > 1e-12*h*h;
                            if (slopeOK && mag(r) <= tol && mag(r/A) <= cTol)
                            {
                                T.c = c - r/A;
                                okAnchor = true;
                                exitCode = 2;
                                break;
                            }
                            scalar cNew;
                            if (bracketed)
                            {
                                const scalar w = hi - lo;
                                const scalar sec = (rhi - rlo)/w;
                                const scalar dlo = (std::isfinite(Alo) && Alo >= 0) ? Alo : sec;
                                const scalar dhi = (std::isfinite(Ahi) && Ahi >= 0) ? Ahi : sec;
                                scalar sRoot = 0.5;
                                const bool ok = hermiteRoot01(rlo, dlo*w, rhi, dhi*w, sRoot);
                                cNew = ok ? lo + sRoot*w : 0.5*(lo + hi);
                                if (!(cNew > lo && cNew < hi) || w > 0.5*wPrev2)
                                {
                                    cNew = 0.5*(lo + hi);
                                }
                                wPrev2 = wPrev; wPrev = w;
                            }
                            else if (slopeOK)
                            {
                                scalar dc = -r/A;
                                if (mag(dc) > 0.5*h) dc = (dc > 0) ? 0.5*h : -0.5*h;
                                cNew = c + dc;
                            }
                            else
                            {
                                // no slope and no bracket: the bracket search
                                break;
                            }
                            T.c = cNew; c = cNew;
                            r = volumeOf() - target;
                            A = slope;
                            if (mag(r) <= 1e-14*Vc)
                            {
                                okAnchor = true;
                                exitCode = 5;
                                break;
                            }
                        }
                        if (!okAnchor) { c0 = c; r0 = r; }
                    }
                    else if (method == "newton")
                    {
                        scalar c = c0, r = r0, A = A0;
                        scalar cPrev = c0, APrev = A0;
                        bool havePrev = false;
                        for (label it = 0; it < 60; ++it)
                        {
                            if (r < 0) { lo = c; rlo = r; haveLo = true; }
                            else { hi = c; rhi = r; haveHi = true; }
                            const bool bracketed = haveLo && haveHi;
                            if (bracketed && hi - lo <= cTol)
                            {
                                T.c = (rlo*hi - rhi*lo)/(rlo - rhi);
                                if (!(T.c > lo && T.c < hi)) T.c = 0.5*(lo + hi);
                                okAnchor = true;
                                exitCode = 3;
                                break;
                            }
                            if (bracketed && mag(r) <= tol && hi - lo <= 1e-6*h)
                            {
                                okAnchor = true;
                                exitCode = 4;
                                break;
                            }
                            const bool slopeOK = std::isfinite(A) && A > 1e-12*h*h;
                            scalar cNew;
                            if (slopeOK)
                            {
                                // the offset criterion on the plain Newton step
                                // r/A, the first-order offset error of the
                                // iterate (the corrected step below can be much
                                // shorter where the slope changes fast)
                                const scalar dcN = -r/A;
                                if (mag(r) <= tol && mag(dcN) <= cTol)
                                {
                                    T.c = c + dcN;
                                    okAnchor = true;
                                    exitCode = 2;
                                    break;
                                }
                                scalar dc = dcN;
                                if (havePrev && c != cPrev)
                                {
                                    // the quadratic model V(c + dc) = V + A dc
                                    // + D dc^2/2 with D from the two latest
                                    // slopes; the correction keeps the sign of
                                    // the step and at most doubles it
                                    const scalar D = (A - APrev)/(c - cPrev);
                                    const scalar den = A + 0.5*D*dc;
                                    if (den > 0.5*A) dc = -r/den;
                                }
                                cNew = c + dc;
                                if (bracketed)
                                {
                                    if (!(cNew > lo && cNew < hi)) cNew = 0.5*(lo + hi);
                                }
                                else if (mag(dc) > 0.5*h)
                                {
                                    cNew = c + ((dc > 0) ? 0.5*h : -0.5*h);
                                }
                            }
                            else if (bracketed)
                            {
                                cNew = 0.5*(lo + hi);
                            }
                            else
                            {
                                // no slope and no bracket: the bracket search
                                break;
                            }
                            cPrev = c; APrev = A; havePrev = true;
                            T.c = cNew; c = cNew;
                            r = volumeOf() - target;
                            A = slope;
                            if (mag(r) <= 1e-14*Vc)
                            {
                                okAnchor = true;
                                exitCode = 5;
                                break;
                            }
                        }
                        // the bracket search and the Illinois rule continue
                        // from the last iterate
                        if (!okAnchor) { c0 = c; r0 = r; }
                    }
                    if (!okAnchor && !(haveLo && haveHi))
                    {
                        scalar step = min(mag(r0)/mMag, 0.25*h);
                        step = max(step, 1e-3*h);
                        if (r0 < 0)
                        {
                            lo = c0; rlo = r0; hi = c0; rhi = r0;
                            for (label k = 0; k < 12 && rhi < 0; ++k)
                            {
                                hi = lo + step; T.c = hi; rhi = volumeOf() - target;
                                if (rhi < 0) { lo = hi; rlo = rhi; }
                                step *= 2;
                            }
                        }
                        else
                        {
                            hi = c0; rhi = r0; lo = c0; rlo = r0;
                            for (label k = 0; k < 12 && rlo > 0; ++k)
                            {
                                lo = hi - step; T.c = lo; rlo = volumeOf() - target;
                                if (rlo > 0) { hi = lo; rhi = rlo; }
                                step *= 2;
                            }
                        }
                    }
                    if (!okAnchor && rlo <= 0 && rhi >= 0)
                    {
                        // Illinois (modified regula falsi), closed on the offset
                        label side = 0;
                        for (label it = 0; it < 80; ++it)
                        {
                            if (mag(rlo) <= 1e-14*Vc) { T.c = lo; okAnchor = true; exitCode = 6; break; }
                            if (mag(rhi) <= 1e-14*Vc) { T.c = hi; okAnchor = true; exitCode = 6; break; }
                            if (hi - lo <= cTol) { T.c = (rlo*hi - rhi*lo)/(rlo - rhi); if (!(T.c > lo && T.c < hi)) T.c = 0.5*(lo + hi); okAnchor = true; exitCode = 7; break; }
                            scalar c1 = (rlo*hi - rhi*lo)/(rlo - rhi);
                            // keep the iterate inside the bracket
                            if (!(c1 > lo && c1 < hi)) c1 = 0.5*(lo + hi);
                            T.c = c1;
                            const scalar r1 = volumeOf() - target;
                            if (mag(r1) <= 1e-14*Vc) { okAnchor = true; exitCode = 8; break; }
                            if (mag(r1) <= tol && hi - lo <= 1e-6*h) { okAnchor = true; exitCode = 9; break; }
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
                        T.c = (mag(rlo) < mag(rhi)) ? lo : hi;
                    }
                }
                return okAnchor;
            };

            label exitCode = 0;
            bool okAnchor = false;
            if (anchorMethod_ != "compare")
            {
                okAnchor = anchorSolve(anchorMethod_, S, exitCode);
            }
            else
            {
                // the diagnostic of the methods (anchorMethod compare): the
                // ladder's illinois and ccs from the same start, the run
                // continues with ccs; the cells whose offsets differ beyond
                // 1e-9 h are counted (all, and those at the edge of the flat
                // interval) and the first 20 per step printed with the
                // volumes and slopes at both offsets
                davof::cellSurface SI = S;
                label exitI = 0;
                const bool okI = anchorSolve("illinois", SI, exitI);
                okAnchor = anchorSolve("ccs", S, exitCode);
                if (mag(SI.c - S.c) > 1e-9*h)
                {
                    ++nCompareDiff_;
                    scalar AI = -1, AN = -1;
                    const scalar vI = volumeRulePtr_->cellLiquidVolumeSlope(celli, SI, anchorTolRel_, AI);
                    const scalar vN = volumeRulePtr_->cellLiquidVolumeSlope(celli, S, anchorTolRel_, AN);
                    // the classes: the edge of the flat interval; a
                    // difference the anchor cannot resolve (the two offsets
                    // enclose less than ten times the residual round-off
                    // 1e-14 V between them, |dc| A <= 1e-13 V: a nearly full
                    // or empty cell whose root is ill-conditioned, the
                    // slope A of the corner sliver small) or an illinois
                    // failure; the rest would be a real disagreement
                    const scalar dV = mag(SI.c - S.c)*max(max(AI, AN), scalar(0));
                    if (exitCode >= 10) ++nCompareDiffEdge_;
                    else if (dV <= 1e-13*Vc || !okI) ++nCompareDiffResolved_;
                    if (nCompare_ < 20)
                    {
                        ++nCompare_;
                        Pout<< "anchor compare: cell " << celli << " alpha " << alpha1In_[celli]
                            << " h " << h << " c0/h " << cStart/h << " r0/V " << rStart/Vc
                            << " A0/h^2 " << AStart/(h*h)
                            << " | illinois ok " << okI << " exit " << exitI
                            << " c/h " << SI.c/h << " r/V " << (vI - target)/Vc << " A/h^2 " << AI/(h*h)
                            << " | ccs ok " << okAnchor << " exit " << exitCode
                            << " c/h " << S.c/h << " r/V " << (vN - target)/Vc << " A/h^2 " << AN/(h*h)
                            << " | (cCcs - cI)/h " << (S.c - SI.c)/h
                            << " |dc| A/V " << dV/Vc << nl;
                    }
                }
            }
            if (exitCode >= 10) ++nAnchorEdge_;
            if (!okAnchor) ++nAnchorFail_;
            // the face fractions the anchored surface cuts against the
            // targets (the face shift): a diagnostic, computed when printed
            // (anchorDiagnostics), with the event print below
            if (anchorDiagnostics_)
            {
                scalarList aF;
                fractionsOf(S, aF);
                scalar dmax = 0;
                forAll(aF, k) dmax = max(dmax, mag(aF[k] - tF[k]));
                anchorShiftMax_ = max(anchorShiftMax_, dmax);
                anchorShiftSum_ += dmax;

            // The event print of the plane test (2026-10-25): a plane-like
            // surface (|b| h below 1e-6) whose face residual before the
            // anchor exceeds 1e-6 is no round-off; the translated plane at
            // N = 160 showed one such cell per run (residual 7.7e-5 and
            // 2.1e-4, the anchored faces off by 4.2e-6 and 1.6e-4, then the
            // transported fractions), with the targets at 1e-10. Printed: the
            // surface, the anchor's iterates, per face the targets, the
            // rule's fractions before and after the anchor and the quadratic
            // face cut, per fan triangle the rule's wetted area against the
            // quadratic cut with the vertex values, and the cell volume of
            // the rule against the tight split quadrature. At most 4 per step.
            const scalar bh = (mag(S.b11) + mag(S.b12) + mag(S.b22))*h;
            static label nEventStep = -1, nEventPrinted = 0;
            if (nEventStep != mesh_.time().timeIndex())
            {
                nEventStep = mesh_.time().timeIndex();
                nEventPrinted = 0;
            }
            if (rmaxPre > 1e-6 && bh < 1e-6 && nEventPrinted < 4)
            {
                ++nEventPrinted;
                auto phiS = [&](const point& x) { return S.value(x); };
                const pointField& P = mesh_.points();
                Pout<< "anchor event: cell " << celli << " h " << h
                    << " alpha " << alpha1In_[celli] << " ok " << okAnchor
                    << " residual before " << rmaxPre << " shift after " << dmax << nl
                    << "  surface o " << S.o << " e3 " << S.e3 << " b " << S.b11
                    << " " << S.b12 << " " << S.b22 << " c before " << Spre.c
                    << " after " << S.c << " (c/h " << (S.c - Spre.c)/h << ")" << nl;
                forAll(anchorIterates, it)
                {
                    Pout<< "  iterate " << it << " c " << anchorIterates[it].first()
                        << " V " << anchorIterates[it].second()
                        << " (V - target)/V " << (anchorIterates[it].second() - target)/Vc << nl;
                }
                {
                    label nl2 = 0;
                    const scalar vKb = volumeRulePtr_->cellLiquidVolume(celli, S, anchorTolRel_);
                    const scalar vTight = davof::cellLiquidVolumeQuadric
                    (
                        mesh_, celli, phiS, 1e-15, 20, nl2
                    );
                    Pout<< "  cell volume: rule " << vKb << " tight split " << vTight
                        << " (rule - tight)/V " << (vKb - vTight)/Vc
                        << " target " << target << nl;
                }
                {
                    const labelList& cp = mesh_.cellPoints()[celli];
                    Pout<< "  vertex phi/h:";
                    forAll(cp, ip) Pout<< " " << S.value(P[cp[ip]])/h;
                    Pout<< nl;
                }
                forAll(cf, k)
                {
                    label nfb = 0;
                    const scalar aq = volumeRulePtr_->davofVolumeRule::faceLiquidFraction(cf[k], vector::zero, S, nfb);
                    Pout<< "  face " << cf[k] << " target " << tF[k]
                        << " rule before " << (k < aFPre.size() ? aFPre[k] : -1)
                        << " after " << aF[k] << " quadratic cut " << aq
                        << " (rule - cut) " << aF[k] - aq << nl;
                    // per fan triangle: the rule's area against the quadratic cut
                    const face& f = mesh_.faces()[cf[k]];
                    const point& xf = mesh_.faceCentres()[cf[k]];
                    forAll(f, ip)
                    {
                        const label ip1 = f.nextLabel(ip);
                        davof::kbFace F;
                        F.p.setSize(3);
                        F.p[0] = xf; F.p[1] = P[f[ip]]; F.p[2] = P[ip1];
                        F.setFrame();
                        davof::kbSums sums;
                        const bool okT = davof::kbIntegrateFace(F, S, 1e-16*F.area(), GREAT, false, sums);
                        label nfq = 0;
                        const scalar aT = davof::triLiquidArea(xf, P[f[ip]], P[ip1], phiS, 6, nfq);
                        Pout<< "    tri " << ip << " area " << F.area() << " wet rule " << sums.lastArea
                            << (okT ? "" : " NOT RESOLVED") << " cut " << aT
                            << " (rule - cut)/area " << (sums.lastArea - aT)/F.area()
                            << " arcs " << sums.nArcs << " phi/h " << S.value(xf)/h
                            << " " << S.value(P[f[ip]])/h << " " << S.value(P[ip1])/h << nl;
                    }
                }
            }
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


void Foam::advection::isoAdvection::anchorTensorCentroids()
{
    addProfilingInFunction(geometricVoF);
    const DynamicField<label>& interfaceLabels = surf_->interfaceLabels();
    const vectorField& usedN = davofPtr_->usedNormals();
    const vectorField& normal = surf_->normal().primitiveField();
    const vectorField& centre = surf_->centre().primitiveField();
    const vectorField& m = davofPtr_->areaNormal().primitiveField();
    const vectorField& C = mesh_.cellCentres();
    const scalarField& V = mesh_.V();
    const pointField& pts = mesh_.points();
    const labelListList& cellPts = mesh_.cellPoints();
    vectorField& xT = davofPtr_->tensorCentres();

    forAll(interfaceLabels, i)
    {
        const label celli = interfaceLabels[i];
        const point o = (mag(normal[celli]) > 0) ? centre[celli] : C[celli];
        xT[celli] = o;
        vector nIn = usedN[i];
        if (!(mag(nIn) > 0))
        {
            const scalar mMag = mag(m[celli]);
            if (!(mMag > 0)) continue;
            nIn = m[celli]/mMag;
        }
        // The plane of the frame through o (c = 0), to be shifted along e3
        davof::cellSurface T;
        T.paraboloid = false;
        T.e3 = -nIn;
        davof::frameOf(T.e3, T.e1, T.e2);
        T.o = o;
        T.c = 0;
        const scalar Vc = V[celli];
        const scalar h = Foam::cbrt(Vc);
        const scalar target = alpha1In_[celli]*Vc;

        // The bracket: the cell's extent along e3 about o, where the
        // volume under the plane is 0 and V (V(c) is non-decreasing)
        scalar lo = 0, hi = 0;
        for (const label pi : cellPts[celli])
        {
            const scalar s = (pts[pi] - o) & T.e3;
            lo = min(lo, s);
            hi = max(hi, s);
        }
        lo -= 1e-9*h;
        hi += 1e-9*h;
        scalar rlo = -target, rhi = Vc - target;
        auto residual = [&](const scalar cc) -> scalar
        {
            T.c = cc;
            return volumeRulePtr_->cellLiquidVolume(celli, T, anchorTolRel_) - target;
        };
        bool ok = false;
        if (rlo >= 0)
        {
            T.c = lo; ok = true;
        }
        else if (rhi <= 0)
        {
            T.c = hi; ok = true;
        }
        else
        {
            // Illinois (modified regula falsi) on the offset, closed to a
            // round-off residual or to 1e-13 h of the offset
            label side = 0;
            for (label it = 0; it < 200; ++it)
            {
                scalar c1 = (rlo*hi - rhi*lo)/(rlo - rhi);
                if (!(c1 > lo && c1 < hi)) c1 = 0.5*(lo + hi);
                const scalar r1 = residual(c1);
                if (mag(r1) <= 1e-15*Vc || hi - lo <= 1e-13*h)
                {
                    T.c = c1; ok = true; break;
                }
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
            if (!ok) T.c = 0.5*(lo + hi);
        }
        if (!ok) ++nCentroidFail_;
        xT[celli] = o + T.c*T.e3;
    }
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

    // the tilt diagnostic's face records of the step
    if (tiltDiagnostics_) tiltDiagBegin();

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
        const point& xf, bool& secondRing, label& surfCell
    ) -> const davof::cellSurface*
    {
        // surfCell: the local cell whose surface is returned (-1 for a
        // remote surface or none), for the tilt diagnostic's record of
        // the choice; the choice itself is unchanged
        secondRing = false;
        surfCell = -1;
        if (up >= 0 && hasSurface_.test(up)) { surfCell = up; return &cellSurfaces_[up]; }
        if (other >= 0 && hasSurface_.test(other)) { surfCell = other; return &cellSurfaces_[other]; }
        if (remoteHas) return remote;
        const davof::cellSurface* best = nullptr;
        label bestCell = -1;
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
        auto considerD = [&](const davof::cellSurface& T, const scalar d, const label c)
        {
            if (!best || d < dMin*(1 - 1e-9))
            {
                dMin = d; best = &T; bestCell = c;
            }
            else if (d <= dMin*(1 + 1e-9) && lexLess(T.o, best->o))
            {
                dMin = min(dMin, d); best = &T; bestCell = c;
            }
        };
        auto consider = [&](const davof::cellSurface& T, const point& o, const label c)
        {
            considerD(T, magSqr(xf - o), c);
        };
        if (remoteRingHas) considerD(*remoteRing, remoteRingDist, -1);
        auto scan = [&](const label c)
        {
            if (c < 0) return;
            for (const label nb : cellCells[c])
            {
                if (!hasSurface_.test(nb)) continue;
                consider(cellSurfaces_[nb], cellSurfaces_[nb].sphere ? mesh_.cellCentres()[nb] : cellSurfaces_[nb].o, nb);
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
                    consider(T, T.o, -1);
                }
            }
        };
        scan(up);
        scan(other);
        secondRing = (best != nullptr);
        surfCell = bestCell;
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
        label surfCell = -1;
        const davof::cellSurface* Sp = surfaceOf
        (
            up, other, remote, remoteHas, remoteRing, remoteRingHas, remoteRingDist,
            Cf[facei], secondRing, surfCell
        );
        if (!Sp)
        {
            ++nCurvedSkipped_;
            if (tiltDiagnostics_)
            {
                tiltDiagRecordFace(facei, -1, UInterp.interpolate(Cf[facei], local)*dt);
            }
            return;
        }
        const davof::cellSurface& S = *Sp;

        // The displacement of the step at the face
        const vector d = UInterp.interpolate(Cf[facei], local)*dt;

        // the tilt diagnostic: the cell whose surface sets this face
        if (tiltDiagnostics_) tiltDiagRecordFace(facei, surfCell, d);

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

        // The face fraction: the pre-image's liquid fraction under the
        // surface; with faceState vector its wetted area vector as well
        // (the pre-image's fan triangles are translates of the face's)
        label nfb = 0;
        scalar af;
        vector aWet(Zero);
        if (wetVecPtr_)
        {
            af = volumeRulePtr_->faceWettedArea(facei, -d, S, aWet, nfb);
        }
        else
        {
            af = volumeRulePtr_->faceLiquidFraction(facei, -d, S, nfb);
        }
        nCurvedFallback_ += nfb;
        if (af < -1e-12 || af > 1 + 1e-12)
        {
            ++nAlphafClipped_;
        }
        af = min(max(af, scalar(0)), scalar(1));
        setFaceValue(*alphafPtr_, facei, af);
        if (wetVecPtr_)
        {
            setFaceValue(*wetVecPtr_, facei, aWet);
        }
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
        syncProcPatches(dVf_, phi_, false, alphafPtr_, wetVecPtr_);
    }
    nCurvedLeaves_ = volumeRulePtr_->nLeaves();
    nRuleFallback_ += volumeRulePtr_->nFallback();
    volumeRulePtr_->resetCounters();

    // the tilt diagnostic: the traces of the exact plane at t^n on the
    // pre-images of the recorded faces, and which cells' surfaces set a
    // traced face in this step
    if (tiltDiagnostics_) tiltDiagFaces(dt);

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
    label nAnchE = nAnchorEvals_, nAnchEdge = nAnchorEdge_;
    label nCmpD = nCompareDiff_, nCmpDE = nCompareDiffEdge_, nCmpDR = nCompareDiffResolved_;
    label nDef = nDefectCells_, nRuleFb = nRuleFallback_ + volumeRulePtr_->nFallback();
    label nCenF = nCentroidFail_;
    scalar fitRMax = fitResMax_, fitRSum = fitResSum_, ashMax = anchorShiftMax_;
    scalar ashSum = anchorShiftSum_, defMax = defectMax_, defSum = defectSum_;
    if (Pstream::parRun())
    {
        reduce(nCenF, sumOp<label>());
        reduce(nFaces, sumOp<label>()); reduce(nSecond, sumOp<label>());
        reduce(nBound, sumOp<label>()); reduce(nSkip, sumOp<label>());
        reduce(nClip, sumOp<label>()); reduce(nLin, sumOp<label>());
        reduce(nLeaves, sumOp<label>()); reduce(nFitC, sumOp<label>());
        reduce(nFitConv, sumOp<label>()); reduce(nFitSt, sumOp<label>());
        reduce(nFitSing, sumOp<label>()); reduce(nFitNoBr, sumOp<label>());
        reduce(nFitNoFr, sumOp<label>()); reduce(nAnchF, sumOp<label>());
        reduce(nAnchL, sumOp<label>()); reduce(nDef, sumOp<label>());
        reduce(nAnchE, sumOp<label>()); reduce(nAnchEdge, sumOp<label>());
        reduce(nCmpD, sumOp<label>()); reduce(nCmpDE, sumOp<label>());
        reduce(nCmpDR, sumOp<label>());
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
            Info<< ", leaves " << nAnchL << ", evaluations " << nAnchE
                << " (per cell " << scalar(nAnchE)/max(nFitC, label(1)) << ")"
                << ", edge " << nAnchEdge;
        }
        if (tensorCentroid_ == "anchored")
        {
            Info<< "; anchored centroids: failed " << nCenF;
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
    if (anchorMethod_ == "compare")
    {
        Info<< "isoAdvection: anchor compare: illinois and ccs differ in "
            << nCmpD << " cells, " << nCmpDE << " of them at the edge of the flat interval, "
            << nCmpDR << " within the anchor's resolution or an illinois failure, "
            << nCmpD - nCmpDE - nCmpDR << " other"
            << endl;
    }
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
