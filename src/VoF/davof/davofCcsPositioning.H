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

Namespace
    Foam::davof

Description
    Consecutive Cubic Spline (CCS) positioning of a PLIC plane in a cell
    (2026-11-04). The kernel hermiteRoot01 is the one the swept pieces'
    volume anchor has used since 2026-10-25 (isoAdvectionCurved.C, moved
    here unchanged so that the anchor and the reconstruction share it).

    Reference:
    \verbatim
        T. Maric,
        Iterative Volume-of-Fluid interface positioning in general
        polyhedrons with Consecutive Cubic Spline interpolation,
        Journal of Computational Physics: X 11 (2021) 100093,
        doi:10.1016/j.jcpx.2021.100093, arXiv:2012.02227
    \endverbatim

    ccsPositionPlane positions the plane {x : (x - x_c).n = p} of a cell so
    that the liquid volume on the side (x - x_c).n > p equals alpha V_c. The
    volume V(p) and the slope dV/dp = -A(p) (A the area of the cut polygon)
    come from cutCellPLIC::calcSubCell, the same cut the vofCutCell
    positioning and the reconstruction's diagnostics use; V(p) is
    non-increasing in p and, between two adjacent vertex levels, exactly a
    cubic (its derivative, the section area, is a quadratic of p). The
    iteration starts from a given p0 (the reconstruction passes the
    explicit Gauss-theorem position pPlane_c, exact on a planar state), and
    keeps the bracket [lo, hi] with r(lo) = V(lo) - alpha V_c >= 0 >= r(hi),
    whose first ends are the vertex extrema (V = V_c and 0, A = 0, known
    without a cut). The iterate is the root of the cubic Hermite
    interpolant of r through the two CONSECUTIVE evaluated points with
    their slopes (the CCS step; exact up to round-off once both lie between
    the same vertex levels as the root, so that the third cut is the
    answer): through the bracket ends once both are evaluated, else through
    the previous and the latest point, extrapolated, taken where the root
    lies inside the bracket; Newton's step r/A from the single first point;
    the bisection point where an iterate leaves the bracket or the bracket
    was not halved over two steps. It stops when |r| <= tolRel V_c (exit
    1), when the residual stops decreasing at a cubic step near the
    tolerance (the round-off floor; exit 2, the best iterate is returned),
    when the bracket has closed to round-off (exit 3) or at maxIter (exit
    4). A target at an end of [0, V_c] (alpha <= 0 or >= 1 within
    round-off) takes the vertex extreme: the full or empty cell, no plane.
    The cut object holds the returned plane (re-cut at the best iterate
    when the last cut was not the best, counted in nEvals).

SourceFiles
    none (inline)

\*---------------------------------------------------------------------------*/

#ifndef davofCcsPositioning_H
#define davofCcsPositioning_H

#include "fvMesh.H"
#include "cutCellPLIC.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace davof
{

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


//- The root of the cubic Hermite interpolant p(s) with p(0) = f0, p(1) = f1
//  and the derivatives p'(0) = d0, p'(1) = d1 (in s), by Newton on the
//  polynomial from s0, anywhere on the line (the CCS step proper: the cubic
//  through the two CONSECUTIVE iterates, extrapolated beyond them where
//  the root lies outside; exact where the cubic is the function). False
//  when Newton on the polynomial fails to settle.
inline bool hermiteRootFrom
(
    const Foam::scalar f0, const Foam::scalar d0,
    const Foam::scalar f1, const Foam::scalar d1,
    const Foam::scalar s0,
    Foam::scalar& s
)
{
    using Foam::scalar; using Foam::label; using Foam::mag;
    const scalar scale = mag(f0) + mag(f1);
    s = s0;
    for (label it = 0; it < 40; ++it)
    {
        const scalar t = s, t2 = t*t, t3 = t2*t;
        const scalar v = (2*t3 - 3*t2 + 1)*f0 + (t3 - 2*t2 + t)*d0 + (-2*t3 + 3*t2)*f1 + (t3 - t2)*d1;
        const scalar dv = (6*t2 - 6*t)*f0 + (3*t2 - 4*t + 1)*d0 + (-6*t2 + 6*t)*f1 + (3*t2 - 2*t)*d1;
        if (mag(v) <= 1e-15*scale) return true;
        if (!(dv != 0) || !std::isfinite(dv)) return false;
        const scalar ds = v/dv;
        s -= ds;
        if (!std::isfinite(s) || mag(s) > 1e6) return false;
        if (mag(ds) <= 1e-14*(1 + mag(s))) return true;
    }
    return true;
}


//- The record of one CCS positioning
struct ccsRecord
{
    label nEvals = 0;        // cuts made (the final re-cut included)
    label exitCode = 0;      // 1 tol, 2 floor (best returned), 3 bracket closed, 4 maxIter
    scalar pStart = 0;       // the start (clamped into the vertex range)
    scalar rStart = 0;       // (V(pStart) - alpha V)/V
    scalar pFinal = 0;
    scalar rFinal = 0;       // (V(pFinal) - alpha V)/V of the returned cut
    label nHermite = 0;      // CCS steps taken
    label nNewton = 0;       // Newton steps taken
    label nBisect = 0;       // bisection steps taken
};


//- Position the plane of celli along the unit normal n by CCS from p0 so
//  that the cut volume equals alpha V_c (see the Description). Returns the
//  status of the final cut (0 cut, -1 fully liquid, +1 fully gas); the cut
//  object holds it.
inline label ccsPositionPlane
(
    const fvMesh& mesh,
    cutCellPLIC& cut,
    const label celli,
    const scalar alpha,
    const vector& n,
    const scalar p0,
    const scalar tolRel,
    const label maxIter,
    ccsRecord& rec
)
{
    const scalar Vc = mesh.V()[celli];
    // the target within [0, V_c]: a volume fraction outside it (a round-off
    // wisp of -3e-15, 2026-11-04) has no plane but the vertex extreme
    const scalar target = min(max(alpha, scalar(0)), scalar(1))*Vc;
    const scalar tol = tolRel*Vc;
    const point& xc = mesh.C()[celli];
    const pointField& P = mesh.points();

    // the vertex extrema: V(pMin) = V_c, V(pMax) = 0, A = 0 at both
    scalar pMin = GREAT, pMax = -GREAT;
    for (const label pi : mesh.cellPoints(celli))
    {
        const scalar f = (P[pi] - xc) & n;
        pMin = min(pMin, f);
        pMax = max(pMax, f);
    }
    const scalar h = Foam::cbrt(Vc);

    rec = ccsRecord();

    label status = 0;
    scalar r = 0, A = 0;
    auto evaluate = [&](const scalar p) -> void
    {
        status = cut.calcSubCell(celli, p, n);
        const scalar V = (status == 0) ? cut.subCellVolume() : ((status < 0) ? Vc : 0);
        A = (status == 0) ? mag(cut.faceArea()) : 0;
        r = V - target;
        ++rec.nEvals;
    };

    // the bracket: lo with r >= 0 (more liquid), hi with r <= 0
    scalar lo = pMin, hi = pMax, rlo = Vc - target, rhi = -target, Alo = 0, Ahi = 0;
    bool loEval = false, hiEval = false;

    // a target at an end of [0, V_c]: the plane at the vertex extreme, exact
    // (the cell is then full or empty: no plane, as vofCutCell returns)
    scalar p = (target <= 0) ? pMax : ((target >= Vc) ? pMin : min(max(p0, pMin), pMax));
    evaluate(p);
    rec.pStart = p;
    rec.rStart = r/Vc;
    scalar pBest = p, rBest = r;
    label statusBest = status;
    scalar wPrev = GREAT, wPrev2 = GREAT;
    // the previous evaluated point (the consecutive cubic)
    scalar pPrev = 0, rPrev = 0, APrev = 0;
    bool havePrev = false;

    for (label it = 0; it < maxIter; ++it)
    {
        if (mag(r) <= tol || target <= 0 || target >= Vc)
        {
            rec.exitCode = 1;
            break;
        }
        if (r > 0) { lo = p; rlo = r; Alo = A; loEval = true; }
        else       { hi = p; rhi = r; Ahi = A; hiEval = true; }
        const scalar w = hi - lo;
        if (w <= 1e-15*max(h, max(mag(lo), mag(hi))))
        {
            rec.exitCode = 3;
            break;
        }

        scalar pNew;
        bool bisect = false, cubic = false;
        if (loEval && hiEval)
        {
            // the CCS step on the bracket: the Hermite cubic of -r
            // (increasing in p) through the bracket ends, in s = (p - lo)/w
            scalar s = 0.5;
            const bool ok = hermiteRoot01(-rlo, Alo*w, -rhi, Ahi*w, s);
            pNew = ok ? lo + s*w : 0.5*(lo + hi);
            if (!(pNew > lo && pNew < hi) || w > 0.5*wPrev2)
            {
                pNew = 0.5*(lo + hi);
                bisect = true;
            }
            else
            {
                ++rec.nHermite;
                cubic = true;
            }
        }
        else
        {
            // the CCS step proper: the cubic through the two consecutive
            // iterates (the previous and this one, with their slopes),
            // its root found from Newton's estimate and taken where it
            // lies inside the bracket; Newton's step from this point where
            // there is no previous one or the cubic gives nothing
            bool done = false;
            if (havePrev && A > 0 && p != pPrev)
            {
                const scalar wc = p - pPrev;
                scalar s = 1 + r/(A*wc);
                if (hermiteRootFrom(rPrev, -APrev*wc, r, -A*wc, s, s))
                {
                    pNew = pPrev + s*wc;
                    if (pNew > lo && pNew < hi)
                    {
                        ++rec.nHermite;
                        cubic = true;
                        done = true;
                    }
                }
            }
            if (!done)
            {
                if (A > 0)
                {
                    // Newton from the latest point: r(p + dp) = r - A dp = 0
                    pNew = p + r/A;
                    if (!(pNew > lo && pNew < hi))
                    {
                        pNew = 0.5*(lo + hi);
                        bisect = true;
                    }
                    else
                    {
                        ++rec.nNewton;
                    }
                }
                else
                {
                    pNew = 0.5*(lo + hi);
                    bisect = true;
                }
            }
        }
        if (bisect) ++rec.nBisect;
        wPrev2 = wPrev; wPrev = w;

        pPrev = p; rPrev = r; APrev = A; havePrev = true;
        p = pNew;
        evaluate(p);
        if (mag(r) < mag(rBest))
        {
            pBest = p; rBest = r; statusBest = status;
        }
        else if (cubic && (mag(rBest) <= 100*tol || (loEval && hiEval)))
        {
            // the residual did not decrease at a cubic step near the
            // tolerance (or with both bracket ends evaluated): the
            // round-off floor; the best cut is returned
            rec.exitCode = 2;
            break;
        }
        if (it == maxIter - 1)
        {
            rec.exitCode = 4;
        }
    }
    if (rec.exitCode == 0)
    {
        rec.exitCode = (mag(r) <= tol) ? 1 : 4;
    }

    if (p != pBest)
    {
        evaluate(pBest);
    }
    rec.pFinal = pBest;
    rec.rFinal = rBest/Vc;
    return statusBest;
}

} // End namespace davof
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
