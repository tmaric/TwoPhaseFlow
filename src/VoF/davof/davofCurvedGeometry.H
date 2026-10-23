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
    davof :: the curved clipping surface of a cell and the geometry of the
    dual advection with it (2026-10-23).

    The surface of a cell is a quadric phi(x), negative in the liquid:
    either the paraboloid of leia's regeneration (davofRegeneration.H,
    copied here with its conventions) in the frame (o; e1, e2, e3), e3 the
    normal pointing OUT of the liquid,

        phi = w - c - a1 u - a2 v - 1/2 (b11 u^2 + 2 b12 u v + b22 v^2),

    u = (x - o).e1, v = (x - o).e2, w = (x - o).e3, or the exact sphere
    phi = |x - xc|^2 - R^2 (the control). Both are quadratic polynomials in
    x, so along any straight segment phi is a quadratic in the parameter
    and its roots are exact (segmentLiquidFraction), and on any planar
    triangle the quadratic Lagrange interpolant of phi reproduces phi, so
    the quadratic face cut of davofQuadraticFaceGeometry.H is the exact
    liquid area up to its Gauss quadrature of the sliver (triLiquidArea).

    The dual advection of a face f over one step moves its pre-image
    f - d (d the displacement of the step, d = U dt for a uniform velocity)
    onto f: the new face fraction is the liquid fraction of the pre-image
    under the surface at t^n (faceCutDisplaced), and the liquid volume flux
    through f is the liquid volume of the prism swept between f - d and f
    (prismLiquidVolume: the liquid fraction of the segment [p - d, p]
    integrated over the face with an adaptive Gauss rule, each segment's
    fraction exact, so the only error is the quadrature's, driven below
    tolRel of the prism volume). Both the volume flux and the face update
    are therefore the moments of ONE surface per cell, which is the
    requirement measured in the sagitta studies of 2026-10-16.

    The six-coefficient fit (fitParaboloidFull) is leia's Gauss-Newton of
    2026-10-12 ported: the paraboloid that cuts every face of the cell at
    the cell's own face fractions, the ring curvature tensor as the prior
    through the Marquardt damping (the system is under-determined when
    fewer than six faces are cut; the step has no component in the null
    space, so those components keep the start, the ring tensor).

SourceFiles
    (header only)

\*---------------------------------------------------------------------------*/

#ifndef davofCurvedGeometry_H
#define davofCurvedGeometry_H

#include "fvMesh.H"
#include "FixedList.H"
#include "davofSimplexGeometry.H"
#include "davofQuadraticFaceGeometry.H"
#include <functional>

namespace Foam
{
namespace davof
{

// * * * * * * * * * * * * * * * The surface * * * * * * * * * * * * * * * //

//- The clipping surface of one cell: a paraboloid in the frame (o; e1, e2,
//  e3) with e3 pointing OUT of the liquid (value < 0 in the liquid), a
//  plane when paraboloid is false, or the sphere |x - o| = R when sphere is
//  true (o the centre, value = |x - o|^2 - R^2: negative inside the ball,
//  the liquid).
struct cellSurface
{
    bool paraboloid = false;
    bool sphere = false;
    point o = Zero;
    vector e1 = Zero, e2 = Zero, e3 = Zero;
    scalar b11 = 0, b12 = 0, b22 = 0, c = 0;
    scalar a1 = 0, a2 = 0;
    scalar R = 0;

    inline scalar value(const point& x) const
    {
        if (sphere)
        {
            return magSqr(x - o) - R*R;
        }
        const vector r = x - o;
        const scalar w = r & e3;
        if (!paraboloid) return w - c;
        const scalar u = r & e1, v = r & e2;
        return w - c - a1*u - a2*v - 0.5*(b11*u*u + 2*b12*u*v + b22*v*v);
    }
};


//- An orthonormal frame (e1, e2) of the plane normal to the unit vector e3.
inline void frameOf(const vector& e3, vector& e1, vector& e2)
{
    vector a = (mag(e3.x()) < 0.9) ? vector(1, 0, 0) : vector(0, 1, 0);
    e1 = a - (a & e3)*e3;
    e1 /= max(mag(e1), VSMALL);
    e2 = e3 ^ e1;
}


// * * * * * * * * * * * * * Segments and triangles  * * * * * * * * * * * //

//- The measure of {s in [0, 1] : q(s) <= 0} for the quadratic
//  q(s) = q0 + B s + C s^2 (exact).
inline scalar quadraticNegativeMeasure(const scalar q0, const scalar B, const scalar C)
{
    auto negLen = [&](const scalar s0, const scalar s1) -> scalar
    {
        // q <= 0 on [s0, s1] when q is of one sign there: test the midpoint
        const scalar sm = 0.5*(s0 + s1);
        return (q0 + B*sm + C*sm*sm <= 0) ? (s1 - s0) : 0;
    };
    const scalar scale = mag(q0) + mag(B) + mag(C);
    if (scale <= VSMALL) return 1;              // phi == 0 along the segment: liquid
    if (mag(C) <= 1e-14*scale)
    {
        // linear
        if (mag(B) <= 1e-14*scale) return (q0 <= 0) ? 1 : 0;
        const scalar s = -q0/B;
        if (s <= 0) return (q0 + B <= 0) ? 1 : 0;
        if (s >= 1) return (q0 <= 0) ? 1 : 0;
        return (B > 0) ? s : (1 - s);
    }
    const scalar disc = B*B - 4*C*q0;
    if (disc <= 0)
    {
        return (C > 0) ? 0 : 1;                 // no real root: the sign of C
    }
    const scalar sq = Foam::sqrt(disc);
    // stable roots
    const scalar qq = -0.5*(B + sign(B)*sq);
    scalar r1 = (qq != 0) ? q0/qq : -B/(2*C);
    scalar r2 = (C != 0) ? qq/C : r1;
    if (r1 > r2) Swap(r1, r2);
    // the segment [0, 1] split by the roots inside it
    scalar s0 = 0, len = 0;
    if (r1 > 0 && r1 < 1) { len += negLen(s0, r1); s0 = r1; }
    if (r2 > 0 && r2 < 1) { len += negLen(s0, r2); s0 = r2; }
    len += negLen(s0, 1);
    return len;
}


//- The liquid fraction of the straight segment from p0 to p1 under the
//  quadric phi: exact from three evaluations (phi is quadratic along it).
template<class PhiAt>
inline scalar segmentLiquidFraction(const point& p0, const point& p1, PhiAt&& phi)
{
    const scalar f0 = phi(p0);
    const scalar f1 = phi(p1);
    const scalar fm = phi(0.5*(p0 + p1));
    const scalar C = 2*(f1 + f0 - 2*fm);
    const scalar B = f1 - f0 - C;
    return quadraticNegativeMeasure(f0, B, C);
}


//- Liquid area of the triangle (a, b, c) under phi: the quadratic cut (six
//  values: vertices and edge midpoints; exact for a quadric); where it
//  falls back (more than two boundary crossings, an arc that is no graph
//  over its chord, a vertex exactly on the surface) the triangle is split
//  into four and the cut retried, down to `levels` more splits, then the
//  linear cut (counted). From leia's davofRegeneration.H.
template<class PhiAt>
inline scalar triLiquidArea
(
    const point& a, const point& b, const point& c,
    PhiAt&& phi,
    const label levels,
    label& nFallback
)
{
    const scalar ph[6] =
    {
        phi(a), phi(b), phi(c),
        phi(0.5*(a + b)), phi(0.5*(b + c)), phi(0.5*(c + a))
    };
    const scalar A = 0.5*mag((b - a) ^ (c - a));
    if (A <= VSMALL) return 0;
    bool fb = false;
    const scalar fr = triQuadraticNegativeFraction(a, b, c, ph, fb);
    if (!fb) return fr*A;
    if (levels <= 0)
    {
        ++nFallback;
        return triNegativeFraction(ph[0], ph[1], ph[2])*A;
    }
    const point mab = 0.5*(a + b), mbc = 0.5*(b + c), mca = 0.5*(c + a);
    return
        triLiquidArea(a, mab, mca, phi, levels - 1, nFallback)
      + triLiquidArea(mab, b, mbc, phi, levels - 1, nFallback)
      + triLiquidArea(mca, mbc, c, phi, levels - 1, nFallback)
      + triLiquidArea(mab, mbc, mca, phi, levels - 1, nFallback);
}


//- The liquid area fraction of the face faceI displaced by d (every vertex
//  and the centre moved by d: the pre-image of the face under a uniform
//  displacement), cut by phi on the fan about the displaced centre.
template<class PhiAt>
inline scalar faceCutDisplaced
(
    const fvMesh& mesh,
    const label faceI,
    const vector& d,
    PhiAt&& phi,
    label& nFallback
)
{
    const pointField& P = mesh.points();
    const face& f = mesh.faces()[faceI];
    const point xf = mesh.faceCentres()[faceI] + d;
    scalar liqA = 0, totA = 0;
    forAll(f, ip)
    {
        const label ip1 = f.nextLabel(ip);
        const point p0 = P[f[ip]] + d;
        const point p1 = P[ip1] + d;
        const scalar triA = 0.5*mag((p0 - xf) ^ (p1 - xf));
        if (triA <= VSMALL) continue;
        liqA += triLiquidArea(xf, p0, p1, phi, 6, nFallback);
        totA += triA;
    }
    return (totA > VSMALL) ? liqA/totA : (phi(xf) < 0 ? scalar(1) : scalar(0));
}


// * * * * * * * The triangle integral split at the quadric traces * * * * * //

//- Gauss-Legendre nodes and weights on [0, 1], n = 12 (cached).
inline const List<Pair<scalar>>& gl12()
{
    return gaussLegendre01(12);
}


//- The roots in (0, 1) of the quadric phi along the segment p0 -> p1
//  (three evaluations: phi is quadratic in the parameter), appended to r.
template<class PhiAt>
inline void segmentRoots(const point& p0, const point& p1, PhiAt&& phi, DynamicList<scalar>& r)
{
    const scalar f0 = phi(p0);
    const scalar f1 = phi(p1);
    const scalar fm = phi(0.5*(p0 + p1));
    const scalar C = 2*(f1 + f0 - 2*fm);
    const scalar B = f1 - f0 - C;
    const scalar scale = mag(f0) + mag(B) + mag(C);
    if (scale <= VSMALL) return;
    if (mag(C) <= 1e-14*scale)
    {
        if (mag(B) > 1e-14*scale)
        {
            const scalar s = -f0/B;
            if (s > 0 && s < 1) r.append(s);
        }
        return;
    }
    const scalar disc = B*B - 4*C*f0;
    if (disc <= 0) return;
    const scalar sq = Foam::sqrt(disc);
    const scalar qq = -0.5*(B + sign(B)*sq);
    scalar r1 = (qq != 0) ? f0/qq : -B/(2*C);
    scalar r2 = (C != 0) ? qq/C : r1;
    if (r1 > r2) Swap(r1, r2);
    if (r1 > 0 && r1 < 1) r.append(r1);
    if (r2 > 0 && r2 < 1 && r2 != r1) r.append(r2);
}


//- The integral over the triangle (a, b, c) of a point weight W(p) whose
//  kinks lie on the zero sets of up to two quadrics (the traces of the
//  surfaces on the triangle): the triangle is swept by the segments from
//  the points P(t) of the edge a -> b to the vertex c, p(t, s) = P(t) +
//  s (c - P(t)), dA = 2 A (1 - s) dt ds; along every such segment the
//  roots of the quadrics are exact and the s-integral is a 12-point Gauss
//  rule on each piece between them (W smooth there); the t-integral is a
//  12-point Gauss rule refined adaptively (a parent against its two halves)
//  to tolAbs, at most maxDepth levels. The breakpoints functor fills the
//  parameters s in (0, 1) at which the segment P -> c crosses a trace.
template<class W, class Breaks>
inline scalar triangleIntegralSplit
(
    const point& a, const point& b, const point& c,
    W&& w,
    Breaks&& breakpoints,
    const scalar tolAbs,
    const label maxDepth,
    label& nLeaves
)
{
    const scalar A2 = mag((b - a) ^ (c - a));      // 2 A
    if (A2 <= VSMALL) return 0;
    const List<Pair<scalar>>& g = gl12();
    DynamicList<scalar> br(8);

    // the inner integral int_0^1 W(p(t, s)) (1 - s) ds for one t
    auto inner = [&](const scalar t) -> scalar
    {
        const point P = a + t*(b - a);
        br.clear();
        breakpoints(P, c, br);
        // sort the breakpoints
        for (label i = 1; i < br.size(); ++i)
        {
            for (label j = i; j > 0 && br[j] < br[j - 1]; --j) Swap(br[j], br[j - 1]);
        }
        scalar s0 = 0, sum = 0;
        for (label k = 0; k <= br.size(); ++k)
        {
            const scalar s1 = (k < br.size()) ? br[k] : scalar(1);
            if (s1 - s0 > 1e-14)
            {
                for (const Pair<scalar>& gw : g)
                {
                    const scalar s = s0 + gw.first()*(s1 - s0);
                    sum += gw.second()*(s1 - s0)*w(P + s*(c - P))*(1 - s);
                }
            }
            s0 = s1;
        }
        return sum;
    };
    auto gaussT = [&](const scalar t0, const scalar t1) -> scalar
    {
        scalar sum = 0;
        for (const Pair<scalar>& gw : g)
        {
            sum += gw.second()*(t1 - t0)*inner(t0 + gw.first()*(t1 - t0));
        }
        return sum;
    };
    // the outer rule is refined at least minDepth times before the
    // parent-against-halves test may accept: the inner integral is not
    // smooth in t where a trace enters or leaves the sweep (a breakpoint
    // appears), and a coarse agreement there is accidental
    const label minDepth = 3;
    std::function<scalar(scalar, scalar, scalar, scalar, label)> rec =
        [&](scalar t0, scalar t1, scalar parent, scalar tol, label depth) -> scalar
        {
            const scalar tm = 0.5*(t0 + t1);
            const scalar c0 = gaussT(t0, tm);
            const scalar c1 = gaussT(tm, t1);
            const scalar children = c0 + c1;
            if ((depth >= minDepth && mag(children - parent) <= tol) || depth >= maxDepth)
            {
                nLeaves += 2;
                return children;
            }
            return rec(t0, tm, c0, 0.5*tol, depth + 1) + rec(tm, t1, c1, 0.5*tol, depth + 1);
        };
    const scalar parent = gaussT(0, 1);
    return A2*rec(0, 1, parent, tolAbs/A2, 0);
}


// * * * * * * * * * * * * * * The prism volume  * * * * * * * * * * * * * //

//- The six-point degree-4 Gauss rule on the unit triangle (Strang and Fix):
//  barycentric coordinates and weights summing to one.
struct triGauss6
{
    static constexpr label n = 6;
    static inline void point(const label k, scalar& l0, scalar& l1, scalar& l2, scalar& w)
    {
        static const scalar a1 = 0.445948490915965, b1 = 0.108103018168070, w1 = 0.223381589678011;
        static const scalar a2 = 0.091576213509771, b2 = 0.816847572980459, w2 = 0.109951743655322;
        switch (k)
        {
            case 0: l0 = b1; l1 = a1; l2 = a1; w = w1; break;
            case 1: l0 = a1; l1 = b1; l2 = a1; w = w1; break;
            case 2: l0 = a1; l1 = a1; l2 = b1; w = w1; break;
            case 3: l0 = b2; l1 = a2; l2 = a2; w = w2; break;
            case 4: l0 = a2; l1 = b2; l2 = a2; w = w2; break;
            default: l0 = a2; l1 = a2; l2 = b2; w = w2; break;
        }
    }
};


//- The liquid volume of the prism swept by the triangle (a, b, c) along the
//  displacement -d, i.e. between the triangle and its pre-image (a - d,
//  b - d, c - d), under the quadric phi: |d . n_T| times the integral over
//  the triangle of the liquid fraction of the segment [p - d, p], by an
//  adaptive six-point Gauss rule (refined where the parent and its four
//  children disagree by more than tolAbs, at most maxDepth levels). The
//  segment fractions are exact, so the quadrature is the only error. nLeaves
//  counts the evaluated leaves (cost).
template<class PhiAt>
inline scalar prismLiquidVolume
(
    const point& a, const point& b, const point& c,
    const vector& d,
    PhiAt&& phi,
    const scalar tolAbs,
    const label maxDepth,
    label& nLeaves
)
{
    const vector S = 0.5*((b - a) ^ (c - a));
    const scalar thickness = mag(d & S)/max(mag(S), VSMALL);   // |d . n_T|
    const scalar A = mag(S);
    if (A <= VSMALL || thickness <= VSMALL) return 0;

    // the segment fraction's kinks: where the face trace or the pre-image's
    // trace crosses the sweep line (both quadrics in the point)
    auto breaks = [&](const point& P, const point& Q, DynamicList<scalar>& r)
    {
        segmentRoots(P, Q, phi, r);
        segmentRoots(P - d, Q - d, phi, r);
    };
    auto weight = [&](const point& p) -> scalar
    {
        return segmentLiquidFraction(p - d, p, phi);
    };
    // tolAbs is a volume: the area integral's tolerance
    return thickness*triangleIntegralSplit(a, b, c, weight, breaks, tolAbs/thickness, maxDepth, nLeaves);
}


//- The liquid volume flux prism of the face faceI (fan of triangles about
//  the face centre) displaced by -d, under phi.
template<class PhiAt>
inline scalar faceSweptLiquidVolume
(
    const fvMesh& mesh,
    const label faceI,
    const vector& d,
    PhiAt&& phi,
    const scalar tolRel,
    const label maxDepth,
    label& nLeaves
)
{
    const pointField& P = mesh.points();
    const face& f = mesh.faces()[faceI];
    const point& xf = mesh.faceCentres()[faceI];
    const vector& Sf = mesh.faceAreas()[faceI];
    const scalar prismVol = mag(d & Sf);
    if (prismVol <= VSMALL) return 0;
    const scalar tolAbs = tolRel*prismVol;
    scalar v = 0;
    forAll(f, ip)
    {
        const label ip1 = f.nextLabel(ip);
        v += prismLiquidVolume(xf, P[f[ip]], P[ip1], d, phi, tolAbs, maxDepth, nLeaves);
    }
    return v;
}


// * * * * * * * * * The cell volume by rays from the cell centre * * * * * //

//- For the ray x(s) = apex + s (p - apex), s in [0, 1]: the integral of
//  3 s^2 over the liquid part {phi(x(s)) <= 0}, i.e. sum (s_b^3 - s_a^3) over
//  the liquid sub-intervals (exact: phi is quadratic in s). The volume of
//  the pyramid (apex, face) under phi is the integral of this weight over
//  the face times (1/3)|(p - apex) . n_f| dA, the cone volume element.
template<class PhiAt>
inline scalar rayLiquidCubicWeight(const point& apex, const point& p, PhiAt&& phi)
{
    const scalar f0 = phi(apex);
    const scalar f1 = phi(p);
    const scalar fm = phi(0.5*(apex + p));
    const scalar C = 2*(f1 + f0 - 2*fm);
    const scalar B = f1 - f0 - C;
    const scalar scale = mag(f0) + mag(B) + mag(C);
    auto q = [&](const scalar s) { return f0 + B*s + C*s*s; };
    auto cubicLen = [&](const scalar s0, const scalar s1) -> scalar
    {
        const scalar sm = 0.5*(s0 + s1);
        return (q(sm) <= 0) ? (s1*s1*s1 - s0*s0*s0) : 0;
    };
    if (scale <= VSMALL) return 1;
    // the roots in (0, 1)
    scalar r[2];
    label nr = 0;
    if (mag(C) <= 1e-14*scale)
    {
        if (mag(B) > 1e-14*scale)
        {
            const scalar s = -f0/B;
            if (s > 0 && s < 1) r[nr++] = s;
        }
    }
    else
    {
        const scalar disc = B*B - 4*C*f0;
        if (disc > 0)
        {
            const scalar sq = Foam::sqrt(disc);
            const scalar qq = -0.5*(B + sign(B)*sq);
            scalar r1 = (qq != 0) ? f0/qq : -B/(2*C);
            scalar r2 = (C != 0) ? qq/C : r1;
            if (r1 > r2) Swap(r1, r2);
            if (r1 > 0 && r1 < 1) r[nr++] = r1;
            if (r2 > 0 && r2 < 1 && r2 != r1) r[nr++] = r2;
        }
    }
    scalar s0 = 0, w = 0;
    for (label k = 0; k < nr; ++k)
    {
        w += cubicLen(s0, r[k]);
        s0 = r[k];
    }
    w += cubicLen(s0, 1);
    return w;
}


//- The liquid volume of the pyramid (apex; a, b, c) under phi by the
//  adaptive six-point rule over the base triangle (as prismLiquidVolume),
//  the ray weights exact.
template<class PhiAt>
inline scalar pyramidLiquidVolume
(
    const point& apex,
    const point& a, const point& b, const point& c,
    PhiAt&& phi,
    const scalar tolAbs,
    const label maxDepth,
    label& nLeaves
)
{
    const vector S = 0.5*((b - a) ^ (c - a));
    const scalar A = mag(S);
    if (A <= VSMALL) return 0;
    const vector nT = S/A;
    const scalar height = mag((a - apex) & nT);
    if (height <= VSMALL) return 0;
    // (1/3) height * integral of the ray weight over the triangle (the
    // height is the same for every point of the planar triangle); the
    // weight's kink is the surface's trace on the triangle (the ray crosses
    // the surface at its end), the silhouette of the surface seen from the
    // apex is left to the adaptive t-refinement
    auto breaks = [&](const point& P, const point& Q, DynamicList<scalar>& r)
    {
        segmentRoots(P, Q, phi, r);
    };
    auto weight = [&](const point& p) -> scalar
    {
        return rayLiquidCubicWeight(apex, p, phi);
    };
    const scalar factor = height/3.0;
    return factor*triangleIntegralSplit(a, b, c, weight, breaks, tolAbs/factor, maxDepth, nLeaves);
}


//- The liquid volume of cellI under the quadric phi: the pyramids from the
//  cell centre over the fan triangles of every face, by the same adaptive
//  quadrature and tolerance (relative to the cell volume) as the flux
//  prisms, so that the volume anchor and the flux are one quadrature.
template<class PhiAt>
inline scalar cellLiquidVolumeQuadric
(
    const fvMesh& mesh,
    const label cellI,
    PhiAt&& phi,
    const scalar tolRel,
    const label maxDepth,
    label& nLeaves
)
{
    const pointField& P = mesh.points();
    const faceList& faces = mesh.faces();
    const vectorField& Cf = mesh.faceCentres();
    const point& xc = mesh.cellCentres()[cellI];
    const labelList& cFaces = mesh.cells()[cellI];
    const scalar tolAbs = tolRel*mesh.V()[cellI];
    scalar v = 0;
    forAll(cFaces, cf)
    {
        const label faceI = cFaces[cf];
        const face& f = faces[faceI];
        const point& xf = Cf[faceI];
        forAll(f, ip)
        {
            const label ip1 = f.nextLabel(ip);
            v += pyramidLiquidVolume(xc, xf, P[f[ip]], P[ip1], phi, tolAbs, maxDepth, nLeaves);
        }
    }
    return v;
}


// * * * * * * * * * * * * * The cell volume (defect) * * * * * * * * * * * //

//- The liquid volume of cellI under phi with every Detrixhe-Aslam tet
//  (x_c, x_f, p_i, p_{i+1}) split `levels` times into eight sub-tets (the
//  linear interpolant on each). From leia's davofRegeneration.H.
template<class PhiAt>
inline scalar cellLiquidVolumeRefined
(
    const fvMesh& mesh,
    const label cellI,
    PhiAt&& phi,
    const label levels
)
{
    std::function<scalar(const point&, const point&, const point&, const point&, label)> rec =
        [&](const point& a, const point& b, const point& c, const point& d, label lev) -> scalar
        {
            const scalar vol = mag(((b - a) ^ (c - a)) & (d - a))/6.0;
            if (vol <= VSMALL) return 0;
            if (lev == 0)
            {
                return tetNegativeFraction(phi(a), phi(b), phi(c), phi(d))*vol;
            }
            const point ab = 0.5*(a + b), ac = 0.5*(a + c), ad = 0.5*(a + d);
            const point bc = 0.5*(b + c), bd = 0.5*(b + d), cd = 0.5*(c + d);
            return
                rec(a, ab, ac, ad, lev - 1) + rec(ab, b, bc, bd, lev - 1)
              + rec(ac, bc, c, cd, lev - 1) + rec(ad, bd, cd, d, lev - 1)
              + rec(ab, ac, ad, bd, lev - 1) + rec(ab, ac, bc, bd, lev - 1)
              + rec(ac, ad, bd, cd, lev - 1) + rec(ac, bc, bd, cd, lev - 1);
        };

    const pointField& P = mesh.points();
    const faceList& faces = mesh.faces();
    const vectorField& Cf = mesh.faceCentres();
    const point& xc = mesh.cellCentres()[cellI];
    const labelList& cFaces = mesh.cells()[cellI];
    scalar liquidVol = 0;
    forAll(cFaces, cf)
    {
        const label faceI = cFaces[cf];
        const face& f = faces[faceI];
        const point& xf = Cf[faceI];
        forAll(f, ip)
        {
            const label ip1 = f.nextLabel(ip);
            liquidVol += rec(xc, xf, P[f[ip]], P[ip1], levels);
        }
    }
    return liquidVol;
}


// * * * * * * * * * * * * * * The offset solve  * * * * * * * * * * * * * //

//- Solve g(c) = target for the offset c by bisection; g is increasing in c
//  (a larger c puts more of the cell into the liquid). The bracket [lo, hi]
//  is widened up to eight times. Returns |g(c) - target|; bracketed is
//  false when no bracket was found. From leia's davofRegeneration.H.
template<class G>
inline scalar solveOffset
(
    G&& g,
    const scalar target,
    scalar lo,
    scalar hi,
    const scalar tol,
    scalar& c,
    bool& bracketed
)
{
    scalar glo = g(lo) - target, ghi = g(hi) - target;
    for (label k = 0; k < 8 && glo > 0; ++k)
    {
        lo -= (hi - lo);
        glo = g(lo) - target;
    }
    for (label k = 0; k < 8 && ghi < 0; ++k)
    {
        hi += (hi - lo);
        ghi = g(hi) - target;
    }
    bracketed = !(glo > 0 || ghi < 0);
    if (!bracketed)
    {
        c = 0.5*(lo + hi);
        return min(mag(glo), mag(ghi));
    }
    scalar mid = 0.5*(lo + hi), gm = 0;
    for (label it = 0; it < 400; ++it)
    {
        mid = 0.5*(lo + hi);
        gm = g(mid) - target;
        if (gm > 0) hi = mid; else lo = mid;
        if (hi - lo <= tol) break;
    }
    c = 0.5*(lo + hi);
    return mag(g(c) - target);
}


// * * * * * * * * * * * * The six-coefficient fit * * * * * * * * * * * * //

struct paraboloidFitResult
{
    bool converged = false;
    bool singular = false;
    bool stalled = false;
    label nIter = 0;
    scalar rmax = 0;        // the max face residual (as a fraction) of the kept surface
};


//- Gauss-Newton on the six dimensionless unknowns y = (c/h, a1, a2, b11 h,
//  b12 h, b22 h) of S so that the paraboloid cuts every face of the cell
//  at the target fractions tF (weights wF = |S_f| / sum |S_f|); the ring
//  tensor and the start offset are the start and, through the Marquardt
//  damping (1 + wPrior^2 on the diagonal plus a floor), the prior of the
//  under-determined directions. fractionsOf(S, aF) fills the fractions of
//  the cell's faces under S in the cell's face order. The best iterate is
//  kept. From leia's leiaTestDavofNormal (regenArm paraboloidFull,
//  2026-10-12).
template<class FractionsOf>
inline paraboloidFitResult fitParaboloidFull
(
    cellSurface& S,
    const scalar h,
    const scalarList& tF,
    const scalarList& wF,
    FractionsOf&& fractionsOf,
    const scalar wPrior,
    const scalar tolFull,
    const label maxIter
)
{
    paraboloidFitResult res;
    const label F = tF.size();
    const scalar eps = 1e-6;
    const label nU = 6;
    FixedList<scalar, 6> y;
    y[0] = S.c/h; y[1] = S.a1; y[2] = S.a2;
    y[3] = S.b11*h; y[4] = S.b12*h; y[5] = S.b22*h;
    auto surfaceOf = [&](const FixedList<scalar, 6>& yy) -> cellSurface
    {
        cellSurface T = S;
        T.paraboloid = true;
        T.c = yy[0]*h;
        T.a1 = yy[1];
        T.a2 = yy[2];
        T.b11 = yy[3]/h;
        T.b12 = yy[4]/h;
        T.b22 = yy[5]/h;
        return T;
    };
    auto residual = [&](const FixedList<scalar, 6>& yy, scalarList& r) -> scalar
    {
        const cellSurface T = surfaceOf(yy);
        scalarList aF;
        fractionsOf(T, aF);
        r.setSize(F);
        scalar rmax = 0;
        forAll(aF, i)
        {
            r[i] = wF[i]*(aF[i] - tF[i]);
            rmax = max(rmax, mag(aF[i] - tF[i]));
        }
        return rmax;
    };
    scalarList r;
    scalar rmax = residual(y, r);
    FixedList<scalar, 6> yBest(y);
    scalar rBest = rmax;
    bool ok = (rmax <= tolFull);
    for (label it = 0; it < maxIter && !ok; ++it)
    {
        ++res.nIter;
        List<scalarList> J(nU);
        for (label dd = 0; dd < nU; ++dd)
        {
            FixedList<scalar, 6> yd(y);
            yd[dd] += eps;
            scalarList rd;
            residual(yd, rd);
            J[dd].setSize(F);
            forAll(rd, i) J[dd][i] = (rd[i] - r[i])/eps;
        }
        FixedList<FixedList<scalar, 6>, 6> A;
        FixedList<scalar, 6> g;
        for (label p = 0; p < nU; ++p)
        {
            for (label q = 0; q < nU; ++q)
            {
                scalar s = 0;
                for (label i = 0; i < F; ++i) s += J[p][i]*J[q][i];
                A[p][q] = s;
            }
            scalar s = 0;
            for (label i = 0; i < F; ++i) s += J[p][i]*r[i];
            g[p] = -s;
        }
        scalar diagMax = 0;
        for (label p = 0; p < nU; ++p) diagMax = max(diagMax, mag(A[p][p]));
        if (diagMax <= VSMALL) { res.singular = true; break; }
        for (label p = 0; p < nU; ++p)
        {
            A[p][p] = A[p][p]*(1 + sqr(wPrior)) + 1e-10*diagMax;
        }
        bool singular = false;
        for (label p = 0; p < nU; ++p)
        {
            label piv = p;
            for (label i = p + 1; i < nU; ++i)
            {
                if (mag(A[i][p]) > mag(A[piv][p])) piv = i;
            }
            if (mag(A[piv][p]) <= 1e-13*diagMax) { singular = true; break; }
            if (piv != p)
            {
                for (label q = 0; q < nU; ++q) Swap(A[p][q], A[piv][q]);
                Swap(g[p], g[piv]);
            }
            for (label i = p + 1; i < nU; ++i)
            {
                const scalar fct = A[i][p]/A[p][p];
                if (fct == 0) continue;
                for (label q = p; q < nU; ++q) A[i][q] -= fct*A[p][q];
                g[i] -= fct*g[p];
            }
        }
        if (singular) { res.singular = true; break; }
        FixedList<scalar, 6> dy;
        for (label p = nU - 1; p >= 0; --p)
        {
            scalar s = g[p];
            for (label q = p + 1; q < nU; ++q) s -= A[p][q]*dy[q];
            dy[p] = s/A[p][p];
        }
        scalar dn = 0;
        forAll(dy, dd) dn += sqr(dy[dd]);
        dn = Foam::sqrt(dn);
        if (dn > 1) forAll(dy, dd) dy[dd] /= dn;
        bool improved = false;
        for (label half = 0; half < 7; ++half)
        {
            FixedList<scalar, 6> yTry(y);
            forAll(yTry, dd) yTry[dd] += dy[dd];
            scalarList rTry;
            const scalar rmaxTry = residual(yTry, rTry);
            if (rmaxTry < rmax)
            {
                y = yTry; r = rTry; rmax = rmaxTry;
                improved = true;
                break;
            }
            forAll(dy, dd) dy[dd] *= 0.5;
        }
        if (!improved) { res.stalled = true; break; }
        if (rmax < rBest) { yBest = y; rBest = rmax; }
        ok = (rmax <= tolFull);
    }
    S = surfaceOf(yBest);
    res.rmax = rBest;
    res.converged = ok;
    return res;
}

} // End namespace davof
} // End namespace Foam

#endif

// ************************************************************************* //
