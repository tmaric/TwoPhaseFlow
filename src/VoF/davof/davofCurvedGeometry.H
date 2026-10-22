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
        liqA += triLiquidArea(xf, p0, p1, phi, 3, nFallback);
        totA += triA;
    }
    return (totA > VSMALL) ? liqA/totA : (phi(xf) < 0 ? scalar(1) : scalar(0));
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

    // the integral of the segment fraction over a triangle by the six-point
    // rule, as an area integral
    auto gauss = [&](const point& p0, const point& p1, const point& p2) -> scalar
    {
        const scalar At = 0.5*mag((p1 - p0) ^ (p2 - p0));
        scalar s = 0;
        for (label k = 0; k < triGauss6::n; ++k)
        {
            scalar l0, l1, l2, w;
            triGauss6::point(k, l0, l1, l2, w);
            const point p = l0*p0 + l1*p1 + l2*p2;
            s += w*segmentLiquidFraction(p - d, p, phi);
        }
        return s*At;
    };

    std::function<scalar(const point&, const point&, const point&, scalar, scalar, label)> rec =
        [&](const point& p0, const point& p1, const point& p2, scalar parent, scalar tol, label depth) -> scalar
        {
            const point m01 = 0.5*(p0 + p1), m12 = 0.5*(p1 + p2), m20 = 0.5*(p2 + p0);
            const scalar c0 = gauss(p0, m01, m20);
            const scalar c1 = gauss(m01, p1, m12);
            const scalar c2 = gauss(m20, m12, p2);
            const scalar c3 = gauss(m01, m12, m20);
            const scalar children = c0 + c1 + c2 + c3;
            if (mag(children - parent) <= tol || depth >= maxDepth)
            {
                nLeaves += 4;
                return children;
            }
            const scalar tc = 0.25*tol;
            return
                rec(p0, m01, m20, c0, tc, depth + 1)
              + rec(m01, p1, m12, c1, tc, depth + 1)
              + rec(m20, m12, p2, c2, tc, depth + 1)
              + rec(m01, m12, m20, c3, tc, depth + 1);
        };

    const scalar parent = gauss(a, b, c);
    // tolAbs is a volume: convert to the area integral's tolerance
    const scalar tolArea = tolAbs/thickness;
    return thickness*rec(a, b, c, parent, tolArea, 0);
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
