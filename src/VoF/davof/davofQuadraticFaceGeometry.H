// copied from leia src/leiaLevelSet/davof, commit f13a529
/*---------------------------------------------------------------------------*\
    leia :: davof :: quadratic face fractions

    The liquid-area fraction of a face triangle (x0, x1, x2) from the EXACT
    signed distance at its three vertices and its three edge midpoints: the
    quadratic Lagrange interpolant q of those six values is cut through the
    triangle and the area of {q <= 0} is integrated exactly up to a
    Gauss-Legendre quadrature of the sliver between the chord of the two
    boundary crossings and the conic arc (the integrand is analytic, 12 points
    are at round-off for a resolved interface).

    Why: the Detrixhe-Aslam formulas are exact for a LINEAR function on each
    simplex, so with exact vertex values they integrate the chord surface,
    whose crossings are O(h^2) off the true interface; that leaves every cut
    face with an O(h) error as a fraction and the Gauss-identity normal first
    order (measured, config/davof/sphereNormal3D.yaml). One polynomial degree
    more moves the crossings to O(h^3), the fractions to O(h^2) and the normal
    to second order, from exact signed distances alone, without iteration and
    without knowing the surface type.

    Scope: one interface crossing per triangle (two boundary crossings of
    q = 0, the arc a graph over its chord). Anything else -- four or six
    crossings, a closed zero set inside the triangle, an arc that folds back
    -- sets `fallback` and the caller keeps the Detrixhe-Aslam value; the
    count is reported.
\*---------------------------------------------------------------------------*/

#ifndef davofQuadraticFaceGeometry_H
#define davofQuadraticFaceGeometry_H

#include "scalar.H"
#include "vector.H"
#include "vector2D.H"
#include "point.H"
#include "List.H"
#include "Pair.H"
#include "mathematicalConstants.H"

namespace Foam
{
namespace davof
{

//- Gauss-Legendre nodes and weights on [0, 1] (Newton on P_n), cached.
inline const List<Pair<scalar>>& gaussLegendre01(const label n = 12)
{
    static List<Pair<scalar>> gl;
    if (gl.size() != n)
    {
        gl.setSize(n);
        for (label i = 0; i < n; ++i)
        {
            scalar x = Foam::cos(constant::mathematical::pi*(i + 0.75)/(n + 0.5));
            scalar dp = 1;
            for (label it = 0; it < 100; ++it)
            {
                scalar p0 = 1, p1 = x;
                for (label k = 2; k <= n; ++k)
                {
                    const scalar p2 = ((2*k - 1)*x*p1 - (k - 1)*p0)/k;
                    p0 = p1; p1 = p2;
                }
                dp = n*(x*p1 - p0)/(x*x - 1);
                const scalar dx = p1/dp;
                x -= dx;
                if (mag(dx) < 1e-15) break;
            }
            gl[i] = Pair<scalar>(0.5*(x + 1), 1.0/((1 - x*x)*dp*dp));
        }
    }
    return gl;
}

inline scalar cross2(const vector2D& a, const vector2D& b)
{
    return a.x()*b.y() - a.y()*b.x();
}

//- Quadratic Lagrange interpolant on a triangle: phi[0..2] at the vertices,
//  phi[3] on edge 01, phi[4] on edge 12, phi[5] on edge 20.
inline scalar quadTri(const scalar l0, const scalar l1, const scalar l2, const scalar phi[6])
{
    return
        phi[0]*l0*(2*l0 - 1) + phi[1]*l1*(2*l1 - 1) + phi[2]*l2*(2*l2 - 1)
      + 4*phi[3]*l0*l1 + 4*phi[4]*l1*l2 + 4*phi[5]*l2*l0;
}

//- Roots in (0, 1) of the quadratic interpolant along an edge with end
//  values a, b and midpoint value m.
inline label edgeRoots(const scalar a, const scalar b, const scalar m, scalar t[2])
{
    const scalar A = 2*a + 2*b - 4*m, B = -3*a - b + 4*m, C = a;
    const scalar tol = 1e-12;
    label n = 0;
    if (mag(A) <= 1e-14*(mag(B) + mag(C) + VSMALL))
    {
        if (mag(B) > VSMALL)
        {
            const scalar r = -C/B;
            if (r > tol && r < 1 - tol) t[n++] = r;
        }
        return n;
    }
    const scalar disc = B*B - 4*A*C;
    if (disc < 0) return 0;
    // The numerically stable pair of roots (no cancellation for a nearly
    // linear function, i.e. a plane, where A is round-off).
    const scalar sq = Foam::sqrt(disc);
    const scalar qq = -0.5*(B + ((B >= 0) ? sq : -sq));
    scalar r1, r2;
    if (mag(qq) > VSMALL)
    {
        r1 = qq/A;
        r2 = C/qq;
    }
    else
    {
        r1 = r2 = 0;
    }
    const scalar lo = min(r1, r2), hi = max(r1, r2);
    if (lo > tol && lo < 1 - tol) t[n++] = lo;
    if (hi > tol && hi < 1 - tol && hi - lo > tol) t[n++] = hi;
    return n;
}

//- Liquid (phi <= 0) area fraction of the triangle (x0, x1, x2) under the
//  quadratic interpolant of phi[6] (vertices 0, 1, 2; midpoints 01, 12, 20).
inline scalar triQuadraticNegativeFraction
(
    const point& x0, const point& x1, const point& x2,
    const scalar phi[6],
    bool& fallback
)
{
    fallback = false;
    bool allNeg = true, allPos = true;
    for (label i = 0; i < 6; ++i)
    {
        if (phi[i] > 0) allNeg = false;
        if (phi[i] < 0) allPos = false;
        if (phi[i] == 0) { fallback = true; return 0; }
    }
    if (allNeg) return 1;
    if (allPos) return 0;

    // The triangle's 2D frame.
    const vector d1 = x1 - x0, d2 = x2 - x0;
    const vector nrm = d1 ^ d2;
    const scalar A2 = mag(nrm);
    if (A2 <= VSMALL) { fallback = true; return 0; }
    const vector e1 = d1/mag(d1);
    const vector e2 = (nrm/A2) ^ e1;
    const vector2D V[3] =
    {
        vector2D(0, 0),
        vector2D(mag(d1), 0),
        vector2D(d2 & e1, d2 & e2)
    };
    const scalar D = cross2(V[1] - V[0], V[2] - V[0]);   // = 2 x area > 0
    auto q = [&](const vector2D& p) -> scalar
    {
        const scalar l1 = cross2(p - V[0], V[2] - V[0])/D;
        const scalar l2 = cross2(V[1] - V[0], p - V[0])/D;
        return quadTri(1 - l1 - l2, l1, l2, phi);
    };

    // Boundary crossings, ordered along the boundary (edge index, then t).
    const label ea[3] = {0, 1, 2}, eb[3] = {1, 2, 0}, em[3] = {3, 4, 5};
    label cEdge[6]; vector2D cP[6];
    label nc = 0;
    for (label e = 0; e < 3; ++e)
    {
        scalar t[2];
        const label nr = edgeRoots(phi[ea[e]], phi[eb[e]], phi[em[e]], t);
        for (label k = 0; k < nr; ++k)
        {
            cEdge[nc] = e;
            cP[nc] = V[ea[e]] + t[k]*(V[eb[e]] - V[ea[e]]);
            ++nc;
        }
    }
    if (nc != 2) { fallback = true; return 0; }

    // The forward polygon: from crossing A along the boundary to crossing B.
    const vector2D& Ap = cP[0];
    const vector2D& Bp = cP[1];
    List<vector2D> fwd;
    fwd.append(Ap);
    vector2D probe;    // a boundary point just after A on the walk
    if (cEdge[0] == cEdge[1])
    {
        probe = 0.5*(Ap + Bp);
    }
    else
    {
        label e = cEdge[0];
        probe = 0.5*(Ap + V[eb[e]]);
        while (e != cEdge[1])
        {
            fwd.append(V[eb[e]]);
            e = (e + 1) % 3;
        }
    }
    fwd.append(Bp);
    scalar areaFwd = 0;
    forAll(fwd, i)
    {
        areaFwd += cross2(fwd[i], fwd[(i + 1) % fwd.size()]);
    }
    areaFwd = 0.5*mag(areaFwd);

    // Chord normal pointing to the forward polygon's side.
    const vector2D d = Bp - Ap;
    const scalar dm = mag(d);
    if (dm <= VSMALL) { fallback = true; return 0; }
    vector2D n(-d.y()/dm, d.x()/dm);
    if (cEdge[0] == cEdge[1])
    {
        // The chord lies on an edge: the forward side is outside the triangle.
        const vector2D toThird = V[(eb[cEdge[0]] + 1) % 3] - Ap;
        if ((n & toThird) > 0) n = -n;
    }
    else
    {
        const vector2D toV = fwd[1] - Ap;
        if ((n & toV) < 0) n = -n;
    }

    // The sliver between chord and arc, signed towards the forward side.
    const List<Pair<scalar>>& gl = gaussLegendre01();
    scalar sliver = 0;
    forAll(gl, g)
    {
        const vector2D c = Ap + gl[g].first()*d;
        const scalar c0 = q(c);
        const scalar qp = q(c + dm*n), qm = q(c - dm*n);
        const scalar b = (qp - qm)/(2*dm);
        const scalar a = (qp + qm - 2*c0)/(2*dm*dm);
        scalar s;
        if (mag(a) <= 1e-14*(mag(b) + VSMALL)/dm)
        {
            if (mag(b) <= VSMALL) { fallback = true; return 0; }
            s = -c0/b;
        }
        else
        {
            const scalar disc = b*b - 4*a*c0;
            if (disc < 0) { fallback = true; return 0; }
            // The root of smaller magnitude, in the stable form c0/qq (it
            // tends to -c0/b as a -> 0 instead of cancelling).
            const scalar sq = Foam::sqrt(disc);
            const scalar qq = -0.5*(b + ((b >= 0) ? sq : -sq));
            if (mag(qq) <= VSMALL) { fallback = true; return 0; }
            s = c0/qq;
        }
        if (mag(s) > dm) { fallback = true; return 0; }   // not a graph over the chord
        sliver += gl[g].second()*s;
    }
    sliver *= dm;

    // Which side is liquid (q <= 0)?
    const scalar triArea = 0.5*D;
    scalar liquid;
    if (q(probe) < 0)
    {
        liquid = areaFwd - sliver;
    }
    else
    {
        liquid = (triArea - areaFwd) + sliver;
    }
    return min(max(liquid/triArea, scalar(0)), scalar(1));
}

} // End namespace davof
} // End namespace Foam

#endif

// ************************************************************************* //
