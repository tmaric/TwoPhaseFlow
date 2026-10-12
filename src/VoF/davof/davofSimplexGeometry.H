// copied from leia src/leiaLevelSet/davof, commit f13a529
/*---------------------------------------------------------------------------*\
    leia :: davof :: simplex geometry

    Closed-form liquid fractions and zero-set patches of a LINEAR field
    phi on a triangle and on a tetrahedron. Liquid is Omega^- = {phi <= 0}.

    The fraction formulas are copies of the ones leia already uses:
      * triNegativeFraction: applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/
        faceAreaFraction.H (exact for a linear phi on a triangle);
      * tetNegativeFraction: src/leiaLevelSet/phaseIndicator/detrixheAslamPhaseIndicator.C
        (Detrixhe & Aslam 2016, Eq. 7, with the corrected two-in-two-out case).
    They are copied, not moved, so that no existing code path changes.

    tetZeroSet is new: it returns the zero set of the linear interpolant of the
    four vertex values as a polygon (triangle for a 1-3 sign split, quad for a
    2-2 split) through its area vector, its area-weighted centroid and its area.
    The area vector points from Omega^- to Omega^+ (along the gradient of the
    interpolant). Summed over the tets of a cell it equals the DAVOF area normal
    m_c = -sum_f alpha_f S_f exactly (up to round-off) when the face fractions
    alpha_f come from the SAME vertex values on the face triangles (the tet
    (xc, xf, p0, p1) restricted to its outer face triangle (xf, p0, p1) IS the
    triangle's linear interpolant, and the inner tet faces cancel pairwise).
\*---------------------------------------------------------------------------*/

#ifndef davofSimplexGeometry_H
#define davofSimplexGeometry_H

#include "scalar.H"
#include "vector.H"
#include "point.H"
#include "label.H"
#include <algorithm>

namespace Foam
{
namespace davof
{

//- Fraction of a triangle in Omega^- = {phi <= 0} from its three vertex values.
inline scalar triNegativeFraction(scalar a, scalar b, scalar c)
{
    scalar p[3] = {a, b, c};
    std::sort(p, p + 3);
    const scalar p1 = p[0], p2 = p[1], p3 = p[2];

    if (p3 <= 0) return 1;                          // all three inside
    if (p1 >  0) return 0;                          // all three outside
    if (p2 <= 0)                                    // two inside, one out (corner p3)
    {
        return 1 - (p3*p3)/stabilise((p3 - p1)*(p3 - p2), VSMALL);
    }
    // one inside (p1), two outside
    return (p1*p1)/stabilise((p2 - p1)*(p3 - p1), VSMALL);
}


//- Fraction of a tetrahedron in Omega^- = {phi <= 0} from its four vertex
//  values (Detrixhe & Aslam 2016, Eq. 7; the middle case corrected as in
//  detrixheAslamPhaseIndicator.C).
inline scalar tetNegativeFraction(scalar a, scalar b, scalar c, scalar d)
{
    scalar p[4] = {a, b, c, d};
    std::sort(p, p + 4);
    const scalar p1 = p[0], p2 = p[1], p3 = p[2], p4 = p[3];

    if (p4 <= 0) return 1;          // all four inside  (Omega^-)
    if (p1 >  0) return 0;          // all four outside (Omega^+)

    if (p3 <= 0)                    // three inside, one outside (corner at p4)
    {
        const scalar den = (p4 - p1)*(p4 - p2)*(p4 - p3);
        return 1 - (p4*p4*p4)/stabilise(den, VSMALL);
    }
    if (p2 <= 0)                    // two inside, two outside (middle)
    {
        const scalar num =
            p3*p4*(p1*p1 + p1*p2 + p2*p2)
          - p1*p2*(p1 + p2)*(p3 + p4)
          + p1*p1*p2*p2;
        const scalar den = (p3 - p1)*(p4 - p1)*(p3 - p2)*(p4 - p2);
        return num/stabilise(den, VSMALL);
    }
    // one inside, three outside (corner at p1)
    const scalar den = (p2 - p1)*(p3 - p1)*(p4 - p1);
    return -(p1*p1*p1)/stabilise(den, VSMALL);
}


//- Point on the edge (xi, xj) where the linear interpolant of (phii, phij)
//  vanishes. The caller guarantees phii and phij have different signs
//  (one <= 0, the other > 0), so the denominator is non-zero.
inline point edgeCrossing
(
    const point& xi, const point& xj, const scalar phii, const scalar phij
)
{
    return xi + (phii/(phii - phij))*(xj - xi);
}


//- The part of the triangle (x0, x1, x2) in Omega^- = {phi <= 0} under the
//  linear interpolant of its vertex values (a, b, c): the area vector Swet
//  (oriented like (x1 - x0) ^ (x2 - x0)) and the area-weighted centroid xwet
//  of that part, built by walking the triangle boundary (wet vertices and
//  the edge crossings). Returns the area fraction; xwet is the triangle
//  centroid when nothing or everything is wet. Convention of
//  triNegativeFraction (phi <= 0 is wet).
inline scalar triNegativePolygon
(
    const point& x0, const point& x1, const point& x2,
    const scalar a, const scalar b, const scalar c,
    vector& Swet, point& xwet
)
{
    const point x[3] = {x0, x1, x2};
    const scalar ph[3] = {a, b, c};
    const vector St = 0.5*((x1 - x0) ^ (x2 - x0));
    const scalar At = mag(St);
    Swet = Zero;
    xwet = (x0 + x1 + x2)/3.0;
    if (At <= VSMALL) return 0;
    const vector nt = St/At;

    point poly[6];
    label n = 0;
    for (label i = 0; i < 3; ++i)
    {
        const label j = (i + 1) % 3;
        if (ph[i] <= 0) poly[n++] = x[i];
        if ((ph[i] <= 0) != (ph[j] <= 0))
        {
            poly[n++] = edgeCrossing(x[i], x[j], ph[i], ph[j]);
        }
    }
    if (n < 3) return 0;

    // Fan about poly[0]: signed areas along the triangle normal.
    vector S(Zero), xs(Zero);
    scalar A = 0;
    for (label k = 1; k + 1 < n; ++k)
    {
        const vector Sk = 0.5*((poly[k] - poly[0]) ^ (poly[k + 1] - poly[0]));
        const scalar Ak = Sk & nt;
        S += Sk;
        A += Ak;
        xs += Ak*((poly[0] + poly[k] + poly[k + 1])/3.0);
    }
    Swet = S;
    if (mag(A) > VSMALL) xwet = xs/A;
    return (S & nt)/At;
}


//- Zero set of the linear interpolant of phi on the tet (x[0..3]).
//  areaVec: area vector oriented from Omega^- (phi <= 0) to Omega^+;
//  centroid: area-weighted centroid of the polygon; area: its area.
//  Returns 0 (no cut: all four values on one side), 3 (triangle) or 4 (quad).
inline label tetZeroSet
(
    const point x[4],
    const scalar phi[4],
    vector& areaVec,
    point& centroid,
    scalar& area
)
{
    areaVec = Zero;
    centroid = Zero;
    area = 0;

    label in[4], out[4];
    label nIn = 0, nOut = 0;
    for (label i = 0; i < 4; ++i)
    {
        if (phi[i] <= 0) in[nIn++] = i;
        else out[nOut++] = i;
    }
    if (nIn == 0 || nOut == 0) return 0;

    // Gradient of the linear interpolant: the orientation reference
    // (points from Omega^- to Omega^+). g & e_i = phi_i - phi_0.
    const vector e1 = x[1] - x[0];
    const vector e2 = x[2] - x[0];
    const vector e3 = x[3] - x[0];
    const scalar sixV = e1 & (e2 ^ e3);
    vector g =
        (phi[1] - phi[0])*(e2 ^ e3)
      + (phi[2] - phi[0])*(e3 ^ e1)
      + (phi[3] - phi[0])*(e1 ^ e2);
    if (mag(sixV) > VSMALL) g /= sixV;

    if (nIn == 1 || nOut == 1)
    {
        // 1-3 split: the lone vertex a against the other three.
        const label a = (nIn == 1) ? in[0] : out[0];
        const label* o = (nIn == 1) ? out : in;
        const point p0 = edgeCrossing(x[a], x[o[0]], phi[a], phi[o[0]]);
        const point p1 = edgeCrossing(x[a], x[o[1]], phi[a], phi[o[1]]);
        const point p2 = edgeCrossing(x[a], x[o[2]], phi[a], phi[o[2]]);
        areaVec = 0.5*((p1 - p0) ^ (p2 - p0));
        if ((areaVec & g) < 0) areaVec = -areaVec;
        area = mag(areaVec);
        centroid = (p0 + p1 + p2)/3.0;
        return 3;
    }

    // 2-2 split: {a, b} on one side, {c, d} on the other. The crossings in the
    // cyclic order (ac, ad, bd, bc) are consecutive on the tet faces acd, abd,
    // bcd, abc, so the quad is simple; it is planar and convex because it is
    // the intersection of a plane with a convex tet.
    const label a = in[0], b = in[1], c = out[0], d = out[1];
    const point pac = edgeCrossing(x[a], x[c], phi[a], phi[c]);
    const point pad = edgeCrossing(x[a], x[d], phi[a], phi[d]);
    const point pbd = edgeCrossing(x[b], x[d], phi[b], phi[d]);
    const point pbc = edgeCrossing(x[b], x[c], phi[b], phi[c]);
    const vector a1 = 0.5*((pad - pac) ^ (pbd - pac));
    const vector a2 = 0.5*((pbd - pac) ^ (pbc - pac));
    areaVec = a1 + a2;
    if ((areaVec & g) < 0) areaVec = -areaVec;
    area = mag(areaVec);
    const scalar A1 = mag(a1), A2 = mag(a2);
    const point c1 = (pac + pad + pbd)/3.0;
    const point c2 = (pac + pbd + pbc)/3.0;
    centroid =
        (A1 + A2 > VSMALL)
      ? (A1*c1 + A2*c2)/(A1 + A2)
      : 0.25*(pac + pad + pbd + pbc);
    return 4;
}

} // End namespace davof
} // End namespace Foam

#endif

// ************************************************************************* //
