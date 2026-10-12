// copied from leia src/leiaLevelSet/davof, commit f13a529
/*---------------------------------------------------------------------------*\
    leia :: davof :: exact sphere geometry

    Closed-form liquid area of a PLANAR polygon face inside a ball (the face
    plane cuts the ball in a circle; the area of a polygon inside a disk is the
    classical sum of signed edge contributions, each a triangle, a circular
    sector or a mix), and the exact liquid volume of an AXIS-ALIGNED box inside
    a ball (adaptive Simpson quadrature over z of the exact rectangle-disk
    cross-section area, to a relative tolerance of 1e-12).

    Together they give the exact DAVOF state (alpha_c, alpha_f) of a sphere on a
    blockMesh hex mesh: the state the reviewers of the DAVOF proposal asked to
    see the Gauss-identity normal recovered from, and the reference against
    which the linear-interpolant state's O(h) face-fraction error is measured.
\*---------------------------------------------------------------------------*/

#ifndef davofSphereGeometry_H
#define davofSphereGeometry_H

#include "scalar.H"
#include "vector.H"
#include "vector2D.H"
#include "point.H"
#include "face.H"
#include "pointField.H"
#include "List.H"

namespace Foam
{
namespace davof
{

inline scalar cross2D(const vector2D& a, const vector2D& b)
{
    return a.x()*b.y() - a.y()*b.x();
}

//- Signed area of the circular sector (radius r, centre at the origin) swept
//  from direction a to direction b.
inline scalar sector2D(const vector2D& a, const vector2D& b, const scalar r)
{
    return 0.5*r*r*Foam::atan2(cross2D(a, b), a & b);
}

//- Signed area of the triangle (origin, a, b).
inline scalar tri2D(const vector2D& a, const vector2D& b)
{
    return 0.5*cross2D(a, b);
}

//- Signed area of the intersection of the triangle (origin, a, b) with the
//  disk of radius r centred at the origin.
inline scalar edgeDiskArea(const vector2D& a, const vector2D& b, const scalar r)
{
    const scalar da = mag(a), db = mag(b);
    if (da <= r && db <= r) return tri2D(a, b);

    const vector2D d = b - a;
    const scalar A = d & d;
    if (A <= VSMALL) return 0;
    const scalar B = 2*(a & d);
    const scalar C = (a & a) - r*r;
    const scalar disc = B*B - 4*A*C;
    if (disc <= 0) return sector2D(a, b, r);           // the edge misses the disk

    const scalar sq = Foam::sqrt(disc);
    const scalar t1 = (-B - sq)/(2*A), t2 = (-B + sq)/(2*A);
    if (da <= r)                                        // a inside, exit at t2
    {
        const vector2D p = a + t2*d;
        return tri2D(a, p) + sector2D(p, b, r);
    }
    if (db <= r)                                        // enter at t1, b inside
    {
        const vector2D p = a + t1*d;
        return sector2D(a, p, r) + tri2D(p, b);
    }
    if (t1 > 0 && t2 < 1 && t1 < t2)                    // both outside, crossing twice
    {
        const vector2D p1 = a + t1*d, p2 = a + t2*d;
        return sector2D(a, p1, r) + tri2D(p1, p2) + sector2D(p2, b, r);
    }
    return sector2D(a, b, r);
}

//- Area of the polygon (2D, any orientation) inside the disk of radius r
//  centred at the origin.
inline scalar polygonDiskArea(const List<vector2D>& poly, const scalar r)
{
    scalar s = 0;
    forAll(poly, i)
    {
        s += edgeDiskArea(poly[i], poly[(i + 1) % poly.size()], r);
    }
    return mag(s);
}

//- Area of the planar face f (centre xf, unit normal nf) inside the ball
//  (centre c, radius R).
inline scalar faceBallArea
(
    const face& f,
    const pointField& points,
    const point& xf,
    const vector& nf,
    const point& c,
    const scalar R
)
{
    const scalar d = (c - xf) & nf;                     // centre to the face plane
    if (mag(d) >= R) return 0;
    const scalar r = Foam::sqrt(R*R - d*d);             // circle in the face plane
    const point cp = c - d*nf;                          // its centre

    vector e1(Zero);
    forAll(f, i)
    {
        vector t = points[f[i]] - xf;
        t -= (t & nf)*nf;
        const scalar tm = mag(t);
        if (tm > VSMALL) { e1 = t/tm; break; }
    }
    if (mag(e1) < VSMALL) return 0;
    const vector e2 = nf ^ e1;

    List<vector2D> poly(f.size());
    forAll(f, i)
    {
        const vector q = points[f[i]] - cp;
        poly[i] = vector2D(q & e1, q & e2);
    }
    return polygonDiskArea(poly, r);
}

//- Area of the rectangle [x0,x1] x [y0,y1] (relative to the disk centre) inside
//  the disk of radius r.
inline scalar rectDiskArea
(
    const scalar x0, const scalar x1, const scalar y0, const scalar y1, const scalar r
)
{
    List<vector2D> poly(4);
    poly[0] = vector2D(x0, y0);
    poly[1] = vector2D(x1, y0);
    poly[2] = vector2D(x1, y1);
    poly[3] = vector2D(x0, y1);
    return polygonDiskArea(poly, r);
}

template<class F>
inline scalar adaptiveSimpsonRec
(
    const F& f, const scalar a, const scalar b,
    const scalar fa, const scalar fm, const scalar fb, const scalar whole,
    const scalar tol, const label depth
)
{
    const scalar m = 0.5*(a + b), lm = 0.5*(a + m), rm = 0.5*(m + b);
    const scalar flm = f(lm), frm = f(rm);
    const scalar left = (m - a)/6.0*(fa + 4*flm + fm);
    const scalar right = (b - m)/6.0*(fm + 4*frm + fb);
    const scalar delta = left + right - whole;
    if (depth <= 0 || mag(delta) <= 15*tol)
    {
        return left + right + delta/15.0;
    }
    return
        adaptiveSimpsonRec(f, a, m, fa, flm, fm, left, 0.5*tol, depth - 1)
      + adaptiveSimpsonRec(f, m, b, fm, frm, fb, right, 0.5*tol, depth - 1);
}

//- Adaptive Simpson quadrature of f on [a, b] to the absolute tolerance tol.
template<class F>
inline scalar adaptiveSimpson
(
    const F& f, const scalar a, const scalar b, const scalar tol,
    const label maxDepth = 44
)
{
    const scalar fa = f(a), fb = f(b), fm = f(0.5*(a + b));
    const scalar whole = (b - a)/6.0*(fa + 4*fm + fb);
    return adaptiveSimpsonRec(f, a, b, fa, fm, fb, whole, tol, maxDepth);
}

//- Volume of the axis-aligned box [lo, hi] inside the ball (centre c, radius R),
//  exact up to the quadrature tolerance relTol x box volume.
inline scalar boxBallVolume
(
    const point& lo, const point& hi, const point& c, const scalar R,
    const scalar relTol = 1e-12
)
{
    const scalar z0 = max(lo.z(), c.z() - R), z1 = min(hi.z(), c.z() + R);
    if (z0 >= z1) return 0;
    const scalar x0 = lo.x() - c.x(), x1 = hi.x() - c.x();
    const scalar y0 = lo.y() - c.y(), y1 = hi.y() - c.y();
    auto A = [&](const scalar z) -> scalar
    {
        const scalar rr = R*R - (z - c.z())*(z - c.z());
        if (rr <= 0) return 0;
        return rectDiskArea(x0, x1, y0, y1, Foam::sqrt(rr));
    };
    const scalar tol = relTol*(hi.x() - lo.x())*(hi.y() - lo.y())*(hi.z() - lo.z());
    return adaptiveSimpson(A, z0, z1, tol);
}

} // End namespace davof
} // End namespace Foam

#endif

// ************************************************************************* //
