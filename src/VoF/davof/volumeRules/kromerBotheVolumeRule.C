/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           |
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

\*---------------------------------------------------------------------------*/

#include "kromerBotheVolumeRule.H"
#include "davofKromerBothe.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace davofVolumeRules
{
    defineTypeNameAndDebug(kromerBothe, 0);
    addToRunTimeSelectionTable(davofVolumeRule, kromerBothe, dictionary);
}
}


namespace
{
    // the self-test (defined below)
    void kbSelfTest(const Foam::scalar tolRel, const Foam::label maxDepth);
    // its vertex-grazing series and the anchor-event cell (defined below)
    void kbVertexGrazingTest
    (
        const Foam::scalar tolRel, const Foam::label maxDepth,
        const Foam::scalar bL, const char* pass, const bool withSplit
    );
    void kbEventCellTest(const char* pass);
    void kbEventCellScan
    (
        const char* pass, const char* name,
        const Foam::label i, const Foam::label j, const Foam::label k,
        const Foam::point& o, const Foam::vector& e3,
        const Foam::scalar b11, const Foam::scalar b12, const Foam::scalar b22,
        const Foam::scalar cAnchor
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::davofVolumeRules::kromerBothe::kromerBothe
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    davofVolumeRule(mesh, dict),
    debug_(dict.getOrDefault<bool>("volumeRuleDebug", false)),
    nDebugPrinted_(0),
    dbgMaxPrism_(0), dbgSumPrism_(0), dbgMaxCell_(0), dbgSumCell_(0),
    dbgNPrism_(0), dbgNCell_(0),
    cacheCell_(-1),
    cacheTimeIndex_(-1),
    cacheFaces_(24),
    prismFaces_(12)
{
    // the debug knobs of davofKromerBothe.H
    davof::kbMinSplit() = dict.getOrDefault<label>("kromerBotheMinSplit", 1);
    // kromerBotheFan: auto (the default: a face is fanned only where it is
    // warped) | true (every face fanned) | false (every face whole)
    davof::kbFanMode() = -1;
    if (const entry* ePtr = dict.findEntry("kromerBotheFan", keyType::LITERAL))
    {
        ITstream& is = ePtr->stream();
        const token tok(is);
        word fan("auto");
        if (tok.isWord()) fan = tok.wordToken();
        else if (tok.isNumber()) fan = (tok.number() != 0) ? "true" : "false";
        if (fan == "auto") davof::kbFanMode() = -1;
        else if (fan == "true" || fan == "on" || fan == "yes") davof::kbFanMode() = 1;
        else if (fan == "false" || fan == "off" || fan == "no") davof::kbFanMode() = 0;
        else
        {
            FatalIOErrorInFunction(dict)
                << "kromerBotheFan '" << fan
                << "' unknown. Valid: auto, true, false." << exit(FatalIOError);
        }
    }
    Info<< "kromerBothe: fan faces "
        << (davof::kbFanMode() < 0 ? "auto (warped faces only)"
          : (davof::kbFanMode() > 0 ? "true (every face)" : "false (whole faces)"))
        << endl;
    if (debug_)
    {
        Info<< "kromerBothe: minSplit " << davof::kbMinSplit()
            << ", fan mode " << davof::kbFanMode() << endl;
        kbSelfTest(tolRel_, maxDepth_);
    }
}


namespace
{

// The self-test of the rule (volumeRuleDebug): the liquid volume of an
// axis-aligned box under the axis-aligned paraboloid z = c + b (x^2 + y^2)/2
// (e3 = z, the liquid below), against a reference that is exact in x (the
// clamped quadratic integrated piecewise between its exact roots) and
// adaptive Gauss-Legendre in y, and against the split quadrature. Thin boxes
// in x reproduce the equatorial prisms of the sphere, whose traces on the
// side faces are parabolas and on the top and bottom circles.
void kbSelfTest(const Foam::scalar tolRel, const Foam::label maxDepth)
{
    using namespace Foam;
    const scalar h = 1.0/32, dx = 0.2*h, b = -7.0, c = 0.0;
    davof::cellSurface S;
    S.paraboloid = true;
    S.o = point::zero;
    S.e1 = vector(1, 0, 0); S.e2 = vector(0, 1, 0); S.e3 = vector(0, 0, 1);
    S.b11 = b; S.b22 = b; S.b12 = 0; S.c = c; S.a1 = 0; S.a2 = 0;
    auto phi = [&](const point& x) { return S.value(x); };
    const List<Pair<scalar>>& g = davof::gl12();

    // the reference
    auto inner = [&](const scalar y, const scalar x0, const scalar x1, const scalar z0) -> scalar
    {
        // clamp(gq(x), 0, h) with gq(x) = (c - z0 + b y^2/2) + b x^2/2
        const scalar k = c - z0 + 0.5*b*y*y;
        DynamicList<scalar> br(6);
        br.append(x0); br.append(x1);
        for (const scalar lvl : {scalar(0), h})
        {
            // b x^2/2 = lvl - k
            const scalar r2 = 2*(lvl - k)/b;
            if (r2 > 0)
            {
                const scalar r = std::sqrt(r2);
                if (r > x0 && r < x1) br.append(r);
                if (-r > x0 && -r < x1) br.append(-r);
            }
        }
        for (label i = 1; i < br.size(); ++i)
            for (label j = i; j > 0 && br[j] < br[j - 1]; --j) Swap(br[j], br[j - 1]);
        scalar sum = 0;
        for (label i = 0; i + 1 < br.size(); ++i)
        {
            const scalar xa = br[i], xb = br[i + 1];
            if (!(xb > xa)) continue;
            const scalar xm = 0.5*(xa + xb);
            const scalar gm = k + 0.5*b*xm*xm;
            if (gm >= h) sum += h*(xb - xa);
            else if (gm > 0) sum += k*(xb - xa) + b*(xb*xb*xb - xa*xa*xa)/6;
        }
        return sum;
    };
    auto reference = [&](const scalar x0, const scalar x1, const scalar y0, const scalar y1, const scalar z0) -> scalar
    {
        std::function<scalar(scalar, scalar, scalar, label)> rec =
            [&](scalar ya, scalar yb, scalar parent, label depth) -> scalar
            {
                auto gauss = [&](scalar a, scalar bq) {
                    scalar s = 0;
                    for (const Pair<scalar>& gw : g) s += gw.second()*(bq - a)*inner(a + gw.first()*(bq - a), x0, x1, z0);
                    return s;
                };
                const scalar ym = 0.5*(ya + yb);
                const scalar c0 = gauss(ya, ym), c1 = gauss(ym, yb);
                if ((depth >= 4 && mag(c0 + c1 - parent) <= 1e-15*h*h*h) || depth >= 40) return c0 + c1;
                return rec(ya, ym, c0, depth + 1) + rec(ym, yb, c1, depth + 1);
            };
        scalar p = 0;
        for (const Pair<scalar>& gw : g) p += gw.second()*(y1 - y0)*inner(y0 + gw.first()*(y1 - y0), x0, x1, z0);
        return rec(y0, y1, p, 0);
    };

    Info<< "kromerBothe self-test: thin boxes (dx = 0.2 h, h = 1/32) under z = "
        << c << " + " << b << " (x^2 + y^2)/2" << nl;
    label nBox = 0;
    for (const scalar x0 : {0.0, 0.02, 0.05, 0.1})
    {
        for (const scalar y0 : {0.0, 0.03, 0.1})
        {
            for (const scalar zShift : {0.5, 0.2, 0.8})
            {
                const scalar x1 = x0 + dx, y1 = y0 + h;
                // the surface height at the box centre, the box placed so that
                // the surface sits at zShift of its height there
                const scalar zc = c + 0.5*b*(sqr(0.5*(x0 + x1)) + sqr(0.5*(y0 + y1)));
                const scalar z0 = zc - zShift*h, z1 = z0 + h;
                const scalar pv = dx*h*h;
                const scalar vRef = reference(x0, x1, y0, y1, z0);
                // the box as a Kromer-Bothe polyhedron (fan faces about an
                // interior point)
                const point xin(0.5*(x0 + x1), 0.5*(y0 + y1), 0.5*(z0 + z1));
                DynamicList<davof::kbFace> faces(12);
                auto quad = [&](point a, point bq, point cq, point d) {
                    List<point> poly(4); poly[0] = a; poly[1] = bq; poly[2] = cq; poly[3] = d;
                    davof::kbFanFaces(poly, xin, faces);
                };
                quad(point(x0,y0,z0), point(x1,y0,z0), point(x1,y1,z0), point(x0,y1,z0));
                quad(point(x0,y0,z1), point(x1,y0,z1), point(x1,y1,z1), point(x0,y1,z1));
                quad(point(x0,y0,z0), point(x1,y0,z0), point(x1,y0,z1), point(x0,y0,z1));
                quad(point(x0,y1,z0), point(x1,y1,z0), point(x1,y1,z1), point(x0,y1,z1));
                quad(point(x0,y0,z0), point(x0,y1,z0), point(x0,y1,z1), point(x0,y0,z1));
                quad(point(x1,y0,z0), point(x1,y1,z0), point(x1,y1,z1), point(x1,y0,z1));
                scalar vKb = 0;
                label nLeavesKb = 0;
                const bool ok = davof::kbPolyhedronLiquidVolume(faces, S, tolRel*pv, vKb, nLeavesKb);
                // the split quadrature: the prism of the face x = x1 displaced by -d
                const vector d(dx, 0, 0);
                const point fc(x1, 0.5*(y0 + y1), 0.5*(z0 + z1));
                const point q0(x1,y0,z0), q1(x1,y1,z0), q2(x1,y1,z1), q3(x1,y0,z1);
                label nls = 0;
                scalar vSplit = 0;
                vSplit += davof::prismLiquidVolume(fc, q0, q1, d, phi, tolRel*pv, maxDepth, nls);
                vSplit += davof::prismLiquidVolume(fc, q1, q2, d, phi, tolRel*pv, maxDepth, nls);
                vSplit += davof::prismLiquidVolume(fc, q2, q3, d, phi, tolRel*pv, maxDepth, nls);
                vSplit += davof::prismLiquidVolume(fc, q3, q0, d, phi, tolRel*pv, maxDepth, nls);
                ++nBox;
                Info<< "  box " << nBox << " x0 " << x0 << " y0 " << y0 << " zShift " << zShift
                    << ": ref " << vRef/pv << " kb " << vKb/pv << (ok ? "" : " (NOT RESOLVED)")
                    << " split " << vSplit/pv
                    << "  (kb - ref)/pv " << (vKb - vRef)/pv
                    << "  (split - ref)/pv " << (vSplit - vRef)/pv << nl;
                if (mag(vKb - vRef) > 1e-9*pv)
                {
                    forAll(faces, k)
                    {
                        davof::kbSums sums;
                        davof::kbIntegrateFace(faces[k], S, 1e-16*pv, 1e-16*pv, true, sums);
                        Info<< "    face " << k << " n " << faces[k].n << " wet " << sums.lastArea
                            << " of " << faces[k].area() << " faceTerm " << sums.faceSum
                            << " line " << sums.lineSum << " arcs " << sums.nArcs << nl;
                    }
                }
            }
        }
    }

    // Exact per-face references for the axis-aligned frame (e3 = z, e1 = x,
    // e2 = y) with the general tensor: the liquid is z < hz(x, y) = c + (b11
    // x^2 + 2 b12 x y + b22 y^2)/2, a quadratic in x for fixed y and in y
    // for fixed x, so the clamped heights integrate exactly between their
    // roots along one direction and adaptively along the other
    struct AxisRef
    {
        scalar c, b11, b12, b22;
        scalar hz(scalar x, scalar y) const { return c + 0.5*(b11*x*x + 2*b12*x*y + b22*y*y); }
        // int_{xa}^{xb} clamp(hz(x, y) - zlo, 0, H) dx, exact
        scalar clampX(scalar y, scalar xa, scalar xb, scalar zlo, scalar H) const
        {
            // hz - zlo = A x^2 + B x + C
            const scalar A = 0.5*b11, B = b12*y, C = c + 0.5*b22*y*y - zlo;
            DynamicList<scalar> br(8);
            br.append(xa); br.append(xb);
            for (const scalar lvl : {scalar(0), H})
            {
                if (mag(A) > VSMALL)
                {
                    const scalar disc = B*B - 4*A*(C - lvl);
                    if (disc > 0)
                    {
                        const scalar sq = std::sqrt(disc);
                        for (const scalar r : {(-B - sq)/(2*A), (-B + sq)/(2*A)})
                            if (r > xa && r < xb) br.append(r);
                    }
                }
                else if (mag(B) > VSMALL)
                {
                    const scalar r = (lvl - C)/B;
                    if (r > xa && r < xb) br.append(r);
                }
            }
            for (label i = 1; i < br.size(); ++i)
                for (label j = i; j > 0 && br[j] < br[j - 1]; --j) Swap(br[j], br[j - 1]);
            scalar sum = 0;
            for (label i = 0; i + 1 < br.size(); ++i)
            {
                const scalar p = br[i], q = br[i + 1];
                if (!(q > p)) continue;
                const scalar m = 0.5*(p + q);
                const scalar gm = A*m*m + B*m + C;
                if (gm >= H) sum += H*(q - p);
                else if (gm > 0) sum += A*(q*q*q - p*p*p)/3 + B*(q*q - p*p)/2 + C*(q - p);
            }
            return sum;
        }
        // int_{xa}^{xb} 1[zlo < hz < zhi] w(x) dx with w(x) = c - (b11 x^2 + 2 b12 x y + b22 y^2)/2 (the patch integrand), exact
        scalar bandX(scalar y, scalar xa, scalar xb, scalar zlo, scalar zhi) const
        {
            const scalar A = 0.5*b11, B = b12*y, C = c + 0.5*b22*y*y;
            DynamicList<scalar> br(8);
            br.append(xa); br.append(xb);
            for (const scalar lvl : {zlo, zhi})
            {
                if (mag(A) > VSMALL)
                {
                    const scalar disc = B*B - 4*A*(C - lvl);
                    if (disc > 0)
                    {
                        const scalar sq = std::sqrt(disc);
                        for (const scalar r : {(-B - sq)/(2*A), (-B + sq)/(2*A)})
                            if (r > xa && r < xb) br.append(r);
                    }
                }
                else if (mag(B) > VSMALL)
                {
                    const scalar r = (lvl - C)/B;
                    if (r > xa && r < xb) br.append(r);
                }
            }
            for (label i = 1; i < br.size(); ++i)
                for (label j = i; j > 0 && br[j] < br[j - 1]; --j) Swap(br[j], br[j - 1]);
            scalar sum = 0;
            for (label i = 0; i + 1 < br.size(); ++i)
            {
                const scalar p = br[i], q = br[i + 1];
                if (!(q > p)) continue;
                const scalar m = 0.5*(p + q);
                const scalar hm = A*m*m + B*m + C;
                if (hm > zlo && hm < zhi)
                {
                    // w = c - A x^2 - B x - (C - c) ... w = c - 0.5 b11 x^2 - b12 x y - 0.5 b22 y^2
                    sum += c*(q - p) - A*(q*q*q - p*p*p)/3 - B*(q*q - p*p)/2 - 0.5*b22*y*y*(q - p);
                }
            }
            return sum;
        }
        // the adaptive outer integral over y of a function f(y)
        scalar outer(const std::function<scalar(scalar)>& f, scalar ya, scalar yb, scalar tol) const
        {
            const List<Pair<scalar>>& g = davof::gl12();
            auto gauss = [&](scalar a, scalar bq) {
                scalar s = 0;
                for (const Pair<scalar>& gw : g) s += gw.second()*(bq - a)*f(a + gw.first()*(bq - a));
                return s;
            };
            std::function<scalar(scalar, scalar, scalar, label)> rec =
                [&](scalar a, scalar bq, scalar parent, label depth) -> scalar
                {
                    const scalar m = 0.5*(a + bq);
                    const scalar c0 = gauss(a, m), c1 = gauss(m, bq);
                    if ((depth >= 4 && mag(c0 + c1 - parent) <= tol) || depth >= 40) return c0 + c1;
                    return rec(a, m, c0, depth + 1) + rec(m, bq, c1, depth + 1);
                };
            return rec(ya, yb, gauss(ya, yb), 0);
        }
    };

    // The second series: the frame tilted as at the sphere's equatorial cells
    // (e3 with a small component along the thin direction x) and a
    // non-diagonal tensor; the split quadrature at a tightened tolerance and
    // depth (1e-15 of the prism, depth 20) is the comparison next to the
    // split quadrature at the run's settings, and for the axis-aligned frame
    // (a = b = 0) the exact per-face references above
    Info<< "kromerBothe self-test, tilted frame: e3 = (sin a, -cos a sin b, cos a cos b), b12 = 0.5" << nl;
    for (const scalar aTilt : {0.0, 0.1, 0.3})
    {
        for (const scalar bTilt : {0.0, 0.33})
        {
            davof::cellSurface T;
            T.paraboloid = true;
            T.o = point::zero;
            T.e3 = vector(std::sin(aTilt), -std::cos(aTilt)*std::sin(bTilt), std::cos(aTilt)*std::cos(bTilt));
            davof::frameOf(T.e3, T.e1, T.e2);
            T.b11 = b; T.b22 = b; T.b12 = 0.5; T.c = 0; T.a1 = 0; T.a2 = 0;
            auto phiT = [&](const point& x) { return T.value(x); };
            for (const scalar t1 : {0.0, 0.02})
            {
                for (const scalar zShift : {0.5, 0.2, 0.8})
                {
                    const scalar t2 = 0.03;
                    // a point of the surface and the box about it
                    const scalar hs = T.c + 0.5*(T.b11*t1*t1 + 2*T.b12*t1*t2 + T.b22*t2*t2);
                    const point ps = T.o + t1*T.e1 + t2*T.e2 + hs*T.e3;
                    const point xin = ps + (0.5 - zShift)*h*vector(0, 0, 1);
                    const scalar x0 = xin.x() - 0.5*dx, x1 = xin.x() + 0.5*dx;
                    const scalar y0 = xin.y() - 0.5*h, y1 = xin.y() + 0.5*h;
                    const scalar z0 = xin.z() - 0.5*h, z1 = xin.z() + 0.5*h;
                    const scalar pv = dx*h*h;
                    DynamicList<davof::kbFace> faces(12);
                    // the box face each kbFace belongs to (one whole face or
                    // four fan triangles per quad, by kbFanMode)
                    DynamicList<label> quadOf(12);
                    label nQuad = 0;
                    auto quad = [&](point a, point bq, point cq, point d) {
                        List<point> poly(4); poly[0] = a; poly[1] = bq; poly[2] = cq; poly[3] = d;
                        davof::kbFanFaces(poly, xin, faces);
                        while (quadOf.size() < faces.size()) quadOf.append(nQuad);
                        ++nQuad;
                    };
                    quad(point(x0,y0,z0), point(x1,y0,z0), point(x1,y1,z0), point(x0,y1,z0));
                    quad(point(x0,y0,z1), point(x1,y0,z1), point(x1,y1,z1), point(x0,y1,z1));
                    quad(point(x0,y0,z0), point(x1,y0,z0), point(x1,y0,z1), point(x0,y0,z1));
                    quad(point(x0,y1,z0), point(x1,y1,z0), point(x1,y1,z1), point(x0,y1,z1));
                    quad(point(x0,y0,z0), point(x0,y1,z0), point(x0,y1,z1), point(x0,y0,z1));
                    quad(point(x1,y0,z0), point(x1,y1,z0), point(x1,y1,z1), point(x1,y0,z1));
                    scalar vKb = 0;
                    label nLeavesKb = 0;
                    const bool ok = davof::kbPolyhedronLiquidVolume(faces, T, tolRel*pv, vKb, nLeavesKb);
                    const vector d(dx, 0, 0);
                    const point fc(x1, 0.5*(y0 + y1), 0.5*(z0 + z1));
                    const point q0(x1,y0,z0), q1(x1,y1,z0), q2(x1,y1,z1), q3(x1,y0,z1);
                    auto splitAt = [&](const scalar tolA, const label depth) -> scalar {
                        label nls = 0;
                        scalar v = 0;
                        v += davof::prismLiquidVolume(fc, q0, q1, d, phiT, tolA, depth, nls);
                        v += davof::prismLiquidVolume(fc, q1, q2, d, phiT, tolA, depth, nls);
                        v += davof::prismLiquidVolume(fc, q2, q3, d, phiT, tolA, depth, nls);
                        v += davof::prismLiquidVolume(fc, q3, q0, d, phiT, tolA, depth, nls);
                        return v;
                    };
                    const scalar vSplit = splitAt(tolRel*pv, maxDepth);
                    const scalar vTight = splitAt(1e-15*pv, 20);
                    Info<< "  tilt a " << aTilt << " b " << bTilt << " t1 " << t1 << " zShift " << zShift
                        << ": kb " << vKb/pv << (ok ? "" : " (NOT RESOLVED)") << " split " << vSplit/pv
                        << " splitTight " << vTight/pv
                        << "  (kb - tight)/pv " << (vKb - vTight)/pv
                        << "  (split - tight)/pv " << (vSplit - vTight)/pv << nl;
                    if (mag(vKb - vTight) > 1e-9*pv)
                    {
                        // the Kromer-Bothe terms per box face (the fan
                        // triangles summed): wetted area, face term, line
                        scalar wetF[6] = {0, 0, 0, 0, 0, 0}, lineF[6] = {0, 0, 0, 0, 0, 0}, termF[6] = {0, 0, 0, 0, 0, 0};
                        forAll(faces, k)
                        {
                            davof::kbSums sums;
                            davof::kbIntegrateFace(faces[k], T, 1e-16*pv, 1e-16*pv, true, sums);
                            const label bf = quadOf[k];   // the box face, in the order quad() was called
                            wetF[bf] += sums.lastArea; lineF[bf] += sums.lineSum; termF[bf] += sums.faceSum;
                        }
                        const char* names[6] = {"z0 (n -z)", "z1 (n +z)", "y0 (n -y)", "y1 (n +y)", "x0 (n -x)", "x1 (n +x)"};
                        scalar lineTot = 0, termTot = 0;
                        for (label bf = 0; bf < 6; ++bf) { lineTot += lineF[bf]; termTot += termF[bf]; }
                        if (aTilt == 0 && bTilt == 0)
                        {
                            // exact references: e3 = z, e1 = x, e2 = y
                            AxisRef R{T.c, T.b11, T.b12, T.b22};
                            const scalar tolA = 1e-16*h*h;
                            // z-faces: the area where hz > z_c, as the clamped height with a unit band
                            auto areaZ = [&](scalar zc) {
                                return R.outer([&](scalar y) {
                                    // length of {x : hz(x, y) > zc} = clamp with H -> the indicator: use a thin band of height 1 scaled
                                    const scalar eps = 1e-9*h;
                                    return R.clampX(y, x0, x1, zc, eps)/eps;
                                }, y0, y1, tolA);
                            };
                            const scalar refZ0 = areaZ(z0), refZ1 = areaZ(z1);
                            const scalar refY0 = R.clampX(y0, x0, x1, z0, h), refY1 = R.clampX(y1, x0, x1, z0, h);
                            // x-faces: clamp along y at fixed x (hz quadratic in y): reuse clampX with the roles swapped
                            AxisRef Rs{T.c, T.b22, T.b12, T.b11};
                            const scalar refX0 = Rs.clampX(x0, y0, y1, z0, h), refX1 = Rs.clampX(x1, y0, y1, z0, h);
                            const scalar refLine = R.outer([&](scalar y) { return R.bandX(y, x0, x1, z0, z1); }, y0, y1, 1e-16*pv);
                            const scalar refWet[6] = {refZ0, refZ1, refY0, refY1, refX0, refX1};
                            // the face terms from the exact areas: (x_k - o).n_k A_k
                            const scalar dist[6] = {-z0, z1, -y0, y1, -x0, x1};
                            scalar refTerm = 0;
                            for (label bf = 0; bf < 6; ++bf)
                            {
                                refTerm += dist[bf]*refWet[bf];
                                Info<< "    " << names[bf] << ": wet kb " << wetF[bf] << " exact " << refWet[bf]
                                    << " diff/area " << (wetF[bf] - refWet[bf])/(bf < 2 ? dx*h : (bf < 4 ? dx*h : h*h))
                                    << "  line kb " << lineF[bf] << nl;
                            }
                            Info<< "    sums: faceTerm kb " << termTot << " exact " << refTerm
                                << "  line kb " << lineTot << " exact (int_S g) " << refLine
                                << "  V kb " << (termTot + lineTot)/3 << " exact " << (refTerm + refLine)/3
                                << " tight " << vTight << nl;
                            // the split quadrature's integrand against the exact liquid
                            // length of the segment [p - d, p] on a grid of the face
                            scalar maxDiff = 0; point worst = point::zero;
                            for (label i = 0; i < 60; ++i)
                            {
                                for (label j = 0; j < 60; ++j)
                                {
                                    const scalar y = y0 + (i + 0.5)/60*(y1 - y0), z = z0 + (j + 0.5)/60*(z1 - z0);
                                    const point p(x1, y, z);
                                    const scalar fSplit = davof::segmentLiquidFraction(p - d, p, phiT)*dx;
                                    const scalar fExact = R.clampX(y, x0, x1, z, 1e30)/1e30*0 + R.clampX(y, x0, x1, z, h*1e-9)/(h*1e-9);
                                    if (mag(fSplit - fExact) > maxDiff) { maxDiff = mag(fSplit - fExact); worst = p; }
                                }
                            }
                            Info<< "    split integrand vs exact segment length on a 60x60 grid: max diff/dx "
                                << maxDiff/dx << " at " << worst << nl;
                            // the sweep direction of the split quadrature's triangles reversed
                            label nlr = 0;
                            scalar vRev = 0;
                            vRev += davof::prismLiquidVolume(q0, q1, fc, d, phiT, 1e-15*pv, 20, nlr);
                            vRev += davof::prismLiquidVolume(q1, q2, fc, d, phiT, 1e-15*pv, 20, nlr);
                            vRev += davof::prismLiquidVolume(q2, q3, fc, d, phiT, 1e-15*pv, 20, nlr);
                            vRev += davof::prismLiquidVolume(q3, q0, fc, d, phiT, 1e-15*pv, 20, nlr);
                            // the split quadrature on the two triangles of the face's diagonal
                            label nld = 0;
                            scalar vDiag = 0;
                            vDiag += davof::prismLiquidVolume(q0, q1, q2, d, phiT, 1e-15*pv, 20, nld);
                            vDiag += davof::prismLiquidVolume(q0, q2, q3, d, phiT, 1e-15*pv, 20, nld);
                            Info<< "    split tight: fan " << vTight/pv << " (sweep reversed " << vRev/pv
                                << ", leaves " << nlr << "), diagonal triangles " << vDiag/pv
                                << " (leaves " << nld << "), exact " << (refTerm + refLine)/3/pv << nl;
                        }
                        else
                        {
                            for (label bf = 0; bf < 6; ++bf)
                                Info<< "    " << names[bf] << ": wet kb " << wetF[bf] << " line kb " << lineF[bf] << nl;
                        }
                    }
                }
            }
        }
    }
    Info<< endl;

    // the vertex-grazing series, before (the treatment of weakly curved
    // traces and the split rule's sweep apex off) and after (both on, the
    // default)
    for (const bool on : {false, true})
    {
        davof::kbGrazingFixOn() = on;
        davof::splitApexRotationOn() = on;
        const char* pass = on ? "after" : "before";
        kbVertexGrazingTest(tolRel, maxDepth, 0, pass, true);
        for (const scalar bL : {1e-12, -1e-10, 1e-8})
        {
            kbVertexGrazingTest(tolRel, maxDepth, bL, pass, false);
        }
        kbEventCellTest(pass);
    }
    davof::kbGrazingFixOn() = true;
    davof::splitApexRotationOn() = true;
}


// The exact liquid area of a planar convex polygon under the plane phi (phi
// linear): the polygon clipped to phi < 0 (Sutherland-Hodgman, the edge
// crossings by linear interpolation, exact for a plane up to round-off), its
// area vector by the shoelace sum about the first vertex.
Foam::vector clippedAreaVector
(
    const Foam::List<Foam::point>& poly,
    const std::function<Foam::scalar(const Foam::point&)>& phi
)
{
    using namespace Foam;
    DynamicList<point> c(poly.size() + 2);
    forAll(poly, i)
    {
        const label j = (i + 1) % poly.size();
        const scalar fi = phi(poly[i]), fj = phi(poly[j]);
        if (fi < 0) c.append(poly[i]);
        if ((fi < 0) != (fj < 0))
        {
            c.append(poly[i] + (fi/(fi - fj))*(poly[j] - poly[i]));
        }
    }
    vector a = Zero;
    for (label i = 1; i + 1 < c.size(); ++i)
    {
        a += 0.5*((c[i] - c[0]) ^ (c[i + 1] - c[0]));
    }
    return a;
}


// The vertex-grazing series (2026-10-25): boxes under a PLANE (b = 0, both
// rules exact) through a point at the signed distance d (relative to the box
// size) from one vertex, on the liquid and on the gas side, for the oblique
// normal of the plane test and axis-near normals, with the grazed vertex the
// highest (a corner cut of size d) or the median vertex along the normal (a
// generic cut passing within d of it). The references are exact: the faces
// clipped to the liquid side (clippedAreaVector) and the divergence theorem
// about a point of the plane, where the cap contributes nothing; the split
// rule is NOT the reference here. Compared: the cell volume (the Kromer-
// Bothe polyhedron of fan faces, the split quadrature's pyramids from the
// centre as in cellLiquidVolumeQuadric), the prism of the face x = x1 swept
// back by an axis and by an oblique displacement (kbPrismFaces' faces, the
// split rule's fan prisms as in faceSweptLiquidVolume) and the face fractions
// (the Kromer-Bothe fan as in kbFaceLiquidFraction, the quadratic cut as in
// faceCutDisplaced), each as |error| relative to the volume or the area.
// The crossing points of a convex polyhedron's edges with the plane phi = 0
// ordered around their centroid in the frame (e1, e2): the cap polygon; and
// its second moments about o in (u, v) = ((x - o).e1, (x - o).e2).
void capMoments
(
    const Foam::List<Foam::List<Foam::point>>& polys,
    const std::function<Foam::scalar(const Foam::point&)>& phi,
    const Foam::point& o, const Foam::vector& e1, const Foam::vector& e2,
    Foam::scalar& Iuu, Foam::scalar& Iuv, Foam::scalar& Ivv
)
{
    using namespace Foam;
    DynamicList<Pair<scalar>> q(16);
    for (const List<point>& poly : polys)
    {
        forAll(poly, i)
        {
            const label j = (i + 1) % poly.size();
            const scalar fi = phi(poly[i]), fj = phi(poly[j]);
            if ((fi < 0) != (fj < 0))
            {
                const point x = poly[i] + (fi/(fi - fj))*(poly[j] - poly[i]);
                q.append(Pair<scalar>((x - o) & e1, (x - o) & e2));
            }
        }
    }
    Iuu = Iuv = Ivv = 0;
    if (q.size() < 3) return;
    scalar cu = 0, cv = 0;
    forAll(q, i) { cu += q[i].first(); cv += q[i].second(); }
    cu /= q.size(); cv /= q.size();
    // sort by the angle about the centroid
    List<scalar> ang(q.size());
    forAll(q, i) ang[i] = std::atan2(q[i].second() - cv, q[i].first() - cu);
    for (label i = 1; i < q.size(); ++i)
        for (label j = i; j > 0 && ang[j] < ang[j - 1]; --j) { Swap(ang[j], ang[j - 1]); Swap(q[j], q[j - 1]); }
    forAll(q, i)
    {
        const label j = (i + 1) % q.size();
        const scalar ui = q[i].first(), vi = q[i].second(), uj = q[j].first(), vj = q[j].second();
        const scalar cr = ui*vj - uj*vi;
        Iuu += cr*(ui*ui + ui*uj + uj*uj)/12;
        Ivv += cr*(vi*vi + vi*vj + vj*vj)/12;
        Iuv += cr*(ui*vj + 2*ui*vi + 2*uj*vj + uj*vi)/24;
    }
}


void kbVertexGrazingTest
(
    const Foam::scalar tolRel,
    const Foam::label maxDepth,
    const Foam::scalar bL,          // the curvature scale times the box size
    const char* pass,               // "before" or "after"
    const bool withSplit            // the split rule on the unit box, oblique normal
)
{
    using namespace Foam;

    struct Box { const char* name; vector L; };
    const Box boxes[2] = {{"unit", vector(1, 1, 1)}, {"thin", vector(0.2, 1, 1)}};
    struct Normal { const char* name; vector n; };
    const Normal normals[4] =
    {
        {"oblique (0.3 0.5 0.8124)", vector(0.3, 0.5, 0.8124)},
        {"diagonal (1 1 1)", vector(1, 1, 1)},
        {"axis-near (1e-3 2e-3 1)", vector(1e-3, 2e-3, 1)},
        {"axis-near (1 1e-7 3e-7)", vector(1, 1e-7, 3e-7)}
    };
    const scalar ds[5] = {1e-6, 1e-8, 1e-10, 1e-12, 1e-14};

    // the maxima per rule and quantity over the series, per d
    scalar mxKbC[5] = {0}, mxSpC[5] = {0}, mxKbP[5] = {0}, mxSpP[5] = {0}, mxKbF[5] = {0}, mxCutF[5] = {0};

    Info<< "kromerBothe self-test, vertex grazing [" << pass << ", b L = " << bL << "]:"
        << " boxes under a surface through a point at d (relative to the box size) from"
        << " the highest or the median vertex along n; |error| against the exact"
        << " clipped polyhedron (with b: plus the first-order cap moment), relative to V"
        << " (cell, prism) or to the face area (faces, the maximum over the six; with b"
        << " against the quadratic cut); cases with an error above 1e-13 listed" << nl
        << "  columns: box normal vertex d side | cell kb split | prism(axis) kb split"
        << " | prism(oblique) kb split | face kb cut" << nl;

    for (const Box& bx : boxes)
    {
        const scalar x0 = 0, y0 = 0, z0 = 0;
        const scalar x1 = bx.L.x(), y1 = bx.L.y(), z1 = bx.L.z();
        const scalar Lref = cmptMax(bx.L);
        const point corners[8] =
        {
            point(x0,y0,z0), point(x1,y0,z0), point(x1,y1,z0), point(x0,y1,z0),
            point(x0,y0,z1), point(x1,y0,z1), point(x1,y1,z1), point(x0,y1,z1)
        };
        // the six faces, outward (counterclockwise seen from outside)
        List<List<point>> boxFaces(6, List<point>(4));
        auto setQ = [&](label k, label a, label b, label c, label d)
        {
            boxFaces[k][0] = corners[a]; boxFaces[k][1] = corners[b];
            boxFaces[k][2] = corners[c]; boxFaces[k][3] = corners[d];
        };
        setQ(0, 0, 3, 2, 1);   // z0
        setQ(1, 4, 5, 6, 7);   // z1
        setQ(2, 0, 1, 5, 4);   // y0
        setQ(3, 3, 7, 6, 2);   // y1
        setQ(4, 0, 4, 7, 3);   // x0
        setQ(5, 1, 2, 6, 5);   // x1
        const point xc(0.5*(x0 + x1), 0.5*(y0 + y1), 0.5*(z0 + z1));
        const scalar Vbox = (x1 - x0)*(y1 - y0)*(z1 - z0);

        for (const Normal& nm : normals)
        {
            const vector n = nm.n/mag(nm.n);
            // the vertices ordered along n: the highest and the median
            label order[8] = {0, 1, 2, 3, 4, 5, 6, 7};
            for (label i = 1; i < 8; ++i)
                for (label j = i; j > 0 && (corners[order[j]] & n) < (corners[order[j - 1]] & n); --j)
                    Swap(order[j], order[j - 1]);
            // the grazed point: the highest and the median vertex, and the
            // centre of the face x = x1 (the apex of its fan triangles in the
            // cell and in the prism's top: the anchor event at step 102 of
            // the Co 0.2 plane case had its dropped sliver there)
            const point fcX1 = 0.25*(corners[1] + corners[2] + corners[5] + corners[6]);
            const point vSel[3] = {corners[order[7]], corners[order[4]], fcX1};
            const char* vName[3] = {"highest", "median", "face centre"};

            for (label iv = 0; iv < 3; ++iv)
            {
                const point& v = vSel[iv];
                for (label id = 0; id < 5; ++id)
                {
                    for (const scalar side : {-1.0, 1.0})
                    {
                        // phi(x) = n.(x - p0), p0 = v + side d L n: phi(v) =
                        // -side d L, side +1 the vertex on the liquid side
                        const scalar dd = ds[id]*Lref;
                        const point p0 = v + side*dd*n;
                        // the paraboloid path; b = 0 a plane, else a mixed
                        // tensor of the scale bL/L (the round-off tensor of a
                        // plane in the solver is 1e-10 to 1e-8 /m at h 2.5e-5,
                        // b h 1e-15 to 1e-13 and above)
                        davof::cellSurface S;
                        S.paraboloid = true;
                        S.o = p0;
                        S.e3 = n;
                        davof::frameOf(S.e3, S.e1, S.e2);
                        S.b11 = -0.4*bL/Lref; S.b12 = 0.2*bL/Lref; S.b22 = bL/Lref;
                        S.c = 0; S.a1 = 0; S.a2 = 0;
                        std::function<scalar(const point&)> phi = [&](const point& x) { return S.value(x); };
                        // the tangent plane at o
                        std::function<scalar(const point&)> phiPlane = [&](const point& x) { return (x - p0) & n; };

                        // the exact liquid volume of a convex polyhedron under
                        // the plane, plus the first-order volume of the
                        // paraboloid's height over the cap (error (b L)^2)
                        auto exactVolume = [&](const List<List<point>>& polys) -> scalar
                        {
                            scalar s = 0;
                            for (const List<point>& poly : polys)
                            {
                                const vector a = clippedAreaVector(poly, phiPlane);
                                s += (poly[0] - p0) & a;
                            }
                            s /= 3;
                            if (bL != 0)
                            {
                                scalar Iuu, Iuv, Ivv;
                                capMoments(polys, phiPlane, S.o, S.e1, S.e2, Iuu, Iuv, Ivv);
                                // the sign of the cap polygon's orientation
                                // (the moments are signed by it): with b11 = b22
                                // = 1 the volume added is positive
                                const scalar sgn = (Iuu + Ivv >= 0) ? 1 : -1;
                                s += sgn*0.5*(S.b11*Iuu + 2*S.b12*Iuv + S.b22*Ivv);
                            }
                            return s;
                        };
                        const bool doSplit = withSplit && iv >= 0 && &bx == &boxes[0] && &nm == &normals[0];

                        // the cell: Kromer-Bothe (fan faces) and split (pyramids)
                        const scalar vRef = exactVolume(boxFaces);
                        DynamicList<davof::kbFace> faces(24);
                        for (const List<point>& poly : boxFaces) davof::kbFanFaces(poly, xc, faces);
                        scalar vKb = 0;
                        label nlk = 0;
                        const bool okC = davof::kbPolyhedronLiquidVolume(faces, S, tolRel*Vbox, vKb, nlk);
                        scalar vSp = 0;
                        label nls = 0;
                        if (doSplit)
                        {
                            for (const List<point>& poly : boxFaces)
                            {
                                point fc = Zero;
                                for (const point& q : poly) fc += q;
                                fc /= poly.size();
                                forAll(poly, i)
                                {
                                    const label j = (i + 1) % poly.size();
                                    vSp += davof::pyramidLiquidVolume(xc, fc, poly[i], poly[j], phi, tolRel*Vbox, maxDepth, nls);
                                }
                            }
                        }
                        const scalar eKbC = okC ? mag(vKb - vRef)/Vbox : GREAT;
                        const scalar eSpC = doSplit ? mag(vSp - vRef)/Vbox : -1;

                        // the prisms of the face x = x1 swept back by dv: the
                        // faces of kbPrismFaces, the fan prisms of the split rule
                        auto prismErrors = [&](const vector& dv, scalar& eKb, scalar& eSp)
                        {
                            const List<point>& top = boxFaces[5];
                            List<point> bottom(top.size());
                            forAll(top, i) bottom[i] = top[i] - dv;
                            point fc = Zero;
                            for (const point& q : top) fc += q;
                            fc /= top.size();
                            const point xin = fc - 0.5*dv;
                            DynamicList<davof::kbFace> pf(12);
                            davof::kbFanFaces(top, xin, pf);
                            davof::kbFanFaces(bottom, xin, pf);
                            forAll(top, i)
                            {
                                const label j = (i + 1) % top.size();
                                davof::kbFace F;
                                F.p.setSize(4);
                                F.p[0] = top[i]; F.p[1] = top[j]; F.p[2] = bottom[j]; F.p[3] = bottom[i];
                                F.setFrame();
                                if (((F.x0 - xin) & F.n) < 0) { F.reverse(); F.setFrame(); }
                                pf.append(F);
                            }
                            // the exact reference from the same outward polygons
                            List<List<point>> polys(pf.size());
                            forAll(pf, k) polys[k] = pf[k].p;
                            const scalar Vp = mag(dv & (0.5*((top[2] - top[0]) ^ (top[3] - top[1]))));
                            const scalar vR = exactVolume(polys);
                            scalar vK = 0;
                            label nl1 = 0;
                            const bool okP = davof::kbPolyhedronLiquidVolume(pf, S, tolRel*Vp, vK, nl1);
                            scalar vS = 0;
                            label nl2 = 0;
                            if (doSplit)
                            {
                                forAll(top, i)
                                {
                                    const label j = (i + 1) % top.size();
                                    vS += davof::prismLiquidVolume(fc, top[i], top[j], dv, phi, tolRel*Vp, maxDepth, nl2);
                                }
                            }
                            eKb = okP ? mag(vK - vR)/Vp : GREAT;
                            eSp = doSplit ? mag(vS - vR)/Vp : -1;
                        };
                        scalar eKbPa, eSpPa, eKbPo, eSpPo;
                        prismErrors(vector(0.37*(x1 - x0), 0, 0), eKbPa, eSpPa);
                        prismErrors(vector(0.37*(x1 - x0), 0.21*(y1 - y0), -0.13*(z1 - z0)), eKbPo, eSpPo);

                        // the face fractions
                        scalar eKbF = 0, eCutF = 0;
                        for (const List<point>& poly : boxFaces)
                        {
                            point fc = Zero;
                            for (const point& q : poly) fc += q;
                            fc /= poly.size();
                            const scalar A = mag(0.5*((poly[2] - poly[0]) ^ (poly[3] - poly[1])));
                            // b = 0: the exact clipped area; else the quadratic
                            // cut (aRef replaced below)
                            scalar aRef = mag(clippedAreaVector(poly, phi));
                            scalar aKb = 0, aCut = 0;
                            bool okF = true;
                            forAll(poly, i)
                            {
                                const label j = (i + 1) % poly.size();
                                davof::kbFace F;
                                F.p.setSize(3);
                                F.p[0] = fc; F.p[1] = poly[i]; F.p[2] = poly[j];
                                F.setFrame();
                                davof::kbSums sums;
                                okF = davof::kbIntegrateFace(F, S, tolRel*A, GREAT, false, sums) && okF;
                                aKb += sums.lastArea;
                                label nfq = 0;
                                aCut += davof::triLiquidArea(fc, poly[i], poly[j], phi, 6, nfq);
                            }
                            if (bL != 0) aRef = aCut;
                            eKbF = max(eKbF, okF ? mag(aKb - aRef)/A : GREAT);
                            eCutF = max(eCutF, mag(aCut - aRef)/A);
                        }

                        mxKbC[id] = max(mxKbC[id], eKbC); mxSpC[id] = max(mxSpC[id], eSpC);
                        mxKbP[id] = max(mxKbP[id], max(eKbPa, eKbPo)); mxSpP[id] = max(mxSpP[id], max(eSpPa, eSpPo));
                        mxKbF[id] = max(mxKbF[id], eKbF); mxCutF[id] = max(mxCutF[id], eCutF);

                        const scalar eMax = max(max(max(eKbC, eSpC), max(max(eKbPa, eSpPa), max(eKbPo, eSpPo))), max(eKbF, eCutF));
                        if (eMax > 1e-13)
                        {
                            Info<< "  [" << pass << ", b L " << bL << "] " << bx.name << " | " << nm.name << " | " << vName[iv]
                                << " | d " << ds[id] << " " << (side > 0 ? "liquid" : "gas")
                                << " | alpha " << vRef/Vbox
                                << " | cell " << eKbC << " " << eSpC
                                << " | prism(axis) " << eKbPa << " " << eSpPa
                                << " | prism(oblique) " << eKbPo << " " << eSpPo
                                << " | face " << eKbF << " " << eCutF << nl;
                        }
                    }
                }
            }
        }
    }
    Info<< "kromerBothe self-test, vertex grazing [" << pass << ", b L = " << bL
        << "]: maxima over boxes, normals, vertices and sides (split: unit box,"
        << " oblique normal; -1 not computed)" << nl;
    for (label id = 0; id < 5; ++id)
    {
        Info<< "  [" << pass << ", b L " << bL << "] d " << ds[id] << ": cell kb " << mxKbC[id] << " split " << mxSpC[id]
            << " | prism kb " << mxKbP[id] << " split " << mxSpP[id]
            << " | face kb " << mxKbF[id] << " cut " << mxCutF[id] << nl;
    }
    Info<< endl;
}


// The cells of the anchor events of the translated plane (N = 160, ring row,
// oblique arm, h = 2.5e-5): Co 0.5 step 12, cell 2815517 = (157, 156, 109),
// a vertex 2e-4 h inside the liquid; Co 0.2 step 102, cell 3427807 =
// (127, 143, 133), a face centre 1.9e-4 h inside. The surface as the solver
// printed it, the offset scanned over 1e-15 m about the anchored value (the
// rule's volume toggled by 1.2e-5 of the cell between offsets 1e-16 m apart):
// the Kromer-Bothe volume against the exact plane volume plus the first-order
// cap moment (b h ~ 1e-12, second order 1e-24), the maximum and the number of
// offsets above 1e-12 of the cell.
void kbEventCellTest(const char* pass)
{
    struct EventCell
    {
        const char* name;
        Foam::label i, j, k;
        Foam::point o;
        Foam::vector e3;
        Foam::scalar b11, b12, b22, c;
    };
    const EventCell cells[2] =
    {
        {
            "Co 0.5 step 12 cell 2815517", 157, 156, 109,
            Foam::point(0.003937499998253168, 0.00391249999825285, 0.002737302906693282),
            Foam::vector(-0.3000009360043708, -0.5000015600073782, -0.8124025346998184),
            -1.182816158405161e-08, 5.056734871100458e-09, 5.603200686295972e-08,
            -8.473360899599996e-11
        },
        {
            "Co 0.2 step 102 cell 3427807", 127, 143, 133,
            Foam::point(0.003189997953676344, 0.003592995493897087, 0.0033456890041218),
            Foam::vector(-0.3000009360039511, -0.5000015600066313, -0.812402534700433),
            -6.987968795410825e-09, 4.442046833528177e-10, -1.159328084991024e-08,
            -4.531655948356755e-10
        }
    };
    for (const EventCell& ec : cells)
    {
        kbEventCellScan(pass, ec.name, ec.i, ec.j, ec.k, ec.o, ec.e3, ec.b11, ec.b12, ec.b22, ec.c);
    }
}


void kbEventCellScan
(
    const char* pass, const char* name,
    const Foam::label i, const Foam::label j, const Foam::label k,
    const Foam::point& o, const Foam::vector& e3,
    const Foam::scalar b11, const Foam::scalar b12, const Foam::scalar b22,
    const Foam::scalar cAnchor
)
{
    using namespace Foam;
    const scalar h = 2.5e-5;
    const scalar x0 = i*h, y0 = j*h, z0 = k*h, x1 = x0 + h, y1 = y0 + h, z1 = z0 + h;
    const point corners[8] =
    {
        point(x0,y0,z0), point(x1,y0,z0), point(x1,y1,z0), point(x0,y1,z0),
        point(x0,y0,z1), point(x1,y0,z1), point(x1,y1,z1), point(x0,y1,z1)
    };
    List<List<point>> boxFaces(6, List<point>(4));
    const label q[6][4] = {{0,3,2,1}, {4,5,6,7}, {0,1,5,4}, {3,7,6,2}, {0,4,7,3}, {1,2,6,5}};
    for (label f = 0; f < 6; ++f) for (label m = 0; m < 4; ++m) boxFaces[f][m] = corners[q[f][m]];
    const point xc(0.5*(x0 + x1), 0.5*(y0 + y1), 0.5*(z0 + z1));
    const scalar Vc = h*h*h;

    davof::cellSurface S;
    S.paraboloid = true;
    S.o = o;
    S.e3 = e3;
    davof::frameOf(S.e3, S.e1, S.e2);
    S.b11 = b11; S.b12 = b12; S.b22 = b22;
    S.a1 = 0; S.a2 = 0;

    DynamicList<davof::kbFace> faces(24);
    for (const List<point>& poly : boxFaces) davof::kbFanFaces(poly, xc, faces);
    scalar eMax = 0, cWorst = cAnchor;
    label nBad = 0, nFail = 0;
    const label nScan = 201;
    for (label s = 0; s < nScan; ++s)
    {
        S.c = cAnchor + (s - nScan/2)*1e-17;
        std::function<scalar(const point&)> phiPlane = [&](const point& x) { return ((x - S.o) & S.e3) - S.c; };
        scalar vRef = 0;
        for (const List<point>& poly : boxFaces)
        {
            const point pO = S.o + S.c*S.e3;   // a point of the plane
            vRef += (poly[0] - pO) & clippedAreaVector(poly, phiPlane);
        }
        vRef /= 3;
        scalar Iuu, Iuv, Ivv;
        capMoments(boxFaces, phiPlane, S.o, S.e1, S.e2, Iuu, Iuv, Ivv);
        const scalar sgn = (Iuu + Ivv >= 0) ? 1 : -1;
        vRef += sgn*0.5*(S.b11*Iuu + 2*S.b12*Iuv + S.b22*Ivv);
        scalar v = 0;
        label nl = 0;
        if (!davof::kbPolyhedronLiquidVolume(faces, S, 1e-11*Vc, v, nl)) { ++nFail; continue; }
        const scalar e = mag(v - vRef)/Vc;
        if (e > 1e-12) ++nBad;
        if (e > eMax) { eMax = e; cWorst = S.c; }
    }
    Info<< "kromerBothe self-test, the anchor-event cell " << name << " [" << pass << "]: " << nScan
        << " offsets within 1e-15 m of " << cAnchor << ": max |kb - exact|/V " << eMax
        << " at c " << cWorst << ", " << nBad << " offsets above 1e-12, " << nFail
        << " unresolved" << endl;
}


// The faces of a polyhedron and their wetted areas under S, for the debug print
void printFaces
(
    const Foam::DynamicList<Foam::davof::kbFace>& faces,
    const Foam::davof::cellSurface& S,
    const Foam::scalar tolAbs
)
{
    using namespace Foam;
    forAll(faces, k)
    {
        const davof::kbFace& F = faces[k];
        davof::kbSums sums;
        const bool ok = davof::kbIntegrateFace(F, S, tolAbs, tolAbs, true, sums);
        Pout<< "    face " << k << " n " << F.n << " x0 " << F.x0 << " Ls " << F.Ls
            << " area " << F.area() << " wet " << sums.lastArea
            << " faceTerm " << sums.faceSum << " line " << sums.lineSum
            << " arcs " << sums.nArcs << (ok ? "" : " NOT RESOLVED") << nl;
        forAll(F.p, i)
        {
            Pout<< "      p " << F.p[i] << " phi " << S.value(F.p[i]) << nl;
        }
    }
}

}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::davofVolumeRules::kromerBothe::cellLiquidVolume
(
    const label celli,
    const davof::cellSurface& S,
    const scalar tolRel
) const
{
    scalar dVdc = -1;
    return cellLiquidVolumeSlope(celli, S, tolRel, dVdc);
}


Foam::scalar Foam::davofVolumeRules::kromerBothe::cellLiquidVolumeSlope
(
    const label celli,
    const davof::cellSurface& S,
    const scalar tolRel,
    scalar& dVdc
) const
{
    dVdc = -1;
    if (S.sphere)
    {
        ++nFallback_;
        return splitCellLiquidVolume(celli, S, tolRel);
    }
    // the cell's faces with their frames depend on the mesh only: built once
    // per cell and time step, reused by the anchor's evaluations
    const label timeIndex = mesh_.time().timeIndex();
    if (celli != cacheCell_ || timeIndex != cacheTimeIndex_)
    {
        davof::kbCellFaces(mesh_, celli, cacheFaces_);
        cacheCell_ = celli;
        cacheTimeIndex_ = timeIndex;
    }
    const DynamicList<davof::kbFace>& faces = cacheFaces_;
    scalar v = 0, patchArea = 0;
    const scalar Vc = mesh_.V()[celli];
    const bool resolved = davof::kbPolyhedronLiquidVolume(faces, S, tolRel*Vc, v, nLeaves_, &patchArea);
    const bool inRange = v >= -1e-9*Vc && v <= (1 + 1e-9)*Vc;
    if (resolved && !inRange) ++davof::kbFailCount(davof::KB_VOLUME);
    if (resolved && inRange)
    {
        dVdc = patchArea;
        if (debug_)
        {
            const scalar vs = splitCellLiquidVolume(celli, S, tolRel);
            const scalar rel = mag(v - vs)/Vc;
            dbgMaxCell_ = max(dbgMaxCell_, rel);
            dbgSumCell_ += rel;
            ++dbgNCell_;
            // the cells whose volume the two rules disagree on beyond
            // round-off (on the translated plane at N = 160 a few cells per
            // step reach 2e-10 of the cell while every prism agrees to 6e-14,
            // 2026-10-25): the per-face breakdown names the face
            if (nDebugPrinted_ < 12 && rel > 1e-12)
            {
                ++nDebugPrinted_;
                Pout<< "kromerBothe cell mismatch: cell " << celli << " V " << Vc
                    << " kb " << v << " split " << vs << " diff/V " << (v - vs)/Vc
                    << " alpha_kb " << v/Vc << " surface o " << S.o << " e3 " << S.e3
                    << " c " << S.c << " b " << S.b11 << " " << S.b12 << " " << S.b22 << nl;
                printFaces(faces, S, tolRel*Vc);
            }
        }
        return min(max(v, scalar(0)), Vc);
    }
    // not resolved (a wild surface, a degenerate conic): the capped split
    // quadrature; clamped to [0, V] as the rule's own value is
    ++nFallback_;
    return min(max(fallbackCellLiquidVolume(celli, S), scalar(0)), Vc);
}


Foam::scalar Foam::davofVolumeRules::kromerBothe::prismLiquidVolume
(
    const label facei,
    const vector& d,
    const davof::cellSurface& S,
    const scalar tolRel
) const
{
    const vector& Sf = mesh_.faceAreas()[facei];
    const scalar prismVol = mag(d & Sf);
    // a face parallel to the displacement sweeps no volume (relative test:
    // the flat prism of a round-off product is not a closed polyhedron and
    // the face sum does not cancel)
    if (prismVol <= 1e-10*mag(d)*mag(Sf)) return 0;
    if (S.sphere)
    {
        ++nFallback_;
        return splitPrismLiquidVolume(facei, d, S, tolRel);
    }
    DynamicList<davof::kbFace>& faces = prismFaces_;
    davof::kbPrismFaces(mesh_, facei, d, faces);
    scalar v = 0;
    const bool resolved = davof::kbPolyhedronLiquidVolume(faces, S, tolRel*prismVol, v, nLeaves_);
    const bool inRange = v >= -1e-9*prismVol && v <= (1 + 1e-9)*prismVol;
    if (resolved && !inRange) ++davof::kbFailCount(davof::KB_VOLUME);
    if (resolved && inRange)
    {
        v = min(max(v, scalar(0)), prismVol);
        if (debug_)
        {
            const scalar vs = splitPrismLiquidVolume(facei, d, S, tolRel);
            const scalar rel = mag(v - vs)/prismVol;
            dbgMaxPrism_ = max(dbgMaxPrism_, rel);
            dbgSumPrism_ += rel;
            ++dbgNPrism_;
            if (nDebugPrinted_ < 12 && mag(v - vs) > 1e-8*prismVol)
            {
                ++nDebugPrinted_;
                Pout<< "kromerBothe prism mismatch: face " << facei << " d " << d
                    << " prismVol " << prismVol << " kb " << v << " split " << vs
                    << " diff/prismVol " << (v - vs)/prismVol << nl;
                // a third value: the brute-force midpoint grid over the face
                // (a parallelogram spanned from its first vertex) with the
                // exact liquid fraction of every segment [p - d, p], at
                // three grid sizes and Richardson-extrapolated (the
                // integrand has kinks on the traces, so the grid converges
                // at second order at best)
                const face& f = mesh_.faces()[facei];
                if (f.size() == 4)
                {
                    const pointField& P = mesh_.points();
                    const point& p0 = P[f[0]];
                    const vector ea = P[f[1]] - p0, eb = P[f[3]] - p0;
                    auto phi = [&](const point& x) { return S.value(x); };
                    scalar vb[3];
                    const label ns[3] = {200, 400, 800};
                    for (label k = 0; k < 3; ++k)
                    {
                        const label n = ns[k];
                        scalar sum = 0;
                        for (label i = 0; i < n; ++i)
                        {
                            for (label j = 0; j < n; ++j)
                            {
                                const point p = p0 + ((i + 0.5)/n)*ea + ((j + 0.5)/n)*eb;
                                sum += davof::segmentLiquidFraction(p - d, p, phi);
                            }
                        }
                        vb[k] = prismVol*sum/(scalar(n)*n);
                    }
                    const scalar rich = vb[2] + (vb[2] - vb[1])/3.0;
                    Pout<< "  brute force (midpoint grid, exact segments) n=200 " << vb[0]
                        << " n=400 " << vb[1] << " n=800 " << vb[2]
                        << " Richardson " << rich << " (kb - bf)/prismVol " << (v - rich)/prismVol
                        << " (split - bf)/prismVol " << (vs - rich)/prismVol << nl;
                }
                printFaces(faces, S, tolRel*prismVol);
            }
        }
        return v;
    }
    ++nFallback_;
    return min(max(fallbackPrismLiquidVolume(facei, d, S), scalar(0)), prismVol);
}


void Foam::davofVolumeRules::kromerBothe::report() const
{
    // the reasons of the step's fallbacks (every processor's own count)
    {
        label n[davof::KB_NFAIL];
        label tot = 0;
        for (label r = 0; r < davof::KB_NFAIL; ++r)
        {
            n[r] = davof::kbFailCount(r);
            davof::kbFailCount(r) = 0;
            if (Pstream::parRun()) reduce(n[r], sumOp<label>());
            tot += n[r];
        }
        if (tot > 0)
        {
            Info<< "kromerBothe fallbacks by reason: conic " << n[davof::KB_CONIC]
                << ", point " << n[davof::KB_AT] << ", root " << n[davof::KB_ROOT]
                << ", gradient " << n[davof::KB_GRAD] << ", area " << n[davof::KB_AREA]
                << ", volume " << n[davof::KB_VOLUME] << endl;
        }
    }
    if (!debug_) return;
    scalar mp = dbgMaxPrism_, sp = dbgSumPrism_, mc = dbgMaxCell_, sc = dbgSumCell_;
    label np = dbgNPrism_, nc = dbgNCell_;
    if (Pstream::parRun())
    {
        reduce(mp, maxOp<scalar>()); reduce(sp, sumOp<scalar>());
        reduce(mc, maxOp<scalar>()); reduce(sc, sumOp<scalar>());
        reduce(np, sumOp<label>()); reduce(nc, sumOp<label>());
    }
    Info<< "kromerBothe debug: |kb - split| over " << np << " prisms max "
        << mp << " mean " << sp/max(np, label(1)) << " of the prism volume; over "
        << nc << " cell volumes max " << mc << " mean " << sc/max(nc, label(1))
        << " of the cell volume" << endl;
    dbgMaxPrism_ = dbgSumPrism_ = dbgMaxCell_ = dbgSumCell_ = 0;
    dbgNPrism_ = dbgNCell_ = 0;
    nDebugPrinted_ = 0;
}


Foam::scalar Foam::davofVolumeRules::kromerBothe::faceLiquidFraction
(
    const label facei,
    const vector& d,
    const davof::cellSurface& S,
    label& nFallback
) const
{
    scalar af = 0;
    if (!S.sphere && davof::kbFaceLiquidFraction(mesh_, facei, d, S, tolRel_, af, nLeaves_))
    {
        if (debug_ && nDebugPrinted_ < 12)
        {
            label nfb = 0;
            const scalar aq = davofVolumeRule::faceLiquidFraction(facei, d, S, nfb);
            if (mag(af - aq) > 1e-8)
            {
                ++nDebugPrinted_;
                Pout<< "kromerBothe face mismatch: face " << facei << " d " << d
                    << " kb " << af << " quadratic cut " << aq << nl;
            }
        }
        return af;
    }
    ++nFallback_;
    return davofVolumeRule::faceLiquidFraction(facei, d, S, nFallback);
}


// ************************************************************************* //
