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

Application
    testCutFaceAdvectEndArea

Description
    Unit check of cutFaceAdvect::subFaceAreaEnd(), the liquid area of a face
    at the end of the sweep of a plane moving at constant normal speed Un0
    along n0 from x0 (the DAVOF face update). On the six unit square faces
    of a single cell (case/), for in-plane, diagonal and oblique sweep
    directions, both signs of Un0 and a range of dt, the stored end area is
    compared with the exact area of the face polygon clipped to the liquid
    side {(p - x0).n0 < Un0 dt} (Sutherland-Hodgman), and the time integral
    timeIntegratedFaceFlux (phi = magSf) with Simpson's rule on every
    sub-interval between the vertex arrival times (the exact area is
    quadratic in t there, so the rule is exact). Prints the maxima and exits
    non-zero above 1e-12 on either.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "cutFaceAdvect.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// The area of the convex polygon pts clipped to {p : (p - x0).n0 < s}
Foam::scalar clippedArea
(
    const Foam::List<Foam::point>& pts,
    const Foam::vector& n0,
    const Foam::point& x0,
    const Foam::scalar s
)
{
    using namespace Foam;
    DynamicList<point> out(8);
    const label n = pts.size();
    for (label i = 0; i < n; ++i)
    {
        const point& a = pts[i];
        const point& b = pts[(i + 1) % n];
        const scalar da = ((a - x0) & n0) - s;
        const scalar db = ((b - x0) & n0) - s;
        if (da <= 0)
        {
            out.append(a);
        }
        if ((da < 0 && db > 0) || (da > 0 && db < 0))
        {
            out.append(a + (da/(da - db))*(b - a));
        }
    }
    if (out.size() < 3)
    {
        return 0;
    }
    vector A(Zero);
    for (label i = 1; i < out.size() - 1; ++i)
    {
        A += 0.5*((out[i] - out[0]) ^ (out[i + 1] - out[0]));
    }
    return mag(A);
}


int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    volScalarField alpha1
    (
        IOobject("alpha1", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE),
        mesh,
        dimensionedScalar(dimless, Zero)
    );

    cutFaceAdvect cf(mesh, alpha1);

    const faceList& faces = mesh.faces();
    const pointField& points = mesh.points();

    scalar maxEndErr = 0, maxIntErr = 0;
    label nChecks = 0;

    // Sweep directions, expressed in the local frame of each face (e1, e2
    // in the face plane, nf the face normal): in-plane axis, in-plane
    // diagonal, oblique out of the plane
    const List<vector> dirs
    ({
        vector(1, 0, 0),
        vector(0, 1, 0),
        vector(1, 1, 0),
        vector(0.6, -0.8, 0),
        vector(0.3, 0.5, 0.8124),
        vector(-0.5, 0.2, 0.7)
    });
    const List<scalar> dts({0.05, 0.2, 0.37, 0.5, 0.75, 1.0, 1.3, 1.9, 2.5});
    const List<scalar> offsets({-0.45, -0.17, 0.0, 0.21, 0.6});

    forAll(faces, facei)
    {
        const face& f = faces[facei];
        List<point> pts(f.points(points));
        const point Cf = mesh.faceCentres()[facei];
        const vector Sf = mesh.faceAreas()[facei];
        const scalar magSf = mag(Sf);
        const vector nf = Sf/magSf;
        vector e1 = pts[1] - pts[0];
        e1 /= mag(e1);
        const vector e2 = nf ^ e1;

        forAll(dirs, di)
        {
            vector n0 = dirs[di].x()*e1 + dirs[di].y()*e2 + dirs[di].z()*nf;
            n0 /= mag(n0);
            forAll(offsets, oi)
            {
                // x0 on the sweep axis through the face centre
                const point x0 = Cf + offsets[oi]*n0 + 0.013*e1 - 0.021*e2;
                for (const scalar sgn : {scalar(1), scalar(-1)})
                {
                    const scalar Un0 = sgn*1.0;
                    forAll(dts, ti)
                    {
                        const scalar dt = dts[ti];
                        const scalar dVf =
                            cf.timeIntegratedFaceFlux
                            (
                                facei, x0, n0, Un0, dt, magSf, magSf
                            );
                        const scalar areaEnd = cf.subFaceAreaEnd();
                        const scalar exactEnd = clippedArea(pts, n0, x0, Un0*dt);

                        // The exact area is quadratic in t between the
                        // vertex arrival times: Simpson on every sub-interval
                        // of [0, dt] between them is exact
                        DynamicList<scalar> knots(8);
                        knots.append(0);
                        forAll(pts, pi)
                        {
                            const scalar tv = ((pts[pi] - x0) & n0)/Un0;
                            if (tv > 0 && tv < dt) knots.append(tv);
                        }
                        knots.append(dt);
                        Foam::sort(knots);
                        scalar integral = 0;
                        for (label k = 0; k + 1 < knots.size(); ++k)
                        {
                            const scalar a = knots[k], b = knots[k + 1];
                            integral +=
                                (b - a)/6.0
                               *(
                                    clippedArea(pts, n0, x0, Un0*a)
                                  + 4*clippedArea(pts, n0, x0, Un0*0.5*(a + b))
                                  + clippedArea(pts, n0, x0, Un0*b)
                                );
                        }

                        const scalar endErr = mag(areaEnd - exactEnd);
                        const scalar intErr = mag(dVf - integral);
                        maxEndErr = max(maxEndErr, endErr);
                        maxIntErr = max(maxIntErr, intErr);
                        ++nChecks;
                        if (endErr > 1e-12 || intErr > 1e-12)
                        {
                            Info<< "MISMATCH face " << facei << " dir " << dirs[di]
                                << " offset " << offsets[oi] << " Un0 " << Un0
                                << " dt " << dt << ": end " << areaEnd
                                << " exact " << exactEnd << " (err " << endErr
                                << "), dVf " << dVf << " exact " << integral
                                << " (err " << intErr << ")" << endl;
                        }
                    }
                }
            }
        }
    }

    Info<< nl << "testCutFaceAdvectEndArea: " << nChecks << " checks on "
        << faces.size() << " faces, max |areaEnd - exact| = " << maxEndErr
        << ", max |dVf - integral(exact)| = " << maxIntErr << nl
        << ((maxEndErr <= 1e-12 && maxIntErr <= 1e-12) ? "PASS" : "FAIL")
        << nl << endl;

    return (maxEndErr <= 1e-12 && maxIntErr <= 1e-12) ? 0 : 1;
}


// ************************************************************************* //
