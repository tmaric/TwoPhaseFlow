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

#include "exactSphereVolumeRule.H"
#include "addToRunTimeSelectionTable.H"

// The reference implementation of Strobl, Formella and Poeschel (2016),
// doi:10.1016/j.jcp.2016.02.003 (MIT; Eigen 3)
#include "overlap/overlap.hpp"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace davofVolumeRules
{
    defineTypeNameAndDebug(exactSphere, 0);
    addToRunTimeSelectionTable(davofVolumeRule, exactSphere, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::davofVolumeRules::exactSphere::exactSphere
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    davofVolumeRule(mesh, dict)
{}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::scalar Foam::davofVolumeRules::exactSphere::tetOverlap
(
    const davof::cellSurface& S,
    const point& a, const point& b, const point& c, const point& d
) const
{
    // a positively oriented tetrahedron (the library asserts the triple
    // product in debug builds)
    const scalar tp = ((b - a) ^ (c - a)) & (d - a);
    if (mag(tp) <= VSMALL) return 0;
    auto V = [](const point& p) { return overlap::Vector(p.x(), p.y(), p.z()); };
    const overlap::Sphere sph(V(S.o), S.R);
    if (tp > 0)
    {
        const overlap::Tetrahedron tet(V(a), V(b), V(c), V(d));
        return overlap::overlap_volume(sph, tet);
    }
    const overlap::Tetrahedron tet(V(a), V(c), V(b), V(d));
    return overlap::overlap_volume(sph, tet);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::davofVolumeRules::exactSphere::cellLiquidVolume
(
    const label celli,
    const davof::cellSurface& S,
    const scalar tolRel
) const
{
    if (!S.sphere)
    {
        ++nFallback_;
        return splitCellLiquidVolume(celli, S, tolRel);
    }
    const pointField& P = mesh_.points();
    const faceList& faces = mesh_.faces();
    const vectorField& Cf = mesh_.faceCentres();
    const point& xc = mesh_.cellCentres()[celli];
    const labelList& cFaces = mesh_.cells()[celli];
    scalar v = 0;
    forAll(cFaces, cf)
    {
        const label faceI = cFaces[cf];
        const face& f = faces[faceI];
        const point& xf = Cf[faceI];
        forAll(f, ip)
        {
            const label ip1 = f.nextLabel(ip);
            v += tetOverlap(S, xc, xf, P[f[ip]], P[ip1]);
            ++nLeaves_;
        }
    }
    return v;
}


Foam::scalar Foam::davofVolumeRules::exactSphere::prismLiquidVolume
(
    const label facei,
    const vector& d,
    const davof::cellSurface& S,
    const scalar tolRel
) const
{
    if (mag(d & mesh_.faceAreas()[facei]) <= VSMALL) return 0;
    if (!S.sphere)
    {
        ++nFallback_;
        return splitPrismLiquidVolume(facei, d, S, tolRel);
    }
    const pointField& P = mesh_.points();
    const face& f = mesh_.faces()[facei];
    const point& xf = mesh_.faceCentres()[facei];
    scalar v = 0;
    forAll(f, ip)
    {
        const label ip1 = f.nextLabel(ip);
        // the triangular prism (A, B, C) -> (A', B', C') = (A, B, C) - d as
        // the three tetrahedra (A, B, C, A'), (B, C, A', B'), (C, A', B', C')
        const point A = xf, B = P[f[ip]], C = P[ip1];
        const point Ap = A - d, Bp = B - d, Cp = C - d;
        v += tetOverlap(S, A, B, C, Ap);
        v += tetOverlap(S, B, C, Ap, Bp);
        v += tetOverlap(S, C, Ap, Bp, Cp);
        nLeaves_ += 3;
    }
    return v;
}


// ************************************************************************* //
