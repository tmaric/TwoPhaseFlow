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
    The position at time t of the origin of an exact body (the sphere of the
    exactNormal dictionaries) that is carried by the uniform velocity U0 and
    turned about the axis through axisPoint with the angular velocity omega
    (rad/s; a zero vector: no rotation), for the solid-body rotation
    benchmarks (2026-10-28):

        x(t) = axisPoint + R(omega t) (x0 - axisPoint) + U0 t

    with R the rotation by |omega| t about omega/|omega|, by the rotation
    formula of O. Rodrigues, Des lois geometriques qui regissent les
    deplacements d'un systeme solide dans l'espace, et de la variation des
    coordonnees provenant de ces deplacements consideres independamment des
    causes qui peuvent les produire, Journal de Mathematiques Pures et
    Appliquees 5 (1840) 380-440 (no DOI; the formula in its vector form,
    R v = v cos(theta) + (k x v) sin(theta) + k (k.v)(1 - cos(theta))).
    Used by the exact normal model, the exactSphere clipping surface and the
    areaFractionError function object, so that the three see the same body.

SourceFiles
    (header only)

\*---------------------------------------------------------------------------*/

#ifndef davofOrbit_H
#define davofOrbit_H

#include "point.H"
#include "vector.H"
#include "dictionary.H"

namespace Foam
{
namespace davof
{

//- The origin at time t (see the file header)
inline point orbitOrigin
(
    const point& x0,
    const vector& U0,
    const vector& omega,
    const point& axisPoint,
    const scalar t
)
{
    const scalar w = mag(omega);
    point x = x0;
    if (w > VSMALL)
    {
        const vector k = omega/w;
        const scalar th = w*t;
        const vector r = x0 - axisPoint;
        x = axisPoint
          + r*Foam::cos(th)
          + (k ^ r)*Foam::sin(th)
          + k*(k & r)*(1 - Foam::cos(th));
    }
    return x + U0*t;
}


//- The rotation keys of an exact-body dictionary: omega (default zero) and
//  axisPoint (default the origin at t = 0, i.e. a rotation about the body's
//  own centre, which does not move it)
inline void readOrbit
(
    const dictionary& dict,
    const point& origin0,
    vector& omega,
    point& axisPoint
)
{
    omega = dict.getOrDefault<vector>("omega", vector::zero);
    axisPoint = dict.getOrDefault<point>("axisPoint", origin0);
}

} // End namespace davof
} // End namespace Foam

#endif

// ************************************************************************* //
