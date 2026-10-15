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

Class
    Foam::davofFaceClippingModels::planeSagitta

Description
    The sagitta correction of the plane clipping: the stored face area is
    the plane-clipped area plus the second-order area between the plane's
    trace on the face (the chord, end points p1, p2, length L, direction
    t = (p2 - p1)/L) and the interface's trace (the arc), from the
    curvature tensor B the davof scheme fits to the ring
    (davof::updateCurvatureTensor, the tensor of the start-of-step planes
    used for the whole step).

    Conventions. n is the plane's unit normal INTO the liquid (isoAdvector's
    n0 = -n points out of it; isoAdvection passes -n0). B is the Hessian of
    the Monge patch z = z_0 + (1/2) x.B.x of the interface over the plane
    (z along n, x in the plane from the plane's centroid c): t.B.t > 0
    where the surface bends TOWARD the liquid (a liquid sphere: B = I_t/R).
    The plane is positioned by the cell's volume fraction, so the
    paraboloid's vertex sits at z_0 = -(1/2) B : M above the plane (M the
    cut polygon's second moments about c): the plane is a secant of the
    surface over the cell, not its tangent plane.

    The correction. On the face with unit normal n_f the trace of the
    surface is displaced from the chord, in the face plane and along the
    chord's in-plane normal nu = (n - (n.n_f) n_f)/cos(phi) that points
    into the liquid, by w/cos(phi), where w(l) = z_0 + (1/2) x(l).B.x(l) is
    the surface's offset from the plane along n at the chord point
    x(l) = (p1 - c') + l t, l in [0, L], c' = c + U dt the centroid at the
    end of the step, and cos(phi) = |n - (n.n_f) n_f| (Meusnier: the
    trace's curvature is the normal curvature over cos(phi)). Where w > 0
    the surface lies on the liquid side of the plane and the true liquid
    area is smaller than the plane's, so

        dA = -(1/cos(phi)) int_0^L w(l) dl
           = -(1/cos(phi)) [ z_0 L + (1/2)( a.B.a L + a.B.t L^2
                                          + t.B.t L^3/3 ) ],   a = p1 - c',

    exact to O(h^4) in the area for a second-order normal and a
    first-order B (the end triangles where chord and arc meet the face's
    edges are O(h^4)); the stored fraction is then second order. On a
    cube cut through its centre by a plane normal to e3 with B = diag(B11,
    B22) this gives dA = -B11 h^3/12 on the faces normal to e1 and
    -B22 h^3/12 on those normal to e2: the plane OVERESTIMATES the wetted
    area of every side face of a convex liquid region, because the secant
    plane lies above the vertex by z_0 = -(B11 + B22) h^2/24.

    sagittaForm segment restricts the correction to the parabolic segment
    between the chord and the arc through the chord's end points,

        dA = + t.B.t L^3 / (12 cos(phi)) = (2/3) L s / cos(phi),
        s = t.B.t L^2 / 8 (the sagitta),

    the original formulation, which assumes the chord's end points lie on
    the arc; it drops the offset z_0 and the a-terms, which are of the
    same order (on the cube above it gives +B22 h^3/12 on the faces normal
    to e1: the opposite sign). Default full.

    Guards: no correction without a chord (the face on one side of the
    plane at the end of the step), in a cell whose tensor fit fell back
    (B = 0), when cos(phi) < sagittaMinCos (the interface nearly parallel
    to the face, the trace ill-conditioned; counted "skipped (cos)"),
    when L < sagittaMinChord sqrt(|S_f|) (the plane grazes a corner;
    "skipped (chord)"); |dA| is capped at sagittaMaxFraction |S_f|
    ("capped") and the stored area is clamped to [0, |S_f|] ("clamped").
    Planar faces: the chord is the plane's trace at the end of the step
    on the face's edges (cutFaceAdvect::chordEnd()); on a non-planar face
    the crossings with the face polygon's edges are used, with the face's
    mean normal, an approximation.

        faceClipping        planeSagitta;
        sagittaMinCos       0.2;
        sagittaMinChord     0.05;
        sagittaMaxFraction  0.25;
        sagittaForm         full;     // full | segment
        tensorSource        ring;     // the ring Monge fit (the only source)

SourceFiles
    planeSagittaClipping.C

\*---------------------------------------------------------------------------*/

#ifndef davofFaceClippingModels_planeSagitta_H
#define davofFaceClippingModels_planeSagitta_H

#include "davofFaceClipping.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
inline namespace twoPhaseFlow
{
namespace davofFaceClippingModels
{

/*---------------------------------------------------------------------------*\
                         Class planeSagitta Declaration
\*---------------------------------------------------------------------------*/

class planeSagitta
:
    public davofFaceClipping
{
    // Private Data

        //- Skip when cos(phi) is below this
        scalar minCos_;

        //- Skip when the chord is shorter than this times sqrt(|S_f|)
        scalar minChord_;

        //- Cap |dA| at this fraction of |S_f|
        scalar maxFraction_;

        //- full | segment
        word form_;

        //- ring (the only source of the tensor in this round)
        word tensorSource_;


public:

    //- Runtime type information
    TypeName("planeSagitta");


    //- Construct from the alpha solver dictionary and the scheme
    planeSagitta
    (
        const dictionary& dict,
        const reconstruction::davof& scheme
    );


    //- Destructor
    virtual ~planeSagitta() = default;


    // Member Functions

        virtual bool corrects() const
        {
            return true;
        }

        virtual bool needsCurvatureTensor() const
        {
            return true;
        }

        virtual scalar correct
        (
            const label facei,
            const label planeCell,
            const vector& nLiq,
            const point& xEnd,
            const UList<point>& chord,
            const scalar areaPlane
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace davofFaceClippingModels
} // End inline namespace twoPhaseFlow
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
