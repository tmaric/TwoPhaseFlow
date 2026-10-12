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
    Foam::davofInit

Description
    Initialisation of the DAVOF (dual area/volume-of-fluid) state, the cell
    volume fractions alpha_c and the face liquid-area fractions alpha_f, from
    an implicitFunction. Port of leia's davofState constructions
    (src/leiaLevelSet/davof/davofState.C: computeFromImplicitSurface,
    computeFromQuadraticFaces, computeFromExactSphere).

    Sign conventions. TwoPhaseFlow's implicitFunction::value() is positive
    inside the liquid (alpha = 1); leia's simplex geometry takes the liquid
    where psi <= 0, so every formula here is fed psi = -value(). The area
    normal recovered from the face fractions,

        m_c = sum_f alpha_f S_f^out,   S_f^out = +S_f (owner), -S_f (neighbour),

    points INTO the liquid, as TwoPhaseFlow's reconstructionSchemes::normal_
    does (leia's m_c = -sum_f alpha_f S_f^out points out of it).

    Three constructions:

      detrixheAslam   psi at the cell centre, the face centres and the mesh
                      points in the Detrixhe-Aslam formulas: the tets
                      (x_c, x_f, p_i, p_{i+1}) of every cell for the volume,
                      the fan triangles (x_f, p_i, p_{i+1}) with the SAME
                      vertex values for the face area. Exact for a plane (the
                      chord surface IS the plane); on a curved surface the
                      chords are O(h^2) off, the face fractions O(h) as
                      fractions and the normal first order.
      quadraticFaces  detrixheAslam for the cell fractions; the face fractions
                      from the quadratic interpolant of psi at the vertices and
                      the edge midpoints of every fan triangle, cut exactly
                      (davofQuadraticFaceGeometry.H): crossings O(h^3), face
                      fractions O(h^2), the normal second order. Triangles the
                      routine cannot handle keep the Detrixhe-Aslam value
                      (nFaceFallback()).
      exactSphere     the exact fractions of a sphere (origin, radius from the
                      dictionary) on axis-aligned hexahedra: planar faces cut by
                      the ball in closed form, cell fractions by adaptive
                      quadrature of the exact cross-section
                      (davofSphereGeometry.H).

    The Detrixhe-Aslam formulas are exact on the chord surface only when psi is
    a signed distance. The dictionary switch `signedDistance true` replaces
    -value() by the signed distance obtained by projecting the point onto the
    zero set of the implicit function (Newton along the gradient, then the
    foot point), for functions whose value is algebraic (the ellipsoid).

    The face field is filled for every polyMesh face: the internal field and,
    patch-wise, the boundary field (empty patches have size zero and are
    skipped).

SourceFiles
    davofInit.C

\*---------------------------------------------------------------------------*/

#ifndef davofInit_H
#define davofInit_H

#include "fvMesh.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "dictionary.H"
#include "implicitFunction.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
inline namespace twoPhaseFlow
{

/*---------------------------------------------------------------------------*\
                          Class davofInit Declaration
\*---------------------------------------------------------------------------*/

class davofInit
{
    // Private Data

        const fvMesh& mesh_;

        const implicitFunction& func_;

        const dictionary& dict_;

        //- Replace -value() by the signed distance (Newton projection)
        const bool signedDistance_;

        //- Largest |raw - clipped| of a cell fraction
        scalar maxAlphaClip_;

        //- quadraticFaces: faces kept at the Detrixhe-Aslam value
        label nFaceFallback_;

        //- Sum of the tet zero-set area vectors per cell, oriented into the
        //  liquid (the reference of the consistency check; the chord
        //  surface's for every construction)
        vectorField mTet_;


    // Private Member Functions

        //- psi = -value() (or minus the signed distance)
        scalar psi(const point& x) const;

        //- Detrixhe-Aslam cell fraction of celli from the vertex values; sets
        //  mTet_[celli]. Returns the clipped fraction.
        scalar cellFromValues
        (
            const label celli,
            const scalarField& psiC,
            const scalarField& psiF,
            const scalarField& psiP
        );

        //- Detrixhe-Aslam fraction of facei (fan about the face centre)
        scalar faceFraction
        (
            const label facei,
            const scalarField& psiF,
            const scalarField& psiP
        ) const;

        //- The Detrixhe-Aslam state: alpha (clipped), one value per polyMesh
        //  face in alphafAll
        void detrixheAslamState
        (
            volScalarField& alpha,
            scalarField& alphafAll,
            scalarField& psiF,
            scalarField& psiP
        );

        //- Copy the per-face values into the surfaceScalarField (internal
        //  field, then patch-wise; empty patches have size zero)
        void fillFaceField
        (
            const scalarField& alphafAll,
            surfaceScalarField& alphaf
        ) const;

        //- No copy construct
        davofInit(const davofInit&) = delete;

        //- No copy assignment
        void operator=(const davofInit&) = delete;


public:

    // Constructors

        //- Construct from mesh, implicit function and dictionary
        davofInit
        (
            const fvMesh& mesh,
            const implicitFunction& func,
            const dictionary& dict
        );


    //- Destructor
    ~davofInit() = default;


    // Member Functions

        //- The Detrixhe-Aslam state
        void detrixheAslam(volScalarField& alpha, surfaceScalarField& alphaf);

        //- Detrixhe-Aslam cell fractions, quadratic face fractions
        void quadraticFaces(volScalarField& alpha, surfaceScalarField& alphaf);

        //- The exact state of the sphere (origin, radius) on axis-aligned
        //  hexahedra; FatalError on any other cell
        void exactSphere(volScalarField& alpha, surfaceScalarField& alphaf);

        //- Interface-cell threshold on |m_c|: wispTol V_c^(2/3)
        static scalar wispThreshold
        (
            const fvMesh& mesh,
            const label celli,
            const scalar wispTol
        );

        //- Distance from x to the zero set of func by Newton projection
        //  (along the gradient onto the surface, then the foot point by
        //  tangential moves); foot receives the surface point.
        static scalar projectedDistance
        (
            const implicitFunction& func,
            const point& x,
            point& foot
        );

        //- Signed distance with the sign of value(): positive inside the
        //  liquid
        static scalar signedDistance
        (
            const implicitFunction& func,
            const point& x
        );

        //- m_c = sum_f alpha_f S_f^out per cell (into the liquid), from the
        //  face field; alphafAll receives the per-polyMesh-face values
        static void areaNormal
        (
            const fvMesh& mesh,
            const surfaceScalarField& alphaf,
            vectorField& m,
            scalarField& alphafAll
        );

        //- max over the interface cells (|m_c| > wispTol V^(2/3)) of
        //  |m_c - mTet_c| / V_c^(2/3), reduced over the ranks
        scalar maxConsistency
        (
            const surfaceScalarField& alphaf,
            const scalar wispTol
        ) const;

        scalar maxAlphaClip() const
        {
            return maxAlphaClip_;
        }

        label nFaceFallback() const
        {
            return nFaceFallback_;
        }

        const vectorField& mTet() const
        {
            return mTet_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End inline namespace twoPhaseFlow
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
