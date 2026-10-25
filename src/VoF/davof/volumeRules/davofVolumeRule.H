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

Class
    Foam::davofVolumeRule

Description
    Run-time selectable rule for the three geometric quantities of the DAVOF
    dual advection with one curved clipping surface per cell
    (isoAdvectionCurved.C): the liquid volume of a cell under its surface
    (the volume anchor and the surface defect), the liquid volume of the
    prism swept between a face's pre-image f - d and the face (the volume
    flux), and the liquid area fraction of a displaced face (the face
    update). The three are the moments of ONE surface, so one rule computes
    all of them.

    Rules (the keyword volumeRule of the isoAdvector dictionary):

        splitQuadrature  the reference (the default, bit-identical to the
                         code before the rules existed): the triangle
                         integrals split at the exact traces, 12-point
                         Gauss-Legendre per piece, adaptive in both sweep
                         directions to curvedTolRel of the volume
                         (davofCurvedGeometry.H)
        kromerBothe      the face-based method of J. Kromer and D. Bothe,
                         J. Comput. Phys. 475 (2023) 111840,
                         doi:10.1016/j.jcp.2022.111840, arXiv:2111.01073:
                         the divergence theorem twice, the paraboloid's
                         patch integral as a line integral along the conic
                         traces on the faces (davofKromerBothe.H); for the
                         paraboloid surface the only error is the
                         one-dimensional quadrature along the arcs. Falls
                         back to splitQuadrature per polyhedron where a
                         configuration is not resolved (counted), and for a
                         sphere surface
        exactSphere      the exact overlap of a sphere with tetrahedra of
                         S. Strobl, A. Formella, T. Poeschel, Exact
                         calculation of the overlap volume of spheres and
                         mesh elements, J. Comput. Phys. 311 (2016) 158-172,
                         doi:10.1016/j.jcp.2016.02.003, through the authors'
                         reference implementation (thirdParty/overlap,
                         MIT): the cell as the tetrahedra (centre, face
                         centre, edge), the prism as three tetrahedra per
                         fan triangle; the face fraction by the quadratic
                         face cut, exact for a circle up to the sliver's
                         Gauss rule. For the sphere translation case; falls
                         back to splitQuadrature for a paraboloid (counted)

    Dictionary:
        volumeRule      splitQuadrature;   // | kromerBothe | exactSphere
        curvedTolRel    1e-11;             // tolerance of the quadratures
        curvedMaxDepth  12;                // splitQuadrature: depth cap

SourceFiles
    davofVolumeRule.C

\*---------------------------------------------------------------------------*/

#ifndef davofVolumeRule_H
#define davofVolumeRule_H

#include "typeInfo.H"
#include "runTimeSelectionTables.H"
#include "autoPtr.H"
#include "dictionary.H"
#include "fvMesh.H"
#include "davofCurvedGeometry.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class davofVolumeRule Declaration
\*---------------------------------------------------------------------------*/

class davofVolumeRule
{
protected:

    // Protected Data

        const fvMesh& mesh_;

        const dictionary dict_;

        //- The tolerance of the quadratures relative to the volume
        const scalar tolRel_;

        //- The depth cap of the split quadrature
        const label maxDepth_;

        //- The split quadrature's tolerance and depth when it is the
        //  FALLBACK of a rule that could not resolve a polyhedron (keys
        //  volumeRuleFallbackTolRel 1e-8, volumeRuleFallbackMaxDepth 6):
        //  unresolved configurations are the wild surfaces of a diverging
        //  loop, where the full-accuracy split quadrature cost hours per step
        const scalar fallbackTolRel_;
        const label fallbackMaxDepth_;

        //- Counters of the step: quadrature leaves (cost) and fallbacks
        mutable label nLeaves_;
        mutable label nFallback_;


    // Protected Member Functions

        //- The capped split quadrature of a fallback
        scalar fallbackCellLiquidVolume
        (
            const label celli,
            const davof::cellSurface& S
        ) const;

        scalar fallbackPrismLiquidVolume
        (
            const label facei,
            const vector& d,
            const davof::cellSurface& S
        ) const;

        //- The reference: the split quadrature of davofCurvedGeometry.H
        scalar splitCellLiquidVolume
        (
            const label celli,
            const davof::cellSurface& S,
            const scalar tolRel
        ) const;

        scalar splitPrismLiquidVolume
        (
            const label facei,
            const vector& d,
            const davof::cellSurface& S,
            const scalar tolRel
        ) const;


public:

    //- Runtime type information
    TypeName("davofVolumeRule");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            davofVolumeRule,
            dictionary,
            (
                const fvMesh& mesh,
                const dictionary& dict
            ),
            (mesh, dict)
        );


    // Constructors

        //- Construct from the isoAdvector dictionary
        davofVolumeRule(const fvMesh& mesh, const dictionary& dict);


    // Selectors

        //- Select by the keyword volumeRule (default splitQuadrature)
        static autoPtr<davofVolumeRule> New
        (
            const fvMesh& mesh,
            const dictionary& dict
        );


    //- Destructor
    virtual ~davofVolumeRule() = default;


    // Member Functions

        //- The liquid volume of celli under S, to tolRel of the cell volume
        virtual scalar cellLiquidVolume
        (
            const label celli,
            const davof::cellSurface& S,
            const scalar tolRel
        ) const = 0;

        //- The liquid volume of celli under S and its slope in the offset
        //  c of the surface, dV/dc = the area of the wetted patch projected
        //  along e3 (the volume anchor's Newton derivative); dVdc is
        //  negative where the rule cannot provide it (the default: the
        //  volume alone, the anchor then brackets)
        virtual scalar cellLiquidVolumeSlope
        (
            const label celli,
            const davof::cellSurface& S,
            const scalar tolRel,
            scalar& dVdc
        ) const
        {
            dVdc = -1;
            return cellLiquidVolume(celli, S, tolRel);
        }

        //- The liquid volume of the prism between facei - d and facei
        virtual scalar prismLiquidVolume
        (
            const label facei,
            const vector& d,
            const davof::cellSurface& S,
            const scalar tolRel
        ) const = 0;

        //- The liquid area fraction of facei displaced by d; nFallback
        //  counts the linear cuts of the quadratic face cut
        virtual scalar faceLiquidFraction
        (
            const label facei,
            const vector& d,
            const davof::cellSurface& S,
            label& nFallback
        ) const;

        //- The wetted area VECTOR of facei displaced by d under S, in the
        //  face's orientation (the sum over the fan triangles about the
        //  face centre of each triangle's wetted area times its unit
        //  normal; alpha_f S_f on a planar face), returned with the scalar
        //  fraction as faceLiquidFraction gives it: the transported face
        //  state on warped faces (faceState vector, 2026-10-25). Default:
        //  the quadratic triangle cut per triangle
        virtual scalar faceWettedArea
        (
            const label facei,
            const vector& d,
            const davof::cellSurface& S,
            vector& aWet,
            label& nFallback
        ) const;

        //- The rule's own line of the step report (nothing by default)
        virtual void report() const {}

        //- The default tolerance of the rule
        scalar tolRel() const { return tolRel_; }

        label maxDepth() const { return maxDepth_; }

        //- The counters of the step
        label nLeaves() const { return nLeaves_; }

        label nFallback() const { return nFallback_; }

        void resetCounters() const { nLeaves_ = 0; nFallback_ = 0; }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
