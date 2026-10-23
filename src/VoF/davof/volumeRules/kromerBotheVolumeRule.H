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
    Foam::davofVolumeRules::kromerBothe

Description
    The face-based volume rule of

        J. Kromer, D. Bothe, Third-order accurate initialization of volume
        fractions on unstructured meshes with arbitrary polyhedral cells,
        Journal of Computational Physics 475 (2023) 111840,
        doi:10.1016/j.jcp.2022.111840, arXiv:2111.01073,

    applied to the DAVOF paraboloid, for which their local paraboloid is the
    surface itself (davofKromerBothe.H): the cell volume and the prism
    volume as face sums plus a line integral along the conic traces, the face
    fraction by the Green integral around the liquid region of the face.
    Where a face is not resolved the polyhedron falls back to the split
    quadrature and is counted; a sphere surface falls back as a whole
    (exactSphere is the rule for it).

    Dictionary (isoAdvector): volumeRule kromerBothe; curvedTolRel (the
    tolerance of the arc quadratures relative to the volume or the area).

SourceFiles
    kromerBotheVolumeRule.C

\*---------------------------------------------------------------------------*/

#ifndef kromerBotheVolumeRule_H
#define kromerBotheVolumeRule_H

#include "davofVolumeRule.H"

namespace Foam
{
namespace davofVolumeRules
{

/*---------------------------------------------------------------------------*\
                         Class kromerBothe Declaration
\*---------------------------------------------------------------------------*/

class kromerBothe
:
    public davofVolumeRule
{
    // Private Data

        //- Compare every prism and cell volume with the split quadrature
        //  and print the first mismatches (key volumeRuleDebug)
        const bool debug_;
        mutable label nDebugPrinted_;


public:

    //- Runtime type information
    TypeName("kromerBothe");


    // Constructors

        kromerBothe(const fvMesh& mesh, const dictionary& dict);


    //- Destructor
    virtual ~kromerBothe() = default;


    // Member Functions

        virtual scalar cellLiquidVolume
        (
            const label celli,
            const davof::cellSurface& S,
            const scalar tolRel
        ) const;

        virtual scalar prismLiquidVolume
        (
            const label facei,
            const vector& d,
            const davof::cellSurface& S,
            const scalar tolRel
        ) const;

        virtual scalar faceLiquidFraction
        (
            const label facei,
            const vector& d,
            const davof::cellSurface& S,
            label& nFallback
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace davofVolumeRules
} // End namespace Foam

#endif

// ************************************************************************* //
