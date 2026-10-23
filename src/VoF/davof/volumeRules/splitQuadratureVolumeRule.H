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
    Foam::davofVolumeRules::splitQuadrature

Description
    The reference volume rule: the triangle integrals of the prism and of
    the pyramids from the cell centre, split at the exact traces of the
    quadric, 12-point Gauss-Legendre per piece, adaptive in both sweep
    directions (davofCurvedGeometry.H). Any quadric surface.

SourceFiles
    splitQuadratureVolumeRule.C

\*---------------------------------------------------------------------------*/

#ifndef splitQuadratureVolumeRule_H
#define splitQuadratureVolumeRule_H

#include "davofVolumeRule.H"

namespace Foam
{
namespace davofVolumeRules
{

/*---------------------------------------------------------------------------*\
                       Class splitQuadrature Declaration
\*---------------------------------------------------------------------------*/

class splitQuadrature
:
    public davofVolumeRule
{
public:

    //- Runtime type information
    TypeName("splitQuadrature");


    // Constructors

        splitQuadrature(const fvMesh& mesh, const dictionary& dict);


    //- Destructor
    virtual ~splitQuadrature() = default;


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
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace davofVolumeRules
} // End namespace Foam

#endif

// ************************************************************************* //
