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
    Foam::davofNormalModels::ringSum

Description
    The normal of "pure" DAVOF (Tomislav's decision, 2026-11-08): the ring
    sum d_c = sum_j m_j over the point-neighbour ring of the cell (the cell
    included) of EVERY ring cell, as the interface-cell rule alphaAndRingSum
    computes it (davof::ringSum()); un-normalised, the scheme's positioning
    normalises it or hands it on as it is (positioningDirection). Under any
    other interface-cell rule the ring sum is zero and the model replaces
    nothing (the cell keeps m_c, counted as kept). Differs from ringArea in
    the cells that contribute: ringArea sums the interface cells' m only.

SourceFiles
    ringSumNormal.C

\*---------------------------------------------------------------------------*/

#ifndef davofNormalModels_ringSum_H
#define davofNormalModels_ringSum_H

#include "davofNormalModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
inline namespace twoPhaseFlow
{
namespace davofNormalModels
{

/*---------------------------------------------------------------------------*\
                           Class ringSum Declaration
\*---------------------------------------------------------------------------*/

class ringSum
:
    public davofNormalModel
{
public:

    //- Runtime type information
    TypeName("ringSum");


    // Constructors

        //- Construct from the davof dictionary and the scheme
        ringSum
        (
            const dictionary& dict,
            const reconstruction::davof& scheme
        );


    //- Destructor
    virtual ~ringSum() = default;


    // Member Functions

        //- The ring sum of every interface cell (zero under another rule)
        virtual void modelNormals(const scalar tState, vectorField& n);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace davofNormalModels
} // End inline namespace twoPhaseFlow
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
