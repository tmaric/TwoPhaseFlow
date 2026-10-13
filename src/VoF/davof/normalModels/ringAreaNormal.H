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
    Foam::davofNormalModels::ringArea

Description
    The area-weighted average of the ring's DAVOF normals,

        n_c = sum_k m_k / |sum_k m_k|,

    the sum over the point-neighbour ring of the cell (the cells sharing at
    least one point with it, the cell included; zoneDistribute's stencil,
    across processor boundaries), where only the interface cells
    contribute (m_k = 0 otherwise). This is the DAVOF normal of the union
    of the ring's cells: exact for a plane on any mesh, no division by a
    small area. Default applyTo all.

SourceFiles
    ringAreaNormal.C

\*---------------------------------------------------------------------------*/

#ifndef davofNormalModels_ringArea_H
#define davofNormalModels_ringArea_H

#include "davofNormalModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
inline namespace twoPhaseFlow
{
namespace davofNormalModels
{

/*---------------------------------------------------------------------------*\
                           Class ringArea Declaration
\*---------------------------------------------------------------------------*/

class ringArea
:
    public davofNormalModel
{
protected:

    // Protected Data

        //- m_c in the interface cells, zero elsewhere (the field exchanged
        //  across processors; its boundary values are zero, so boundary
        //  faces in the stencil contribute nothing)
        volVectorField mInterface_;


    // Protected Member Functions

        //- Fill mInterface_ from the scheme's area normal
        void fillInterfaceNormals();


public:

    //- Runtime type information
    TypeName("ringArea");


    //- Construct from the davof dictionary and the scheme
    ringArea
    (
        const dictionary& dict,
        const reconstruction::davof& scheme
    );


    //- Destructor
    virtual ~ringArea() = default;


    // Member Functions

        virtual void modelNormals(const scalar tState, vectorField& n);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace davofNormalModels
} // End inline namespace twoPhaseFlow
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
