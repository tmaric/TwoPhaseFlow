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
    Foam::davofNormalModels::cellLocal

Description
    The raw DAVOF normal m_c/|m_c| in every interface cell: the model
    replaces nothing (the plane positioned along m_c stays).

SourceFiles
    cellLocalNormal.C

\*---------------------------------------------------------------------------*/

#ifndef davofNormalModels_cellLocal_H
#define davofNormalModels_cellLocal_H

#include "davofNormalModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
inline namespace twoPhaseFlow
{
namespace davofNormalModels
{

/*---------------------------------------------------------------------------*\
                          Class cellLocal Declaration
\*---------------------------------------------------------------------------*/

class cellLocal
:
    public davofNormalModel
{
public:

    //- Runtime type information
    TypeName("cellLocal");


    //- Construct from the davof dictionary and the scheme
    cellLocal
    (
        const dictionary& dict,
        const reconstruction::davof& scheme
    );


    //- Destructor
    virtual ~cellLocal() = default;


    // Member Functions

        virtual void modelNormals(const scalar tState, vectorField& n);

        virtual bool replaces() const
        {
            return false;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace davofNormalModels
} // End inline namespace twoPhaseFlow
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
