/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
     \\/     M anipulation  |
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
    Foam::gravityReconImprove

Description

SourceFiles
    gravityReconImprove.C
    newPhaseChangeModel.C

\*---------------------------------------------------------------------------*/

#ifndef gravityReconImprove_H
#define gravityReconImprove_H

#include "typeInfo.H"
#include "volFields.H"
#include "dimensionedScalar.H"
#include "autoPtr.H"
#include "zoneDistribute.H"
#include "accelerationForceModel.H"
#include "markInterfaceRegion.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class gravityReconImprove Declaration
\*---------------------------------------------------------------------------*/

class gravityReconImprove
    : public accelerationForceModel
{

private:

    //-
        //- Stabilisation for normalisation of the interface normal
        const dictionary& gravityDict_;
        dimensionedVector g_;

    // Private Member Functions

        //- Disallow copy construct
        gravityReconImprove(const gravityReconImprove&);

        //- Disallow default bitwise assignment
        void operator=(const gravityReconImprove&);

        vector closestDist(const point p, const vector n ,const vector center);


        protected:

        //- Re-calculate the acceleration
        virtual void calculateAcc();

public:

    //- Runtime type information
    TypeName("gravityReconImprove");


    // Constructors

        //- Construct from components
        gravityReconImprove
        (
            const dictionary& dict,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~gravityReconImprove()
    {}

    virtual tmp<surfaceScalarField> accelerationForce();


};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
