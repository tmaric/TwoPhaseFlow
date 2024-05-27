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
    Foam::gravityDirac

Description

SourceFiles
    gravityDirac.C
    newPhaseChangeModel.C

\*---------------------------------------------------------------------------*/

#ifndef gravityDirac_H
#define gravityDirac_H

#include "typeInfo.H"
#include "volFields.H"
#include "dimensionedScalar.H"
#include "autoPtr.H"
#include "zoneDistribute.H"
#include "accelerationForceModel.H"
#include "markInterfaceRegion.H"
// #include "immiscibleIncompressibleTwoPhaseMixture.H"
#include "incompressibleTwoPhaseMixture.H"
#include "cutFacePLIC.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class gravityDirac Declaration
\*---------------------------------------------------------------------------*/

class gravityDirac
    : public accelerationForceModel
{

private:

    //-
        //- Stabilisation for normalisation of the interface normal
        const dictionary& gravityDict_;
        dimensionedVector g_;
//        surfaceScalarField snGradAlpha_;

        //- Get density jump
        incompressibleTwoPhaseMixture mixture_;
        surfaceScalarField snGradAlpha_;
        //- Cut face to get cut edge
        cutFacePLIC  cutFace_;
    // Private Member Functions

        //- Disallow copy construct
        gravityDirac(const gravityDirac&);

        //- Disallow default bitwise assignment
        void operator=(const gravityDirac&);

//        vector closestDist(const point p, const vector n ,const vector center);


        protected:

        scalar calculateDiracInteg (
            const label faceID,
            scalar accf,
            const scalar snGradAlpha,
            const scalar rmagSf,
            const vector Sf,
            const vector faceNormal,
            const vector faceCenter,
            const vector cellFaceNormal
        );
        //- Re-calculate the acceleration
        virtual void calculateAcc();

public:

    //- Runtime type information
    TypeName("gravityDirac");


    // Constructors

        //- Construct from components
        gravityDirac
        (
            const dictionary& dict,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~gravityDirac()
    {}

    virtual tmp<surfaceScalarField> accelerationForce();


};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
