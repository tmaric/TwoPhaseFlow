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
    Foam::gravityVDirac

Description

SourceFiles
    gravityVDirac.C
    newPhaseChangeModel.C

\*---------------------------------------------------------------------------*/

#ifndef gravityVDirac_H
#define gravityVDirac_H

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
                Class gravityVDirac Declaration
\*---------------------------------------------------------------------------*/

class gravityVDirac
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
//        surfaceScalarField snGradAlpha_;
        volVectorField diracV_;
        surfaceVectorField diracVf_;
        //- Cut face to get cut edge
//        cutFacePLIC  cutFace_;
    // Private Member Functions

        //- Disallow copy construct
        gravityVDirac(const gravityVDirac&);

        //- Disallow default bitwise assignment
        void operator=(const gravityVDirac&);

//        vector closestDist(const point p, const vector n ,const vector center);


        protected:

        //- Re-calculate the acceleration
        virtual void calculateAcc();

public:

    //- Runtime type information
    TypeName("gravityVDirac");


    // Constructors

        //- Construct from components
        gravityVDirac
        (
            const dictionary& dict,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~gravityVDirac()
    {}

    virtual tmp<surfaceScalarField> accelerationForce();


};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
