/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022 Tomislav Maric TU Darmstadt 
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
    Foam::reconstruction::areaFraction

Description


SourceFiles
    areaFraction.C

\*---------------------------------------------------------------------------*/

#ifndef areaFraction_H
#define areaFraction_H

#include "gradAlpha.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace reconstruction
{

/*---------------------------------------------------------------------------*\
                          Class areaFraction Declaration
\*---------------------------------------------------------------------------*/

class areaFraction
:
    public gradAlpha 
{

protected: 

    const fvMesh& mesh_;

    const Time& runTime_; 
    
    //- FIXME(TM): this should be protected in gradAlpha, or const-access to 
    //  it should be provided. 
    scalar surfCellTol_;
    
    surfaceScalarField alphaf_; 

public:

    //- Runtime type information
    TypeName("areaFraction");

    //- Construct from components
    areaFraction
    (
        volScalarField& alpha1,
        const surfaceScalarField& phi,
        const volVectorField& U,
        const dictionary& dict
    );

    //- Destructor
    virtual ~areaFraction() = default;

    // Member Functions
    
        //- Compute area fractions given (extrenal) data. 
        virtual void calcAreaFractions
        (
            const volScalarField& alpha1,
            const volVectorField& centres,
            const volVectorField& normals, 
            const boolList& interfaceCell
        );

        //- Compute area fractions using internal data.
        virtual void calcAreaFractions(); 

        const surfaceScalarField& areaFractions() const
        {
            return alphaf_;
        }
        
        //- Reconstruct interface
        //- TODO(TM): Youngs + area fractions from Youngs + normal calculation.
        virtual void reconstruct(bool forceUpdate = true);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace reconstruction
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
