/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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
    Foam::dynamicAlphaContactAngleVoxVoinovFvPatchScalarField

Description
    A dynamic alphaContactAngle scalar boundary condition
    (alphaContactAngleTwoPhaseFvPatchScalarField)

SourceFiles
    dynamicAlphaContactAngleVoxVoinovFvPatchScalarField.C



\*---------------------------------------------------------------------------*/

#ifndef dynamicAlphaContactAngleVoxVoinovFvPatchScalarField_H
#define dynamicAlphaContactAngleVoxVoinovFvPatchScalarField_H

#include "alphaContactAngleTwoPhaseFvPatchScalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class dynamicAlphaContactAngleVoxVoinovFvPatch Declaration
\*---------------------------------------------------------------------------*/

class dynamicAlphaContactAngleVoxVoinovFvPatchScalarField
:
    public alphaContactAngleTwoPhaseFvPatchScalarField
{
    // Private data

        //- Equilibrium contact angle
        scalar theta0_;

        //- Coeffizient
        scalar ct_;


public:

    //- Runtime type information
    TypeName("dynamicAlphaContactAngleVoxVoinov");


    // Constructors

        //- Construct from patch and internal field
        dynamicAlphaContactAngleVoxVoinovFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        dynamicAlphaContactAngleVoxVoinovFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  dynamicAlphaContactAngleVoxVoinovFvPatchScalarField
        //  onto a new patch
        dynamicAlphaContactAngleVoxVoinovFvPatchScalarField
        (
            const dynamicAlphaContactAngleVoxVoinovFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        dynamicAlphaContactAngleVoxVoinovFvPatchScalarField
        (
            const dynamicAlphaContactAngleVoxVoinovFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new dynamicAlphaContactAngleVoxVoinovFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        dynamicAlphaContactAngleVoxVoinovFvPatchScalarField
        (
            const dynamicAlphaContactAngleVoxVoinovFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new dynamicAlphaContactAngleVoxVoinovFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

        //- Evaluate and return dynamic contact-angle
        virtual tmp<scalarField> theta
        (
            const fvPatchVectorField& Up,
            const fvsPatchVectorField& nHat
        ) const;

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
