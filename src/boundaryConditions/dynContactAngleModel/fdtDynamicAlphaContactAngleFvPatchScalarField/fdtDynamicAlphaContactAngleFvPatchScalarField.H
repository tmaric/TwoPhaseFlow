/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2013 OpenFOAM Foundation
    Copyright (C) 2019 OpenCFD Ltd.
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
    Foam::fdtDynamicAlphaContactAngleFvPatchScalarField

Description

    A dynamic contact angle model that implements hysteresis using the 

    Feedback Deceleration Technique 

    Krämer, Veronika, Beawer Barwari, Sebastian Burgmann, Martin Rohde, Simon
    Rentschler, Christopher Holzknecht, Christoph Gmelin, and Uwe Janoske.
    "Numerical analysis of an adhering droplet applying an adapted feedback
    deceleration technique." International Journal of Multiphase Flow 145
    (2021): 103808., https://doi.org/10.1016/j.ijmultiphaseflow.2021.103808

SourceFiles
    fdtDynamicAlphaContactAngleFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef fdtDynamicAlphaContactAngleFvPatchScalarField_H
#define fdtDynamicAlphaContactAngleFvPatchScalarField_H

#include "alphaContactAngleTwoPhaseFvPatchScalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
               Class dynamicAlphaContactAngleFvPatch Declaration
\*---------------------------------------------------------------------------*/

class fdtDynamicAlphaContactAngleFvPatchScalarField
:
    public alphaContactAngleTwoPhaseFvPatchScalarField
{
    // Private data

        //- Equilibrium contact angle
        scalar theta0_;

        //- Limiting advancing contact angle
        scalar thetaA_;

        //- Limiting receding contact angle
        scalar thetaR_;

        //- Scaling factor
        scalar dxdy_;


   // Private member functions
       
       //- Determine whether a face contains a contact line
       bool hasContactLine(label faceI) const;

public:

    //- Runtime type information
    TypeName("fdtHysteresis");


    // Constructors

        //- Construct from patch and internal field
        fdtDynamicAlphaContactAngleFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        fdtDynamicAlphaContactAngleFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  fdtDynamicAlphaContactAngleFvPatchScalarField
        //  onto a new patch
        fdtDynamicAlphaContactAngleFvPatchScalarField
        (
            const fdtDynamicAlphaContactAngleFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        fdtDynamicAlphaContactAngleFvPatchScalarField
        (
            const fdtDynamicAlphaContactAngleFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new fdtDynamicAlphaContactAngleFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        fdtDynamicAlphaContactAngleFvPatchScalarField
        (
            const fdtDynamicAlphaContactAngleFvPatchScalarField&,
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
                new fdtDynamicAlphaContactAngleFvPatchScalarField(*this, iF)
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
