/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2021 OpenCFD Ltd.
    Copyright (C) YEAR AUTHOR, AFFILIATION
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

\*---------------------------------------------------------------------------*/

#include "fixedValueFvPatchFieldTemplate.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "unitConversion.H"
#include "PatchFunction1.H"

//{{{ begin codeInclude

//}}} end codeInclude


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

//{{{ begin localCode

//}}} end localCode


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

// dynamicCode:
// SHA1 = aa5abfeb8f0c91250404a6f90d7332a5cfa350ca
//
// unique function name that can be checked if the correct library version
// has been loaded
extern "C" void codedPatchBC_aa5abfeb8f0c91250404a6f90d7332a5cfa350ca(bool load)
{
    if (load)
    {
        // Code that can be explicitly executed after loading
    }
    else
    {
        // Code that can be explicitly executed before unloading
    }
}

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

makeRemovablePatchTypeField
(
    fvPatchVectorField,
    codedPatchBCFixedValueFvPatchVectorField
);

} // End namespace Foam


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::
codedPatchBCFixedValueFvPatchVectorField::
codedPatchBCFixedValueFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    parent_bctype(p, iF)
{
    if (false)
    {
        printMessage("Construct codedPatchBC : patch/DimensionedField");
    }
}


Foam::
codedPatchBCFixedValueFvPatchVectorField::
codedPatchBCFixedValueFvPatchVectorField
(
    const codedPatchBCFixedValueFvPatchVectorField& rhs,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    parent_bctype(rhs, p, iF, mapper)
{
    if (false)
    {
        printMessage("Construct codedPatchBC : patch/DimensionedField/mapper");
    }
}


Foam::
codedPatchBCFixedValueFvPatchVectorField::
codedPatchBCFixedValueFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    parent_bctype(p, iF, dict)
{
    if (false)
    {
        printMessage("Construct codedPatchBC : patch/dictionary");
    }
}


Foam::
codedPatchBCFixedValueFvPatchVectorField::
codedPatchBCFixedValueFvPatchVectorField
(
    const codedPatchBCFixedValueFvPatchVectorField& rhs
)
:
    parent_bctype(rhs),
    dictionaryContent(rhs)
{
    if (false)
    {
        printMessage("Copy construct codedPatchBC");
    }
}


Foam::
codedPatchBCFixedValueFvPatchVectorField::
codedPatchBCFixedValueFvPatchVectorField
(
    const codedPatchBCFixedValueFvPatchVectorField& rhs,
    const DimensionedField<vector, volMesh>& iF
)
:
    parent_bctype(rhs, iF)
{
    if (false)
    {
        printMessage("Construct codedPatchBC : copy/DimensionedField");
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::
codedPatchBCFixedValueFvPatchVectorField::
~codedPatchBCFixedValueFvPatchVectorField()
{
    if (false)
    {
        printMessage("Destroy codedPatchBC");
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void
Foam::
codedPatchBCFixedValueFvPatchVectorField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    if (false)
    {
        printMessage("updateCoeffs codedPatchBC");
    }

//{{{ begin code
    #line 40 "/work/groups/da_mma_b/jun/TwoPhaseFlow/run/rhoVoF/jetCrossFlow/0/U.boundaryField.jetInlet"
const fvPatch& patch = this->patch();
            const vectorField& cf = patch.Cf();
            vectorField& Ufield = *this;
            const vector jetCenter = vector(0.0, 0.0, 0.0);

            forAll(cf, i)
            {
                const scalar dist = mag(cf[i]-jetCenter);
                Ufield[i] = vector(0, -21.434*pow(dist/0.001,3)+15.512*pow(dist/0.001,2)-2.5722*(dist/0.001)+8.6504, 0);
            }
//}}} end code

    this->parent_bctype::updateCoeffs();
}


// ************************************************************************* //

