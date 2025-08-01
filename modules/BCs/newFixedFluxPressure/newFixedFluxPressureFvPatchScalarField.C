/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2015-2020 OpenCFD Ltd.
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

#include "newFixedFluxPressureFvPatchScalarField.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::newFixedFluxPressureFvPatchScalarField::newFixedFluxPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(p, iF),
    curTimeIndex_(-1)
{}


Foam::newFixedFluxPressureFvPatchScalarField::newFixedFluxPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedGradientFvPatchScalarField(p, iF),
    curTimeIndex_(-1)
{
    patchType() = dict.getOrDefault<word>("patchType", word::null);
    if (dict.found("value") && dict.found("gradient"))
    {
        fvPatchField<scalar>::operator=
        (
            scalarField("value", dict, p.size())
        );
        gradient() = scalarField("gradient", dict, p.size());
    }
    else
    {
        fvPatchField<scalar>::operator=(patchInternalField());
        gradient() = 0.0;
    }
}


Foam::newFixedFluxPressureFvPatchScalarField::newFixedFluxPressureFvPatchScalarField
(
    const newFixedFluxPressureFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedGradientFvPatchScalarField(p, iF),
    curTimeIndex_(-1)
{
    patchType() = ptf.patchType();

    // Map gradient. Set unmapped values and overwrite with mapped ptf
    gradient() = 0.0;
    gradient().map(ptf.gradient(), mapper);

    // Evaluate the value field from the gradient if the internal field is valid
    if (notNull(iF))
    {
        if (iF.size())
        {
            // Note: cannot ask for nf() if zero faces

            scalarField::operator=
            (
                //patchInternalField() + gradient()/patch().deltaCoeffs()
                // ***HGW Hack to avoid the construction of mesh.deltaCoeffs
                // which fails for AMI patches for some mapping operations
                patchInternalField()
              + gradient()*(patch().nf() & patch().delta())
            );
        }
    }
    else
    {
        // Enforce mapping of values so we have a valid starting value
        this->map(ptf, mapper);
    }
}


Foam::newFixedFluxPressureFvPatchScalarField::newFixedFluxPressureFvPatchScalarField
(
    const newFixedFluxPressureFvPatchScalarField& wbppsf
)
:
    fixedGradientFvPatchScalarField(wbppsf),
    curTimeIndex_(-1)
{}


Foam::newFixedFluxPressureFvPatchScalarField::newFixedFluxPressureFvPatchScalarField
(
    const newFixedFluxPressureFvPatchScalarField& wbppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(wbppsf, iF),
    curTimeIndex_(-1)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::newFixedFluxPressureFvPatchScalarField::updateSnGrad
(
    const scalarField& snGradp
)
{
    if (updated())
    {
        return;
    }

    curTimeIndex_ = this->db().time().timeIndex();

    gradient() = snGradp;
    fixedGradientFvPatchScalarField::updateCoeffs();
}


void Foam::newFixedFluxPressureFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    if (curTimeIndex_ != this->db().time().timeIndex())
    {
       /* FatalErrorInFunction
            << "updateCoeffs(const scalarField& snGradp) MUST be called before"
               " updateCoeffs() or evaluate() to set the boundary gradient."
            << exit(FatalError);*/
        Info << "##curTimeIndex_: " << curTimeIndex_ << endl;
        Info << "##current time index: " << this->db().time().timeIndex() << endl;
    }
}


void Foam::newFixedFluxPressureFvPatchScalarField::write(Ostream& os) const
{
    fixedGradientFvPatchScalarField::write(os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        newFixedFluxPressureFvPatchScalarField
    );
}


// ************************************************************************* //
