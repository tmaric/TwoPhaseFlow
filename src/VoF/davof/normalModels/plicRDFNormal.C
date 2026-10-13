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

\*---------------------------------------------------------------------------*/

#include "plicRDFNormal.H"
#include "davof.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
inline namespace twoPhaseFlow
{
namespace davofNormalModels
{
    defineTypeNameAndDebug(plicRDF, 0);
    addToRunTimeSelectionTable(davofNormalModel, plicRDF, dictionary);
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::davofNormalModels::plicRDF::plicRDF
(
    const dictionary& dict,
    const reconstruction::davof& scheme
)
:
    davofNormalModel(typeName, dict, scheme, "check"),
    helperPtr_(nullptr)
{
    Info<< "davof: constructing the helper plicRDF reconstruction" << endl;
    helperPtr_.reset
    (
        new reconstruction::plicRDF
        (
            const_cast<volScalarField&>(scheme.alpha1()),
            scheme.phi(),
            scheme.U(),
            scheme.modelDict()
        )
    );
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::davofNormalModels::plicRDF::modelNormals
(
    const scalar,
    vectorField& n
)
{
    helperPtr_->reconstruct(true);
    const volVectorField& normal = helperPtr_->normal();
    const DynamicField<label>& labels = scheme_.interfaceLabels();

    forAll(labels, i)
    {
        n[i] = normal[labels[i]];
    }
}


// ************************************************************************* //
