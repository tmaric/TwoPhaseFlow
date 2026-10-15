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

#include "exactNormal.H"
#include "davof.H"
#include "implicitFunction.H"
#include "addToRunTimeSelectionTable.H"
#include "profiling.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
inline namespace twoPhaseFlow
{
namespace davofNormalModels
{
    defineTypeNameAndDebug(exact, 0);
    addToRunTimeSelectionTable(davofNormalModel, exact, dictionary);
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::davofNormalModels::exact::exact
(
    const dictionary& dict,
    const reconstruction::davof& scheme
)
:
    davofNormalModel(typeName, dict, scheme, "all"),
    exactNormalDict_(dict.subDict("exactNormal")),
    functionType_(exactNormalDict_.get<word>("functionType")),
    origin0_(exactNormalDict_.get<vector>("origin")),
    U0_(exactNormalDict_.get<vector>("U0"))
{
    Info<< "davof: exact normal (" << functionType_ << ", origin "
        << origin0_ << ", U0 " << U0_ << ")" << endl;
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::davofNormalModels::exact::modelNormals
(
    const scalar tState,
    vectorField& n
)
{
    addProfilingInFunction(geometricVoF);
    dictionary funcDict(exactNormalDict_);
    funcDict.set("type", functionType_);
    funcDict.set("origin", origin0_ + U0_*tState);
    autoPtr<implicitFunction> func
    (
        implicitFunction::New(functionType_, funcDict)
    );

    const DynamicField<label>& labels = scheme_.interfaceLabels();
    const volVectorField& normal = scheme_.normal();
    const volVectorField& centre = scheme_.centre();
    const vectorField& C = mesh_.cellCentres();

    forAll(labels, i)
    {
        const label celli = labels[i];

        // The exact normal at the centroid of the plane positioned by
        // alpha along m_c, or at the cell centre when there is none.
        point x = C[celli];
        if (mag(normal[celli]) != 0)
        {
            x = centre[celli];
        }
        n[i] = func().grad(x);
    }
}


// ************************************************************************* //
