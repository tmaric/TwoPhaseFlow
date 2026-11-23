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

#include "ringSumRDFNormal.H"
#include "davof.H"
#include "zoneDistribute.H"
#include "leastSquareGrad.H"
#include "addToRunTimeSelectionTable.H"
#include "profiling.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
inline namespace twoPhaseFlow
{
namespace davofNormalModels
{
    defineTypeNameAndDebug(ringSumRDF, 0);
    addToRunTimeSelectionTable(davofNormalModel, ringSumRDF, dictionary);
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::davofNormalModels::ringSumRDF::ringSumRDF
(
    const dictionary& dict,
    const reconstruction::davof& scheme
)
:
    davofNormalModel(typeName, dict, scheme, "all"),
    rdfIterations_(dict.getOrDefault<label>("rdfIterations", 1)),
    RDF_(reconstructedDistanceFunction::New(scheme.mesh()))
{
    if (rdfIterations_ < 0)
    {
        FatalIOErrorInFunction(dict)
            << "normalModel ringSumRDF: rdfIterations " << rdfIterations_
            << " < 0." << exit(FatalIOError);
    }
    if (scheme.interfaceCellRule() != "alphaAndRingSum")
    {
        WarningInFunction
            << "normalModel ringSumRDF with interfaceCellRule "
            << scheme.interfaceCellRule()
            << ": the ring sum is computed by the rule alphaAndRingSum only"
            << endl;
    }
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::davofNormalModels::ringSumRDF::positionPlanes(const vectorField& n)
{
    // davof is a friend of this class's access: positionPlane is private
    reconstruction::davof& scheme = const_cast<reconstruction::davof&>(scheme_);
    const DynamicField<label>& labels = scheme_.interfaceLabels();
    forAll(labels, i)
    {
        scheme.positionPlane(labels[i], n[i]);
    }
    const_cast<volVectorField&>(scheme_.normal()).correctBoundaryConditions();
    const_cast<volVectorField&>(scheme_.centre()).correctBoundaryConditions();
}


void Foam::davofNormalModels::ringSumRDF::gradSurf
(
    const volScalarField& phi,
    vectorField& g
)
{
    // The body of plicRDF::gradSurf(), over the davof scheme's interface cells
    leastSquareGrad<scalar> lsGrad("polyDegree1", mesh_.geometricD());
    zoneDistribute& exchangeFields = zoneDistribute::New(mesh_);
    const boolList& zone = scheme_.interfaceCell();
    exchangeFields.setUpCommforZone(zone, false);

    Map<vector> mapCC
    (
        exchangeFields.getDatafromOtherProc(zone, mesh_.C())
    );
    Map<scalar> mapPhi
    (
        exchangeFields.getDatafromOtherProc(zone, phi)
    );

    DynamicField<vector> cellCentre(100);
    DynamicField<scalar> phiValues(100);
    const labelListList& stencil = exchangeFields.getStencil();
    const DynamicField<label>& labels = scheme_.interfaceLabels();
    g.resize(labels.size());

    forAll(labels, i)
    {
        const label celli = labels[i];
        cellCentre.clear();
        phiValues.clear();
        for (const label gblIdx : stencil[celli])
        {
            cellCentre.append
            (
                exchangeFields.getValue(mesh_.C(), mapCC, gblIdx)
            );
            phiValues.append
            (
                exchangeFields.getValue(phi, mapPhi, gblIdx)
            );
        }
        cellCentre -= mesh_.C()[celli];
        g[i] = lsGrad.grad(cellCentre, phiValues);
    }
}


void Foam::davofNormalModels::ringSumRDF::modelNormals
(
    const scalar,
    vectorField& n
)
{
    addProfilingInFunction(geometricVoF);
    const vectorField& d = scheme_.ringSum().primitiveField();
    const DynamicField<label>& labels = scheme_.interfaceLabels();

    // Step 1: the ring sums
    n.resize(labels.size());
    forAll(labels, i)
    {
        n[i] = d[labels[i]];
    }
    if (rdfIterations_ == 0)
    {
        return;
    }

    // The planes positioned by alpha along the ring sums
    positionPlanes(n);

    // The stencil of the cells next to the interface cells, as
    // plicRDF::setInitNormals() sets it
    zoneDistribute& exchangeFields = zoneDistribute::New(mesh_);
    RDF_.markCellsNearSurf(scheme_.interfaceCell(), 1);
    const boolList& nextToInterface = RDF_.nextToInterface();
    exchangeFields.updateStencil(nextToInterface);

    vectorField g;
    for (label iter = 0; iter < rdfIterations_; ++iter)
    {
        RDF_.constructRDF
        (
            nextToInterface,
            scheme_.centre(),
            scheme_.normal(),
            exchangeFields,
            false
        );
        gradSurf(RDF_, g);
        forAll(g, i)
        {
            const scalar gMag = mag(g[i]);
            n[i] = (gMag > 0) ? g[i]/gMag : vector(Zero);
        }
        if (iter + 1 < rdfIterations_)
        {
            positionPlanes(n);
        }
    }
}


// ************************************************************************* //
