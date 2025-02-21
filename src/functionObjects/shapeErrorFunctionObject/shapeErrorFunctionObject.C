/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2024 AUTHOR,AFFILIATION
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

#include "shapeErrorFunctionObject.H"
#include "Time.H"
#include "fvMesh.H"
#include "addToRunTimeSelectionTable.H"

#include "implicitFunction.H"
#include "cutCellIso.H"
#include "cutFaceIso.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(shapeErrorFunctionObject, 0);
    addToRunTimeSelectionTable(functionObject, shapeErrorFunctionObject, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::shapeErrorFunctionObject::shapeErrorFunctionObject
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    setAlphaFieldDict_
    (
        IOobject
        (
            "setAlphaFieldDict", 
            mesh_.time().system(), 
            mesh_, 
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),
	alpha_(mesh_.lookupObject<volScalarField>( setAlphaFieldDict_.get<word>("field") )),
    alphaExa_ 
    (
        IOobject
        (
            "alphaExa_", 
            mesh_.time().timeName(), 
            mesh_, 
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("alphaExa_", dimless, pTraits<scalar>::zero)
    ),
	vel_(Foam::vector(0,0,0)),
    origin_(setAlphaFieldDict_.get<vector>("origin"))//(Foam::vector(0,0,0))
{
    read(dict);
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::functionObjects::shapeErrorFunctionObject::setAlpha
(
    volScalarField& alpha1,
    scalarField& f
)
{
    const fvMesh& mesh = alpha1.mesh();
    cutCellIso cutCell(mesh, f);
    cutFaceIso cutFace(mesh, f);

    forAll(alpha1, cellI)
    {
        cutCell.calcSubCell(cellI, 0.0);

        alpha1[cellI] = max(min(cutCell.VolumeOfFluid(), 1), 0);
    }

    // Setting boundary alpha1 values
    forAll(mesh.boundary(), patchi)
    {
        if (mesh.boundary()[patchi].size() > 0)
        {
            const label start = mesh.boundary()[patchi].patch().start();
            scalarField& alphap = alpha1.boundaryFieldRef()[patchi];
            const scalarField& magSfp = mesh.magSf().boundaryField()[patchi];

            forAll(alphap, patchFacei)
            {
             	const label facei = patchFacei + start;
                cutFace.calcSubFace(facei, 0.0);
                alphap[patchFacei] =
                    mag(cutFace.subFaceArea())/magSfp[patchFacei];
            }
		}
    }
}


bool Foam::functionObjects::shapeErrorFunctionObject::read(const dictionary& dict)
{
	vel_ = dict.get<vector>("velocity");
    //origin_ = setAlphaFieldDict_.get<vector>("origin");
    return true;
}

bool Foam::functionObjects::shapeErrorFunctionObject::execute()
{

	//vector origin = setAlphaFieldDict_.get<vector>("origin");
	vector origin = origin_ + vel_*mesh_.time().value();
	setAlphaFieldDict_.add("origin", origin, true);
    autoPtr<implicitFunction> func = implicitFunction::New
     (
         setAlphaFieldDict_.get<word>("type"),
         setAlphaFieldDict_
     );

     scalarField f(mesh_.nPoints(), Zero);

     forAll(f, pi)
     {
         f[pi] = func->value(mesh_.points()[pi]);
     };

	setAlpha(alphaExa_, f);

	volScalarField alphaDiff = mag(alphaExa_ - alpha_);

	scalar shapeErr = gSum((alphaDiff*mesh_.V()).ref())/gSum((alphaExa_*mesh_.V()).ref());
	//alphaExa_.write();
	Info << "### shape err: " << shapeErr << endl;
    return true;
}


bool Foam::functionObjects::shapeErrorFunctionObject::end()
{
    return true;
}


bool Foam::functionObjects::shapeErrorFunctionObject::write()
{
	alphaExa_.write();
    return true;
}


// ************************************************************************* //
