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

#include "functionObjectTemplate.H"
#define namespaceFoam  // Suppress <using namespace Foam;>
#include "fvCFD.H"
#include "unitConversion.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(visRDFFunctionObject, 0);

addRemovableToRunTimeSelectionTable
(
    functionObject,
    visRDFFunctionObject,
    dictionary
);


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

// dynamicCode:
// SHA1 = dab14549aa9cd73d878d5d81815aad306b5bea11
//
// unique function name that can be checked if the correct library version
// has been loaded
extern "C" void visRDF_dab14549aa9cd73d878d5d81815aad306b5bea11(bool load)
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


// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

//{{{ begin localCode

//}}} end localCode

} // End namespace Foam


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

const Foam::fvMesh&
Foam::visRDFFunctionObject::mesh() const
{
    return refCast<const fvMesh>(obr_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::
visRDFFunctionObject::
visRDFFunctionObject
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    functionObjects::regionFunctionObject(name, runTime, dict)
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::
visRDFFunctionObject::
~visRDFFunctionObject()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool
Foam::
visRDFFunctionObject::read(const dictionary& dict)
{
    if (false)
    {
        printMessage("read visRDF");
    }

//{{{ begin code
    
//}}} end code

    return true;
}


bool
Foam::
visRDFFunctionObject::execute()
{
    if (false)
    {
        printMessage("execute visRDF");
    }

//{{{ begin code
    
//}}} end code

    return true;
}


bool
Foam::
visRDFFunctionObject::write()
{
    if (false)
    {
        printMessage("write visRDF");
    }

//{{{ begin code
    #line 37 "/work/groups/da_mma_b/jun/TwoPhaseFlow/run/rhoVoF/jetCrossFlow/system/controlDict.functions.visRDF"
Info<< "###Reading/reassigning RDF field \n" << endl;

            volScalarField alpha1 = mesh().lookupObject<volScalarField>("alpha.water");
            const volVectorField& U = mesh().lookupObject<volVectorField>("U");
            const surfaceScalarField& phi = mesh().lookupObject<surfaceScalarField>("phi");
            Foam::advection::isoAdvection advector(alpha1, phi, U);
            word reconScheme(advector.surf().modelDict().lookup("reconstructionScheme"));
            scalar surfCellTol = advector.surf().modelDict().lookupOrDefault("surfCellTol", 1e-8);
            const word plicRDF("plicRDF");
            if(reconScheme == plicRDF)
            {
                const volScalarField& RDF = mesh().lookupObject<volScalarField>("reconstructedDistanceFunction");
                volScalarField visRDF
                (
                    IOobject
                    (
                        "visRDF",
                        mesh().time().timeName(),
                        mesh(),
                        IOobject::NO_READ,
                        IOobject::AUTO_WRITE 
                     ),
                    RDF
                );
                const volScalarField& cellDistLevel = mesh().lookupObject<volScalarField>("cellDistLevel_");

                forAll(cellDistLevel, celli)
                {
                    if (cellDistLevel[celli] == -1)
                    {
                        if(alpha1[celli] > 1 - surfCellTol)
                        {
                            visRDF[celli] = 100;
                        } 
                        else if (alpha1[celli] < surfCellTol) 
                        {
                            visRDF[celli] = -100;
                        }
                    }
                }
    
                forAll(visRDF.boundaryField(), patchI)
                {
                    fvPatchScalarField& pRDF = visRDF.boundaryFieldRef()[patchI];
                    if (isA<calculatedFvPatchScalarField>(pRDF))
                    {
                        const polyPatch& pp = pRDF.patch().patch();
                        forAll(pRDF, i)
                        {
                            const label pCellI = pp.faceCells()[i];
                            if(cellDistLevel[pCellI] == -1)
                            {
                                pRDF[i] = visRDF[pCellI];
                            }
                        }
                    }
                }
                visRDF.correctBoundaryConditions();
                visRDF.write();
            }
//}}} end code

    return true;
}


bool
Foam::
visRDFFunctionObject::end()
{
    if (false)
    {
        printMessage("end visRDF");
    }

//{{{ begin code
    
//}}} end code

    return true;
}


// ************************************************************************* //

