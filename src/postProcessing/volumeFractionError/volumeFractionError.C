/*---------------------------------------------------------------------------*\
            Copyright (c) 2017-2019, German Aerospace Center (DLR)
-------------------------------------------------------------------------------
License
    This file is part of the VoFLibrary source code library, which is an
    unofficial extension to OpenFOAM.

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

#include "volumeFractionError.H"
#include "addToRunTimeSelectionTable.H"

#include "advectionSchemes.H"
#include "cutCellImpFunc.H"
#include "cutCellIso.H"
#include "profiling.H"
#include "davofOrbit.H"
#include "davofInit.H"
#include "surfaceFields.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(volumeFractionError, 0);
    addToRunTimeSelectionTable(functionObject, volumeFractionError, dictionary);
}
}


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void Foam::functionObjects::volumeFractionError::writeFileHeader(Ostream& os) const
{
    // Add headers to output data (the exact state named when it is the
    // exact sphere's, 2026-11-02; the former header otherwise)
    writeHeader
    (
        os,
        (exactState_ == "exactSphere")
      ? "volumeFractionError (exactState exactSphere)"
      : "volumeFractionError "
    );
    writeCommented(os, "Time");
    writeTabbed(os, "Eshape");
    writeTabbed(os, "Emass");
    writeTabbed(os, "Ebound");
    writeTabbed(os, "tRec");
    writeTabbed(os, "Tadv");
    writeTabbed(os, "celltoCellDist");
    writeTabbed(os, "Ncells");
    writeTabbed(os, "avg(CellVol)");
    os  << endl;
}



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::volumeFractionError::volumeFractionError
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    writeFile(obr_, name, typeName, dict),
    Field_(),
    dict_(dict),
    initMass_(0),
    initCentre_(dict.get<vector>("origin")),
    exactState_("cutCellIso")
{



    read(dict);

    writeFileHeader(file());
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::volumeFractionError::~volumeFractionError()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::volumeFractionError::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);
    writeFile::read(dict);

    Field_ = dict.get<word>("field");

    // The exact state (key exactState, 2026-11-02; the header file): auto
    // takes the exact sphere where davofInit::exactSphere applies (a sphere
    // on axis-aligned boxes), the former cutCellIso construction otherwise
    {
        const word fType(dict.get<word>("functionType"));
        const word key(dict.getOrDefault<word>("exactState", "auto"));
        if (key != "auto" && key != "exactSphere" && key != "cutCellIso")
        {
            FatalIOErrorInFunction(dict)
                << "exactState '" << key << "' unknown. Valid: auto,"
                << " exactSphere, cutCellIso." << exit(FatalIOError);
        }
        bool boxes = true;
        {
            const vectorField& Sf = mesh_.faceAreas();
            forAll(Sf, facei)
            {
                const scalar A = mag(Sf[facei]);
                if (A > VSMALL && cmptMax(cmptMag(Sf[facei]/A)) < 1 - 1e-10)
                {
                    boxes = false;
                    break;
                }
            }
            if (boxes)
            {
                const pointField& P = mesh_.points();
                const labelListList& cellPoints = mesh_.cellPoints();
                const scalarField& V = mesh_.V();
                forAll(V, celli)
                {
                    const boundBox bb(pointField(P, cellPoints[celli]), false);
                    const vector s = bb.span();
                    if (mag(s.x()*s.y()*s.z() - V[celli]) > 1e-8*V[celli])
                    {
                        boxes = false;
                        break;
                    }
                }
            }
            reduce(boxes, andOp<bool>());
        }
        if (key == "auto")
        {
            exactState_ = (fType == "sphere" && boxes) ? "exactSphere" : "cutCellIso";
        }
        else
        {
            exactState_ = key;
        }
        if (exactState_ == "exactSphere" && (fType != "sphere" || !boxes))
        {
            FatalIOErrorInFunction(dict)
                << "exactState exactSphere needs functionType sphere on"
                << " axis-aligned hexahedra" << exit(FatalIOError);
        }
        Info<< type() << " " << name() << ": exactState " << exactState_
            << " (key " << key << ")" << endl;
    }


    // set initial mass
    Foam::autoPtr<Foam::implicitFunction> func
    (
        implicitFunction::New
        (
            dict.get<word>("functionType"),
            dict
        )
    );

    scalarField f(mesh_.nPoints(),0.0);


    forAll(f,pI)
    {
        f[pI] =  func->value(mesh_.points()[pI]);
    }


    cutCellIso cutCell(mesh_,f);

    volScalarField alphaExact
    (
          IOobject
          (
              "alphaExact",
              mesh_.time().timeName(),
              mesh_,
              IOobject::NO_READ,
              IOobject::NO_WRITE
          ),
          mesh_,
          dimensionedScalar("0", dimless, 0),
          "calculated"
    );


    forAll(alphaExact,cellI)
    {
        label cellStatus = cutCell.calcSubCell(cellI,0.0);

        if(cellStatus == -1)
        {
            alphaExact[cellI] = 1;
        }
        else if(cellStatus == 1)
        {
            alphaExact[cellI] = 0;
        }
        else if(cellStatus == 0)
        {
            alphaExact[cellI]= cutCell.VolumeOfFluid();
        }

    }
    const volScalarField& alpha = lookupObject<volScalarField>(Field_);

    initMass_ = fvc::domainIntegrate(alpha.internalField()).value();

    return true;
}


bool Foam::functionObjects::volumeFractionError::execute()
{
    // do nothing

    return true;
}


bool Foam::functionObjects::volumeFractionError::write()
{
    addProfilingInFunction(geometricVoF);

    const volScalarField& alpha = lookupObject<volScalarField>(Field_);

    // the rotation about an axis (keys omega, axisPoint; davofOrbit.H,
    // 2026-11-01), zero by default: the exact body at its position at t
    vector omega; point axisPoint;
    davof::readOrbit(dict_, initCentre_, omega, axisPoint);

    // The translation velocity U0 of the exact body (2026-11-01):
    // - the key U0 when the dictionary gives it (as in areaFractionError);
    // - without U0 and with a rotation (omega nonzero): zero. The orbit
    //   already moves the body; a uniform U under a rotation is the circular
    //   translation without spin (setFlow mode function), whose exact body is
    //   the orbit itself, and U t added to it displaced the sphere (Eshape of
    //   the orbit runs read the whole sphere volume, 1.414e-2);
    // - without U0 and without a rotation: as before, U of the first cell
    //   when |U| is uniform over the cells and the patches (within SMALL),
    //   else zero. This fallback is valid for a constant uniform U only
    //   (origin + U0 t); the translation templates that give no U0 rely on it
    //   and are unchanged.
    vector U0(Zero);
    if (dict_.found("U0"))
    {
        U0 = dict_.get<vector>("U0");
    }
    else if (mag(omega) <= VSMALL)
    {
        const volVectorField& U = mesh_.lookupObject<volVectorField>("U");

        volScalarField magU ("magU",mag(U));

        scalar maxMagU = max(magU).value();
        scalar minMagU = min(magU).value();

        reduce(maxMagU,maxOp<scalar>());
        reduce(minMagU,minOp<scalar>());

        if (U.size() && (maxMagU - minMagU) <= SMALL)
        {
            U0 = U[0];
        }
        reduce(U0, maxMagSqrOp<vector>());
    }

    vector centrePos = davof::orbitOrigin(initCentre_, U0, omega, axisPoint, mesh_.time().value());






    dict_.set<vector>("origin",centrePos);

    Foam::autoPtr<Foam::implicitFunction> func =  implicitFunction::New
    (
           word(dict_.lookup("functionType")),
           dict_
    );



    volScalarField alphaExact
    (
          IOobject
          (
              "alphaExact",
              mesh_.time().timeName(),
              mesh_,
              IOobject::NO_READ,
              IOobject::NO_WRITE
          ),
          mesh_,
          dimensionedScalar("0", dimless, 0),
          "calculated"
    );

    if (exactState_ == "exactSphere")
    {
        // the exact sphere at origin(t) (2026-11-02): the cell fractions of
        // davofInit::exactSphere, the state of initDavofFields and the
        // reference of areaFractionError
        surfaceScalarField alphafExact
        (
            IOobject
            (
                "volumeFractionError:alphafExact",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensionedScalar(dimless, Zero),
            calculatedFvsPatchScalarField::typeName
        );
        davofInit init(mesh_, func(), dict_);
        init.exactSphere(alphaExact, alphafExact);
    }
    else
    {
        // the former construction: the cell cut by the isosurface of the
        // vertex values
        scalarField f(mesh_.nPoints(),0.0);

        forAll(f,pI)
        {
            f[pI] =  func->value(mesh_.points()[pI]);
        }

        cutCellIso cutCell(mesh_,f); // changes f

        forAll(alphaExact,cellI)
        {
            label cellStatus = cutCell.calcSubCell(cellI,0.0);

            if(cellStatus == -1)
            {
                alphaExact[cellI] = 1;
            }
            else if(cellStatus == 1)
            {
                alphaExact[cellI] = 0;
            }
            else if(cellStatus == 0)
            {
                alphaExact[cellI]= cutCell.VolumeOfFluid();
            }
        }
    }


    scalar E1 = fvc::domainIntegrate(mag(alpha.internalField() - alphaExact.internalField())).value();
    scalar Emass = mag(fvc::domainIntegrate(alpha.internalField()).value() - initMass_);

    scalarField fluidVol = alpha.primitiveField();
    scalar Ebound = mag(max(-min(fluidVol),max(fluidVol-1)));

    reduce(Ebound,maxOp<scalar>());





    const advectionSchemes& advection = mesh_.lookupObjectRef<advectionSchemes>("advectionSchemes");

    scalar reconAvgTimePerIter = advection.runTime().x() / mesh_.time().timeIndex();
    scalar advectionAvgTimePerIter = advection.runTime().y() / mesh_.time().timeIndex();

    scalar length = average(mag(mesh_.delta())()).value();
    scalar nCells = mesh_.nCells();
    scalar avgVol = average(mesh_.V()).value();


//    Log << type() << " " << volumeFractionError.name() << " is " << intValue << endl;


    if (Pstream::master())
    {
        writeCurrentTime(file());

        file()
            << token::TAB << E1
            << token::TAB << Emass
            << token::TAB << Ebound
            << token::TAB << reconAvgTimePerIter
            << token::TAB << advectionAvgTimePerIter
            << token::TAB << length
            << token::TAB << nCells
            << token::TAB << avgVol
            << endl;
    }



    return true;
}


// ************************************************************************* //
