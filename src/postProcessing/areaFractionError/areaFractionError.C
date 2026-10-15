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

#include "areaFractionError.H"
#include "addToRunTimeSelectionTable.H"

#include "advectionSchemes.H"
#include "reconstructionSchemes.H"
#include "implicitFunction.H"
#include "davofInit.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(areaFractionError, 0);
    addToRunTimeSelectionTable(functionObject, areaFractionError, dictionary);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::areaFractionError::writeFileHeader(Ostream& os) const
{
    writeHeader(os, "areaFractionError ");
    writeCommented(os, "Time");
    writeTabbed(os, "L1");
    writeTabbed(os, "L2");
    writeTabbed(os, "Linf");
    writeTabbed(os, "nFaces");
    writeTabbed(os, "nInterfaceCells");
    writeTabbed(os, "consistency");
    os  << endl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::areaFractionError::areaFractionError
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    writeFile(obr_, name, typeName, dict),
    dict_(dict),
    fieldName_(),
    functionType_(),
    origin0_(Zero),
    U0_(Zero),
    U0FromField_(true),
    davofMethod_("quadraticFaces"),
    wispTol_(1e-3),
    consistency_(true),
    reconstructionError_(true),
    reconstructionErrorEveryStep_(false),
    recErrPtr_(nullptr)
{
    read(dict);

    writeFileHeader(file());
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::areaFractionError::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);
    writeFile::read(dict);

    fieldName_ = dict.get<word>("field");
    functionType_ = dict.get<word>("functionType");
    origin0_ = dict.get<vector>("origin");
    U0FromField_ = !dict.found("U0");
    U0_ = dict.getOrDefault<vector>("U0", vector::zero);
    davofMethod_ = dict.getOrDefault<word>("davofMethod", "quadraticFaces");
    wispTol_ = dict.getOrDefault<scalar>("wispTol", 1e-3);
    consistency_ = dict.getOrDefault<Switch>("consistency", Switch(true));
    reconstructionError_ = dict.getOrDefault<bool>("reconstructionError", true);
    reconstructionErrorEveryStep_ =
        dict.getOrDefault<bool>("reconstructionErrorEveryStep", false);

    if (reconstructionError_)
    {
        recErrPtr_.reset(new Foam::reconstructionError(mesh_, obr_, dict));
    }

    return true;
}


bool Foam::functionObjects::areaFractionError::execute()
{
    // Nothing: the exact state, the error and the extra reconstructions
    // are built in write() only, so a write control that fires at the end
    // (onEnd) computes nothing before it
    return true;
}


bool Foam::functionObjects::areaFractionError::write()
{
    const Time& runTime = mesh_.time();
    const scalar t = runTime.value();

    vector U0 = U0_;
    if (U0FromField_)
    {
        const volVectorField& U = lookupObject<volVectorField>("U");
        const scalarField magU(mag(U.primitiveField()));
        if (U.size() && (gMax(magU) - gMin(magU)) <= SMALL)
        {
            U0 = U[0];
        }
        reduce(U0, maxMagSqrOp<vector>());
    }

    // The implicit function translated to time t
    dictionary funcDict(dict_);
    funcDict.set("type", functionType_);
    funcDict.set("origin", origin0_ + U0*t);

    autoPtr<implicitFunction> func
    (
        implicitFunction::New(functionType_, funcDict)
    );

    const volScalarField& alpha = lookupObject<volScalarField>(fieldName_);
    const word afName(IOobject::groupName("alphaf", alpha.group()));

    // The exact DAVOF state of the translated function
    volScalarField alphaEx
    (
        IOobject
        (
            "areaFractionError:alphaExact",
            runTime.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        mesh_,
        dimensionedScalar(dimless, Zero),
        calculatedFvPatchScalarField::typeName
    );
    surfaceScalarField alphafEx
    (
        IOobject
        (
            "areaFractionError:alphafExact",
            runTime.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        mesh_,
        dimensionedScalar(dimless, Zero),
        calculatedFvsPatchScalarField::typeName
    );

    davofInit init(mesh_, func(), funcDict);
    if (davofMethod_ == "exactSphere")
    {
        init.exactSphere(alphaEx, alphafEx);
    }
    else if (davofMethod_ == "quadraticFaces")
    {
        init.quadraticFaces(alphaEx, alphafEx);
    }
    else if (davofMethod_ == "detrixheAslam" || davofMethod_ == "linearInterpolant")
    {
        init.detrixheAslam(alphaEx, alphafEx);
    }
    else
    {
        FatalErrorInFunction
            << "Unknown davofMethod '" << davofMethod_
            << "'. Valid: exactSphere, quadraticFaces, detrixheAslam."
            << exit(FatalError);
    }

    scalar L1 = -1, L2 = -1, Linf = -1, consistency = -1;
    label nFaces = 0, nInterface = 0;

    advectionSchemes* advPtr = nullptr;
    if (mesh_.foundObject<advectionSchemes>("advectionSchemes"))
    {
        advPtr = &mesh_.lookupObjectRef<advectionSchemes>("advectionSchemes");
    }

    if (foundObject<surfaceScalarField>(afName))
    {
        const surfaceScalarField& af = lookupObject<surfaceScalarField>(afName);

        vectorField mCur, mEx;
        scalarField afCur, afEx;
        davofInit::areaNormal(mesh_, af, mCur, afCur);
        davofInit::areaNormal(mesh_, alphafEx, mEx, afEx);

        const labelUList& own = mesh_.faceOwner();
        const labelUList& nei = mesh_.faceNeighbour();
        const label nInt = mesh_.nInternalFaces();

        boolList isInterface(mesh_.nCells(), false);
        label nIntLocal = 0;
        forAll(isInterface, celli)
        {
            const scalar thr = davofInit::wispThreshold(mesh_, celli, wispTol_);
            isInterface[celli] = (mag(mCur[celli]) > thr || mag(mEx[celli]) > thr);
            if (isInterface[celli]) ++nIntLocal;
        }

        boolList hasValue(mesh_.nFaces(), false);
        for (label facei = 0; facei < nInt; ++facei)
        {
            hasValue[facei] = true;
        }
        const surfaceScalarField::Boundary& bf = af.boundaryField();
        forAll(bf, patchi)
        {
            const polyPatch& pp = mesh_.boundaryMesh()[patchi];
            if (bf[patchi].size() == pp.size())
            {
                forAll(pp, i)
                {
                    hasValue[pp.start() + i] = true;
                }
            }
        }

        scalar sum = 0, sum2 = 0, worst = 0;
        label nLocal = 0;
        forAll(hasValue, facei)
        {
            if (!hasValue[facei]) continue;
            const bool adjacent =
                isInterface[own[facei]]
             || (facei < nInt && isInterface[nei[facei]]);
            if (!adjacent) continue;
            const scalar e = mag(afCur[facei] - afEx[facei]);
            sum += e;
            sum2 += e*e;
            worst = max(worst, e);
            ++nLocal;
        }

        nFaces = returnReduce(nLocal, sumOp<label>());
        nInterface = returnReduce(nIntLocal, sumOp<label>());
        sum = returnReduce(sum, sumOp<scalar>());
        sum2 = returnReduce(sum2, sumOp<scalar>());
        L1 = sum/max(nFaces, label(1));
        L2 = Foam::sqrt(sum2/max(nFaces, label(1)));
        Linf = returnReduce(worst, maxOp<scalar>());

        // The consistency of the transported state with the plane the
        // davof scheme positions in it (one extra reconstruction; -1 when
        // switched off)
        if (consistency_ && advPtr && advPtr->surf().type() == "davof")
        {
            reconstructionSchemes& surf = advPtr->surf();
            surf.reconstruct(true);
            const volVectorField& n = surf.normal();
            const scalarField& V = mesh_.V();
            scalar worstC = 0;
            for (const label celli : surf.interfaceLabels())
            {
                worstC = max
                (
                    worstC,
                    mag(mCur[celli] - n[celli])/Foam::pow(V[celli], 2.0/3.0)
                );
            }
            consistency = returnReduce(worstC, maxOp<scalar>());
        }
    }

    if (Pstream::master())
    {
        writeCurrentTime(file());
        file()
            << token::TAB << L1
            << token::TAB << L2
            << token::TAB << Linf
            << token::TAB << nFaces
            << token::TAB << nInterface
            << token::TAB << consistency
            << endl;
    }

    Log << type() << " " << name() << ": alphaf error L1 = " << L1
        << ", L2 = " << L2 << ", Linf = " << Linf << " over " << nFaces
        << " faces of " << nInterface << " interface cells, consistency "
        << consistency << endl;

    // The reconstruction error against the translated function
    if (recErrPtr_ && advPtr)
    {
        const bool last =
            (t + 0.5*runTime.deltaTValue() >= runTime.endTime().value());
        if (reconstructionErrorEveryStep_ || last)
        {
            reconstructionSchemes& surf = advPtr->surf();
            surf.reconstruct(true);
            volScalarField curv
            (
                IOobject
                (
                    "areaFractionError:curv",
                    runTime.timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    false
                ),
                mesh_,
                dimensionedScalar(dimless/dimLength, Zero)
            );
            recErrPtr_->calcError(funcDict, surf.centre(), surf.normal(), curv, true);
            recErrPtr_->write();
        }
    }

    return true;
}


// ************************************************************************* //
