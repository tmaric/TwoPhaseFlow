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

#include "reconstructionError.H"
#include "addToRunTimeSelectionTable.H"

#include "reconstructionSchemes.H"
#include "implicitFunction.H"
#include "mathematicalConstants.H"
#include "davofInit.H"



// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void Foam::reconstructionError::writeFileHeader(Ostream& os) const
{
    // Add headers to output data
    writeHeader(os, "reconstructionError ");
    writeCommented(os, "Time");
    writeTabbed(os, "LNormalInf_");
    writeTabbed(os, "LNormal1_");
    writeTabbed(os, "LCentreInf_");
    writeTabbed(os, "LCentre1_");
    writeTabbed(os, "LCurvInf_");
    writeTabbed(os, "LCurv1_");
    writeTabbed(os, "k1");
    writeTabbed(os, "k2");
    writeTabbed(os, "length");
    writeTabbed(os, "nCells");
    writeTabbed(os, "avgVol");
    // appended (davof): existing consumers index the columns above by position
    writeTabbed(os, "LNormalDiff1");
    writeTabbed(os, "LNormalDiff2");
    writeTabbed(os, "LNormalDiffInf");
    writeTabbed(os, "areaSum");
    writeTabbed(os, "areaExact");
    writeTabbed(os, "nInterface");
    // appended (davof, 2026-10-12): the position error of centre()
    writeTabbed(os, "LCentreDiff1");
    writeTabbed(os, "LCentreDiff2");
    writeTabbed(os, "LCentreDiffInf");
    os  << endl;
}



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::reconstructionError::reconstructionError
(
    const fvMesh& mesh,
    const objectRegistry& obr,
    const dictionary& dict
)
:
    writeFile(obr, "reconstructionError", "error", dict),
    mesh_(mesh),
    LNormalInf_(0),
    LNormal1_(0),
    LCentreInf_(0),
    LCentre1_(0),
    LCurvInf_(0),
    LCurv1_(0),
    avgLNormal_(100),
    avgLCentre_(100),
    avgLCurv_(100),
    LNormalDiff1_(0),
    LNormalDiff2_(0),
    LNormalDiffInf_(0),
    areaSum_(0),
    areaExact_(0),
    nInterface_(0),
    LCentreDiff1_(0),
    LCentreDiff2_(0),
    LCentreDiffInf_(0)
{


    writeFileHeader(file());
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::reconstructionError::~reconstructionError()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::reconstructionError::calcError
(
    const dictionary& dict,
    const volVectorField& centre,
    const volVectorField& normal,
    const volScalarField& curv,
    bool randomMode
)
{

    if(randomMode)
    {
        avgLCentre_.clear();
        avgLNormal_.clear();
        avgLCurv_.clear();
        LNormalInf_ = 0;
        LNormal1_ = 0;
        LCurv1_ = 0;
        LCurvInf_ = 0;
    }

    Foam::autoPtr<Foam::implicitFunction> func
    (
        implicitFunction::New
        (
            word(dict.lookup("type")),
            dict
        )
    );

    word type (dict.lookup("type"));
    scalar exactCurv = 0;
    if(type == "sphere")
    {
        exactCurv = 2/(readScalar(dict.lookup("radius")));
    }
    else if(type == "cylinder")
    {
        exactCurv = 1/(readScalar(dict.lookup("radius")));
    }
    else if(type == "plane")
    {
        exactCurv = 0;
    }
    else
    {
        notImplemented("only sphere, cylinder and plane are supported")
    }

    // The position error (davof): the exact distance for the sphere, the
    // Newton-projected distance to the implicit surface otherwise.
    const bool isSphere = (type == "sphere");
    const point sphereOrigin =
        isSphere ? dict.get<point>("origin") : point(Zero);
    const scalar sphereRadius =
        isSphere ? dict.get<scalar>("radius") : scalar(0);
    auto positionError = [&](const point& x) -> scalar
    {
        if (isSphere)
        {
            return mag(sphereRadius - mag(x - sphereOrigin));
        }
        point foot;
        return davofInit::projectedDistance(func(), x, foot);
    };



    DynamicList<scalar> LNormal(1000);
    DynamicList<scalar> LCentre(1000);
    DynamicList<scalar> LCurv(1000);

    // appended (davof): |n - n_exact| moments, total area, interface count
    scalar sumD = 0, sumD2 = 0, maxD = 0, areaLocal = 0;
    scalar sumP = 0, sumP2 = 0, maxP = 0;
    label nLocal = 0;

    forAll(normal,cellI)
    {
        if(mag(normal[cellI]) != 0)
        {
           vector centreCellI = centre[cellI];
           vector normalCellI = normal[cellI];
           normalCellI /= mag(normalCellI);

           vector exactNormal = func->grad(centreCellI);
           exactNormal /= mag(exactNormal);

           LNormal.append(1-(normalCellI & exactNormal));
           LCentre.append(mag(func->distanceToSurfaces(centreCellI)));
           LCurv.append(mag(curv[cellI]-exactCurv));

           const scalar diff = mag(normalCellI - exactNormal);
           sumD += diff;
           sumD2 += diff*diff;
           maxD = max(maxD, diff);
           areaLocal += mag(normal[cellI]);
           ++nLocal;

           const scalar dp = positionError(centreCellI);
           sumP += dp;
           sumP2 += dp*dp;
           maxP = max(maxP, dp);
        }
    }

    nInterface_ = returnReduce(nLocal, sumOp<label>());
    const scalar sumDg = returnReduce(sumD, sumOp<scalar>());
    const scalar sumD2g = returnReduce(sumD2, sumOp<scalar>());
    LNormalDiff1_ = sumDg/max(nInterface_, label(1));
    LNormalDiff2_ = Foam::sqrt(sumD2g/max(nInterface_, label(1)));
    LNormalDiffInf_ = returnReduce(maxD, maxOp<scalar>());
    areaSum_ = returnReduce(areaLocal, sumOp<scalar>());
    const scalar sumPg = returnReduce(sumP, sumOp<scalar>());
    const scalar sumP2g = returnReduce(sumP2, sumOp<scalar>());
    LCentreDiff1_ = sumPg/max(nInterface_, label(1));
    LCentreDiff2_ = Foam::sqrt(sumP2g/max(nInterface_, label(1)));
    LCentreDiffInf_ = returnReduce(maxP, maxOp<scalar>());
    if (type == "sphere")
    {
        const scalar R = readScalar(dict.lookup("radius"));
        areaExact_ = 4.0*constant::mathematical::pi*R*R;
    }
    else
    {
        areaExact_ = 0;
    }


    // normal errors
    scalar Lmax = max(LNormal);
    scalar L1 = average(LNormal);

    // centre errors
    scalar LCmax = max(LCentre);
    scalar LC1 = average(LCentre);

    // curv errors
    scalar LCurvmax = max(LCurv); // /exactCurv;
    scalar LCurv1 = average(LCurv); // /exactCurv;

    avgLNormal_.append(L1);
    avgLCentre_.append(LC1);
    avgLCurv_.append(LCurv1);

    if(Lmax > LNormalInf_)
    {
        LNormalInf_ = Lmax;
    }

    if(LCmax > LCentreInf_)
    {
        LCentreInf_ = LCmax;
    }

    if(LCurvmax > LCurvInf_)
    {
        LCurvInf_ = LCurvmax;
    }

    LNormal1_ = average(avgLNormal_);
    LCentre1_ = average(avgLCentre_);
    LCurv1_ = average(avgLCurv_);


}


void Foam::reconstructionError::write(const scalar k1, const scalar k2)
{
    scalar length = average(mag(mesh_.delta())()).value();
    scalar nCells = mesh_.nCells();
    scalar avgVol = average(mesh_.V()).value();


    if (Pstream::master())
    {
        writeCurrentTime(file());

        file()
            << token::TAB << LNormalInf_
            << token::TAB << LNormal1_
            << token::TAB << LCentreInf_
            << token::TAB << LCentre1_
            << token::TAB << LCurvInf_
            << token::TAB << LCurv1_
            << token::TAB << k1
            << token::TAB << k2
            << token::TAB << length
            << token::TAB << nCells
            << token::TAB << avgVol
            << token::TAB << LNormalDiff1_
            << token::TAB << LNormalDiff2_
            << token::TAB << LNormalDiffInf_
            << token::TAB << areaSum_
            << token::TAB << areaExact_
            << token::TAB << nInterface_
            << token::TAB << LCentreDiff1_
            << token::TAB << LCentreDiff2_
            << token::TAB << LCentreDiffInf_
            << endl;
    }


}

void Foam::reconstructionError::write()
{
    scalar length = average(mag(mesh_.delta())()).value();
    scalar nCells = mesh_.nCells();
    scalar avgVol = average(mesh_.V()).value();


    if (Pstream::master())
    {
        writeCurrentTime(file());

        file()
            << token::TAB << LNormalInf_
            << token::TAB << LNormal1_
            << token::TAB << LCentreInf_
            << token::TAB << LCentre1_
            << token::TAB << LCurvInf_
            << token::TAB << LCurv1_
            << token::TAB << length
            << token::TAB << nCells
            << token::TAB << avgVol
            << token::TAB << LNormalDiff1_
            << token::TAB << LNormalDiff2_
            << token::TAB << LNormalDiffInf_
            << token::TAB << areaSum_
            << token::TAB << areaExact_
            << token::TAB << nInterface_
            << token::TAB << LCentreDiff1_
            << token::TAB << LCentreDiff2_
            << token::TAB << LCentreDiffInf_
            << endl;
    }


}


// ************************************************************************* //
