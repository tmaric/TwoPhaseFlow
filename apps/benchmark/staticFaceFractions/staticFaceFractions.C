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

Application
    staticFaceFractions

Description
    The static face fractions of a transported volume-fraction field
    (2026-11-02): at every selected time the face fractions reconstructed
    from alpha alone, scored with the L1(alphaf) of the areaFractionError
    function object against the exact DAVOF state of the sphere at its
    position at that time, next to the transported face fractions alphaf
    when the time directory holds them, and the cell-level error of alpha.
    Serial; a post-processing step that writes no field.

    The exact body and the scoring are those of the areaFractionError
    entry of system/controlDict's functions (option -entry, default
    areaFractionError): functionType sphere, origin, radius, U0, omega,
    axisPoint (fork davofOrbit.H), davofMethod, wispTol. The counted faces
    are those of areaFractionError: every face with a value adjacent to a
    cell whose area vector m_c = sum_f alpha_f S_f^out (of the scored state
    or of the exact one) exceeds wispTol V_c^(2/3).

    The static state from a field a (alpha, or the exact alpha_exact):
    - interface cells: surfCellTol < a < 1 - surfCellTol (isoAdvector's
      surface cells; option -surfCellTol, default 1e-8);
    - the plane of an interface cell: the exact normal n = (o - x)/|o - x|
      (into the liquid ball of centre o at time t) positioned by a, the
      volume under it a V (surfaceIteratorPLIC::vofCutCell, tolerance
      -positionTol in alpha, default 1e-11, the loop's anchorTolRel), first
      with x the cell centre, then repositioned with x the centroid of the
      first plane (the exact normal model takes the normal at the centroid
      of the positioned plane);
    - a face of a cell: the area of the face on the liquid side of the
      cell's plane over |S_f| (cutFacePLIC); a cell without a plane gives
      1 for a >= 1/2 and 0 otherwise (full, empty, or a cell whose
      positioning failed, counted);
    - the face fraction, rule upwind: the value of the face's donor cell,
      the owner for phi_f >= 0 and the neighbour otherwise (boundary faces:
      the owner), the limit dt -> 0 of the loop's face update
      (sweptFaceCells pieces with the donorSide coplanar rule: the
      pre-image is the face itself, given to the donor); phi from the time
      directory, else U0.S_f; rule mean: the mean of the two cells' values.

    The cell-level error over the interface cells of the exact state
    (surfCellTol < alpha_exact < 1 - surfCellTol, nEx cells):
        Ec   = sum_{exact interface} |alpha - alpha_exact| / nEx,
        Eout = sum_{other cells} |alpha - alpha_exact| / nEx,
        Evol = sum_all |alpha - alpha_exact| V (volumeFractionError's
               Eshape numerator, for the check).

    The same static rule serves a comparator's alpha (plicRDF + isoAdvector:
    J. Roenby, H. Bredmose, H. Jasak, A computational method for sharp
    interface advection, R. Soc. Open Sci. 3 (2016) 160405,
    doi:10.1098/rsos.160405, arXiv:1601.05392; H. Scheufler, J. Roenby,
    Accurate and efficient surface reconstruction from volume fraction data
    on general meshes, J. Comput. Phys. 383 (2019) 1-23,
    doi:10.1016/j.jcp.2019.01.009, arXiv:1801.05382), whose time directories
    hold no alphaf.

    Written to postProcessing/staticFaceFractions/staticFaceFractions.dat,
    one row per selected time:
        Time L1loop nFacesLoop L1static nFacesStatic L1staticMean
        L1staticExact nPlanes nPlaneFail Ec Eout Evol nExInterface
    (L1loop -1 without alphaf in the time directory; L1staticExact the
    upwind static state of alpha_exact, the reconstruction floor).

        staticFaceFractions -time '0:1'      (the written times, 0 included)
        staticFaceFractions -time '0.1:0.45'

    The re-scoring of a shape error (2026-11-03): -entry volumeFractionError
    -davofMethod exactSphere -evolOnly -out <file> computes, at the selected
    times, Evol = sum |alpha - alpha_exact| V against the exact sphere (the
    Eshape of volumeFractionError with exactState exactSphere) and Ec, Eout,
    without the face fractions and without writing into the case.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "timeSelector.H"
#include "implicitFunction.H"
#include "davofInit.H"
#include "davofOrbit.H"
#include "surfaceIteratorPLIC.H"
#include "cutFacePLIC.H"
#include "OFstream.H"
#include "OSspecific.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace
{

using namespace Foam;

// The static face fractions of the field a (see the description): the
// upwind rule into afU, the mean rule into afM
void staticFractions
(
    const fvMesh& mesh,
    const scalarField& a,
    const point& o,
    const scalarField& phiF,
    const scalar surfCellTol,
    const scalar positionTol,
    surfaceScalarField& afU,
    surfaceScalarField& afM,
    label& nPlanes,
    label& nFail
)
{
    const label nCells = mesh.nCells();
    const label nInt = mesh.nInternalFaces();
    const vectorField& C = mesh.cellCentres();
    const vectorField& Sf = mesh.faceAreas();
    const labelUList& own = mesh.faceOwner();
    const labelUList& nei = mesh.faceNeighbour();

    surfaceIteratorPLIC sIter(mesh, surfCellTol);
    cutFacePLIC cutF(mesh);

    vectorField nC(nCells, Zero);
    pointField bC(nCells, Zero);
    bitSet hasPlane(nCells);
    nPlanes = 0;
    nFail = 0;

    forAll(a, celli)
    {
        if (!(a[celli] > surfCellTol && a[celli] < 1 - surfCellTol))
        {
            continue;
        }
        bool ok = false;
        const vector n1 = o - C[celli];
        if (mag(n1) > 0)
        {
            const label s1 = sIter.vofCutCell(celli, a[celli], positionTol, 100, n1);
            if (s1 == 0 && mag(sIter.surfaceArea()) > 0)
            {
                const vector n2 = o - sIter.surfaceCentre();
                if (mag(n2) > 0)
                {
                    const label s2 =
                        sIter.vofCutCell(celli, a[celli], positionTol, 100, n2);
                    if (s2 == 0 && mag(sIter.surfaceArea()) > 0)
                    {
                        nC[celli] = n2/mag(n2);
                        bC[celli] = sIter.surfaceCentre();
                        hasPlane.set(celli);
                        ok = true;
                    }
                }
            }
        }
        if (ok) ++nPlanes; else ++nFail;
    }

    auto cellValue = [&](const label celli, const label facei) -> scalar
    {
        if (hasPlane.test(celli))
        {
            cutF.calcSubFace(facei, nC[celli], bC[celli]);
            return min(max(mag(cutF.subFaceArea())/mag(Sf[facei]), scalar(0)), scalar(1));
        }
        return (a[celli] >= 0.5) ? scalar(1) : scalar(0);
    };

    scalarField& u = afU.primitiveFieldRef();
    scalarField& m = afM.primitiveFieldRef();
    label nBoth = 0, nDiffer = 0;
    scalar maxDiff = 0;
    for (label facei = 0; facei < nInt; ++facei)
    {
        const scalar vo = cellValue(own[facei], facei);
        const scalar vn = cellValue(nei[facei], facei);
        u[facei] = (phiF[facei] >= 0) ? vo : vn;
        m[facei] = 0.5*(vo + vn);
        if (hasPlane.test(own[facei]) && hasPlane.test(nei[facei]))
        {
            ++nBoth;
            if (vo != vn) ++nDiffer;
            maxDiff = max(maxDiff, mag(vo - vn));
        }
    }
    Info<< "    staticFractions: " << nBoth << " faces between two planes, "
        << nDiffer << " with different clips, max |difference| " << maxDiff
        << endl;
    const polyBoundaryMesh& pbm = mesh.boundaryMesh();
    forAll(pbm, patchi)
    {
        fvsPatchScalarField& pu = afU.boundaryFieldRef()[patchi];
        fvsPatchScalarField& pm = afM.boundaryFieldRef()[patchi];
        if (pu.size() != pbm[patchi].size()) continue;
        forAll(pu, i)
        {
            const label facei = pbm[patchi].start() + i;
            pu[i] = cellValue(own[facei], facei);
            pm[i] = pu[i];
        }
    }
}


// The L1(alphaf) of areaFractionError (serial): the face fractions af
// against the exact ones over the faces adjacent to an interface cell of
// either state
scalar scoreL1
(
    const fvMesh& mesh,
    const surfaceScalarField& af,
    const scalarField& afEx,
    const vectorField& mEx,
    const scalar wispTol,
    label& nFaces
)
{
    vectorField mCur;
    scalarField afCur;
    davofInit::areaNormal(mesh, af, mCur, afCur);

    const labelUList& own = mesh.faceOwner();
    const labelUList& nei = mesh.faceNeighbour();
    const label nInt = mesh.nInternalFaces();

    boolList isInterface(mesh.nCells(), false);
    forAll(isInterface, celli)
    {
        const scalar thr = davofInit::wispThreshold(mesh, celli, wispTol);
        isInterface[celli] = (mag(mCur[celli]) > thr || mag(mEx[celli]) > thr);
    }

    boolList hasValue(mesh.nFaces(), false);
    for (label facei = 0; facei < nInt; ++facei)
    {
        hasValue[facei] = true;
    }
    forAll(af.boundaryField(), patchi)
    {
        const polyPatch& pp = mesh.boundaryMesh()[patchi];
        if (af.boundaryField()[patchi].size() == pp.size())
        {
            forAll(pp, i)
            {
                hasValue[pp.start() + i] = true;
            }
        }
    }

    scalar sum = 0;
    nFaces = 0;
    forAll(hasValue, facei)
    {
        if (!hasValue[facei]) continue;
        const bool adjacent =
            isInterface[own[facei]] || (facei < nInt && isInterface[nei[facei]]);
        if (!adjacent) continue;
        sum += mag(afCur[facei] - afEx[facei]);
        ++nFaces;
    }
    return sum/max(nFaces, label(1));
}

} // End anonymous namespace


int main(int argc, char *argv[])
{
    argList::addNote
    (
        "The static face fractions of the transported alpha (exact-normal"
        " planes positioned by alpha), scored with areaFractionError's"
        " L1(alphaf) next to the transported alphaf; serial"
    );
    timeSelector::addOptions();
    argList::addOption
    (
        "entry", "word",
        "the entry of system/controlDict's functions that holds the exact"
        " body and the scoring keys (default areaFractionError)"
    );
    argList::addOption("surfCellTol", "scalar", "interface cells (default 1e-8)");
    argList::addOption("positionTol", "scalar", "plane positioning (default 1e-11)");
    argList::addOption
    (
        "davofMethod", "word",
        "the exact state's construction, overriding the entry's davofMethod"
        " (an entry without the key, volumeFractionError's, defaults to"
        " quadraticFaces)"
    );
    argList::addBoolOption
    (
        "evolOnly",
        "the cell-level errors Ec, Eout, Evol only (no face fractions, no phi;"
        " the L1 columns and the plane counts are written as -1)"
    );
    argList::addOption
    (
        "out", "file",
        "the output file (default postProcessing/staticFaceFractions/"
        "staticFaceFractions.dat of the case)"
    );
    argList::noParallel();

    #include "setRootCase.H"
    #include "createTime.H"
    instantList timeDirs = timeSelector::select0(runTime, args);
    #include "createMesh.H"

    const word entry(args.getOrDefault<word>("entry", "areaFractionError"));
    const scalar surfCellTol = args.getOrDefault<scalar>("surfCellTol", 1e-8);
    const scalar positionTol = args.getOrDefault<scalar>("positionTol", 1e-11);

    const dictionary& foDict =
        runTime.controlDict().subDict("functions").subDict(entry);
    const word fieldName(foDict.get<word>("field"));
    const word functionType(foDict.get<word>("functionType"));
    if (functionType != "sphere")
    {
        FatalIOErrorInFunction(foDict)
            << "staticFaceFractions: the exact normal is the sphere's;"
            << " functionType " << functionType << " is not supported"
            << exit(FatalIOError);
    }
    const point origin0(foDict.get<point>("origin"));
    vector omega;
    point axisPoint;
    davof::readOrbit(foDict, origin0, omega, axisPoint);
    const bool U0FromField = !foDict.found("U0") && mag(omega) <= VSMALL;
    const word davofMethod
    (
        args.found("davofMethod")
      ? args.get<word>("davofMethod")
      : foDict.getOrDefault<word>("davofMethod", "quadraticFaces")
    );
    const bool evolOnly = args.found("evolOnly");
    const scalar wispTol = foDict.getOrDefault<scalar>("wispTol", 1e-3);

    fileName outFile
    (
        runTime.globalPath()/"postProcessing"/"staticFaceFractions"
       /"staticFaceFractions.dat"
    );
    if (args.found("out"))
    {
        outFile = args.get<fileName>("out");
    }
    mkDir(outFile.path());
    OFstream os(outFile);
    os.precision(12);
    os  << "# staticFaceFractions (entry " << entry << ", surfCellTol "
        << surfCellTol << ", positionTol " << positionTol << ", wispTol "
        << wispTol << ", davofMethod " << davofMethod << ")" << nl
        << "# Time\tL1loop\tnFacesLoop\tL1static\tnFacesStatic\tL1staticMean"
        << "\tL1staticExact\tnPlanes\tnPlaneFail\tEc\tEout\tEvol\tnExInterface"
        << nl;

    const vectorField& Sf = mesh.faceAreas();
    const scalarField& V = mesh.V();
    const label nInt = mesh.nInternalFaces();

    forAll(timeDirs, timei)
    {
        runTime.setTime(timeDirs[timei], timei);
        const scalar t = runTime.value();
        Info<< "Time = " << runTime.timeName() << endl;

        volScalarField alpha
        (
            IOobject
            (
                fieldName, runTime.timeName(), mesh,
                IOobject::MUST_READ, IOobject::NO_WRITE, false
            ),
            mesh
        );

        // the exact body's translation (as areaFractionError)
        vector U0 = foDict.getOrDefault<vector>("U0", vector::zero);
        if (U0FromField)
        {
            volVectorField U
            (
                IOobject
                (
                    "U", runTime.timeName(), mesh,
                    IOobject::MUST_READ, IOobject::NO_WRITE, false
                ),
                mesh
            );
            const scalarField magU(mag(U.primitiveField()));
            if (U.size() && (max(magU) - min(magU)) <= SMALL)
            {
                U0 = U[0];
            }
        }
        const point o = davof::orbitOrigin(origin0, U0, omega, axisPoint, t);

        // phi of the time directory, else U0.S_f
        scalarField phiF(evolOnly ? 0 : nInt);
        if (!evolOnly)
        {
            IOobject phiIO
            (
                "phi", runTime.timeName(), mesh,
                IOobject::MUST_READ, IOobject::NO_WRITE, false
            );
            if (phiIO.typeHeaderOk<surfaceScalarField>(true))
            {
                const surfaceScalarField phi(phiIO, mesh);
                phiF = phi.primitiveField();
            }
            else
            {
                Info<< "    no phi in " << runTime.timeName()
                    << ": the upwind direction from U0.S_f" << endl;
                forAll(phiF, facei)
                {
                    phiF[facei] = U0 & Sf[facei];
                }
            }
        }

        // the exact state (areaFractionError's construction)
        dictionary funcDict(foDict);
        funcDict.set("type", functionType);
        funcDict.set("origin", o);
        autoPtr<implicitFunction> func(implicitFunction::New(functionType, funcDict));
        volScalarField alphaEx
        (
            IOobject("staticFaceFractions:alphaExact", runTime.timeName(), mesh,
                IOobject::NO_READ, IOobject::NO_WRITE, false),
            mesh, dimensionedScalar(dimless, Zero), calculatedFvPatchScalarField::typeName
        );
        surfaceScalarField alphafEx
        (
            IOobject("staticFaceFractions:alphafExact", runTime.timeName(), mesh,
                IOobject::NO_READ, IOobject::NO_WRITE, false),
            mesh, dimensionedScalar(dimless, Zero), calculatedFvsPatchScalarField::typeName
        );
        davofInit init(mesh, func(), funcDict);
        if (davofMethod == "exactSphere")
        {
            init.exactSphere(alphaEx, alphafEx);
        }
        else if (davofMethod == "quadraticFaces")
        {
            init.quadraticFaces(alphaEx, alphafEx);
        }
        else
        {
            init.detrixheAslam(alphaEx, alphafEx);
        }
        scalar L1loop = -1, L1static = -1, L1mean = -1, L1exact = -1;
        label nLoop = 0, nStatic = 0;
        label nPlanes = -1, nFail = -1, nPlanesEx = -1, nFailEx = -1;
        if (!evolOnly)
        {
        vectorField mEx;
        scalarField afEx;
        davofInit::areaNormal(mesh, alphafEx, mEx, afEx);

        // the transported face fractions, when written
        {
            IOobject afIO
            (
                IOobject::groupName("alphaf", alpha.group()), runTime.timeName(), mesh,
                IOobject::MUST_READ, IOobject::NO_WRITE, false
            );
            if (afIO.typeHeaderOk<surfaceScalarField>(true))
            {
                const surfaceScalarField afLoop(afIO, mesh);
                L1loop = scoreL1(mesh, afLoop, afEx, mEx, wispTol, nLoop);
            }
        }

        // the static states of alpha and of alpha_exact
        auto makeSurf = [&](const word& nm)
        {
            return surfaceScalarField
            (
                IOobject(nm, runTime.timeName(), mesh,
                    IOobject::NO_READ, IOobject::NO_WRITE, false),
                mesh, dimensionedScalar(dimless, Zero),
                calculatedFvsPatchScalarField::typeName
            );
        };
        surfaceScalarField afU(makeSurf("staticFaceFractions:upwind"));
        surfaceScalarField afM(makeSurf("staticFaceFractions:mean"));
        surfaceScalarField afUEx(makeSurf("staticFaceFractions:upwindExact"));
        surfaceScalarField afMEx(makeSurf("staticFaceFractions:meanExact"));
        staticFractions
        (
            mesh, alpha.primitiveField(), o, phiF, surfCellTol, positionTol,
            afU, afM, nPlanes, nFail
        );
        staticFractions
        (
            mesh, alphaEx.primitiveField(), o, phiF, surfCellTol, positionTol,
            afUEx, afMEx, nPlanesEx, nFailEx
        );
        label nDummy = 0;
        L1static = scoreL1(mesh, afU, afEx, mEx, wispTol, nStatic);
        L1mean = scoreL1(mesh, afM, afEx, mEx, wispTol, nDummy);
        L1exact = scoreL1(mesh, afUEx, afEx, mEx, wispTol, nDummy);
        }

        // the cell-level error
        const scalarField& a = alpha.primitiveField();
        const scalarField& aEx = alphaEx.primitiveField();
        scalar sIn = 0, sOut = 0, sVol = 0;
        label nEx = 0;
        forAll(a, celli)
        {
            const scalar e = mag(a[celli] - aEx[celli]);
            sVol += e*V[celli];
            if (aEx[celli] > surfCellTol && aEx[celli] < 1 - surfCellTol)
            {
                sIn += e;
                ++nEx;
            }
            else
            {
                sOut += e;
            }
        }
        const scalar Ec = sIn/max(nEx, label(1));
        const scalar Eout = sOut/max(nEx, label(1));

        Info<< "    origin " << o << ", L1 loop " << L1loop << " (" << nLoop
            << " faces), static " << L1static << " (" << nStatic << " faces, "
            << nPlanes << " planes, " << nFail << " failed), static mean "
            << L1mean << ", static of alpha_exact " << L1exact << " ("
            << nPlanesEx << " planes, " << nFailEx << " failed), Ec " << Ec
            << ", Eout " << Eout << ", Evol " << sVol << " over " << nEx
            << " exact interface cells" << endl;

        os  << t << '\t' << L1loop << '\t' << nLoop << '\t' << L1static << '\t'
            << nStatic << '\t' << L1mean << '\t' << L1exact << '\t' << nPlanes
            << '\t' << nFail << '\t' << Ec << '\t' << Eout << '\t' << sVol
            << '\t' << nEx << nl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
