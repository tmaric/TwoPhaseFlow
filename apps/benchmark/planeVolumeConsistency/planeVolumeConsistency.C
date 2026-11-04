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
    planeVolumeConsistency

Description
    The volume consistency of the PLIC plane positioning the DAVOF
    reconstruction uses (2026-11-04, Tomislav Maric: "after positioning is
    finished, clip the mesh with the new PLIC interface and sum the
    phase-specific total volume; it should be the same up to machine epsilon
    times the number of interface cells").

    For a written volume fraction alpha and the face fractions alphaf of a
    time directory, in every interface cell of the DAVOF reconstruction
    (|m_c| > wispTol V^(2/3), m_c = sum_f alpha_f S_f^out, the all-equal
    face rule of davof::reconstruct), the plane is positioned exactly as
    davof::positionPlane positions it: surfaceIteratorPLIC::vofCutCell
    (celli, alpha, positionTol, 100, n), the normal n one of
      areaNormal  m_c itself (the DAVOF identity; the first pass of the ring
                  models and the only pass of normalModel none),
      exact       the exact normal of the sphere of the controlDict entry
                  (-entry, default areaFractionError; origin, U0, omega,
                  axisPoint, fork davofOrbit.H) at the cell centre, then
                  re-positioned with the normal at the first plane's
                  centroid (the exact normal model; staticFaceFractions),
      field       the written interfaceNormal (the plane normal the run
                  used, e.g. the ring normal),
    with the production tolerance positionTol = isoFaceTol = 1e-8 (option
    -positionTol). The returned plane (cutValue p along the unit normal from
    the cell centre, liquid where (x - x_c - p n).n > 0) is then clipped
    against the cell twice, independently of the iteration:
      (1) cutCellPLIC::calcSubCell on a fresh object (the pyramid sum of the
          production clip), V_plic;
      (2) the closed-form volume of the plane-cut polyhedron by the
          divergence theorem, V = (1/3) sum_F x_F . S_F over the closed
          boundary (the wet polygon of every cell face, cut by the plane
          with its own edge crossings, plus the cap, whose area vector
          closes the boundary, S_cap = -sum_F S_F^wet), V_ind.
    Reported per time and normal, over the interface cells with
    positionTol < alpha < 1 - positionTol and status 0 (converged):
      - d_c = (V_ind - alpha_c V_c)/V_c: max |d|, mean |d|, the decade
        histogram of |d|, the count above 1e-14 and above positionTol;
      - the signed and the absolute sums of (V_ind - alpha_c V_c), against
        the total liquid volume sum alpha V and against
        N_interface x eps x mean V_c (eps = 2.22e-16);
      - |V_plic - V_ind|/V_c max (the two clips against each other), and
        |S_cap + A_plic|/V^(2/3) max (the cap area vector of (2), out of
        the liquid, against the production cut polygon's faceArea_, which
        cutCell::calcGeomDataCutFace turns into the sub-cell);
      - the residual the iteration returned, |VolumeOfFluid() - alpha|, max,
        against the clip's;
      - the status counts, split by alpha inside/outside [tol, 1 - tol];
        every inside cell whose status is not 0 is printed (cell, alpha,
        normal, plane);
      - the branch of the iteration, from a verbatim copy of vofCutCell with
        the iterates recorded (vofCutCellTraced below; its cutValue is
        checked bitwise against the production call): the bisection hit,
        the coincident vertex levels, the flat bracket (|a1 - a2| < tol: the
        midpoint rule), Newton on the cubic, the secant, the secant
        exhausted (maxIter) or stalled (g1 == g2), and the cells whose
        secant left with a residual above tol but status 0;
      - the orientation: the cells where V_ind is closer to (1 - alpha) V
        than to alpha V (a flip), and on a sphere the cells whose liquid
        sub-cell centroid is farther from the sphere centre than the gas
        sub-cell centroid.

    PRE-REGISTERED EXPECTATION (Tomislav, 2026-11-04, written before the
    runs): per cell |V_cut - alpha_c V_c| < 1e-14 V_c (E1); the sum over the
    interface cells |sum (V_cut - alpha_c V_c)| < N_interface x 1e-15 x V_c
    (E2); every status 0 except the cells with alpha outside
    [tol, 1 - tol] (E3). Each is reported HELD or FALSIFIED with the
    numbers. The agent's reading of the code before the runs: the
    iteration stops at |alpha(p) - alpha| <= positionTol = 1e-8 (the Newton
    loop on the cubic model and the final check), so residuals between
    1e-14 and 1e-8 are admitted by construction; E1 and E2 hold only if the
    cubic model (exact between adjacent vertex levels, up to the LU
    round-off) leaves the first Newton exit at round-off in every cell.

    Option -positioning ccs (2026-11-04): the plane is positioned by the
    Consecutive Cubic Spline positioning of davofCcsPositioning.H instead,
    exactly as davof::positionPlane with the key positioning ccs: from the
    explicit position pPlane_c = (sum_f alpha_f (x_f - x_c).S_f^out -
    3 alpha_c V_c)/|m_c| along the selected normal, on the same cutCellPLIC
    cut, to -ccsTol of V_c (default 1e-14 = the key's default) or the
    round-off floor, at most -ccsMaxIter cuts (default 20); the cells with
    alpha within tol of 0 or 1 are positioned too. Reported in addition:
    the histogram of cuts per cell, the exit codes, the CCS / Newton /
    bisection steps, the end-of-interval cells with and without a plane,
    and the wall time per positioning (the positioning call alone, the
    same for both methods, on the same node).

    PRE-REGISTERED for positioning ccs (2026-11-04, before the runs;
    Tomislav's expectations, the agent's reading in brackets):
      C1 per cell |V_cut - alpha_c V_c| < 1e-14 V_c everywhere [HELD at
         N <= 64; at N = 128 the floor measured with vofCutCell at 1e-14
         reached 3e-14 in some cells, so C1 may fail there by a factor 3];
      C2 |sum| of the order N x eps x V_c [HELD within a factor 20, as the
         1e-14 sweep showed];
      C3 cuts per cell 1 to 2 from the closed-form start [2 to 3: the first
         step is Newton's from the single evaluated point, the second the
         exact Hermite cubic through the two; 1 only where the state is a
         plane within round-off];
      C4 the end-of-interval cells get a plane where |m_c| passes the wisp
         threshold [HELD except where alpha V_c is itself below the floor];
      C5 the cost per positioning is not above vofCutCell's [HELD:
         vofCutCell makes 5 to 6 cuts per cell].

    Runs serial or parallel (-parallel on a decomposed case: the sums and
    maxima are reduced, the per-cell lists are per rank in Pout). Writes
    one row per time and normal to -out (default postProcessing/
    planeVolumeConsistency/planeVolumeConsistency.dat of the case) and
    nothing into the time directories.

    The positioning: H. Scheufler, J. Roenby, Accurate and efficient
    surface reconstruction from volume fraction data on general meshes,
    J. Comput. Phys. 383 (2019) 1-23, doi:10.1016/j.jcp.2019.01.009,
    arXiv:1801.05382 (surfaceIteratorPLIC, cutCellPLIC, cutFacePLIC of the
    TwoPhaseFlow library, DLR). The divergence-theorem volume of a
    polyhedron: B. Mirtich, Fast and accurate computation of polyhedral
    mass properties, J. Graphics Tools 1 (1996) 31-50,
    doi:10.1080/10867651.1996.10487458.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "timeSelector.H"
#include "davofInit.H"
#include "davofOrbit.H"
#include "surfaceIteratorPLIC.H"
#include "cutCellPLIC.H"
#include "davofCcsPositioning.H"
#include "scalarMatrices.H"
#include "OFstream.H"
#include "OSspecific.H"
#include <chrono>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace
{

using namespace Foam;

//- The record of one positioning by the traced copy
struct traceRec
{
    word branch;
    label nNewton = 0;
    label nSecant = 0;
    scalar f1 = 0, f2 = 0, a1 = 0, a2 = 0;   // the bracket after bisection
    scalar fMin = 0, fMax = 0;               // the vertex extrema
    scalar fRet = 0;                         // cutValue of the returned cut
    scalar resRet = 0;                       // |VOF - alpha| of the returned cut
    scalar resNewton = -1;                   // the cut residual after Newton
    DynamicList<Pair<scalar>> it;            // (f, alpha(f)) of every cut
};


//- A verbatim copy of surfaceIteratorPLIC::vofCutCell (TwoPhaseFlow fork,
//  src/VoF/surfaceIterators/surfaceIteratorPLIC.C, as of 2026-11-04) with
//  the branch taken and the iterates recorded. Every cut goes through the
//  same cutCellPLIC::calcSubCell, so the returned cutValue must equal the
//  production call's bitwise; the check is reported.
label vofCutCellTraced
(
    const fvMesh& mesh,
    cutCellPLIC& cutCell,
    const label celli,
    const scalar alpha1,
    const scalar tol,
    const label maxIter,
    vector normal,
    traceRec& T
)
{
    auto cut = [&](const scalar f) -> label
    {
        const label s = cutCell.calcSubCell(celli, f, normal);
        T.it.append(Pair<scalar>(f, cutCell.VolumeOfFluid()));
        return s;
    };

    if (mag(normal) == 0)
    {
        T.branch = "zeroNormal";
        return sign(alpha1 - 0.5);
    }

    normal.normalise();

    const labelList& pLabels = mesh.cellPoints(celli);
    scalarField fvert(pLabels.size());
    forAll(pLabels, pi)
    {
        fvert[pi] = (mesh.points()[pLabels[pi]] - mesh.C()[celli]) & (normal);
    }

    labelList order(fvert.size());
    sortedOrder(fvert, order);
    scalar f1 = fvert[order.first()];
    scalar f2 = fvert[order.last()];
    T.fMin = f1;
    T.fMax = f2;

    if (alpha1 < tol)
    {
        T.branch = "belowTol";
        return -1;
    }
    else if (1 - alpha1 < tol)
    {
        T.branch = "aboveTol";
        return 1;
    }

    label L1 = 0;
    label L2 = fvert.size() - 1;
    scalar a1 = 1;
    scalar a2 = 0;
    scalar L3, f3, a3;

    while (L2 - L1 > 1)
    {
        L3 = round(0.5*(L1 + L2));
        f3 = fvert[order[L3]];
        cut(f3);
        a3 = cutCell.VolumeOfFluid();
        if (a3 > alpha1)
        {
            L1 = L3; f1 = f3; a1 = a3;
        }
        else if (a3 < alpha1)
        {
            L2 = L3; f2 = f3; a2 = a3;
        }
        else
        {
            T.branch = "bisectionHit";
            const label s = cut(f3);
            T.fRet = cutCell.cutValue();
            T.resRet = mag(cutCell.VolumeOfFluid() - alpha1);
            return s;
        }
    }
    T.f1 = f1; T.f2 = f2; T.a1 = a1; T.a2 = a2;

    if (mag(f1 - f2) < 10*SMALL)
    {
        T.branch = "coincidentLevels";
        const label s = cut(f1);
        T.fRet = cutCell.cutValue();
        T.resRet = mag(cutCell.VolumeOfFluid() - alpha1);
        return s;
    }

    if (mag(a1 - a2) < tol)
    {
        T.branch = "flatBracket";
        const label s = cut(0.5*(f1 + f2));
        T.fRet = cutCell.cutValue();
        T.resRet = mag(cutCell.VolumeOfFluid() - alpha1);
        return s;
    }

    f3 = f1 + (f2 - f1)/3;
    cut(f3);
    a3 = cutCell.VolumeOfFluid();

    scalar f4 = f1 + 2*(f2 - f1)/3;
    cut(f4);
    scalar a4 = cutCell.VolumeOfFluid();

    scalarField a(4), f(4), C(4), fOld(4);
    {
        a[0] = a1, a[1] = a3, a[2] = a4, a[3] = a2;
        fOld[0] = f1, fOld[1] = f3, fOld[2] = f4, fOld[3] = f2;
        f[0] = 0, f[1] = (f3-f1)/(f2-f1), f[2] = (f4-f1)/(f2-f1), f[3] = 1;
        scalarSquareMatrix M(4);
        forAll(f, i)
        {
            forAll(f, j)
            {
                M[i][j] = pow(f[i], 3 - j);
            }
        }
        C = a;
        LUsolve(M, C);
    }

    f3 = f[1]; a3 = a[1];
    label nIter = 0;
    scalar res = mag(a3 - alpha1);
    while (res > tol && nIter < 10*maxIter)
    {
        f3 -= (C[0]*pow3(f3) + C[1]*sqr(f3) + C[2]*f3 + C[3] - alpha1)
            /(3*C[0]*sqr(f3) + 2*C[1]*f3 + C[2]);
        a3 = C[0]*pow3(f3) + C[1]*sqr(f3) + C[2]*f3 + C[3];
        res = mag(a3 - alpha1);
        nIter++;
    }
    T.nNewton = nIter;
    f3 = f3*(f2 - f1) + f1;

    label status = cut(f3);
    const scalar VOF = cutCell.VolumeOfFluid();
    res = mag(VOF - alpha1);
    T.resNewton = res;

    if (res > tol)
    {
    }
    else
    {
        T.branch = "newton";
        T.fRet = cutCell.cutValue();
        T.resRet = res;
        return status;
    }

    scalar x2 = f3;
    scalar g2 = VOF - alpha1;
    scalar x1 = max(1e-3*(f2 - f1), 100*SMALL);
    x1 = max(x1, f1);
    x1 = min(x1, f2);
    cut(x1);
    scalar g1 = cutCell.VolumeOfFluid() - alpha1;

    nIter = 0;
    scalar g0(0), x0(0);
    while (res > tol && nIter < maxIter && g1 != g2)
    {
        x0 = (x2*g1 - x1*g2)/(g1 - g2);
        status = cut(x0);
        g0 = cutCell.VolumeOfFluid() - alpha1;
        res = mag(g0);
        x2 = x1; g2 = g1;
        x1 = x0; g1 = g0;
        nIter++;
    }
    T.nSecant = nIter;
    T.fRet = cutCell.cutValue();
    T.resRet = mag(cutCell.VolumeOfFluid() - alpha1);
    if (res <= tol)
    {
        T.branch = "secant";
    }
    else if (nIter == 0)
    {
        // g1 == g2 at entry: the loop never ran, the cut object holds the
        // cut at x1 while status is the cut's at f3
        T.branch = "secantStalledAtEntry";
    }
    else if (nIter >= maxIter)
    {
        T.branch = "secantExhausted";
    }
    else
    {
        T.branch = "secantStalled";
    }

    return status;
}


//- The liquid volume of celli on the side (x - base).n > 0 of the plane,
//  by the divergence theorem over the closed boundary of the clipped cell:
//  V = (1/3) [sum_F x_F . S_F + base . S_cap], S_cap = -sum_F S_F, with F
//  the wet polygon of every cell face (the vertices with (p - base).n >= 0
//  and the crossings of the edges whose end values differ in sign), its
//  area vector by the fan about its first vertex, oriented out of the
//  cell. Independent of cutFacePLIC and cutCellPLIC.
scalar clippedVolume
(
    const fvMesh& mesh,
    const label celli,
    const vector& n,
    const point& base,
    vector& Scap
)
{
    const pointField& P = mesh.points();
    const faceList& faces = mesh.faces();
    const labelUList& own = mesh.faceOwner();

    vector sumS(Zero);
    scalar sumXS = 0;
    DynamicList<point> poly(16);

    for (const label facei : mesh.cells()[celli])
    {
        const face& f = faces[facei];
        poly.clear();
        forAll(f, i)
        {
            const point& pi = P[f[i]];
            const point& pj = P[f.nextLabel(i)];
            const scalar di = (pi - base) & n;
            const scalar dj = (pj - base) & n;
            if (di >= 0)
            {
                poly.append(pi);
            }
            if ((di > 0 && dj < 0) || (di < 0 && dj > 0))
            {
                poly.append(pi + (pj - pi)*(di/(di - dj)));
            }
        }
        if (poly.size() < 3) continue;

        vector S(Zero);
        for (label k = 1; k + 1 < poly.size(); ++k)
        {
            S += 0.5*((poly[k] - poly[0]) ^ (poly[k + 1] - poly[0]));
        }
        if (own[facei] != celli)
        {
            S = -S;
        }
        sumS += S;
        sumXS += poly[0] & S;
    }

    Scap = -sumS;
    return (sumXS + (base & Scap))/3.0;
}


//- The decade bin of |d|: 0 for |d| < 1e-16, k for 1e-(16-k) <= |d| < 1e-(15-k), 9 for >= 1e-8
label decadeBin(const scalar d)
{
    if (!(d >= 1e-16)) return 0;
    const label b = label(std::floor(Foam::log10(d))) + 17;   // 1e-16 -> 1
    return min(max(b, label(1)), label(9));
}

} // End anonymous namespace


int main(int argc, char *argv[])
{
    argList::addNote
    (
        "The volume consistency of the DAVOF plane positioning: every"
        " interface cell's plane positioned as davof::positionPlane does,"
        " then clipped independently; the mismatch per cell and summed"
    );
    timeSelector::addOptions();
    argList::addOption
    (
        "entry", "word",
        "the entry of system/controlDict's functions that holds the exact"
        " sphere (default areaFractionError)"
    );
    argList::addOption
    (
        "normals", "list",
        "the normals to test, a subset of (areaNormal exact field)"
        " (default (areaNormal exact))"
    );
    argList::addOption("positionTol", "scalar", "the positioning tolerance in alpha (default 1e-8 = isoFaceTol)");
    argList::addOption("positioning", "word", "vofCutCell (default) | ccs");
    argList::addOption("ccsTol", "scalar", "ccs: the tolerance of V_c (default 1e-14)");
    argList::addOption("ccsMaxIter", "label", "ccs: the most cuts (default 20)");
    argList::addOption("wispTol", "scalar", "the interface-cell threshold (default 1e-3)");
    argList::addOption("nPrint", "label", "the worst cells printed per rank (default 10)");
    argList::addOption
    (
        "out", "file",
        "the output file (default postProcessing/planeVolumeConsistency/"
        "planeVolumeConsistency.dat of the case)"
    );

    #include "setRootCase.H"
    #include "createTime.H"
    instantList timeDirs = timeSelector::select0(runTime, args);
    #include "createMesh.H"

    const word entry(args.getOrDefault<word>("entry", "areaFractionError"));
    const scalar positionTol = args.getOrDefault<scalar>("positionTol", 1e-8);
    const word positioning = args.getOrDefault<word>("positioning", "vofCutCell");
    const scalar ccsTol = args.getOrDefault<scalar>("ccsTol", 1e-14);
    const label ccsMaxIter = args.getOrDefault<label>("ccsMaxIter", 30);
    const bool ccs = (positioning == "ccs");
    if (!ccs && positioning != "vofCutCell")
    {
        FatalErrorInFunction
            << "positioning " << positioning << " unknown: vofCutCell, ccs"
            << exit(FatalError);
    }
    const scalar wispTol = args.getOrDefault<scalar>("wispTol", 1e-3);
    const label nPrint = args.getOrDefault<label>("nPrint", 10);
    wordList normals({"areaNormal", "exact"});
    if (args.found("normals"))
    {
        normals = args.getList<word>("normals");
    }

    const dictionary& foDict =
        runTime.controlDict().subDict("functions").subDict(entry);
    const word fieldName(foDict.get<word>("field"));
    const point origin0(foDict.get<point>("origin"));
    const vector U0 = foDict.getOrDefault<vector>("U0", vector::zero);
    vector omega;
    point axisPoint;
    davof::readOrbit(foDict, origin0, omega, axisPoint);

    fileName outFile
    (
        runTime.globalPath()/"postProcessing"/"planeVolumeConsistency"
       /"planeVolumeConsistency.dat"
    );
    if (args.found("out"))
    {
        outFile = args.get<fileName>("out");
    }
    autoPtr<OFstream> osPtr;
    if (Pstream::master())
    {
        mkDir(outFile.path());
        osPtr.reset(new OFstream(outFile));
        osPtr().precision(12);
        osPtr()
            << "# planeVolumeConsistency (entry " << entry << ", positionTol "
            << positionTol << ", wispTol " << wispTol << ", positioning "
            << positioning << ", ccsTol " << ccsTol << ", ranks "
            << Pstream::nProcs() << ")" << nl
            << "# Time\tnormal\tnInterface\tnInside\tnConverged\tnInsideNotConv"
            << "\tnOutside\tmaxAbsD\tmeanAbsD\tnAbove1e-14\tnAboveTol"
            << "\tsumSigned\tsumAbs\tVliquid\tNepsV\tmaxPlicVsInd\tmaxCapArea"
            << "\tmaxResIter\tnNewton\tnSecant\tnSecantBad\tnFlat\tnCoinc"
            << "\tnBisectHit\tnCopyMismatch\tnFlip\tnSphereOrient\tE1\tE2\tE3"
            << "\tpositioning\tnsPerPos\tcutsPerCell\tmaxCuts\tnFloor\tnEnd\tnEndNoPlane"
            << nl;
    }

    const vectorField& C = mesh.cellCentres();
    const scalarField& V = mesh.V();
    const cellList& cells = mesh.cells();
    const label nInt = mesh.nInternalFaces();
    const scalar eps = 2.220446049250313e-16;

    forAll(timeDirs, timei)
    {
        runTime.setTime(timeDirs[timei], timei);
        const scalar t = runTime.value();
        Info<< "Time = " << runTime.timeName() << endl;

        volScalarField alpha
        (
            IOobject(fieldName, runTime.timeName(), mesh, IOobject::MUST_READ, IOobject::NO_WRITE, false),
            mesh
        );
        surfaceScalarField alphaf
        (
            IOobject
            (
                IOobject::groupName("alphaf", alpha.group()), runTime.timeName(), mesh,
                IOobject::MUST_READ, IOobject::NO_WRITE, false
            ),
            mesh
        );
        autoPtr<volVectorField> fieldNormalPtr;
        if (normals.found("field"))
        {
            fieldNormalPtr.reset
            (
                new volVectorField
                (
                    IOobject
                    (
                        IOobject::groupName("interfaceNormal", alpha.group()),
                        runTime.timeName(), mesh, IOobject::MUST_READ, IOobject::NO_WRITE, false
                    ),
                    mesh
                )
            );
        }
        const scalarField& a = alpha.primitiveField();
        const point o = davof::orbitOrigin(origin0, U0, omega, axisPoint, t);

        // The interface cells as davof::reconstruct selects them
        vectorField m;
        scalarField afAll;
        davofInit::areaNormal(mesh, alphaf, m, afAll);
        boolList hasValue(mesh.nFaces(), false);
        for (label facei = 0; facei < nInt; ++facei) hasValue[facei] = true;
        forAll(alphaf.boundaryField(), patchi)
        {
            const polyPatch& pp = mesh.boundaryMesh()[patchi];
            if (alphaf.boundaryField()[patchi].size() == pp.size())
            {
                forAll(pp, i) hasValue[pp.start() + i] = true;
            }
        }
        DynamicList<label> interfaceLabels(mesh.nCells()/10 + 16);
        forAll(m, celli)
        {
            if (m[celli] == vector::zero) continue;
            scalar lo = GREAT, hi = -GREAT;
            for (const label facei : cells[celli])
            {
                if (!hasValue[facei]) continue;
                lo = min(lo, afAll[facei]);
                hi = max(hi, afAll[facei]);
            }
            if (lo == hi) continue;
            if (mag(m[celli]) > davofInit::wispThreshold(mesh, celli, wispTol))
            {
                interfaceLabels.append(celli);
            }
        }
        const label nInterface = returnReduce(interfaceLabels.size(), sumOp<label>());
        scalar Vliquid = 0;
        forAll(a, celli) Vliquid += a[celli]*V[celli];
        Vliquid = returnReduce(Vliquid, sumOp<scalar>());
        Info<< "    " << nInterface << " interface cells, liquid volume "
            << Vliquid << ", sphere centre " << o << endl;

        for (const word& normalName : normals)
        {
            surfaceIteratorPLIC sIter(mesh, positionTol);
            surfaceIteratorPLIC sIterPre(mesh, positionTol);
            cutCellPLIC ccTrace(mesh);
            cutCellPLIC ccClip(mesh);

            // the counters
            label nInside = 0, nOutside = 0, nConverged = 0, nInsideNotConv = 0;
            label nNoNormal = 0, nOutsideWrongStatus = 0;
            label nAbove14 = 0, nAboveTol = 0, nCopyMismatch = 0;
            label nFlip = 0, nSphereOrient = 0, nSphereTested = 0;
            label nNewton = 0, nSecant = 0, nSecantBad = 0, nFlat = 0, nCoinc = 0, nBisect = 0;
            label nSecantStalledEntry = 0, nOutsideBracket = 0;
            label maxNewtonIter = 0, maxSecantIter = 0;
            labelList hist(10, Zero);
            // ccs: the cuts per cell (1..7, 8+), the exit codes (1..4), the
            // step kinds, the end-of-interval cells, the timing
            cutCellPLIC ccPos(mesh);
            labelList evalHist(9, Zero);
            labelList exitHist(5, Zero);
            label nCcsHermite = 0, nCcsNewton = 0, nCcsBisect = 0, maxEvals = 0, sumEvals = 0;
            label nEndCcs = 0, nEndNoPlane = 0;
            scalar maxStartRes = 0, sumStartRes = 0;
            std::chrono::steady_clock::duration tPos(0);
            label nTimed = 0;
            DynamicList<label> sortedFaces(16);
            auto explicitStart = [&](const label celli) -> scalar
            {
                // pPlane_c as davof::reconstruct forms it (the faces in
                // ascending order)
                sortedFaces.clear();
                sortedFaces.append(cells[celli]);
                Foam::sort(sortedFaces);
                scalar faceMoment = 0;
                for (const label facei : sortedFaces)
                {
                    if (!hasValue[facei]) continue;
                    const scalar af = afAll[facei];
                    if (af == 0) continue;
                    const scalar term = af*((mesh.faceCentres()[facei] - C[celli]) & mesh.faceAreas()[facei]);
                    faceMoment += (mesh.faceOwner()[facei] == celli) ? term : -term;
                }
                return (faceMoment - 3.0*a[celli]*V[celli])/mag(m[celli]);
            };
            scalar maxAbsD = 0, sumAbsD = 0, sumSigned = 0, sumAbs = 0;
            scalar maxPlicVsInd = 0, maxCapArea = 0, maxResIter = 0, maxResIterVsClip = 0;
            scalar sumVint = 0;
            scalar maxOutsideRes = 0;
            label worstCell = -1;
            scalar worstD = 0;
            // the worst cells of this rank: (|d|, celli) and the normal used
            DynamicList<Pair<scalar>> worst(interfaceLabels.size());
            DynamicList<vector> worstN(interfaceLabels.size());

            for (const label celli : interfaceLabels)
            {
                const scalar ac = a[celli];
                const scalar Vc = V[celli];
                const scalar h = Foam::cbrt(Vc);

                vector nRaw(Zero);
                if (normalName == "areaNormal")
                {
                    nRaw = m[celli];
                }
                else if (normalName == "exact")
                {
                    nRaw = o - C[celli];
                    if (mag(nRaw) > 0 && ac >= positionTol && 1 - ac >= positionTol)
                    {
                        const label s1 = sIterPre.vofCutCell(celli, ac, positionTol, 100, nRaw);
                        if (s1 == 0 && mag(sIterPre.surfaceArea()) > 0)
                        {
                            const vector n2 = o - sIterPre.surfaceCentre();
                            if (mag(n2) > 0) nRaw = n2;
                        }
                    }
                }
                else if (normalName == "field")
                {
                    nRaw = fieldNormalPtr()[celli];
                }
                else
                {
                    FatalErrorInFunction
                        << "normal " << normalName << " unknown: areaNormal, exact, field"
                        << exit(FatalError);
                }
                if (mag(nRaw) == 0)
                {
                    ++nNoNormal;
                    continue;
                }

                // (1) the production positioning, timed (the call alone)
                vector nUnit = nRaw;
                nUnit.normalise();
                const bool atEnd = (ac < positionTol) || (1 - ac < positionTol);
                label status = 0;
                scalar p = 0, vofRet = 0, Vl = 0;
                point xl(Zero);
                Foam::davof::ccsRecord rec;
                if (ccs)
                {
                    const scalar p0 = explicitStart(celli);
                    const auto t0 = std::chrono::steady_clock::now();
                    status = Foam::davof::ccsPositionPlane
                    (
                        mesh, ccPos, celli, ac, nUnit, p0, ccsTol, ccsMaxIter, rec
                    );
                    tPos += std::chrono::steady_clock::now() - t0;
                    ++nTimed;
                    p = ccPos.cutValue();
                    vofRet = (status == 0) ? ccPos.VolumeOfFluid() : ((status < 0) ? 1 : 0);
                    Vl = (status == 0) ? ccPos.subCellVolume() : ((status < 0) ? Vc : 0);
                    xl = ccPos.subCellCentre();
                    ++evalHist[min(rec.nEvals, evalHist.size() - 1)];
                    ++exitHist[min(max(rec.exitCode, label(0)), exitHist.size() - 1)];
                    nCcsHermite += rec.nHermite; nCcsNewton += rec.nNewton; nCcsBisect += rec.nBisect;
                    maxEvals = max(maxEvals, rec.nEvals);
                    sumEvals += rec.nEvals;
                    maxStartRes = max(maxStartRes, mag(rec.rStart));
                    sumStartRes += mag(rec.rStart);
                    if (atEnd)
                    {
                        ++nEndCcs;
                        if (status != 0 || mag(ccPos.faceArea()) == 0)
                        {
                            ++nEndNoPlane;
                            Pout<< "END CELL WITHOUT PLANE: cell " << celli << " alpha " << ac
                                << " |m|/h^2 " << mag(m[celli])/(h*h) << " status " << status
                                << " cuts " << rec.nEvals << " exit " << rec.exitCode
                                << " rStart " << rec.rStart << " rFinal " << rec.rFinal << nl;
                            continue;
                        }
                    }
                    ++nInside;
                    if (status != 0)
                    {
                        ++nInsideNotConv;
                        Pout<< "NOT CONVERGED (ccs, status " << status << "): cell " << celli
                            << " alpha " << ac << " n " << nRaw << " cuts " << rec.nEvals
                            << " exit " << rec.exitCode << " pStart/h " << rec.pStart/h
                            << " rStart " << rec.rStart << " pFinal/h " << rec.pFinal/h
                            << " rFinal " << rec.rFinal << nl;
                        continue;
                    }
                    ++nConverged;
                    maxResIter = max(maxResIter, mag(rec.rFinal));
                }
                else
                {
                    const auto t0 = std::chrono::steady_clock::now();
                    status = sIter.vofCutCell(celli, ac, positionTol, 100, nRaw);
                    tPos += std::chrono::steady_clock::now() - t0;
                    ++nTimed;
                    const bool inside = !atEnd;

                    if (!inside)
                    {
                        ++nOutside;
                        const label expected = (ac < positionTol) ? -1 : 1;
                        if (status != expected) ++nOutsideWrongStatus;
                        maxOutsideRes = max(maxOutsideRes, (ac < positionTol) ? ac : 1 - ac);
                        continue;
                    }
                    ++nInside;

                    // the traced copy, for the branch and the iterates
                    traceRec T;
                    const label st2 = vofCutCellTraced(mesh, ccTrace, celli, ac, positionTol, 100, nRaw, T);
                    if (st2 != status || (status == 0 && T.fRet != sIter.cutValue()))
                    {
                        ++nCopyMismatch;
                    }
                    if (T.branch == "newton") ++nNewton;
                    else if (T.branch == "secant") ++nSecant;
                    else if (T.branch == "flatBracket") ++nFlat;
                    else if (T.branch == "coincidentLevels") ++nCoinc;
                    else if (T.branch == "bisectionHit") ++nBisect;
                    else if (T.branch == "secantStalledAtEntry") { ++nSecantBad; ++nSecantStalledEntry; }
                    else ++nSecantBad;
                    maxNewtonIter = max(maxNewtonIter, T.nNewton);
                    maxSecantIter = max(maxSecantIter, T.nSecant);

                    if (status != 0)
                    {
                        ++nInsideNotConv;
                        Pout<< "NOT CONVERGED (status " << status << "): cell " << celli
                            << " alpha " << ac << " n " << nRaw << " branch " << T.branch
                            << " nNewton " << T.nNewton << " nSecant " << T.nSecant
                            << " bracket f [" << T.f1 << ", " << T.f2 << "] a [" << T.a2 << ", " << T.a1
                            << "] vertex f [" << T.fMin << ", " << T.fMax << "]" << nl;
                        forAll(T.it, k)
                        {
                            Pout<< "    cut " << k << " f/h " << T.it[k].first()/h
                                << " alpha(f) " << T.it[k].second()
                                << " r " << T.it[k].second() - ac << nl;
                        }
                        continue;
                    }
                    ++nConverged;

                    p = sIter.cutValue();
                    vofRet = sIter.VolumeOfFluid();
                    Vl = sIter.subCellVolume();
                    xl = sIter.subCellCentre();
                    const scalar resIter = mag(vofRet - ac);
                    maxResIter = max(maxResIter, resIter);
                    if (p < T.fMin || p > T.fMax) ++nOutsideBracket;
                    if (T.branch != "newton" && T.branch != "secant" && T.resRet > positionTol)
                    {
                        // a non-converged plane returned with status 0
                        Pout<< "RESIDUAL ABOVE TOL WITH STATUS 0: cell " << celli
                            << " alpha " << ac << " n " << nRaw << " branch " << T.branch
                            << " res " << T.resRet << " nNewton " << T.nNewton
                            << " nSecant " << T.nSecant << nl;
                    }
                }

                // (2) the production clip on a fresh object
                const label st3 = ccClip.calcSubCell(celli, p, nUnit);
                const scalar Vplic = (st3 == 0) ? ccClip.subCellVolume() : ((st3 < 0) ? Vc : 0);
                const vector Aplic = (st3 == 0) ? ccClip.faceArea() : vector(Zero);
                maxResIterVsClip = max(maxResIterVsClip, mag(Vplic/Vc - vofRet));

                // (3) the independent clip
                vector Scap;
                const point base = C[celli] + nUnit*p;
                const scalar Vind = clippedVolume(mesh, celli, nUnit, base, Scap);

                const scalar d = (Vind - ac*Vc)/Vc;
                const scalar ad = mag(d);
                sumVint += Vc;
                sumSigned += (Vind - ac*Vc);
                sumAbs += mag(Vind - ac*Vc);
                sumAbsD += ad;
                if (ad > maxAbsD) { maxAbsD = ad; worstCell = celli; worstD = d; }
                if (ad >= 1e-14) ++nAbove14;
                if (ad > positionTol) ++nAboveTol;
                ++hist[decadeBin(ad)];
                maxPlicVsInd = max(maxPlicVsInd, mag(Vplic - Vind)/Vc);
                // cutCell::calcGeomDataCutFace turns faceArea_ into the
                // sub-cell (into the liquid); S_cap points out of it
                maxCapArea = max(maxCapArea, mag(Scap + Aplic)/(h*h));
                worst.append(Pair<scalar>(ad, scalar(celli)));
                worstN.append(nRaw);

                // the orientation
                if (mag(Vind/Vc - (1 - ac)) < ad) ++nFlip;
                if (Vl > 0.01*Vc && Vc - Vl > 0.01*Vc)
                {
                    const point xg = (C[celli]*Vc - xl*Vl)/(Vc - Vl);
                    ++nSphereTested;
                    if (mag(xl - o) > mag(xg - o)) ++nSphereOrient;
                }
            }

            // the worst cells of this rank
            {
                labelList order(worst.size());
                scalarField key(worst.size());
                forAll(worst, k) key[k] = -worst[k].first();
                sortedOrder(key, order);
                const label nShow = min(nPrint, order.size());
                for (label k = 0; k < nShow; ++k)
                {
                    const label celli = label(worst[order[k]].second());
                    const scalar ac = a[celli];
                    const vector nRaw = worstN[order[k]];
                    if (ccs)
                    {
                        vector nUnit = nRaw;
                        nUnit.normalise();
                        Foam::davof::ccsRecord rec;
                        const scalar p0 = explicitStart(celli);
                        const scalar hh = Foam::cbrt(V[celli]);
                        Foam::davof::ccsPositionPlane(mesh, ccPos, celli, ac, nUnit, p0, ccsTol, ccsMaxIter, rec);
                        Pout<< "worst " << k << " (" << normalName << ", ccs): cell " << celli
                            << " alpha " << ac << " |d| " << worst[order[k]].first()
                            << " cuts " << rec.nEvals << " exit " << rec.exitCode
                            << " hermite " << rec.nHermite << " newton " << rec.nNewton
                            << " bisect " << rec.nBisect << " pStart/h " << rec.pStart/hh
                            << " rStart " << rec.rStart << " pFinal/h " << rec.pFinal/hh
                            << " rFinal " << rec.rFinal << nl;
                        continue;
                    }
                    traceRec T;
                    vofCutCellTraced(mesh, ccTrace, celli, ac, positionTol, 100, nRaw, T);
                    Pout<< "worst " << k << " (" << normalName << "): cell " << celli
                        << " alpha " << ac << " |d| " << worst[order[k]].first()
                        << " branch " << T.branch << " nNewton " << T.nNewton
                        << " nSecant " << T.nSecant << " resNewton " << T.resNewton
                        << " resRet " << T.resRet << " nCuts " << T.it.size()
                        << " bracket a [" << T.a2 << ", " << T.a1 << "] f/h ["
                        << T.f1/Foam::cbrt(V[celli]) << ", " << T.f2/Foam::cbrt(V[celli]) << "]" << nl;
                }
            }

            // the reductions
            reduce(nInside, sumOp<label>());
            reduce(nOutside, sumOp<label>());
            reduce(nConverged, sumOp<label>());
            reduce(nInsideNotConv, sumOp<label>());
            reduce(nNoNormal, sumOp<label>());
            reduce(nOutsideWrongStatus, sumOp<label>());
            reduce(nAbove14, sumOp<label>());
            reduce(nAboveTol, sumOp<label>());
            reduce(nCopyMismatch, sumOp<label>());
            reduce(nFlip, sumOp<label>());
            reduce(nSphereOrient, sumOp<label>());
            reduce(nSphereTested, sumOp<label>());
            reduce(nNewton, sumOp<label>());
            reduce(nSecant, sumOp<label>());
            reduce(nSecantBad, sumOp<label>());
            reduce(nSecantStalledEntry, sumOp<label>());
            reduce(nFlat, sumOp<label>());
            reduce(nCoinc, sumOp<label>());
            reduce(nBisect, sumOp<label>());
            reduce(nOutsideBracket, sumOp<label>());
            reduce(maxNewtonIter, maxOp<label>());
            reduce(maxSecantIter, maxOp<label>());
            forAll(hist, b) reduce(hist[b], sumOp<label>());
            reduce(sumAbsD, sumOp<scalar>());
            reduce(sumSigned, sumOp<scalar>());
            reduce(sumAbs, sumOp<scalar>());
            reduce(sumVint, sumOp<scalar>());
            const scalar maxAbsDAll = returnReduce(maxAbsD, maxOp<scalar>());
            reduce(maxPlicVsInd, maxOp<scalar>());
            reduce(maxCapArea, maxOp<scalar>());
            reduce(maxResIter, maxOp<scalar>());
            reduce(maxResIterVsClip, maxOp<scalar>());
            reduce(maxOutsideRes, maxOp<scalar>());
            // the ccs counters and the timing (per rank mean, the slowest rank)
            forAll(evalHist, b) reduce(evalHist[b], sumOp<label>());
            forAll(exitHist, b) reduce(exitHist[b], sumOp<label>());
            reduce(nCcsHermite, sumOp<label>()); reduce(nCcsNewton, sumOp<label>());
            reduce(nCcsBisect, sumOp<label>()); reduce(maxEvals, maxOp<label>());
            reduce(sumEvals, sumOp<label>()); reduce(nEndCcs, sumOp<label>());
            reduce(nEndNoPlane, sumOp<label>());
            reduce(maxStartRes, maxOp<scalar>()); reduce(sumStartRes, sumOp<scalar>());
            scalar nsPerPos =
                nTimed > 0
              ? scalar(std::chrono::duration_cast<std::chrono::nanoseconds>(tPos).count())/nTimed
              : 0;
            scalar nsPerPosMax = nsPerPos;
            reduce(nsPerPosMax, maxOp<scalar>());
            label nTimedAll = returnReduce(nTimed, sumOp<label>());
            scalar nsSum = nsPerPos*nTimed;
            reduce(nsSum, sumOp<scalar>());
            const scalar nsPerPosMean = nsSum/max(nTimedAll, label(1));
            if (worstCell >= 0 && maxAbsD == maxAbsDAll)
            {
                Pout<< "global worst cell " << worstCell << " alpha " << a[worstCell]
                    << " d " << worstD << " V " << V[worstCell] << " C " << C[worstCell] << nl;
            }

            const scalar meanAbsD = sumAbsD/max(nConverged, label(1));
            const scalar meanV = sumVint/max(nConverged, label(1));
            const scalar NepsV = nConverged*eps*meanV;
            const bool E1 = (nConverged > 0) && (maxAbsDAll < 1e-14);
            const bool E2 = (nConverged > 0) && (mag(sumSigned) < nConverged*1e-15*meanV);
            const bool E3 = (nInsideNotConv == 0) && (nSecantBad == 0) && (nOutsideWrongStatus == 0);

            Info<< nl << "  normal " << normalName << ": " << nInterface
                << " interface cells, " << nInside << " with alpha inside [tol, 1 - tol], "
                << nOutside << " outside (" << nOutsideWrongStatus
                << " with a status other than the expected -1/+1, max alpha residual "
                << maxOutsideRes << "), " << nNoNormal << " without a normal" << nl
                << "    converged (status 0): " << nConverged << ", not converged: "
                << nInsideNotConv << ", traced copy mismatches: " << nCopyMismatch << nl
                << "    branches: newton " << nNewton << " (max Newton iterations " << maxNewtonIter
                << "), secant " << nSecant << " (max secant iterations " << maxSecantIter
                << "), secant left above tol or stalled " << nSecantBad
                << " (stalled at entry " << nSecantStalledEntry << "), flat bracket " << nFlat
                << ", coincident levels " << nCoinc << ", bisection hit " << nBisect
                << ", plane outside the vertex range " << nOutsideBracket << nl
                << "    per cell d = (V_ind - alpha V)/V: max |d| " << maxAbsDAll
                << ", mean |d| " << meanAbsD << ", cells >= 1e-14: " << nAbove14
                << ", cells > tol: " << nAboveTol << nl
                << "    histogram of |d| [<1e-16";
            for (label b = 1; b < 10; ++b)
            {
                Info<< ", 1e-" << (17 - b) << ".." << ((b == 9) ? "inf" : "");
                if (b < 9) Info<< "1e-" << (16 - b);
            }
            Info<< "]: " << hist << nl
                << "    sums over the converged cells: signed " << sumSigned
                << ", absolute " << sumAbs << ", liquid volume " << Vliquid
                << " (signed/Vliquid " << sumSigned/max(Vliquid, VSMALL)
                << "), N eps V_mean " << NepsV << " (signed/NepsV "
                << sumSigned/max(NepsV, VSMALL) << ", absolute/NepsV "
                << sumAbs/max(NepsV, VSMALL) << ")" << nl
                << "    clips: max |V_plic - V_ind|/V " << maxPlicVsInd
                << ", max |S_cap - A_plic|/h^2 " << maxCapArea
                << ", max residual returned by the iteration " << maxResIter
                << ", max |V_plic/V - VOF_returned| " << maxResIterVsClip << nl
                << "    orientation: flips " << nFlip << ", sphere-side violations "
                << nSphereOrient << " of " << nSphereTested << " tested" << nl
                << "    E1 (per cell |d| < 1e-14): " << (E1 ? "HELD" : "FALSIFIED")
                << " (max |d| " << maxAbsDAll << ", " << nAbove14 << " cells at or above)" << nl
                << "    E2 (|sum| < N 1e-15 V): " << (E2 ? "HELD" : "FALSIFIED")
                << " (|sum| " << mag(sumSigned) << " against " << nConverged*1e-15*meanV << ")" << nl
                << "    E3 (every inside status 0): " << (E3 ? "HELD" : "FALSIFIED")
                << " (" << nInsideNotConv << " not converged, " << nSecantBad
                << " secant bad)" << nl
                << "    timing (" << positioning << "): " << nsPerPosMean
                << " ns per positioning (mean over " << nTimedAll
                << " calls; slowest rank " << nsPerPosMax << ")" << endl;
            if (ccs)
            {
                const scalar meanEvals = scalar(sumEvals)/max(nInside + nEndNoPlane, label(1));
                const bool C1 = E1;
                const bool C2 = (nConverged > 0) && (mag(sumSigned) <= 20*NepsV);
                const bool C3 = (evalHist[3] + evalHist[4] + evalHist[5] + evalHist[6]
                               + evalHist[7] + evalHist[8] == 0);
                const bool C3b = (evalHist[4] + evalHist[5] + evalHist[6] + evalHist[7] + evalHist[8] == 0);
                const bool C4 = (nEndNoPlane == 0);
                Info<< "    ccs: cuts per cell mean " << meanEvals << " max " << maxEvals
                    << ", histogram [1..7, 8+] " << evalHist
                    << ", exits [0, tol, floor, closed, maxIter] " << exitHist
                    << ", steps hermite " << nCcsHermite << " newton " << nCcsNewton
                    << " bisect " << nCcsBisect << nl
                    << "    ccs: start residual |V(p0) - alpha V|/V max " << maxStartRes
                    << " mean " << sumStartRes/max(nInside + nEndNoPlane, label(1))
                    << "; end-of-interval cells " << nEndCcs << " (" << nEndNoPlane
                    << " without a plane)" << nl
                    << "    C1 (per cell |d| < 1e-14): " << (C1 ? "HELD" : "FALSIFIED")
                    << " (max |d| " << maxAbsDAll << ")" << nl
                    << "    C2 (|sum| <= 20 N eps V): " << (C2 ? "HELD" : "FALSIFIED")
                    << " (|sum|/(N eps V) " << mag(sumSigned)/max(NepsV, VSMALL) << ")" << nl
                    << "    C3 (cuts per cell 1 to 2): " << (C3 ? "HELD" : "FALSIFIED")
                    << "; agent's reading (1 to 3): " << (C3b ? "HELD" : "FALSIFIED") << nl
                    << "    C4 (every end-of-interval cell has a plane): "
                    << (C4 ? "HELD" : "FALSIFIED") << " (" << nEndNoPlane << " of "
                    << nEndCcs << " without)" << endl;
            }

            if (osPtr)
            {
                osPtr()
                    << t << '\t' << normalName << '\t' << nInterface << '\t' << nInside
                    << '\t' << nConverged << '\t' << nInsideNotConv << '\t' << nOutside
                    << '\t' << maxAbsDAll << '\t' << meanAbsD << '\t' << nAbove14
                    << '\t' << nAboveTol << '\t' << sumSigned << '\t' << sumAbs
                    << '\t' << Vliquid << '\t' << NepsV << '\t' << maxPlicVsInd
                    << '\t' << maxCapArea << '\t' << maxResIter << '\t' << nNewton
                    << '\t' << nSecant << '\t' << nSecantBad << '\t' << nFlat
                    << '\t' << nCoinc << '\t' << nBisect << '\t' << nCopyMismatch
                    << '\t' << nFlip << '\t' << nSphereOrient
                    << '\t' << (E1 ? "HELD" : "FALSIFIED")
                    << '\t' << (E2 ? "HELD" : "FALSIFIED")
                    << '\t' << (E3 ? "HELD" : "FALSIFIED")
                    << '\t' << positioning << '\t' << nsPerPosMean
                    << '\t' << scalar(sumEvals)/max(nInside + nEndNoPlane, label(1))
                    << '\t' << maxEvals << '\t' << exitHist[2] << '\t' << nEndCcs
                    << '\t' << nEndNoPlane << nl;
            }
        }
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
