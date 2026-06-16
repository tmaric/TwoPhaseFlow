/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  |
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
    acousticDropEquilibriumFoam

Author
    Chuanchao Xu, MMA, TU Darmstadt
    Email: xu@mma.tu-darmstadt.de

Description
    Pseudo-time acoustic drop levitation equilibrium solver.

    Pseudo-solver following Andrade et al. Fig. 2: for each fixed drop shape,
    solve the frequency-domain acoustic field and iterate the vertical drop
    position until force balance is reached; then apply one Young-Laplace shape
    correction and repeat the fixed-shape position routine for the updated
    interface.

    Base solver: MPI-capable block-coupled frequency-domain acoustic solver.
    Uses the same pressure block structure as acousticHelmholtzSerialFoam:

        [ A  -(B1 + B2) ] [Pim] = [bPim]
        [ (B1 + B2)  A ] [Pre]   [bPre]

    but assembles A and coupling contributions on decomposed subdomains,
    explicitly including processor-interface couplings in the PETSc matrix.
    PETSc then solves the global distributed linear system (default:
    preonly+lu+mumps).

    Practical difference to acousticHelmholtzSerialFoam:
    - acousticHelmholtzSerialFoam: serial OpenFOAM assembly only (reference).
    - acousticHelmholtzFoam: distributed OpenFOAM assembly + distributed PETSc solve.
\*---------------------------------------------------------------------------*/

#include <petscksp.h>
#include "fvCFD.H"
#include "fvOptions.H"
#include "simpleControl.H"
#include "emptyPolyPatch.H"
#include "wedgePolyPatch.H"
#include "isoAdvection.H"
#include "cutFaceAdvect.H"
#include "surfaceIteratorPLIC.H"
#include "reconstructionSchemes.H"
#include "upwind.H"
#include "processorPolyPatch.H"
#include "processorLduInterface.H"
#include "processorBC.H"
#include "pointMesh.H"
#include "pointFields.H"
#include "fixedValuePointPatchFields.H"
#include "wedgePointPatch.H"
#include "processorPointPatch.H"
#include "displacementLaplacianFvMotionSolver.H"
#include "DynamicList.H"
#include "HashSet.H"
#include <algorithm>
#include <vector>

static inline scalar twoPi() { return constant::mathematical::twoPi; }

#include "diagnosticsHelpers.H"
#include "petscBlockAssembly.H"
#include "petscBlockSolve.H"

struct ShapeResidualStats
{
    scalar h;
    scalar rms;
    scalar maxMag;
    scalar area;
};

struct AxisymProfile
{
    List<scalar> theta;
    List<scalar> radius;
    List<scalar> y;
    List<scalar> kappa;
    scalar horizontalRadius;
    scalar verticalRadius;
};

struct ShapeMotionProfile
{
    List<scalar> theta;
    List<scalar> deltaN;
    scalar volumeCorrection;
    scalar minDeltaN;
    scalar maxDeltaN;
    scalar rmsDeltaN;
};

static inline scalar clampScalar
(
    const scalar x,
    const scalar lo,
    const scalar hi
)
{
    return min(max(x, lo), hi);
}

static scalar smoothWeight
(
    const point& p,
    const point& centre,
    const scalar dropRadius,
    const scalar motionRadius
)
{
    const scalar r = mag(p - centre);

    if (r <= dropRadius)
    {
        return 1.0;
    }

    if (r >= motionRadius)
    {
        return 0.0;
    }

    const scalar s = (r - dropRadius)/(motionRadius - dropRadius + VSMALL);
    return 0.5*(1.0 + Foam::cos(constant::mathematical::pi*s));
}

static void collectPatchPointLabels
(
    const polyMesh& mesh,
    const label patchID,
    labelHashSet& pointLabels
)
{
    const polyPatch& pp = mesh.boundaryMesh()[patchID];
    const label start = pp.start();

    forAll(pp, facei)
    {
        const face& f = mesh.faces()[start + facei];
        forAll(f, fp)
        {
            pointLabels.insert(f[fp]);
        }
    }
}

static vector dropForce
(
    const fvMesh& mesh,
    const volScalarField& pr,
    const label patchID,
    const scalar axisymFactor
)
{
    const auto& Sf = mesh.Sf().boundaryField()[patchID];
    const auto& pB = pr.boundaryField()[patchID];

    // Sf on dropWall points from the gas domain into the drop. The radiation
    // pressure traction on the drop follows this area-vector orientation.
    return axisymFactor*gSum(pB*Sf);
}

static point dropPatchCentre
(
    const fvMesh& mesh,
    const label patchID,
    const scalar axisymFactor
)
{
    labelHashSet patchPointLabels;
    collectPatchPointLabels(mesh, patchID, patchPointLabels);

    scalar minY = GREAT;
    scalar maxY = -GREAT;
    label nLocal = 0;

    forAllConstIter(labelHashSet, patchPointLabels, iter)
    {
        const point& p = mesh.points()[iter.key()];
        minY = min(minY, p.y());
        maxY = max(maxY, p.y());
        ++nLocal;
    }

    reduce(minY, minOp<scalar>());
    reduce(maxY, maxOp<scalar>());
    reduce(nLocal, sumOp<label>());

    if (nLocal == 0)
    {
        FatalErrorInFunction
            << "dropWall patch has no points." << exit(FatalError);
    }

    // Axisymmetric profile coordinates use the rotation axis as radial origin.
    return point(0.0, 0.5*(minY + maxY), 0.0);
}

static scalar profileTheta(const point& p, const point& centre)
{
    const scalar radial = Foam::sqrt(sqr(p.x()) + sqr(p.z()));
    return Foam::atan2(radial, p.y() - centre.y());
}

static inline scalar cross2D
(
    const scalar ax,
    const scalar ay,
    const scalar bx,
    const scalar by
)
{
    return ax*by - ay*bx;
}

static void computeProfileCurvature
(
    AxisymProfile& profile,
    const scalar fallbackRadius
)
{
    const label n = profile.theta.size();
    if (n < 3)
    {
        forAll(profile.kappa, profI)
        {
            profile.kappa[profI] = 2.0/(fallbackRadius + VSMALL);
        }
        return;
    }

    forAll(profile.theta, profI)
    {
        const label prevI = max(profI - 1, 0);
        const label nextI = min(profI + 1, n - 1);

        label leftI = prevI;
        label rightI = nextI;

        if (profI == 0)
        {
            leftI = 0;
            rightI = 2;
        }
        else if (profI == n - 1)
        {
            leftI = n - 3;
            rightI = n - 1;
        }

        const scalar r0 = profile.radius[leftI];
        const scalar y0 = profile.y[leftI];
        const scalar r1 = profile.radius[profI];
        const scalar y1 = profile.y[profI];
        const scalar r2 = profile.radius[rightI];
        const scalar y2 = profile.y[rightI];

        const scalar e01r = r1 - r0;
        const scalar e01y = y1 - y0;
        const scalar e12r = r2 - r1;
        const scalar e12y = y2 - y1;
        const scalar e02r = r2 - r0;
        const scalar e02y = y2 - y0;

        const scalar a = Foam::sqrt(sqr(e01r) + sqr(e01y));
        const scalar b = Foam::sqrt(sqr(e12r) + sqr(e12y));
        const scalar c = Foam::sqrt(sqr(e02r) + sqr(e02y));

        scalar meridional = 1.0/(fallbackRadius + VSMALL);
        if (a > SMALL && b > SMALL && c > SMALL)
        {
            // theta increases from the top pole to the bottom pole. The
            // negative sign gives positive curvature for a convex sphere.
            meridional = -2.0*cross2D(e01r, e01y, e12r, e12y)/(a*b*c);
        }

        const scalar tr = e02r;
        const scalar ty = e02y;
        const scalar tMag = Foam::sqrt(sqr(tr) + sqr(ty));

        scalar azimuthal = meridional;
        if (profile.radius[profI] > SMALL && tMag > SMALL)
        {
            const scalar outwardNormalR = -ty/tMag;
            azimuthal = outwardNormalR/(profile.radius[profI] + VSMALL);
        }

        profile.kappa[profI] = meridional + azimuthal;
    }
}

static AxisymProfile buildAxisymProfile
(
    const fvMesh& mesh,
    const label patchID,
    const point& centre,
    const scalar fallbackRadius
)
{
    const auto& Cf = mesh.Cf().boundaryField()[patchID];

    List<List<point>> allProcPoints(Pstream::nProcs());
    allProcPoints[Pstream::myProcNo()] = Cf;

    if (Pstream::parRun())
    {
        Pstream::gatherList(allProcPoints);
        Pstream::scatterList(allProcPoints);
    }

    label nProfile = 0;
    forAll(allProcPoints, proci)
    {
        nProfile += allProcPoints[proci].size();
    }

    List<point> rawProfilePoints(nProfile);
    label count = 0;
    forAll(allProcPoints, proci)
    {
        forAll(allProcPoints[proci], pointi)
        {
            rawProfilePoints[count++] = allProcPoints[proci][pointi];
        }
    }

    std::vector<label> order(nProfile);
    for (label i = 0; i < nProfile; ++i)
    {
        order[i] = i;
    }

    std::sort
    (
        order.begin(),
        order.end(),
        [&](const label a, const label b)
        {
            return profileTheta(rawProfilePoints[a], centre)
                 < profileTheta(rawProfilePoints[b], centre);
        }
    );

    DynamicList<scalar> thetaValues(nProfile);
    DynamicList<scalar> radiusValues(nProfile);
    DynamicList<scalar> yValues(nProfile);

    const scalar thetaTol = 1e-8;
    label i = 0;
    while (i < nProfile)
    {
        const scalar theta0 = profileTheta(rawProfilePoints[order[i]], centre);
        scalar thetaSum = 0.0;
        scalar radiusSum = 0.0;
        scalar ySum = 0.0;
        label nBin = 0;

        while
        (
            i < nProfile
         && mag(profileTheta(rawProfilePoints[order[i]], centre) - theta0)
          <= thetaTol
        )
        {
            const point& p = rawProfilePoints[order[i]];

            thetaSum += profileTheta(p, centre);
            radiusSum += Foam::sqrt(sqr(p.x()) + sqr(p.z()));
            ySum += p.y() - centre.y();
            ++nBin;
            ++i;
        }

        thetaValues.append(thetaSum/nBin);
        radiusValues.append(radiusSum/nBin);
        yValues.append(ySum/nBin);
    }

    nProfile = thetaValues.size();

    AxisymProfile profile;
    profile.theta.setSize(nProfile);
    profile.radius.setSize(nProfile);
    profile.y.setSize(nProfile);
    profile.kappa.setSize(nProfile);
    profile.horizontalRadius = fallbackRadius;
    profile.verticalRadius = fallbackRadius;

    for (label profI = 0; profI < nProfile; ++profI)
    {
        profile.theta[profI] = thetaValues[profI];
        profile.radius[profI] = radiusValues[profI];
        profile.y[profI] = yValues[profI];
        profile.kappa[profI] = 2.0/(fallbackRadius + VSMALL);
    }

    if (nProfile == 0)
    {
        return profile;
    }

    scalar horizontalRadius = 0.0;
    scalar verticalRadius = 0.0;
    forAll(profile.theta, profI)
    {
        horizontalRadius = max(horizontalRadius, profile.radius[profI]);
        verticalRadius = max(verticalRadius, mag(profile.y[profI]));
    }

    profile.horizontalRadius = max(horizontalRadius, fallbackRadius);
    profile.verticalRadius = max(verticalRadius, fallbackRadius);

    computeProfileCurvature(profile, fallbackRadius);

    return profile;
}

static scalar curvatureAt
(
    const AxisymProfile& profile,
    const point& p,
    const point& centre,
    const scalar fallbackRadius
)
{
    const label n = profile.theta.size();
    if (n == 0)
    {
        return 2.0/(fallbackRadius + VSMALL);
    }

    const scalar theta = profileTheta(p, centre);

    label best = 0;
    scalar bestDist = GREAT;
    forAll(profile.theta, i)
    {
        const scalar dist = mag(profile.theta[i] - theta);
        if (dist < bestDist)
        {
            best = i;
            bestDist = dist;
        }
    }

    return profile.kappa[best];
}

static scalar nearestProfileValue
(
    const List<scalar>& thetaValues,
    const List<scalar>& values,
    const scalar theta,
    const scalar fallback
)
{
    if (thetaValues.empty())
    {
        return fallback;
    }

    label best = 0;
    scalar bestDist = GREAT;
    forAll(thetaValues, i)
    {
        const scalar dist = mag(thetaValues[i] - theta);
        if (dist < bestDist)
        {
            best = i;
            bestDist = dist;
        }
    }

    return values[best];
}

static vector axisymOutwardNormal
(
    const point& p,
    const point& centre,
    const AxisymProfile& profile
)
{
    const scalar radial = Foam::sqrt(sqr(p.x()) + sqr(p.z()));

    vector radialHat(vector::zero);
    if (radial > SMALL)
    {
        radialHat = vector(p.x()/radial, 0.0, p.z()/radial);
    }

    const label nProfile = profile.theta.size();
    if (nProfile < 2)
    {
        vector n = p - centre;
        if (mag(n) <= VSMALL)
        {
            n = vector(0, p.y() >= centre.y() ? 1 : -1, 0);
        }
        return n/(mag(n) + VSMALL);
    }

    const scalar theta = profileTheta(p, centre);
    label best = 0;
    scalar bestDist = GREAT;
    forAll(profile.theta, profI)
    {
        const scalar dist = mag(profile.theta[profI] - theta);
        if (dist < bestDist)
        {
            best = profI;
            bestDist = dist;
        }
    }

    label leftI = max(best - 1, 0);
    label rightI = min(best + 1, nProfile - 1);

    if (best == 0)
    {
        rightI = 1;
    }
    else if (best == nProfile - 1)
    {
        leftI = nProfile - 2;
    }

    const scalar tr = profile.radius[rightI] - profile.radius[leftI];
    const scalar ty = profile.y[rightI] - profile.y[leftI];
    const scalar tMag = Foam::sqrt(sqr(tr) + sqr(ty));

    if (tMag <= VSMALL)
    {
        vector n = p - centre;
        if (mag(n) <= VSMALL)
        {
            n = vector(0, p.y() >= centre.y() ? 1 : -1, 0);
        }
        return n/(mag(n) + VSMALL);
    }

    const scalar normalR = -ty/tMag;
    const scalar normalY = tr/tMag;

    vector n
    (
        normalR*radialHat.x(),
        normalY,
        normalR*radialHat.z()
    );

    if (mag(n) <= VSMALL)
    {
        n = vector(0, normalY >= 0 ? 1 : -1, 0);
    }

    return n/(mag(n) + VSMALL);
}

static wordList pointDisplacementPatchFieldTypes(const fvMesh& mesh)
{
    const pointBoundaryMesh& pbm = pointMesh::New(mesh).boundary();

    wordList patchFieldTypes(pbm.size());
    forAll(pbm, patchi)
    {
        const word patchType = pbm[patchi].type();
        if
        (
            patchType == wedgePointPatch::typeName
         || patchType == processorPointPatch::typeName
        )
        {
            patchFieldTypes[patchi] = patchType;
        }
        else
        {
            patchFieldTypes[patchi] =
                fixedValuePointPatchField<vector>::typeName;
        }
    }

    return patchFieldTypes;
}

static wordList pointDisplacementPatchTypes(const fvMesh& mesh)
{
    const pointBoundaryMesh& pbm = pointMesh::New(mesh).boundary();

    wordList patchTypes(pbm.size());
    forAll(pbm, patchi)
    {
        patchTypes[patchi] = pbm[patchi].type();
    }

    return patchTypes;
}

static scalar axisymVolumeFromPatchPoints
(
    const fvMesh& mesh,
    const label patchID,
    const pointField& points
)
{
    labelHashSet patchPointLabels;
    collectPatchPointLabels(mesh, patchID, patchPointLabels);

    DynamicList<vector> localEntries(patchPointLabels.size());

    forAllConstIter(labelHashSet, patchPointLabels, iter)
    {
        const point& p = points[iter.key()];
        localEntries.append
        (
            vector(p.y(), Foam::sqrt(sqr(p.x()) + sqr(p.z())), 0.0)
        );
    }

    List<vector> localList(localEntries);
    List<List<vector>> allProcEntries(Pstream::nProcs());
    allProcEntries[Pstream::myProcNo()] = localList;

    if (Pstream::parRun())
    {
        Pstream::gatherList(allProcEntries);
        Pstream::scatterList(allProcEntries);
    }

    label nEntries = 0;
    forAll(allProcEntries, proci)
    {
        nEntries += allProcEntries[proci].size();
    }

    if (nEntries < 2)
    {
        return 0.0;
    }

    List<vector> entries(nEntries);
    label count = 0;
    forAll(allProcEntries, proci)
    {
        forAll(allProcEntries[proci], entryi)
        {
            entries[count++] = allProcEntries[proci][entryi];
        }
    }

    std::vector<label> order(nEntries);
    for (label i = 0; i < nEntries; ++i)
    {
        order[i] = i;
    }

    std::sort
    (
        order.begin(),
        order.end(),
        [&](const label a, const label b)
        {
            return entries[a].x() < entries[b].x();
        }
    );

    const scalar minY = entries[order.front()].x();
    const scalar maxY = entries[order.back()].x();
    const scalar yTol = max(1e-12, 1e-8*(maxY - minY));

    DynamicList<vector> profile(nEntries);
    label i = 0;
    while (i < nEntries)
    {
        const scalar y0 = entries[order[i]].x();
        scalar ySum = 0.0;
        scalar maxR = 0.0;
        label nBin = 0;

        while
        (
            i < nEntries
         && mag(entries[order[i]].x() - y0) <= yTol
        )
        {
            ySum += entries[order[i]].x();
            maxR = max(maxR, entries[order[i]].y());
            ++nBin;
            ++i;
        }

        profile.append(vector(ySum/nBin, maxR, 0.0));
    }

    scalar volume = 0.0;
    for (label profI = 1; profI < profile.size(); ++profI)
    {
        const scalar y0 = profile[profI - 1].x();
        const scalar y1 = profile[profI].x();
        const scalar r0 = profile[profI - 1].y();
        const scalar r1 = profile[profI].y();
        const scalar dy = max(y1 - y0, 0.0);

        volume += constant::mathematical::pi*dy
            *(sqr(r0) + r0*r1 + sqr(r1))/3.0;
    }

    return volume;
}

static ShapeMotionProfile buildShapeMotionProfile
(
    const fvMesh& mesh,
    const volScalarField& pr,
    const label patchID,
    const scalar axisymFactor,
    const scalar liquidDensity,
    const scalar gravityMagnitude,
    const scalar surfaceTension,
    const AxisymProfile& axisymProfile,
    const point& centre,
    const scalar fallbackRadius,
    const ShapeResidualStats& shapeStats,
    const scalar shapeRelaxation,
    const scalar maxShapeDelta
)
{
    const auto& Cf = mesh.Cf().boundaryField()[patchID];
    const auto& magSf = mesh.magSf().boundaryField()[patchID];
    const auto& pB = pr.boundaryField()[patchID];

    DynamicList<vector> localEntries(pB.size());
    scalar area = 0.0;
    scalar dV = 0.0;

    forAll(pB, facei)
    {
        const scalar kappa =
            curvatureAt(axisymProfile, Cf[facei], centre, fallbackRadius);
        const scalar q =
            pB[facei]
          + liquidDensity*gravityMagnitude*Cf[facei].y()
          + surfaceTension*kappa;
        const scalar residual = q - shapeStats.h;

        const scalar rawDeltaN =
            clampScalar(-shapeRelaxation*residual, -maxShapeDelta, maxShapeDelta);

        const scalar dA = axisymFactor*magSf[facei];
        area += dA;
        dV += rawDeltaN*dA;

        localEntries.append(vector(profileTheta(Cf[facei], centre), rawDeltaN, 0));
    }

    reduce(area, sumOp<scalar>());
    reduce(dV, sumOp<scalar>());

    const scalar volumeCorrection = -dV/(area + VSMALL);

    scalar minDeltaN = VGREAT;
    scalar maxDeltaN = -VGREAT;
    scalar sumSqrDeltaN = 0.0;
    label nDeltaN = 0;

    forAll(localEntries, i)
    {
        localEntries[i].y() =
            clampScalar
            (
                localEntries[i].y() + volumeCorrection,
               -maxShapeDelta,
                maxShapeDelta
            );

        minDeltaN = min(minDeltaN, localEntries[i].y());
        maxDeltaN = max(maxDeltaN, localEntries[i].y());
        sumSqrDeltaN += sqr(localEntries[i].y());
        ++nDeltaN;
    }

    reduce(minDeltaN, minOp<scalar>());
    reduce(maxDeltaN, maxOp<scalar>());
    reduce(sumSqrDeltaN, sumOp<scalar>());
    reduce(nDeltaN, sumOp<label>());

    if (nDeltaN == 0)
    {
        minDeltaN = 0.0;
        maxDeltaN = 0.0;
    }

    List<vector> localList(localEntries);
    List<List<vector>> allProcEntries(Pstream::nProcs());
    allProcEntries[Pstream::myProcNo()] = localList;

    if (Pstream::parRun())
    {
        Pstream::gatherList(allProcEntries);
        Pstream::scatterList(allProcEntries);
    }

    label nEntries = 0;
    forAll(allProcEntries, proci)
    {
        nEntries += allProcEntries[proci].size();
    }

    List<vector> rawEntries(nEntries);
    label count = 0;
    forAll(allProcEntries, proci)
    {
        forAll(allProcEntries[proci], entryi)
        {
            rawEntries[count++] = allProcEntries[proci][entryi];
        }
    }

    std::vector<label> order(nEntries);
    for (label i = 0; i < nEntries; ++i)
    {
        order[i] = i;
    }

    std::sort
    (
        order.begin(),
        order.end(),
        [&](const label a, const label b)
        {
            return rawEntries[a].x() < rawEntries[b].x();
        }
    );

    DynamicList<scalar> thetaValues(nEntries);
    DynamicList<scalar> deltaNValues(nEntries);

    const scalar thetaTol = 1e-8;
    label i = 0;
    while (i < nEntries)
    {
        const scalar theta0 = rawEntries[order[i]].x();
        scalar thetaSum = 0.0;
        scalar deltaNSum = 0.0;
        label nBin = 0;

        while
        (
            i < nEntries
         && mag(rawEntries[order[i]].x() - theta0) <= thetaTol
        )
        {
            thetaSum += rawEntries[order[i]].x();
            deltaNSum += rawEntries[order[i]].y();
            ++nBin;
            ++i;
        }

        thetaValues.append(thetaSum/nBin);
        deltaNValues.append(deltaNSum/nBin);
    }

    ShapeMotionProfile profile;
    profile.theta.setSize(thetaValues.size());
    profile.deltaN.setSize(deltaNValues.size());
    profile.volumeCorrection = volumeCorrection;
    profile.minDeltaN = minDeltaN;
    profile.maxDeltaN = maxDeltaN;
    profile.rmsDeltaN = Foam::sqrt(sumSqrDeltaN/(scalar(nDeltaN) + VSMALL));

    forAll(profile.theta, profI)
    {
        profile.theta[profI] = thetaValues[profI];
        profile.deltaN[profI] = deltaNValues[profI];
    }

    return profile;
}

static void writeEquilibriumHistory
(
    const Time& runTime,
    const label iter,
    const label shapeIter,
    const label positionIter,
    const word& phase,
    const scalar centreY,
    const vector& force,
    const scalar weight,
    const scalar residual,
    const scalar deltaY,
    const ShapeResidualStats& shapeStats,
    const scalar shapeVolumeCorrection,
    const scalar shapeDeltaNMin,
    const scalar shapeDeltaNMax,
    const scalar shapeDeltaNRms
)
{
    if (!Pstream::master())
    {
        return;
    }

    const fileName outDir
    (
        runTime.globalPath()/"postProcessing"/"dropEquilibrium"
    );
    mkDir(outDir);

    const fileName outFile(outDir/"history.dat");
    const bool writeHeader = !exists(outFile);

    OFstream os(outFile, IOstreamOption(), IOstreamOption::APPEND);
    if (writeHeader)
    {
        os  << "# iter"
            << tab << "shapeIter"
            << tab << "positionIter"
            << tab << "phase"
            << tab << "time"
            << tab << "centerY"
            << tab << "Fx"
            << tab << "Fy"
            << tab << "Fz"
            << tab << "weight"
            << tab << "forceResidual"
            << tab << "deltaY"
            << tab << "shapeResidualRms"
            << tab << "shapeResidualMax"
            << tab << "laplaceH"
            << tab << "dropArea"
            << tab << "shapeVolumeCorrection"
            << tab << "shapeDeltaNMin"
            << tab << "shapeDeltaNMax"
            << tab << "shapeDeltaNRms"
            << nl;
    }

    os  << iter
        << tab << shapeIter
        << tab << positionIter
        << tab << phase
        << tab << runTime.timeName()
        << tab << centreY
        << tab << force.x()
        << tab << force.y()
        << tab << force.z()
        << tab << weight
        << tab << residual
        << tab << deltaY
        << tab << shapeStats.rms
        << tab << shapeStats.maxMag
        << tab << shapeStats.h
        << tab << shapeStats.area
        << tab << shapeVolumeCorrection
        << tab << shapeDeltaNMin
        << tab << shapeDeltaNMax
        << tab << shapeDeltaNRms
        << nl;
}

static ShapeResidualStats computeShapeResidualStats
(
    const fvMesh& mesh,
    const volScalarField& pr,
    const label patchID,
    const scalar axisymFactor,
    const scalar liquidDensity,
    const scalar gravityMagnitude,
    const scalar surfaceTension,
    const AxisymProfile& profile,
    const point& centre,
    const scalar fallbackRadius
)
{
    const auto& Cf = mesh.Cf().boundaryField()[patchID];
    const auto& magSf = mesh.magSf().boundaryField()[patchID];
    const auto& pB = pr.boundaryField()[patchID];

    scalar area = 0.0;
    scalar qIntegral = 0.0;

    forAll(pB, facei)
    {
        const scalar kappa =
            curvatureAt(profile, Cf[facei], centre, fallbackRadius);
        const scalar dA = axisymFactor*magSf[facei];
        const scalar q =
            pB[facei]
          + liquidDensity*gravityMagnitude*Cf[facei].y()
          + surfaceTension*kappa;

        area += dA;
        qIntegral += dA*q;
    }

    reduce(area, sumOp<scalar>());
    reduce(qIntegral, sumOp<scalar>());

    const scalar h = qIntegral/(area + VSMALL);

    scalar rmsIntegral = 0.0;
    scalar maxMagResidual = 0.0;

    forAll(pB, facei)
    {
        const scalar kappa =
            curvatureAt(profile, Cf[facei], centre, fallbackRadius);
        const scalar dA = axisymFactor*magSf[facei];
        const scalar q =
            pB[facei]
          + liquidDensity*gravityMagnitude*Cf[facei].y()
          + surfaceTension*kappa;
        const scalar residual = q - h;

        rmsIntegral += dA*sqr(residual);
        maxMagResidual = max(maxMagResidual, mag(residual));
    }

    reduce(rmsIntegral, sumOp<scalar>());
    reduce(maxMagResidual, maxOp<scalar>());

    ShapeResidualStats stats;
    stats.h = h;
    stats.rms = Foam::sqrt(rmsIntegral/(area + VSMALL));
    stats.maxMag = maxMagResidual;
    stats.area = area;
    return stats;
}

static void writeShapeResidualProfile
(
    const Time& runTime,
    const label iter,
    const fvMesh& mesh,
    const volScalarField& pr,
    const label patchID,
    const point& centre,
    const scalar liquidDensity,
    const scalar gravityMagnitude,
    const scalar surfaceTension,
    const AxisymProfile& profile,
    const scalar fallbackRadius,
    const ShapeResidualStats& shapeStats
)
{
    const auto& Cf = mesh.Cf().boundaryField()[patchID];
    const auto& magSf = mesh.magSf().boundaryField()[patchID];
    const auto& pB = pr.boundaryField()[patchID];

    const fileName outDir
    (
        runTime.globalPath()/"postProcessing"/"dropEquilibrium"
    );
    mkDir(outDir);

    const fileName outFile
    (
        outDir/
        (
            "profile_"
          + Foam::name(iter)
          + "_proc"
          + Foam::name(Pstream::myProcNo())
          + ".dat"
        )
    );

    OFstream os(outFile);
    os  << "# theta"
        << tab << "radius"
        << tab << "x"
        << tab << "y"
        << tab << "z"
        << tab << "pr"
        << tab << "kappa"
        << tab << "q"
        << tab << "h"
        << tab << "residual"
        << tab << "area"
        << nl;

    forAll(pB, facei)
    {
        const vector rel = Cf[facei] - centre;
        const scalar radial = Foam::sqrt(sqr(rel.x()) + sqr(rel.z()));
        const scalar theta = Foam::atan2(radial, rel.y());
        const scalar radius = mag(rel);
        const scalar kappa =
            curvatureAt(profile, Cf[facei], centre, fallbackRadius);
        const scalar q =
            pB[facei]
          + liquidDensity*gravityMagnitude*Cf[facei].y()
          + surfaceTension*kappa;
        const scalar residual = q - shapeStats.h;

        os  << theta
            << tab << radius
            << tab << Cf[facei].x()
            << tab << Cf[facei].y()
            << tab << Cf[facei].z()
            << tab << pB[facei]
            << tab << kappa
            << tab << q
            << tab << shapeStats.h
            << tab << residual
            << tab << magSf[facei]
            << nl;
    }
}

static void moveDropVertically
(
    fvMesh& mesh,
    const label dropPatchID,
    const point& centre,
    const scalar dropRadius,
    const scalar motionRadius,
    const scalar deltaY,
    const bool enableShapeUpdate,
    const scalar targetDropVolume,
    const AxisymProfile& axisymProfile,
    const ShapeMotionProfile& shapeMotionProfile
)
{
    pointField newPoints(mesh.points());

    labelHashSet dropPatchPoints;
    collectPatchPointLabels(mesh, dropPatchID, dropPatchPoints);

    const polyBoundaryMesh& patches = mesh.boundaryMesh();
    labelHashSet fixedBoundaryPoints;
    labelHashSet axisBoundaryPoints;

    forAll(patches, patchi)
    {
        const polyPatch& pp = patches[patchi];

        if (pp.name() == "axis")
        {
            const label start = pp.start();
            forAll(pp, facei)
            {
                const face& f = mesh.faces()[start + facei];
                forAll(f, fp)
                {
                    axisBoundaryPoints.insert(f[fp]);
                }
            }
            break;
        }
    }

    forAll(patches, patchi)
    {
        if
        (
            patchi == dropPatchID
         || patches[patchi].coupled()
         || isA<wedgePolyPatch>(patches[patchi])
         || isA<emptyPolyPatch>(patches[patchi])
        )
        {
            continue;
        }

        const polyPatch& pp = patches[patchi];
        const label start = pp.start();

        if (pp.name() == "axis")
        {
            continue;
        }

        const bool anchorsAxis =
            pp.name() == "transducer1"
         || pp.name() == "transducer2"
         || pp.name() == "reflector1"
         || pp.name() == "reflector2";

        forAll(pp, facei)
        {
            const face& f = mesh.faces()[start + facei];
            forAll(f, fp)
            {
                if (!axisBoundaryPoints.found(f[fp]) || anchorsAxis)
                {
                    fixedBoundaryPoints.insert(f[fp]);
                }
            }
        }
    }

    pointField prescribedDisplacement(newPoints.size(), vector::zero);

    forAll(newPoints, pointi)
    {
        scalar w = 0.0;

        if (dropPatchPoints.found(pointi))
        {
            w = 1.0;
        }
        else if (!fixedBoundaryPoints.found(pointi))
        {
            w = smoothWeight(newPoints[pointi], centre, dropRadius, motionRadius);
        }

        vector displacement(0, deltaY, 0);

        if (enableShapeUpdate && w > VSMALL)
        {
            const scalar theta = profileTheta(newPoints[pointi], centre);
            const scalar deltaN =
                nearestProfileValue
                (
                    shapeMotionProfile.theta,
                    shapeMotionProfile.deltaN,
                    theta,
                    0.0
                );
            displacement += deltaN*axisymOutwardNormal
            (
                newPoints[pointi],
                centre,
                axisymProfile
            );
        }

        if (axisBoundaryPoints.found(pointi))
        {
            displacement = vector(0, displacement.y(), 0);
        }

        prescribedDisplacement[pointi] = w*displacement;
    }

    if (enableShapeUpdate && targetDropVolume > VSMALL)
    {
        pointField correctedPoints(newPoints);
        forAll(correctedPoints, pointi)
        {
            correctedPoints[pointi] += prescribedDisplacement[pointi];
        }

        const scalar uncorrectedVolume =
            axisymVolumeFromPatchPoints(mesh, dropPatchID, correctedPoints);

        const scalar volumeError = uncorrectedVolume - targetDropVolume;

        if (mag(volumeError) > SMALL*targetDropVolume)
        {
            const scalar correctionLimit = max(0.1*dropRadius, 10.0*SMALL);
            scalar lo = -correctionLimit;
            scalar hi = correctionLimit;

            auto volumeWithOffset = [&](const scalar offset)
            {
                pointField candidatePoints(newPoints);
                forAll(candidatePoints, pointi)
                {
                    candidatePoints[pointi] += prescribedDisplacement[pointi];

                    if (dropPatchPoints.found(pointi))
                    {
                        vector correction =
                            offset*axisymOutwardNormal
                            (
                                newPoints[pointi],
                                centre,
                                axisymProfile
                            );

                        if (axisBoundaryPoints.found(pointi))
                        {
                            correction = vector(0.0, correction.y(), 0.0);
                        }

                        candidatePoints[pointi] += correction;
                    }
                }

                return axisymVolumeFromPatchPoints
                (
                    mesh,
                    dropPatchID,
                    candidatePoints
                );
            };

            scalar fLo = volumeWithOffset(lo) - targetDropVolume;
            scalar fHi = volumeWithOffset(hi) - targetDropVolume;

            if (fLo*fHi <= 0.0)
            {
                for (label iter = 0; iter < 40; ++iter)
                {
                    const scalar mid = 0.5*(lo + hi);
                    const scalar fMid = volumeWithOffset(mid) - targetDropVolume;

                    if (fLo*fMid <= 0.0)
                    {
                        hi = mid;
                        fHi = fMid;
                    }
                    else
                    {
                        lo = mid;
                        fLo = fMid;
                    }
                }

                const scalar normalOffset = 0.5*(lo + hi);
                forAll(prescribedDisplacement, pointi)
                {
                    if (dropPatchPoints.found(pointi))
                    {
                        vector correction =
                            normalOffset*axisymOutwardNormal
                            (
                                newPoints[pointi],
                                centre,
                                axisymProfile
                            );

                        if (axisBoundaryPoints.found(pointi))
                        {
                            correction = vector(0.0, correction.y(), 0.0);
                        }

                        prescribedDisplacement[pointi] += correction;
                    }
                }

                Info<< "Exact volume projection: uncorrectedVolume="
                    << uncorrectedVolume
                    << " targetVolume=" << targetDropVolume
                    << " volumeError=" << volumeError
                    << " normalOffset=" << normalOffset
                    << nl;
            }
            else
            {
                WarningInFunction
                    << "Could not bracket volume correction. "
                    << "uncorrectedVolume=" << uncorrectedVolume
                    << " targetVolume=" << targetDropVolume
                    << " volumeError=" << volumeError
                    << " fLo=" << fLo
                    << " fHi=" << fHi
                    << nl;
            }
        }
    }

    IOdictionary motionDict
    (
        IOobject
        (
            "dynamicMeshDict",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );

    pointIOField points0
    (
        IOobject
        (
            "points0",
            mesh.time().timeName(),
            polyMesh::meshSubDir,
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        ),
        mesh.points()
    );

    pointVectorField pointDisplacement
    (
        IOobject
        (
            "pointDisplacement",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        pointMesh::New(mesh),
        dimensionedVector(dimLength, vector::zero),
        pointDisplacementPatchFieldTypes(mesh),
        pointDisplacementPatchTypes(mesh)
    );

    pointDisplacement.primitiveFieldRef() = prescribedDisplacement;

    forAll(pointDisplacement.boundaryField(), patchi)
    {
        pointPatchField<vector>& ppf =
            pointDisplacement.boundaryFieldRef()[patchi];

        if (isA<fixedValuePointPatchField<vector>>(ppf))
        {
            fixedValuePointPatchField<vector>& fixedPpf =
                refCast<fixedValuePointPatchField<vector>>(ppf);
            const labelUList& meshPoints = ppf.patch().meshPoints();
            forAll(meshPoints, pointPatchI)
            {
                fixedPpf[pointPatchI] =
                    prescribedDisplacement[meshPoints[pointPatchI]];
            }
        }
    }

    pointDisplacement.correctBoundaryConditions();

    displacementLaplacianFvMotionSolver motionSolver
    (
        mesh,
        motionDict,
        pointDisplacement,
        points0
    );

    motionSolver.solve();
    mesh.movePoints(motionSolver.curPoints());
}

struct DropEquilibriumState
{
    vector force;
    scalar forceResidual;
    point centre;
    AxisymProfile axisymProfile;
    ShapeResidualStats shapeStats;

    DropEquilibriumState()
    :
        force(vector::zero),
        forceResidual(GREAT),
        centre(vector::zero),
        axisymProfile(),
        shapeStats()
    {
        shapeStats.h = 0.0;
        shapeStats.rms = GREAT;
        shapeStats.maxMag = GREAT;
        shapeStats.area = 0.0;
    }
};

static void solveAcousticAndRadiationFields
(
    simpleControl& simple,
    const globalIndex& globalCells,
    const PetscInt N,
    Mat& M,
    Vec& x,
    Vec& b,
    KSP& ksp,
    bool& dumpedMatrixStats,
    bool& dumpedOpStats,
    volScalarField& rho,
    const surfaceScalarField& invRhof,
    const volScalarField& k2,
    volScalarField& Pre,
    volScalarField& Pim,
    volScalarField& SC0,
    volScalarField& SC1,
    volTensorField& T0,
    volTensorField& T1,
    volVectorField& Ure,
    volVectorField& Uim,
    volScalarField& pa,
    volScalarField& pr,
    volTensorField& momFlux,
    const volScalarField& alpha1,
    const dimensionedScalar& f,
    const dimensionedScalar& kl,
    const dimensionedScalar& kg
)
{
    while (simple.correctNonOrthogonal())
    {
        MatZeroEntries(M);
        VecSet(b, 0.0);

        fvScalarMatrix AopPre
        (
          rho*fvm::laplacian(invRhof, Pre)
          + fvm::laplacian(T0, Pre)
          + fvm::Sp(k2 - SC0, Pre)
        );

        fvScalarMatrix AopPim
        (
          rho*fvm::laplacian(invRhof, Pim)
          + fvm::laplacian(T0, Pim)
          + fvm::Sp(k2 - SC0, Pim)
        );

        // Keep B1 (laplacian) and B2 (Sp) as separate operators.
        // Their off-block diagonal handling differs in assembly.
        fvScalarMatrix couplingLaplPre(fvm::laplacian(T1, Pre));  // B1
        fvScalarMatrix couplingMassPre(fvm::Sp(SC1, Pre));         // B2

        fvScalarMatrix couplingLaplPim(fvm::laplacian(T1, Pim));  // B1
        fvScalarMatrix couplingMassPim(fvm::Sp(SC1, Pim));         // B2

        if (!dumpedOpStats)
        {
            reportFvMatrixBreakdown("AopPre beforeBoundaryManipulate", AopPre);
            reportFvMatrixBreakdown("AopPim beforeBoundaryManipulate", AopPim);
            reportFvMatrixBreakdown("couplingLaplPre beforeBoundaryManipulate", couplingLaplPre);
            reportFvMatrixBreakdown("couplingMassPre beforeBoundaryManipulate", couplingMassPre);
        }

        assembleBlockSystem
        (
            M, globalCells, N,
            AopPim, AopPre,
            couplingLaplPre, couplingMassPre,
            couplingLaplPim, couplingMassPim
        );

        if (!dumpedOpStats)
        {
            reportFvMatrixBreakdown("AopPre afterBoundaryManipulate", AopPre);
            reportFvMatrixBreakdown("AopPim afterBoundaryManipulate", AopPim);
            reportFvMatrixBreakdown("couplingLaplPre afterBoundaryManipulate", couplingLaplPre);
            reportFvMatrixBreakdown("couplingMassPre afterBoundaryManipulate", couplingMassPre);
            dumpedOpStats = true;
        }

        scalarField bPim;
        scalarField bPre;
        buildRhs
        (
            AopPim,
            AopPre,
            couplingLaplPre,
            couplingMassPre,
            couplingLaplPim,
            couplingMassPim,
            bPim,
            bPre
        );

        setBlockRhs(b, globalCells, N, bPim, bPre);

        MatAssemblyBegin(M, MAT_FINAL_ASSEMBLY);
        MatAssemblyEnd(M, MAT_FINAL_ASSEMBLY);
        VecAssemblyBegin(b);
        VecAssemblyEnd(b);

        if (!dumpedMatrixStats)
        {
            dumpedMatrixStats = true;
            reportMatrixStats(M, b, N);
            dumpProcessorInterfaceRows(M, AopPim, globalCells, N);
        }

        KSPSolve(ksp, b, x);
        reportKspStats(ksp, b);
        scatterBlockSolution(x, globalCells, N, Pim, Pre);

        Pre.correctBoundaryConditions();
        Pim.correctBoundaryConditions();
    }

    Ure == 1/(2*constant::mathematical::pi*f*rho) * fvc::grad(Pim);
    Uim == -1/(2*constant::mathematical::pi*f*rho) * fvc::grad(Pre);
    pa == Foam::sqrt(Pim*Pim + Pre*Pre);
    pr == 0.25*(kl*alpha1 + kg*(1-alpha1))*(Pre*Pre + Pim*Pim)
        - 0.25*rho*((Ure&Ure) + (Uim&Uim));
    momFlux == 0.5*rho*(Ure*Ure + Uim*Uim);
}

static DropEquilibriumState evaluateDropEquilibrium
(
    fvMesh& mesh,
    const volScalarField& pr,
    const label dropPatchID,
    const scalar axisymFactor,
    const scalar weight,
    const scalar liquidDensity,
    const scalar gravityMagnitude,
    const scalar surfaceTension,
    const scalar curvatureFallbackRadius
)
{
    DropEquilibriumState state;

    state.force = dropForce(mesh, pr, dropPatchID, axisymFactor);
    state.forceResidual = state.force.y() - weight;
    state.centre = dropPatchCentre(mesh, dropPatchID, axisymFactor);
    state.axisymProfile =
        buildAxisymProfile
        (
            mesh,
            dropPatchID,
            state.centre,
            curvatureFallbackRadius
        );
    state.shapeStats =
        computeShapeResidualStats
        (
            mesh,
            pr,
            dropPatchID,
            axisymFactor,
            liquidDensity,
            gravityMagnitude,
            surfaceTension,
            state.axisymProfile,
            state.centre,
            curvatureFallbackRadius
        );

    return state;
}

static bool runFixedShapePositionUpdate
(
    simpleControl& simple,
    Time& runTime,
    fvMesh& mesh,
    const globalIndex& globalCells,
    const PetscInt N,
    Mat& M,
    Vec& x,
    Vec& b,
    KSP& ksp,
    bool& dumpedMatrixStats,
    bool& dumpedOpStats,
    volScalarField& rho,
    const surfaceScalarField& invRhof,
    const volScalarField& k2,
    volScalarField& Pre,
    volScalarField& Pim,
    volScalarField& SC0,
    volScalarField& SC1,
    volTensorField& T0,
    volTensorField& T1,
    volVectorField& Ure,
    volVectorField& Uim,
    volScalarField& pa,
    volScalarField& pr,
    volTensorField& momFlux,
    const volScalarField& alpha1,
    const dimensionedScalar& f,
    const dimensionedScalar& kl,
    const dimensionedScalar& kg,
    const label dropPatchID,
    const scalar axisymFactor,
    const scalar weight,
    const scalar liquidDensity,
    const scalar gravityMagnitude,
    const scalar surfaceTension,
    const scalar curvatureFallbackRadius,
    const scalar dropRadius,
    const scalar motionRadius,
    const scalar displacementRelaxation,
    const scalar maxDeltaY,
    const Switch enablePositionUpdate,
    const scalar forceTolerance,
    const label maxPositionIter,
    const label maxEquilibriumIter,
    const label shapeIter,
    label& eqIter,
    const Switch writeShapeProfiles,
    const scalar maxAllowedShapeResidual,
    const scalar maxAllowedForceMagnitude,
    DropEquilibriumState& state
)
{
    for
    (
        label positionIter = 0;
        positionIter < maxPositionIter && eqIter < maxEquilibriumIter;
        ++positionIter
    )
    {
        if (!simple.loop())
        {
            return false;
        }

        Info<< "Time = " << runTime.timeName()
            << " shapeIter=" << shapeIter
            << " positionIter=" << positionIter
            << nl << endl;

        solveAcousticAndRadiationFields
        (
            simple,
            globalCells,
            N,
            M,
            x,
            b,
            ksp,
            dumpedMatrixStats,
            dumpedOpStats,
            rho,
            invRhof,
            k2,
            Pre,
            Pim,
            SC0,
            SC1,
            T0,
            T1,
            Ure,
            Uim,
            pa,
            pr,
            momFlux,
            alpha1,
            f,
            kl,
            kg
        );

        state =
            evaluateDropEquilibrium
            (
                mesh,
                pr,
                dropPatchID,
                axisymFactor,
                weight,
                liquidDensity,
                gravityMagnitude,
                surfaceTension,
                curvatureFallbackRadius
            );

        scalar deltaY =
            displacementRelaxation*dropRadius*state.forceResidual/(weight + VSMALL);
        deltaY = clampScalar(deltaY, -maxDeltaY, maxDeltaY);

        writeEquilibriumHistory
        (
            runTime,
            eqIter,
            shapeIter,
            positionIter,
            word("position"),
            state.centre.y(),
            state.force,
            weight,
            state.forceResidual,
            deltaY,
            state.shapeStats,
            0.0,
            0.0,
            0.0,
            0.0
        );

        if (writeShapeProfiles)
        {
            writeShapeResidualProfile
            (
                runTime,
                eqIter,
                mesh,
                pr,
                dropPatchID,
                state.centre,
                liquidDensity,
                gravityMagnitude,
                surfaceTension,
                state.axisymProfile,
                curvatureFallbackRadius,
                state.shapeStats
            );
        }

        Info<< "Position iter " << positionIter
            << " in shapeIter " << shapeIter
            << ": centerY=" << state.centre.y()
            << " Fy=" << state.force.y()
            << " weight=" << weight
            << " residual=" << state.forceResidual
            << " deltaY=" << deltaY
            << " shapeResidualRms=" << state.shapeStats.rms
            << nl;

        runTime.write();
        ++eqIter;

        if
        (
            state.shapeStats.rms > maxAllowedShapeResidual
         || mag(state.force) > maxAllowedForceMagnitude
        )
        {
            WarningInFunction
                << "Stopping drop equilibrium loop because stability limits "
                << "were exceeded: shapeResidualRms=" << state.shapeStats.rms
                << " maxAllowedShapeResidual=" << maxAllowedShapeResidual
                << " |force|=" << mag(state.force)
                << " maxAllowedForceMagnitude=" << maxAllowedForceMagnitude
                << nl;
            return false;
        }

        if (!enablePositionUpdate)
        {
            Info<< "Fixed-shape vertical position update is disabled." << nl;
            return true;
        }

        if (mag(state.forceResidual) <= forceTolerance)
        {
            Info<< "Fixed-shape vertical equilibrium reached." << nl;
            return true;
        }

        moveDropVertically
        (
            mesh,
            dropPatchID,
            state.centre,
            dropRadius,
            motionRadius,
            deltaY,
            false,
            0.0,
            state.axisymProfile,
            ShapeMotionProfile()
        );

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    WarningInFunction
        << "Fixed-shape position equilibrium was not reached within "
        << "maxPositionIter=" << maxPositionIter
        << " or maxEquilibriumIter=" << maxEquilibriumIter << nl;

    return false;
}

static void applyShapeUpdate
(
    Time& runTime,
    fvMesh& mesh,
    const volScalarField& pr,
    const label dropPatchID,
    const scalar axisymFactor,
    const scalar weight,
    const scalar liquidDensity,
    const scalar gravityMagnitude,
    const scalar surfaceTension,
    const scalar curvatureFallbackRadius,
    const scalar dropRadius,
    const scalar dropVolume,
    const scalar motionRadius,
    const scalar shapeRelaxation,
    const scalar maxShapeDelta,
    const label eqIter,
    const label shapeIter,
    const DropEquilibriumState& state
)
{
    ShapeMotionProfile shapeMotionProfile =
        buildShapeMotionProfile
        (
            mesh,
            pr,
            dropPatchID,
            axisymFactor,
            liquidDensity,
            gravityMagnitude,
            surfaceTension,
            state.axisymProfile,
            state.centre,
            curvatureFallbackRadius,
            state.shapeStats,
            shapeRelaxation,
            maxShapeDelta
        );

    writeEquilibriumHistory
    (
        runTime,
        eqIter,
        shapeIter,
        -1,
        word("shapeUpdate"),
        state.centre.y(),
        state.force,
        weight,
        state.forceResidual,
        0.0,
        state.shapeStats,
        shapeMotionProfile.volumeCorrection,
        shapeMotionProfile.minDeltaN,
        shapeMotionProfile.maxDeltaN,
        shapeMotionProfile.rmsDeltaN
    );

    Info<< "Shape update " << shapeIter
        << ": shapeResidualRms=" << state.shapeStats.rms
        << " shapeVolumeCorrection=" << shapeMotionProfile.volumeCorrection
        << " shapeDeltaN[min,max,rms]=("
        << shapeMotionProfile.minDeltaN
        << ", "
        << shapeMotionProfile.maxDeltaN
        << ", "
        << shapeMotionProfile.rmsDeltaN
        << ")"
        << nl;

    moveDropVertically
    (
        mesh,
        dropPatchID,
        state.centre,
        dropRadius,
        motionRadius,
        0.0,
        true,
        dropVolume,
        state.axisymProfile,
        shapeMotionProfile
    );
}

int main(int argc, char *argv[])
{
    #include "postProcess.H"
    #include "addCheckCaseOptions.H"
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createFields.H"
    #include "computePMLCoefs.H"
    #include "computeAlphaf.H"

    simpleControl simple(mesh);

    PetscInitialize(&argc, &argv, nullptr, nullptr);

    // Global indexing for block system
    globalIndex globalCells(mesh.nCells());
    const PetscInt N      = (PetscInt)globalCells.size();
    const PetscInt nLocal = (PetscInt)mesh.nCells();

    Mat M;
    Vec x, b;
    KSP ksp;
    bool dumpedMatrixStats = false;
    bool dumpedOpStats = false;

    initializePetscSystem(nLocal, N, M, x, b, ksp);

    IOdictionary dropEquilibriumDict
    (
        IOobject
        (
            "dropEquilibriumDict",
            runTime.system(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        )
    );

    const word dropPatchName =
        dropEquilibriumDict.lookupOrDefault<word>("dropPatch", "dropWall");
    const label dropPatchID = mesh.boundaryMesh().findPatchID(dropPatchName);
    if (dropPatchID < 0)
    {
        FatalErrorInFunction
            << "Patch '" << dropPatchName << "' not found."
            << exit(FatalError);
    }

    const scalar axisymFactor =
        dropEquilibriumDict.lookupOrDefault<scalar>("axisymFactor", 360.0);
    const scalar dropRadius =
        dropEquilibriumDict.lookupOrDefault<scalar>("dropRadius", 0.001);
    const scalar dropVolume =
        dropEquilibriumDict.lookupOrDefault<scalar>
        (
            "dropVolume",
            4.0*constant::mathematical::pi*pow3(dropRadius)/3.0
        );
    const scalar liquidDensity =
        dropEquilibriumDict.lookupOrDefault<scalar>("liquidDensity", rhol.value());
    const scalar gravityMagnitude =
        dropEquilibriumDict.lookupOrDefault<scalar>("gravityMagnitude", 9.81);
    const scalar surfaceTension =
        dropEquilibriumDict.lookupOrDefault<scalar>("surfaceTension", 0.072);
    const scalar weight = liquidDensity*dropVolume*gravityMagnitude;
    const scalar forceTolerance =
        dropEquilibriumDict.lookupOrDefault<scalar>("forceTolerance", 1e-7);
    const scalar displacementRelaxation =
        dropEquilibriumDict.lookupOrDefault<scalar>
        (
            "displacementRelaxation",
            0.05
        );
    const scalar maxDeltaY =
        dropEquilibriumDict.lookupOrDefault<scalar>("maxDeltaY", 5e-5);
    const Switch enablePositionUpdate =
        dropEquilibriumDict.lookupOrDefault<Switch>("enablePositionUpdate", true);
    const scalar motionRadius =
        dropEquilibriumDict.lookupOrDefault<scalar>
        (
            "motionRadius",
            4.0*dropRadius
        );
    const label maxEquilibriumIter =
        dropEquilibriumDict.lookupOrDefault<label>("maxEquilibriumIter", 50);
    const label maxShapeIter =
        dropEquilibriumDict.lookupOrDefault<label>("maxShapeIter", 20);
    const label maxPositionIter =
        dropEquilibriumDict.lookupOrDefault<label>("maxPositionIter", 50);
    const Switch writeShapeProfiles =
        dropEquilibriumDict.lookupOrDefault<Switch>("writeShapeProfiles", true);
    const Switch enableShapeUpdate =
        dropEquilibriumDict.lookupOrDefault<Switch>("enableShapeUpdate", false);
    const scalar shapeRelaxation =
        dropEquilibriumDict.lookupOrDefault<scalar>("shapeRelaxation", 1e-9);
    const scalar maxShapeDelta =
        dropEquilibriumDict.lookupOrDefault<scalar>("maxShapeDelta", 2e-6);
    const scalar shapeTolerance =
        dropEquilibriumDict.lookupOrDefault<scalar>("shapeTolerance", 1e-2);
    const scalar maxAllowedShapeResidual =
        dropEquilibriumDict.lookupOrDefault<scalar>
        (
            "maxAllowedShapeResidual",
            500.0
        );
    const scalar maxAllowedForceMagnitude =
        dropEquilibriumDict.lookupOrDefault<scalar>
        (
            "maxAllowedForceMagnitude",
            1e-2
        );
    const scalar curvatureFallbackRadius =
        dropEquilibriumDict.lookupOrDefault<scalar>
        (
            "curvatureFallbackRadius",
            dropRadius
        );
    if (motionRadius <= dropRadius)
    {
        FatalErrorInFunction
            << "motionRadius must be larger than dropRadius." << exit(FatalError);
    }
    if (curvatureFallbackRadius <= VSMALL)
    {
        FatalErrorInFunction
            << "curvatureFallbackRadius must be positive." << exit(FatalError);
    }

    Info<< "\nStarting acoustic drop equilibrium loop\n"
        << "dropPatch=" << dropPatchName
        << " weight=" << weight
        << " surfaceTension=" << surfaceTension
        << " curvatureFallbackRadius=" << curvatureFallbackRadius
        << " forceTolerance=" << forceTolerance
        << " displacementRelaxation=" << displacementRelaxation
        << " maxDeltaY=" << maxDeltaY
        << " enablePositionUpdate=" << enablePositionUpdate
        << " motionRadius=" << motionRadius
        << " maxEquilibriumIter=" << maxEquilibriumIter
        << " maxShapeIter=" << maxShapeIter
        << " maxPositionIter=" << maxPositionIter
        << " enableShapeUpdate=" << enableShapeUpdate
        << " shapeRelaxation=" << shapeRelaxation
        << " maxShapeDelta=" << maxShapeDelta
        << " shapeTolerance=" << shapeTolerance
        << " maxAllowedShapeResidual=" << maxAllowedShapeResidual
        << " maxAllowedForceMagnitude=" << maxAllowedForceMagnitude
        << endl;

    rho = 1/(alpha1/rhol + (1 - alpha1)/rhog);
    invRhof = alphaf/rhol + (1 - alphaf)/rhog;

    volScalarField k2
    (
        IOobject
        (
            "k2",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        sqr(twoPi()*f)*(alpha1/sqr(cl) + (1 - alpha1)/sqr(cg))
    );

    label eqIter = 0;
    label shapeIter = 0;
    DropEquilibriumState state;

    bool positionConverged =
        runFixedShapePositionUpdate
        (
            simple,
            runTime,
            mesh,
            globalCells,
            N,
            M,
            x,
            b,
            ksp,
            dumpedMatrixStats,
            dumpedOpStats,
            rho,
            invRhof,
            k2,
            Pre,
            Pim,
            SC0,
            SC1,
            T0,
            T1,
            Ure,
            Uim,
            pa,
            pr,
            momFlux,
            alpha1,
            f,
            kl,
            kg,
            dropPatchID,
            axisymFactor,
            weight,
            liquidDensity,
            gravityMagnitude,
            surfaceTension,
            curvatureFallbackRadius,
            dropRadius,
            motionRadius,
            displacementRelaxation,
            maxDeltaY,
            enablePositionUpdate,
            forceTolerance,
            maxPositionIter,
            maxEquilibriumIter,
            shapeIter,
            eqIter,
            writeShapeProfiles,
            maxAllowedShapeResidual,
            maxAllowedForceMagnitude,
            state
        );

    while
    (
        positionConverged
     && enableShapeUpdate
     && state.shapeStats.rms > shapeTolerance
     && shapeIter < maxShapeIter
     && eqIter < maxEquilibriumIter
    )
    {
        ++shapeIter;

        applyShapeUpdate
        (
            runTime,
            mesh,
            pr,
            dropPatchID,
            axisymFactor,
            weight,
            liquidDensity,
            gravityMagnitude,
            surfaceTension,
            curvatureFallbackRadius,
            dropRadius,
            dropVolume,
            motionRadius,
            shapeRelaxation,
            maxShapeDelta,
            eqIter,
            shapeIter,
            state
        );

        // Andrade et al. Fig. 2: after each shape update, recompute the
        // levitation position for the updated fixed shape before testing the
        // next shape residual.
        positionConverged =
            runFixedShapePositionUpdate
            (
                simple,
                runTime,
                mesh,
                globalCells,
                N,
                M,
                x,
                b,
                ksp,
                dumpedMatrixStats,
                dumpedOpStats,
                rho,
                invRhof,
                k2,
                Pre,
                Pim,
                SC0,
                SC1,
                T0,
                T1,
                Ure,
                Uim,
                pa,
                pr,
                momFlux,
                alpha1,
                f,
                kl,
                kg,
                dropPatchID,
                axisymFactor,
                weight,
                liquidDensity,
                gravityMagnitude,
                surfaceTension,
                curvatureFallbackRadius,
                dropRadius,
                motionRadius,
                displacementRelaxation,
                maxDeltaY,
                enablePositionUpdate,
                forceTolerance,
                maxPositionIter,
                maxEquilibriumIter,
                shapeIter,
                eqIter,
                writeShapeProfiles,
                maxAllowedShapeResidual,
                maxAllowedForceMagnitude,
                state
            );
    }

    if
    (
        positionConverged
     && (!enableShapeUpdate || state.shapeStats.rms <= shapeTolerance)
    )
    {
        Info<< "Drop force/shape equilibrium reached." << nl;
    }

    KSPDestroy(&ksp);
    VecDestroy(&x);
    VecDestroy(&b);
    MatDestroy(&M);
    PetscFinalize();

    Info<< "End\n" << endl;
    return 0;
}

// ************************************************************************* //
