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

Description
    The DAVOF dual advection with ONE curved clipping surface per cell
    (2026-10-23; dictionary key clippingSurface: none | exactSphere |
    paraboloid): the volume flux through a face is the liquid volume of the
    prism swept between the face's pre-image f - U dt and the face under the
    surface of the upwind cell, and the new face fraction is the liquid
    fraction of that pre-image under the same surface, so that alpha_k and
    alpha_f are the volume and the boundary moments of one surface at every
    step (the requirement of the sagitta studies, 2026-10-16). Replaces the
    plane sweep of timeIntegratedFlux() when active; serial runs only.

    exactSphere: the sphere of the exactNormal dictionary (origin + U0 t at
    the state's time) is every cell's surface: the machinery gate, exact to
    the quadrature tolerance under uniform translation by construction.

    paraboloid: the six-coefficient paraboloid of each interface cell in the
    frame of the normal the scheme used (offset, two slopes, the curvature
    tensor), all six fitted to the cell's own face fractions by the
    Gauss-Newton of leia's regeneration (2026-10-12) with the ring curvature
    tensor as the prior, from the one-scalar anchoring as the start; cells
    whose fit does not converge keep the best iterate (counted).

SourceFiles
    isoAdvectionCurved.C

\*---------------------------------------------------------------------------*/

#include "isoAdvection.H"
#include "davofCurvedGeometry.H"
#include "interpolationCellPoint.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::advection::isoAdvection::buildCellSurfaces(const scalar tState)
{
    addProfilingInFunction(geometricVoF);
    const label nCells = mesh_.nCells();
    if (cellSurfaces_.size() != nCells)
    {
        cellSurfaces_.setSize(nCells);
    }
    hasSurface_.reset();
    hasSurface_.resize(nCells);

    nFitCells_ = nFitConverged_ = nFitStalled_ = nFitSingular_ = nFitNoBracket_ = 0;
    nFitNoFrame_ = 0;
    fitResMax_ = 0;
    fitResSum_ = 0;

    const DynamicField<label>& interfaceLabels = surf_->interfaceLabels();

    if (clippingSurface_ == "exactSphere")
    {
        const dictionary& ed = modelDict().subDict("exactNormal");
        const point origin0 = ed.get<point>("origin");
        const vector U0 = ed.getOrDefault<vector>("U0", vector::zero);
        const scalar R = ed.get<scalar>("radius");
        davof::cellSurface S;
        S.sphere = true;
        S.o = origin0 + U0*tState;
        S.R = R;
        forAll(interfaceLabels, i)
        {
            const label celli = interfaceLabels[i];
            cellSurfaces_[celli] = S;
            hasSurface_.set(celli);
        }
        return;
    }

    // paraboloid: the ring curvature tensor of the start-of-step state
    davofPtr_->updateCurvatureTensor();
    const symmTensorField& B = davofPtr_->curvatureTensor().primitiveField();
    const vectorField& usedN = davofPtr_->usedNormals();
    const vectorField& normal = surf_->normal().primitiveField();
    const vectorField& centre = surf_->centre().primitiveField();
    const vectorField& m = davofPtr_->areaNormal().primitiveField();
    const vectorField& C = mesh_.cellCentres();
    const scalarField& V = mesh_.V();
    const vectorField& Sf = mesh_.faceAreas();
    const cellList& cells = mesh_.cells();

    if (writeCurvedFields_ && !fitResidualPtr_)
    {
        fitResidualPtr_.reset
        (
            new volScalarField
            (
                IOobject
                (
                    IOobject::groupName("fitResidualDavof", alpha1_.group()),
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh_,
                dimensionedScalar(dimless, Zero)
            )
        );
    }
    if (fitResidualPtr_)
    {
        fitResidualPtr_->primitiveFieldRef() = 0;
    }

    forAll(interfaceLabels, i)
    {
        const label celli = interfaceLabels[i];

        // The frame: e3 points out of the liquid (the used normal points in)
        vector nIn = usedN[i];
        if (!(mag(nIn) > 0))
        {
            const scalar mMag = mag(m[celli]);
            if (!(mMag > 0))
            {
                ++nFitNoFrame_;
                continue;
            }
            nIn = m[celli]/mMag;
        }
        davof::cellSurface S;
        S.e3 = -nIn;
        davof::frameOf(S.e3, S.e1, S.e2);
        S.o = (mag(normal[celli]) > 0) ? centre[celli] : C[celli];
        S.c = 0;
        const symmTensor& Bc = B[celli];
        S.b11 = -(S.e1 & (Bc & S.e1));
        S.b12 = -(S.e1 & (Bc & S.e2));
        S.b22 = -(S.e2 & (Bc & S.e2));
        S.paraboloid = true;
        const scalar h = Foam::cbrt(V[celli]);

        // The targets: the cell's own face fractions, weighted by the areas
        const labelList& cf = cells[celli];
        const label F = cf.size();
        scalarList tF(F), wF(F);
        scalar sumS = 0, target = 0;
        forAll(cf, k) sumS += mag(Sf[cf[k]]);
        forAll(cf, k)
        {
            const label f = cf[k];
            tF[k] = faceValue(*alphafPtr_, f);
            wF[k] = mag(Sf[f])/max(sumS, VSMALL);
            target += mag(Sf[f])*tF[k];
        }

        // The fractions of the cell's faces under a surface, in face order
        auto fractionsOf = [&](const davof::cellSurface& T, scalarList& aF)
        {
            aF.setSize(F);
            label nfb = 0;
            forAll(cf, k)
            {
                aF[k] = davof::faceCutDisplaced
                (
                    mesh_, cf[k], vector::zero,
                    [&](const point& x) { return T.value(x); }, nfb
                );
            }
        };

        // The start: the offset that reproduces the total wetted area with
        // the ring tensor and zero slopes (leia's paraboloidFaces anchoring)
        {
            auto g = [&](const scalar cc) -> scalar
            {
                davof::cellSurface T = S;
                T.c = cc;
                scalarList aF;
                fractionsOf(T, aF);
                scalar a = 0;
                forAll(cf, k) a += aF[k]*mag(Sf[cf[k]]);
                return a;
            };
            scalar cc = 0;
            bool bracketed = false;
            davof::solveOffset(g, target, -h, h, 1e-13*h, cc, bracketed);
            if (bracketed)
            {
                S.c = cc;
            }
            else
            {
                ++nFitNoBracket_;
            }
        }

        // All six coefficients
        ++nFitCells_;
        const davof::paraboloidFitResult res = davof::fitParaboloidFull
        (
            S, h, tF, wF, fractionsOf, fitPriorWeight_, fitTol_, fitMaxIter_
        );
        if (res.converged) ++nFitConverged_;
        if (res.stalled) ++nFitStalled_;
        if (res.singular) ++nFitSingular_;
        fitResMax_ = max(fitResMax_, res.rmax);
        fitResSum_ += res.rmax;
        if (fitResidualPtr_)
        {
            fitResidualPtr_->primitiveFieldRef()[celli] = res.rmax;
        }

        cellSurfaces_[celli] = S;
        hasSurface_.set(celli);
    }
}


void Foam::advection::isoAdvection::curvedFlux(const scalar dt)
{
    addProfilingInFunction(geometricVoF);

    if (Pstream::parRun())
    {
        FatalErrorInFunction
            << "clippingSurface " << clippingSurface_
            << " runs in serial only (no exchange of the surfaces across"
            << " processor boundaries)" << exit(FatalError);
    }

    // The records of the face fractions, as timeIntegratedFlux() keeps them
    clearIsoFaceData();
    alphafSetFace_.reset();
    alphafSetFace_.resize(mesh_.nFaces());
    alphafSetFaces_.clear();
    alphafFallbackFaces_.clear();
    if
    (
        mesh_.topoChanging()
     || alphaDefault_.size() != mesh_.nCells()
     || phiPos_.size() != mesh_.nInternalFaces()
     || wetFace_.size() != mesh_.nFaces()
    )
    {
        defaultsValid_ = false;
        wetFace_.resize(mesh_.nFaces());
    }
    nAlphafFromPlane_ = 0;
    nAlphafSecondRing_ = 0;
    nAlphafFallback_ = 0;
    nAlphafClipped_ = 0;
    nFluxSecondRing_ = 0;
    nContentBoundFaces_ = 0;
    nContentBoundCells_ = 0;
    nCurvedFaces_ = nCurvedSecondRing_ = nCurvedBoundary_ = nCurvedSkipped_ = 0;
    nCurvedFallback_ = 0;
    nCurvedLeaves_ = 0;
    defectMax_ = 0;
    defectSum_ = 0;
    nDefectCells_ = 0;

    // The surfaces of the state at t^n: the time has been advanced to
    // t^{n+1} when the advection runs (the scheme's stateTime() counts its
    // calls and cannot be asked again here)
    buildCellSurfaces(mesh_.time().value() - dt);

    interpolationCellPoint<vector> UInterp(U_);

    const labelList& own = mesh_.faceOwner();
    const labelList& nei = mesh_.faceNeighbour();
    const label nInt = mesh_.nInternalFaces();
    const labelListList& cellCells = mesh_.cellCells();
    const vectorField& Cf = mesh_.faceCentres();
    const vectorField& Sf = mesh_.faceAreas();
    const scalarField& phiIn = phi_.primitiveField();
    scalarField& dVfIn = dVf_.primitiveFieldRef();
    const DynamicField<label>& interfaceLabels = surf_->interfaceLabels();
    const polyBoundaryMesh& pbm = mesh_.boundaryMesh();

    // The candidate cells: the interface cells and their face neighbours
    // (the surface moves less than one cell per step); the interface cells
    // are the surface cells of the bounding (limitFluxes marks the cells
    // next to them from surfCells_)
    bitSet cand(mesh_.nCells());
    forAll(interfaceLabels, i)
    {
        const label celli = interfaceLabels[i];
        surfCells_.append(celli);
        cand.set(celli);
        for (const label nb : cellCells[celli])
        {
            cand.set(nb);
        }
    }

    // The surface of a face: the upwind cell's, else the other cell's, else
    // the closest surface cell among the face neighbours of both (second
    // ring); -1 when none
    auto surfaceCellOf = [&](const label up, const label other, const point& xf, bool& secondRing) -> label
    {
        secondRing = false;
        if (hasSurface_.test(up)) return up;
        if (other >= 0 && hasSurface_.test(other)) return other;
        label best = -1;
        scalar dMin = GREAT;
        auto scan = [&](const label c)
        {
            if (c < 0) return;
            for (const label nb : cellCells[c])
            {
                if (!hasSurface_.test(nb)) continue;
                const point& o = cellSurfaces_[nb].sphere ? mesh_.cellCentres()[nb] : cellSurfaces_[nb].o;
                const scalar d = magSqr(xf - o);
                if (d < dMin) { dMin = d; best = nb; }
            }
        };
        scan(up);
        scan(other);
        secondRing = (best >= 0);
        return best;
    };

    auto processFace = [&](const label facei, const label up, const label other, const scalar phiF, const bool boundary)
    {
        bool secondRing = false;
        const label sc = surfaceCellOf(up, other, Cf[facei], secondRing);
        if (sc < 0)
        {
            ++nCurvedSkipped_;
            return;
        }
        const davof::cellSurface& S = cellSurfaces_[sc];
        auto phi = [&](const point& x) { return S.value(x); };

        // The displacement of the step at the face
        const vector d = UInterp.interpolate(Cf[facei], up)*dt;

        // The volume flux: the liquid volume of the prism between f - d and f
        if (!(boundary && phiF < 0))
        {
            scalar vLiq = 0;
            if (mag(d & Sf[facei]) > VSMALL)
            {
                vLiq = davof::faceSweptLiquidVolume
                (
                    mesh_, facei, d, phi, curvedTolRel_, curvedMaxDepth_, nCurvedLeaves_
                );
            }
            const scalar dVfNew = (phiF >= 0 ? vLiq : -vLiq);
            if (boundary)
            {
                setFaceValue(dVf_, facei, dVfNew);
            }
            else
            {
                dVfIn[facei] = dVfNew;
            }
        }

        // The face fraction: the pre-image's liquid fraction under the surface
        label nfb = 0;
        scalar af = davof::faceCutDisplaced(mesh_, facei, -d, phi, nfb);
        nCurvedFallback_ += nfb;
        if (af < -1e-12 || af > 1 + 1e-12)
        {
            ++nAlphafClipped_;
        }
        af = min(max(af, scalar(0)), scalar(1));
        setFaceValue(*alphafPtr_, facei, af);
        wetFace_.set(facei, (af != 0) ? 1u : 0u);
        if (!alphafSetFace_.test(facei))
        {
            alphafSetFaces_.append(facei);
        }
        alphafSetFace_.set(facei);
        ++nCurvedFaces_;
        if (secondRing) ++nCurvedSecondRing_;
        if (boundary) ++nCurvedBoundary_;
    };

    for (label facei = 0; facei < nInt; ++facei)
    {
        const label o = own[facei];
        const label n = nei[facei];
        if (!cand.test(o) && !cand.test(n)) continue;
        const scalar phiF = phiIn[facei];
        const label up = (phiF >= 0) ? o : n;
        const label other = (phiF >= 0) ? n : o;
        processFace(facei, up, other, phiF, false);
    }
    forAll(pbm, patchi)
    {
        const polyPatch& pp = pbm[patchi];
        if (alphafPtr_->boundaryField()[patchi].size() != pp.size()) continue;
        if (pp.coupled()) continue;
        forAll(pp, i)
        {
            const label facei = pp.start() + i;
            const label o = own[facei];
            if (!cand.test(o)) continue;
            processFace(facei, o, -1, faceValue(phi_, facei), true);
        }
    }

    // The paraboloid defect: alpha_k against the volume under the cell's
    // surface, at the start of the step (diagnostic)
    if (defectLevels_ >= 0)
    {
        if (writeCurvedFields_ && !paraboloidDefectPtr_)
        {
            paraboloidDefectPtr_.reset
            (
                new volScalarField
                (
                    IOobject
                    (
                        IOobject::groupName("paraboloidDefectDavof", alpha1_.group()),
                        mesh_.time().timeName(),
                        mesh_,
                        IOobject::NO_READ,
                        IOobject::AUTO_WRITE
                    ),
                    mesh_,
                    dimensionedScalar(dimless, Zero)
                )
            );
        }
        if (paraboloidDefectPtr_)
        {
            paraboloidDefectPtr_->primitiveFieldRef() = 0;
        }
        const scalarField& V = mesh_.V();
        forAll(interfaceLabels, i)
        {
            const label celli = interfaceLabels[i];
            if (!hasSurface_.test(celli)) continue;
            const davof::cellSurface& S = cellSurfaces_[celli];
            const scalar vLiq = davof::cellLiquidVolumeRefined
            (
                mesh_, celli, [&](const point& x) { return S.value(x); }, defectLevels_
            );
            const scalar defect = alpha1In_[celli] - vLiq/V[celli];
            defectMax_ = max(defectMax_, mag(defect));
            defectSum_ += mag(defect);
            ++nDefectCells_;
            if (paraboloidDefectPtr_)
            {
                paraboloidDefectPtr_->primitiveFieldRef()[celli] = defect;
            }
        }
    }
}


void Foam::advection::isoAdvection::reportCurved() const
{
    Info<< "isoAdvection: curved clipping (" << clippingSurface_
        << "): alphaf and dVf set on " << nCurvedFaces_ << " faces ("
        << nCurvedSecondRing_ << " second ring, " << nCurvedBoundary_
        << " boundary), " << nCurvedSkipped_ << " skipped, "
        << nAlphafClipped_ << " clips, " << nCurvedFallback_
        << " linear cuts, " << nCurvedLeaves_ << " quadrature leaves";
    if (clippingSurface_ == "paraboloid")
    {
        Info<< "; fit: cells " << nFitCells_ << ", converged " << nFitConverged_
            << ", stalled " << nFitStalled_ << ", singular " << nFitSingular_
            << ", no bracket " << nFitNoBracket_ << ", no frame " << nFitNoFrame_
            << ", residual max " << fitResMax_
            << " mean " << fitResSum_/max(nFitCells_, label(1));
    }
    if (nDefectCells_ > 0)
    {
        Info<< "; surface defect max " << defectMax_
            << " L1 " << defectSum_/nDefectCells_
            << " over " << nDefectCells_ << " cells";
    }
    Info<< endl;
}


// ************************************************************************* //
