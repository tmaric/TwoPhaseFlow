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
    The tilt diagnostic of the curved dual advection on the translated
    plane (key tiltDiagnostics, default false; serial; 2026-10-26). The
    question: the own normal of the loop fails the plane gate by a slow
    geometric growth of the face-fraction error although every operation
    is exact for a plane. The claimed mechanism: a_k = -sum_f alpha_f S_f
    is exact for the fractions of ONE surface, but the fractions of cell
    k's faces are set by DIFFERENT surfaces (the surface of the cell the
    code chose for the face: the upwind cell, the lexicographically smaller
    centre on a zero-flux face, a second-ring surface). The linearised
    prediction of the turn of a_k^{n+1} away from the exact normal n (out
    of the liquid):

        tau_pred = (1/A_k) sum_f ell_f d_{u(f)}(m_f) t_f

    over the faces f of k whose pre-image f - U dt carries a trace of the
    exact plane Pi^n (length ell_f, midpoint m_f), u(f) the cell whose
    surface set alpha_f^{n+1} (recorded where the choice is made, in
    curvedFlux), d_u(x) = (x - x_u).n_u the signed distance from the
    tangent plane of that surface at its origin (n_u = e3 out of the
    liquid, x_u = o + c e3), t_f the unit transverse part of the outward
    face normal, A_k the wetted patch area of the exact plane Pi^{n+1} in
    cell k (cutCellPLIC). Measured: tau_meas = the transverse part of a_k
    over |a_k| from the stored fractions at the end of the step. Weights
    w_f = ell_f |m_f - x*_u|/A_k, x*_u the centroid of the wetted patch of
    the surface's tangent plane in cell u (cutCellPLIC), the neighbour's
    tilt e_u = angle(n_u, n) and offset p_u = (x*_u - x_Pi).n.

    The sign of the exact normal: the exactNormal dictionary's plane
    (alphaFieldFunctions: value(p) = (p - origin).normal, the gradient
    points INTO the liquid as the scheme's normals do) is checked against
    alpha at the first step and the normal used here points OUT of the
    liquid; the outcome is printed.

    Output (postProcessing/tiltDiagnostics/): cells_<step>.csv, one line
    per cell with 0 < alpha < 1 at the end of the step (format in
    etc/README-davof.md), and ownSetsTrace.csv (step, celli, flag: 1 when a
    face set by the cell's surface in that step carries a trace), appended.

SourceFiles
    isoAdvectionTiltDiagnostics.C

\*---------------------------------------------------------------------------*/

#include "isoAdvection.H"
#include "cutCellPLIC.H"
#include "OFstream.H"
#include "OSspecific.H"
#include <cmath>

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace
{

// The trace of the plane (xP, nP) on the polygon f shifted by `shift`:
// the two crossings of the plane with the polygon's edges (a vertex on
// the plane counts as non-negative); false when the plane misses the
// polygon or grazes a vertex only.
bool planeTraceOnPolygon
(
    const Foam::face& f,
    const Foam::pointField& P,
    const Foam::vector& shift,
    const Foam::point& xP,
    const Foam::vector& nP,
    Foam::point& m,
    Foam::scalar& ell
)
{
    using namespace Foam;
    point q[2];
    label nq = 0;
    forAll(f, i)
    {
        const point pi = P[f[i]] + shift;
        const point pj = P[f.nextLabel(i)] + shift;
        const scalar vi = (pi - xP) & nP;
        const scalar vj = (pj - xP) & nP;
        if ((vi < 0) != (vj < 0))
        {
            const point x = pi + (vi/(vi - vj))*(pj - pi);
            if (nq < 2)
            {
                q[nq] = x;
            }
            ++nq;
        }
    }
    if (nq != 2) return false;
    ell = mag(q[1] - q[0]);
    m = 0.5*(q[0] + q[1]);
    return ell > 0;
}

} // End anonymous namespace


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::advection::isoAdvection::tiltDiagSetup()
{
    if (!modelDict().found("exactNormal"))
    {
        FatalIOErrorInFunction(modelDict())
            << "tiltDiagnostics needs the exactNormal sub-dictionary"
            << " (functionType plane, origin, normal, U0)"
            << exit(FatalIOError);
    }
    const dictionary& ed = modelDict().subDict("exactNormal");
    const word ft = ed.getOrDefault<word>("functionType", "plane");
    if (ft != "plane")
    {
        FatalIOErrorInFunction(modelDict())
            << "tiltDiagnostics is the plane test's diagnostic; functionType "
            << ft << " is not a plane" << exit(FatalIOError);
    }
    tiltPlaneP0_ = ed.get<point>("origin");
    const vector n = ed.get<vector>("normal");
    tiltPlaneN_ = n/mag(n);
    tiltPlaneU_ = ed.getOrDefault<vector>("U0", vector::zero);
    tiltSignChecked_ = false;
    tiltDir_ = mesh_.time().path()/"postProcessing"/"tiltDiagnostics";
    mkDir(tiltDir_);
    // the per-step join list, started afresh
    {
        OFstream os(tiltDir_/"ownSetsTrace.csv");
        os  << "step,celli,ownSetsTrace" << nl;
    }
    Info<< "isoAdvection: tiltDiagnostics on (the per-cell tilt prediction"
        << " of the translated plane, postProcessing/tiltDiagnostics/)"
        << endl;
}


void Foam::advection::isoAdvection::tiltDiagBegin()
{
    const label nFaces = mesh_.nFaces();
    if (tiltFaceRec_.size() != nFaces)
    {
        tiltFaceRec_.setSize(nFaces);
        tiltFaceRec_ = -1;
    }
    else
    {
        forAll(tiltRecFace_, r)
        {
            tiltFaceRec_[tiltRecFace_[r]] = -1;
        }
    }
    tiltRecFace_.clear();
    tiltRecCell_.clear();
    tiltRecD_.clear();
    tiltRecEll_.clear();
    tiltRecM_.clear();
    tiltRecTrace_.clear();
}


void Foam::advection::isoAdvection::tiltDiagRecordFace
(
    const label facei,
    const label surfCell,
    const vector& d
)
{
    tiltFaceRec_[facei] = tiltRecFace_.size();
    tiltRecFace_.append(facei);
    tiltRecCell_.append(surfCell);
    tiltRecD_.append(d);
}


void Foam::advection::isoAdvection::tiltDiagFaces(const scalar dt)
{
    if (Pstream::parRun())
    {
        WarningInFunction
            << "tiltDiagnostics is serial only; nothing written" << endl;
        return;
    }
    const scalar tState = mesh_.time().value() - dt;
    const point xP = tiltPlaneP0_ + tiltPlaneU_*tState;
    const vectorField& C = mesh_.cellCentres();

    // The sign of the exact normal against alpha at the start of the step
    // (alpha1In_ is still alpha^n here): the normal used points OUT of the
    // liquid
    if (!tiltSignChecked_)
    {
        const DynamicField<label>& labels = surf_->interfaceLabels();
        scalar s = 0;
        forAll(labels, i)
        {
            const label celli = labels[i];
            s += ((C[celli] - xP) & tiltPlaneN_)*(alpha1In_[celli] - 0.5);
        }
        if (s > 0)
        {
            tiltPlaneN_ = -tiltPlaneN_;
            Info<< "isoAdvection: tiltDiagnostics: the liquid lies on the +normal"
                << " side of the exactNormal plane; the outward normal is "
                << tiltPlaneN_ << endl;
        }
        else
        {
            Info<< "isoAdvection: tiltDiagnostics: the liquid lies on the -normal"
                << " side of the exactNormal plane; the outward normal is "
                << tiltPlaneN_ << endl;
        }
        tiltSignChecked_ = true;
    }

    // The traces of Pi^n on the pre-images of the processed faces
    const pointField& P = mesh_.points();
    const faceList& faces = mesh_.faces();
    const label nRec = tiltRecFace_.size();
    tiltRecEll_.setSize(nRec);
    tiltRecM_.setSize(nRec);
    tiltRecTrace_.setSize(nRec);
    bitSet ownSets(mesh_.nCells());
    for (label r = 0; r < nRec; ++r)
    {
        point m(Zero);
        scalar ell = 0;
        const bool has = planeTraceOnPolygon
        (
            faces[tiltRecFace_[r]], P, -tiltRecD_[r], xP, tiltPlaneN_, m, ell
        );
        tiltRecTrace_[r] = has;
        tiltRecEll_[r] = has ? ell : 0;
        tiltRecM_[r] = m;
        if (has && tiltRecCell_[r] >= 0)
        {
            ownSets.set(tiltRecCell_[r]);
        }
    }

    // The join list: every cell with a surface in this step
    {
        OFstream os(tiltDir_/"ownSetsTrace.csv", IOstreamOption(), IOstreamOption::APPEND_ATE);
        const label step = mesh_.time().timeIndex();
        const DynamicField<label>& labels = surf_->interfaceLabels();
        forAll(labels, i)
        {
            const label celli = labels[i];
            if (!hasSurface_.test(celli)) continue;
            os  << step << ',' << celli << ',' << (ownSets.test(celli) ? 1 : 0) << nl;
        }
    }
}


void Foam::advection::isoAdvection::tiltDiagCells()
{
    if (Pstream::parRun()) return;
    if (tiltRecFace_.empty()) return;

    const label step = mesh_.time().timeIndex();
    const scalar dt = mesh_.time().deltaTValue();
    const scalar t1 = mesh_.time().value();
    const point xP1 = tiltPlaneP0_ + tiltPlaneU_*t1;          // Pi^{n+1}
    const point xP0 = tiltPlaneP0_ + tiltPlaneU_*(t1 - dt);   // Pi^n
    const vector& nP = tiltPlaneN_;

    const cellList& cells = mesh_.cells();
    const labelList& own = mesh_.faceOwner();
    const vectorField& Sf = mesh_.faceAreas();
    const vectorField& C = mesh_.cellCentres();
    const scalarField& V = mesh_.V();
    const label nCells = mesh_.nCells();
    const label nInt = mesh_.nInternalFaces();
    const polyBoundaryMesh& pbm = mesh_.boundaryMesh();

    cutCellPLIC cutter(mesh_);

    // The per-neighbour cache of the step: the centroid x*_u of the wetted
    // patch of the tangent plane of u's surface in cell u, its tilt e_u and
    // offset p_u from Pi^n
    if (tiltStamp_.size() != nCells)
    {
        tiltStamp_.setSize(nCells);
        tiltStamp_ = -1;
        tiltXStar_.setSize(nCells);
        tiltTilt_.setSize(nCells);
        tiltOffset_.setSize(nCells);
        tiltCutOK_.setSize(nCells);
    }
    auto neighbourData = [&](const label u)
    {
        if (tiltStamp_[u] == step) return;
        tiltStamp_[u] = step;
        const davof::cellSurface& S = cellSurfaces_[u];
        const point xu = S.o + S.c*S.e3;
        const label status = cutter.calcSubCell(u, (xu - C[u]) & S.e3, S.e3);
        if (status == 0 && mag(cutter.faceArea()) > 0)
        {
            tiltXStar_[u] = cutter.faceCentre();
            tiltCutOK_[u] = true;
        }
        else
        {
            tiltXStar_[u] = C[u] + ((xu - C[u]) & S.e3)*S.e3;
            tiltCutOK_[u] = false;
        }
        // the angle through atan2 of the sine and the cosine: acos of a
        // dot product rounded to 1 - 1e-16 floors the angle at 1.5e-8,
        // above every tilt of the first steps of the plane test
        tiltTilt_[u] = std::atan2(mag(S.e3 ^ nP), S.e3 & nP);
        tiltOffset_[u] = (tiltXStar_[u] - xP0) & nP;
    };

    // The face fraction of any face (faces of patches without values: 0)
    auto alphaOf = [&](const label facei) -> scalar
    {
        if (facei < nInt) return alphafPtr_->primitiveField()[facei];
        const label patchi = pbm.whichPatch(facei);
        const fvsPatchScalarField& pf = alphafPtr_->boundaryField()[patchi];
        if (pf.size() != pbm[patchi].size()) return 0;
        return pf[facei - pbm[patchi].start()];
    };

    OFstream os(tiltDir_/("cells_" + Foam::name(step) + ".csv"));
    os.precision(12);
    os  << "step,celli,alpha,A_over_h2,a_over_h2,nFaces,nTraced,nForeign,nNoSurface,"
        << "sumAbsW,maxW,tauMeas,tauPred,tauDiff,cosine,maxNeighbourTilt,"
        << "maxNeighbourOffset,maxTraceDist,selfTilt,selfOffset,exactCut" << nl;

    label nWritten = 0;
    for (label celli = 0; celli < nCells; ++celli)
    {
        const scalar a = alpha1In_[celli];
        if (!(a > 1e-12 && a < 1 - 1e-12)) continue;

        const labelList& cf = cells[celli];
        const scalar h = Foam::cbrt(V[celli]);
        const scalar h2 = h*h;

        // a_k = -sum_f alpha_f S_f^out, the measured turn
        vector ak(Zero);
        forAll(cf, k)
        {
            const label f = cf[k];
            const vector SfOut = (own[f] == celli) ? Sf[f] : -Sf[f];
            ak -= alphaOf(f)*SfOut;
        }
        const scalar magA = mag(ak);
        if (!(magA > 0)) continue;
        const vector tauMeas = (ak - (ak & nP)*nP)/magA;

        // A_k: the wetted patch of Pi^{n+1} in the cell
        scalar A = 0;
        label exactCut = 0;
        {
            const label status = cutter.calcSubCell(celli, (xP1 - C[celli]) & nP, nP);
            if (status == 0)
            {
                A = mag(cutter.faceArea());
                exactCut = 1;
            }
        }

        // The prediction over the traced faces
        vector tauPred(Zero);
        label nTraced = 0, nForeign = 0, nNoSurface = 0;
        scalar sumAbsW = 0, maxW = 0, maxTilt = 0, maxOffset = 0, maxTraceDist = 0;
        forAll(cf, k)
        {
            const label f = cf[k];
            const label r = tiltFaceRec_[f];
            if (r < 0) continue;                   // not processed: no surface reached it
            if (!tiltRecTrace_[r]) continue;
            ++nTraced;
            const label u = tiltRecCell_[r];
            if (u < 0)
            {
                ++nNoSurface;
                continue;
            }
            if (u != celli) ++nForeign;
            neighbourData(u);
            const davof::cellSurface& S = cellSurfaces_[u];
            const point xu = S.o + S.c*S.e3;
            const point& m = tiltRecM_[r];
            const scalar ell = tiltRecEll_[r];
            const scalar du = (m - xu) & S.e3;
            const vector SfOut = (own[f] == celli) ? Sf[f] : -Sf[f];
            const vector s = SfOut/mag(SfOut);
            vector sPerp = s - (s & nP)*nP;
            const scalar magPerp = mag(sPerp);
            if (magPerp > 1e-8)
            {
                const vector tf = sPerp/magPerp;
                if (A > 0) tauPred += (ell*du/A)*tf;
            }
            const scalar w = (A > 0) ? ell*mag(m - tiltXStar_[u])/A : 0;
            sumAbsW += mag(w);
            maxW = max(maxW, mag(w));
            maxTilt = max(maxTilt, tiltTilt_[u]);
            maxOffset = max(maxOffset, mag(tiltOffset_[u])/h);
            maxTraceDist = max(maxTraceDist, mag(du)/h);
        }

        scalar selfTilt = -1, selfOffset = 0;
        if (hasSurface_.test(celli))
        {
            neighbourData(celli);
            selfTilt = tiltTilt_[celli];
            selfOffset = tiltOffset_[celli]/h;
        }

        const scalar magM = mag(tauMeas), magP = mag(tauPred);
        const scalar cosine =
            (magM > 0 && magP > 0) ? (tauMeas & tauPred)/(magM*magP) : 0;

        os  << step << ',' << celli << ',' << a << ',' << A/h2 << ','
            << magA/h2 << ',' << cf.size() << ',' << nTraced << ','
            << nForeign << ',' << nNoSurface << ',' << sumAbsW << ','
            << maxW << ',' << magM << ',' << magP << ','
            << mag(tauMeas - tauPred) << ',' << cosine << ','
            << maxTilt << ',' << maxOffset << ',' << maxTraceDist << ','
            << selfTilt << ',' << selfOffset << ',' << exactCut << nl;
        ++nWritten;
    }
    Info<< "isoAdvection: tiltDiagnostics: " << nWritten
        << " cells written to " << tiltDir_.name()/("cells_" + Foam::name(step) + ".csv")
        << endl;
}


// ************************************************************************* //
