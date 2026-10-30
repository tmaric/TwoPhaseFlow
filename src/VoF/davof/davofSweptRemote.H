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
    The cells of another rank in the vertex ring of a donor cell of the
    swept-face fluxes (fluxIntegration sweptFace and sweptFaceCells,
    isoAdvectionSwept.C), 2026-10-30.

    The engineering problem: the swept volume of a face is cut into pieces in
    the donor and in every cell sharing a vertex with it (one layer), and the
    backward trace of a vertex interpolates the velocity in the cell that
    holds the trajectory's midpoint. Both read the neighbour cells by label
    from the local mesh, so under domain decomposition a neighbour across a
    processor boundary was not seen (its part of the swept volume counted as
    beyond the layer, the midpoint fell back to the Euler step). The
    neighbour set is the cell-point-cell stencil of the donor, which
    OpenFOAM's zoneDistribute / zoneCPCStencil provides across the processor
    boundaries with global cell indices (introduced for the plicRDF
    reconstruction by H. Scheufler, J. Roenby, Accurate and efficient surface
    reconstruction from volume fraction data on general meshes, J. Comput.
    Phys. 383 (2019) 1-23, doi:10.1016/j.jcp.2019.01.009). The stencil's
    remote cells are fetched once per step through a mapDistributeBase built
    from their global indices; each remote cell travels as one packed
    scalarList (below), because the face planes and the tetrahedra of a cell
    are lists that a volume field cannot carry.

    A remote cell carries exactly the data the local path reads, bit for bit
    as the owning rank holds it:
      - the state of t^n: alpha^n, the surface flag and the surface (the 21
        scalars of the processor-patch exchange of isoAdvectionCurved.C);
      - the bounding box of its points (the ring prefilter);
      - per face the outward area vector and the face centre (the
        half-spaces of the pieces and the point-in-cell test, which
        normalise them on arrival with the same arithmetic as the local
        path);
      - the tetrahedra of its cell-point decomposition (the cell centre and
        the three points of each face triangle, OpenFOAM's tetIndices with
        the face's tet base point) with the volPointInterpolation velocity
        at the three points, the cell velocity and the cell volume, so that
        the velocity at a point in the cell is OpenFOAM's
        interpolationCellPoint value (the tetrahedron search of
        cellPointWeight::findTetrahedron, replicated below with the same
        tolerance and the same order of the sums).

SourceFiles
    davofSweptRemote.H

\*---------------------------------------------------------------------------*/

#ifndef davofSweptRemote_H
#define davofSweptRemote_H

#include "davofCurvedGeometry.H"
#include "davofSweptPieces.H"
#include "tetrahedron.H"
#include "barycentric.H"
#include "cellPointWeight.H"
#include "boundBox.H"
#include "FixedList.H"
#include "List.H"

namespace Foam
{
namespace davof
{

//- A cell of another rank in the vertex ring of a donor
struct remoteCell
{
    //- The global index (zoneDistribute's numbering)
    label gid = -1;

    point C = Zero;
    scalar V = 0;
    scalar alpha0 = 0;
    bool hasSurface = false;
    cellSurface S;
    point bbMin = Zero, bbMax = Zero;
    vector Ucell = Zero;

    //- Per face: the outward area vector and the face centre
    List<vector> Sout;
    List<point> Cf;

    //- Per tetrahedron (cell centre, a, b, c): the points a, b, c and the
    //  volPointInterpolation velocity at them
    List<FixedList<point, 3>> tetP;
    List<FixedList<vector, 3>> tetU;
};


//- The 21 scalars of a surface (the packing of the processor-patch exchange)
inline void packCellSurface(const cellSurface& S, DynamicList<scalar>& v)
{
    v.append(S.paraboloid ? 1 : 0);
    v.append(S.sphere ? 1 : 0);
    for (direction k = 0; k < 3; ++k) v.append(S.o[k]);
    for (direction k = 0; k < 3; ++k) v.append(S.e1[k]);
    for (direction k = 0; k < 3; ++k) v.append(S.e2[k]);
    for (direction k = 0; k < 3; ++k) v.append(S.e3[k]);
    v.append(S.b11); v.append(S.b12); v.append(S.b22);
    v.append(S.c); v.append(S.a1); v.append(S.a2); v.append(S.R);
}


//- The packed cell: [1, C, V, alpha0, has, S(21), bbMin, bbMax, Ucell,
//  nFaces, (Sout, Cf) per face, nTets, (a, b, c, Ua, Ub, Uc) per tet]; a
//  stencil entry that is a boundary face packs [0]
inline void appendVector(const vector& x, DynamicList<scalar>& v)
{
    v.append(x.x()); v.append(x.y()); v.append(x.z());
}

inline vector readVector(const UList<scalar>& v, label& i)
{
    const vector x(v[i], v[i + 1], v[i + 2]);
    i += 3;
    return x;
}


//- Unpack a cell; false for a boundary-face entry of the stencil
inline bool unpackRemoteCell(const UList<scalar>& v, remoteCell& r)
{
    if (v.empty() || v[0] < 0.5) return false;
    label i = 1;
    r.C = readVector(v, i);
    r.V = v[i++];
    r.alpha0 = v[i++];
    r.hasSurface = (v[i++] > 0.5);
    cellSurface& S = r.S;
    S.paraboloid = (v[i++] > 0.5);
    S.sphere = (v[i++] > 0.5);
    S.o = readVector(v, i);
    S.e1 = readVector(v, i);
    S.e2 = readVector(v, i);
    S.e3 = readVector(v, i);
    S.b11 = v[i++]; S.b12 = v[i++]; S.b22 = v[i++];
    S.c = v[i++]; S.a1 = v[i++]; S.a2 = v[i++]; S.R = v[i++];
    r.bbMin = readVector(v, i);
    r.bbMax = readVector(v, i);
    r.Ucell = readVector(v, i);
    const label nF = label(v[i++] + 0.5);
    r.Sout.setSize(nF);
    r.Cf.setSize(nF);
    for (label j = 0; j < nF; ++j)
    {
        r.Sout[j] = readVector(v, i);
        r.Cf[j] = readVector(v, i);
    }
    const label nT = label(v[i++] + 0.5);
    r.tetP.setSize(nT);
    r.tetU.setSize(nT);
    for (label j = 0; j < nT; ++j)
    {
        for (direction k = 0; k < 3; ++k) r.tetP[j][k] = readVector(v, i);
        for (direction k = 0; k < 3; ++k) r.tetU[j][k] = readVector(v, i);
    }
    return true;
}


//- The point-in-cell test by the face planes with a tolerance in units of
//  h (the local inCellTol of isoAdvectionSwept.C on the packed planes)
inline bool inRemoteCellTol
(
    const remoteCell& r, const point& x, const scalar tolRel, const scalar h
)
{
    forAll(r.Sout, j)
    {
        const vector& nOut = r.Sout[j];
        if (((x - r.Cf[j]) & nOut) > tolRel*h*mag(nOut))
        {
            return false;
        }
    }
    return true;
}


//- The half-spaces of the cell (the local path's normalisation)
inline void remoteCellPlanes(const remoteCell& r, List<halfSpace>& planes)
{
    planes.setSize(r.Sout.size());
    forAll(r.Sout, j)
    {
        const vector& nOut = r.Sout[j];
        planes[j].n = nOut/max(mag(nOut), VSMALL);
        planes[j].x0 = r.Cf[j];
    }
}


//- The cell-point interpolate of the velocity at x in the remote cell:
//  cellPointWeight::findTetrahedron (the first tetrahedron whose
//  barycentric coordinates are inside within cellPointWeight::tol, else the
//  nearest) and interpolationCellPoint::interpolate (cell value times the
//  first weight plus the point values times the others, in that order)
inline vector interpolateRemote(const remoteCell& r, const point& x)
{
    const scalar tol = cellPointWeight::tol;
    barycentric w;
    label found = -1;
    forAll(r.tetP, k)
    {
        const tetPointRef t(r.C, r.tetP[k][0], r.tetP[k][1], r.tetP[k][2]);
        const scalar det = t.pointToBarycentric(x, w);
        if (mag(det/r.V) > tol)
        {
            const scalar& u = w[0];
            const scalar& v = w[1];
            const scalar& ww = w[2];
            if ((u + tol > 0) && (v + tol > 0) && (ww + tol > 0) && (u + v + ww < 1 + tol))
            {
                found = k;
                break;
            }
        }
    }
    if (found < 0)
    {
        scalar minNearDist = VGREAT;
        forAll(r.tetP, k)
        {
            const tetPointRef t(r.C, r.tetP[k][0], r.tetP[k][1], r.tetP[k][2]);
            const scalar d = t.nearestPoint(x).distance();
            if (d < minNearDist)
            {
                minNearDist = d;
                found = k;
            }
        }
        if (found < 0) return r.Ucell;
        const tetPointRef t(r.C, r.tetP[found][0], r.tetP[found][1], r.tetP[found][2]);
        w = t.pointToBarycentric(x);
    }
    vector U = r.Ucell*w[0];
    U += r.tetU[found][0]*w[1];
    U += r.tetU[found][1]*w[2];
    U += r.tetU[found][2]*w[3];
    return U;
}

} // End namespace davof
} // End namespace Foam

#endif

// ************************************************************************* //
