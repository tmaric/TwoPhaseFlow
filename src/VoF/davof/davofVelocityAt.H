/*---------------------------------------------------------------------------*\
    The velocity at a point of a cell for the swept-face fluxes (fork
    davof-swept, 2026-11-05, round 4): interpolationCellPoint's cell-point
    interpolate (icp), or the same arithmetic on the cell's tetrahedra cached
    per cell (key sweptTetCache): the tetrahedra of
    polyMeshTetDecomposition::cellTetIndices as the point triples of
    tetIndices::faceTriIs, built once per cell on first use and kept while
    the mesh stands, the search and the weights of
    cellPointWeight::findTetrahedron and interpolationCellPoint::interpolate
    repeated operation by operation (the first tetrahedron whose barycentric
    coordinates are inside within cellPointWeight::tol, else the nearest;
    the cell value times the first weight plus the point values times the
    others, in that order), so that the result is bit for bit the library's
    without the per-lookup tetrahedron list. The point velocities Up are
    asked of the lazy evaluator first where one is given (ensure; key
    sweptPointVelocityLazy).
\*---------------------------------------------------------------------------*/

#ifndef davofVelocityAt_H
#define davofVelocityAt_H

#include "polyMesh.H"
#include "vectorField.H"
#include "interpolationCellPoint.H"
#include "polyMeshTetDecomposition.H"
#include "tetIndices.H"
#include "cellPointWeight.H"
#include "triFace.H"
#include <functional>

namespace Foam
{
namespace davof
{

class velocityAt
{
    const polyMesh& mesh_;
    const vectorField& U_;
    const vectorField& Up_;
    const interpolationCellPoint<vector>* icp_;
    List<List<triFace>>* cellTets_;
    std::function<void(label)> ensure_;

public:

    velocityAt
    (
        const polyMesh& mesh,
        const vectorField& U,
        const vectorField& Up,
        const interpolationCellPoint<vector>* icp,
        List<List<triFace>>* cellTets,
        std::function<void(label)> ensure
    )
    :
        mesh_(mesh), U_(U), Up_(Up), icp_(icp), cellTets_(cellTets), ensure_(ensure)
    {}

    //- The cell's tetrahedra as point triples, cached
    const List<triFace>& tets(const label celli) const
    {
        List<triFace>& t = (*cellTets_)[celli];
        if (t.empty())
        {
            const List<tetIndices> ti(polyMeshTetDecomposition::cellTetIndices(mesh_, celli));
            t.setSize(ti.size());
            forAll(ti, k) t[k] = ti[k].faceTriIs(mesh_);
        }
        return t;
    }

    //- The velocity at x in cell celli
    vector interpolate(const point& x, const label celli) const
    {
        if (icp_) return icp_->interpolate(x, celli);
        const List<triFace>& t = tets(celli);
        const pointField& P = mesh_.points();
        const point& C = mesh_.cellCentres()[celli];
        const scalar cellVolume = mesh_.cellVolumes()[celli];
        const scalar tol = cellPointWeight::tol;
        barycentric w;
        label found = -1;
        forAll(t, k)
        {
            const tetPointRef tet(C, P[t[k][0]], P[t[k][1]], P[t[k][2]]);
            const scalar det = tet.pointToBarycentric(x, w);
            if (mag(det/cellVolume) > tol)
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
            forAll(t, k)
            {
                const tetPointRef tet(C, P[t[k][0]], P[t[k][1]], P[t[k][2]]);
                const scalar d = tet.nearestPoint(x).distance();
                if (d < minNearDist)
                {
                    minNearDist = d;
                    found = k;
                }
            }
            const tetPointRef tet(C, P[t[found][0]], P[t[found][1]], P[t[found][2]]);
            w = tet.pointToBarycentric(x);
        }
        const triFace& tri = t[found];
        if (ensure_)
        {
            ensure_(tri[0]); ensure_(tri[1]); ensure_(tri[2]);
        }
        vector r = U_[celli]*w[0];
        r += Up_[tri[0]]*w[1];
        r += Up_[tri[1]]*w[2];
        r += Up_[tri[2]]*w[3];
        return r;
    }
};

} // End namespace davof
} // End namespace Foam

#endif
