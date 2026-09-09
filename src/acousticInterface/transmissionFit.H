#ifndef acousticTransmissionFit_H
#define acousticTransmissionFit_H
#include "fvMesh.H"
#include "scalarMatrices.H"
#include "cutFaceAdvect.H"

namespace Foam { namespace acoustic
{
struct CutArea
{
    scalar fraction = 0;
    vector liquid = vector::zero;
    vector gas = vector::zero;
};
// Canonical triangle fan, identical after reversing/reordering face vertices.
CutArea cutArea
(
    const face&, const pointField&, const vector& Sf,
    const point& planeCentre, const vector& planeNormal, cutFaceAdvect&
);
struct FitResult
{
    bool valid = false;
    scalar condition = GREAT;
    scalarField weights;
};
// points[0:2] are the exact pressure anchors. activeBasis is orthonormal.
FitResult transmissionFit
(
    const pointField& points, const point& faceCentre,
    const point& planeCentre, const vector& planeNormal,
    const vector& liquidArea, const vector& gasArea,
    const vectorField& activeBasis, scalar rhoL, scalar rhoG,
    scalar svdTolerance, scalar maxCondition
);
}}
#endif
