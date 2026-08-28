/*---------------------------------------------------------------------------*\
Class
    Foam::movingFrameFlow

Description
    Prescribed, exactly periodic, non-retracing advection test field built by
    superposing a T-periodic Euclidean frame onto a forward-backward base
    vortex, following

        Bothe & Maric, "A Recipe for Periodic, Non-Retracing Deformation Test
        Cases" (2026), Sections 3 and 8.

    With frame E_t(y) = y0(t) + c + Q(t)(y - c) the laboratory-frame field is

        w(t,x) = y0'(t) + th'(t) r^perp + Q(t) v(t, c + Q(t)^T r),
        r      = x - y0(t) - c,     v(t,y) = alpha(t) u0(y).

    Rather than evaluating w and interpolating it to faces, the *stream
    function* of w is assembled and the volumetric face fluxes are obtained as
    the discrete circulation of the vector potential Psi*e_y around each face
    loop.  That keeps sum_f phi_f identically zero (to round-off) for any
    frame, which is what the machine-precision volume conservation of
    isoAdvection relies on.  The lab-frame stream function is

        Psi(t,x) = y0'_z x - y0'_x z            (translation)
                 + 0.5 th'(t) |r|^2             (rigid rotation)
                 + alpha(t) psi0(c + Q^T r)     (base vortex, unchanged)

    the last term using that the 90-degree rotation operator commutes with
    Q in SO(2), so a rotated planar field keeps the same scalar stream
    function evaluated at the pulled-back point.

    The base stream function psi0 is extended by zero outside the unit square;
    psi0 and grad(psi0) both vanish on its boundary, so this extension is C1.

SourceFiles
    movingFrameFlow.C

\*---------------------------------------------------------------------------*/

#ifndef movingFrameFlow_H
#define movingFrameFlow_H

#include "fvMesh.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "dictionary.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class movingFrameFlow Declaration
\*---------------------------------------------------------------------------*/

class movingFrameFlow
{
    // Private data

        const fvMesh& mesh_;
        volVectorField& U_;
        surfaceScalarField& phi_;

        //- Rotation centre (x and z components used)
        vector centre_;

        //- Frame revolutions performed over endTime_
        scalar revolutions_;

        //- Recipe final time T; frame is periodic on [0,T]
        scalar endTime_;

        //- Base amplitude period: alpha(t) = cos(2 pi t/period_)
        scalar period_;

        //- Scaling of the base velocity field (0 -> frame motion only)
        scalar baseAmplitude_;

        //- Scaling of the frame translation offset (0 -> pure rotation)
        scalar translationAmplitude_;

        //- Stream function at mesh points
        scalarField Psi_;

        //- max|div(u)| from the most recent update
        scalar maxMagDivPhi_;


    // Private Member Functions

        scalar theta(const scalar t) const;
        scalar thetaDot() const;
        vector yOff(const scalar t) const;
        vector yOffDot(const scalar t) const;
        scalar amplitude(const scalar t) const;
        scalar psi0(const scalar x, const scalar z) const;
        vector u0(const scalar x, const scalar z) const;
        scalar faceCirculation(const label facei) const;

        //- No copy construct / assignment
        movingFrameFlow(const movingFrameFlow&) = delete;
        void operator=(const movingFrameFlow&) = delete;


public:

    // Constructors

        movingFrameFlow
        (
            const fvMesh& mesh,
            volVectorField& U,
            surfaceScalarField& phi,
            const dictionary& dict
        );


    // Member Functions

        //- Set phi and U for time t (pass the step mid-point)
        void update(const scalar t);

        //- max|div(u)| from the most recent update
        scalar maxMagDivPhi() const
        {
            return maxMagDivPhi_;
        }
};


} // End namespace Foam

#endif

// ************************************************************************* //
