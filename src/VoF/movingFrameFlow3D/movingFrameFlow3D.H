/*---------------------------------------------------------------------------*\
Class
    Foam::movingFrameFlow3D

Description
    Three-dimensional counterpart of Foam::movingFrameFlow: a T-periodic
    Euclidean frame superposed on the spiralling deformation field of the
    \c test-deformationSphere / \c test-Leveque benchmarks, following

        Bothe & Maric, "A Recipe for Periodic, Non-Retracing Deformation Test
        Cases" (2026), Sections 3 and 6.

    The base field is

        u = sin(2 pi y) sin^2(pi x)
        v = -sin(2 pi x) sin^2(pi y)
        w = (1 - 2 rho)^2,   rho = |(x,y) - s|,  s the swirl axis,

    a horizontal, two-dimensionally solenoidal part plus a purely axial part.
    Both admit a closed-form vector potential, so the whole field does:

        A0 = psi0(x,y) e_z + f(rho) ( e_z x (p_h - s) ),
        psi0 = (1/pi) sin^2(pi x) sin^2(pi y),
        f(rho) = 1/2 - 4 rho/3 + rho^2,

    since (curl A0)_z = (1/rho) d/drho (rho^2/2 - 4 rho^3/3 + rho^4)
                      = 1 - 4 rho + 4 rho^2 = (1 - 2 rho)^2,
    and the 1/rho cancels, so A0 is regular on the axis.

    THIS REPLACES A PRESSURE PROJECTION. The reference generator builds the
    flux as a face-centre quadrature U_f & S_f, which is NOT discretely
    solenoidal, and then repairs it with CorrectPhi -- a Poisson solve that
    also perturbs the field away from the analytic one. Here the fluxes are the
    discrete circulation of A around each face loop,

        phi_f = sum_edges 1/2 (A_p1 + A_p2) . (x_p2 - x_p1),

    which telescopes to zero over the faces of any closed cell. Discrete
    solenoidality therefore holds to round-off by construction, with no linear
    solve and with the advected field exactly the analytic one -- which is what
    makes a moving frame affordable at all, since under a frame the flux
    changes shape every step rather than merely being rescaled.

    The frame rotates about an axis parallel to e_z through \c rotationCentre.
    That choice is deliberate: the horizontal motion is confined to level sets
    of psi0, so the material support is a streamline band exactly as in two
    dimensions, while a z-parallel rotation leaves the axial drift untouched.

SourceFiles
    movingFrameFlow3D.C

\*---------------------------------------------------------------------------*/

#ifndef movingFrameFlow3D_H
#define movingFrameFlow3D_H

#include "fvMesh.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "dictionary.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class movingFrameFlow3D Declaration
\*---------------------------------------------------------------------------*/

class movingFrameFlow3D
{
    // Private data

        const fvMesh& mesh_;
        volVectorField& U_;
        surfaceScalarField& phi_;

        //- Which base field: "spiralling" (2D vortex + axial) or "leveque"
        //  (the fully three-dimensional deformation field)
        word baseField_;

        //- Point on the frame rotation axis (axis is parallel to e_z)
        vector centre_;

        //- Swirl axis of the base field
        vector swirl_;

        //- Frame revolutions performed over endTime_
        scalar revolutions_;

        //- Recipe final time T
        scalar endTime_;

        //- Base amplitude period: alpha(t) = cos(2 pi t/period_).
        //  The reference case is spelt "spirallingFlow S", i.e.
        //  alpha = cos(pi t/S), so period_ = 2 S.
        scalar period_;

        //- Scaling of the base field (0 -> frame motion only)
        scalar baseAmplitude_;

        //- Scaling of the frame translation offset
        scalar translationAmplitude_;

        //- Vector potential at the mesh points
        vectorField A_;

        //- max|div(u)| from the most recent update
        scalar maxMagDivPhi_;


    // Private Member Functions

        scalar theta(const scalar t) const;
        scalar thetaDot() const;
        vector yOff(const scalar t) const;
        vector yOffDot(const scalar t) const;
        scalar amplitude(const scalar t) const;

        //- Base vector potential at a point
        vector A0(const vector& xi) const;

        //- Base velocity at a point
        vector u0(const vector& xi) const;

        //- Discrete circulation of A_ around one face loop
        scalar faceCirculation(const label facei) const;

        //- No copy construct / assignment
        movingFrameFlow3D(const movingFrameFlow3D&) = delete;
        void operator=(const movingFrameFlow3D&) = delete;


public:

    // Constructors

        movingFrameFlow3D
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
