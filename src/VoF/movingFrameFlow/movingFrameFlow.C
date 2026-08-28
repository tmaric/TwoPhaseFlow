/*---------------------------------------------------------------------------*\
    Implementation of Foam::movingFrameFlow.  See movingFrameFlow.H.
\*---------------------------------------------------------------------------*/

#include "movingFrameFlow.H"
#include "fvcSurfaceIntegrate.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::movingFrameFlow::movingFrameFlow
(
    const fvMesh& mesh,
    volVectorField& U,
    surfaceScalarField& phi,
    const dictionary& dict
)
:
    mesh_(mesh),
    U_(U),
    phi_(phi),
    centre_(dict.lookupOrDefault<vector>("rotationCentre", vector(0.5, 0, 0.5))),
    revolutions_(dict.lookupOrDefault<scalar>("revolutions", 1.0)),
    endTime_
    (
        dict.lookupOrDefault<scalar>("endTime", mesh.time().endTime().value())
    ),
    period_(dict.lookupOrDefault<scalar>("period", 0.0)),
    baseAmplitude_(dict.lookupOrDefault<scalar>("baseAmplitude", 1.0)),
    translationAmplitude_
    (
        dict.lookupOrDefault<scalar>("translationAmplitude", 0.0)
    ),
    Psi_(mesh.points().size(), 0.0),
    maxMagDivPhi_(0.0)
{
    if (endTime_ <= 0.0)
    {
        FatalErrorInFunction
            << "movingFrameFlow: endTime must be positive, got " << endTime_
            << exit(FatalError);
    }

    Info<< "movingFrameFlow:" << nl
        << "    rotationCentre        " << centre_ << nl
        << "    revolutions           " << revolutions_ << nl
        << "    endTime (T)           " << endTime_ << nl
        << "    period                " << period_ << nl
        << "    baseAmplitude         " << baseAmplitude_ << nl
        << "    translationAmplitude  " << translationAmplitude_ << endl;
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::scalar Foam::movingFrameFlow::theta(const scalar t) const
{
    return 2.0*constant::mathematical::pi*revolutions_*t/endTime_;
}


Foam::scalar Foam::movingFrameFlow::thetaDot() const
{
    return 2.0*constant::mathematical::pi*revolutions_/endTime_;
}


Foam::vector Foam::movingFrameFlow::yOff(const scalar t) const
{
    if (translationAmplitude_ == 0.0)
    {
        return vector(0, 0, 0);
    }

    // Recipe Eq. (8) with u(s) = (s, s^2), v(s) = (sin s, sin 2s)
    const scalar s = t/endTime_;
    const scalar w = 1.0 - s;

    return translationAmplitude_*vector
    (
        (1.0 - s)*s     + s*Foam::sin(w),
        0,
        (1.0 - s)*s*s   + s*Foam::sin(2.0*w)
    );
}


Foam::vector Foam::movingFrameFlow::yOffDot(const scalar t) const
{
    if (translationAmplitude_ == 0.0)
    {
        return vector(0, 0, 0);
    }

    const scalar s = t/endTime_;
    const scalar w = 1.0 - s;

    return (translationAmplitude_/endTime_)*vector
    (
        -s + (1.0 - s)           + Foam::sin(w)        - s*Foam::cos(w),
        0,
        -s*s + (1.0 - s)*2.0*s   + Foam::sin(2.0*w)    - 2.0*s*Foam::cos(2.0*w)
    );
}


Foam::scalar Foam::movingFrameFlow::amplitude(const scalar t) const
{
    if (period_ > 0.0)
    {
        return baseAmplitude_
            *Foam::cos(2.0*constant::mathematical::pi*t/period_);
    }

    return baseAmplitude_;
}


Foam::scalar Foam::movingFrameFlow::psi0(const scalar x, const scalar z) const
{
    // C1 extension by zero: psi0 and grad(psi0) vanish on the unit square
    if (x < 0.0 || x > 1.0 || z < 0.0 || z > 1.0)
    {
        return 0.0;
    }

    const scalar sx = Foam::sin(constant::mathematical::pi*x);
    const scalar sz = Foam::sin(constant::mathematical::pi*z);

    return (1.0/constant::mathematical::pi)*sx*sx*sz*sz;
}


Foam::vector Foam::movingFrameFlow::u0(const scalar x, const scalar z) const
{
    if (x < 0.0 || x > 1.0 || z < 0.0 || z > 1.0)
    {
        return vector(0, 0, 0);
    }

    const scalar sx = Foam::sin(constant::mathematical::pi*x);
    const scalar sz = Foam::sin(constant::mathematical::pi*z);

    return vector
    (
        -sx*sx*Foam::sin(2.0*constant::mathematical::pi*z),
        0,
        Foam::sin(2.0*constant::mathematical::pi*x)*sz*sz
    );
}


Foam::scalar Foam::movingFrameFlow::faceCirculation(const label facei) const
{
    const pointField& pts = mesh_.points();
    const face& f = mesh_.faces()[facei];
    const label nPoints = f.size();

    scalar sum = 0.0;

    forAll(f, i)
    {
        const label p1 = f[i];
        const label p2 = f[(i + 1) % nPoints];

        sum += 0.5*(Psi_[p1] + Psi_[p2])*(pts[p2].y() - pts[p1].y());
    }

    return sum;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::movingFrameFlow::update(const scalar t)
{
    const scalar th  = theta(t);
    const scalar ct  = Foam::cos(th);
    const scalar st  = Foam::sin(th);
    const scalar thd = thetaDot();
    const vector yo  = yOff(t);
    const vector yd  = yOffDot(t);
    const scalar a   = amplitude(t);

    // 1) Lab-frame stream function at the mesh points
    const pointField& pts = mesh_.points();
    Psi_.setSize(pts.size());

    forAll(pts, pointi)
    {
        const point& p = pts[pointi];

        const scalar rx = p.x() - yo.x() - centre_.x();
        const scalar rz = p.z() - yo.z() - centre_.z();

        // xi = c + Q^T r
        const scalar xix = centre_.x() + ( rx*ct + rz*st);
        const scalar xiz = centre_.z() + (-rx*st + rz*ct);

        Psi_[pointi] =
            yd.z()*p.x() - yd.x()*p.z()     // translation
          + 0.5*thd*(rx*rx + rz*rz)         // rigid rotation
          + a*psi0(xix, xiz);               // base vortex in frame coords
    }

    // 2) Fluxes as the discrete circulation of Psi*e_y around each face loop.
    //    Every edge is traversed twice with opposite sense when a cell is
    //    closed, so sum_f phi_f vanishes to round-off for any Psi.
    scalarField& phiIn = phi_.primitiveFieldRef();
    forAll(phiIn, facei)
    {
        phiIn[facei] = faceCirculation(facei);
    }

    surfaceScalarField::Boundary& phiBf = phi_.boundaryFieldRef();
    forAll(phiBf, patchi)
    {
        fvsPatchScalarField& pphi = phiBf[patchi];
        const label start = mesh_.boundaryMesh()[patchi].start();

        forAll(pphi, i)
        {
            pphi[i] = faceCirculation(start + i);
        }
    }

    // 3) Cell-centred velocity, evaluated analytically (Courant number/output
    //    only -- the advection itself uses phi)
    vectorField& Uin = U_.primitiveFieldRef();
    forAll(Uin, celli)
    {
        const point& p = mesh_.C()[celli];

        const scalar rx = p.x() - yo.x() - centre_.x();
        const scalar rz = p.z() - yo.z() - centre_.z();

        const scalar xix = centre_.x() + ( rx*ct + rz*st);
        const scalar xiz = centre_.z() + (-rx*st + rz*ct);

        const vector ub = a*u0(xix, xiz);
        const vector Qub(ub.x()*ct - ub.z()*st, 0, ub.x()*st + ub.z()*ct);

        Uin[celli] = yd + thd*vector(-rz, 0, rx) + Qub;
    }

    U_.correctBoundaryConditions();

    // 4) Discrete solenoidality diagnostic
    const volScalarField divU(fvc::surfaceIntegrate(phi_));
    maxMagDivPhi_ = gMax(mag(divU.primitiveField()));
}


// ************************************************************************* //
