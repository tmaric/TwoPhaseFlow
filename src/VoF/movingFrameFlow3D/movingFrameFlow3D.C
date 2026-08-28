/*---------------------------------------------------------------------------*\
    Implementation of Foam::movingFrameFlow3D.  See movingFrameFlow3D.H.
\*---------------------------------------------------------------------------*/

#include "movingFrameFlow3D.H"
#include "fvcSurfaceIntegrate.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::movingFrameFlow3D::movingFrameFlow3D
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
    centre_(dict.lookupOrDefault<vector>("rotationCentre", vector(0.5, 0.5, 0))),
    swirl_(dict.lookupOrDefault<vector>("swirlCentre", vector(0.5, 0.5, 0))),
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
    A_(mesh.points().size(), vector(0, 0, 0)),
    maxMagDivPhi_(0.0)
{
    if (endTime_ <= 0.0)
    {
        FatalErrorInFunction
            << "movingFrameFlow3D: endTime must be positive, got " << endTime_
            << exit(FatalError);
    }

    Info<< "movingFrameFlow3D:" << nl
        << "    rotationCentre        " << centre_ << nl
        << "    swirlCentre           " << swirl_ << nl
        << "    revolutions           " << revolutions_ << nl
        << "    endTime (T)           " << endTime_ << nl
        << "    period                " << period_ << nl
        << "    baseAmplitude         " << baseAmplitude_ << nl
        << "    translationAmplitude  " << translationAmplitude_ << endl;
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::scalar Foam::movingFrameFlow3D::theta(const scalar t) const
{
    return 2.0*constant::mathematical::pi*revolutions_*t/endTime_;
}


Foam::scalar Foam::movingFrameFlow3D::thetaDot() const
{
    return 2.0*constant::mathematical::pi*revolutions_/endTime_;
}


Foam::vector Foam::movingFrameFlow3D::yOff(const scalar t) const
{
    if (translationAmplitude_ == 0.0)
    {
        return vector(0, 0, 0);
    }

    // Recipe Eq. (8) with u(s) = (s, s^2), v(s) = (sin s, sin 2s), applied in
    // the horizontal plane so that the axial drift is left untouched.
    const scalar s = t/endTime_;
    const scalar w = 1.0 - s;

    return translationAmplitude_*vector
    (
        (1.0 - s)*s   + s*Foam::sin(w),
        (1.0 - s)*s*s + s*Foam::sin(2.0*w),
        0
    );
}


Foam::vector Foam::movingFrameFlow3D::yOffDot(const scalar t) const
{
    if (translationAmplitude_ == 0.0)
    {
        return vector(0, 0, 0);
    }

    const scalar s = t/endTime_;
    const scalar w = 1.0 - s;

    return (translationAmplitude_/endTime_)*vector
    (
        -s + (1.0 - s)         + Foam::sin(w)     - s*Foam::cos(w),
        -s*s + (1.0 - s)*2.0*s + Foam::sin(2.0*w) - 2.0*s*Foam::cos(2.0*w),
        0
    );
}


Foam::scalar Foam::movingFrameFlow3D::amplitude(const scalar t) const
{
    if (period_ > 0.0)
    {
        return baseAmplitude_
            *Foam::cos(2.0*constant::mathematical::pi*t/period_);
    }

    return baseAmplitude_;
}


Foam::vector Foam::movingFrameFlow3D::A0(const vector& xi) const
{
    const scalar pi = constant::mathematical::pi;

    // horizontal part: the same stream function as the two-dimensional test,
    // extended by zero outside the unit square (psi0 and grad(psi0) both
    // vanish on its boundary, so the extension is C1)
    scalar psi0 = 0.0;
    if (xi.x() >= 0.0 && xi.x() <= 1.0 && xi.y() >= 0.0 && xi.y() <= 1.0)
    {
        const scalar sx = Foam::sin(pi*xi.x());
        const scalar sy = Foam::sin(pi*xi.y());
        psi0 = (1.0/pi)*sx*sx*sy*sy;
    }

    // axial part: A_theta(rho) e_theta with A_theta = rho f(rho); the 1/rho of
    // e_theta cancels, so this is regular on the axis
    const scalar dx = xi.x() - swirl_.x();
    const scalar dy = xi.y() - swirl_.y();
    const scalar rho = Foam::sqrt(dx*dx + dy*dy);
    const scalar f = 0.5 - (4.0/3.0)*rho + rho*rho;

    return vector(-f*dy, f*dx, psi0);
}


Foam::vector Foam::movingFrameFlow3D::u0(const vector& xi) const
{
    const scalar pi = constant::mathematical::pi;

    scalar ux = 0.0;
    scalar uy = 0.0;
    if (xi.x() >= 0.0 && xi.x() <= 1.0 && xi.y() >= 0.0 && xi.y() <= 1.0)
    {
        const scalar sx = Foam::sin(pi*xi.x());
        const scalar sy = Foam::sin(pi*xi.y());
        ux =  Foam::sin(2.0*pi*xi.y())*sx*sx;
        uy = -Foam::sin(2.0*pi*xi.x())*sy*sy;
    }

    const scalar dx = xi.x() - swirl_.x();
    const scalar dy = xi.y() - swirl_.y();
    const scalar rho = Foam::sqrt(dx*dx + dy*dy);
    const scalar uz = Foam::sqr(1.0 - 2.0*rho);

    return vector(ux, uy, uz);
}


Foam::scalar Foam::movingFrameFlow3D::faceCirculation(const label facei) const
{
    const pointField& pts = mesh_.points();
    const face& f = mesh_.faces()[facei];
    const label nPoints = f.size();

    scalar sum = 0.0;

    forAll(f, i)
    {
        const label p1 = f[i];
        const label p2 = f[(i + 1) % nPoints];

        sum += 0.5*((A_[p1] + A_[p2]) & (pts[p2] - pts[p1]));
    }

    return sum;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::movingFrameFlow3D::update(const scalar t)
{
    const scalar th  = theta(t);
    const scalar ct  = Foam::cos(th);
    const scalar st  = Foam::sin(th);
    const scalar thd = thetaDot();
    const vector yo  = yOff(t);
    const vector yd  = yOffDot(t);
    const scalar a   = amplitude(t);

    // 1) lab-frame vector potential at the mesh points
    const pointField& pts = mesh_.points();
    A_.setSize(pts.size());

    forAll(pts, pointi)
    {
        const point& p = pts[pointi];

        const scalar rx = p.x() - yo.x() - centre_.x();
        const scalar ry = p.y() - yo.y() - centre_.y();

        // xi = c + Q^T r, the rotation being about the axis parallel to e_z
        const vector xi
        (
            centre_.x() + ( rx*ct + ry*st),
            centre_.y() + (-rx*st + ry*ct),
            p.z() - yo.z()
        );

        // Q A0(xi): rotate the horizontal components, leave e_z invariant
        const vector Ab = A0(xi);
        const vector QAb
        (
            Ab.x()*ct - Ab.y()*st,
            Ab.x()*st + Ab.y()*ct,
            Ab.z()
        );

        // rigid rotation: A = -1/2 |r_perp|^2 omega, omega = thetaDot e_z
        const vector Arot(0, 0, -0.5*thd*(rx*rx + ry*ry));

        // translation: A = 1/2 (ydot x p)
        const vector Atr = 0.5*(yd ^ p);

        A_[pointi] = a*QAb + Arot + Atr;
    }

    // 2) fluxes as the discrete circulation of A around each face loop; every
    //    edge of a cell is traversed twice in opposite senses, so the discrete
    //    divergence vanishes to round-off for any A
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

    // 3) cell-centred velocity, evaluated analytically (Courant/output only)
    vectorField& Uin = U_.primitiveFieldRef();
    forAll(Uin, celli)
    {
        const point& p = mesh_.C()[celli];

        const scalar rx = p.x() - yo.x() - centre_.x();
        const scalar ry = p.y() - yo.y() - centre_.y();

        const vector xi
        (
            centre_.x() + ( rx*ct + ry*st),
            centre_.y() + (-rx*st + ry*ct),
            p.z() - yo.z()
        );

        const vector ub = a*u0(xi);
        const vector Qub
        (
            ub.x()*ct - ub.y()*st,
            ub.x()*st + ub.y()*ct,
            ub.z()
        );

        Uin[celli] = yd + thd*vector(-ry, rx, 0) + Qub;
    }

    U_.correctBoundaryConditions();

    // 4) discrete solenoidality diagnostic
    const volScalarField divU(fvc::surfaceIntegrate(phi_));
    maxMagDivPhi_ = gMax(mag(divU.primitiveField()));
}


// ************************************************************************* //
