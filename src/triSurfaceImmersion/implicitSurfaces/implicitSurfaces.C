/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020 2020 Tomislav Maric, TU Darmstadt
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

\*---------------------------------------------------------------------------*/

#include "implicitSurfaces.H"

#include "Ostream.H"
#include "addToRunTimeSelectionTable.H"
#include "error.H"
#include "mathematicalConstants.H"
#include <quaternion.H>

namespace Foam::TriSurfaceImmersion
{
// * * * * * * * * * * * * Class implicitSurface  * * * * * * * * * * * //

// * * * * * * * * * * * * * * Static Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(implicitSurface, false);
defineRunTimeSelectionTable(implicitSurface, ITstream);
defineRunTimeSelectionTable(implicitSurface, Dictionary);

autoPtr<implicitSurface> implicitSurface::New(const word& name, ITstream& is)
{
    auto* ctorPtr = ITstreamConstructorTable(name);

    if (!ctorPtr)
    {
        FatalIOErrorInLookup
        (
            is,
            "implicitSurface",
            name,
            *ITstreamConstructorTablePtr_
        ) << exit(FatalIOError);
    }

    return autoPtr<implicitSurface>(ctorPtr(is));
}


autoPtr<implicitSurface> implicitSurface::New(const dictionary& configDict)
{
    const auto name = configDict.get<word>("surface");

    auto* ctorPtr = DictionaryConstructorTable(name);

    if (!ctorPtr)
    {
        FatalIOErrorInLookup
        (
            configDict,
            "implicitSurface",
            name,
            *DictionaryConstructorTablePtr_
        ) << exit(FatalIOError);
    }

    return autoPtr<implicitSurface>(ctorPtr(configDict));
}


scalar implicitSurface::orientation(const word& keyword)
{
    if (keyword == "outward")
    {
        return 1.0;
    }
    else if (keyword == "inward")
    {
        return -1.0;
    }
    else
    {
        FatalErrorIn("Foam::TriSurfaceImmersion::implicitSurface::"
                             "orientation(const word& keyword)")
            << "Unknown keyword " << keyword << nl << nl
            << "Valid keyword are 'outward' (no change in orientation) "
            << "and 'inward' (flip orientation)"
            << exit(FatalError);

        // Disable compiler warning, statement never reached (TT)
        return 0.0;
    }
}

// Only makes sense for closed surfaces, but needs to be supported by the
// interface class (TT)
scalar implicitSurface::volume() const
{
    return -SMALL;
}


// * * * * * * * * * * * * Class plane  * * * * * * * * * * * //

// * * * * * * * * * * * * * * Static Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(plane, false);
addToRunTimeSelectionTable(implicitSurface, plane, ITstream);
addToRunTimeSelectionTable(implicitSurface, plane, Dictionary);

// * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

plane::plane(vector position, vector normal)
    : position_{position}, normal_{normal}
{
    normal_ /= Foam::mag(normal_);
}


plane::plane(ITstream& is)
{
    is >> position_;
    is >> normal_;
}


plane::plane(const dictionary& configDict)
    : position_{configDict.get<vector>("position")}, normal_{
                                                         configDict.get<vector>(
                                                             "normal")}
{
    normal_ /= Foam::mag(normal_);
}


// * * * * * * * * * * * * * Member Functions * * * * * * * * * * * * * //

scalar plane::value(const vector& x) const
{
    return Foam::dot(x - position_, normal_);
}


scalar plane::operator()(const vector& x) const
{
    return value(x);
}


vector plane::grad(const vector&) const
{
    return normal_;
}


vector plane::position() const
{
    return position_;
}


vector plane::normal() const
{
    return normal_;
}

// * * * * * * * * * * * * Class sphere * * * * * * * * * * * //

// * * * * * * * * * * * * * * Static Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(sphere, false);
addToRunTimeSelectionTable(implicitSurface, sphere, ITstream);
addToRunTimeSelectionTable(implicitSurface, sphere, Dictionary);


// * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

sphere::sphere(vector center, scalar radius, scalar orientation)
    : center_{center}, radius_{radius}, orientation_{orientation}
{
}


sphere::sphere(ITstream& is)
{
    is >> center_;
    is >> radius_;
    is >> orientation_;
}


sphere::sphere(const dictionary& configDict)
    : center_{configDict.get<vector>("center")}, radius_{configDict.get<scalar>(
                                                     "radius")}
{
    auto keyword = configDict.getOrDefault<word>("orientation", "outward");
    orientation_ = orientation(keyword);
}


// * * * * * * * * * * * * * Member Functions * * * * * * * * * * * * * //

scalar sphere::value(const vector& x) const
{
    return orientation_*(magSqr(x - center_) - radius_*radius_);
}


scalar sphere::operator()(const vector& x) const
{
    return value(x);
}


vector sphere::grad(const vector& x) const
{
    return orientation_*2.0*(x - center_);
}


scalar sphere::volume() const
{
    return 4.0 / 3.0 * Foam::constant::mathematical::pi * pow(radius_, 3.0);
}


vector sphere::center() const
{
    return center_;
}


scalar sphere::radius() const
{
    return radius_;
}


// * * * * * * * * * * * * Class cylinder * * * * * * * * * * * //

// * * * * * * * * * * * * * * Static Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(cylinder, false);
addToRunTimeSelectionTable(implicitSurface, cylinder, ITstream);
addToRunTimeSelectionTable(implicitSurface, cylinder, Dictionary);


// * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

cylinder::cylinder(vector axis, vector pointOnAxis, scalar radius, scalar height, scalar orientation)
    : axis_{axis}, pointOnAxis_{pointOnAxis},radius_{radius}, height_{height}, orientation_{orientation}
{
    initializeRotation();
}


cylinder::cylinder(ITstream& is)
{
    is >> axis_;
    is >> pointOnAxis_;
    is >> radius_;
    is >> height_;
    is >> orientation_;

    initializeRotation();
}


cylinder::cylinder(const dictionary& configDict)
    : axis_{configDict.get<vector>("axis")},
      pointOnAxis_{configDict.get<vector>("pointOnAxis")},
      radius_{configDict.get<scalar>("radius")},
      height_{configDict.get<scalar>("height")}
{
    auto keyword = configDict.getOrDefault<word>("orientation", "outward");
    orientation_ = orientation(keyword);
    axis_ /= mag(axis_);
    initializeRotation();
}


// * * * * * * * * * * * * * Member Functions * * * * * * * * * * * * * //
void cylinder::initializeRotation()
{
    // In the cylinder aligned coordinate system the z-axis aligns with
    // the cylinder centre axis
    auto rotationAxis = axis_ ^ vector{0,0,1};

    // Catch case when the cylinder axis is parallel to the z-direction
    if (mag(rotationAxis) > SMALL)
    {
        rotationAxis /= mag(rotationAxis);
        auto rotationAngle = Foam::acos(vector{0,0,1} & axis_);
        rotation_ = quaternion{rotationAxis, rotationAngle};
    }
    else
    {
        rotation_ = quaternion{vector{0,0,1}, 0.0};
    }
}


vector cylinder::transformToCylinderKOS(const vector& x) const
{
    // The cylinder with infinite extension in z-direction is essentially
    // a 2D problem. Thus, the z-components of the transformed x is
    // disregarded.
    auto localX = rotation_.transform(x - pointOnAxis_);
    localX.z() = 0.0;
    return localX;
}


scalar cylinder::value(const vector& x) const
{
    auto localX = transformToCylinderKOS(x);
    return orientation_*(magSqr(localX) - radius_*radius_);
}


scalar cylinder::operator()(const vector& x) const
{
    return value(x);
}


vector cylinder::grad(const vector& x) const
{
    auto localX = transformToCylinderKOS(x);
    return rotation_.invTransform(orientation_*2.0*localX);
}


scalar cylinder::volume() const
{
    return Foam::constant::mathematical::pi * pow(radius_, 2.0) * height_;
}


vector cylinder::axis() const
{
    return axis_;
}


vector cylinder::pointOnAxis() const
{
    return pointOnAxis_;
}


scalar cylinder::radius() const
{
    return radius_;
}


scalar cylinder::height() const
{
    return height_;
}


// * * * * * * * * * * * * Class ellipsoid * * * * * * * * * * * //

// * * * * * * * * * * * * * * Static Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(ellipsoid, false);
addToRunTimeSelectionTable(implicitSurface, ellipsoid, ITstream);
addToRunTimeSelectionTable(implicitSurface, ellipsoid, Dictionary);


// * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

ellipsoid::ellipsoid(vector center, vector axes, scalar orientation)
    : center_{center}, axes_{axes}, orientation_{orientation}
{
    setAxesSqr(axes);
}


ellipsoid::ellipsoid(ITstream& is)
{
    is >> center_;
    is >> axes_;
    is >> orientation_;
    setAxesSqr(axes_);
}


ellipsoid::ellipsoid(const dictionary& configDict)
    : center_{configDict.get<vector>("center")}, axes_{configDict.get<vector>(
                                                     "axes")}
{
    auto keyword = configDict.getOrDefault<word>("orientation", "outward");
    orientation_ = orientation(keyword);
    setAxesSqr(axes_);
}


// * * * * * * * * * * * * * Member Functions * * * * * * * * * * * * * //

void ellipsoid::setAxesSqr(const vector& axes)
{
    axesSqr_[0] = Foam::sqr(axes[0]);
    axesSqr_[1] = Foam::sqr(axes[1]);
    axesSqr_[2] = Foam::sqr(axes[2]);
}


scalar ellipsoid::value(const vector& x) const
{
    return orientation_*(Foam::sqr(x[0] - center_[0]) / axesSqr_[0] +
        Foam::sqr(x[1] - center_[1]) / axesSqr_[1] +
        Foam::sqr(x[2] - center_[2]) / axesSqr_[2] - 1);
}


scalar ellipsoid::operator()(const vector& x) const
{
    return value(x);
}


vector ellipsoid::grad(const vector& x) const
{
    return orientation_*2*
        vector((x[0] - center_[0]) / axesSqr_[0],
            (x[1] - center_[1]) / axesSqr_[1],
            (x[2] - center_[2]) / axesSqr_[2]);
}


scalar ellipsoid::volume() const
{
    return 4.0 / 3.0 * Foam::constant::mathematical::pi * axes_[0] * axes_[1] *
        axes_[2];
}


vector ellipsoid::center() const
{
    return center_;
}


vector ellipsoid::axes() const
{
    return axes_;
}

// * * * * * * * * * * * * Class sinc * * * * * * * * * * * //

// * * * * * * * * * * * * * * Static Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(sinc, false);
addToRunTimeSelectionTable(implicitSurface, sinc, ITstream);
addToRunTimeSelectionTable(implicitSurface, sinc, Dictionary);


// * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

sinc::sinc(vector origin, scalar amplitude, scalar omega, scalar orientation)
    : origin_{origin}, amplitude_{amplitude}, omega_{omega},
      orientation_{orientation}
{
}


sinc::sinc(ITstream& is)
{
    is >> origin_;
    is >> amplitude_;
    is >> omega_;
    is >> orientation_;
}


sinc::sinc(const dictionary& configDict)
    : origin_{configDict.get<vector>("origin")},
      amplitude_{configDict.get<scalar>("amplitude")},
      omega_{configDict.get<scalar>("omega")}
{
    auto keyword = configDict.getOrDefault<word>("orientation", "outward");
    orientation_ = orientation(keyword);
}


// * * * * * * * * * * * * * Member Functions * * * * * * * * * * * * * //

scalar sinc::value(const vector& x) const
{
    double r =
        Foam::sqrt(Foam::sqr(x[0] - origin_[0]) + Foam::sqr(x[1] - origin_[1]));

    if (r < std::numeric_limits<double>::min())
    {
        return amplitude_;
    }

    return orientation_*
        (x[2] - origin_[2] - amplitude_ * sin(omega_ * r) / (omega_ * r));
}


scalar sinc::operator()(const vector& x) const
{
    return value(x);
}


vector sinc::grad(const vector& x) const
{
    const scalar& A = amplitude_;
    const scalar& O0 = origin_[0];
    const scalar& O1 = origin_[1];

    const scalar& x0 = x[0];
    const scalar& x1 = x[1];

    return orientation_*vector // Expression calculated in sympy.
        (A * (O0 - x0) *
                (omega_ * pow(pow(O0 - x0, 2) + pow(O1 - x1, 2), 3.0 / 2.0) *
                        cos(omega_ * sqrt(pow(O0 - x0, 2) + pow(O1 - x1, 2))) -
                    (pow(O0 - x0, 2) + pow(O1 - x1, 2)) *
                        sin(omega_ * sqrt(pow(O0 - x0, 2) + pow(O1 - x1, 2)))) /
                (omega_ * pow(pow(O0 - x0, 2) + pow(O1 - x1, 2), 5.0 / 2.0)),
            A * (O1 - x1) *
                (omega_ * pow(pow(O0 - x0, 2) + pow(O1 - x1, 2), 3.0 / 2.0) *
                        cos(omega_ * sqrt(pow(O0 - x0, 2) + pow(O1 - x1, 2))) -
                    (pow(O0 - x0, 2) + pow(O1 - x1, 2)) *
                        sin(omega_ * sqrt(pow(O0 - x0, 2) + pow(O1 - x1, 2)))) /
                (omega_ * pow(pow(O0 - x0, 2) + pow(O1 - x1, 2), 5.0 / 2.0)),
            1);
}


vector sinc::origin() const
{
    return origin_;
}


scalar sinc::amplitude() const
{
    return amplitude_;
}


scalar sinc::omega() const
{
    return omega_;
}

/*
// * * * * * * * * * * * * Class sincScaled * * * * * * * * * * * //

// * * * * * * * * * * * * * * Static Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(sincScaled, false);
addToRunTimeSelectionTable(implicitSurface, sincScaled, ITstream);

// * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

sincScaled::sincScaled(vector origin, scalar amplitude, scalar omega)
    : origin_{origin}, amplitude_{amplitude}, omega_{omega}
{
}

sincScaled::sincScaled(ITstream& is)
{
    is >> origin_;
    is >> amplitude_;
    is >> omega_;
}

sincScaled::sincScaled(const dictionary& configDict)
    : origin_{configDict.get<vector>("origin")},
      amplitude_{configDict.get<scalar>("amplitude")},
      omega_{configDict.get<scalar>("omega")}
{
}

// * * * * * * * * * * * * * Member Functions * * * * * * * * * * * * * //

scalar sincScaled::value(const vector& x) const // TODO (TM) scale the amplitude
{
    double r =
        Foam::sqrt(Foam::sqr(x[0] - origin_[0]) + Foam::sqr(x[1] - origin_[1]));

    if (r < std::numeric_limits<double>::min())
    {
        return amplitude_;
    }

    return x[2] - origin_[2] - amplitude_ * sin(omega_ * r) / (omega_ * r);
}

scalar sincScaled::operator()(const vector& x) const
{
    return value(x);
}

vector sincScaled::grad(const vector& x) const
{
    // const scalar& A = amplitude_;
    ////const scalar& omega = omega_;
    // const scalar& O0 = origin_[0];
    // const scalar& O1 = origin_[1];

    // const scalar& x0 = x[0];
    // const scalar& x1 = x[1];

    // FIXME: insert sympy expression.
    return vector(GREAT, GREAT, GREAT);
}

vector sincScaled::origin() const
{
    return origin_;
}

scalar sincScaled::amplitude() const
{
    return amplitude_;
}

scalar sincScaled::omega() const
{
    return omega_;
}
*/

} // End namespace Foam::TriSurfaceImmersion
