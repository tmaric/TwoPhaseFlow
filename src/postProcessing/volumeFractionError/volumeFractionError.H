/*---------------------------------------------------------------------------*\
            Copyright (c) 2017-2019, German Aerospace Center (DLR)
-------------------------------------------------------------------------------
License
    This file is part of the VoFLibrary source code library, which is an
    unofficial extension to OpenFOAM.

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

Class
    Foam::functionObjects::volumeFractionError


Description
    Eshape = sum_c |alpha_c - alpha_exact,c| V_c, the mass and the bound
    error of a volume fraction against an implicit function (functionType,
    origin at t = 0 and the function's keys) moved to its position at t by
    davof::orbitOrigin (davofOrbit.H): origin(t) = axisPoint
    + R(omega t)(origin - axisPoint) + U0 t. Keys omega (default zero) and
    axisPoint (default origin); U0 (2026-11-01): when given it is used; when
    absent it is zero with a rotation (omega nonzero) and, without a
    rotation, U of the first cell when |U| is uniform (a constant uniform
    translation), else zero.

Author
    Henning Scheufler, DLR, all rights reserved.
    Johan Roenby, DHI, all rights reserved.

SourceFiles
    volumeFractionError.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_volumeFractionError_H
#define functionObjects_volumeFractionError_H

#include "fvMeshFunctionObject.H"
#include "writeFile.H"
#include "volFieldsFwd.H"
#include "HashSet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                       Class volumeFractionError Declaration
\*---------------------------------------------------------------------------*/

class volumeFractionError
:
    public fvMeshFunctionObject,
    public writeFile
{

protected:

    // Protected data

        word Field_;

        dictionary dict_; //

        scalar initMass_;

        vector initCentre_;


    // Protected Member Functions

        //- File header information
        virtual void writeFileHeader(Ostream& os) const;



private:

    // Private member functions

        //- Disallow default bitwise copy construct
        volumeFractionError(const volumeFractionError&);

        //- Disallow default bitwise assignment
        void operator=(const volumeFractionError&);


public:

    //- Runtime type information
    TypeName("volumeFractionError");


    // Constructors

        //- Construct from Time and dictionary
        volumeFractionError
        (
            const word& name,
            const Time& runTime,
            const dictionary&
        );


    //- Destructor
    virtual ~volumeFractionError();


    // Member Functions

        //- Read the volumeFractionError data
        virtual bool read(const dictionary&);

        //- Calculate the wall heat-flux
        virtual bool execute();

        //- Write the wall heat-flux
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
