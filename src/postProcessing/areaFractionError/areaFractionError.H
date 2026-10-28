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

Class
    Foam::functionObjects::areaFractionError

Description
    The error of the transported DAVOF face fractions alphaf.<group> against
    the exact face fractions of an implicit function translated uniformly,
    origin(t) = origin + U0 t, and the reconstruction error of the interface
    (normal and centre) against the same translated function.

    Every write: the exact DAVOF state of the translated function is built
    with davofInit (davofMethod: exactSphere on axis-aligned hexahedra,
    quadraticFaces or detrixheAslam) and, when alphaf.<group> is registered
    (the davof reconstruction scheme with isoAdvection), the face error
    |alphaf - alphaf_exact| is taken over the faces of the interface cells of
    either state (a cell is an interface cell when |m_c| > wispTol V^(2/3)
    with m_c = sum_f alphaf S_f^out, in the current or the exact state):
    mean (L1), root mean square (L2), maximum (Linf), the face count and the
    interface-cell count. With the davof scheme and consistency true (the
    default) the consistency of the transported state is written as
    max |m_c - A_c n_c| / V_c^(2/3) over the interface cells, A_c n_c the
    area vector of the plane the scheme positions in it, which costs one
    extra reconstruction of the scheme per write; consistency false skips
    that reconstruction and writes -1. Without alphaf the face columns
    are -1.

    Nothing is computed outside write(): execute() does nothing, and the
    exact state, the error sums and the extra reconstructions are built
    in write() only, so writeControl onEnd (or timeStep with a
    writeInterval that the end step meets) costs nothing before the write.
    A timeStep interval the end step does not divide writes NOTHING (the
    function object has no end() of its own): for the end only use
    writeControl onEnd. The skipped consistency reconstruction leaves the
    transport bit-identical for the models that reconstruct from alphaf
    and alpha alone (cellLocal, youngs, exact, ringArea, ringFit) and
    changes the plicRDF model's path (its helper starts its iteration from
    the previous normals, see reconstructionErrorEveryStep below).

    The reconstruction error (reconstructionError, the columns of the
    reconstruction benchmarks: LNormalDiff1/2/Inf, LCentreDiff1/2/Inf, ...)
    is written to postProcessing/reconstructionError/<time>/error.dat for
    the scheme reconstructed from the end-of-step state at the end time
    only (reconstructionErrorEveryStep true: at every write; plicRDF's
    iteration starts from the previous normals, so an extra reconstruction
    per step changes its path slightly).

    Dictionary:
        type        areaFractionError;
        libs        ("libpostProcess.so");
        field       alpha.water;
        functionType sphere;            // the implicitFunction type
        origin      (0.5 0.5 0.5);      // at t = 0
        radius      0.15;               // and the function's other keys
        U0          (1 0 0);            // default: U when uniform
        davofMethod exactSphere;        // exactSphere | quadraticFaces | detrixheAslam
        wispTol     1e-3;
        consistency true;               // the extra reconstruction per write; false writes -1
        reconstructionError true;
        reconstructionErrorEveryStep false;
        writeControl timeStep;          // onEnd: the end step only
        writeInterval 1;

    Output columns (postProcessing/<name>/<time>/areaFractionError.dat):
        Time L1 L2 Linf nFaces nInterfaceCells consistency

SourceFiles
    areaFractionError.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_areaFractionError_H
#define functionObjects_areaFractionError_H

#include "fvMeshFunctionObject.H"
#include "writeFile.H"
#include "volFieldsFwd.H"
#include "Switch.H"
#include "reconstructionError.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                      Class areaFractionError Declaration
\*---------------------------------------------------------------------------*/

class areaFractionError
:
    public fvMeshFunctionObject,
    public writeFile
{
    // Private Data

        //- The dictionary (the implicit function is built from a copy with
        //  the translated origin)
        dictionary dict_;

        word fieldName_;

        word functionType_;

        vector origin0_;

        vector U0_;

        //- U0 not given: taken from the velocity field when uniform
        bool U0FromField_;

        //- The rotation about an axis (keys omega, rad/s, default zero;
        //  axisPoint, default origin): the body's origin at t is
        //  davof::orbitOrigin (davofOrbit.H, 2026-10-28)
        vector omega_;
        point axisPoint_;

        word davofMethod_;

        scalar wispTol_;

        //- The consistency measure (one extra reconstruction per write)
        Switch consistency_;

        bool reconstructionError_;

        bool reconstructionErrorEveryStep_;

        //- Write the signed per-face error alphaf - alphaf_exact on the
        //  faces of the interface cells as the surface field
        //  alphafError.<group> at every write of this object (diagnostic,
        //  2026-10-16; default false)
        bool writeFaceError_;

        autoPtr<Foam::reconstructionError> recErrPtr_;


    // Private Member Functions

        //- File header information
        virtual void writeFileHeader(Ostream& os) const;

        //- No copy construct
        areaFractionError(const areaFractionError&) = delete;

        //- No copy assignment
        void operator=(const areaFractionError&) = delete;


public:

    //- Runtime type information
    TypeName("areaFractionError");


    // Constructors

        //- Construct from Time and dictionary
        areaFractionError
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~areaFractionError() = default;


    // Member Functions

        //- Read the areaFractionError data
        virtual bool read(const dictionary&);

        //- Do nothing
        virtual bool execute();

        //- Calculate and write the errors
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
