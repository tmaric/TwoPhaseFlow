/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2024 AUTHOR, AFFILIATION
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
    Foam::functionObjects::vofForcesFunctionObject

Group
    grpFieldFunctionObjects

Description
    <minimal description of the function object>

    <equation>
    \f[
        x = x_{ref}^x + \rho \omega
    \f]

    <variable-explanation table>
    where
    \vartable
      \rho        | <explanation> [units, e.g. kg/m3]
      \omega      | \f$ \nabla \cdot \vec U \f$
      ...         | ...
    \endvartable

    <inline equation>
    where \f$ x_k \f$ is ...

    <input-output table>
    \table
      Operand      | Type                 | Location
      input        | {vol,surface}\<Type\>Field(s) <!--
               --> |$FOAM_CASE/\<time\>/\<inpField\>s
      output file  | dat                  <!--
              -->  | $FOAM_CASE/postProcessing/\<FO\>/\<time\>/\<file\>
      output field | volScalarField       | $FOAM_CASE/\<time\>/\<outField\>
    \endtable

Usage
    Minimal example by using \c system/controlDict.functions:
    \verbatim
    vofForcesFunctionObject1
    {
        // Mandatory entries (unmodifiable)
        type           vofForcesFunctionObject;
        libs           (vofForcesFunctionObjectFunctionObject);

        // Mandatory entries (runtime modifiable)
        ...

        // Mandatory (inherited) entries (unmodifiable)
        ...

        // Mandatory (inherited) entries (runtime unmodifiable)
        ...

        // Optional entries (unmodifiable)
        ...

        // Optional entries (runtime modifiable)
        boolData       <bool>;
        labelData      <label>;
        wordData       <word>;
        scalarData     <scalar>;

        // Optional (inherited) entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property   | Description                      | Type | Req'd | Dflt
      type       | Type name: vofForcesFunctionObject        | word |  yes  | -
      libs       | Library name: vofForcesFunctionObjectFunctionObject  <!--
             --> | word  | yes | -
      boolData   | <explanation>                    | bool |  yes  | -
      labelData  | <explanation>                    | label |  yes | -
      wordData   | <explanation>                    | word |  yes  | -
      scalarData | <explanation>                    | scalar |  no | 1.0
      wordListData  | <explanation>                 | wordList | yes | -
    \endtable

    Options for the \c ENTRY entry:
    \verbatim
        <option1>
        <option2> | <explanation>
      ...
    \endverbatim

    The inherited entries are elaborated in:
     - \link functionObject.H \endlink
     - \link fieldExpression.H \endlink
     - \link fieldsExpression.H \endlink
     - \link writeFile.H \endlink
     ...

    <if \c postProcess is applicable>
    Minimal example by using the \c postProcess utility:
    \verbatim
        postProcess -func vofForcesFunctionObject
    \endverbatim

    <if \c postProcess is not applicable>
    Usage by the \c postProcess utility is not available.

Note
    - <note1>
    - <note2>
    ...

See also
    - Foam::functionObject
    - Foam::functionObjects::fvMeshFunctionObject
    - ExtendedCodeGuide::functionObjects::field::vofForcesFunctionObject
    ...

SourceFiles
    vofForcesFunctionObject.C
    vofForcesFunctionObjectTEMPLATES.C
    ...

\*---------------------------------------------------------------------------*/

#ifndef vofForcesFunctionObject_H
#define vofForcesFunctionObject_H

#include "fvMeshFunctionObject.H"
#include "volFields.H"
#include "immiscibleIncompressibleTwoPhaseMixture.H"
#include "surfaceForces.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                   Class vofForcesFunctionObject Declaration
\*---------------------------------------------------------------------------*/

class vofForcesFunctionObject
:
    public fvMeshFunctionObject
{
    const volScalarField& rho_; 
    const volVectorField& U_; 

    const volScalarField& p_rgh_;
    const surfaceScalarField& phi_;
    const surfaceScalarField& rhoPhi_;


    // Interface Properties
    immiscibleIncompressibleTwoPhaseMixture mixture_; 

    const volScalarField& alpha1_;

    volScalarField mu_;

    // Surface tension force
    surfaceForces surfForces_; 

    // Terms from the momentum equation 

    volVectorField ddtRhoU_; 

    volVectorField inertialForce_; 

    volVectorField viscousForce_;  

    volVectorField surfaceTensionForce_;

    volVectorField accelerationForce_;

    volVectorField pressureForce_;

public:

    //- Runtime type information
    TypeName("vofForcesFunctionObject");


    // Constructors

        //- Construct from Time and dictionary
        vofForcesFunctionObject
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- No copy construct
        vofForcesFunctionObject(const vofForcesFunctionObject&) = delete;

        //- No copy assignment
        void operator=(const vofForcesFunctionObject&) = delete;


    //- Destructor
    virtual ~vofForcesFunctionObject() = default;


    // Member Functions

        //- Read the vofForcesFunctionObject data
        virtual bool read(const dictionary& dict);

        //- Execute, currently does nothing
        virtual bool execute();

        //- Execute at the final time-loop, currently does nothing
        virtual bool end();

        //- Write the vofForcesFunctionObject
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
