/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
     \\/     M anipulation  |
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
    Foam::laminarBLInletFvPatchVectorField

Group
    grpInletBoundaryConditions

Description
    Laminar Boundary Layer Velocity profile based on the inlet flow speed
    and boundary layer height. A parabolic profile is used to prescribe
    velocity inside the boundary layer.

    \heading Patch usage

    \table
        Property     | Description             | Required    | Default value
        flowSpeed    | inlet velocity          | yes         |
        delta        | inlet boundary layer height | yes     |
    \endtable

    Example of the boundary condition specification:
    \verbatim
    myPatch
    {
        type            laminarBLInlet;
        flowSpeed       (0.1 0 0);          // free-stream velocity
        delta           0.1;                // BL height
        value           uniform (0 0 0);    // initial value
    }
    \endverbatim

SeeAlso
    Foam::fixedValueFvPatchField

SourceFiles
    laminarBLInletFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef laminarBLInletFvPatchVectorField_H
#define laminarBLInletFvPatchVectorField_H

#include "fvPatchFields.H"
#include "fixedValueFvPatchFields.H"
#include "Switch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
               Class SRFVelocityFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class laminarBLInletFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{
    // Private data

        //- Boundary Layer Height at the inlet section
        scalar delta_;

        //- Inlet value [m/s]
        vector flowSpeed_;

        //- Y-Direction. 0 => X axis; 1 => Y axis and 2 => Z axis.
        int yDir_ = 1;


public:

    //- Runtime type information
    TypeName("laminarBLInlet");


    // Constructors

        //- Construct from patch and internal field
        laminarBLInletFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        laminarBLInletFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given SRFVelocityFvPatchVectorField
        //  onto a new patch
        laminarBLInletFvPatchVectorField
        (
            const laminarBLInletFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        laminarBLInletFvPatchVectorField
        (
            const laminarBLInletFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new laminarBLInletFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        laminarBLInletFvPatchVectorField
        (
            const laminarBLInletFvPatchVectorField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchVectorField> clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<fvPatchVectorField>
            (
                new laminarBLInletFvPatchVectorField(*this, iF)
            );
        }


    // Member functions

        // Access

            //- Return const access to the relative flag
            const scalar& delta() const
            {
                return delta_;
            }

        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();

        // I-O

            //- Write
            virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
