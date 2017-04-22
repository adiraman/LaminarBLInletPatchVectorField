/*---------------------------------------------------------------------------* \
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "laminarBLInletFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::laminarBLInletFvPatchVectorField::laminarBLInletFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF),
    delta_(0),
    flowSpeed_(0,0,0)
{}


Foam::laminarBLInletFvPatchVectorField::laminarBLInletFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF)
{
    dict.lookup("delta") >> delta_;
    dict.lookup("flowSpeed") >> flowSpeed_;
}


Foam::laminarBLInletFvPatchVectorField::laminarBLInletFvPatchVectorField
(
    const laminarBLInletFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper),
    delta_(ptf.delta_),
    flowSpeed_(ptf.flowSpeed_)
{}




Foam::laminarBLInletFvPatchVectorField::laminarBLInletFvPatchVectorField
(
    const laminarBLInletFvPatchVectorField& srfvpvf
)
:
    fixedValueFvPatchVectorField(srfvpvf),
    delta_(srfvpvf.delta_),
    flowSpeed_(srfvpvf.flowSpeed_)
{}


Foam::laminarBLInletFvPatchVectorField::laminarBLInletFvPatchVectorField
(
    const laminarBLInletFvPatchVectorField& srfvpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(srfvpvf, iF),
    delta_(srfvpvf.delta_),
    flowSpeed_(srfvpvf.flowSpeed_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::laminarBLInletFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // Define a container to hold the data
    vectorField Uin = patch().Cf()*0.;

    // iterate over all faces in the patch
    forAll(patch().Cf(), facei)
    {
        scalar yD = patch().Cf()[facei][yDir_]/delta_;

        if (yD > 1.0)
        {
            Uin[facei] = flowSpeed_;
        }
        else
        {
            Uin[facei] = flowSpeed_*(2*yD - pow(yD, 2.0));
        }
    }

    fvPatchVectorField::operator=(Uin);

    fixedValueFvPatchVectorField::updateCoeffs();
}


void Foam::laminarBLInletFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    os.writeKeyword("delta") << delta_ << token::END_STATEMENT << nl;
    os.writeKeyword("flowSpeed") << flowSpeed_ << token::END_STATEMENT << nl;
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        laminarBLInletFvPatchVectorField
    );
}

// ************************************************************************* //
