/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2008 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "theta_bc.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "fvCFD.H"
#include "fvcSnGrad.H"
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::theta_bc::
theta_bc
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(p, iF),

    D_os_(1.0),
    S_TPB_(1.0),
    kox_TPB_(1.0),
    kred_TPB_(1.0),
    tao_(1.0),
    kdes_(1.0),
    kads_(1.0),
    T_(1.0),
    d_phi_(1.0)
{
gradient() = 0.0;
}

Foam::theta_bc::
theta_bc
(
    const theta_bc& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedGradientFvPatchScalarField(ptf, p, iF, mapper),

    D_os_(ptf.D_os_),
    S_TPB_(ptf.S_TPB_),
    kox_TPB_(ptf.kox_TPB_),
    kred_TPB_(ptf.kred_TPB_),
    tao_(ptf.tao_),
    kdes_(ptf.kdes_),
    kads_(ptf.kads_),
    T_(ptf.T_),
    d_phi_(ptf.d_phi_)
{}


Foam::theta_bc::
theta_bc
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedGradientFvPatchScalarField(p, iF, dict)
{
D_os_=readScalar(dict.lookup("D_os"));
S_TPB_=readScalar(dict.lookup("S_TPB"));
kox_TPB_=readScalar(dict.lookup("kox_TPB"));
kred_TPB_=readScalar(dict.lookup("kred_TPB"));
tao_=readScalar(dict.lookup("tao"));
kdes_=readScalar(dict.lookup("kdes"));
kads_=readScalar(dict.lookup("kads"));
T_=readScalar(dict.lookup("T"));
d_phi_=readScalar(dict.lookup("d_phi"));

gradient()=0.0;

   if (dict.found("gradient"))
    {
        fvPatchScalarField::operator=
        (
         scalarField("gradient", dict, p.size())
        );
     }
    else
    {
       fvPatchScalarField::operator=(gradient());
    }

}

Foam::theta_bc::
theta_bc
(
    const theta_bc& tppsf
)
:
    fixedGradientFvPatchScalarField(tppsf),

    D_os_(tppsf.D_os_),
    S_TPB_(tppsf.S_TPB_),
    kox_TPB_(tppsf.kox_TPB_),
    kred_TPB_(tppsf.kred_TPB_),
    tao_(tppsf.tao_),
    kdes_(tppsf.kdes_),
    kads_(tppsf.kads_),
    T_(tppsf.T_),
    d_phi_(tppsf.d_phi_)
{}


Foam::theta_bc::
theta_bc
(
    const theta_bc& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(tppsf, iF),
    D_os_(tppsf.D_os_),
    S_TPB_(tppsf.S_TPB_),
    kox_TPB_(tppsf.kox_TPB_),
    kred_TPB_(tppsf.kred_TPB_),
    tao_(tppsf.tao_),
    kdes_(tppsf.kdes_),
    kads_(tppsf.kads_),
    T_(tppsf.T_),
    d_phi_(tppsf.d_phi_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::theta_bc::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedGradientFvPatchScalarField::autoMap(m);
}


void Foam::theta_bc::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    fixedGradientFvPatchScalarField::rmap(ptf, addr);
}


void Foam::theta_bc::updateCoeffs()
{
    if (updated())
    {
        return;
    }

   const fvPatchScalarField& Cw = *this;

forAll(Cw, counter)
	{
if (Cw[counter]>0  && Cw[counter]<50)
{
	gradient()= -S_TPB_/D_os_*tao_*(kox_TPB_*Foam::exp(96485*d_phi_/8.31/T_)*(Cw/50)*0.975-kred_TPB_*Foam::exp(-96485*d_phi_/8.31/T_)*(1-Cw/50)*0.025);
}
else if (Cw[counter]<0)
{
	gradient()= -S_TPB_/D_os_*tao_*(-kred_TPB_*Foam::exp(-96485*d_phi_/8.31/T_)*0.025);
}
else 
{
	gradient()= -S_TPB_/D_os_*tao_*(kox_TPB_*Foam::exp(96485*d_phi_/8.31/T_)*0.975);
}
	}
    //gradient()= -S_TPB_/D_os_*tao_*(kox_TPB_*Foam::exp(96485*d_phi_/8.31/T_)*(Cw/50)*0.975-kred_TPB_*Foam::exp(-96485*d_phi_/8.31/T_)*(1-Cw/50)*0.025);

    fixedGradientFvPatchScalarField::updateCoeffs(); 
}


void Foam::theta_bc::write
(
    Ostream& os
) const
{
    fixedGradientFvPatchScalarField::write(os);
    os.writeKeyword("D_os") << D_os_ << token::END_STATEMENT << nl;
    os.writeKeyword("S_TPB") << S_TPB_ << token::END_STATEMENT << nl;
    os.writeKeyword("kox_TPB") << kox_TPB_ << token::END_STATEMENT << nl;
    os.writeKeyword("kred_TPB") << kred_TPB_ << token::END_STATEMENT << nl;
    os.writeKeyword("tao") << tao_ << token::END_STATEMENT << nl;
    os.writeKeyword("kdes") << kdes_ << token::END_STATEMENT << nl;
    os.writeKeyword("kads") << kads_ << token::END_STATEMENT << nl;
    os.writeKeyword("T") << T_ << token::END_STATEMENT << nl;
    os.writeKeyword("d_phi") << d_phi_ << token::END_STATEMENT << nl;

    os.writeKeyword("i") << D_os_*gradient()*2*96485 << token::END_STATEMENT << nl;
writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        theta_bc
    );
}

// ************************************************************************* //
