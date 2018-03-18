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

#include "vacancy_bc.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "fvCFD.H"
#include "fvcSnGrad.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::vacancy_bc::
vacancy_bc
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(p, iF),

    D_v_(1.0),
    Sp_CGO_MIEC_(1.0),
    kplus0_(1.0),
    kplus1_(1.0),
    kplus2_(1.0),
    kplus3_(1.0),
    kplus4_(1.0),
    kplus5_(1.0),
    kminus0_(1.0),
    kminus1_(1.0),
    kminus2_(1.0),
    kminus3_(1.0),
    kminus4_(1.0),
    kminus5_(1.0),
    T_(1.0),
    d_phi_(1.0),
    C_omax_(1.0),
    alpha_1_(1.0),
    alpha_2_(1.0),
	var_energy_(1.0),
	var_diff_(1.0)
{
gradient() = 0.0;
}

Foam::vacancy_bc::
vacancy_bc
(
    const vacancy_bc& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedGradientFvPatchScalarField(ptf, p, iF, mapper),

    D_v_(ptf.D_v_),
    Sp_CGO_MIEC_(ptf.Sp_CGO_MIEC_),
    kplus0_(ptf.kplus0_),
    kplus1_(ptf.kplus1_),
    kplus2_(ptf.kplus2_),
    kplus3_(ptf.kplus3_),
    kplus4_(ptf.kplus4_),
    kplus5_(ptf.kplus5_),
    kminus0_(ptf.kminus0_),
    kminus1_(ptf.kminus1_),
    kminus2_(ptf.kminus2_),
    kminus3_(ptf.kminus3_),
    kminus4_(ptf.kminus4_),
    kminus5_(ptf.kminus5_),
    T_(ptf.T_),
    d_phi_(ptf.d_phi_),
    C_omax_(ptf.C_omax_),
    alpha_1_(ptf.alpha_1_),
    alpha_2_(ptf.alpha_2_),
    var_energy_(ptf.var_energy_),
    var_diff_(ptf.var_diff_)
{}


Foam::vacancy_bc::
vacancy_bc
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedGradientFvPatchScalarField(p, iF, dict)
{
D_v_=readScalar(dict.lookup("D_v"));
Sp_CGO_MIEC_=readScalar(dict.lookup("Sp_CGO_MIEC"));
kplus0_=readScalar(dict.lookup("kplus0"));
kplus1_=readScalar(dict.lookup("kplus1"));
kplus2_=readScalar(dict.lookup("kplus2"));
kplus3_=readScalar(dict.lookup("kplus3"));
kplus4_=readScalar(dict.lookup("kplus4"));
kplus5_=readScalar(dict.lookup("kplus5"));
kminus0_=readScalar(dict.lookup("kminus0"));
kminus1_=readScalar(dict.lookup("kminus1"));
kminus2_=readScalar(dict.lookup("kminus2"));
kminus3_=readScalar(dict.lookup("kminus3"));
kminus4_=readScalar(dict.lookup("kminus4"));
kminus5_=readScalar(dict.lookup("kminus5"));
T_=readScalar(dict.lookup("T"));
d_phi_=readScalar(dict.lookup("d_phi"));
C_omax_=readScalar(dict.lookup("C_omax"));
alpha_1_=readScalar(dict.lookup("alpha_1"));
alpha_2_=readScalar(dict.lookup("alpha_2"));
var_energy_=readScalar(dict.lookup("var_energy"));
var_diff_=readScalar(dict.lookup("var_diff"));

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

Foam::vacancy_bc::
vacancy_bc
(
    const vacancy_bc& tppsf
)
:
    fixedGradientFvPatchScalarField(tppsf),
    D_v_(tppsf.D_v_),
    Sp_CGO_MIEC_(tppsf.Sp_CGO_MIEC_),
    kplus0_(tppsf.kplus0_),
    kplus1_(tppsf.kplus1_),
    kplus2_(tppsf.kplus2_),
    kplus3_(tppsf.kplus3_),
    kplus4_(tppsf.kplus4_),
    kplus5_(tppsf.kplus5_),
    kminus0_(tppsf.kminus0_),
    kminus1_(tppsf.kminus1_),
    kminus2_(tppsf.kminus2_),
    kminus3_(tppsf.kminus3_),
    kminus4_(tppsf.kminus4_),
    kminus5_(tppsf.kminus5_),
    T_(tppsf.T_),
    d_phi_(tppsf.d_phi_),
    C_omax_(tppsf.C_omax_),
    alpha_1_(tppsf.alpha_1_),
    alpha_2_(tppsf.alpha_2_),
    var_energy_(tppsf.var_energy_),
    var_diff_(tppsf.var_diff_)
{}


Foam::vacancy_bc::
vacancy_bc
(
    const vacancy_bc& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(tppsf, iF),
    D_v_(tppsf.D_v_),
    Sp_CGO_MIEC_(tppsf.Sp_CGO_MIEC_),
    kplus0_(tppsf.kplus0_),
    kplus1_(tppsf.kplus1_),
    kplus2_(tppsf.kplus2_),
    kplus3_(tppsf.kplus3_),
    kplus4_(tppsf.kplus4_),
    kplus5_(tppsf.kplus5_),
    kminus0_(tppsf.kminus0_),
    kminus1_(tppsf.kminus1_),
    kminus2_(tppsf.kminus2_),
    kminus3_(tppsf.kminus3_),
    kminus4_(tppsf.kminus4_),
    kminus5_(tppsf.kminus5_),
    T_(tppsf.T_),
    d_phi_(tppsf.d_phi_),
    C_omax_(tppsf.C_omax_),
    alpha_1_(tppsf.alpha_1_),
    alpha_2_(tppsf.alpha_2_),
    var_energy_(tppsf.var_energy_),
    var_diff_(tppsf.var_diff_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::vacancy_bc::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedGradientFvPatchScalarField::autoMap(m);
}


void Foam::vacancy_bc::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    fixedGradientFvPatchScalarField::rmap(ptf, addr);
}


void Foam::vacancy_bc::updateCoeffs()
{
    if (updated())
    {
        return;
    }

	const fvPatchScalarField& Cw = *this;
	int ite = 0;                          //iteration counter

	scalar v0_;
	scalar v1_;
	scalar v2_;
	scalar v3_;
	scalar v4_;
	scalar v5_=2078.675;
	scalar kplus0_ref = kplus0_;
	scalar kplus1_ref = kplus1_;
	scalar kplus2_ref = kplus2_;
	scalar kplus3_ref = kplus3_;
	scalar kplus4_ref = kplus4_;
	scalar kplus5_ref = kplus5_;
	scalar kminus0_ref = kminus0_;
	scalar kminus1_ref = kminus1_;
	scalar kminus2_ref = kminus2_;
	scalar kminus3_ref = kminus3_;
	scalar kminus4_ref = kminus4_;
	scalar kminus5_ref = kminus5_;
	scalar capacitance=20;
	scalar over_p;

	forAll(Cw, counter)
	{
	kplus0_ref = kplus0_*Foam::exp(var_energy_*Cw[counter]/C_omax_)*Foam::exp(-var_energy_*(1154-554)/C_omax_);
	kminus0_ref = kminus0_;
	kplus1_ref = kplus1_;
	kminus1_ref = kminus1_;
	kplus2_ref = kplus2_;//*Foam::exp(-var_diff_*Cw[counter]/C_omax_);///**Foam::exp(52.8*1154.82/C_omax_)*/*Foam::exp(96485*d_phi_*alpha_1_/8.31/T_);
	kminus2_ref = kminus2_;//*Foam::exp(-var_diff_*Cw[counter]/C_omax_);///**Foam::exp(52.8*1154.82/C_omax_)*/*Foam::exp(-96485*d_phi_*alpha_2_/8.31/T_);
	kplus3_ref = kplus3_;
	kminus3_ref = kminus3_*Foam::exp(-96485*d_phi_*alpha_1_/8.31/T_);//*Foam::exp(-var_energy_*Cw[counter]/C_omax_)*Foam::exp(var_energy_*600/C_omax_);
	kplus4_ref = kplus4_*Foam::exp(96485*d_phi_*alpha_2_/8.31/T_);
	kminus4_ref = kminus4_;
	kplus5_ref = kplus5_;
	kminus5_ref = kminus5_;
	}

	v4_ = C_omax_*kminus5_ref*v5_/(kplus5_ref*(C_omax_-v5_)+kminus5_ref*v5_);    // initial values for v0, v1, v2, v3, v4. Reaction rates are 0
	v3_ = C_omax_*kminus4_ref*v4_/(kplus4_ref*(C_omax_-v4_)+kminus4_ref*v4_);    
	v2_ = C_omax_*kminus3_ref*v3_/(kplus3_ref*(C_omax_-v3_)+kminus3_ref*v3_);
	v1_ = C_omax_*kminus2_ref*v2_/(kplus2_ref*(C_omax_-v2_)+kminus2_ref*v2_);
	v0_ = C_omax_*kminus1_ref*v1_/(kplus1_ref*(C_omax_-v1_)+kminus1_ref*v1_);

	for (ite = 0;ite < 5;ite ++)
	{
		forAll(Cw, counter)
		{
			v0_ = C_omax_*(kplus0_ref*Cw[counter]+kminus1_ref*v1_)/(kplus0_ref*Cw[counter]+kminus0_ref*(C_omax_-Cw[counter])+kplus1_ref*(C_omax_-v1_)+kminus1_ref*v1_);
			v1_ = C_omax_*(kplus1_ref*v0_+kminus2_ref*v2_)/(kplus1_ref*v0_+kminus1_ref*(C_omax_-v0_)+kplus2_ref*(C_omax_-v2_)+kminus2_ref*v2_);
			v2_ = C_omax_*(kplus2_ref*v1_+kminus3_ref*v3_)/(kplus2_ref*v1_+kminus2_ref*(C_omax_-v1_)+kplus3_ref*(C_omax_-v3_)+kminus3_ref*v3_);
			v3_ = C_omax_*(kplus3_ref*v2_+kminus4_ref*v4_)/(kplus3_ref*v2_+kminus3_ref*(C_omax_-v2_)+kplus4_ref*(C_omax_-v4_)+kminus4_ref*v4_);
			v4_ = C_omax_*(kplus4_ref*v3_+kminus5_ref*v5_)/(kplus4_ref*v3_+kminus4_ref*(C_omax_-v3_)+kplus5_ref*(C_omax_-v5_)+kminus5_ref*v5_);
		}
	}
//***************************************************************************************************************//	
	gradient()= -Sp_CGO_MIEC_*Foam::exp(-var_diff_*Cw/C_omax_)*Foam::exp(var_diff_*600/C_omax_)/(D_v_*Foam::exp(-var_diff_*Cw/C_omax_)*Foam::exp(var_diff_*600/C_omax_))*1e-5*(kplus0_ref*(Cw/C_omax_)*(1-v0_/C_omax_)-kminus0_ref*(1-Cw/C_omax_)*v0_/C_omax_);
	//gradient()= -Sp_CGO_MIEC_*Foam::exp(-var_diff_*Cw/C_omax_)/(D_v_*Foam::exp(-var_diff_*Cw/C_omax_)*Foam::exp(var_diff_*600/C_omax_))*1e-5*(kplus1_ref*(v0_/C_omax_)*(1-v1_/C_omax_)-kminus1_ref*(1-v0_/C_omax_)*v1_/C_omax_);

//D_v_ = D_v_*Foam::exp(-149.72*Cw[0]/C_omax_);
//kplus_ = kplus_*Foam::exp(-149.72*Cw/C_omax_)
//kminus_ = kminus_*Foam::exp(-149.72*Cw/C_omax_)

	fixedGradientFvPatchScalarField::updateCoeffs(); 
}


void Foam::vacancy_bc::write
(
    Ostream& os
) const
{
	const fvPatchScalarField& Cw = *this;
    fixedGradientFvPatchScalarField::write(os);
    os.writeKeyword("D_v") << D_v_*Foam::exp(-var_diff_*Cw/C_omax_)*Foam::exp(var_diff_*600/C_omax_) << token::END_STATEMENT << nl;
    os.writeKeyword("Sp_CGO_MIEC") << Sp_CGO_MIEC_ << token::END_STATEMENT << nl;
    os.writeKeyword("kplus0") << kplus0_ << token::END_STATEMENT << nl;
    os.writeKeyword("kplus1") << kplus1_ << token::END_STATEMENT << nl;
    os.writeKeyword("kplus2") << kplus2_ << token::END_STATEMENT << nl;
    os.writeKeyword("kplus3") << kplus3_ << token::END_STATEMENT << nl;
    os.writeKeyword("kplus4") << kplus4_ << token::END_STATEMENT << nl;
    os.writeKeyword("kplus5") << kplus5_ << token::END_STATEMENT << nl;
    os.writeKeyword("kminus0") << kminus0_ << token::END_STATEMENT << nl;
    os.writeKeyword("kminus1") << kminus1_ << token::END_STATEMENT << nl;
    os.writeKeyword("kminus2") << kminus2_ << token::END_STATEMENT << nl;
    os.writeKeyword("kminus3") << kminus3_ << token::END_STATEMENT << nl;
    os.writeKeyword("kminus4") << kminus4_ << token::END_STATEMENT << nl;
    os.writeKeyword("kminus5") << kminus5_ << token::END_STATEMENT << nl;
    os.writeKeyword("T") << T_ << token::END_STATEMENT << nl;
    os.writeKeyword("d_phi") << d_phi_ << token::END_STATEMENT << nl;
    os.writeKeyword("C_omax") << C_omax_ << token::END_STATEMENT << nl;
    os.writeKeyword("alpha_1") << alpha_1_ << token::END_STATEMENT << nl;
    os.writeKeyword("alpha_2") << alpha_2_ << token::END_STATEMENT << nl;
    os.writeKeyword("var_energy") << var_energy_ << token::END_STATEMENT << nl;
    os.writeKeyword("var_diff") << var_diff_ << token::END_STATEMENT << nl;
    os.writeKeyword("i") << D_v_*Foam::exp(-var_diff_*Cw/C_omax_)*Foam::exp(var_diff_*600/C_omax_)*gradient()*2*96485 << token::END_STATEMENT << nl;
writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        vacancy_bc
    );
}

// ************************************************************************* //
