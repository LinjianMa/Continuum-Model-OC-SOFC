/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
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

Application
    Sr_o2p

Description
    Transient solver for incompressible, laminar flow of Newtonian fluids.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "fixedGradientFvPatchFields.H"
#include "vacancy_bc.H"
#include "theta_bc.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createFields.H"
    #include "readSolverControls.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;
	double theta_vbulk,theta_obulk,theta_osuf,theta_vsuf,theta_o2suf;
	double m,n,q;                                //the expression of m,n,q can be found in the notes
	double kads_ref = kads;
	double kdes_ref = kdes;
	double kox_ref = kox_MIEC_Gas;               // store the referance value of diffusivity
	double kred_ref = kred_MIEC_Gas;             // store the referance value of diffusivity
	double kox_MIEC_Gas_2_3ref = kox_MIEC_Gas_2_3;
	int    ite = 0;                              // record how many iterations have finished
	int    ite_eqi = 0;                          // used to do the iteration for the equilibrium
	double dphi_const;                           // used in phi calculation
	double v_3, v_4, o_3, o_4;                   // used in surface/subsurface reaction 
	double over_surf;
	// y_o2 = 11.327
//***************************************************************************************************************//	
if (runTime.timeName() == "0")
{
		forAll(mesh.C(),counter)
	{
		C_ads.internalField()[counter] = kads_phys/kdes_phys*y_o2.internalField()[counter];
		n = kads*C_ads.internalField()[counter]*kred_MIEC_Gas*kred_MIEC_Gas_2*kred_MIEC_Gas_2/(kdes*kox_MIEC_Gas*kox_MIEC_Gas_2*kox_MIEC_Gas_2);
		//11.327 is the concentration of O2
        n = Foam::sqrt(n);  //the expression of n can be found in the notes
		C_v.internalField()[counter] = 1/(n+1)*C_omax;
		C_o.internalField()[counter] = C_omax - C_v.internalField()[counter];
		m = kox_MIEC_Gas*(kads*C_ads.internalField()[counter]+kdes)*(kads*C_ads.internalField()[counter]+kdes)/(kdes*kred_MIEC_Gas*kads*C_ads.internalField()[counter]);
		m = Foam::sqrt(m);  //the expression of m can be found in the notes
		C_osuf.internalField()[counter] = 1/(m+1)*C_omax_suf;
		q = kads*C_ads.internalField()[counter]/(kads*C_ads.internalField()[counter]+kdes);  //the expression of q can be found in the notes
		C_o2suf.internalField()[counter] = q*(C_omax_suf-C_osuf.internalField()[counter]);
		C_vsuf.internalField()[counter] = C_omax_suf-C_o2suf.internalField()[counter]-C_osuf.internalField()[counter];

		//initialize the reaction rates
		V_3.internalField()[counter] = 0;
		V_2.internalField()[counter] = 0;
		V_1.internalField()[counter] = 0;
	}
}
//D_o2.value() = 10000*2/3*r_pore*porosity*tortuosity*Foam::sqrt(8*R*T/Pi/M_o2)*pressure/R/T;
//*****************************************************************************************************************//
    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

	forAll(mesh.C(),counter)
	{
	//Iterate on k, details can be found in the notes
	//over_surf = (C_osuf.internalField()[counter]-1/(m+1)*C_omax_suf)/5e6*96485*2/over_ratio;
	over_surf = (C_osuf.internalField()[counter]-33.404)/5e6*96485*2/over_ratio;
    //kads = kads_ref*Foam::exp(-96485*alpha_surf*over_surf/8.31/1073);
    //kdes = kdes_ref*Foam::exp(96485*(2-alpha_surf)*over_surf/8.31/1073);
	//kads = kads_ref * Foam::exp(var_energy*C_v.internalField()[counter]/C_omax) * Foam::exp(-var_energy*600/C_omax);
	kox_MIEC_Gas = kox_ref*Foam::exp(-var_k*C_v.internalField()[counter]/C_omax)*Foam::exp(96485*(2-alpha_surf)*over_surf/8.31/1073);//*Foam::exp(52.8*1154.82/C_omax);
	kred_MIEC_Gas = kred_ref*Foam::exp(-var_k*C_v.internalField()[counter]/C_omax)*Foam::exp(-96485*alpha_surf*over_surf/8.31/1073);//*Foam::exp(52.8*1154.82/C_omax);
	kox_MIEC_Gas_2_3 = kox_MIEC_Gas_2_3ref * Foam::exp(-var_energy*C_v.internalField()[counter]/C_omax) * Foam::exp(var_energy*(1154-554)/C_omax);
    /////////////////////////////
	C_o.internalField()[counter] = C_omax - C_v.internalField()[counter];
	C_osuf.internalField()[counter] = C_omax_suf - C_vsuf.internalField()[counter] - C_o2suf.internalField()[counter];
	theta_vbulk = C_v.internalField()[counter]/C_omax;
    theta_obulk = C_o.internalField()[counter]/C_omax;
	theta_osuf = C_osuf.internalField()[counter]/C_omax_suf;
	theta_vsuf = C_vsuf.internalField()[counter]/C_omax_suf;
	theta_o2suf = C_o2suf.internalField()[counter]/C_omax_suf;
	//Iterate on surface/subsurface reaction 
	v_3 = kox_MIEC_Gas_2_1*C_omax_suf*C_vsuf.internalField()[counter]/(kred_MIEC_Gas_2_1*C_osuf.internalField()[counter] + kox_MIEC_Gas_2_1*C_vsuf.internalField()[counter]);
	o_3 = C_omax_suf - v_3;
	v_4 = kred_MIEC_Gas_2_3*C_omax_suf*C_v.internalField()[counter]/(kred_MIEC_Gas_2_3*C_v.internalField()[counter] + kox_MIEC_Gas_2_3*C_o.internalField()[counter]);
	o_4 = C_omax_suf - v_4;
	for (ite_eqi=0;ite_eqi<4;ite_eqi++)
	{
		v_3 = (kred_MIEC_Gas_2_2*v_4 + kox_MIEC_Gas_2_1*C_vsuf.internalField()[counter])*C_omax_suf/(kred_MIEC_Gas_2_1*C_osuf.internalField()[counter] + kox_MIEC_Gas_2_1*C_vsuf.internalField()[counter] + kred_MIEC_Gas_2_2*v_4 + kox_MIEC_Gas_2_2*o_4);
		v_4 = (kred_MIEC_Gas_2_3*C_v.internalField()[counter]/C_omax + kox_MIEC_Gas_2_2*v_3/C_omax_suf)*C_omax_suf/(kred_MIEC_Gas_2_2*o_3/C_omax_suf + kox_MIEC_Gas_2_2*v_3/C_omax_suf + kred_MIEC_Gas_2_3*C_v.internalField()[counter]/C_omax + kox_MIEC_Gas_2_3*C_o.internalField()[counter]/C_omax);
		o_3 = C_omax_suf - v_3;
		o_4 = C_omax_suf - v_4;
	}
    /////////////////////////////
	V_3.internalField()[counter] = 1e-6*per_surf*(kox_MIEC_Gas_2_3*theta_obulk*v_4-kred_MIEC_Gas_2_3*theta_vbulk*o_4);
	V_2.internalField()[counter] = 1e-6*per_surf*(kox_MIEC_Gas*C_osuf.internalField()[counter]*theta_osuf - kred_MIEC_Gas*C_o2suf.internalField()[counter]*theta_vsuf);
	V_1.internalField()[counter] = 1e-6*per_surf*(kdes*C_o2suf.internalField()[counter] - kads*C_ads.internalField()[counter]*C_vsuf.internalField()[counter]);
	V_0.internalField()[counter] = 1e-6*per_surf*(kdes_phys*C_ads.internalField()[counter] - kads_phys*y_o2.internalField()[counter]);
	}

	// Iterations on the phi 
	//solve
	//(
	//	fvm::ddt(grad_phi_o)
	//);
	dphi = -2*F*D_v/(Sigma+4*F*F*C_o*D_v/(R*T))*fvc::grad(C_v);
	dphi_const = dphi.internalField()[2].y();
	forAll(mesh.C(),counter)
	{
		if (ite <= 2000)
		{
			dphi.internalField()[counter].y() = 0;
		}
		else
		{
			dphi.internalField()[counter].y() = dphi.internalField()[counter].y() - dphi_const;
		}
		grad_phi_o.internalField()[counter] = 2*F/(R*T)*D_v.internalField()[counter]*dphi.internalField()[counter];
		grad_phi_vsuf.internalField()[counter] = 2*F/(R*T)*D_vsuf.value()*dphi.internalField()[counter];
		grad_phi_o2suf.internalField()[counter] = 2*F/(R*T)*D_o2suf.value()*dphi.internalField()[counter];
	}
	phi_flux_o = linearInterpolate(grad_phi_o) & mesh.Sf();
	phi_flux_vsuf = linearInterpolate(grad_phi_vsuf) & mesh.Sf();
	phi_flux_o2suf = linearInterpolate(grad_phi_o2suf) & mesh.Sf();
	// end for phi

	//if (ite <= 2000000000000)
	//{
	C_v.storePrevIter();
	solve
		(
			fvm::ddt(C_v)
			- fvm::laplacian(D_v,C_v) - fvm::div(phi_flux_o, C_v)
			==
			V_3-C_omax*fvc::div(phi_flux_o)
		);
	C_v.relax();

/*for (ite_eqi=0;ite_eqi<4;ite_eqi++)
		{
			forAll(mesh.C(),counter)
			{
				C_o.internalField()[counter] = C_omax - C_v.internalField()[counter];
				theta_vbulk = C_v.internalField()[counter]/C_omax;
				theta_obulk = C_o.internalField()[counter]/C_omax;
				C_osuf.internalField()[counter] = C_omax_suf - C_vsuf.internalField()[counter] - C_o2suf.internalField()[counter];
				theta_osuf = C_osuf.internalField()[counter]/C_omax_suf;
				theta_vsuf = C_vsuf.internalField()[counter]/C_omax_suf;
				theta_o2suf = C_o2suf.internalField()[counter]/C_omax_suf;
				//C_vsuf.internalField()[counter]   =  (-2*kox_MIEC_Gas*C_osuf.internalField()[counter]*theta_osuf+kox_MIEC_Gas_2_3*theta_obulk*v_4-kred_MIEC_Gas_2_3*theta_vbulk*o_4)/(- 2*kred_MIEC_Gas*C_o2suf.internalField()[counter]/C_omax_suf);
				C_o2suf.internalField()[counter]  =   (kads*C_ads.internalField()[counter]*C_vsuf.internalField()[counter]+kox_MIEC_Gas*C_osuf.internalField()[counter]*theta_osuf)/(kdes+kred_MIEC_Gas*theta_vsuf);
			}
		}
for (ite_eqi=0;ite_eqi<4;ite_eqi++)
		{
			forAll(mesh.C(),counter)
			{
				C_o.internalField()[counter] = C_omax - C_v.internalField()[counter];
				theta_vbulk = C_v.internalField()[counter]/C_omax;
				theta_obulk = C_o.internalField()[counter]/C_omax;
				C_osuf.internalField()[counter] = C_omax_suf - C_vsuf.internalField()[counter] - C_o2suf.internalField()[counter];
				theta_osuf = C_osuf.internalField()[counter]/C_omax_suf;
				theta_vsuf = C_vsuf.internalField()[counter]/C_omax_suf;
				theta_o2suf = C_o2suf.internalField()[counter]/C_omax_suf;
				//C_vsuf.internalField()[counter]   =  (-kox_MIEC_Gas*C_osuf.internalField()[counter]*theta_osuf-kdes*C_o2suf.internalField()[counter]+kox_MIEC_Gas_2_3*theta_obulk*v_4-kred_MIEC_Gas_2_3*theta_vbulk*o_4)/(-kads*C_ads.internalField()[counter]- kred_MIEC_Gas*C_o2suf.internalField()[counter]/C_omax_suf);
				//C_o2suf.internalField()[counter]  =   (kads*C_ads.internalField()[counter]*C_vsuf.internalField()[counter]+kox_MIEC_Gas*C_osuf.internalField()[counter]*theta_osuf)/(kdes+kred_MIEC_Gas*theta_vsuf);
				C_vsuf.internalField()[counter]   =  (-2*kox_MIEC_Gas*C_osuf.internalField()[counter]*theta_osuf+kox_MIEC_Gas_2_3*theta_obulk*v_4-kred_MIEC_Gas_2_3*theta_vbulk*o_4)/(- 2*kred_MIEC_Gas*C_o2suf.internalField()[counter]/C_omax_suf);
			}
		}*/
	C_vsuf.storePrevIter();
	solve
		(
			time_vsuf*fvm::ddt(C_vsuf)
			- fvm::laplacian(D_vsuf, C_vsuf) - fvm::div(phi_flux_vsuf, C_vsuf)
			==
			V_1+V_2-V_3-C_omax_suf*fvc::div(grad_phi_vsuf)
		);
	C_vsuf.relax();

	C_o2suf.storePrevIter();
	solve
		(
			time_o2suf*fvm::ddt(C_o2suf)
			- fvm::laplacian(D_o2suf, C_o2suf) + fvm::div(phi_flux_o2suf, C_o2suf)
			==
			-V_1+V_2
		);
	C_o2suf.relax();

	C_ads.storePrevIter();
	solve
		(
			1e9*fvm::ddt(C_ads)
			- fvm::laplacian(D_ads, C_ads)
			==
			V_1-V_0
		);
	C_ads.relax();

	y_o2.storePrevIter();
	solve
		(
			1e9*fvm::ddt(y_o2)
			- fvm::laplacian(D_o2, y_o2)
			==
			V_0
		);
	y_o2.relax();

	//}
	/*else 
	{
		// C_v iteration
		C_v.storePrevIter();
		solve
			(
				fvm::ddt(C_v)
				- fvm::laplacian(D_v,C_v) - fvm::div(phi_flux_o, C_v)
				==
				V_3-C_omax*fvc::div(phi_flux_o)
			);
		C_v.relax();
		// end for C_v
		for (ite_eqi=0;ite_eqi<4;ite_eqi++)
		{
			forAll(mesh.C(),counter)
			{
				C_o.internalField()[counter] = C_omax - C_v.internalField()[counter];

				kox_MIEC_Gas = kox_ref*Foam::exp(-52.8*C_v.internalField()[counter]/C_omax)*Foam::exp(-120*C_v.internalField()[counter]/C_omax)*Foam::exp(120*600/C_omax);
				kred_MIEC_Gas = kred_ref*Foam::exp(-52.8*C_v.internalField()[counter]/C_omax);
				A = kred_MIEC_Gas/C_omax_suf;
				B = -(kred_MIEC_Gas/C_omax_suf*(C_omax_suf-C_osuf.internalField()[counter])+0.5*kox_MIEC_Gas_2/C_omax*C_o.internalField()[counter]);
				C = kox_MIEC_Gas/C_omax_suf*C_osuf.internalField()[counter]*C_osuf.internalField()[counter]+0.5*kred_MIEC_Gas_2*C_v.internalField()[counter]/C_omax*C_osuf.internalField()[counter];
				C_vsuf.internalField()[counter] = (-B-Foam::sqrt(B*B-4*A*C))/(2*A);
				D = 0.5*kox_MIEC_Gas_2*C_o.internalField()[counter]/C_omax+kdes+11.327*kads;
				E = -kdes*C_omax_suf;
				f = 0.5*kred_MIEC_Gas_2*C_v.internalField()[counter]/C_omax-kdes;
				C_osuf.internalField()[counter] = (D*C_vsuf.internalField()[counter]+E)/f;
				C_o2suf.internalField()[counter] = C_omax_suf - C_vsuf.internalField()[counter] - C_osuf.internalField()[counter];
			}
		}
	}*/
	//D_v iteration
	solve
		(
			fvm::ddt(D_v)
		);
	forAll(mesh.C(),counter)
	{
		//Iterate on diffusivity, details can be found in the notes
		D_v.internalField()[counter] = D_v_ref*Foam::exp(-var_diff*C_v.internalField()[counter]/C_omax)*Foam::exp(var_diff*600/C_omax);
	}
	/*forAll(mesh.C(),counter)
	{
		if (C_vsuf.internalField()[counter] > 50)
		{
			C_vsuf.internalField()[counter] = 50;
		}
		else if (C_vsuf.internalField()[counter] < 0)
		{
			C_vsuf.internalField()[counter] = 0;
		}
	}*/
	// end for D_v

	Flux_v = -2*F*(D_v*fvc::grad(C_v)+grad_phi_o*C_o);
	Flux_o2suf = -2*F*(-D_o2suf*fvc::grad(C_o2suf)+grad_phi_o2suf*C_o2suf);
	Flux_vsuf = -2*F*(D_vsuf*fvc::grad(C_vsuf)+grad_phi_vsuf*C_vsuf);
	Flux_o2 = -D_o2*fvc::grad(y_o2);
	//electricfield = -Sigma_h*fvc::grad(Phi);
	
		ite = ite + 1;
        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End" << endl;

    return 0;
}


// ************************************************************************* //
