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

Class
    zeroIonFluxFvPatchScalarField

Description
    ion flux boundary condition

    Example usage:

        myWallPatch
        {
            type            ionNormalFlux;
            electricField   electricField;           // Name of electric field
	    Zi              Zi;
        }

SourceFiles
    zeroIonFluxFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef theta_bc_H
#define theta_bc_H

#include "fixedGradientFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
      Class zeroIonFluxFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class theta_bc
:
public fixedGradientFvPatchScalarField 
{
    // Private data

        //- Name of potential field
        scalar D_os_;
	scalar S_TPB_;
        scalar kox_TPB_;
        scalar kred_TPB_;
	scalar tao_;
	scalar kdes_;
	scalar kads_;
	scalar T_;
        scalar d_phi_;

public:

    //- Runtime type information
    TypeName("theta_bc");


    // Constructors

        //- Construct from patch and internal field
        theta_bc
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        theta_bc
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        // solidWallHeatFluxTemperatureFvPatchScalarField
        // onto a new patch
        theta_bc
        (
            const theta_bc&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        theta_bc
        (
            const theta_bc&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new theta_bc(*this)
            );
        }

        //- Construct as copy setting internal field reference
        theta_bc
        (
            const theta_bc&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new theta_bc(*this, iF)
            );
        }


    // Member functions

        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();

//Access

//- Return the electric field name
//const word& 

        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const fvPatchFieldMapper&
            );

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap
            (
                const fvPatchScalarField&,
                const labelList&
            );


        // I-O

            //- Write
            void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
