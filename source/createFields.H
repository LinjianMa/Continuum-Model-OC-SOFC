Info<< "Reading transportProperties\n" << endl;

IOdictionary transportProperties
(
    IOobject
    (
        "transportProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

Info<< "Reading diffusivity D\n" << endl;

dimensionedScalar D_o2
(
    transportProperties.lookup("D_o2")
);

dimensionedScalar D_ads
(
    transportProperties.lookup("D_ads")
);

dimensionedScalar D_vsuf
(
    transportProperties.lookup("D_vsuf")
);

dimensionedScalar D_o2suf
(
    transportProperties.lookup("D_o2suf")
);

dimensionedScalar Sigma
(
    transportProperties.lookup("Sigma")
);

dimensionedScalar tao
(
    transportProperties.lookup("tao")
);
/***************************************************************/

Info<< "Reading field C\n" << endl;

volScalarField C_v
(
    IOobject
    (
        "C_v",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField C_ads
(
    IOobject
    (
        "C_ads",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField D_v
(
    IOobject
    (
        "D_v",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField y_o2
(
    IOobject
    (
        "y_o2",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField C_o
(
    IOobject
    (
        "C_o",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField C_o2suf
(
    IOobject
    (
        "C_o2suf",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField C_vsuf
(
    IOobject
    (
        "C_vsuf",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField C_osuf
(
    IOobject
    (
        "C_osuf",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading source term V\n" << endl;
volScalarField V_0
(
    IOobject
    (
        "V_0",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField V_1
(
    IOobject
    (
        "V_1",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField V_2
(
    IOobject
    (
        "V_2",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField V_3
(
    IOobject
    (
        "V_3",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading/calculating flux \n" << endl;
    volVectorField Flux_v
    (
        IOobject
        (
            "Flux_v",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
       -D_v*fvc::grad(C_v)
    );

    volVectorField Flux_o2
    (
        IOobject
        (
            "Flux_o2",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
       -D_o2*fvc::grad(y_o2)
    );

    volVectorField Flux_o2suf
    (
        IOobject
        (
            "Flux_o2suf",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
       -D_o2suf*fvc::grad(C_o2suf)
    );

    volVectorField Flux_vsuf
    (
        IOobject
        (
            "Flux_vsuf",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
       -D_vsuf*fvc::grad(C_vsuf)
    );

    volVectorField dphi
    (
        IOobject
        (
            "dphi",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
       -fvc::grad(C_v)
    );

    volVectorField grad_phi_o
    (
        IOobject
        (
            "grad_phi_o",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
       -D_v*fvc::grad(C_v)
    );

	surfaceScalarField phi_flux_o
	(
		IOobject
		(
			"phi_flux_o",
			runTime.timeName(),
			mesh,
			IOobject::READ_IF_PRESENT,
			IOobject::AUTO_WRITE
		),
	linearInterpolate(grad_phi_o) & mesh.Sf()
	);

	    volVectorField grad_phi_o2suf
    (
        IOobject
        (
            "grad_phi_o2suf",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
       -D_v*fvc::grad(C_v)
    );

	surfaceScalarField phi_flux_o2suf
	(
		IOobject
		(
			"phi_flux_o2suf",
			runTime.timeName(),
			mesh,
			IOobject::READ_IF_PRESENT,
			IOobject::AUTO_WRITE
		),
	linearInterpolate(grad_phi_o2suf) & mesh.Sf()
	);

	    volVectorField grad_phi_vsuf
    (
        IOobject
        (
            "grad_phi_vsuf",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
       -D_v*fvc::grad(C_vsuf)
    );

	surfaceScalarField phi_flux_vsuf
	(
		IOobject
		(
			"phi_flux_vsuf",
			runTime.timeName(),
			mesh,
			IOobject::READ_IF_PRESENT,
			IOobject::AUTO_WRITE
		),
	linearInterpolate(grad_phi_vsuf) & mesh.Sf()
	);