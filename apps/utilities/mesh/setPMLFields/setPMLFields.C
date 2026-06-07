/*---------------------------------------------------------------------------*\
Application
    setPMLFields

Description
    Computes and writes the scalar and tensor PML coefficient fields used by the
    frequency-domain acoustic solvers.
\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "PMLFieldSetup.H"

int main(int argc, char *argv[])
{
    #include "addRegionOption.H"
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createNamedMesh.H"

    Info<< "Reading transportProperties\n" << endl;

    IOdictionary transportProperties
    (
        IOobject
        (
            "transportProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    Info<< "Reading alpha.water\n" << endl;

    volScalarField alpha
    (
        IOobject
        (
            "alpha.water",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("alpha.water", dimless, Zero)
    );

    volScalarField SC0
    (
        IOobject
        (
            "SC0",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("SC0", dimensionSet(0, -2, 0, 0, 0, 0, 0), Zero)
    );

    volScalarField SC1
    (
        IOobject
        (
            "SC1",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("SC1", dimensionSet(0, -2, 0, 0, 0, 0, 0), Zero)
    );

    volScalarField TC1
    (
        IOobject
        (
            "TC1",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("TC1", dimless, Zero)
    );

    volTensorField T0
    (
        IOobject
        (
            "T0",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedTensor("T0", dimless, tensor::zero)
    );

    volTensorField T1
    (
        IOobject
        (
            "T1",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedTensor("T1", dimless, tensor::zero)
    );

    volScalarField sigma
    (
        IOobject
        (
            "sigma",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("sigma", dimless/dimTime, Zero)
    );

    const acousticPML::Config cfg = acousticPML::readConfig(transportProperties);
    acousticPML::updateFields(mesh, cfg, alpha, SC0, SC1, TC1, T0, T1, sigma);

    SC0.write();
    SC1.write();
    TC1.write();
    T0.write();
    T1.write();
    sigma.write();

    Info<< "Wrote PML coefficient fields at time " << runTime.timeName() << nl
        << "    SC0, SC1, TC1, T0, T1, sigma" << endl;

    return 0;
}
