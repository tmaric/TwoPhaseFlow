/*---------------------------------------------------------------------------*\
Application
    setPMLFields

Description
    Computes and writes the directional PML damping tensor sigma. Each acoustic
    solver derives its formulation-specific PML coefficients from this field.
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

    volTensorField sigma
    (
        IOobject
        (
            "sigma",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedTensor("sigma", dimless/dimTime, tensor::zero)
    );

    const acousticPML::Config cfg = acousticPML::readConfig(transportProperties);
    acousticPML::updateSigma(mesh, cfg, sigma);
    sigma.write();

    Info<< "Wrote PML damping field sigma at time "
        << runTime.timeName() << endl;

    return 0;
}
