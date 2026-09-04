if (radiationWindowCount == radiationWindowSamples)
{
    const volScalarField& oldPr = prWindow[radiationWindowIndex];
    const volTensorField& oldMomFlux =
        momFluxWindow[radiationWindowIndex];
    prWindowSum -= oldPr;
    momFluxWindowSum -= oldMomFlux;
    prWindow[radiationWindowIndex] == prT;
    momFluxWindow[radiationWindowIndex] == momFluxT;
}
else
{
    prWindow.set
    (
        radiationWindowIndex,
        new volScalarField
        (
            IOobject
            (
                "prWindow_" + Foam::name(radiationWindowIndex),
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh,
            dimensionedScalar("prWindow", dimPressure, Zero),
            calculatedFvPatchScalarField::typeName
        )
    );

    momFluxWindow.set
    (
        radiationWindowIndex,
        new volTensorField
        (
            IOobject
            (
                "momFluxWindow_" + Foam::name(radiationWindowIndex),
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh,
            dimensionedTensor
            (
                "momFluxWindow",
                dimDensity*dimVelocity*dimVelocity,
                tensor::zero
            ),
            calculatedFvPatchTensorField::typeName
        )
    );

    prWindow[radiationWindowIndex] == prT;
    momFluxWindow[radiationWindowIndex] == momFluxT;

    ++radiationWindowCount;
}

prWindowSum += prT;
momFluxWindowSum += momFluxT;

radiationWindowIndex =
    (radiationWindowIndex + 1) % radiationWindowSamples;

const scalar invWindowCount = 1.0/radiationWindowCount;
pr == invWindowCount*prWindowSum;
momFlux == invWindowCount*momFluxWindowSum;
