FoamFile
{
    version 2.0;
    format ascii;
    class volScalarField;
    object Pre;
}
dimensions [1 -1 -2 0 0 0 0];
internalField uniform 0;
boundaryField
{
    left { type zeroGradient; }
    right { type zeroGradient; }
    yMin { type empty; }
    yMax { type empty; }
    zMin { type empty; }
    zMax { type empty; }
}
