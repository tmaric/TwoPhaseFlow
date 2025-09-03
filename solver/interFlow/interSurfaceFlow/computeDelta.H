forAll(delta, cellI)
{
    const scalar rdfCell = RDF[cellI];
    const scalar rdfCellAbs = std::abs(rdfCell);
    if ((rdfCellAbs <= support) && (rdfCellAbs > SMALL)) 
    {
        delta[cellI] = (1.0/(sigma*Foam::sqrt(2.0*M_PI)))
                                * Foam::exp(-0.5*sqr(rdfCell/sigma));
    }
    else
    {
        delta[cellI] = 0.0;
    }
}
