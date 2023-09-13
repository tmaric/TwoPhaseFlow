{
// Update alpha1
#include "alphaSuSp.H"
advector->advect(Sp,(Su + divU*min(alpha1(), scalar(1)))());
}

Info<< "Phase-1 volume fraction = "
    << alpha1.weightedAverage(mesh.Vsc()).value()
    << "  Min(" << alpha1.name() << ") = " << min(alpha1).value()
    << "  Max(" << alpha1.name() << ") = " << max(alpha1).value()
    << endl;
