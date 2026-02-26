static inline scalar twoPi() { return constant::mathematical::twoPi; }

// NOTE:
// This is the MPI acoustic backend used by interMPIALFlow.
// Assembly/solve is split into helper headers to keep the main coupling file
// readable. Block ordering matches freBCMFoam and interFALFlow:
//   rows [0, N)   -> Pim equation
//   rows [N, 2N)  -> Pre equation
//   cols [0, N)   -> Pim unknowns
//   cols [N, 2N)  -> Pre unknowns

#include "acousticAssemblyHelpers.H"
#include "acousticSolveHelpers.H"
