# overlap (vendored)

`overlap.hpp` is the reference implementation of

S. Strobl, A. Formella, T. Pöschel, Exact calculation of the overlap volume of
spheres and mesh elements, Journal of Computational Physics 311 (2016) 158-172,
doi:10.1016/j.jcp.2016.02.003,

taken unchanged from https://github.com/severinstrobl/overlap (MIT licence, see
LICENSE; CITATION.cff as shipped) on 2026-10-24. It is header-only and needs
Eigen 3 (`-I$(EIGEN_INC)` in `src/VoF/Make/options`; `etc/davof-env.sh` sets
`EIGEN_INC` when it is not set). Used by the volume rule `exactSphere`
(`src/VoF/davof/volumeRules/exactSphereVolumeRule.C`).
