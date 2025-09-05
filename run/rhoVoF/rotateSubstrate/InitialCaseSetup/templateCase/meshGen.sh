#!/bin/sh
# NOTE: it is not possible to prescribe a specific shell for interpretation
# of this script since all comments are lost when the template parameters
# are replaced.
# So this script should be interpretable by all POSIX compatible shell

#gmsh constant/geometry/squareSubstrate.geo -parse_and_exit
#surfaceFeatureEdges constant/geometry/squareSubstrate.stl constant/geometry/squareSubstrate.fms
#cartesianMesh
blockMesh
# topoSet
# refineHexMesh -overwrite interfaceRegion
