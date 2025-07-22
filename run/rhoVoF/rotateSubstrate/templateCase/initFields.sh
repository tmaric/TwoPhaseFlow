#!/bin/sh
# NOTE: it is not possible to prescribe a specific shell for interpretation
# of this script since all comments are lost when the template parameters
# are replaced.
# So this script should be interpretable by all POSIX compatible shell
# cd ${0%/*} || exit 1    # run from this directory

# Source tutorial run functions
. $WM_PROJECT_DIR/bin/tools/RunFunctions

./meshGen.sh

runApplication setAlphaField
runApplication setExprFields

runApplication decomposePar
