# Python toolchain for the plotting job on Lichtenberg. Sourced, not executed.
#
# Deliberately separate from env-lichtenberg.sh: that one starts with
# `module purge`, which removes the python module carrying numpy, matplotlib
# and foamlib. The snapshot rule needs those and does NOT need OpenFOAM.
module purge
module load python          # python/3.11.14 -> numpy, matplotlib, foamlib
