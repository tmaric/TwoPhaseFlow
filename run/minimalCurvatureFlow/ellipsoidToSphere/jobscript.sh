# job submission
subopenfoam -n 1 \
            -N 1 \
            -r 5:30 \
            -T $HOME/OpenFOAM/OpenFOAM-v2212 \
            -j squareCapillary -V v2212 interFlow 



