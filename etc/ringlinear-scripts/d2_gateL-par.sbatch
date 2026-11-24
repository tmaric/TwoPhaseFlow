#!/bin/bash
#SBATCH -A special00004
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --mem-per-cpu=3600
#SBATCH -t 01:00:00
#SBATCH -o %x.%j.out
#SBATCH -e %x.%j.err
# 1 rank against 8 ranks, 3D sphere N 64 standard centre, ringLinear 1 pass (2026-11-24): the case of the serial row E1 of gateL, decomposed (scotch), the same step.
W=/work/scratch/tm83tomy/davof/agents/forkdefault
TPF=$W/TwoPhaseFlow
module purge; module load gcc/11.5.0-z7mc openmpi/4.1.8-6xzv
source "$HOME/OpenFOAM/OpenFOAM-v2606/etc/bashrc"
. $TPF/etc/davof-env.sh
S=$W/gateL/3D_N64_p0/E1
P=$W/gateL/par8_N64
rm -rf $P; mkdir -p $P
cp -r $S/system $S/constant $S/0 $P/
cd $P
foamDictionary -entry numberOfSubdomains -set 8 system/decomposeParDict > /dev/null
foamDictionary -entry method -set scotch system/decomposeParDict > /dev/null
decomposePar -force > log.decomposePar 2>&1
echo "decompose rc=$?"
mpirun -np 8 advectorVoF -parallel > log.advectorVoF 2>&1
echo "run rc=$? end=$(grep -c '^End' log.advectorVoF)"
printf "binary %s\njob %s\n" "$(sha256sum $(command -v advectorVoF) | cut -c1-16)" "$SLURM_JOB_ID" > .stamp
echo "serial : $(grep -v '^#' $S/postProcessing/reconstructionError/*/error.dat | tail -1 | cut -f1-3,13-15)"
echo "8 ranks: $(grep -v '^#' $P/postProcessing/reconstructionError/*/error.dat | tail -1 | cut -f1-3,13-15)"
grep -h "ringLinear pass" $S/log.advectorVoF | tail -1
grep -h "ringLinear pass" $P/log.advectorVoF | tail -1
