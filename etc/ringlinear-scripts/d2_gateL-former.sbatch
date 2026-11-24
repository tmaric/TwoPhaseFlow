#!/bin/bash
#SBATCH -A special00004
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --mem-per-cpu=3600
#SBATCH -t 01:00:00
#SBATCH -o %x.%j.out
#SBATCH -e %x.%j.err
# 8 ranks, row ringSum (B) of the gate, 3D sphere N 64 standard centre, with davofProcFaceRule former (all other keys unchanged); same decomposition as job 57168520.
W=/work/scratch/tm83tomy/davof/agents/forkdefault
TPF=$W/TwoPhaseFlow
module purge; module load gcc/11.5.0-z7mc openmpi/4.1.8-6xzv
source "$HOME/OpenFOAM/OpenFOAM-v2606/etc/bashrc"
. $TPF/etc/davof-env.sh
echo "binary $(sha256sum $(command -v advectorVoF) | cut -c1-16) libVoF $(sha256sum $TPF/platforms/$WM_OPTIONS/lib/libVoF.so | cut -c1-16) job $SLURM_JOB_ID"
G=$W/gateL/3D_N64_p0
P=$W/gateL/ctl8_Bformer; rm -rf $P; mkdir -p $P
cp -r $G/B/system $G/B/constant $G/B/0 $P/
grep -c "davofProcFaceRule symmetric" $P/system/fvSolution
sed -i 's/davofProcFaceRule symmetric;/davofProcFaceRule former;/' $P/system/fvSolution
grep -o "davofProcFaceRule [a-z]*;" $P/system/fvSolution
cd $P
foamDictionary -entry numberOfSubdomains -set 8 system/decomposeParDict > /dev/null
foamDictionary -entry method -set scotch system/decomposeParDict > /dev/null
decomposePar -force > log.decomposePar 2>&1; echo "decompose rc=$?"
mpirun -np 8 advectorVoF -parallel > log.advectorVoF 2>&1; echo "run rc=$? end=$(grep -c '^End' log.advectorVoF)"
grep -v '^#' postProcessing/reconstructionError/*/error.dat | tail -1 | cut -f13-15,18
