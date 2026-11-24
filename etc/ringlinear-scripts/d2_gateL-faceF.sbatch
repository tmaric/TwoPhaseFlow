#!/bin/bash
#SBATCH -A special00004
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --mem-per-cpu=3600
#SBATCH -t 01:00:00
#SBATCH -o %x.%j.out
#SBATCH -e %x.%j.err
# Part F (2026-11-24): row ringSum of the gate (3D sphere N 64, standard centre, one step), 1 rank and 8 ranks (same decomposition as job 57168520), with the print-only key
# davofFaceSourceLog true (the registered field alphafSource.water is written with the fields); the record of the processor faces whose end-time alphaf differs from serial.
W=/work/scratch/tm83tomy/davof/agents/forkdefault
TPF=$W/TwoPhaseFlow
module purge; module load gcc/11.5.0-z7mc openmpi/4.1.8-6xzv
source "$HOME/OpenFOAM/OpenFOAM-v2606/etc/bashrc"
. $TPF/etc/davof-env.sh
echo "binary $(sha256sum $(command -v advectorVoF) | cut -c1-16) libVoF $(sha256sum $TPF/platforms/$WM_OPTIONS/lib/libVoF.so | cut -c1-16) job $SLURM_JOB_ID"
G=$W/gateL/3D_N64_p0
for r in 1 8; do
  P=$W/gateL/ctlF_$r; rm -rf $P; mkdir -p $P
  cp -r $G/B/system $G/B/constant $G/B/0 $P/
  sed -i 's/davofProcFaceRule symmetric;/davofProcFaceRule symmetric; davofFaceSourceLog true;/' $P/system/fvSolution
  grep -c "davofFaceSourceLog true" $P/system/fvSolution
  cd $P
  if [ $r = 1 ]; then
    advectorVoF > log.advectorVoF 2>&1; echo "1 rank run rc=$? end=$(grep -c '^End' log.advectorVoF)"
  else
    foamDictionary -entry numberOfSubdomains -set 8 system/decomposeParDict > /dev/null
    foamDictionary -entry method -set scotch system/decomposeParDict > /dev/null
    decomposePar -force > log.decomposePar 2>&1; echo "decompose rc=$?"
    mpirun -np 8 advectorVoF -parallel > log.advectorVoF 2>&1; echo "8 ranks run rc=$? end=$(grep -c '^End' log.advectorVoF)"
  fi
done
ls $W/gateL/ctlF_1/5.20833e-06 | tr '\n' ' '
/usr/bin/python3 $W/d2_faceRecord.py $W/gateL/ctlF_1 $W/gateL/ctlF_8 5.20833e-06 > $W/faceRecord.md 2> $W/faceRecord.err
echo "record rc=$?"; tail -3 $W/faceRecord.err
cat $W/faceRecord.md
