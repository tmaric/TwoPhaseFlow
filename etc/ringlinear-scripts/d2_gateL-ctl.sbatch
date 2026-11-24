#!/bin/bash
#SBATCH -A special00004
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --mem-per-cpu=3600
#SBATCH -t 01:00:00
#SBATCH -o %x.%j.out
#SBATCH -e %x.%j.err
# The controls of the 1-rank against 8-rank difference of ringLinear (2026-11-24): 8 ranks x {plicRDF A, ringSum B, ringLinear correction off E0}, at 3D sphere N 64,
# standard centre, one step of deltaT/1000; the cases are those of the serial rows of gateL (A, B, E1) with the same decomposition as job 57155296 (scotch, 8),
# E0 = the E1 case with ringLinearCorrection false. Also: the decomposition and the cells of the interface that have a face on a processor boundary.
W=/work/scratch/tm83tomy/davof/agents/forkdefault
TPF=$W/TwoPhaseFlow
module purge; module load gcc/11.5.0-z7mc openmpi/4.1.8-6xzv
source "$HOME/OpenFOAM/OpenFOAM-v2606/etc/bashrc"
. $TPF/etc/davof-env.sh
echo "binary $(sha256sum $(command -v advectorVoF) | cut -c1-16) libVoF $(sha256sum $TPF/platforms/$WM_OPTIONS/lib/libVoF.so | cut -c1-16) job $SLURM_JOB_ID"
G=$W/gateL/3D_N64_p0
for r in A B E0; do
  P=$W/gateL/ctl8_$r; rm -rf $P; mkdir -p $P
  src=$G/$r; [ $r = E0 ] && src=$G/E1
  cp -r $src/system $src/constant $src/0 $P/
  if [ $r = E0 ]; then sed -i 's/ringLinearCorrection true;/ringLinearCorrection false;/' $P/system/fvSolution; grep -c "ringLinearCorrection false" $P/system/fvSolution; fi
  cd $P
  foamDictionary -entry numberOfSubdomains -set 8 system/decomposeParDict > /dev/null
  foamDictionary -entry method -set scotch system/decomposeParDict > /dev/null
  decomposePar -force > log.decomposePar 2>&1; echo "$r decompose rc=$?"
  mpirun -np 8 advectorVoF -parallel > log.advectorVoF 2>&1; echo "$r run rc=$? end=$(grep -c '^End' log.advectorVoF)"
done
# the decomposition and the interface
P=$W/gateL/ctl8_A
cd $P
cat system/decomposeParDict | grep -v "^\s*//" | sed -n 8,40p
cp -r $W/gateL/3D_N64_p0/A/system $W/gateL/3D_N64_p0/A/constant $W/gateL/3D_N64_p0/A/0 $W/gateL/ctl8_cell/ 2>/dev/null
rm -rf $W/gateL/ctl8_cell; mkdir -p $W/gateL/ctl8_cell; cp -r $W/gateL/3D_N64_p0/A/system $W/gateL/3D_N64_p0/A/constant $W/gateL/3D_N64_p0/A/0 $W/gateL/ctl8_cell/
cd $W/gateL/ctl8_cell
foamDictionary -entry numberOfSubdomains -set 8 system/decomposeParDict > /dev/null
foamDictionary -entry method -set scotch system/decomposeParDict > /dev/null
decomposePar -force -cellDist > log.decomposePar 2>&1; echo "cellDist decompose rc=$?"
ls 0 | head; ls constant/polyMesh | head
python3 $W/d2_ifcross.py $W/gateL/ctl8_cell
