#!/bin/bash
#SBATCH -A special00004
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -o %x.%j.out
#SBATCH -e %x.%j.err
# The static normal gate (2026-11-23, branch davof-isoAdvector-normal): one state (DIM 3|2, N, POS 0|1 = the standard or the shifted centre),
# serial, every row after the other. The case: system and constant from the prepared study case of the row's family (the shear3D / shear2D study
# of the 2026-11-09 matrix, target Co 0.5; N 512 of the 2D circle from the N256 study with n_cells 512), blockMesh, the exact initial state by
# initDavofFields (with the shifted centre in POS 1), one step of deltaT/1000, reconstructionError at that time. Rows: A plicRDF, B ringSum (PURE),
# C ringFit (the ring quadric), D0 / D1 / D2 ringSumRDF with rdfIterations 0 / 1 / 2 (ROWS selects). Every davof key explicit.
#   sbatch -J gate-3-32-0 --mem=8G -t 02:00:00 --export=ALL,DIM=3,N=32,POS=0,ROWS="A B C D1 D2" d2_gate.sbatch
W=/work/scratch/tm83tomy/davof/agents/forkdefault
S=/work/scratch/tm83tomy/davof/agents/swept/iso/leia/studies/davof
TPF=${TPF_DIR:-$W/TwoPhaseFlow}
OUT=${OUT_W:-$W/gate}
: ${DIM:?} ${N:?} ${POS:?}
ROWS=${ROWS:-"A B C D1 D2"}
module purge; module load gcc/11.5.0-z7mc openmpi/4.1.8-6xzv
source "$HOME/OpenFOAM/OpenFOAM-v2606/etc/bashrc"
. $TPF/etc/davof-env.sh
echo "binary $(command -v advectorVoF); $(git -C $TPF log --oneline -1 | cut -c1-80); job $SLURM_JOB_ID; $(date)"
BIN_SHA=$(sha256sum $(command -v advectorVoF) | cut -c1-16); LIB_SHA=$(sha256sum $TPF/platforms/$WM_OPTIONS/lib/libVoF.so | cut -c1-16)
if [ $DIM = 3 ]; then ST=isoShear3D; NS=$N; OLD="[(]0\.5 0\.75 0\.25[)]"; NEW="(0.513 0.757 0.253)"
else ST=isoShear2D; NS=$N; [ $N = 512 ] && NS=256; OLD="[(]0\.5 0 0\.75[)]"; NEW="(0.513 0.007 0.75)"; fi
pick() { /usr/bin/python3 - "$S/${ST}_hexN$NS" "$1" <<'EOF'
import json, glob, os, sys
st, row = sys.argv[1], sys.argv[2]
for c in sorted(glob.glob(st + "/*_0*")):
    if "." in os.path.basename(c) or not os.path.exists(c + "/case_params.json"): continue
    t = json.load(open(c + "/case_params.json"))["tokens"]
    if t["TPF_ROW"] == row and float(t["TARGET_ADVECTIVE_CO"]) == 0.5: print(c); break
EOF
}
SRC_RDF=$(pick plicRDF); SRC_RING=$(pick davof-ringArea)
echo "sources: $SRC_RDF $SRC_RING"
B=$OUT/${DIM}D_N${N}_p${POS}
rm -rf $B; mkdir -p $B
# the base: the RING case (system with the davof dictionary), mesh, the exact initial state
base=$B/base; mkdir -p $base/constant; cp -r $SRC_RING/system $base/; cp $SRC_RING/constant/transportProperties $base/constant/
cp -r $SRC_RING/0.org $base/0
if [ $N != $NS ]; then sed -i -E "s/^n_cells +[0-9]+;/n_cells $N;/" $base/system/blockMeshDict; fi
if [ $POS = 1 ]; then
  for f in system/initAlphaFieldDict system/fvSolution system/controlDict; do sed -i "s/$OLD/$NEW/g" $base/$f; done
fi
grep -n "origin\|axisPoint" $base/system/initAlphaFieldDict $base/system/controlDict | head
(cd $base && blockMesh > log.blockMesh 2>&1 && initDavofFields > log.initDavofFields 2>&1; echo "init rc=$?"; tail -3 log.initDavofFields | cut -c1-150)
dt0=$(foamDictionary -entry deltaT -value $SRC_RING/system/controlDict)
if [ $N != $NS ]; then dt0=$(python3 -c "print(repr($dt0*$NS/$N))"); fi
dtn=$(python3 -c "print(repr($dt0/1000.0))")
echo "deltaT0 $dt0, one step $dtn"
PURE_KEYS="davofFluxExtension none; davofFaceFractionSource sweepOnly; boundSecondLayer false; sweptBoundedFaceFraction scale; interfaceCellRule alphaAndRingSum; wispTol 0; davofProcFaceRule symmetric;"
for r in $ROWS; do
  c=$B/$r; mkdir -p $c; cp -r $base/system $base/constant $base/0 $c/
  # the step: one step of deltaT/1000, written at that time
  foamDictionary -entry deltaT -set $dtn $c/system/controlDict > /dev/null
  foamDictionary -entry endTime -set $dtn $c/system/controlDict > /dev/null
  foamDictionary -entry writeInterval -set $dtn $c/system/controlDict > /dev/null
  foamDictionary -entry functions/areaFractionError/writeInterval -set $dtn $c/system/controlDict > /dev/null
  foamDictionary -entry functions/volumeFractionError/writeInterval -set $dtn $c/system/controlDict > /dev/null
  sed -i "s/^\(\s*\)reconstructionError true;/&\n\1normalCellDiag true;/" $c/system/controlDict
  case $r in
    A) cp $SRC_RDF/system/fvSolution $c/system/fvSolution
       [ $POS = 1 ] && sed -i "s/$OLD/$NEW/g" $c/system/fvSolution
       # print-only: the plicRDF iterations of the step
       printf '\nDebugSwitches\n{\n    plicRDF 1;\n}\n' >> $c/system/controlDict ;;
    *) case $r in
         B) NM="normalModel ringSum;" ;;
         C) NM="normalModel ringFit; applyTo all; weight area; alignTol 0.2; minCells 6; ringFitDiag 1;" ;;
         E0) NM="normalModel ringLinear; ringLinearPasses 1; ringLinearCorrection false; ringLinearDiag true;" ;;
         E1) NM="normalModel ringLinear; ringLinearPasses 1; ringLinearCorrection true; ringLinearDiag true;" ;;
         E2) NM="normalModel ringLinear; ringLinearPasses 2; ringLinearCorrection true; ringLinearDiag true;" ;;
         D0) NM="normalModel ringSumRDF; rdfIterations 0;" ;;
         D1) NM="normalModel ringSumRDF; rdfIterations 1;" ;;
         D2) NM="normalModel ringSumRDF; rdfIterations 2;" ;;
       esac
       sed -i -E "s/^(\s*)(normalModel\s+[A-Za-z]+;)/\1\/\/ \2 ($r)/; s/^(\s*)(wispTol\s+[0-9.e-]+;)/\1\/\/ \2 ($r)/" $c/system/fvSolution
       sed -i "0,/^\(\s*\)reconstructionScheme davof;/s||&\n\1$PURE_KEYS $NM|" $c/system/fvSolution ;;
  esac
  ( cd $c && advectorVoF > log.advectorVoF 2>&1; rc=$?
    printf "fork_commit=%s\nadvectorVoF_sha256_16=%s\nlibVoF_sha256_16=%s\njob=%s\n" "$(git -C $TPF rev-parse --short HEAD)" "$BIN_SHA" "$LIB_SHA" "$SLURM_JOB_ID" > .stamp
    echo "GATE $DIM D N$N p$POS $r rc=$rc steps $(grep -c '^Time = ' log.advectorVoF) end=$(grep -c '^End' log.advectorVoF) $(grep -o 'Courant Number mean: [^ ]* max: [^ ]*' log.advectorVoF | tail -1)" )
done
echo "GATE-DONE $(date)"
