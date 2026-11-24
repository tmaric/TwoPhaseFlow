#!/bin/bash
#SBATCH -A special00004
#SBATCH -N 1
#SBATCH --exclusive
#SBATCH --constraint=i01
#SBATCH --mem-per-cpu=3000
#SBATCH -o %x.%j.out
#SBATCH -e %x.%j.err
# The timed verification of the isoAdvector package (fork branch davof-isoAdvector,
# 2026-11-06; the pre-registration F-1 .. F-5 in config/davof/iso*3D_hexN32.yaml). One
# exclusive i01 node per (FLOW, N, CO index): REPS repetitions of the rows RDF (plicRDF,
# vofCutCell), RDFccs (plicRDF, ccs), RING (davof ring normal, ccs), RINGvcc (davof ring
# normal, vofCutCell) and ISO (isoAlpha), the order rotated per repetition, each solve one
# srun step on NP ranks (--cpu-bind=cores); the case copies come from the preparation
# study iso<FLOW>3D_hexN<N> (rendered, meshed, initialised, decomposed), endTime set to the
# flow's T (the rotation 1, the shear and the deformation 3), the fields and the
# diagnostics written at the quarter times. After the repetitions ONE profiled run per
# row (OpenFOAM profiling active: the reconstruct()/advect() shares), not timed.
#   sbatch -J davof-iso-<FLOW>-N<N>-co<CO> -n <NP> -t <limit> --export=ALL,FLOW=Rotation|Shear|Deformation,N=..,CO=0|1|2 iso-time.sbatch
# REP0 / REP1 (optional): only the repetitions REP0..REP1; METHODS (optional): only these
# rows; NOPROF=1: no profiled runs. ISO_W: the workspace (default the swept agent's iso/).
# MESHCLASS (default hex): hex | perturbed | poly | tet, the preparation study
# iso<FLOW>3D_<MESHCLASS>N<N> and the output directory time/<FLOW>_<MESHCLASS>_N<N>_co<CO>.
W=${ISO_W:-/work/scratch/tm83tomy/davof/agents/swept/iso}
S=${STUDIES:-$W/leia/studies/davof}
OW=${OUT_W:-/work/scratch/tm83tomy/davof/agents/forkdefault/iso}   # the outputs of this task (ringsum matrix, 2026-11-23); W only supplies the prepared cases and the scripts, read only
module purge
module load gcc/11.5.0-z7mc openmpi/4.1.8-6xzv
source "$HOME/OpenFOAM/OpenFOAM-v2606/etc/bashrc"
# the fork worktree whose binaries run (TPF_DIR; default iso/TwoPhaseFlow): leia-env.sh
# appends its platforms to PATH and LD_LIBRARY_PATH
export LEIA_TPF_DIR=${TPF_DIR:-$W/TwoPhaseFlow}
cd $W/leia && . etc/leia-env.sh || exit 1
: ${FLOW:?} ${N:?} ${CO:?}
MESHCLASS=${MESHCLASS:-hex}
REPS=${REPS:-3}
NP=$SLURM_NTASKS
case $FLOW in Rotation) END=1; WI=0.25;; Shear|Deformation) END=3; WI=0.75;; Shear2D) END=8; WI=2;; *) echo "FLOW $FLOW?"; exit 1;; esac
# the preparation study: iso<FLOW>3D_<class>N<N>, or iso<FLOW>_... for a 2D flow
STB=iso${FLOW}3D; case $FLOW in *2D) STB=iso$FLOW;; esac
echo "fork $(git -C $LEIA_TPF_DIR log --oneline -1 | cut -c1-90) ($LEIA_TPF_DIR); leia $(git -C $W/leia log --oneline -1 | cut -c1-90)"
echo "node $(hostname) $(lscpu | grep -m1 'Model name' | sed 's/  */ /g'), cores $(nproc --all), job $SLURM_JOB_ID, np $NP, FLOW $FLOW MESHCLASS $MESHCLASS N $N CO index $CO, $(date)"
echo "advectorVoF $(command -v advectorVoF)"
# the stamp of every run (.stamp): the fork commit and worktree, the sha256 of the binary and
# of libVoF, the leia commit
FORK_SHA=$(git -C $LEIA_TPF_DIR rev-parse --short HEAD); LEIA_SHA=$(git -C $W/leia rev-parse --short HEAD)
BIN_SHA=$(sha256sum $(command -v advectorVoF) | cut -c1-16); LIB_SHA=$(sha256sum $LEIA_TPF_DIR/platforms/$WM_OPTIONS/lib/libVoF.so | cut -c1-16)
T=$OW/time/${FLOW}_${MESHCLASS}_N${N}_co${CO}
# METHODS alone (REP0 empty) re-runs those rows of an existing cell and keeps the others
if [ -z "$REP0" ] && [ -z "$METHODS" ]; then rm -rf $T
elif [ -n "$METHODS" ]; then for m in $METHODS; do rm -rf $T/${m}_r${REP0:-*}; done
else rm -rf $T/*_r$REP0; fi
mkdir -p $T
SUF=${REP0:+_r$REP0}${METHODS:+_${METHODS// /-}}
pick() { /usr/bin/python3 - "$S/${STB}_${MESHCLASS}N$N" "$1" "$CO" <<'EOF'
import json, glob, os, sys
st, row, co = sys.argv[1], sys.argv[2], int(sys.argv[3])
cs = []
for c in sorted(glob.glob(st + "/*_0*")):
    if "." in os.path.basename(c) or not os.path.exists(c + "/case_params.json"): continue
    t = json.load(open(c + "/case_params.json"))["tokens"]
    if t["TPF_ROW"] == row: cs.append((float(t["TARGET_ADVECTIVE_CO"]), c))
cs.sort()
sys.stderr.write("pick %s: Co %s -> %s\n" % (row, [x[0] for x in cs], cs[co][0]))
print(cs[co][1])
EOF
}
declare -A SRC
SRC[RDF]=$(pick plicRDF); SRC[RDFccs]=$(pick plicRDF-ccs); SRC[RING]=$(pick davof-ringArea)
SRC[RINGvcc]=$(pick davof-ringArea-vofCutCell); SRC[ISO]=$(pick isoAlpha)
ALL="RDF PURE ONERING ENTERING"
# the row PURE ("pure" DAVOF, Tomislav 2026-11-08; opt-in through METHODS until it is the default):
# the RING case copy with the keys appended to the rendered fvSolution by sed (a foamDictionary
# rewrite perturbs a run at round-off); needs a build with the keys (TPF_DIR, 7646055 or later)
SRC[PURE]=${SRC[RING]}
PURE_KEYS='davofFluxExtension none; davofFaceFractionSource sweepOnly; boundSecondLayer false; sweptBoundedFaceFraction scale; interfaceCellRule alphaAndRingSum; normalModel ringSum; wispTol 0; davofProcFaceRule symmetric;'
# the row ONERING (Tomislav's hypothesis of 2026-11-09: the third face-fraction source interfaceCellFaces,
# every face of every interface cell swept by the moving plane, no second ring; the other keys PURE's;
# opt-in through METHODS; needs a build with the value, the fork after 7646055)
SRC[ONERING]=${SRC[RING]}
ONERING_KEYS='davofFluxExtension none; davofFaceFractionSource interfaceCellFaces; boundSecondLayer false; sweptBoundedFaceFraction scale; interfaceCellRule alphaAndRingSum; normalModel ringSum; wispTol 0; davofProcFaceRule symmetric;'
# the row ENTERING (Tomislav's realizability rule of 2026-11-09: the fourth face-fraction source enteringFaces,
# the sweeps of interfaceCellFaces plus the second-ring sweep restricted to the entering cells and the guard
# that dries the faces of a dry cell and fills those of a full one; the other keys PURE's; opt-in through
# METHODS; needs the fork after 294891d)
SRC[ENTERING]=${SRC[RING]}
SRC[CTRL]=${SRC[RING]}
SRC[HYB1]=${SRC[RING]}
SRC[QUAD]=${SRC[RING]}
SRC[LIN]=${SRC[RING]}
CTRL_KEYS='davofFluxExtension none; davofFaceFractionSource sweepOnly; boundSecondLayer false; sweptBoundedFaceFraction scale; interfaceCellRule alphaAndRingSum; normalModel plicRDF; wispTol 0; davofProcFaceRule symmetric;'
HYB1_KEYS="${PURE_KEYS/normalModel ringSum;/normalModel ringSumRDF; rdfIterations 1;}"
QUAD_KEYS="${PURE_KEYS/normalModel ringSum;/normalModel ringFit; applyTo all; weight area; alignTol 0.2; minCells 6;}"
LIN_KEYS="${PURE_KEYS/normalModel ringSum;/normalModel ringLinear; ringLinearPasses 1; ringLinearCorrection true; ringLinearDiag true;}"
ENTERING_KEYS='davofFluxExtension none; davofFaceFractionSource enteringFaces; boundSecondLayer false; sweptBoundedFaceFraction scale; interfaceCellRule alphaAndRingSum; normalModel ringSum; wispTol 0; davofProcFaceRule symmetric;'
for m in $ALL; do echo "source $m ${SRC[$m]}"; done
prep() { # method rep [prof] -> case dir
  local m=$1 r=$2 src=${SRC[$1]} c=$T/${1}_r$2
  mkdir -p $c; cp -r $src/system $src/constant $c/
  if [ $NP -gt 1 ]; then
    for p in $src/processor*; do mkdir -p $c/$(basename $p); cp -r $p/0 $p/constant $c/$(basename $p)/; done
  else
    cp -r $src/0 $c/
  fi
  foamDictionary -entry endTime -set $END $c/system/controlDict > /dev/null
  foamDictionary -entry writeInterval -set $WI $c/system/controlDict > /dev/null
  foamDictionary -entry functions/areaFractionError/writeInterval -set $WI $c/system/controlDict > /dev/null
  foamDictionary -entry functions/volumeFractionError/writeInterval -set $WI $c/system/controlDict > /dev/null
  [ "$3" = prof ] && foamDictionary -entry profiling/active -set true $c/system/controlDict > /dev/null
  if [ "$m" = PURE ] || [ "$m" = ONERING ] || [ "$m" = ENTERING ] || [ "$m" = CTRL ] || [ "$m" = HYB1 ] || [ "$m" = QUAD ] || [ "$m" = LIN ]; then
    local keys=$PURE_KEYS; [ "$m" = ONERING ] && keys=$ONERING_KEYS; [ "$m" = ENTERING ] && keys=$ENTERING_KEYS; [ "$m" = CTRL ] && keys=$CTRL_KEYS; [ "$m" = HYB1 ] && keys=$HYB1_KEYS; [ "$m" = QUAD ] && keys=$QUAD_KEYS; [ "$m" = LIN ] && keys=$LIN_KEYS
    sed -i -E "s/^(\s*)(normalModel\s+[A-Za-z]+;)/\1\/\/ \2 ($m)/; s/^(\s*)(wispTol\s+[0-9.e-]+;)/\1\/\/ \2 ($m)/" $c/system/fvSolution
    sed -i "0,/^\(\s*\)reconstructionScheme davof;/s||&\n\1$keys|" $c/system/fvSolution
  fi
  if [ -n "$EXTRA" ]; then sed -i "0,/^\(\s*\)reconstructionScheme [A-Za-z]*;/s||&\n\1$EXTRA|" $c/system/fvSolution; fi
  if [ -n "$NAB" ]; then sed -i -E "s/^(\s*)nAlphaBounds\s+[0-9]+;/\1nAlphaBounds $NAB;/" $c/system/fvSolution; fi
  if [ -n "$D1" ]; then
    local dt0=$(foamDictionary -entry deltaT -value $c/system/controlDict); local dtn=$(python3 -c "print(repr($dt0/1000.0))")
    foamDictionary -entry deltaT -set $dtn $c/system/controlDict > /dev/null
    foamDictionary -entry endTime -set $dtn $c/system/controlDict > /dev/null
    foamDictionary -entry writeInterval -set $dtn $c/system/controlDict > /dev/null
    foamDictionary -entry functions/areaFractionError/writeInterval -set $dtn $c/system/controlDict > /dev/null
    foamDictionary -entry functions/volumeFractionError/writeInterval -set $dtn $c/system/controlDict > /dev/null
  fi
  echo "$src" > $c/.source; echo $c
}
csv=$T/timing$SUF.csv; echo "flow,N,np,co_index,method,rep,steps,dt,end_time,clock_s,exec_s,wall_s,median_step_s,node" > $csv
par=""; [ $NP -gt 1 ] && par="-parallel"
run1() { # method rep [prof]
  local m=$1 r=$2 c; c=$(prep $m $r $3)
  local t0=$(date +%s.%N)
  (cd $c && srun --ntasks=$NP --cpu-bind=cores --kill-on-bad-exit=1 --job-name=${m}_r$r advectorVoF $par > log.advectorVoF 2>&1); local rc=$?
  local t1=$(date +%s.%N)
  printf "fork_commit=%s\nfork_dir=%s\nadvectorVoF_sha256_16=%s\nlibVoF_sha256_16=%s\nleia_commit=%s\njob=%s\n" "$FORK_SHA" "$LEIA_TPF_DIR" "$BIN_SHA" "$LIB_SHA" "$LEIA_SHA" "$SLURM_JOB_ID" > $c/.stamp
  (cd $c && python3 $W/leia/workflow/scripts/tpf_advection_metrics.py --case . --out advectorVoF.csv > log.metrics 2>&1)
  /usr/bin/python3 - $c $m $r $rc $t0 $t1 "$FLOW" "$N" "$NP" "$CO" >> $csv <<'EOF'
import sys, re, os, statistics
c, m, r, rc, t0, t1, flow, N, np_, co = sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4]), float(sys.argv[5]), float(sys.argv[6]), sys.argv[7], sys.argv[8], sys.argv[9], sys.argv[10]
log = open(c + "/log.advectorVoF").read()
ct = re.findall(r"ExecutionTime = (\S+) s\s+ClockTime = (\S+) s(?:\s+timeStepTime = (\S+) s)?", log)
steps = len(re.findall(r"\nTime = ", log)); st = [float(x[2]) for x in ct if x[2]]
cd = open(c + "/system/controlDict").read()
dt = re.search(r"\ndeltaT\s+([^;]+);", cd).group(1); et = re.search(r"\nendTime\s+([^;]+);", cd).group(1)
print(",".join(str(x) for x in (flow, N, np_, co, m, r, steps, dt, et, ct[-1][1] if ct else "", ct[-1][0] if ct else "", "%.2f" % (t1 - t0),
      "%.4f" % statistics.median(st) if st else "", os.uname().nodename)) + ("" if rc == 0 else ",rc=%d" % rc))
EOF
  echo "$m r$r rc=$rc wall $(python3 -c "print('%.1f' % ($t1-$t0))") s, $(grep -c '^Time = ' $c/log.advectorVoF) steps; $(grep -o 'ABORT at step [0-9]*, t = [^:]*: [a-zA-Z() ]*' $c/log.advectorVoF | head -1)"
}
if [ -z "$REP0" ] || [ -z "$METHODS" ]; then METHODS=${METHODS:-$ALL}; fi
# REP0=prof: the profiled runs only (one per row, after the timed repetitions of another job)
if [ "$REP0" = prof ]; then
  for m in $METHODS; do run1 $m prof prof; done
  sleep 5; sacct -j $SLURM_JOB_ID -n -P -o JobID,JobName,Elapsed,ElapsedRaw,NNodes,NTasks,State > $T/sacct_steps_prof.txt
  echo "ISO-TIME-DONE $FLOW $MESHCLASS N$N co$CO (profiled runs) $(date)"; exit 0
fi
for r in $(seq ${REP0:-1} ${REP1:-$REPS}); do
  # the order rotated per repetition
  set -- $METHODS; k=$(( (r - 1) % $# )); order="${@:$((k+1))} ${@:1:$k}"
  for m in $order; do
    # one repetition per aborted cell (Tomislav, 2026-11-06): a method whose first
    # repetition ended by the abort guard is not repeated
    if [ "$r" != 1 ] && [ "$r" != prof ] && grep -q 'advectorVoF: ABORT' $T/${m}_r1/log.advectorVoF 2>/dev/null; then
      echo "$m r$r skipped (r1 aborted: $(grep -o 'ABORT at step [0-9]*' $T/${m}_r1/log.advectorVoF | head -1))"; continue
    fi
    run1 $m $r
  done
done
if [ -z "$NOPROF" ] && [ -z "$REP0" ]; then
  for m in $METHODS; do run1 $m prof prof; done
fi
sleep 5; sacct -j $SLURM_JOB_ID -n -P -o JobID,JobName,Elapsed,ElapsedRaw,NNodes,NTasks,State > $T/sacct_steps$SUF.txt
echo "ISO-TIME-DONE $FLOW $MESHCLASS N$N co$CO $(date)"
