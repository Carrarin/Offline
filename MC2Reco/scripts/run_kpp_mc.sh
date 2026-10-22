#!/bin/bash
# run_kpp_mc.sh -- digitization, reconstruction and ntuple of the KPP hall simulation (gm2ita), file by file as the Geant4 jobs finish:
# every dts.mc2.<NAME>.<subrun>.art of SIM whose job log ends with art status 0 and that has not changed for 2 minutes goes through
#   mc2_kpp_digi_reco.fcl -> products/kpp_mc_reco/<NAME>/<ver>/rec_<subrun>.art
#   kpp_ntuple.fcl        -> products/kpp_mc_ntuple/<NAME>/<ver>/nts_<subrun>.root
#   kpp_dump_events.py    -> products/kpp_mc_events/<NAME>/<ver>/ev_<subrun>.json
# at most NJOBS at a time (SCHED_IDLE, nice 19), with a MANIFEST.json per product. Runs until STOPFILE exists and nothing is left.
# usage: run_kpp_mc.sh <muse area> <sim dir> <name, e.g. CosmicCRYHallKPP> <version, e.g. v001>
#   env: NJOBS (8), CHECK (120 s), WORK (/data1/mu2e/carrarin/kpp_mc_<name>_<version>), STOPFILE ($WORK/STOP)
M=$(readlink -f $1); SIM=$(readlink -f $2); NAME=$3; VER=$4
[ -d "$M/Offline/MC2Reco" ] && [ -d "$SIM" ] && [ -n "$VER" ] || { echo "usage: $0 <muse area> <sim dir> <name> <version>"; exit 1; }
NJOBS=${NJOBS:-8}; CHECK=${CHECK:-120}; WORK=${WORK:-/data1/mu2e/carrarin/kpp_mc_${NAME}_$VER}; STOPFILE=${STOPFILE:-$WORK/STOP}
P=/data1/mu2e/mc2/products; IMG=/cvmfs/singularity.opensciencegrid.org/fermilab/fnal-wn-el9:latest
COMMIT=$(git -C $M/Offline rev-parse --short HEAD); DIRTY=$(git -C $M/Offline status --porcelain --untracked-files=no | wc -l)
R=$P/kpp_mc_reco/$NAME/$VER; N=$P/kpp_mc_ntuple/$NAME/$VER; E=$P/kpp_mc_events/$NAME/$VER
mkdir -p $WORK $R $N $E; cd $WORK
log() { echo "$(date '+%F %T') $*" | tee -a $WORK/run.log; }
manifest() {
  for d in $R $N $E; do
    python3 - "$d" <<EOF
import json, glob, os, sys, datetime
d = sys.argv[1]; done = sorted(os.path.basename(f)[:-5] for f in glob.glob("$WORK/*.done"))
m = dict(product=d.replace("$P/", ""), description="KPP hall simulation $NAME: " + {"$R": "digitization and reconstruction (tracker and CRV as data, MC2 calorimeter response)",
         "$N": "EventNtuple", "$E": "tracks with in-time calorimeter crystals and CRV coincidences (kpp_dump_events.py)"}[d],
         fork="github.com/Carrarin/Offline mc2-calo $COMMIT" + (" (local changes)" if int("$DIRTY") else ""),
         environment="muse setup $M (backing DataJob 1.11.0) in fnal-wn-el9", input="$SIM (dts.mc2.$NAME.<subrun>.art, see its MANIFEST.json)",
         fcl=["Offline/MC2Reco/fcl/sim/mc2_kpp_digi_reco.fcl", "Offline/MC2Reco/fcl/kpp/kpp_ntuple.fcl"],
         geometry="Offline/MC2Reco/data/geom/geom_extracted_mc2kpp.txt", command="Offline/MC2Reco/scripts/run_kpp_mc.sh $M $SIM $NAME $VER",
         subruns=done, n_files=len(done), updated=datetime.datetime.now().astimezone().isoformat(timespec="seconds"), author="$(id -un)")
json.dump(m, open(os.path.join(d, "MANIFEST.json"), "w"), indent=1)
EOF
  done
}
job() {   # $1 = subrun (3 digits)
  local s=$1 f=$SIM/dts.mc2.$NAME.$1.art
  chrt --idle 0 nice -n 19 apptainer exec -B /cvmfs,/data1 $IMG bash -c "source /cvmfs/mu2e.opensciencegrid.org/setupmu2e-art.sh >/dev/null 2>&1; cd $M; muse setup > /dev/null 2>&1
    export MU2E_SEARCH_PATH=\$MU2E_SEARCH_PATH:$M:/data1/mu2e/mc2
    mu2e -c Offline/MC2Reco/fcl/sim/mc2_kpp_digi_reco.fcl -s $f -o $R/rec_$s.art > $WORK/rec_$s.log 2>&1 || exit 1
    mu2e -c Offline/MC2Reco/fcl/kpp/kpp_ntuple.fcl -s $R/rec_$s.art -T $N/nts_$s.root > $WORK/nts_$s.log 2>&1 || exit 2
    python3 Offline/MC2Reco/scripts/kpp_dump_events.py $E/ev_$s.json $N/nts_$s.root > $WORK/ev_$s.log 2>&1 || exit 3" < /dev/null
  local st=$?
  if [ $st = 0 ]; then touch $WORK/$s.done; log "subrun $s done ($(tail -1 $WORK/ev_$s.log))"; manifest
  else touch $WORK/$s.failed; log "subrun $s FAILED at step $st, see $WORK"; fi
}
log "start: $SIM ($NAME), version $VER, fork $COMMIT (local changes: $DIRTY), NJOBS $NJOBS"
declare -A PID
while true; do
  for s in "${!PID[@]}"; do kill -0 ${PID[$s]} 2>/dev/null || unset PID[$s]; done
  for f in $(ls $SIM/dts.mc2.$NAME.*.art 2>/dev/null); do
    s=$(basename $f .art); s=${s##*.}; n=$((10#$s))
    [ -e $WORK/$s.done ] || [ -e $WORK/$s.failed ] || [ -n "${PID[$s]}" ] && continue
    [ ${#PID[@]} -ge $NJOBS ] && break
    grep -q "Art has completed and will exit with status 0" $SIM/job_$n.log 2>/dev/null || continue
    [ $(( $(date +%s) - $(stat -c %Y $f) )) -gt 120 ] || continue
    job $s & PID[$s]=$!; log "subrun $s started"
  done
  [ -e $STOPFILE ] && [ ${#PID[@]} = 0 ] && break
  sleep $CHECK
done
log "stopped"
