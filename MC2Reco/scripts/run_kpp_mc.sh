#!/bin/bash
#
# run_kpp_mc.sh: digitization, reconstruction and ntuple of the KPP hall simulation on gm2ita, file by file as the
# Geant4 jobs finish.
#
# Every dts.mc2.<name>.<subrun>.art of <sim dir> whose job log (job_<n>.log, n = subrun without leading zeros) ends
# with art status 0 and that has not changed for 2 minutes goes through, in the EL9 container:
#
#   mc2_kpp_digi_reco.fcl -> products/kpp_mc_reco/<name>/<version>/rec_<subrun>.art
#   kpp_ntuple.fcl        -> products/kpp_mc_ntuple/<name>/<version>/nts_<subrun>.root
#   kpp_dump_events.py    -> products/kpp_mc_events/<name>/<version>/ev_<subrun>.json
#
# at most NJOBS at a time (SCHED_IDLE, nice 19). After every subrun the MANIFEST.json of the three products is
# rewritten. The script keeps waiting for new files until STOPFILE exists and no job is running.
#
# Usage (on gm2ita, outside the container; meant to be started detached, e.g. with setsid nohup ... &):
#
#   run_kpp_mc.sh <muse area> <sim dir> <name> <version>
#
#   e.g. Offline/MC2Reco/scripts/run_kpp_mc.sh /data1/mu2e/carrarin/mc2offline \
#          /data1/mu2e/mc2/sim/CosmicCRYHallKPP/v001 CosmicCRYHallKPP v001
#   and, when the simulation is over:  touch /data1/mu2e/carrarin/kpp_mc_CosmicCRYHallKPP_v001/STOP
#
# Arguments:
#
#   <muse area>   muse area with Offline/MC2Reco (built) and EventNtuple
#   <sim dir>     directory of the simulation: dts.mc2.<name>.<subrun>.art and job_<n>.log
#   <name>        name of the production in the file names, e.g. CosmicCRYHallKPP
#   <version>     version of the products, e.g. v001
#
# Options (environment variables):
#
#   NJOBS=8      maximum number of subruns at a time
#   CHECK=120    seconds between checks
#   WORK=/data1/mu2e/carrarin/kpp_mc_<name>_<version>   logs and end markers (give it absolute)
#   STOPFILE=<WORK>/STOP                                when it exists, the script ends once no job is running
#
# Output:
#
#   /data1/mu2e/mc2/products/{kpp_mc_reco,kpp_mc_ntuple,kpp_mc_events}/<name>/<version>/
#                                                  the files above + MANIFEST.json
#   <WORK>/run.log                                 log of the whole production
#   <WORK>/rec_<subrun>.log, nts_<subrun>.log, ev_<subrun>.log   logs of the three steps
#   <WORK>/<subrun>.done / <subrun>.failed         end markers (a failed subrun is not retried)
#

MUSE_AREA=$(readlink -f $1)
SIM=$(readlink -f $2)
NAME=$3
VER=$4

if ! { [ -d "$MUSE_AREA/Offline/MC2Reco" ] && [ -d "$SIM" ] && [ -n "$VER" ]; }; then
  echo "usage: $0 <muse area> <sim dir> <name> <version>"
  exit 1
fi

NJOBS=${NJOBS:-8}
CHECK=${CHECK:-120}
WORK=${WORK:-/data1/mu2e/carrarin/kpp_mc_${NAME}_$VER}
STOPFILE=${STOPFILE:-$WORK/STOP}
PRODUCTS=/data1/mu2e/mc2/products
IMAGE=/cvmfs/singularity.opensciencegrid.org/fermilab/fnal-wn-el9:latest

# Commit of the fork and number of modified tracked files, for the log and the manifests
COMMIT=$(git -C $MUSE_AREA/Offline rev-parse --short HEAD)
DIRTY=$(git -C $MUSE_AREA/Offline status --porcelain --untracked-files=no | wc -l)

# Output directories of the three products
RECO_DIR=$PRODUCTS/kpp_mc_reco/$NAME/$VER
NTUPLE_DIR=$PRODUCTS/kpp_mc_ntuple/$NAME/$VER
EVENTS_DIR=$PRODUCTS/kpp_mc_events/$NAME/$VER

mkdir -p $WORK $RECO_DIR $NTUPLE_DIR $EVENTS_DIR
cd $WORK


# ---- helpers ---------------------------------------------------------------------------------------------------------

log() {
  echo "$(date '+%F %T') $*" | tee -a $WORK/run.log
}

# One MANIFEST.json per product, listing the subruns done so far
manifest() {
  local product_dir
  for product_dir in $RECO_DIR $NTUPLE_DIR $EVENTS_DIR; do
    python3 - "$product_dir" "$PRODUCTS" "$RECO_DIR" "$NTUPLE_DIR" "$EVENTS_DIR" "$WORK" "$NAME" "$VER" "$COMMIT" \
      "$DIRTY" "$MUSE_AREA" "$SIM" "$(id -un)" <<'EOF'
import datetime
import glob
import json
import os
import sys

(product_dir, products, reco_dir, ntuple_dir, events_dir, work, name, version, commit,
 dirty, muse_area, sim, author) = sys.argv[1:14]

# subruns done so far: the end markers <subrun>.done in the work directory
done = sorted(os.path.basename(f)[:-5] for f in glob.glob(work + "/*.done"))

descriptions = {
    reco_dir: "digitization and reconstruction (tracker and CRV as data, MC2 calorimeter response)",
    ntuple_dir: "EventNtuple",
    events_dir: "tracks with in-time calorimeter crystals and CRV coincidences (kpp_dump_events.py)",
}
manifest = dict(
    product=product_dir.replace(products + "/", ""),
    description=f"KPP hall simulation {name}: " + descriptions[product_dir],
    fork="github.com/Carrarin/Offline mc2-calo " + commit + (" (local changes)" if int(dirty) else ""),
    environment=f"muse setup {muse_area} (backing DataJob 1.11.0) in fnal-wn-el9",
    input=f"{sim} (dts.mc2.{name}.<subrun>.art, see its MANIFEST.json)",
    fcl=["Offline/MC2Reco/fcl/sim/mc2_kpp_digi_reco.fcl", "Offline/MC2Reco/fcl/kpp/kpp_ntuple.fcl"],
    geometry="Offline/MC2Reco/data/geom/geom_extracted_mc2kpp.txt",
    command=f"Offline/MC2Reco/scripts/run_kpp_mc.sh {muse_area} {sim} {name} {version}",
    subruns=done,
    n_files=len(done),
    updated=datetime.datetime.now().astimezone().isoformat(timespec="seconds"),
    author=author,
)
json.dump(manifest, open(os.path.join(product_dir, "MANIFEST.json"), "w"), indent=1)
EOF
  done
}

# Digitization + reconstruction, ntuple and event dump of one subrun, in the container at idle priority:
# job <subrun (3 digits)>
job() {
  local subrun=$1
  local dts_file=$SIM/dts.mc2.$NAME.$1.art
  local reco_file=$RECO_DIR/rec_$subrun.art
  local ntuple_file=$NTUPLE_DIR/nts_$subrun.root
  local events_file=$EVENTS_DIR/ev_$subrun.json

  # the three steps; the exit status says which one failed (1 reco, 2 ntuple, 3 event dump)
  local commands="source /cvmfs/mu2e.opensciencegrid.org/setupmu2e-art.sh >/dev/null 2>&1
    cd $MUSE_AREA
    muse setup > /dev/null 2>&1
    export MU2E_SEARCH_PATH=\$MU2E_SEARCH_PATH:$MUSE_AREA:/data1/mu2e/mc2
    mu2e -c Offline/MC2Reco/fcl/sim/mc2_kpp_digi_reco.fcl -s $dts_file -o $reco_file > $WORK/rec_$subrun.log 2>&1 \
|| exit 1
    mu2e -c Offline/MC2Reco/fcl/kpp/kpp_ntuple.fcl -s $reco_file -T $ntuple_file > $WORK/nts_$subrun.log 2>&1 || exit 2
    python3 Offline/MC2Reco/scripts/kpp_dump_events.py $events_file $ntuple_file > $WORK/ev_$subrun.log 2>&1 || exit 3"
  chrt --idle 0 nice -n 19 apptainer exec -B /cvmfs,/data1 $IMAGE bash -c "$commands" < /dev/null
  local status=$?

  if [ $status = 0 ]; then
    touch $WORK/$subrun.done
    log "subrun $subrun done ($(tail -1 $WORK/ev_$subrun.log))"
    manifest
  else
    touch $WORK/$subrun.failed
    log "subrun $subrun FAILED at step $status, see $WORK"
  fi
}


# ---- main loop: every CHECK seconds, until STOPFILE exists and no job is running ------------------------------------

log "start: $SIM ($NAME), version $VER, fork $COMMIT (local changes: $DIRTY), NJOBS $NJOBS"
declare -A PID   # subrun -> pid of the running job
while true; do

  # forget the jobs that ended (their .done/.failed marker is written by the job itself)
  for subrun in "${!PID[@]}"; do
    if ! kill -0 ${PID[$subrun]} 2>/dev/null; then
      unset PID[$subrun]
    fi
  done

  # start the finished simulation files, in name order, while there is room
  for dts_file in $(ls $SIM/dts.mc2.$NAME.*.art 2>/dev/null); do
    subrun=$(basename $dts_file .art)
    subrun=${subrun##*.}
    subrun_number=$((10#$subrun))   # the simulation job log has no leading zeros
    if [ -e $WORK/$subrun.done ] || [ -e $WORK/$subrun.failed ] || [ -n "${PID[$subrun]}" ]; then
      continue   # done, failed or running
    fi
    if [ ${#PID[@]} -ge $NJOBS ]; then
      break      # no room
    fi
    if ! grep -q "Art has completed and will exit with status 0" $SIM/job_$subrun_number.log 2>/dev/null; then
      continue   # simulation not finished (or failed)
    fi
    if ! [ $(( $(date +%s) - $(stat -c %Y $dts_file) )) -gt 120 ]; then
      continue   # changed less than 2 minutes ago
    fi
    job $subrun &
    PID[$subrun]=$!
    log "subrun $subrun started"
  done

  if [ -e $STOPFILE ] && [ ${#PID[@]} = 0 ]; then
    break
  fi
  sleep $CHECK
done

log "stopped"
