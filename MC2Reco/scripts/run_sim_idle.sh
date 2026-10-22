#!/bin/bash
# run_sim_idle.sh -- a Mu2eG4 hall production on the idle cycles of a shared machine.
# Every job runs in the EL9 container under SCHED_IDLE (chrt --idle 0: it gets a CPU only when nobody else wants one) and nice 19.
# At most MAXJOBS jobs; a new one starts only if MemAvailable stays above MEMFREE GB after its JOBMEM GB; if MemAvailable drops below
# MEMMIN GB, or the other users take more than OTHERCPU cores, the newest job is stopped and its subrun goes back to the queue, so
# that people who log in get the machine back. A job that ends with art status 0 is copied with rsync to DEST (host:dir, optional)
# and its art file removed here; a copy that fails is retried at the next check. Failed jobs are logged, not rerun.
# Job subrun s -> seed 7100+s, as the other productions: give each machine its own subrun range.
# Only Offline/MC2Reco/fcl/sim and Offline/MC2Reco/data/geom of the area are copied to <work dir>/area and put in the search paths:
# a whole Offline there would shadow the fcl of SimJob (an older Offline) with ones its code does not accept. The copy also freezes
# the production configuration.
# usage: run_sim_idle.sh <area with Offline/MC2Reco> <work dir> <first subrun> <last subrun> <events per job> [DEST]
#   env: FCL (cosmic_cry_hall_calo.fcl), NAME (CosmicCRYHallCalo), MAXJOBS (16), MEMFREE (8), MEMMIN (4), JOBMEM (2.1),
#        OTHERCPU (2 cores used by other users), CHECK (60 s between checks)
AREA=$(readlink -f $1); WORK=$2; FIRST=$3; LAST=$4; NEV=$5; DEST=$6
[ -d "$AREA/Offline/MC2Reco/fcl/sim" ] && [ -n "$WORK" ] && [ -n "$LAST" ] && [ -n "$NEV" ] || {
  echo "usage: $0 <area with Offline/MC2Reco> <work dir> <first subrun> <last subrun> <events per job> [DEST]"; exit 1; }
FCL=${FCL:-cosmic_cry_hall_calo.fcl}; NAME=${NAME:-CosmicCRYHallCalo}
MAXJOBS=${MAXJOBS:-16}; MEMFREE=${MEMFREE:-8}; MEMMIN=${MEMMIN:-4}; JOBMEM=${JOBMEM:-2.1}; OTHERCPU=${OTHERCPU:-2}; CHECK=${CHECK:-60}
IMG=/cvmfs/singularity.opensciencegrid.org/fermilab/fnal-wn-el9:latest
mkdir -p $WORK; cd $WORK; WORK=$(pwd)
for d in fcl/sim data/geom; do
  mkdir -p area/Offline/MC2Reco/$(dirname $d); [ -d area/Offline/MC2Reco/$d ] || [ ! -d $AREA/Offline/MC2Reco/$d ] || cp -r $AREA/Offline/MC2Reco/$d area/Offline/MC2Reco/$d
done
AREA=$WORK/area
log() { echo "$(date '+%F %T') $*" | tee -a $WORK/idle.log; }
memavail() { awk '/MemAvailable/ {printf "%.1f", $2 / 1e6}' /proc/meminfo; }   # GB
othercpu() { ps -eo user=,pcpu= | awk -v me=$(id -un) '$1 != me && $1 != "root" {s += $2} END {printf "%.1f", s / 100}'; }   # cores
declare -A PID
QUEUE=(); for s in $(seq $FIRST $LAST); do [ -e job_$s.done ] || [ -e job_$s.failed ] || QUEUE+=($s); done
log "start: $FCL ($NAME), subruns $FIRST-$LAST, ${#QUEUE[@]} to do, $NEV events each, MAXJOBS $MAXJOBS, MEMFREE $MEMFREE GB, OTHERCPU $OTHERCPU, dest ${DEST:-none}"

launch() {
  local s=$1 f=dts.mc2.$NAME.$(printf %03d $1).art
  rm -f $f
  printf '#include "Offline/MC2Reco/fcl/sim/%s"\nsource.firstSubRun : %d\nservices.SeedService.baseSeed : %d\noutputs.PrimaryOutput.fileName : "%s"\n' $FCL $s $((7100+s)) $f > job_$s.fcl
  setsid chrt --idle 0 nice -n 19 apptainer exec -B /cvmfs,$WORK $IMG bash -c "source /cvmfs/mu2e.opensciencegrid.org/setupmu2e-art.sh >/dev/null 2>&1; muse setup SimJob >/dev/null 2>&1; export MU2E_SEARCH_PATH=$AREA:\$MU2E_SEARCH_PATH FHICL_FILE_PATH=$AREA:\$FHICL_FILE_PATH; cd $WORK; mu2e -c job_$s.fcl -n $NEV" > job_$s.log 2>&1 < /dev/null &
  PID[$s]=$!; log "subrun $s started (pid ${PID[$s]}, MemAvailable $(memavail) GB)"
}
stop_newest() {   # $1 = reason
  local s=$(for k in "${!PID[@]}"; do echo $k; done | sort -n | tail -1)
  kill -- -${PID[$s]} 2>/dev/null; unset PID[$s]; QUEUE=($s "${QUEUE[@]}"); log "$1: subrun $s stopped, back in the queue"
}
ship() {   # copy a finished subrun to DEST, remove the art file here when the copy has the same size
  local s=$1 f=dts.mc2.$NAME.$(printf %03d $1).art
  [ -e $f ] || return 0; [ -z "$DEST" ] && return 0
  rsync -a -e "ssh -x -o BatchMode=yes" job_$s.fcl job_$s.log $f $DEST/ 2>> $WORK/idle.log || { log "subrun $s: copy to $DEST failed, retry later"; return 1; }
  local here=$(stat -c %s $f) there=$(ssh -x -o BatchMode=yes ${DEST%%:*} "stat -c %s ${DEST#*:}/$f" 2>/dev/null)
  [ "$here" = "$there" ] && { rm -f $f; log "subrun $s copied to $DEST ($here bytes)"; } || log "subrun $s: size here $here, there ${there:-none}, kept here"
}

while [ ${#QUEUE[@]} -gt 0 ] || [ ${#PID[@]} -gt 0 ]; do
  for s in "${!PID[@]}"; do
    kill -0 ${PID[$s]} 2>/dev/null && continue
    unset PID[$s]
    if grep -q "Art has completed and will exit with status 0" job_$s.log; then touch job_$s.done; log "subrun $s done"
    else touch job_$s.failed; log "subrun $s FAILED, see job_$s.log"; fi
  done
  for d in job_*.done; do [ -e "$d" ] && { s=${d#job_}; ship ${s%.done}; }; done
  if [ ${#PID[@]} -gt 0 ] && [ $(echo "$(memavail) < $MEMMIN" | bc) = 1 ]; then stop_newest "MemAvailable $(memavail) GB < $MEMMIN"; fi
  others=$(othercpu); limit=$MAXJOBS
  if [ $(echo "$others > $OTHERCPU" | bc) = 1 ]; then
    limit=$(( $(nproc) - ${others%.*} - 2 )); [ $limit -gt $MAXJOBS ] && limit=$MAXJOBS; [ $limit -lt 1 ] && limit=1
    while [ ${#PID[@]} -gt $limit ]; do stop_newest "other users use $others cores"; done
  fi
  while [ ${#QUEUE[@]} -gt 0 ] && [ ${#PID[@]} -lt $limit ] && [ $(echo "$(memavail) - $JOBMEM > $MEMFREE" | bc) = 1 ]; do
    launch ${QUEUE[0]}; QUEUE=("${QUEUE[@]:1}"); sleep 20   # let the job allocate its memory before the next check
  done
  sleep $CHECK
done
log "all subruns $FIRST-$LAST processed"
