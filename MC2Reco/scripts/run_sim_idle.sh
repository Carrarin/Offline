#!/bin/bash
#
# run_sim_idle.sh: a Mu2eG4 hall production on the idle cycles of a shared machine.
#
# Every job runs in the EL9 container under SCHED_IDLE (chrt --idle 0: it gets a CPU only when nobody else wants one)
# and nice 19, in its own session (setsid), so that it can be stopped with its whole process group.
# At most MAXJOBS jobs at a time. A new job starts only if MemAvailable stays above MEMFREE GB after its JOBMEM GB.
# If MemAvailable drops below MEMMIN GB, or the other users take more than OTHERCPU cores, the job with the highest
# subrun (normally the newest) is stopped and its subrun goes back to the queue, so that people who log in get the
# machine back. A job that ends with art status 0 is copied with rsync to DEST (host:dir, optional) and its art file is
# removed here once the copy has the same size; a copy that fails is retried at the next check. Failed jobs are logged,
# not rerun.
# Subrun s -> seed 7100 + s, as in the other productions: give each machine its own subrun range.
#
# Only Offline/MC2Reco/fcl/sim and Offline/MC2Reco/data/geom of the area are copied to <work dir>/area and put in the
# search paths: a whole Offline there would shadow the fcl of SimJob (an older Offline) with ones its code does not
# accept. The copy also freezes the configuration of the production (it is not updated if it already exists).
#
# Usage (outside the container; meant to be started detached, e.g. with setsid nohup ... &):
#
#   run_sim_idle.sh <area> <work dir> <first subrun> <last subrun> <events per job> [DEST]
#
#   e.g. Offline/MC2Reco/scripts/run_sim_idle.sh /data1/mu2e/carrarin/mc2offline /data1/mu2e/carrarin/sim_idle \
#          200 299 20000 gm2ita:/data1/mu2e/mc2/sim/CosmicCRYHallCalo/v002
#
# Arguments:
#
#   <area>             directory with Offline/MC2Reco (only its fcl/sim and data/geom are used)
#   <work dir>         where the jobs run and write their files, created if missing
#   <first subrun>, <last subrun>   subruns to produce (one job each); subruns already .done or .failed are skipped
#   <events per job>   events of each job (mu2e -n)
#   [DEST]             host:dir where the finished art files are copied (rsync over ssh, BatchMode: needs a key)
#
# Options (environment variables):
#
#   FCL=cosmic_cry_hall_calo.fcl   fcl of Offline/MC2Reco/fcl/sim
#   NAME=CosmicCRYHallCalo         name in the output files dts.mc2.<NAME>.<sss>.art
#   MAXJOBS=16                     maximum number of jobs at a time
#   MEMFREE=8                      GB of MemAvailable that must remain after starting a job
#   MEMMIN=4                       GB of MemAvailable below which a job is stopped
#   JOBMEM=2.1                     GB of memory of one job
#   OTHERCPU=2                     cores used by other users above which jobs are stopped
#   CHECK=60                       seconds between checks
#
# Output in <work dir>:
#
#   idle.log                                  log of the production (starts, stops, ends, copies)
#   area/Offline/MC2Reco/{fcl/sim,data/geom}  frozen copy of the configuration
#   job_<s>.fcl, job_<s>.log                  job fcl and log of subrun s
#   job_<s>.done / job_<s>.failed             end markers (art status 0 / anything else)
#   dts.mc2.<NAME>.<sss>.art                  simulated events (removed after a verified copy to DEST)
#
# Where: any machine with apptainer and cvmfs (gm2ita, mu2ecalo01), outside the container.
#

AREA=$(readlink -f $1)
WORK=$2
FIRST_SUBRUN=$3
LAST_SUBRUN=$4
EVENTS_PER_JOB=$5
DEST=$6

if ! { [ -d "$AREA/Offline/MC2Reco/fcl/sim" ] && [ -n "$WORK" ] && [ -n "$LAST_SUBRUN" ] && [ -n "$EVENTS_PER_JOB" ]; }
then
  echo "usage: $0 <area with Offline/MC2Reco> <work dir> <first subrun> <last subrun> <events per job> [DEST]"
  exit 1
fi

FCL=${FCL:-cosmic_cry_hall_calo.fcl}
NAME=${NAME:-CosmicCRYHallCalo}
MAXJOBS=${MAXJOBS:-16}
MEMFREE=${MEMFREE:-8}
MEMMIN=${MEMMIN:-4}
JOBMEM=${JOBMEM:-2.1}
OTHERCPU=${OTHERCPU:-2}
CHECK=${CHECK:-60}
IMAGE=/cvmfs/singularity.opensciencegrid.org/fermilab/fnal-wn-el9:latest
SSH="ssh -x -o BatchMode=yes"

mkdir -p $WORK
cd $WORK
WORK=$(pwd)


# ---- frozen copy of the configuration --------------------------------------------------------------------------------
# Each directory is copied only if it is not there yet and exists in the area.

for config_dir in fcl/sim data/geom; do
  mkdir -p area/Offline/MC2Reco/$(dirname $config_dir)
  if [ ! -d area/Offline/MC2Reco/$config_dir ] && [ -d $AREA/Offline/MC2Reco/$config_dir ]; then
    cp -r $AREA/Offline/MC2Reco/$config_dir area/Offline/MC2Reco/$config_dir
  fi
done
AREA=$WORK/area


# ---- helpers ---------------------------------------------------------------------------------------------------------

log() {
  echo "$(date '+%F %T') $*" | tee -a $WORK/idle.log
}

# MemAvailable of the machine, in GB with one decimal
memavail() {
  awk '/MemAvailable/ {printf "%.1f", $2 / 1e6}' /proc/meminfo
}

# CPU used by the other users (not me, not root), in cores with one decimal
othercpu() {
  ps -eo user=,pcpu= | awk -v me=$(id -un) '$1 != me && $1 != "root" {s += $2} END {printf "%.1f", s / 100}'
}

# Floating-point comparison (bc): is_true "<expression>" succeeds if the expression is true
is_true() {
  [ $(echo "$1" | bc) = 1 ]
}

# Name of the art file of a subrun (3 digits)
art_file_name() {
  echo dts.mc2.$NAME.$(printf %03d $1).art
}

# Start the job of a subrun, detached and at idle priority: launch <subrun>
launch() {
  local subrun=$1
  local art_file=$(art_file_name $subrun)
  rm -f $art_file
  {
    printf '#include "Offline/MC2Reco/fcl/sim/%s"\n' $FCL
    printf 'source.firstSubRun : %d\n' $subrun
    printf 'services.SeedService.baseSeed : %d\n' $((7100 + subrun))
    printf 'outputs.PrimaryOutput.fileName : "%s"\n' $art_file
  } > job_$subrun.fcl
  local commands="source /cvmfs/mu2e.opensciencegrid.org/setupmu2e-art.sh >/dev/null 2>&1"
  commands+="; muse setup SimJob >/dev/null 2>&1"
  commands+="; export MU2E_SEARCH_PATH=$AREA:\$MU2E_SEARCH_PATH FHICL_FILE_PATH=$AREA:\$FHICL_FILE_PATH"
  commands+="; cd $WORK"
  commands+="; mu2e -c job_$subrun.fcl -n $EVENTS_PER_JOB"
  setsid chrt --idle 0 nice -n 19 apptainer exec -B /cvmfs,$WORK $IMAGE bash -c "$commands" \
    > job_$subrun.log 2>&1 < /dev/null &
  PID[$subrun]=$!
  log "subrun $subrun started (pid ${PID[$subrun]}, MemAvailable $(memavail) GB)"
}

# Stop the running job with the highest subrun (its whole process group) and put the subrun back at the head of the
# queue: stop_newest <reason>
stop_newest() {
  local subrun=$(for running in "${!PID[@]}"; do echo $running; done | sort -n | tail -1)
  kill -- -${PID[$subrun]} 2>/dev/null
  unset PID[$subrun]
  QUEUE=($subrun "${QUEUE[@]}")
  log "$1: subrun $subrun stopped, back in the queue"
}

# Copy a finished subrun to DEST and remove the art file here when the copy has the same size: ship <subrun>
ship() {
  local subrun=$1
  local art_file=$(art_file_name $subrun)
  if [ ! -e $art_file ]; then
    return 0   # already copied
  fi
  if [ -z "$DEST" ]; then
    return 0   # no destination: the file stays here
  fi
  if ! rsync -a -e "$SSH" job_$subrun.fcl job_$subrun.log $art_file $DEST/ 2>> $WORK/idle.log; then
    log "subrun $subrun: copy to $DEST failed, retry later"
    return 1
  fi
  local size_here=$(stat -c %s $art_file)
  local size_there=$($SSH ${DEST%%:*} "stat -c %s ${DEST#*:}/$art_file" 2>/dev/null)
  if [ "$size_here" = "$size_there" ]; then
    rm -f $art_file
    log "subrun $subrun copied to $DEST ($size_here bytes)"
  else
    log "subrun $subrun: size here $size_here, there ${size_there:-none}, kept here"
  fi
}


# ---- queue of the subruns to do --------------------------------------------------------------------------------------

declare -A PID   # subrun -> pid of the running job (process group leader)
QUEUE=()
for subrun in $(seq $FIRST_SUBRUN $LAST_SUBRUN); do
  if [ ! -e job_$subrun.done ] && [ ! -e job_$subrun.failed ]; then
    QUEUE+=($subrun)
  fi
done
log "start: $FCL ($NAME), subruns $FIRST_SUBRUN-$LAST_SUBRUN, ${#QUEUE[@]} to do, $EVENTS_PER_JOB events each," \
  "MAXJOBS $MAXJOBS, MEMFREE $MEMFREE GB, OTHERCPU $OTHERCPU, dest ${DEST:-none}"


# ---- main loop: every CHECK seconds, until the queue is empty and no job is running ----------------------------------

while [ ${#QUEUE[@]} -gt 0 ] || [ ${#PID[@]} -gt 0 ]; do

  # 1. jobs that ended: done (art status 0) or failed
  for subrun in "${!PID[@]}"; do
    if kill -0 ${PID[$subrun]} 2>/dev/null; then
      continue
    fi
    unset PID[$subrun]
    if grep -q "Art has completed and will exit with status 0" job_$subrun.log; then
      touch job_$subrun.done
      log "subrun $subrun done"
    else
      touch job_$subrun.failed
      log "subrun $subrun FAILED, see job_$subrun.log"
    fi
  done

  # 2. copy the finished subruns to DEST (those already copied have no art file left and are skipped)
  for done_file in job_*.done; do
    if [ -e "$done_file" ]; then
      subrun=${done_file#job_}
      ship ${subrun%.done}
    fi
  done

  # 3. too little memory: stop one job
  if [ ${#PID[@]} -gt 0 ] && is_true "$(memavail) < $MEMMIN"; then
    stop_newest "MemAvailable $(memavail) GB < $MEMMIN"
  fi

  # 4. other users active: limit = free cores (all minus theirs minus 2), between 1 and MAXJOBS; stop jobs above it
  others=$(othercpu)
  limit=$MAXJOBS
  if is_true "$others > $OTHERCPU"; then
    limit=$(( $(nproc) - ${others%.*} - 2 ))
    if [ $limit -gt $MAXJOBS ]; then
      limit=$MAXJOBS
    fi
    if [ $limit -lt 1 ]; then
      limit=1
    fi
    while [ ${#PID[@]} -gt $limit ]; do
      stop_newest "other users use $others cores"
    done
  fi

  # 5. start new jobs while there is room and memory
  while [ ${#QUEUE[@]} -gt 0 ] && [ ${#PID[@]} -lt $limit ] && is_true "$(memavail) - $JOBMEM > $MEMFREE"; do
    launch ${QUEUE[0]}
    QUEUE=("${QUEUE[@]:1}")
    sleep 20   # let the job allocate its memory before the next check
  done

  sleep $CHECK
done

log "all subruns $FIRST_SUBRUN-$LAST_SUBRUN processed"
