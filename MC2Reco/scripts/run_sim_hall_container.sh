#!/bin/bash
#
# run_sim_hall_container.sh: the Mu2eG4 hall production of run_sim_hall.sh on a machine without an EL9 Mu2e
# environment (e.g. gm2ita, AlmaLinux 10): each job runs in the Fermilab EL9 container from cvmfs.
# Job i writes subrun FIRST + i with seed 7100 + FIRST + i: productions with different FIRST do not overlap.
# The jobs are launched detached (setsid nohup) and the script returns at once, without waiting for them.
#
# Usage:
#
#   run_sim_hall_container.sh <area> <out dir> <n jobs> <events per job> <first subrun>
#
#   e.g. Offline/MC2Reco/scripts/run_sim_hall_container.sh /data1/mu2e/carrarin/mc2offline \
#          /data1/mu2e/carrarin/sim_hall 16 20000 100
#
# Arguments:
#
#   <area>             directory with Offline/MC2Reco/fcl/sim (put first in the search paths in the container)
#   <out dir>          output directory, created if missing; the jobs run in it
#   <n jobs>           number of jobs
#   <events per job>   events of each job (mu2e -n)
#   <first subrun>     subrun of the first job
#
# Output in <out dir>, for every subrun s:
#
#   job_<s>.fcl                            job fcl (subrun, seed, output file name)
#   job_<s>.log                            job log
#   dts.mc2.CosmicCRYHallCalo.<sss>.art    simulated events (sss = s with 3 digits)
#
# Where: any machine with apptainer and cvmfs (gm2ita), outside the container. Give <area> and <out dir> as absolute
# paths: they are bound into the container and used there. The jobs run at normal priority (no nice).
#

AREA=$1
OUT=$2
N_JOBS=$3
EVENTS_PER_JOB=$4
FIRST_SUBRUN=$5
IMAGE=/cvmfs/singularity.opensciencegrid.org/fermilab/fnal-wn-el9:latest

mkdir -p $OUT
cd $OUT


# Job fcl of one subrun: write_job_fcl <subrun>
write_job_fcl() {
  local subrun=$1
  {
    printf '#include "Offline/MC2Reco/fcl/sim/cosmic_cry_hall_calo.fcl"\n'
    printf 'source.firstSubRun : %d\n' $subrun
    printf 'services.SeedService.baseSeed : %d\n' $((7100 + subrun))
    printf 'outputs.PrimaryOutput.fileName : "dts.mc2.CosmicCRYHallCalo.%03d.art"\n' $subrun
  } > job_$subrun.fcl
}

# One detached job in the container: Mu2e environment, SimJob, the fcl of <area> first in the search paths
launch_job() {
  local subrun=$1
  local commands="source /cvmfs/mu2e.opensciencegrid.org/setupmu2e-art.sh >/dev/null 2>&1"
  commands+="; muse setup SimJob >/dev/null 2>&1"
  commands+="; export MU2E_SEARCH_PATH=$AREA:\$MU2E_SEARCH_PATH FHICL_FILE_PATH=$AREA:\$FHICL_FILE_PATH"
  commands+="; cd $OUT"
  commands+="; mu2e -c job_$subrun.fcl -n $EVENTS_PER_JOB"
  setsid nohup apptainer exec -B /cvmfs,$AREA,$OUT $IMAGE bash -c "$commands" > job_$subrun.log 2>&1 < /dev/null &
}


for i in $(seq 0 $((N_JOBS - 1))); do
  subrun=$((FIRST_SUBRUN + i))
  write_job_fcl $subrun
  launch_job $subrun
done

echo "launched $N_JOBS jobs in $OUT"
