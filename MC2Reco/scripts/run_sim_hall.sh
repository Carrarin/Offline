#!/bin/bash
#
# run_sim_hall.sh: Mu2eG4 production of cosmic rays (CRY) in the Mu2e hall for the calorimeter.
# Runs <n jobs> jobs of <events per job> events of Offline/MC2Reco/fcl/sim/cosmic_cry_hall_calo.fcl, all at the same
# time; job i writes subrun i with seed 7100 + i.
#
# Usage (from the muse area, the directory with Offline/; the script does "muse setup SimJob" itself):
#
#   Offline/MC2Reco/scripts/run_sim_hall.sh <out dir> <n jobs> <events per job>
#
#   e.g. Offline/MC2Reco/scripts/run_sim_hall.sh /data1/mu2e/carrarin/sim_hall 16 20000
#
# Arguments:
#
#   <out dir>          output directory, created if missing; the jobs run in it
#   <n jobs>           number of jobs (subruns 0 .. n-1)
#   <events per job>   events of each job (mu2e -n)
#
# Output in <out dir>, for every job i:
#
#   job_<i>.fcl                            job fcl (subrun, seed, output file name)
#   job_<i>.log                            job log
#   dts.mc2.CosmicCRYHallCalo.<iii>.art    simulated events (iii = i with 3 digits)
#
# Where: a machine with the EL9 Mu2e environment (in the fnal-wn-el9 container: the script does not start it; for a
# machine without it, e.g. gm2ita outside the container, see run_sim_hall_container.sh). The jobs run with nohup at
# normal priority; the script waits for all of them.
#

OUT=$1
N_JOBS=$2
EVENTS_PER_JOB=$3
AREA=$(pwd)

mkdir -p $OUT

# Mu2e environment; SimJob is set up from <out dir>, and the fcl of the area are put first in the search paths
source /cvmfs/mu2e.opensciencegrid.org/setupmu2e-art.sh >/dev/null 2>&1
cd $OUT && muse setup SimJob >/dev/null 2>&1
export MU2E_SEARCH_PATH=$AREA:$MU2E_SEARCH_PATH FHICL_FILE_PATH=$AREA:$FHICL_FILE_PATH


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


# Start all the jobs in the background, then wait for them
for i in $(seq 0 $((N_JOBS - 1))); do
  write_job_fcl $i
  nohup mu2e -c job_$i.fcl -n $EVENTS_PER_JOB > job_$i.log 2>&1 &
done

wait
echo done
