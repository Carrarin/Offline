#!/bin/bash
# The Mu2eG4 hall production (as run_sim_hall.sh) on a machine without an EL9 Mu2e environment (e.g. gm2ita, AlmaLinux 10): each job
# in the Fermilab EL9 container from cvmfs. Job i -> subrun FIRST+i, seed 7100+FIRST+i (no overlap with other productions if FIRST differs).
# usage: run_sim_hall_container.sh <area with Offline/MC2Reco/fcl/sim> <out dir> <n jobs> <events per job> <first subrun>
AREA=$1; OUT=$2; NJ=$3; NEV=$4; FIRST=$5
IMG=/cvmfs/singularity.opensciencegrid.org/fermilab/fnal-wn-el9:latest
mkdir -p $OUT; cd $OUT
for i in $(seq 0 $((NJ-1))); do
  s=$((FIRST+i))
  printf '#include "Offline/MC2Reco/fcl/sim/cosmic_cry_hall_calo.fcl"\nsource.firstSubRun : %d\nservices.SeedService.baseSeed : %d\noutputs.PrimaryOutput.fileName : "dts.mc2.CosmicCRYHallCalo.%03d.art"\n' $s $((7100+s)) $s > job_$s.fcl
  setsid nohup apptainer exec -B /cvmfs,$AREA,$OUT $IMG bash -c "source /cvmfs/mu2e.opensciencegrid.org/setupmu2e-art.sh >/dev/null 2>&1; muse setup SimJob >/dev/null 2>&1; export MU2E_SEARCH_PATH=$AREA:\$MU2E_SEARCH_PATH FHICL_FILE_PATH=$AREA:\$FHICL_FILE_PATH; cd $OUT; mu2e -c job_$s.fcl -n $NEV" > job_$s.log 2>&1 < /dev/null &
done
echo "launched $NJ jobs in $OUT"
