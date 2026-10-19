#!/bin/bash
# Mu2eG4 cosmic production in the Mu2e hall for the calorimeter: N jobs of NEV events, job i -> subrun i, seed 7100+i.
# usage: Offline/MC2Reco/scripts/run_sim_hall.sh <out dir> <n jobs> <events per job>     (from the mc2reco area)
OUT=$1; NJ=$2; NEV=$3; AREA=$(pwd); mkdir -p $OUT
source /cvmfs/mu2e.opensciencegrid.org/setupmu2e-art.sh >/dev/null 2>&1
cd $OUT && muse setup SimJob >/dev/null 2>&1
export MU2E_SEARCH_PATH=$AREA:$MU2E_SEARCH_PATH FHICL_FILE_PATH=$AREA:$FHICL_FILE_PATH
for i in $(seq 0 $((NJ-1))); do
  printf '#include "Offline/MC2Reco/fcl/sim/cosmic_cry_hall_calo.fcl"\nsource.firstSubRun : %d\nservices.SeedService.baseSeed : %d\noutputs.PrimaryOutput.fileName : "dts.mc2.CosmicCRYHallCalo.%03d.art"\n' $i $((7100+i)) $i > job_$i.fcl
  nohup mu2e -c job_$i.fcl -n $NEV > job_$i.log 2>&1 &
done
wait
echo done
