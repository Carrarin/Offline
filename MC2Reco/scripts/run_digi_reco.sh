#!/bin/bash
# Digitization + reconstruction (Offline/MC2Reco/fcl/sim/mc2_digi_reco.fcl) of every dts file of a hall production, one job per file
# (seed 8100 + file number). usage (from the mc2reco area): Offline/MC2Reco/scripts/run_digi_reco.sh <production dir> <out dir> [fcl, default mc2_digi_reco.fcl]
IN=$1; OUT=$2; FCL=${3:-mc2_digi_reco.fcl}; AREA=$(pwd); mkdir -p $OUT
source /cvmfs/mu2e.opensciencegrid.org/setupmu2e-art.sh >/dev/null 2>&1; muse setup >/dev/null 2>&1
for f in $IN/dts.*.art; do
  n=$(basename $f .art | awk -F. '{print $NF}')
  printf '#include "Offline/MC2Reco/fcl/sim/%s"\nservices.SeedService.baseSeed : %d\nphysics.analyzers.stopped.summaryFile : "%s/stopped_%s.txt"\n' $FCL $((8100+10#$n)) $OUT $n > $OUT/digi_$n.fcl
  ( cd $OUT && mu2e -c digi_$n.fcl -s $f -T digi_$n.root > digi_$n.log 2>&1 ) &
done
wait; echo done
