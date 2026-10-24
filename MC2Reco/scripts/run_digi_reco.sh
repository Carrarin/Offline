#!/bin/bash
#
# run_digi_reco.sh: digitization and reconstruction of a hall simulation production, one job per dts file.
# Every dts.*.art file of the production goes through Offline/MC2Reco/fcl/sim/<fcl> (by default mc2_digi_reco.fcl);
# the job of file number n uses the seed 8100 + n. All the jobs run at the same time, in the background.
#
# Usage (from the muse area, the directory with Offline/; the script does "muse setup" itself):
#
#   Offline/MC2Reco/scripts/run_digi_reco.sh <production dir> <out dir> [fcl]
#
#   e.g. Offline/MC2Reco/scripts/run_digi_reco.sh sim/CosmicCRYHallCalo/v001 digi/CosmicCRYHallCalo
#
# Arguments:
#
#   <production dir>   directory with the dts.<...>.<n>.art files of the simulation (n = file number, last field)
#   <out dir>          output directory, created if missing
#   [fcl]              fcl of Offline/MC2Reco/fcl/sim to include (default mc2_digi_reco.fcl)
#
# Output in <out dir>, for every file number n:
#
#   digi_<n>.fcl     job fcl (the base fcl + seed + summary file of the stopped muons)
#   digi_<n>.root    TFileService output (dump/tracks, stopped/..., used by sim_vs_data.py)
#   digi_<n>.log     job log
#   stopped_<n>.txt  summary of MC2StoppedMuons
#
# Where: from the muse area, in a shell with the EL9 Mu2e environment available (inside the fnal-wn-el9 container on
# gm2ita or mu2ecalo01: the script does not start the container). No nice: the jobs run at normal priority.
# Give <out dir> and <production dir> as absolute paths: the jobs run inside <out dir>.
#

PRODUCTION_DIR=$1
OUT=$2
FCL=${3:-mc2_digi_reco.fcl}

mkdir -p $OUT

# Mu2e environment and the muse area of the current directory
source /cvmfs/mu2e.opensciencegrid.org/setupmu2e-art.sh >/dev/null 2>&1
muse setup >/dev/null 2>&1


# Job fcl of one file: write_job_fcl <file number>
write_job_fcl() {
  local file_number=$1
  {
    printf '#include "Offline/MC2Reco/fcl/sim/%s"\n' $FCL
    printf 'services.SeedService.baseSeed : %d\n' $((8100 + 10#$file_number))
    printf 'physics.analyzers.stopped.summaryFile : "%s/stopped_%s.txt"\n' $OUT $file_number
  } > $OUT/digi_$file_number.fcl
}


# One background job per dts file; the file number is the last dot-separated field of the name
for dts_file in $PRODUCTION_DIR/dts.*.art; do
  file_number=$(basename $dts_file .art | awk -F. '{print $NF}')
  write_job_fcl $file_number
  ( cd $OUT && mu2e -c digi_$file_number.fcl -s $dts_file -T digi_$file_number.root > digi_$file_number.log 2>&1 ) &
done

wait
echo done
