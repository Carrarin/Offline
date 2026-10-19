#!/bin/bash
# run_machine.sh -- the calorimeter machine on one run: calibration passes until they converge, then the physics with the final tables.
#
# usage: Offline/MC2Reco/scripts/run_machine.sh <run> <raw list> <out dir> [start CalTimeCalib table] [start CalEnergyCalib table]
#   (from the muse area: the directory with MC2Reco, EventNtuple and backing; muse setup done; the tables must be below it:
#   DbService looks for its text files in the search path, so they are written relative to it)
#
# Start tables:
#   time   : any table that exists (e.g. the last converged table of a nearby run, or the calorimeter group's); without one, pass 0
#            (no tables: all crystals of an event in one cluster, L and R of a crystal merged within 2 us) makes the first table
#   energy : by default one ADC2MeV for all the SiPMs (0.0625, crystals equalised); the machine measures the scale from the MIP
#   NEVER the laser for the times of cosmics (fibre lengths), NEVER the laser for photoelectrons (different wavelength, different PDE).
# Passes:
#   pass 1   : start tables -> calibration only (dump + calib) -> tables_1; WIDE (all crystals of an event in one cluster, L and R
#              merged within 2 us): a SiPM with a start offset wrong by more than the 6 ns L/R merging window would never join a track
#   pass k+1 : tables_k as conditions -> calibration again; converged when the time corrections have |median| < 0.1 ns and
#              MAD sigma < 0.3 ns and the energy ratio has |median - 1| < 0.005 (at most MAXPASS passes)
#   physics  : with the converged tables: tracks, stopped muons (Michel per disk, lifetime), captures (examples, injection)
# Output: <out dir>/pass<k>.{fcl,log,root}, CalTimeCalib_pass<k>.txt, CalEnergyCalib_pass<k>.txt, CalSimCrystals/CalBaselines/response,
#         final/ (the converged tables + README with the provenance), physics.root, stopped.txt, capture/ (analysis figures)
set -e
set -o pipefail   # a check piped into tee must keep its own exit status
RUN=$1; LIST=$(readlink -f $2); OUT=$3; TSTART=$4; ESTART=$5; MAXPASS=${MAXPASS:-4}
[ -z "$RUN" ] || [ ! -s "$LIST" ] || [ -z "$OUT" ] && { echo "usage: $0 <run> <raw list> <out dir> [start time table] [start energy table]"; exit 1; }
mkdir -p $OUT; OUT=$(readlink -f $OUT); BASE=Offline/MC2Reco/fcl/mc2_calo_machine.fcl
log() { echo "[run_machine $RUN] $*" | tee -a $OUT/machine.log; }
uniform() { python3 - "$1" "$2" "$3" <<'EOF'
import sys
f, name, val = sys.argv[1], sys.argv[2], float(sys.argv[3])
with open(f, "w") as o:
    o.write(f"TABLE {name}\n#roid,{'timeOffset' if name == 'CalTimeCalib' else 'ADC2MeV'}  start: one value for all SiPMs\n")
    for i in range(2720): o.write(f"{i},{val}\n")
EOF
}
[ -n "$ESTART" ] || { ESTART=$OUT/CalEnergyCalib_start_uniform.txt; uniform $ESTART CalEnergyCalib 0.0625; log "energy start: uniform 0.0625 MeV/ADC"; }
writefcl() {   # $1 pass name, $2 time table, $3 energy table, $4 analyzers, $5 extra lines
  local F=$OUT/$1.fcl; cp $BASE $F
  cat >> $F <<EOF
# $1 of run $RUN (run_machine.sh)
services.ProditionsService.calCalib.useDb : true
services.DbService.textFile : [ "$(realpath --relative-to=$PWD $3)", "$(realpath --relative-to=$PWD $2)" ]
physics.ana : [ $4 ]
physics.end_paths : [ ana ]
physics.analyzers.calib.timeTableFile : "$OUT/CalTimeCalib_$1.txt"
physics.analyzers.calib.energyTableFile : "$OUT/CalEnergyCalib_$1.txt"
physics.analyzers.calib.statsFile : "$OUT/calib_$1_stats.csv"
physics.analyzers.calib.simTableFile : "$OUT/CalSimCrystals_$1.txt"
physics.analyzers.calib.baselineTableFile : "$OUT/CalBaselines_$1.txt"
physics.analyzers.calib.responseFile : "$OUT/response_$1.csv"
physics.analyzers.stopped.summaryFile : "$OUT/stopped_$1.txt"
$5
EOF
}
runpass() { log "running $1"; mu2e -c $OUT/$1.fcl -S $LIST -T $OUT/$1.root > $OUT/$1.log 2>&1 || { log "FAILED $1, see $OUT/$1.log"; exit 1; }
  grep -E "MC2CaloCalibFromTracks:|MC2StoppedMuons|lifetime|Michel disk|MC2MuCapture:" $OUT/$1.log | tee -a $OUT/machine.log; }
converged() { python3 - "$OUT/calib_$1_stats.csv" <<'EOF'
import sys, numpy as np
s = np.genfromtxt(sys.argv[1], delimiter=",", names=True); t = s[s["timeFromTracks"] == 1]; e = s[s["energyFromTracks"] == 1]
mt = np.median(t["timeCorrection_ns"]); st = 1.4826 * np.median(np.abs(t["timeCorrection_ns"] - mt)); re = np.median(e["newADC2MeV"] / e["startADC2MeV"])
print(f"time correction median {mt:.3f} ns, MAD sigma {st:.3f} ns ({len(t)} SiPMs); energy ratio median {re:.4f} ({len(e)} SiPMs)")
sys.exit(0 if abs(mt) < 0.1 and st < 0.3 and abs(re - 1) < 0.005 else 1)
EOF
}
# ---- pass 0 if there is no time table ----
if [ -z "$TSTART" ]; then
  uniform $OUT/CalTimeCalib_start_zero.txt CalTimeCalib 0; TSTART=$OUT/CalTimeCalib_start_zero.txt
  writefcl pass0 $TSTART $ESTART "dump, calib" "physics.producers.MC2CaloCosmicTrackFinder.maxGap : 1e9
physics.producers.MC2CaloHitMaker.time4Merge : 2000."
  runpass pass0; TSTART=$OUT/CalTimeCalib_pass0.txt; ESTART=$OUT/CalEnergyCalib_pass0.txt
fi
log "time start: $TSTART ; energy start: $ESTART"
# ---- calibration passes ----
T=$TSTART; E=$ESTART; k=1
while [ $k -le $MAXPASS ]; do
  WIDE=""; [ $k -eq 1 ] && WIDE="physics.producers.MC2CaloCosmicTrackFinder.maxGap : 1e9
physics.producers.MC2CaloHitMaker.time4Merge : 2000."   # first pass always wide: SiPMs with a wrong start offset still join the tracks
  writefcl pass$k $T $E "dump, calib" "$WIDE"; runpass pass$k
  T=$OUT/CalTimeCalib_pass$k.txt; E=$OUT/CalEnergyCalib_pass$k.txt
  if [ $k -gt 1 ]; then
    if converged pass$k > $OUT/conv_$k.txt; then cat $OUT/conv_$k.txt | tee -a $OUT/machine.log; log "converged at pass $k"; break; fi
    cat $OUT/conv_$k.txt | tee -a $OUT/machine.log
  else converged pass1 > $OUT/conv_1.txt || true; cat $OUT/conv_1.txt | tee -a $OUT/machine.log; fi
  k=$((k + 1))
done
[ $k -gt $MAXPASS ] && log "NOT converged after $MAXPASS passes: the tables of the last pass are used, check $OUT/machine.log"
# the converged tables are the input of the last pass (the last pass only confirms them)
mkdir -p $OUT/final; TF=$OUT/CalTimeCalib_pass$((k > MAXPASS ? MAXPASS : k - 1)).txt; EF=$OUT/CalEnergyCalib_pass$((k > MAXPASS ? MAXPASS : k - 1)).txt
[ $k -eq 1 ] && { TF=$T; EF=$E; }
cp $TF $OUT/final/CalTimeCalib_$RUN.txt; cp $EF $OUT/final/CalEnergyCalib_$RUN.txt
for x in CalSimCrystals CalBaselines response; do f=$(ls $OUT/${x}_pass$((k > MAXPASS ? MAXPASS : k))* 2>/dev/null | head -1); [ -n "$f" ] && cp $f $OUT/final/${x}_$RUN.${f##*.}; done
cat > $OUT/final/README.txt <<EOF
Calorimeter tables of run $RUN made by Offline/MC2Reco/scripts/run_machine.sh on $(date +%Y-%m-%d) (raw list $LIST).
Start: time $TSTART ; energy $ESTART. Passes and convergence: $OUT/machine.log.
CalTimeCalib/CalEnergyCalib: Offline condition format (DbService.textFile + ProditionsService.calCalib.useDb : true;
MC2CaloRecoDigiMaker.useCalCalibTimeOffset : true). CalSimCrystals/CalBaselines/response: detector response for the simulation
(thresholds in ADC and MeV, baselines, noise; pePerMeV from L-R photostatistics is NOT reliable yet).
EOF
log "final tables in $OUT/final"
# ---- physics with the final tables ----
writefcl physics $OUT/final/CalTimeCalib_$RUN.txt $OUT/final/CalEnergyCalib_$RUN.txt "dump, stopped, capture" ""
runpass physics
PY="env -u PYTHONPATH -u PYTHONHOME /usr/bin/python3"   # the analysis scripts need matplotlib (not in the spack python of the Mu2e environment)
STATS=$(ls $OUT/calib_pass*_stats.csv | sort -V | tail -1)
root -l -b -q "Offline/MC2Reco/scripts/capture_to_csv.C(\"$OUT/physics.root\",\"$OUT/physics\")" > /dev/null 2>&1 && $PY Offline/MC2Reco/scripts/mu_capture_analysis.py $OUT/physics $OUT/capture > $OUT/capture_analysis.log 2>&1 || log "capture analysis failed, see $OUT/capture_analysis.log"
root -l -b -q "Offline/MC2Reco/scripts/delayed_to_csv.C(\"$OUT/physics.root\",\"$OUT/delayed.csv\")" > /dev/null 2>&1 && $PY Offline/MC2Reco/scripts/neighbour_delayed_analysis.py $OUT/delayed.csv $OUT/physics_capture.csv $STATS $OUT/neighbour > $OUT/neighbour_analysis.log 2>&1 || log "neighbour analysis failed, see $OUT/neighbour_analysis.log"
log "done"
