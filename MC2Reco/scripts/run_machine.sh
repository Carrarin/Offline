#!/bin/bash
#
# run_machine.sh: the calorimeter "machine" on one run.
# Calibration passes on the cosmic tracks until the tables converge, then the physics job with the final tables.
#
# Usage (from the muse area, after "muse setup"):
#
#   Offline/MC2Reco/scripts/run_machine.sh <run> <raw list> <out dir> [time start table] [energy start table]
#
#   <raw list>   text file with the raw .art files of the run, one per line
#   <out dir>    must be below the muse area: DbService looks for its text tables in the search path, so the tables are
#                passed to it relative to the area
#
# Options (environment variables):
#
#   MAXPASS=4        maximum number of calibration passes
#   FIXED_ENERGY=1   keep the energy start table (required) in all passes and in the physics: only the times are
#                    calibrated (use it when the energy scale of this HV configuration is already measured)
#
# Start tables:
#
#   time    any existing CalTimeCalib table (e.g. the converged table of a nearby run). Without one, a pass 0 with all
#           offsets at 0 and wide clustering (all crystals of an event in one cluster, L and R merged within 2 us)
#           makes the first table.
#   energy  by default one ADC2MeV for all SiPMs (0.0625 MeV/ADC, crystals equalised); the machine measures the scale
#           from the MIP peak.
#   Never the laser for the times of cosmics (fibre lengths), never the laser for photoelectrons (other wavelength).
#
# Steps:
#
#   1. pass 0 (only without a time start table)
#   2. pass 1: start tables -> calibration only (analyzers dump + calib) -> tables of pass 1. Always WIDE (track finder
#      maxGap 1e9, L/R merging 2 us): a SiPM whose start offset is wrong by more than the 6 ns L/R merging window would
#      otherwise never join a track.
#   3. pass k+1: tables of pass k as conditions -> calibration again, until converged: time corrections with
#      |median| < 0.1 ns and MAD sigma < 0.3 ns, and energy ratio with |median - 1| < 0.005 (times only with
#      FIXED_ENERGY). At most MAXPASS passes.
#   4. final tables: the input tables of the converged pass, copied to <out dir>/final with a README
#   5. physics with the final tables: tracks (dump), stopped muons (Michel spectra per disk, lifetime), mu- captures
#   6. analysis of the captures and of the delayed deposits in the neighbours (figures in capture/ and neighbour/)
#
# Output in <out dir>:
#
#   machine.log                       log of the whole procedure (passes, convergence, summaries)
#   pass<k>.fcl/.log/.root            job of each pass
#   CalTimeCalib_pass<k>.txt, CalEnergyCalib_pass<k>.txt, CalSimCrystals_pass<k>.txt, CalBaselines_pass<k>.txt,
#   response_pass<k>.csv, calib_pass<k>_stats.csv, conv_<k>.txt
#   final/                            converged tables of the run + README.txt with the provenance
#   physics.fcl/.log/.root, stopped_physics.txt, capture/, neighbour/
#
set -e
set -o pipefail   # a check piped into tee must keep its own exit status

RUN=$1
LIST=$(readlink -f "$2")
OUT=$3
TIME_START=$4
ENERGY_START=$5
MAXPASS=${MAXPASS:-4}
FIXED_ENERGY=${FIXED_ENERGY:-0}

if [ -z "$RUN" ] || [ ! -s "$LIST" ] || [ -z "$OUT" ]; then
  echo "usage: $0 <run> <raw list> <out dir> [time start table] [energy start table]"
  exit 1
fi
if [ "$FIXED_ENERGY" = 1 ] && [ ! -s "$ENERGY_START" ]; then
  echo "FIXED_ENERGY=1 needs the energy start table"
  exit 1
fi

mkdir -p "$OUT"
OUT=$(readlink -f "$OUT")
BASE_FCL=Offline/MC2Reco/fcl/mc2_calo_machine.fcl
WIDE_CLUSTERING="physics.producers.MC2CaloCosmicTrackFinder.maxGap : 1e9
physics.producers.MC2CaloHitMaker.time4Merge : 2000."


# ---- helpers -------------------------------------------------------------------------------------------------------

log() {
  echo "[run_machine $RUN] $*" | tee -a "$OUT/machine.log"
}

# Table with the same value for every SiPM: uniform_table <file> <CalTimeCalib|CalEnergyCalib> <value>
uniform_table() {
  python3 - "$1" "$2" "$3" <<'EOF'
import sys
file_name, table, value = sys.argv[1], sys.argv[2], float(sys.argv[3])
column = "timeOffset" if table == "CalTimeCalib" else "ADC2MeV"
with open(file_name, "w") as out:
    out.write(f"TABLE {table}\n#roid,{column}  start: one value for all SiPMs\n")
    for roid in range(2720):
        out.write(f"{roid},{value}\n")
EOF
}

# fcl of one pass: write_fcl <pass name> <time table> <energy table> <analyzers> <extra fcl lines>
write_fcl() {
  local name=$1 time_table=$2 energy_table=$3 analyzers=$4 extra=$5
  local fcl=$OUT/$name.fcl
  cp "$BASE_FCL" "$fcl"
  cat >> "$fcl" <<EOF
# $name of run $RUN (run_machine.sh)
services.ProditionsService.calCalib.useDb : true
services.DbService.textFile : [ "$(realpath --relative-to="$PWD" "$energy_table")", "$(realpath --relative-to="$PWD" "$time_table")" ]
physics.ana : [ $analyzers ]
physics.end_paths : [ ana ]
physics.analyzers.calib.timeTableFile : "$OUT/CalTimeCalib_$name.txt"
physics.analyzers.calib.energyTableFile : "$OUT/CalEnergyCalib_$name.txt"
physics.analyzers.calib.statsFile : "$OUT/calib_${name}_stats.csv"
physics.analyzers.calib.simTableFile : "$OUT/CalSimCrystals_$name.txt"
physics.analyzers.calib.baselineTableFile : "$OUT/CalBaselines_$name.txt"
physics.analyzers.calib.responseFile : "$OUT/response_$name.csv"
physics.analyzers.stopped.summaryFile : "$OUT/stopped_$name.txt"
$extra
EOF
}

# Run one pass and copy its summary lines to the log: run_pass <pass name>
run_pass() {
  local name=$1
  log "running $name"
  if ! mu2e -c "$OUT/$name.fcl" -S "$LIST" -T "$OUT/$name.root" > "$OUT/$name.log" 2>&1; then
    log "FAILED $name, see $OUT/$name.log"
    exit 1
  fi
  grep -E "MC2CaloCalibFromTracks:|MC2StoppedMuons|lifetime|Michel disk|MC2MuCapture:" "$OUT/$name.log" \
    | tee -a "$OUT/machine.log"
}

# Convergence of a pass from its stats csv (exit status 0 = converged); prints one summary line.
converged() {
  python3 - "$OUT/calib_$1_stats.csv" "$FIXED_ENERGY" <<'EOF'
import sys
import numpy as np
stats = np.genfromtxt(sys.argv[1], delimiter=",", names=True)
fixed_energy = sys.argv[2] == "1"
time_rows = stats[stats["timeFromTracks"] == 1]
energy_rows = stats[stats["energyFromTracks"] == 1]
time_median = np.median(time_rows["timeCorrection_ns"])
time_sigma = 1.4826 * np.median(np.abs(time_rows["timeCorrection_ns"] - time_median))
energy_ratio = np.median(energy_rows["newADC2MeV"] / energy_rows["startADC2MeV"])
print(f"time correction median {time_median:.3f} ns, MAD sigma {time_sigma:.3f} ns ({len(time_rows)} SiPMs); "
      f"energy ratio median {energy_ratio:.4f} ({len(energy_rows)} SiPMs)")
times_ok = abs(time_median) < 0.1 and time_sigma < 0.3
energy_ok = fixed_energy or abs(energy_ratio - 1) < 0.005
sys.exit(0 if times_ok and energy_ok else 1)
EOF
}


# ---- start tables --------------------------------------------------------------------------------------------------

if [ -z "$ENERGY_START" ]; then
  ENERGY_START=$OUT/CalEnergyCalib_start_uniform.txt
  uniform_table "$ENERGY_START" CalEnergyCalib 0.0625
  log "energy start: uniform 0.0625 MeV/ADC"
fi


# ---- step 1: pass 0, only without a time start table ---------------------------------------------------------------

if [ -z "$TIME_START" ]; then
  TIME_START=$OUT/CalTimeCalib_start_zero.txt
  uniform_table "$TIME_START" CalTimeCalib 0
  write_fcl pass0 "$TIME_START" "$ENERGY_START" "dump, calib" "$WIDE_CLUSTERING"
  run_pass pass0
  TIME_START=$OUT/CalTimeCalib_pass0.txt
  if [ "$FIXED_ENERGY" != 1 ]; then
    ENERGY_START=$OUT/CalEnergyCalib_pass0.txt
  fi
fi
if [ "$FIXED_ENERGY" = 1 ]; then
  log "time start: $TIME_START ; energy start: $ENERGY_START (fixed)"
else
  log "time start: $TIME_START ; energy start: $ENERGY_START"
fi


# ---- steps 2 and 3: calibration passes until converged -------------------------------------------------------------

time_table=$TIME_START
energy_table=$ENERGY_START
pass=1
while [ $pass -le "$MAXPASS" ]; do
  extra=""
  if [ $pass -eq 1 ]; then
    extra=$WIDE_CLUSTERING
  fi
  write_fcl pass$pass "$time_table" "$energy_table" "dump, calib" "$extra"
  run_pass pass$pass

  time_table=$OUT/CalTimeCalib_pass$pass.txt
  if [ "$FIXED_ENERGY" != 1 ]; then
    energy_table=$OUT/CalEnergyCalib_pass$pass.txt
  fi

  # pass 1 always continues (its corrections are large by construction); later passes stop when converged
  if converged pass$pass > "$OUT/conv_$pass.txt"; then
    tee -a "$OUT/machine.log" < "$OUT/conv_$pass.txt"
    if [ $pass -gt 1 ]; then
      log "converged at pass $pass"
      break
    fi
  else
    tee -a "$OUT/machine.log" < "$OUT/conv_$pass.txt"
  fi
  pass=$((pass + 1))
done
if [ $pass -gt "$MAXPASS" ]; then
  log "NOT converged after $MAXPASS passes: the tables of the last pass are used, check $OUT/machine.log"
fi


# ---- step 4: final tables ------------------------------------------------------------------------------------------
# Converged at pass k: the tables of pass k-1 were its input and pass k only confirmed them.
# Not converged: the tables of the last pass.

if [ $pass -gt "$MAXPASS" ]; then
  table_pass=$MAXPASS
  response_pass=$MAXPASS
else
  table_pass=$((pass - 1))
  response_pass=$pass
fi
mkdir -p "$OUT/final"
final_time=$OUT/CalTimeCalib_pass$table_pass.txt
final_energy=$OUT/CalEnergyCalib_pass$table_pass.txt
if [ $pass -eq 1 ]; then   # MAXPASS=0: no calibration pass, the start tables
  final_time=$time_table
  final_energy=$energy_table
fi
if [ "$FIXED_ENERGY" = 1 ]; then
  final_energy=$ENERGY_START
fi
cp "$final_time" "$OUT/final/CalTimeCalib_$RUN.txt"
cp "$final_energy" "$OUT/final/CalEnergyCalib_$RUN.txt"
for table in CalSimCrystals CalBaselines response; do
  file=$(ls "$OUT/${table}_pass$response_pass"* 2>/dev/null | head -1)
  if [ -n "$file" ]; then
    cp "$file" "$OUT/final/${table}_$RUN.${file##*.}"
  fi
done
cat > "$OUT/final/README.txt" <<EOF
Calorimeter tables of run $RUN made by Offline/MC2Reco/scripts/run_machine.sh on $(date +%Y-%m-%d) (raw list $LIST).
Start: time $TIME_START ; energy $ENERGY_START$([ "$FIXED_ENERGY" = 1 ] && echo " (fixed)"). Passes and convergence: $OUT/machine.log.
CalTimeCalib/CalEnergyCalib: Offline condition format (DbService.textFile + ProditionsService.calCalib.useDb : true;
MC2CaloRecoDigiMaker.useCalCalibTimeOffset : true). CalSimCrystals/CalBaselines/response: detector response for the simulation
(thresholds in ADC and MeV, baselines, noise; pePerMeV from L-R photostatistics is NOT reliable yet).
EOF
log "final tables in $OUT/final"


# ---- step 5: physics with the final tables -------------------------------------------------------------------------

write_fcl physics "$OUT/final/CalTimeCalib_$RUN.txt" "$OUT/final/CalEnergyCalib_$RUN.txt" "dump, stopped, capture" ""
run_pass physics


# ---- step 6: analysis of the captures and of the delayed deposits --------------------------------------------------
# The analysis scripts need matplotlib, which is not in the python of the Mu2e environment: system python.

PYTHON="env -u PYTHONPATH -u PYTHONHOME /usr/bin/python3"
last_stats=$(ls "$OUT"/calib_pass*_stats.csv | sort -V | tail -1)

if root -l -b -q "Offline/MC2Reco/scripts/capture_to_csv.C(\"$OUT/physics.root\",\"$OUT/physics\")" > /dev/null 2>&1 \
   && $PYTHON Offline/MC2Reco/scripts/mu_capture_analysis.py "$OUT/physics" "$OUT/capture" \
        > "$OUT/capture_analysis.log" 2>&1; then
  :
else
  log "capture analysis failed, see $OUT/capture_analysis.log"
fi

if root -l -b -q "Offline/MC2Reco/scripts/delayed_to_csv.C(\"$OUT/physics.root\",\"$OUT/delayed.csv\")" > /dev/null 2>&1 \
   && $PYTHON Offline/MC2Reco/scripts/neighbour_delayed_analysis.py "$OUT/delayed.csv" "$OUT/physics_capture.csv" \
        "$last_stats" "$OUT/neighbour" > "$OUT/neighbour_analysis.log" 2>&1; then
  :
else
  log "neighbour analysis failed, see $OUT/neighbour_analysis.log"
fi

log "done"
