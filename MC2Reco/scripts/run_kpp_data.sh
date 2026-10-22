#!/bin/bash
# run_kpp_data.sh -- reconstruction of the KPP combined runs from the archive (gm2ita): for every raw file of a list, as soon as it is in
# /data1/mu2e/mc2/archive/raw/kpp/<run>/ with its md5 in md5sums.txt, run in the EL9 container
#   kpp_combined.fcl (+ conditions set) -> products/kpp_reco/<run>/<ver>/rec_<run>_<file>.art
#   kpp_ntuple.fcl                       -> products/kpp_ntuple/<run>/<ver>/nts_<run>_<file>.root
#   kpp_dump_events.py                   -> products/kpp_events/<run>/<ver>/ev_<run>_<file>.json
# with at most NJOBS files at a time (priority nice NICE), and keep a MANIFEST.json per run and product.
# usage: run_kpp_data.sh <muse area> <list of raw file names (first column)> <version, e.g. v001> <conditions set, e.g. conditions/sets/kpp_calo_v001.txt>
#   env: NJOBS (16), NICE (5), CHECK (60 s), WORK (/data1/mu2e/carrarin/kpp_data_<version>: logs and job fcl)
M=$(readlink -f $1); LIST=$(readlink -f $2); VER=$3; SET=$4
[ -d "$M/Offline/MC2Reco" ] && [ -s "$LIST" ] && [ -n "$VER" ] && [ -n "$SET" ] || { echo "usage: $0 <muse area> <list> <version> <conditions set>"; exit 1; }
NJOBS=${NJOBS:-16}; NICE=${NICE:-5}; CHECK=${CHECK:-60}; WORK=${WORK:-/data1/mu2e/carrarin/kpp_data_$VER}
A=/data1/mu2e/mc2/archive/raw/kpp; P=/data1/mu2e/mc2/products; IMG=/cvmfs/singularity.opensciencegrid.org/fermilab/fnal-wn-el9:latest
COMMIT=$(git -C $M/Offline rev-parse --short HEAD); DIRTY=$(git -C $M/Offline status --porcelain --untracked-files=no | wc -l)
mkdir -p $WORK; cd $WORK
log() { echo "$(date '+%F %T') $*" | tee -a $WORK/run.log; }
manifest() {   # $1 = run: one MANIFEST.json per product version, listing the files done so far with their input md5
  local run=$1
  for prod in kpp_reco kpp_ntuple kpp_events; do
    local d=$P/$prod/$run/$VER; [ -d $d ] || continue
    python3 - "$d" "$prod" "$run" <<EOF
import json, glob, os, sys, datetime
d, prod, run = sys.argv[1:4]
md5 = dict(l.split()[::-1] for l in open("$A/" + run + "/md5sums.txt")) if os.path.exists("$A/" + run + "/md5sums.txt") else {}
done = sorted(os.path.basename(f)[:-5] for f in glob.glob("$WORK/" + run + "_*.done"))
inputs = [dict(file=f"raw.mu2e.cosmics.kpp.{b}-001.art", md5=md5.get(f"raw.mu2e.cosmics.kpp.{b}-001.art")) for b in done]
m = dict(product=f"{prod}/{run}/$VER", description={"kpp_reco": "reconstruction of the KPP combined run: tracker KKLine and CRV (PassN), MC2 calorimeter chain",
         "kpp_ntuple": "EventNtuple of kpp_reco", "kpp_events": "tracks with in-time calorimeter crystals and CRV coincidences (kpp_dump_events.py)"}[prod],
         fork="github.com/Carrarin/Offline mc2-calo $COMMIT" + (" (local changes)" if int("$DIRTY") else ""),
         environment="muse setup $M (backing DataJob 1.11.0, PassN Mu2e/PassN $(git -C $M/PassN rev-parse --short HEAD 2>/dev/null)) in fnal-wn-el9",
         fcl=["Offline/MC2Reco/fcl/kpp/kpp_combined.fcl", "Offline/MC2Reco/fcl/kpp/kpp_ntuple.fcl"], conditions="$SET",
         geometry="Offline/MC2Reco/data/geom/geom_extracted_mc2kpp.txt", command="Offline/MC2Reco/scripts/run_kpp_data.sh $M $LIST $VER $SET",
         inputs=inputs, n_files=len(inputs), updated=datetime.datetime.now().astimezone().isoformat(timespec="seconds"), author="$(id -un)")
json.dump(m, open(os.path.join(d, "MANIFEST.json"), "w"), indent=1)
EOF
  done
}
job() {   # $1 = raw file name
  local f=$1 b=$(echo $1 | sed -E 's/raw\.mu2e\.cosmics\.kpp\.(.*)-001\.art/\1/'); local run=${b%%_*}
  local R=$P/kpp_reco/$run/$VER N=$P/kpp_ntuple/$run/$VER E=$P/kpp_events/$run/$VER; mkdir -p $R $N $E
  printf '#include "Offline/MC2Reco/fcl/kpp/kpp_combined.fcl"\nservices.DbService.textFile : [ @sequence::services.DbService.textFile, "%s" ]\n' $SET > $WORK/job_$b.fcl
  nice -n $NICE apptainer exec -B /cvmfs,/data1 $IMG bash -c "source /cvmfs/mu2e.opensciencegrid.org/setupmu2e-art.sh >/dev/null 2>&1; cd $M; muse setup > /dev/null 2>&1
    export MU2E_SEARCH_PATH=\$MU2E_SEARCH_PATH:$M:/data1/mu2e/mc2
    mu2e -c $WORK/job_$b.fcl -s $A/$run/$f -o $R/rec_$b.art -T /dev/null > $WORK/rec_$b.log 2>&1 || exit 1
    mu2e -c Offline/MC2Reco/fcl/kpp/kpp_ntuple.fcl -s $R/rec_$b.art -T $N/nts_$b.root > $WORK/nts_$b.log 2>&1 || exit 2
    python3 Offline/MC2Reco/scripts/kpp_dump_events.py $E/ev_$b.json $N/nts_$b.root > $WORK/ev_$b.log 2>&1 || exit 3" < /dev/null
  local s=$?
  if [ $s = 0 ]; then touch $WORK/$b.done; cp $WORK/rec_$b.log $R/; log "$b done ($(tail -1 $WORK/ev_$b.log))"; manifest $run
  else touch $WORK/$b.failed; log "$b FAILED at step $s, see $WORK"; fi
}
log "start: $(wc -l < $LIST) files, version $VER, conditions $SET, fork $COMMIT (local changes: $DIRTY), NJOBS $NJOBS"
declare -A PID
while true; do
  for b in "${!PID[@]}"; do kill -0 ${PID[$b]} 2>/dev/null || unset PID[$b]; done
  todo=0
  for f in $(awk '{print $1}' $LIST); do
    b=$(echo $f | sed -E 's/raw\.mu2e\.cosmics\.kpp\.(.*)-001\.art/\1/'); run=${b%%_*}
    [ -e $WORK/$b.done ] || [ -e $WORK/$b.failed ] && continue
    todo=$((todo + 1)); [ -n "${PID[$b]}" ] && continue
    [ ${#PID[@]} -ge $NJOBS ] && continue
    grep -q " $f$" $A/$run/md5sums.txt 2>/dev/null || continue   # not copied (and checked) yet
    job $f & PID[$b]=$!; log "$b started"
  done
  [ $todo = 0 ] && break
  sleep $CHECK
done
log "all files of $LIST processed"
