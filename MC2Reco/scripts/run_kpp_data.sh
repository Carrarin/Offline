#!/bin/bash
#
# run_kpp_data.sh: reconstruction of the KPP combined runs (tracker + CRV + MC2 calorimeter) from the archive on gm2ita.
#
# For every raw file of a list, as soon as it is in /data1/mu2e/mc2/archive/raw/kpp/<run>/ with its md5 in
# md5sums.txt (i.e. copied and checked), three steps run in the EL9 container:
#
#   kpp_combined.fcl (+ conditions set) -> products/kpp_reco/<run>/<version>/rec_<run>_<file>.art
#   kpp_ntuple.fcl                      -> products/kpp_ntuple/<run>/<version>/nts_<run>_<file>.root
#   kpp_dump_events.py                  -> products/kpp_events/<run>/<version>/ev_<run>_<file>.json
#
# with at most NJOBS files at a time (priority nice NICE). After every file the MANIFEST.json of each product version
# of its run is rewritten with the files done so far and the md5 of their inputs. The script ends when every file of
# the list is done or failed; files not yet in the archive are waited for (it checks again every CHECK seconds).
#
# Usage (on gm2ita, outside the container; meant to be started detached, e.g. with setsid nohup ... &):
#
#   run_kpp_data.sh <muse area> <list> <version> <conditions set>
#
#   e.g. Offline/MC2Reco/scripts/run_kpp_data.sh /data1/mu2e/carrarin/mc2offline kpp_files.txt v001 \
#          conditions/sets/kpp_calo_v001.txt
#
# Arguments:
#
#   <muse area>        muse area with Offline/MC2Reco (built), PassN, EventNtuple
#   <list>             text file with the raw file names raw.mu2e.cosmics.kpp.<run>_<file>-001.art (first column)
#   <version>          version of the products, e.g. v001 (directory <version> of each product)
#   <conditions set>   DbService text file added to the tables of kpp_combined.fcl, relative to /data1/mu2e/mc2
#                      (e.g. conditions/sets/kpp_calo_v001.txt)
#
# Options (environment variables):
#
#   NJOBS=16    maximum number of files at a time
#   NICE=5      nice of the jobs
#   CHECK=60    seconds between checks
#   WORK=/data1/mu2e/carrarin/kpp_data_<version>   logs, job fcl and end markers (give it absolute)
#
# Output:
#
#   /data1/mu2e/mc2/products/{kpp_reco,kpp_ntuple,kpp_events}/<run>/<version>/   the files above + MANIFEST.json
#   (kpp_reco also gets a copy of the reconstruction log rec_<run>_<file>.log)
#   <WORK>/run.log                                         log of the whole production
#   <WORK>/job_<b>.fcl, rec_<b>.log, nts_<b>.log, ev_<b>.log   job fcl and logs of file b = <run>_<file>
#   <WORK>/<b>.done / <b>.failed                           end markers (a failed file is not retried)
#

MUSE_AREA=$(readlink -f $1)
LIST=$(readlink -f $2)
VER=$3
CONDITIONS_SET=$4

if ! { [ -d "$MUSE_AREA/Offline/MC2Reco" ] && [ -s "$LIST" ] && [ -n "$VER" ] && [ -n "$CONDITIONS_SET" ]; }; then
  echo "usage: $0 <muse area> <list> <version> <conditions set>"
  exit 1
fi

NJOBS=${NJOBS:-16}
NICE=${NICE:-5}
CHECK=${CHECK:-60}
WORK=${WORK:-/data1/mu2e/carrarin/kpp_data_$VER}
ARCHIVE=/data1/mu2e/mc2/archive/raw/kpp
PRODUCTS=/data1/mu2e/mc2/products
IMAGE=/cvmfs/singularity.opensciencegrid.org/fermilab/fnal-wn-el9:latest

# Commit of the fork and number of modified tracked files, for the log and the manifests
COMMIT=$(git -C $MUSE_AREA/Offline rev-parse --short HEAD)
DIRTY=$(git -C $MUSE_AREA/Offline status --porcelain --untracked-files=no | wc -l)

mkdir -p $WORK
cd $WORK


# ---- helpers ---------------------------------------------------------------------------------------------------------

log() {
  echo "$(date '+%F %T') $*" | tee -a $WORK/run.log
}

# raw.mu2e.cosmics.kpp.<run>_<file>-001.art -> <run>_<file>
file_tag() {
  echo $1 | sed -E 's/raw\.mu2e\.cosmics\.kpp\.(.*)-001\.art/\1/'
}

# One MANIFEST.json per product version of a run, listing the files done so far with the md5 of their input:
# manifest <run>
manifest() {
  local run=$1
  local passn_commit=$(git -C $MUSE_AREA/PassN rev-parse --short HEAD 2>/dev/null)
  local product
  for product in kpp_reco kpp_ntuple kpp_events; do
    local product_dir=$PRODUCTS/$product/$run/$VER
    if [ ! -d $product_dir ]; then
      continue
    fi
    python3 - "$product_dir" "$product" "$run" "$ARCHIVE" "$WORK" "$VER" "$COMMIT" "$DIRTY" "$MUSE_AREA" \
      "$passn_commit" "$CONDITIONS_SET" "$LIST" "$(id -un)" <<'EOF'
import datetime
import glob
import json
import os
import sys

(product_dir, product, run, archive, work, version, commit, dirty, muse_area,
 passn_commit, conditions_set, file_list, author) = sys.argv[1:14]

# md5 of the raw files of the run: lines "<md5>  <file name>"
md5_file = archive + "/" + run + "/md5sums.txt"
md5 = dict(line.split()[::-1] for line in open(md5_file)) if os.path.exists(md5_file) else {}

# files done so far: the end markers <run>_<file>.done in the work directory
done = sorted(os.path.basename(f)[:-5] for f in glob.glob(work + "/" + run + "_*.done"))
inputs = [dict(file=f"raw.mu2e.cosmics.kpp.{tag}-001.art", md5=md5.get(f"raw.mu2e.cosmics.kpp.{tag}-001.art"))
          for tag in done]

descriptions = {
    "kpp_reco": "reconstruction of the KPP combined run: tracker KKLine and CRV (PassN), MC2 calorimeter chain",
    "kpp_ntuple": "EventNtuple of kpp_reco",
    "kpp_events": "tracks with in-time calorimeter crystals and CRV coincidences (kpp_dump_events.py)",
}
manifest = dict(
    product=f"{product}/{run}/{version}",
    description=descriptions[product],
    fork="github.com/Carrarin/Offline mc2-calo " + commit + (" (local changes)" if int(dirty) else ""),
    environment=f"muse setup {muse_area} (backing DataJob 1.11.0, PassN Mu2e/PassN {passn_commit}) in fnal-wn-el9",
    fcl=["Offline/MC2Reco/fcl/kpp/kpp_combined.fcl", "Offline/MC2Reco/fcl/kpp/kpp_ntuple.fcl"],
    conditions=conditions_set,
    geometry="Offline/MC2Reco/data/geom/geom_extracted_mc2kpp.txt",
    command=f"Offline/MC2Reco/scripts/run_kpp_data.sh {muse_area} {file_list} {version} {conditions_set}",
    inputs=inputs,
    n_files=len(inputs),
    updated=datetime.datetime.now().astimezone().isoformat(timespec="seconds"),
    author=author,
)
json.dump(manifest, open(os.path.join(product_dir, "MANIFEST.json"), "w"), indent=1)
EOF
  done
}

# Reconstruction, ntuple and event dump of one raw file, in the container: job <raw file name>
job() {
  local raw_file=$1
  local tag=$(file_tag $1)
  local run=${tag%%_*}
  local reco_dir=$PRODUCTS/kpp_reco/$run/$VER
  local ntuple_dir=$PRODUCTS/kpp_ntuple/$run/$VER
  local events_dir=$PRODUCTS/kpp_events/$run/$VER
  mkdir -p $reco_dir $ntuple_dir $events_dir

  # job fcl: kpp_combined.fcl with the conditions set added to its DbService tables
  {
    printf '#include "Offline/MC2Reco/fcl/kpp/kpp_combined.fcl"\n'
    printf 'services.DbService.textFile : [ @sequence::services.DbService.textFile, "%s" ]\n' $CONDITIONS_SET
  } > $WORK/job_$tag.fcl

  # the three steps in the container; the exit status says which one failed (1 reco, 2 ntuple, 3 event dump)
  local reco_file=$reco_dir/rec_$tag.art
  local ntuple_file=$ntuple_dir/nts_$tag.root
  local events_file=$events_dir/ev_$tag.json
  local commands="source /cvmfs/mu2e.opensciencegrid.org/setupmu2e-art.sh >/dev/null 2>&1
    cd $MUSE_AREA
    muse setup > /dev/null 2>&1
    export MU2E_SEARCH_PATH=\$MU2E_SEARCH_PATH:$MUSE_AREA:/data1/mu2e/mc2
    mu2e -c $WORK/job_$tag.fcl -s $ARCHIVE/$run/$raw_file -o $reco_file -T /dev/null > $WORK/rec_$tag.log 2>&1 || exit 1
    mu2e -c Offline/MC2Reco/fcl/kpp/kpp_ntuple.fcl -s $reco_file -T $ntuple_file > $WORK/nts_$tag.log 2>&1 || exit 2
    python3 Offline/MC2Reco/scripts/kpp_dump_events.py $events_file $ntuple_file > $WORK/ev_$tag.log 2>&1 || exit 3"
  nice -n $NICE apptainer exec -B /cvmfs,/data1 $IMAGE bash -c "$commands" < /dev/null
  local status=$?

  if [ $status = 0 ]; then
    touch $WORK/$tag.done
    cp $WORK/rec_$tag.log $reco_dir/
    log "$tag done ($(tail -1 $WORK/ev_$tag.log))"
    manifest $run
  else
    touch $WORK/$tag.failed
    log "$tag FAILED at step $status, see $WORK"
  fi
}


# ---- main loop: every CHECK seconds, until every file of the list is done or failed ----------------------------------

log "start: $(wc -l < $LIST) files, version $VER, conditions $CONDITIONS_SET, fork $COMMIT (local changes: $DIRTY)," \
  "NJOBS $NJOBS"
declare -A PID   # <run>_<file> -> pid of the running job
while true; do

  # forget the jobs that ended (their .done/.failed marker is written by the job itself)
  for tag in "${!PID[@]}"; do
    if ! kill -0 ${PID[$tag]} 2>/dev/null; then
      unset PID[$tag]
    fi
  done

  # start the files that are in the archive with their md5, while there is room; count those still to do
  todo=0
  for raw_file in $(awk '{print $1}' $LIST); do
    tag=$(file_tag $raw_file)
    run=${tag%%_*}
    if [ -e $WORK/$tag.done ] || [ -e $WORK/$tag.failed ]; then
      continue
    fi
    todo=$((todo + 1))
    if [ -n "${PID[$tag]}" ]; then
      continue   # running
    fi
    if [ ${#PID[@]} -ge $NJOBS ]; then
      continue   # no room
    fi
    if ! grep -q " $raw_file$" $ARCHIVE/$run/md5sums.txt 2>/dev/null; then
      continue   # not copied (and checked) yet
    fi
    job $raw_file &
    PID[$tag]=$!
    log "$tag started"
  done

  if [ $todo = 0 ]; then
    break
  fi
  sleep $CHECK
done

log "all files of $LIST processed"
