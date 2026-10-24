#!/usr/bin/env python3
"""Each reconstructed near pair of the hall Monte Carlo matched to the truth of the same event: which stopped muon, if
any, has a decay / capture at that delay.

Pairs (tree stopped/pairs of MC2StoppedMuons), cuts of book chapter 3: cls 1 (near pair), nD <= 6 crystals and
aD <= 3000 ADC in the delayed cluster, dt in [1, 20) us.
Truth (tree truth/stops of MC2StopTruth, file truth_<n>.root of the same job number as digi_<n>.root): for every stopped
muon the first decay (Decay or DIO) and the first nuclear capture time. A pair is a "mu+ decay" or "mu- DIO" (in or
outside the crystals) if a stopped muon of the event decays within 0.2 us of the pair dt (first match wins); otherwise
the event has stopped muons without a decay at that delay (flagged when a mu- is captured in the event), or none.

Inputs   <digi dir>/digi_*.root, <truth dir>/truth_*.root
Output   JSON (also printed): number of pairs; per class [pairs, fraction, mean dt in us], most common first;
         number of events with 1, 2, 3, >= 4 pairs.

Usage
  PYTHONPATH=$ROOTSYS/lib python3 Offline/MC2Reco/scripts/pairs_vs_truth.py <digi dir> <truth dir> <out json>
  e.g. ... pairs_vs_truth.py /data1/mu2e/mc2/sim/hall/v001/digi /data1/mu2e/mc2/sim/hall/v001/truth pairs_truth.json
"""
import argparse
import glob
import json
import pathlib
from collections import Counter, defaultdict

import numpy as np
import ROOT

DECAY, DIO, NCAP = 14, 114, 115       # ProcessCode ids (Offline MCDataProducts/inc/ProcessCode.hh)
DELAY_MATCH_US = 0.2                  # us, |truth decay time - pair dt| to match a pair to a decay
MAX_PAIRS_PER_EVENT_BIN = 4           # events with 4 or more pairs counted together


def parse_args():
    """Command line: digi directory, truth directory, output JSON."""
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("digi_dir", help="directory with the digi_*.root files (stopped/pairs)")
    parser.add_argument("truth_dir", help="directory with the truth_*.root files (truth/stops)")
    parser.add_argument("output", help="output JSON file")
    return parser.parse_args()


def first_time(codes, times, wanted):
    """Earliest time of a daughter whose process code is in wanted, -1 if none."""
    return min([t for c, t in zip(codes, times) if c in wanted], default=-1.0)


def read_truth(truth_file):
    """Stopped muons by event: (charge, in a crystal, first decay time [ns], first capture time [ns], energy [MeV])."""
    stops_by_event = defaultdict(list)
    f = ROOT.TFile.Open(str(truth_file))
    for e in f.Get("truth/stops"):
        codes, times = list(e.dCode), list(e.dT)
        decay_time = first_time(codes, times, (DECAY, DIO))
        capture_time = first_time(codes, times, (NCAP,))
        stops_by_event[e.event].append((e.q, bool(e.inCry), decay_time, capture_time, e.eTot))
    f.Close()
    return stops_by_event


def is_selected(pair):
    """Near pair with the cuts of book chapter 3."""
    return pair.cls == 1 and pair.nD <= 6 and pair.aD <= 3000 and 1 <= pair.dt < 20


def classify(pair_dt, stops):
    """Truth class of a pair with delay pair_dt [us] given the stopped muons of its event."""
    label = "no stopped muon in the event"
    for charge, in_crystal, decay_time, capture_time, energy in stops:
        if decay_time >= 0 and abs(decay_time / 1000 - pair_dt) < DELAY_MATCH_US:
            label = (("mu+ decay" if charge > 0 else "mu- DIO")
                     + (" (in crystal)" if in_crystal else " (outside crystals)"))
            break
    else:
        if stops:
            captured = any(s[0] < 0 and s[3] >= 0 for s in stops)
            label = "stopped muon(s), no decay at that delay" + (" [mu- capture in event]" if captured else "")
    return label


def main():
    """Classify all selected pairs, write and print the summary."""
    args = parse_args()
    digi_dir, truth_dir, output = pathlib.Path(args.digi_dir), pathlib.Path(args.truth_dir), pathlib.Path(args.output)
    pairs_per_class = Counter()
    pairs_per_event = Counter()
    dt_per_class = defaultdict(list)
    for digi_file in sorted(glob.glob(str(digi_dir / "digi_*.root"))):
        job = digi_file.split("_")[-1].split(".")[0]
        stops_by_event = read_truth(truth_dir / f"truth_{job}.root")
        f = ROOT.TFile.Open(digi_file)
        for pair in f.Get("stopped/pairs"):
            if not is_selected(pair):
                continue
            event = int(pair.entry)          # entry = art event number
            pairs_per_event[(job, event)] += 1
            label = classify(pair.dt, stops_by_event.get(event, []))
            pairs_per_class[label] += 1
            dt_per_class[label].append(pair.dt)
        f.Close()

    total = sum(pairs_per_class.values())
    multiplicity = Counter(min(v, MAX_PAIRS_PER_EVENT_BIN) for v in pairs_per_event.values())
    result = {"pairs": total,
              "classes": {k: [v, round(v / total, 3), round(float(np.mean(dt_per_class[k])), 2)]
                          for k, v in pairs_per_class.most_common()},
              "pairs_per_event": {str(k): v for k, v in sorted(multiplicity.items())}}
    output.write_text(json.dumps(result, indent=1))
    print(json.dumps(result, indent=1))


if __name__ == "__main__":
    main()
