#!/usr/bin/env python3
"""Hall Monte Carlo (mc2_digi_reco_capture.fcl) with the cuts of the 125367 data tables of the notebook, STOP candidates
split by truth class, compared with the data.

Candidates (tree capture/capture of MC2MuCapture): kind 0 STOP, kind 1 THROUGH.
Truth class of a candidate (tree truth/stops of MC2StopTruth): the first stopped muon of the same event whose stop
crystal is the candidate crystal or within 50 mm of it in the same disk: mu- DIO (a DIO daughter), mu- capture
(a nuclear-capture daughter), mu+ decay; otherwise "no stop". mu- without DIO or capture daughters are not used.
Observables
  second pulse: maxS > 5, A2 / 23.4 ADC/MeV in the bins 5-10, 10-20, 20-40, > 40 MeV;
  delayed clusters (tree capture/cluster) with E >= 10 MeV in 30-300 ns and in 1-10 us (each candidate counted once
  per window).
MC counts scaled to the 79709 STOP and 157770 THROUGH candidates of 125367; the data numbers are those of the notebook.

Inputs   <digicap dir>/digi_*.root, the channel map DMAP.txt (crystal positions)
Output   JSON: per observable the MC STOP total, by class, MC THROUGH and the data [STOP, THROUGH].
Printed  the same table.

Usage
  PYTHONPATH=$ROOTSYS/lib python3 Offline/MC2Reco/scripts/capture_vs_data_125367.py <digicap dir> <DMAP.txt> <out json>
"""
import argparse
import glob
import json
import math
from collections import Counter, defaultdict

import ROOT

DIO, NCAP = 114, 115                  # ProcessCode ids (Offline MCDataProducts/inc/ProcessCode.hh)
CLASSES = ["mu- capture", "mu- DIO", "mu+ decay", "no stop"]
MATCH_DISTANCE = 50                   # mm, stop crystal - candidate crystal distance for a truth match
ADC_PER_MEV = 23.4                    # energy scale of the second-pulse amplitude used in the notebook tables
CLUSTER_EMIN = 10                     # MeV, delayed cluster energy
STOP_CANDIDATES_125367 = 79709
THROUGH_CANDIDATES_125367 = 157770
# data of run 125367 (notebook): [STOP, THROUGH] candidates with each observable
DATA_125367 = {"5-10": [277, 54], "10-20": [94, 2], "20-40": [12, 0], ">40": [0, 0],
               "cl>=10 30-300ns": [85, 11], "cl>=10 1-10us": [623, 69]}


def parse_args():
    """Command line: digicap directory, channel map, output JSON."""
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("digicap_dir", help="directory with the digi_*.root files of mc2_digi_reco_capture.fcl")
    parser.add_argument("dmap", help="channel map DMAP.txt (crystal positions)")
    parser.add_argument("output", help="output JSON file")
    return parser.parse_args()


def read_crystal_positions(dmap_path):
    """Crystal positions from the channel map: {crystalId: (disk, x, y)} for the CAL rows of disks 0 and 1."""
    positions = {}
    for line in open(dmap_path).read().splitlines()[1:]:
        fields = line.split()
        if len(fields) > 17 and fields[17] == "CAL" and fields[2] in ("0", "1"):
            positions[int(fields[14])] = (int(fields[2]), float(fields[12]), float(fields[13]))
    return positions


def is_near(positions, crystal_a, crystal_b):
    """Both crystals known, in the same disk and within MATCH_DISTANCE."""
    a, b = positions.get(crystal_a), positions.get(crystal_b)
    return bool(a and b and a[0] == b[0] and math.hypot(a[1] - b[1], a[2] - b[2]) <= MATCH_DISTANCE)


def truth_class(stop):
    """Class of a stopped muon from its charge and daughter process codes; None for a mu- without DIO or capture."""
    codes = list(stop.dCode)
    if stop.q < 0:
        return "mu- DIO" if DIO in codes else ("mu- capture" if NCAP in codes else None)
    return "mu+ decay"


def amplitude_bin(energy):
    """Second-pulse energy bin [MeV], None below 5 MeV."""
    if 5 <= energy < 10:
        return "5-10"
    if 10 <= energy < 20:
        return "10-20"
    if 20 <= energy < 40:
        return "20-40"
    if energy >= 40:
        return ">40"
    return None


def cluster_window(dt):
    """Delay window of a delayed cluster [ns], None outside both."""
    if 30 <= dt <= 300:
        return "cl>=10 30-300ns"
    if 1000 <= dt <= 10000:
        return "cl>=10 1-10us"
    return None


def count_file(filename, positions, counts, n_candidates):
    """Add the candidates and observables of one file to counts[kind][(observable, class)] and n_candidates[kind]."""
    f = ROOT.TFile.Open(filename)
    stops_by_event = defaultdict(list)
    for stop in f.Get("truth/stops"):
        label = truth_class(stop)
        if label:
            stops_by_event[stop.event].append((label, stop.cry))
    label_of = {}
    for c in f.Get("capture/capture"):
        label = next((cls for cls, cry in stops_by_event.get(c.event, [])
                      if cry == c.cry or is_near(positions, cry, c.cry)), "no stop")
        label_of[(c.event, c.cry, c.kind)] = label
        n_candidates[c.kind] += 1
        if c.maxS > 5:
            energy_bin = amplitude_bin(c.A2 / ADC_PER_MEV)
            if energy_bin:
                counts[c.kind][(energy_bin, label)] += 1
    # each candidate counted at most once per window
    seen = set()
    for c in f.Get("capture/cluster"):
        if c.E < CLUSTER_EMIN:
            continue
        window = cluster_window(c.dt)
        key = (c.event, c.cry, c.kind, window)
        if window is None or key in seen:
            continue
        seen.add(key)
        counts[c.kind][(window, label_of.get((c.event, c.cry, c.kind), "no stop"))] += 1
    f.Close()


def main():
    """Count all files, scale to 125367, print and write the comparison."""
    args = parse_args()
    positions = read_crystal_positions(args.dmap)
    counts = {0: Counter(), 1: Counter()}
    n_candidates = Counter()
    for filename in sorted(glob.glob(args.digicap_dir + "/digi_*.root")):
        count_file(filename, positions, counts, n_candidates)

    scale_stop = STOP_CANDIDATES_125367 / n_candidates[0]
    scale_through = THROUGH_CANDIDATES_125367 / n_candidates[1]
    result = {}
    print("%-16s %8s %7s %6s %6s %7s %7s | %s"
          % ("observable", "MC STOP", "capture", "DIO", "mu+", "no stop", "MC THR", "data STOP, THROUGH"))
    for observable in DATA_125367:
        parts = [counts[0][(observable, label)] * scale_stop for label in CLASSES]
        through = sum(counts[1][(observable, label)] for label in CLASSES) * scale_through
        result[observable] = {"mc_stop": round(sum(parts)), "by_class": dict(zip(CLASSES, [round(x) for x in parts])),
                              "mc_through": round(through), "data_stop_through": DATA_125367[observable]}
        print("%-16s %8.0f %7.0f %6.0f %6.0f %7.0f %7.0f | %s"
              % (observable, sum(parts), *parts, through, DATA_125367[observable]))
    json.dump(result, open(args.output, "w"), indent=1)


if __name__ == "__main__":
    main()
