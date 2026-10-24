#!/usr/bin/env python3
"""Signature of mu- captures and mu- decays in orbit (DIO) in the hall Monte Carlo, with the selection used on the data
(MC2MuCapture) and the truth of the same events (MC2StopTruth), from the files of mc2_digi_reco_capture.fcl.

Candidate (tree capture/capture): kind 0 STOP (last crystal of a one-disk track, next crystal empty), kind 1 THROUGH
(middle crystal).
Truth class of a candidate: a stopped muon of the same event whose stop crystal (or, if it stopped outside the
crystals, the last crystal it crossed) is within 50 mm of the candidate crystal in the same disk: mu- capture, mu- DIO,
mu+ decay; otherwise "no stop".
Observables (as on the data):
  second pulse in the candidate crystal: found = maxS > 5 and A2 > 30 ADC (MC2MuCapture example thresholds);
      A2 in MeV = A2 / 22.7;
  delayed cluster (tree capture/cluster, crystals within 80 mm + second pulse): dt in [30, 300] ns with E > 5 MeV /
      > 15 MeV;
  Michel: nMichel >= 1 (CaloHit within 80 mm, 1-20 us later).
Scaled to the data: per STOP candidate, and to the 79709 STOP / 157770 THROUGH of run 125367 (raw files 901-931,
notebook).

Inputs
  <digicap dir>  digi_*.root of mc2_digi_reco_capture.fcl (trees truth/stops, capture/capture, capture/cluster)
  <DMAP.txt>     channel map with the crystal positions (columns: 2 disk, 12 x, 13 y, 14 crystal id, 17 "CAL")
Outputs (in <out dir>)
  capture_signature.json  candidates per truth class, signature of each class, scaled to 125367 (also printed)
  capture_signature.png   second-pulse amplitude, delayed-cluster energy and second-pulse delay of the STOP
                          candidates, by truth class

Usage
  PYTHONPATH=$ROOTSYS/lib python3 Offline/MC2Reco/scripts/capture_signature.py <digicap dir> <DMAP.txt> <out dir>
"""
import argparse
import glob
import json
import math
import pathlib
from collections import Counter, defaultdict

import numpy as np
import ROOT

# ProcessCode ids (Offline MCDataProducts/inc/ProcessCode.hh)
DIO = 114
NCAP = 115

TRUTH_CLASSES = ["mu- capture", "mu- DIO", "mu+ decay", "no stop"]
# truth match: stop crystal within this distance of the candidate crystal, same disk [mm]
MATCH_RADIUS = 50.
# second pulse found (MC2MuCapture example thresholds)
SECOND_PULSE_MAXS = 5
SECOND_PULSE_A2_ADC = 30
# ADC per MeV of the second-pulse amplitude A2
A2_ADC_PER_MEV = 22.7
# STOP and THROUGH candidates of run 125367 (raw files 901-931, notebook)
DATA_STOP_125367 = 79709
DATA_THROUGH_125367 = 157770


def parse_args():
    """Command line: digicap dir, DMAP file, output dir."""
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("digicap_dir", type=pathlib.Path, help="directory with digi_*.root of the capture job")
    parser.add_argument("dmap", type=pathlib.Path, help="DMAP.txt with the crystal positions")
    parser.add_argument("out_dir", type=pathlib.Path, help="output directory")
    return parser.parse_args()


def read_crystal_positions(dmap_path):
    """{crystal id: (disk, x, y)} of the calorimeter crystals of the channel map (header line skipped)."""
    positions = {}
    for line in dmap_path.read_text().splitlines()[1:]:
        fields = line.split()
        if len(fields) > 17 and fields[17] == "CAL" and fields[2] in ("0", "1"):
            positions[int(fields[14])] = (int(fields[2]), float(fields[12]), float(fields[13]))
    return positions


def is_near(positions, crystal1, crystal2, radius=MATCH_RADIUS):
    """True if both crystals are known, on the same disk and within radius [mm]."""
    a, b = positions.get(crystal1), positions.get(crystal2)
    return a is not None and b is not None and a[0] == b[0] and math.hypot(a[1] - b[1], a[2] - b[2]) <= radius


def read_truth(root_file):
    """{event: [(class, crystal, inCry, eTot, eWin, eStopWin)]} of the stopped muons with a truth class."""
    truth = defaultdict(list)
    for stop in root_file.Get("truth/stops"):
        codes, times = list(stop.dCode), list(stop.dT)
        if stop.q < 0:
            t_capture = min([t for c, t in zip(codes, times) if c == NCAP], default=-1.0)
            t_dio = min([t for c, t in zip(codes, times) if c == DIO], default=-1.0)
            truth_class = "mu- DIO" if t_dio >= 0 else ("mu- capture" if t_capture >= 0 else None)
        else:
            truth_class = "mu+ decay"
        if truth_class:
            truth[stop.event].append((truth_class, stop.cry, bool(stop.inCry), stop.eTot, list(stop.eWin),
                                      list(stop.eStopWin)))
    return truth


def read_clusters(root_file):
    """{(event, crystal, kind): [(dt, E, Estop)]} of the delayed clusters of each candidate."""
    clusters = defaultdict(list)
    for cluster in root_file.Get("capture/cluster"):
        clusters[(cluster.event, cluster.cry, cluster.kind)].append((cluster.dt, cluster.E, cluster.Estop))
    return clusters


def candidate_row(candidate, truth, clusters, positions):
    """One candidate: its truth class (first matching stopped muon) and its observables."""
    label, in_crystal, e_true, e_true_win, e_true_stop = "no stop", False, 0., [0] * 4, [0] * 4
    for truth_class, crystal, inc, e_tot, e_win, e_stop_win in truth.get(candidate.event, []):
        if crystal == candidate.cry or is_near(positions, crystal, candidate.cry):
            label, in_crystal, e_true, e_true_win, e_true_stop = truth_class, inc, e_tot, e_win, e_stop_win
            break
    candidate_clusters = clusters.get((candidate.event, candidate.cry, candidate.kind), [])
    # delayed clusters 30-300 ns after the candidate, above 5 and 15 MeV; largest energy in the window
    return dict(kind=candidate.kind, cls=label, inCry=in_crystal, eTrue=e_true, eTrueWin=e_true_win,
                eTrueStop=e_true_stop, maxS=candidate.maxS, A2=candidate.A2, tau=candidate.tau, Q=candidate.Q,
                amp=candidate.amp, nMichel=candidate.nMichel,
                cl5=any(30 <= x[0] <= 300 and x[1] > 5 for x in candidate_clusters),
                cl15=any(30 <= x[0] <= 300 and x[1] > 15 for x in candidate_clusters),
                clE=max([x[1] for x in candidate_clusters if 30 <= x[0] <= 300], default=0.))


def read_candidates(digicap_dir, positions):
    """Rows of all the candidates of all the files, file by file in name order."""
    rows = []
    for file_name in sorted(glob.glob(str(digicap_dir / "digi_*.root"))):
        root_file = ROOT.TFile.Open(file_name)
        truth = read_truth(root_file)
        clusters = read_clusters(root_file)
        for candidate in root_file.Get("capture/capture"):
            rows.append(candidate_row(candidate, truth, clusters, positions))
        root_file.Close()
    return rows


def class_signature(class_rows):
    """Fractions of candidates of one truth class with each observable, and medians / means of the truth."""
    with_pulse = [r for r in class_rows
                  if r["maxS"] > SECOND_PULSE_MAXS and r["A2"] > SECOND_PULSE_A2_ADC]
    a2_median = (round(float(np.median([r["A2"] / A2_ADC_PER_MEV for r in with_pulse])), 1)
                 if with_pulse else None)
    tau_median = round(float(np.median([r["tau"] for r in with_pulse])), 0) if with_pulse else None
    return {"n": len(class_rows), "in_crystal": round(np.mean([r["inCry"] for r in class_rows]), 3),
            "second_pulse": round(len(with_pulse) / len(class_rows), 3),
            "second_pulse_A2_MeV_median": a2_median,
            "second_pulse_tau_ns_median": tau_median,
            "cluster_30_300ns_gt5MeV": round(np.mean([r["cl5"] for r in class_rows]), 3),
            "cluster_30_300ns_gt15MeV": round(np.mean([r["cl15"] for r in class_rows]), 3),
            "michel_1_20us": round(np.mean([r["nMichel"] >= 1 for r in class_rows]), 3),
            "true_E_products_median_MeV": round(float(np.median([r["eTrue"] for r in class_rows])), 1),
            # truth energy in the windows 0 and 1 (0-50 ns, 50-500 ns)
            "true_E_0_500ns_mean_MeV": round(float(np.mean([r["eTrueWin"][0] + r["eTrueWin"][1]
                                                            for r in class_rows])), 2)}


def summarise_kind(rows, kind):
    """Candidates of one kind (0 STOP, 1 THROUGH): counts and fractions per truth class, signature per class."""
    kind_rows = [r for r in rows if r["kind"] == kind]
    n = len(kind_rows)
    counts = Counter(r["cls"] for r in kind_rows)
    out = {"candidates": n, "classes": {c: counts.get(c, 0) for c in TRUTH_CLASSES},
           "fraction": {c: round(counts.get(c, 0) / n, 4) for c in TRUTH_CLASSES}, "signature": {}}
    for c in TRUTH_CLASSES:
        class_rows = [r for r in kind_rows if r["cls"] == c]
        if not class_rows:
            continue
        out["signature"][c] = class_signature(class_rows)
    return out


def summarise(rows):
    """The content of capture_signature.json."""
    results = {}
    for kind, kind_name in ((0, "STOP"), (1, "THROUGH")):
        results[kind_name] = summarise_kind(rows, kind)
    n_stop, n_through = results["STOP"]["candidates"], results["THROUGH"]["candidates"]
    results["scaled_to_125367"] = {
        "STOP": {c: round(DATA_STOP_125367 * results["STOP"]["classes"][c] / n_stop) for c in TRUTH_CLASSES},
        "THROUGH": {c: round(DATA_THROUGH_125367 * results["THROUGH"]["classes"][c] / n_through)
                    for c in TRUTH_CLASSES}}
    expected = {}
    for c in TRUTH_CLASSES:
        signature = results["STOP"]["signature"].get(c)
        if signature:
            expected[c] = {k: round(results["scaled_to_125367"]["STOP"][c] * signature[k])
                           for k in ("second_pulse", "cluster_30_300ns_gt5MeV", "cluster_30_300ns_gt15MeV",
                                     "michel_1_20us")}
    results["expected_in_125367_STOP"] = expected
    return results


def plot(rows, out_dir):
    """capture_signature.png: STOP candidates by truth class, three observables (normalised histograms)."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.4))
    colours = {"mu- capture": "C0", "mu- DIO": "C1", "mu+ decay": "C3", "no stop": "0.45"}
    stop_rows = [r for r in rows if r["kind"] == 0]
    for c in TRUTH_CLASSES:
        class_rows = [r for r in stop_rows if r["cls"] == c]
        if not class_rows:
            continue
        ax[0].hist([r["A2"] / A2_ADC_PER_MEV for r in class_rows if r["maxS"] > SECOND_PULSE_MAXS],
                   bins=np.linspace(0, 40, 41), histtype="step", color=colours[c], density=True,
                   label=f"{c} ({len(class_rows)})")
        ax[1].hist([r["clE"] for r in class_rows if r["clE"] > 0], bins=np.linspace(0, 60, 61), histtype="step",
                   color=colours[c], density=True, label=c)
        ax[2].hist([r["tau"] for r in class_rows
                    if r["maxS"] > SECOND_PULSE_MAXS and r["A2"] > SECOND_PULSE_A2_ADC],
                   bins=np.linspace(0, 450, 46), histtype="step", color=colours[c], density=True, label=c)
    ax[0].set_xlabel("second pulse A2 [MeV] (maxS > 5)")
    ax[0].set_ylabel("normalised")
    ax[0].set_yscale("log")
    ax[0].legend(fontsize=7.5)
    ax[0].set_title("second pulse in the stop crystal", fontsize=10)
    ax[1].set_xlabel("delayed cluster energy, 30–300 ns [MeV]")
    ax[1].set_yscale("log")
    ax[1].legend(fontsize=7.5)
    ax[1].set_title("cluster in the crystals around (and the stop crystal)", fontsize=10)
    ax[2].set_xlabel("delay of the second pulse [ns]")
    ax[2].legend(fontsize=7.5)
    ax[2].set_title("when the second pulse comes", fontsize=10)
    fig.suptitle("Hall Monte Carlo, STOP candidates of MC2MuCapture by truth class", fontsize=11)
    fig.tight_layout()
    fig.savefig(out_dir / "capture_signature.png", dpi=125)


def main():
    args = parse_args()
    ROOT.gROOT.SetBatch(True)
    out_dir = args.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)
    positions = read_crystal_positions(args.dmap)
    rows = read_candidates(args.digicap_dir, positions)
    results = summarise(rows)
    (out_dir / "capture_signature.json").write_text(json.dumps(results, indent=1))
    print(json.dumps(results, indent=1))
    plot(rows, out_dir)


if __name__ == "__main__":
    main()
