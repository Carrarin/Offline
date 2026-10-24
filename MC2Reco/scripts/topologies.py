#!/usr/bin/env python3
"""Topologies of the stopped muons in the hall Monte Carlo, in simple terms, and two ways of finding the mu- capture:
  A. stop crystal: second pulse of MC2MuCapture (maxS > 5, A2 > 30 ADC)
  B. neighbours: a crystal within 80 mm of the stop crystal (same disk), not on the muon track (no hit within +-30 ns
     of the stop), with a CaloHit of both SiPMs 30-500 ns after the stop hit; energy = the largest such crystal and
     the sum.
Truth class of each STOP candidate as capture_signature.py (stopped muon of the same event whose crystal is the
candidate crystal or within 50 mm of it, same disk). Event displays: one example per topology, crystals around the
stop crystal in three time windows (prompt -30..+30 ns, 30-500 ns, 1-20 us).

Inputs
  <digicap dir>  digi_*.root of mc2_digi_reco_capture_hits.fcl (trees truth/stops, hits/hits, capture/capture)
  <DMAP.txt>     channel map with the crystal positions (columns: 2 disk, 12 x, 13 y, 14 crystal id, 17 "CAL")
Outputs (in <out dir>)
  topologies.json  per truth class: fraction of STOP candidates found by A, by B (>= 3, >= 5 MeV), A or B, B without
                   a late electron, late electron (CaloHit > 5 MeV within 80 mm, 1-20 us); truth energies 0-500 ns;
                   purity of the mu- capture and counts scaled to the 79709 STOP candidates of run 125367 (printed)
  topologies.png   event displays, one row per topology

Usage
  PYTHONPATH=$ROOTSYS/lib python3 Offline/MC2Reco/scripts/topologies.py <digicap dir> <DMAP.txt> <out dir>
"""
import argparse
import glob
import json
import math
import pathlib
from collections import defaultdict

import numpy as np
import ROOT

# ProcessCode ids (Offline MCDataProducts/inc/ProcessCode.hh)
DIO = 114
NCAP = 115

TRUTH_CLASSES = ["mu- capture", "mu- DIO", "mu+ decay", "no stop"]
# STOP candidates of run 125367 (raw files 901-931, notebook)
DATA_STOP_125367 = 79709
# ADC per MeV of the second-pulse amplitude A2
A2_ADC_PER_MEV = 22.7
# examples kept per truth class for the event displays
MAX_EXAMPLES = 300
# half size of a crystal square in the displays, half width of the display around the stop crystal [mm]
HALF_CRYSTAL = 17
HALF_VIEW = 260


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


def distance(positions, crystal1, crystal2):
    """Distance [mm] of two crystals on the same disk; 1e9 if on different disks or unknown."""
    a, b = positions.get(crystal1), positions.get(crystal2)
    return 1e9 if (a is None or b is None or a[0] != b[0]) else math.hypot(a[1] - b[1], a[2] - b[2])


def read_truth(root_file):
    """{event: [(class, crystal, inCry, eStopWin, eNearWin, eFarWin)]} of the stopped muons with a truth class."""
    truth = defaultdict(list)
    for stop in root_file.Get("truth/stops"):
        codes = list(stop.dCode)
        truth_class = ("mu- DIO" if DIO in codes else ("mu- capture" if NCAP in codes else None)) \
            if stop.q < 0 else "mu+ decay"
        if truth_class:
            truth[stop.event].append((truth_class, stop.cry, bool(stop.inCry), list(stop.eStopWin),
                                      list(stop.eNearWin), list(stop.eFarWin)))
    return truth


def read_event_hits(root_file):
    """{event: (crystals, times, energies, nSiPM)} of all the CaloHits of each event (last entry if repeated)."""
    hits = {}
    for event_hits in root_file.Get("hits/hits"):
        hits[event_hits.event] = (list(event_hits.cry), list(event_hits.t), list(event_hits.e),
                                  list(event_hits.nSiPM))
    return hits


def candidate_row(candidate, truth, event_hits, positions):
    """Truth class and topology of one STOP candidate; (row, t0), or None if the stop crystal has no hit."""
    truth_class, in_crystal, truth_windows = "no stop", False, None
    for cl, crystal, inc, stop_win, near_win, far_win in truth.get(candidate.event, []):
        if crystal == candidate.cry or distance(positions, crystal, candidate.cry) <= 50:
            truth_class, in_crystal, truth_windows = cl, inc, (stop_win, near_win, far_win)
            break
    hit_crystals, hit_times, hit_energies, hit_nsipm = event_hits
    stop_times = [t for cr, t, en in zip(hit_crystals, hit_times, hit_energies) if cr == candidate.cry]
    if not stop_times:
        return None
    # time of the stop hit: the hit of the stop crystal closest to t1 of the candidate (or the first one)
    t0 = min(stop_times, key=lambda t: abs(t - candidate.t1)) if candidate.t1 > 0 else min(stop_times)
    # crystals on the muon track: a hit within 30 ns of the stop
    prompt = {cr for cr, t in zip(hit_crystals, hit_times) if abs(t - t0) < 30}
    # B: delayed energy per neighbour crystal (within 80 mm, off the track, both SiPMs, 30-500 ns)
    neighbour_energy = defaultdict(float)
    for cr, t, en, ns in zip(hit_crystals, hit_times, hit_energies, hit_nsipm):
        if cr == candidate.cry or cr in prompt or ns < 2 or distance(positions, cr, candidate.cry) > 80:
            continue
        if 30 <= t - t0 <= 500:
            neighbour_energy[cr] += en
    # late electron: a hit above 5 MeV within 80 mm, 1-20 us after the stop
    late = any(1000 <= t - t0 <= 20000 and distance(positions, cr, candidate.cry) <= 80 and en > 5
               for cr, t, en in zip(hit_crystals, hit_times, hit_energies))
    row = dict(cls=truth_class, inCry=in_crystal, A=(candidate.maxS > 5 and candidate.A2 > 30),
               A2=candidate.A2 / A2_ADC_PER_MEV, nbMax=max(neighbour_energy.values(), default=0.),
               nbSum=sum(neighbour_energy.values()), late=late, truthWin=truth_windows)
    return row, t0


def read_candidates(digicap_dir, positions):
    """Rows of the STOP candidates of all the files, and up to MAX_EXAMPLES examples per truth class
    (file, event, crystal, t0, row) for the displays."""
    rows = []
    examples = defaultdict(list)
    for file_name in sorted(glob.glob(str(digicap_dir / "digi_*.root"))):
        root_file = ROOT.TFile.Open(file_name)
        truth = read_truth(root_file)
        hits = read_event_hits(root_file)
        for candidate in root_file.Get("capture/capture"):
            if candidate.kind != 0 or candidate.event not in hits:
                continue
            result = candidate_row(candidate, truth, hits[candidate.event], positions)
            if result is None:
                continue
            row, t0 = result
            rows.append(row)
            if len(examples[row["cls"]]) < MAX_EXAMPLES:
                examples[row["cls"]].append((file_name, candidate.event, candidate.cry, t0, row))
        root_file.Close()
    return rows, examples


def summarise(rows):
    """The content of topologies.json."""
    results = {}
    for cl in TRUTH_CLASSES:
        class_rows = [r for r in rows if r["cls"] == cl]
        n = len(class_rows)
        if not n:
            continue
        results[cl] = {"n": n,
                       "A_stop_crystal_second_pulse": round(np.mean([r["A"] for r in class_rows]), 4),
                       "B_neighbour_ge3MeV": round(np.mean([r["nbMax"] >= 3 for r in class_rows]), 4),
                       "B_neighbour_ge5MeV": round(np.mean([r["nbMax"] >= 5 for r in class_rows]), 4),
                       "A_or_B3": round(np.mean([r["A"] or r["nbMax"] >= 3 for r in class_rows]), 4),
                       "B3_and_no_late_electron": round(np.mean([(r["nbMax"] >= 3) and not r["late"]
                                                                 for r in class_rows]), 4),
                       "late_electron_1_20us": round(np.mean([r["late"] for r in class_rows]), 4)}
        truth_windows = [r["truthWin"] for r in class_rows if r["truthWin"]]
        if truth_windows:
            a = np.array(truth_windows)   # [n, 3 (stop, near, far), 4 windows]
            # windows 0 and 1: 0-500 ns
            results[cl]["truth_mean_MeV_0_500ns"] = {"stop": round(float((a[:, 0, 0] + a[:, 0, 1]).mean()), 2),
                                                     "near": round(float((a[:, 1, 0] + a[:, 1, 1]).mean()), 2),
                                                     "far": round(float((a[:, 2, 0] + a[:, 2, 1]).mean()), 2)}
    present = [c for c in TRUTH_CLASSES if c in results]
    n_stop = sum(results[c]["n"] for c in present)
    for key in ("A_stop_crystal_second_pulse", "B_neighbour_ge3MeV", "B3_and_no_late_electron"):
        total = sum(results[c]["n"] * results[c][key] for c in present)
        captures = results["mu- capture"]["n"] * results["mu- capture"][key]
        results[f"purity_{key}"] = round(captures / total, 3) if total else None
        results[f"per_79709_STOP_{key}"] = {c: round(DATA_STOP_125367 / n_stop * results[c]["n"] * results[c][key])
                                            for c in present}
    return results


def first_example(examples, cl, condition):
    """First example of class cl whose row satisfies condition, None if there is none."""
    for example in examples[cl]:
        if condition(example[4]):
            return example


def pick_examples(examples):
    """One example per topology (display title -> example or None)."""
    picks = {}
    picks["μ⁻ capture (stopped in the crystal)"] = first_example(
        examples, "mu- capture", lambda r: r["inCry"] and r["nbMax"] >= 3)
    picks["μ⁻ decay in orbit (DIO)"] = (first_example(examples, "mu- DIO", lambda r: r["nbMax"] >= 5)
                                        or first_example(examples, "mu- DIO", lambda r: True))
    picks["μ⁺ decay (Michel electron after µs)"] = first_example(
        examples, "mu+ decay", lambda r: r["inCry"] and r["late"])
    picks["muon stopped outside the crystals"] = first_example(
        examples, "mu+ decay", lambda r: (not r["inCry"]) and r["late"])
    picks["no stop: the muon leaves the disk"] = first_example(examples, "no stop", lambda r: True)
    return picks


def draw_window(a, plt, Rectangle, positions, hits, t0, stop_crystal, lo, hi):
    """Crystals around the stop crystal with their energy in [lo, hi) ns after t0, stop crystal in red."""
    hit_crystals, hit_times, hit_energies = hits
    disk0, x0, y0 = positions[stop_crystal]
    energy = defaultdict(float)
    for cr, t, en in zip(hit_crystals, hit_times, hit_energies):
        if lo <= t - t0 < hi and cr in positions and positions[cr][0] == disk0:
            energy[cr] += en
    for cr, (disk, xx, yy) in positions.items():
        if disk != disk0 or abs(xx - x0) > HALF_VIEW or abs(yy - y0) > HALF_VIEW:
            continue
        en = energy.get(cr, 0.)
        face = plt.cm.viridis(min(en / 40, 1)) if en > 0.5 else "white"
        a.add_patch(Rectangle((xx - HALF_CRYSTAL, yy - HALF_CRYSTAL), 2 * HALF_CRYSTAL, 2 * HALF_CRYSTAL,
                              facecolor=face, edgecolor="0.8", lw=0.5))
        if en > 0.5:
            a.text(xx, yy, f"{en:.0f}", ha="center", va="center", fontsize=6, color="w" if en < 25 else "k")
    a.add_patch(Rectangle((x0 - HALF_CRYSTAL, y0 - HALF_CRYSTAL), 2 * HALF_CRYSTAL, 2 * HALF_CRYSTAL, fill=False,
                          edgecolor="C3", lw=1.8))
    a.set_xlim(x0 - HALF_VIEW, x0 + HALF_VIEW)
    a.set_ylim(y0 - HALF_VIEW, y0 + HALF_VIEW)
    a.set_aspect("equal")
    a.set_xticks([])
    a.set_yticks([])


def plot_event_displays(examples, positions, out_dir):
    """topologies.png: one row per topology, three time windows."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle
    picks = pick_examples(examples)
    windows = [("prompt −30…+30 ns: the muon", -30, 30), ("30–500 ns: capture / DIO", 30, 500),
               ("1–20 µs: Michel electron", 1000, 20000)]
    fig, ax = plt.subplots(len(picks), 3, figsize=(13, 4.2 * len(picks)))
    for i, (title, example) in enumerate(picks.items()):
        if example is None:
            for j in range(3):
                ax[i, j].axis("off")
            continue
        file_name, event, stop_crystal, t0, _ = example
        # all the hits of the event, read again from its file
        root_file = ROOT.TFile.Open(file_name)
        event_hits = next(hh for hh in root_file.Get("hits/hits") if hh.event == event)
        hits = (list(event_hits.cry), list(event_hits.t), list(event_hits.e))
        root_file.Close()
        for j, (window_label, lo, hi) in enumerate(windows):
            a = ax[i, j]
            draw_window(a, plt, Rectangle, positions, hits, t0, stop_crystal, lo, hi)
            a.set_title(window_label if i == 0 else "", fontsize=10)
            if j == 0:
                a.set_ylabel(title, fontsize=10)
    fig.suptitle("Hall Monte Carlo: one event per topology (red square = crystal where the track ends; "
                 "numbers = MeV)", fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.975))
    fig.savefig(out_dir / "topologies.png", dpi=110)


def main():
    args = parse_args()
    out_dir = args.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)
    positions = read_crystal_positions(args.dmap)
    rows, examples = read_candidates(args.digicap_dir, positions)
    results = summarise(rows)
    (out_dir / "topologies.json").write_text(json.dumps(results, indent=1))
    print(json.dumps(results, indent=1))
    plot_event_displays(examples, positions, out_dir)


if __name__ == "__main__":
    main()
