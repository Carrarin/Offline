#!/usr/bin/env python3
"""Checks on the Monte Carlo truth of the stopped muons (MC2StopTruth trees of the hall production): where they stop,
what their decay / capture products leave in the crystals and when, and whether the counts close against the decays
counted by the reconstruction of the same production (sim_vs_data.py).

Inputs
  <truth dir>            truth_*.root with the tree truth/stops of MC2StopTruth (one entry per stopped muon)
  <livetime s>           CRY livetime of the production [s] (as in sim_vs_data.json, mc.livetime_s)
  <reco-counted decays>  decays counted by the reconstruction (sim_vs_data.json, mc.decays_counted); only copied to
                         the output, to compare with the predicted number

Checks
  - stops of mu+ / mu- in the crystals and outside, all and with more than 5 MeV of products in the crystals; rates
  - mu+ stopped in a crystal: energy of the products
  - decays at 1-20 us with products > 5 MeV, and the same times 0.785 (window factor of the prediction), to compare
    with the decays counted by the reconstruction
  - lifetimes from an unbinned maximum likelihood: mu+ (decays, 0-20 us) and mu- stopped in a crystal (capture or DIO
    time, 0-2 us)
  - mu- nuclear captures in a crystal: energy of the products, by time window, in the stop crystal and around

Outputs (in <out dir>)
  stop_truth_check.json  the numbers above (also printed)
  stop_truth.npz         the per-muon columns (q, inCry, endX/Y/Z, eTot, eStop, eNear, eFar, ePrompt, tdec, tcap,
                         nCryMu) and the energies per time window (win: all crystals, swin: stop crystal)
  stop_truth.png         mu- lifetime in CsI, product energies, when the energy of a capture arrives

Usage (plain ROOT 6.36 Python, not inside the muse environment)
  PYTHONPATH=$ROOTSYS/lib python3 Offline/MC2Reco/scripts/stop_truth_check.py \\
      <truth dir> <livetime s> <reco-counted decays> <out dir>
  e.g. (livetime and decays from sim_vs_data.json) ... stop_truth_check.py <production>/truth 36000 1234.5 out
"""
import argparse
import glob
import json
import math
import pathlib

import numpy as np
import ROOT

# ProcessCode ids (Offline MCDataProducts/inc/ProcessCode.hh)
DECAY = 14
DIO = 114
NCAP = 115

# per-muon scalar branches of truth/stops copied to the columns, in this order
SCALAR_BRANCHES = ("q", "inCry", "endX", "endY", "endZ", "eTot", "eStop", "eNear", "eFar", "ePrompt", "nCryMu")
# all the columns (the order is the order of the arrays in stop_truth.npz)
COLUMNS = ("q", "inCry", "endX", "endY", "endZ", "eTot", "eStop", "eNear", "eFar", "ePrompt", "tdec", "tcap", "nCryMu")

# factor from the decays at 1-20 us with visible products to the decays observed by the reconstruction
# (json key predicted_observed_x_window0785)
WINDOW_FRACTION = 0.785


def parse_args():
    """Command line: truth dir, livetime [s], reco-counted decays, output dir."""
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("truth_dir", type=pathlib.Path, help="directory with the truth_*.root of MC2StopTruth")
    parser.add_argument("livetime", type=float, help="CRY livetime of the production [s]")
    parser.add_argument("reco_counted", type=float, help="decays counted by the reconstruction (sim_vs_data.py)")
    parser.add_argument("out_dir", type=pathlib.Path, help="output directory")
    return parser.parse_args()


def read_stops(truth_dir):
    """Per-muon columns (numpy arrays) and the energies per time window (all crystals, stop crystal)."""
    chain = ROOT.TChain("truth/stops")
    for file_name in sorted(glob.glob(str(truth_dir / "truth_*.root"))):
        chain.Add(file_name)
    columns = {k: [] for k in COLUMNS}
    energy_windows = []
    stop_crystal_windows = []
    for stop in chain:
        codes, times = list(stop.dCode), list(stop.dT)
        for k in SCALAR_BRANCHES:
            columns[k].append(getattr(stop, k))
        # time of the decay (free decay or DIO) and of the nuclear capture, -1 if not there
        columns["tdec"].append(min([t for c, t in zip(codes, times) if c in (DECAY, DIO)], default=-1.0))
        columns["tcap"].append(min([t for c, t in zip(codes, times) if c == NCAP], default=-1.0))
        energy_windows.append(list(stop.eWin))
        stop_crystal_windows.append(list(stop.eStopWin))
    arrays = {k: np.array(v) for k, v in columns.items()}
    return arrays, np.array(energy_windows), np.array(stop_crystal_windows)


def tau_ml(times, lo, hi):
    """Lifetime from the maximum of the likelihood of an exponential truncated to [lo, hi), scanned on a grid.

    Returns (tau, tau / sqrt(n), n)."""
    times = times[(times >= lo) & (times < hi)]
    grid = np.linspace(max(1e-3, 0.002 * hi), 5 * hi, 20000)
    log_likelihood = [-len(times) * math.log(x * (math.exp(-lo / x) - math.exp(-hi / x))) - times.sum() / x
                      for x in grid]
    tau = grid[int(np.argmax(log_likelihood))]
    return float(tau), float(tau / math.sqrt(len(times))), int(len(times))


def summarise(d, energy_windows, stop_crystal_windows, live, reco_counted):
    """All the numbers of stop_truth_check.json, and the arrays needed by the plots."""
    q, in_crystal, e_total = d["q"], d["inCry"].astype(bool), d["eTot"]
    charges_and_places = [(qq, nm, lab, m) for qq, nm in ((1, "mu+"), (-1, "mu-"))
                          for lab, m in (("inCry", in_crystal), ("outside", ~in_crystal))]
    results = {}
    results["stops"] = {f"{nm}_{lab}": int(((q == qq) & m).sum()) for qq, nm, lab, m in charges_and_places}
    results["stops_visible5MeV"] = {f"{nm}_{lab}": int(((q == qq) & m & (e_total > 5)).sum())
                                    for qq, nm, lab, m in charges_and_places}
    results["rates_Hz"] = {k: v / live for k, v in results["stops"].items()}

    mu_plus_in_crystal = (q > 0) & in_crystal
    results["mu+_inCry_products"] = {
        "median_MeV": float(np.median(e_total[mu_plus_in_crystal])),
        "frac_gt20": float(np.mean(e_total[mu_plus_in_crystal] > 20)),
        "mean_stop_near_far": [float(d[k][mu_plus_in_crystal].mean()) for k in ("eStop", "eNear", "eFar")]}

    # decays in the window of the reconstruction (1-20 us) with visible products
    late = (d["tdec"] >= 1000) & (d["tdec"] < 20000)
    results["decays_1_20us_products_gt5MeV"] = {
        lab: int((late & mm & (e_total > 5)).sum())
        for lab, mm in (("inCry", in_crystal), ("all", np.ones_like(in_crystal)))}
    results["predicted_observed_x_window0785"] = {
        k: WINDOW_FRACTION * v for k, v in results["decays_1_20us_products_gt5MeV"].items()}
    results["reco_counted"] = reco_counted

    # lifetimes: mu+ decay times; mu- in a crystal: capture time if captured, else DIO time
    mu_plus_decay_times = d["tdec"][(q > 0) & (d["tdec"] >= 0)]
    captured = (q < 0) & in_crystal & (d["tcap"] >= 0)
    mu_minus_times = np.where(d["tcap"] >= 0, d["tcap"], d["tdec"])[
        (q < 0) & in_crystal & ((d["tcap"] >= 0) | (d["tdec"] >= 0))]
    results["tau_mu+_ns"] = tau_ml(mu_plus_decay_times, 0, 20000)
    results["tau_mu-_inCry_ns"] = tau_ml(mu_minus_times, 0, 2000)

    results["mu-_inCry_capture"] = {
        "n": int(captured.sum()),
        "of_mu-_inCry": int(((q < 0) & in_crystal).sum()),
        "dio": int(((q < 0) & in_crystal & (d["tdec"] >= 0)).sum()),
        "E_frac_gt1_gt5_gt10": [float(np.mean(e_total[captured] > x)) for x in (1, 5, 10)],
        "E_median_84_99": [float(np.percentile(e_total[captured], p)) for p in (50, 84, 99)],
        "E_by_window_mean": energy_windows[captured].mean(0).tolist(),
        "E_stop_crystal_by_window_mean": stop_crystal_windows[captured].mean(0).tolist(),
        "E_stop_near_far_mean": [float(d[k][captured].mean()) for k in ("eStop", "eNear", "eFar")],
        # windows 0 and 1 are 0-50 ns and 50-500 ns
        "frac_gt5MeV_stop_crystal_within_500ns": float(np.mean(
            stop_crystal_windows[captured][:, 0] + stop_crystal_windows[captured][:, 1] > 5))}
    return results, mu_minus_times, captured, in_crystal, e_total


def plot(results, d, energy_windows, stop_crystal_windows, mu_minus_times, captured, in_crystal, e_total, out_dir):
    """stop_truth.png: mu- lifetime in CsI, product energies, mean energy per capture by time window."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    q = d["q"]
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.4))
    tau, tau_err = results["tau_mu-_inCry_ns"][0], results["tau_mu-_inCry_ns"][1]
    ax[0].hist(mu_minus_times, bins=np.linspace(0, 600, 61), histtype="step", color="k",
               label=f"μ⁻ stopped in a crystal: τ = {tau:.0f} ± {tau_err:.0f} ns")
    ax[0].set_yscale("log")
    ax[0].set_xlabel("capture or DIO time after the stop [ns]")
    ax[0].set_ylabel("μ⁻ / 10 ns")
    ax[0].legend(fontsize=7.5)
    ax[0].set_title("Geant4 (ShieldingM): μ⁻ lifetime in CsI", fontsize=10)

    bins = np.linspace(0, 80, 81)
    ax[1].hist(e_total[captured], bins=bins, histtype="step", color="C0", label="μ⁻ nuclear capture")
    ax[1].hist(e_total[(q > 0) & in_crystal & (d["tdec"] >= 0)], bins=bins, histtype="step", color="C3",
               label="μ⁺ decay")
    ax[1].set_yscale("log")
    ax[1].set_xlabel("energy of the products in the crystals [MeV, Geant4]")
    ax[1].set_ylabel("stopped muons / MeV")
    ax[1].legend(fontsize=7.5)
    ax[1].set_title("what the products leave (no detector response)", fontsize=10)

    labels = ["0–50 ns", "50–500 ns", "0.5–20 µs", "> 20 µs"]
    xpos = np.arange(4)
    ax[2].bar(xpos - 0.2, energy_windows[captured].mean(0), 0.4, color="C0", label="all crystals")
    ax[2].bar(xpos + 0.2, stop_crystal_windows[captured].mean(0), 0.4, color="C2", label="stop crystal")
    ax[2].set_xticks(xpos)
    ax[2].set_xticklabels(labels)
    ax[2].set_ylabel("mean energy per capture [MeV]")
    ax[2].legend(fontsize=7.5)
    ax[2].set_title("μ⁻ capture: when the energy arrives", fontsize=10)
    fig.suptitle("Monte Carlo truth of the stopped muons, hall production", fontsize=11)
    fig.tight_layout()
    fig.savefig(out_dir / "stop_truth.png", dpi=125)


def main():
    args = parse_args()
    ROOT.gROOT.SetBatch(True)
    out_dir = args.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    d, energy_windows, stop_crystal_windows = read_stops(args.truth_dir)
    results, mu_minus_times, captured, in_crystal, e_total = summarise(
        d, energy_windows, stop_crystal_windows, args.livetime, args.reco_counted)
    (out_dir / "stop_truth_check.json").write_text(json.dumps(results, indent=1))
    print(json.dumps(results, indent=1))
    np.savez(out_dir / "stop_truth.npz", **d, win=energy_windows, swin=stop_crystal_windows)

    plot(results, d, energy_windows, stop_crystal_windows, mu_minus_times, captured, in_crystal, e_total, out_dir)


if __name__ == "__main__":
    main()
