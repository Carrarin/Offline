#!/usr/bin/env python3
"""Hall Monte Carlo (Mu2eG4 + CRY, digitized and reconstructed with mc2_digi_reco.fcl) against the cosmic data of run
125329: track rates and observed decays of stopped muons, with no free parameter (normalisation = CRY livetime).

Inputs
  MC   : digi_*.root of run_digi_reco.sh: dump/tracks (MC2CaloCosmicTrackDump; type 1, 2 one-disk, 3 two-disk) and
         stopped/pairs (MC2StoppedMuons). Livetime: sum of "Cosmic livetime" of the Geant4 job logs (job_*.log of the
         production directory).
  Data : 125329 from mc2-ana (output/book/ch05, ch03 results.json, output/cosmic_tracks_run125329/summary.txt): tracks
         of the mc2tools tagger and pairs of cosmic_delayed_coincidence.C (same selection; the data macro also asks
         |tL - tR| < 10 ns per crystal).

Decays observed: near pairs (cls 1) with nD <= 6 and aD <= 3000 ADC, as book chapter 3 / 12; binned Poisson likelihood
fit of N0 exp(-t/tau) + B on dt in [1, 20] us, 0.25 us bins; observed decays = integral of the exponential on
[1, 20] us. MC rates are multiplied by the laser live fraction of the data (1 - ch05.dead_laser_frac).

Outputs (in <out dir>)
  sim_vs_data.json   MC and data rates, their ratio, decays per 1000 one-disk tracks (also printed)
  sim_vs_data.png    decay-time fit of the MC near pairs; data / MC ratios

Usage
  PYTHONPATH=$ROOTSYS/lib python3 Offline/MC2Reco/scripts/sim_vs_data.py \\
      <production dir> <digi dir> <mc2-ana dir> <out dir>
"""
import argparse
import glob
import json
import math
import pathlib
import re

import numpy as np
import ROOT

# decay-time fit range and binning [us]
TMIN = 1.0
TMAX = 20.0
TBIN = 0.25


def parse_args():
    """Command line: production dir, digi dir, mc2-ana dir, output dir."""
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("production_dir", type=pathlib.Path, help="hall production with the Geant4 job_*.log")
    parser.add_argument("digi_dir", type=pathlib.Path, help="directory with the digi_*.root of run_digi_reco.sh")
    parser.add_argument("ana_dir", type=pathlib.Path, help="mc2-ana directory (data of run 125329)")
    parser.add_argument("out_dir", type=pathlib.Path, help="output directory")
    return parser.parse_args()


def read_livetime(production_dir):
    """Sum of the CRY livetimes [s] printed in the Geant4 job logs."""
    return sum(float(m.group(1)) for f in production_dir.glob("job_*.log")
               for m in re.finditer(r"Cosmic livetime: ([0-9.eE+-]+) s", f.read_text()))


def read_mc_tracks_and_pairs(files):
    """Tracks per type (dump/tracks) and the delays [us] of the selected near pairs (stopped/pairs)."""
    tracks_per_type = {1: 0, 2: 0, 3: 0}
    pair_delays = []
    for file_name in files:
        root_file = ROOT.TFile.Open(file_name)
        for track in root_file.Get("dump/tracks"):
            tracks_per_type[track.type] = tracks_per_type.get(track.type, 0) + 1
        # near pairs: cluster class 1, at most 6 delayed crystals, delayed amplitude <= 3000 ADC (book chapter 3)
        for pair in root_file.Get("stopped/pairs"):
            if pair.cls == 1 and pair.nD <= 6 and pair.aD <= 3000:
                pair_delays.append(pair.dt)
        root_file.Close()
    return tracks_per_type, np.array(pair_delays)


def fit_decay_curve(pair_delays):
    """Histogram of the delays in [TMIN, TMAX) and Poisson-likelihood fit N0 exp(-t/tau) + B."""
    hist = ROOT.TH1D("h", "", int((TMAX - TMIN) / TBIN), TMIN, TMAX)
    hist.SetDirectory(0)
    for delay in pair_delays[(pair_delays >= TMIN) & (pair_delays < TMAX)]:
        hist.Fill(delay)
    func = ROOT.TF1("fx", "[0]*exp(-x/[1])+[2]", TMIN, TMAX)
    func.SetParameters(max(1, hist.GetBinContent(1)) * 1.5, 2.2, 0.1)
    func.SetParLimits(1, 0.3, 20)
    func.SetParLimits(2, 0, 1e6)
    fit_result = hist.Fit(func, "LSQN0R")
    return hist, func, fit_result


def decays_in_window(n0, tau):
    """Integral of N0 exp(-t/tau) on [TMIN, TMAX], in pairs (N0 is per bin of TBIN)."""
    return n0 * tau / TBIN * (math.exp(-TMIN / tau) - math.exp(-TMAX / tau))


def decays_error(n0, tau, signal, fit_result):
    """Error on the counted decays from the covariance of (N0, tau, B), with a numerical derivative in tau."""
    cov = fit_result.GetCovarianceMatrix()
    dlog_signal_dtau = (decays_in_window(n0, tau * 1.001) - signal) / (signal * 0.001 * tau)
    gradient = np.array([signal / n0, dlog_signal_dtau * signal, 0.0])
    cov_matrix = np.array([[cov(i, j) for j in range(3)] for i in range(3)])
    return float(np.sqrt(gradient @ cov_matrix @ gradient))


def read_data_125329(ana_dir):
    """Rates of run 125329 from the mc2-ana outputs, and the laser live fraction of the run."""
    ch05 = json.loads((ana_dir / "output/book/ch05/results.json").read_text())
    ch03 = json.loads((ana_dir / "output/book/ch03/results.json").read_text())
    track_summary = (ana_dir / "output/cosmic_tracks_run125329/summary.txt").read_text()
    match = re.search(r"disk 0 only (\d+), disk 1 only (\d+), both disks (\d+)", track_summary)
    run_seconds = ch05["ch05.run_s"]
    eps_laser = 1 - ch05["ch05.dead_laser_frac"]
    data = {"one_disk_Hz": (int(match.group(1)) + int(match.group(2))) / run_seconds,
            "two_disk_Hz": int(match.group(3)) / run_seconds,
            "decays_Hz": ch03["ch03.main.signal"] / ch03["ch03.run_s"],
            "decays_err_Hz": ch03["ch03.main.signal"] ** 0.5 / ch03["ch03.run_s"],
            "run_s": run_seconds}
    return data, eps_laser


def plot(hist, func, live, tau, tau_err, mc, data, out_dir):
    """Left: MC near-pair delays with the fit. Right: data / MC for the three rates."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 2, figsize=(15, 4.6))
    centres = np.array([hist.GetBinCenter(i) for i in range(1, hist.GetNbinsX() + 1)])
    contents = np.array([hist.GetBinContent(i) for i in range(1, hist.GetNbinsX() + 1)])
    ax[0].errorbar(centres, contents, np.sqrt(np.clip(contents, 1, None)), fmt="o", ms=4, color="k",
                   label=f"hall MC, near pairs ({live:.0f} s live)")
    ax[0].plot(centres, [func.Eval(x) for x in centres], color="C0", lw=0.9,
               label=f"fit: τ = {tau:.2f} ± {tau_err:.2f} µs")
    ax[0].set_yscale("log")
    ax[0].set_xlabel("delay of the delayed deposit [µs]")
    ax[0].set_ylabel("pairs / 0.25 µs")
    ax[0].legend(fontsize=7.5)
    ax[0].set_title("stopped muons in the hall Monte Carlo", fontsize=10)

    labels = ["one-disk tracks", "two-disk tracks", "decays 1–20 µs"]
    keys = ["one_disk_Hz", "two_disk_Hz", "decays_Hz"]
    mc_values = [mc[k] for k in keys]
    mc_errors = [mc["one_disk_err_Hz"], mc["two_disk_err_Hz"], mc["decays_err_Hz"]]
    data_values = [data[k] for k in keys]
    xpos = np.arange(3)
    ratio = np.array(data_values) / np.array(mc_values)
    # error bar: relative MC error only
    ratio_err = np.array(mc_errors) / np.array(mc_values) * np.array(data_values) / np.array(mc_values)
    ax[1].errorbar(xpos - 0.1, ratio, ratio_err, fmt="o", ms=5, color="C3", label="data 125329 / hall MC")
    ax[1].axhline(1, color="0.45", lw=0.9)
    ax[1].set_xticks(xpos)
    ax[1].set_xticklabels(labels)
    ax[1].set_ylabel("data / MC")
    ax[1].legend(fontsize=7.5)
    ax[1].set_ylim(0, 1.6)
    ax[1].set_title("no free parameter: normalisation = CRY livetime", fontsize=10)
    fig.suptitle("Hall Monte Carlo (Mu2eG4 + CRY) against run 125329", fontsize=11)
    fig.tight_layout()
    fig.savefig(out_dir / "sim_vs_data.png", dpi=125)


def main():
    args = parse_args()
    ROOT.gROOT.SetBatch(True)
    ROOT.gErrorIgnoreLevel = ROOT.kError
    out_dir = args.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    # Monte Carlo: livetime, tracks, near pairs and their decay-time fit
    live = read_livetime(args.production_dir)
    files = sorted(glob.glob(str(args.digi_dir / "digi_*.root")))
    tracks_per_type, pair_delays = read_mc_tracks_and_pairs(files)
    hist, func, fit_result = fit_decay_curve(pair_delays)
    n0, tau, background = (func.GetParameter(i) for i in range(3))
    tau_err = func.GetParError(1)
    signal = decays_in_window(n0, tau)
    signal_err = decays_error(n0, tau, signal, fit_result)

    # data 125329
    data, eps_laser = read_data_125329(args.ana_dir)

    n_one_disk = tracks_per_type[1] + tracks_per_type[2]
    mc = {"livetime_s": live, "files": len(files), "one_disk_Hz": n_one_disk / live * eps_laser,
          "one_disk_err_Hz": n_one_disk ** 0.5 / live,
          "two_disk_Hz": tracks_per_type[3] / live * eps_laser, "two_disk_err_Hz": tracks_per_type[3] ** 0.5 / live,
          "tracks": tracks_per_type,
          "near_pairs_1_20us": int(hist.Integral()), "tau_us": tau, "tau_err_us": tau_err, "B_per_bin": background,
          "decays_Hz": signal / live * eps_laser, "decays_err_Hz": signal_err / live, "decays_counted": signal}
    results = {"mc": mc, "data": data,
               "ratio_data_mc": {k: data[k] / mc[k] for k in ("one_disk_Hz", "two_disk_Hz", "decays_Hz")},
               "decays_per_1000_onedisk": {"mc": 1000 * mc["decays_Hz"] / mc["one_disk_Hz"],
                                           "data": 1000 * data["decays_Hz"] / data["one_disk_Hz"]}}
    (out_dir / "sim_vs_data.json").write_text(json.dumps(results, indent=1))
    print(json.dumps(results, indent=1))

    plot(hist, func, live, tau, tau_err, mc, data, out_dir)


if __name__ == "__main__":
    main()
