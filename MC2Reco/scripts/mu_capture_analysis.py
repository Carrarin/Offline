"""mu- capture from the output of MC2MuCapture (capture_to_csv.C): second pulses in the STOP crystals over the THROUGH
controls, corrected for the detection efficiency measured by injection, and the decay time.

Samples (table <PREFIX>_capture.csv, one row per candidate crystal)
  STOP     kind 0: the crystal where a one-disk track ends; also with the mu+ veto (nMichel == 0: no Michel electron
           after the stop, which removes most mu+ and keeps the mu- that are captured)
  THROUGH  kind 1: a crystal crossed by the track (control: no muon stops there)
Selection
  second pulse found if maxS > SCUT and A2 > A2MIN ADC (and tau >= 0); delay = tau.
  Controls: THROUGH reweighted to the first-pulse amplitude distribution of the STOP sample (the fake rate depends on
  the amplitude), 20 quantile bins of the amplitude of STOP + THROUGH.
Excess and efficiency
  Excess per STOP crystal: f_stop(t) - f_through(t) (fractions per 20 ns delay bin). Efficiency eff(t) = fraction of
  the injected pulses (table <PREFIX>_inject.csv, injected A2 above A2MIN) found with |tau - injDelay| < 15 ns.
  excess / eff fitted with N exp(-t/tauC) + C on [TLO, THI] ns, for all STOP and for STOP with the mu+ veto.

Outputs (in OUTDIR)
  c1_maxS.png   maxS of STOP, STOP with veto, THROUGH (reweighted)
  c2_delay.png  delay of the second pulses, STOP against THROUGH
  c3_eff.png    efficiency from injection against the injected delay
  c4_fit.png    efficiency-corrected excess with the exponential fit
  summary.json  sample sizes, cuts, fraction found, fit results (also printed)

Usage (run_machine.sh calls it with PREFIX and OUTDIR only)
  python3 mu_capture_analysis.py PREFIX OUTDIR [SCUT=5] [A2MIN=30] [TLO=20] [THI=600]
  e.g. python3 Offline/MC2Reco/scripts/mu_capture_analysis.py out/physics out/capture 5 30 20 600
"""
import argparse
import json
import pathlib
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, str(pathlib.Path(__file__).parent))
from crv_calo_coinc import curve_fit  # noqa: E402

BIN_WIDTH = 20.               # ns, delay bins
MAX_DELAY = 1500              # ns, end of the delay histograms
N_AMPLITUDE_BINS = 20         # quantile bins of the first-pulse amplitude for the reweighting
INJECTION_MATCH = 15          # ns, |tau - injected delay| for an injected pulse to count as found
TAU_START = 80.               # ns, starting value of the capture time in the fit


def parse_args():
    """Command line: table prefix, output directory and the optional cuts as positional numbers."""
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("prefix",
                        help="prefix of the capture_to_csv.C tables (<prefix>_capture.csv, <prefix>_inject.csv)")
    parser.add_argument("outdir", help="output directory (created if missing)")
    parser.add_argument("scut", nargs="?", type=float, default=5., help="cut on maxS (default 5)")
    parser.add_argument("a2min", nargs="?", type=float, default=30., help="minimum A2 [ADC] (default 30)")
    parser.add_argument("tlo", nargs="?", type=float, default=20., help="start of the fit [ns] (default 20)")
    parser.add_argument("thi", nargs="?", type=float, default=600., help="end of the fit [ns] (default 600)")
    return parser.parse_args()


def amplitude_weights(reference, control, amplitude_edges):
    """Per-row weights of the control sample that give it the first-pulse amplitude distribution of the reference.

    Normalised so that the weights of the control sum (about) to its number of rows.
    """
    hist_reference, _ = np.histogram(reference["amp"], amplitude_edges)
    hist_control, _ = np.histogram(control["amp"], amplitude_edges)
    ratio = (np.where(hist_control > 0, hist_reference / np.maximum(hist_control, 1), 0) * len(control)
             / max(len(reference), 1))
    return ratio[np.clip(np.digitize(control["amp"], amplitude_edges) - 1, 0, len(ratio) - 1)]


def found(sample, scut, a2min):
    """Second pulse found: matched-filter maximum above scut, amplitude above a2min, valid delay."""
    return (sample["maxS"] > scut) & (sample["A2"] > a2min) & (sample["tau"] >= 0)


def delay_fractions(sample, delay_edges, scut, a2min, weights=None):
    """Fraction of the candidates with a second pulse in each delay bin, and its error (sum of weights squared)."""
    is_found = found(sample, scut, a2min)
    row_weights = (weights if weights is not None else np.ones(len(sample)))
    hist, _ = np.histogram(sample["tau"][is_found], delay_edges, weights=row_weights[is_found])
    hist_w2, _ = np.histogram(sample["tau"][is_found], delay_edges, weights=row_weights[is_found] ** 2)
    total = row_weights.sum()
    return hist / total, np.sqrt(hist_w2) / total


def plot_max_s(stop, stop_veto, through, through_weights, scut, outdir):
    """c1: maxS distributions (fraction of crystals per bin)."""
    fig, ax = plt.subplots(figsize=(10, 5))
    bins = np.linspace(-5, 40, 91)
    for sample, label, weights in ((stop, "STOP", None), (stop_veto, "STOP, mu+ veto", None),
                                   (through, "THROUGH (reweighted)", through_weights)):
        normalised = (weights if weights is not None else np.ones(len(sample))) / max(len(sample), 1)
        ax.hist(sample["maxS"], bins, weights=normalised, histtype="step", lw=1.3, label=f"{label}: {len(sample)}")
    ax.axvline(scut, color="r", ls="--", label=f"cut maxS > {scut:g}")
    ax.set_yscale("log")
    ax.set_xlabel("max of the matched filter on the residual, S")
    ax.set_ylabel("fraction of crystals / bin")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(outdir / "c1_maxS.png", dpi=110)
    plt.close(fig)


def plot_delays(centres, pairs, outdir):
    """c2: delay of the second pulses, STOP (left: all, right: with the veto) against THROUGH (reweighted)."""
    fig, axs = plt.subplots(1, 2, figsize=(14, 5))
    for ax, (frac_stop, err_stop, frac_through, err_through, label) in zip(axs, pairs):
        ax.errorbar(centres, frac_stop, err_stop, fmt=".", color="k", ms=3, label=label)
        ax.errorbar(centres, frac_through, err_through, fmt=".", color="0.6", ms=3, label="THROUGH (reweighted)")
        ax.set_yscale("log")
        ax.set_xlabel("delay of the second pulse [ns]")
        ax.set_ylabel(f"fraction of crystals / {BIN_WIDTH:g} ns")
        ax.grid(alpha=0.3)
        ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(outdir / "c2_delay.png", dpi=110)
    plt.close(fig)


def injection_efficiency(inject, delay_edges, scut, a2min):
    """Efficiency per delay bin from the injected pulses (injected A2 above a2min), with its binomial error.

    NaN in the bins without injected pulses.
    """
    injected = inject[inject["injA2"] > a2min]
    recovered = ((injected["maxS"] > scut) & (injected["A2"] > a2min)
                 & (np.abs(injected["tau"] - injected["injDelay"]) < INJECTION_MATCH))
    hist_injected, _ = np.histogram(injected["injDelay"], delay_edges)
    hist_recovered, _ = np.histogram(injected["injDelay"][recovered], delay_edges)
    efficiency = np.where(hist_injected > 0, hist_recovered / np.maximum(hist_injected, 1), np.nan)
    efficiency_err = np.sqrt(np.maximum(efficiency * (1 - efficiency), 0) / np.maximum(hist_injected, 1))
    return efficiency, efficiency_err


def plot_efficiency(centres, efficiency, efficiency_err, scut, a2min, outdir):
    """c3: efficiency against the injected delay."""
    fig, ax = plt.subplots(figsize=(10, 4.5))
    ax.errorbar(centres, efficiency, efficiency_err, fmt="o", ms=3, color="k")
    ax.set_ylim(0, 1.05)
    ax.set_xlim(0, 650)
    ax.grid(alpha=0.3)
    ax.set_xlabel("injected delay [ns]")
    ax.set_ylabel("efficiency")
    ax.set_title(f"injected pulses A2 {a2min:g}-400 ADC in THROUGH crystals: found with maxS > {scut:g}, "
                 f"A2 > {a2min:g}, |tau - delay| < 15 ns", fontsize=9)
    fig.tight_layout()
    fig.savefig(outdir / "c3_eff.png", dpi=110)
    plt.close(fig)


def capture_model(t, n0, tau, c):
    """Exponential decay of the captured mu- on a constant: N exp(-t/tau) + C."""
    return n0 * np.exp(-t / tau) + c


def fit_excess(ax, centres, excess, excess_err, in_fit, tlo, thi):
    """Fit N exp(-t/tau) + C to the excess in the fit range, draw it, return the result for the summary.

    Starting values: N = 1.5 x the first point, tau = 80 ns, C = mean of the last 5 points (at least 1e-7).
    """
    x_fit = centres[in_fit]
    y_fit = excess[in_fit]
    p0 = [y_fit[0] * 1.5, TAU_START, max(y_fit[-5:].mean(), 1e-7)]
    params, covariance = curve_fit(capture_model, x_fit, y_fit, p0=p0, sigma=excess_err[in_fit])
    param_errors = np.sqrt(np.diag(covariance))
    chi2 = float(np.sum(((y_fit - capture_model(x_fit, *params)) / excess_err[in_fit]) ** 2))
    ndf = int(in_fit.sum() - 3)
    curve_t = np.linspace(tlo, thi, 300)
    ax.plot(curve_t, capture_model(curve_t, *params), "r-", label="N exp(-t/tau) + C")
    ax.text(0.35, 0.95,
            f"tau = {params[1]:.1f} ± {param_errors[1]:.1f} ns\nN = {params[0]:.2e} ± {param_errors[0]:.1e}"
            f"\nC = {params[2]:.2e} ± {param_errors[2]:.1e}\nchi2/ndf = {chi2:.1f}/{ndf}\nfit {tlo:g}-{thi:g} ns",
            transform=ax.transAxes, va="top", fontsize=9, bbox=dict(fc="w", ec="0.7"))
    return dict(tau=float(params[1]), tau_err=float(param_errors[1]), N=float(params[0]), C=float(params[2]),
                C_err=float(param_errors[2]), chi2=chi2, ndf=ndf)


def plot_excess_fit(centres, pairs, efficiency, tlo, thi, summary, outdir):
    """c4: (STOP - THROUGH) / efficiency in the fit range with the exponential fit, all STOP and with the veto."""
    fig, axs = plt.subplots(1, 2, figsize=(14, 5))
    summary["fit"] = {}
    for ax, (frac_stop, err_stop, frac_through, err_through, label) in zip(axs, pairs):
        excess = (frac_stop - frac_through) / efficiency
        excess_err = np.hypot(err_stop, err_through) / efficiency
        in_fit = (centres >= tlo) & (centres <= thi) & np.isfinite(excess) & (excess_err > 0)
        ax.errorbar(centres[in_fit], excess[in_fit], excess_err[in_fit], fmt="o", ms=3, color="k",
                    label="(STOP - THROUGH) / efficiency")
        try:
            summary["fit"][label] = fit_excess(ax, centres, excess, excess_err, in_fit, tlo, thi)
        except Exception as exn:
            ax.text(0.4, 0.9, f"fit failed: {exn}", transform=ax.transAxes)
        ax.set_title(label, fontsize=10)
        ax.set_xlabel("delay of the second pulse [ns]")
        ax.set_ylabel(f"excess per STOP crystal / {BIN_WIDTH:g} ns (efficiency corrected)")
        ax.grid(alpha=0.3)
        ax.legend(fontsize=8, loc="upper right")
    fig.tight_layout()
    fig.savefig(outdir / "c4_fit.png", dpi=110)
    plt.close(fig)


def main():
    """Samples, reweighting, delay fractions, efficiency, excess fit, summary."""
    args = parse_args()
    outdir = pathlib.Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    scut, a2min, tlo, thi = args.scut, args.a2min, args.tlo, args.thi

    candidates = np.genfromtxt(f"{args.prefix}_capture.csv", delimiter=",", names=True)
    inject = np.genfromtxt(f"{args.prefix}_inject.csv", delimiter=",", names=True)
    stop = candidates[candidates["kind"] == 0]
    through = candidates[candidates["kind"] == 1]
    stop_veto = stop[stop["nMichel"] == 0]
    summary = dict(stop=int(len(stop)), stop_veto=int(len(stop_veto)), through=int(len(through)),
                   inject=int(len(inject)), scut=scut, a2min=a2min, tlo=tlo, thi=thi)

    # amplitude reweighting of the controls: quantile bins of the first-pulse amplitude of STOP + THROUGH
    amplitude_edges = np.quantile(np.r_[stop["amp"], through["amp"]], np.linspace(0, 1, N_AMPLITUDE_BINS + 1))
    amplitude_edges[-1] += 1
    through_weights_stop = amplitude_weights(stop, through, amplitude_edges)
    through_weights_veto = amplitude_weights(stop_veto, through, amplitude_edges)

    plot_max_s(stop, stop_veto, through, through_weights_stop, scut, outdir)

    delay_edges = np.arange(0, MAX_DELAY + BIN_WIDTH, BIN_WIDTH)
    centres = 0.5 * (delay_edges[1:] + delay_edges[:-1])
    frac_stop, err_stop = delay_fractions(stop, delay_edges, scut, a2min)
    frac_veto, err_veto = delay_fractions(stop_veto, delay_edges, scut, a2min)
    frac_through, err_through = delay_fractions(through, delay_edges, scut, a2min, through_weights_stop)
    frac_through_veto, err_through_veto = delay_fractions(through, delay_edges, scut, a2min, through_weights_veto)
    summary["frac_found"] = dict(stop=float(found(stop, scut, a2min).mean()),
                                 stop_veto=float(found(stop_veto, scut, a2min).mean()),
                                 through=float(np.average(found(through, scut, a2min), weights=through_weights_stop)))
    plot_delays(centres, ((frac_stop, err_stop, frac_through, err_through, "STOP (all)"),
                          (frac_veto, err_veto, frac_through_veto, err_through_veto, "STOP with the mu+ veto")),
                outdir)

    efficiency, efficiency_err = injection_efficiency(inject, delay_edges, scut, a2min)
    plot_efficiency(centres, efficiency, efficiency_err, scut, a2min, outdir)

    plot_excess_fit(centres, ((frac_stop, err_stop, frac_through, err_through, "all STOP"),
                              (frac_veto, err_veto, frac_through_veto, err_through_veto, "STOP with the mu+ veto")),
                    efficiency, tlo, thi, summary, outdir)

    json.dump(summary, open(outdir / "summary.json", "w"), indent=1)
    print(json.dumps(summary, indent=1))


if __name__ == "__main__":
    main()
