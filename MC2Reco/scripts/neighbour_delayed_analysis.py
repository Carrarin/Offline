"""mu- capture with the delayed deposits in the crystals around the stopping crystal (tree "delayed" of MC2MuCapture).

Deposits
  One entry per neighbouring crystal deposit (SiPM 0 of the crystal only, so that L and R are not counted twice),
  ring 1 = within 80 mm of the candidate crystal (not on the track), ring 2 = 250-400 mm (accidentals).
  SiPMs without a time offset from the tracks (calibration stats file, timeFromTracks != 1) are left out: their
  deposits can be displaced by board-level offsets. Optional energy threshold EMIN.
Samples
  all STOP, STOP with a mu+ decay tag (nMichel > 0, control), THROUGH; rates per candidate per ns. The flat level
  before the muon (dt in [-1000, -50) ns) is shown as a dotted line.
Fit
  The ring-1 spectrum of all STOP is fitted with A exp(-t/tauC) + B exp(-t/2197 ns) + C (capture, mu+ decay, flat)
  from several start times (20, 50, 100 ns): the stability of tauC with the start tells a real fast component from the
  tail of the prompt deposits. The same fit without the capture term gives chi2_no_capture.

Inputs
  DELAYED.csv      deposits (delayed_to_csv.C): kind, nMichel, sipm, e, ring, dt
  CAPTURE.csv      candidates (capture_to_csv.C, <prefix>_capture.csv): kind, nMichel
  CALIB_STATS.csv  calibration stats of the last calibration pass: roid, timeFromTracks
Outputs (in OUTDIR)
  n1_dt.png     rate of the deposits against dt, ring 1 and ring 2, for the three samples
  n2_fit.png    ring-1 rate of all STOP (and THROUGH) with the fits from the three start times
  summary.json  threshold, SiPMs used, candidates per sample, fit results (also printed)

Usage
  python3 neighbour_delayed_analysis.py DELAYED.csv CAPTURE.csv CALIB_STATS.csv OUTDIR [EMIN_MeV=0]
  e.g. python3 Offline/MC2Reco/scripts/neighbour_delayed_analysis.py out/delayed.csv out/physics_capture.csv \\
           out/calib_pass3_stats.csv out/neighbour
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

TAU_MU_PLUS = 2197.             # ns, free mu+ lifetime
FIT_STARTS = ((20, "r"), (50, "m"), (100, "b"))   # start of the fit [ns] and colour of its curve
# variable dt bins [ns]: coarse before the muon, 10 ns up to 500 ns, then wider
DT_EDGES = np.r_[np.arange(-1000, 0, 50), np.arange(0, 500, 10), np.arange(500, 2000, 100),
                 np.arange(2000, 10001, 500)]
SAMPLE_COLOURS = ("k", "C3", "0.5")


def parse_args():
    """Command line: the three tables, the output directory and the optional energy threshold."""
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("delayed", help="table of the delayed deposits (delayed_to_csv.C)")
    parser.add_argument("capture", help="table of the candidates (<prefix>_capture.csv of capture_to_csv.C)")
    parser.add_argument("calib_stats", help="calibration stats of the SiPMs (roid, timeFromTracks)")
    parser.add_argument("outdir", help="output directory (created if missing)")
    parser.add_argument("emin", nargs="?", type=float, default=0., help="deposit energy threshold [MeV] (default 0)")
    return parser.parse_args()


def select_deposits(deposits, calib_stats, emin):
    """SiPM 0 of each crystal, above emin, on SiPMs with a time offset from the tracks; returns (deposits, SiPMs)."""
    good_sipms = set(calib_stats["roid"][calib_stats["timeFromTracks"] == 1].astype(int))
    selected = deposits[(deposits["sipm"] % 2 == 0) & (deposits["e"] > emin)
                        & np.isin(deposits["sipm"].astype(int), np.array(sorted(good_sipms)))]
    return selected, good_sipms


def sample_masks(candidates, deposits):
    """The three samples: (candidate mask, deposit mask) by name, in plotting order."""
    return {"STOP (all)": (candidates["kind"] == 0, deposits["kind"] == 0),
            "STOP with a mu+ decay tag (control)": ((candidates["kind"] == 0) & (candidates["nMichel"] > 0),
                                                    (deposits["kind"] == 0) & (deposits["nMichel"] > 0)),
            "THROUGH": (candidates["kind"] == 1, deposits["kind"] == 1)}


def plot_rates(deposits, samples, centres, bin_widths, emin, summary, outdir):
    """n1: rate of deposits per candidate per ns against dt, ring 1 (left) and ring 2 (right).

    Returns {(sample, ring): (rate, error, flat level before the muon, counts, candidates)}.
    """
    rates = {}
    fig, axs = plt.subplots(1, 2, figsize=(15, 5.5))
    for (name, (candidate_mask, deposit_mask)), colour in zip(samples.items(), SAMPLE_COLOURS):
        n_candidates = int(candidate_mask.sum())
        for ring, ax in ((1, axs[0]), (2, axs[1])):
            dt = deposits[deposit_mask & (deposits["ring"] == ring)]["dt"]
            counts, _ = np.histogram(dt, DT_EDGES)
            rate = counts / bin_widths / max(n_candidates, 1)
            rate_err = np.sqrt(np.maximum(counts, 1)) / bin_widths / max(n_candidates, 1)
            # flat level from the 950 ns before the muon
            flat = ((dt >= -1000) & (dt < -50)).sum() / 950. / max(n_candidates, 1)
            rates[(name, ring)] = (rate, rate_err, flat, counts, n_candidates)
            non_empty = counts > 0
            ax.errorbar(centres[non_empty], rate[non_empty], rate_err[non_empty], fmt=".", ms=3, color=colour,
                        label=f"{name}: {n_candidates} candidates, before the muon {flat:.1e} /ns")
            ax.axhline(flat, color=colour, ls=":", lw=0.8)
        summary[name] = dict(candidates=n_candidates)
    for ax, title in zip(axs, ("ring 1 (within 80 mm, not on the track)", "ring 2 (250-400 mm, accidentals)")):
        ax.set_xscale("symlog", linthresh=100)
        ax.set_yscale("log")
        ax.set_ylim(1e-9, 0.2)
        ax.set_xlabel("dt = t(deposit) - t(candidate crystal) [ns]")
        ax.set_ylabel("deposits per candidate per ns")
        ax.grid(alpha=0.3)
        ax.legend(fontsize=8)
        ax.set_title(title, fontsize=10)
    fig.suptitle("Delayed deposits around the candidate crystals (SiPM 0 of each crystal, calibrated SiPMs only, "
                 f"E > {emin:g} MeV)", fontsize=10.5)
    fig.tight_layout()
    fig.savefig(outdir / "n1_dt.png", dpi=110)
    plt.close(fig)
    return rates


def capture_model(t, A, tauC, B, Cc):
    """Capture (A, tauC) + mu+ decay (B, 2197 ns) + flat (Cc)."""
    return A * np.exp(-t / tauC) + B * np.exp(-t / TAU_MU_PLUS) + Cc


def no_capture_model(t, B, Cc):
    """mu+ decay (B, 2197 ns) + flat (Cc), without the capture term."""
    return B * np.exp(-t / TAU_MU_PLUS) + Cc


def fit_from(ax, centres, rate, rate_err, counts, fit_start, colour, summary):
    """Fit the ring-1 rate from fit_start with and without the capture term; draw it, store it, return its text line.

    Starting values: A = 3e-5, tauC = 80 ns, B = 1e-5, C = 1e-8 per candidate per ns.
    """
    in_fit = (centres >= fit_start) & (counts > 0)
    params, covariance = curve_fit(capture_model, centres[in_fit], rate[in_fit], p0=[3e-5, 80., 1e-5, 1e-8],
                                   sigma=rate_err[in_fit])
    param_errors = np.sqrt(np.diag(covariance))
    chi2 = float(np.sum(((rate[in_fit] - capture_model(centres[in_fit], *params)) / rate_err[in_fit]) ** 2))
    params_no_capture, _ = curve_fit(no_capture_model, centres[in_fit], rate[in_fit], p0=[1e-5, 1e-8],
                                     sigma=rate_err[in_fit])
    chi2_no_capture = float(np.sum(((rate[in_fit] - no_capture_model(centres[in_fit], *params_no_capture))
                                    / rate_err[in_fit]) ** 2))
    curve_t = np.linspace(fit_start, 10000, 2000)
    ax.plot(curve_t, capture_model(curve_t, *params), "-", color=colour, lw=1, label=f"fit from {fit_start} ns")
    summary["fits"][str(fit_start)] = dict(tauC=float(params[1]), tauC_err=float(param_errors[1]),
                                           A=float(params[0]), A_err=float(param_errors[0]),
                                           B=float(params[2]), B_err=float(param_errors[2]),
                                           C=float(params[3]), chi2=chi2, ndf=int(in_fit.sum() - 4),
                                           chi2_no_capture=chi2_no_capture)
    return (f"from {fit_start:3d} ns: tauC = {params[1]:.0f} ± {param_errors[1]:.0f} ns, "
            f"A = {params[0]:.1e} ± {param_errors[0]:.0e}, chi2/ndf {chi2:.0f}/{in_fit.sum()-4} "
            f"(no capture {chi2_no_capture:.0f})")


def plot_fits(rates, centres, summary, outdir):
    """n2: ring-1 rate of all STOP (and THROUGH for comparison) from 20 ns, with the fits from each start time."""
    rate, rate_err, _, counts, _ = rates[("STOP (all)", 1)]
    fig, ax = plt.subplots(figsize=(10.5, 5.8))
    shown = (centres >= 20) & (counts > 0)
    ax.errorbar(centres[shown], rate[shown], rate_err[shown], fmt="o", ms=3, color="k", label="STOP (all), ring 1")
    rate_through, rate_err_through, _, _, _ = rates[("THROUGH", 1)]
    ax.errorbar(centres[shown], rate_through[shown], rate_err_through[shown], fmt=".", ms=3, color="0.6",
                label="THROUGH, ring 1")
    summary["fits"] = {}
    lines = []
    for fit_start, colour in FIT_STARTS:
        try:
            lines.append(fit_from(ax, centres, rate, rate_err, counts, fit_start, colour, summary))
        except Exception as exn:
            lines.append(f"from {fit_start} ns: fit failed ({exn})")
    ax.text(0.02, 0.03, "A exp(-t/tauC) + B exp(-t/2197 ns) + C\n" + "\n".join(lines), transform=ax.transAxes,
            fontsize=7.5, va="bottom", bbox=dict(fc="w", ec="0.7"))
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("dt [ns]")
    ax.set_ylabel("deposits per candidate per ns")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8, loc="upper right")
    fig.tight_layout()
    fig.savefig(outdir / "n2_fit.png", dpi=110)
    plt.close(fig)


def main():
    """Selection of the deposits, rates per sample and ring, fits of the ring-1 STOP rate, summary."""
    args = parse_args()
    deposits = np.genfromtxt(args.delayed, delimiter=",", names=True)
    candidates = np.genfromtxt(args.capture, delimiter=",", names=True)
    calib_stats = np.genfromtxt(args.calib_stats, delimiter=",", names=True)
    outdir = pathlib.Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    emin = args.emin

    deposits, good_sipms = select_deposits(deposits, calib_stats, emin)
    centres = 0.5 * (DT_EDGES[1:] + DT_EDGES[:-1])
    bin_widths = np.diff(DT_EDGES)
    samples = sample_masks(candidates, deposits)
    summary = dict(emin=emin, sipms_used=len(good_sipms))

    rates = plot_rates(deposits, samples, centres, bin_widths, emin, summary, outdir)
    plot_fits(rates, centres, summary, outdir)

    json.dump(summary, open(outdir / "summary.json", "w"), indent=1)
    print(json.dumps(summary, indent=1))


if __name__ == "__main__":
    main()
