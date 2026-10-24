"""Time coincidences between CRV coincidence clusters and calorimeter hits (runs with the CRV, KPP extracted position).

For every event with CRV coincidences, every pair (CRV cluster, calorimeter hit with eDep > EMIN) gives
dt = t_calo - t_CRV (raw times, no offsets).

Inputs
  the tables of crv_calo_flatten.C, one or more prefixes:
    <prefix>_crv.csv   CRV coincidence clusters (columns ev, time, sector, ...)
    <prefix>_calo.csv  calorimeter hits (columns ev, time, eDep, ...)

Outputs (in OUTDIR)
  c1_dt_wide.png  dt of all pairs on a wide range, 20 ns bins (to find the coincidence peak)
  c2_dt_peak.png  dt around the peak (+-150 ns, 2 ns bins), one panel per CRV sector type, with a Gaussian + flat fit
                  (parameters on the plot)
  summary.json    runs, number of CRV clusters and pairs, the highest bin of c1, the fit of every sector
  The summary is also printed.

This file also provides the fit helpers used by the other analysis scripts (crv_calo_spatial.py, mu_capture_analysis.py,
neighbour_delayed_analysis.py): curve_fit, gflat, EMIN, SECTOR.

Usage
  python3 crv_calo_coinc.py OUTDIR PREFIX [PREFIX ...]
  e.g. python3 Offline/MC2Reco/scripts/crv_calo_coinc.py out/coinc_125400 flat/r125400
"""
import argparse
import json
import pathlib

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

EMIN = 5.0                                           # MeV, calorimeter hit threshold (Offline energy scale)
SECTOR = {1: "EX", 2: "T1", 3: "T2", 4: "M1-M8"}     # sectorType as in crv_counters_extracted_v04.txt

WIDE_HALF_RANGE = 5000      # ns, half range of the wide dt histogram
WIDE_BIN = 20               # ns, bin of the wide dt histogram
PEAK_HALF_RANGE = 150       # ns, half range of the histogram around the peak
PEAK_BIN = 2                # ns, bin of the histogram around the peak
MIN_PAIRS_FOR_FIT = 30      # pairs in the peak window needed to fit


# ---------------------------------------------------------------------------------------------------------------------
# Fit helpers (shared with the other scripts)
# ---------------------------------------------------------------------------------------------------------------------

def curve_fit(f, x, y, p0, sigma, absolute_sigma=True, nit=200):
    """Weighted least squares fit of y(x) with the model f(x, *p), Levenberg-Marquardt with numerical derivatives
    (numpy only, no scipy in the Mu2e environment).

    sigma are the absolute errors of y (the covariance is not rescaled by chi2/ndf; absolute_sigma is accepted for
    compatibility with scipy but not used). At most nit iterations; stops when the chi2 improves by less than
    1e-8 (relative). Returns (best parameters, covariance matrix).
    """
    params = np.array(p0, float)
    weights = 1.0 / np.asarray(sigma, float)
    damping = 1e-3

    def residuals(q):
        """Weighted residuals (y - f) / sigma."""
        return (y - f(x, *q)) * weights

    def jacobian(q):
        """Numerical derivatives of the weighted model with respect to each parameter."""
        jac = np.empty((len(x), len(q)))
        for j in range(len(q)):
            step_size = 1e-6 * max(abs(q[j]), 1e-3)
            shifted = q.copy()
            shifted[j] += step_size
            jac[:, j] = (f(x, *shifted) - f(x, *q)) * weights / step_size
        return jac

    chi2_current = np.sum(residuals(params) ** 2)
    for _ in range(nit):
        jac = jacobian(params)
        curvature = jac.T @ jac
        gradient = jac.T @ residuals(params)
        step = np.linalg.solve(curvature + damping * np.diag(np.diag(curvature)), gradient)
        chi2_trial = np.sum(residuals(params + step) ** 2)
        if chi2_trial < chi2_current:
            converged = abs(chi2_current - chi2_trial) < 1e-8 * max(chi2_current, 1)
            params = params + step
            damping = damping / 10
            chi2_current = chi2_trial
            if converged:
                break
        else:
            damping *= 10
    jac = jacobian(params)
    return params, np.linalg.inv(jac.T @ jac)


def gflat(x, a, mu, s, b):
    """Gaussian (amplitude a, mean mu, sigma s) on a flat background b."""
    return a * np.exp(-0.5 * ((x - mu) / s) ** 2) + b


# ---------------------------------------------------------------------------------------------------------------------
# Steps
# ---------------------------------------------------------------------------------------------------------------------

def parse_args():
    """Command line: output directory and one or more table prefixes."""
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("outdir", help="output directory (created if missing)")
    parser.add_argument("prefixes", nargs="+", metavar="PREFIX",
                        help="prefix of the crv_calo_flatten.C tables (<prefix>_crv.csv, <prefix>_calo.csv)")
    return parser.parse_args()


def load(prefix):
    """Read the CRV-cluster and calorimeter-hit tables of one prefix."""
    crv = np.genfromtxt(f"{prefix}_crv.csv", delimiter=",", names=True)
    calo = np.genfromtxt(f"{prefix}_calo.csv", delimiter=",", names=True)
    return np.atleast_1d(crv), np.atleast_1d(calo)


def pairs(crv, calo):
    """All (CRV cluster, calorimeter hit above EMIN) pairs of the same event.

    Returns dt = t_calo - t_CRV, the CRV sector type of the pair and the index of the CRV cluster.
    """
    calo = calo[calo["eDep"] > EMIN]
    order = np.argsort(calo["ev"], kind="stable")
    calo = calo[order]
    calo_event = calo["ev"]
    # hits of the event of each CRV cluster: calo[first[i]:last[i]]
    first = np.searchsorted(calo_event, crv["ev"], "left")
    last = np.searchsorted(calo_event, crv["ev"], "right")
    dt_list, sector_list, index_list = [], [], []
    for i in range(len(crv)):
        if last[i] > first[i]:
            dt = calo["time"][first[i]:last[i]] - crv["time"][i]
            dt_list.append(dt)
            sector_list.append(np.full(len(dt), crv["sector"][i]))
            index_list.append(np.full(len(dt), i))
    if dt_list:
        return np.concatenate(dt_list), np.concatenate(sector_list).astype(int), np.concatenate(index_list)
    return (np.array([]),) * 3


def text_box(ax, x, y, text):
    """White box with the fit / count results in axes coordinates."""
    ax.text(x, y, text, transform=ax.transAxes, va="top", fontsize=9, bbox=dict(fc="w", ec="0.7"))


def plot_wide(dt, runs, summary, outdir):
    """c1: dt of all pairs on +-5 us; the highest bin is the coincidence peak used by c2."""
    fig, ax = plt.subplots(figsize=(11, 5))
    bin_edges = np.arange(-WIDE_HALF_RANGE, WIDE_HALF_RANGE + WIDE_BIN, WIDE_BIN)
    counts, edges = np.histogram(dt, bin_edges)
    centres = 0.5 * (edges[1:] + edges[:-1])
    ax.errorbar(centres, counts, yerr=np.sqrt(counts), fmt=".", ms=3, lw=0.6, color="k")
    peak = centres[np.argmax(counts)]
    summary["wide_peak_ns"] = float(peak)
    ax.set_xlabel("t_calo - t_CRV [ns]  (calorimeter hit time - CRV coincidence time, raw, no offsets)")
    ax.set_ylabel("pairs / 20 ns")
    ax.set_title(f"CRV-calorimeter time difference, all pairs (calo hits > {EMIN:g} MeV), runs {', '.join(runs)}",
                 fontsize=10)
    text_box(ax, 0.02, 0.92, f"pairs: {len(dt)}\nhighest bin at {peak:.0f} ns")
    ax.set_yscale("log")
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(outdir / "c1_dt_wide.png", dpi=110)
    plt.close(fig)
    return peak


def fit_sector_peak(ax, centres, counts, sector_name, summary):
    """Gaussian + flat fit of the peak of one sector: draws the curve, stores the parameters, returns the text lines.

    Starting values: amplitude = max - median, mean at the highest bin, sigma 5 ns, flat = median (at least 0.1).
    Errors sqrt(n), at least 1. The signal is the Gaussian integral divided by the 2 ns bin width.
    """
    p0 = [counts.max() - np.median(counts), centres[np.argmax(counts)], 5, max(np.median(counts), 0.1)]
    errors = np.sqrt(np.maximum(counts, 1))
    params, covariance = curve_fit(gflat, centres, counts, p0=p0, sigma=errors, absolute_sigma=True)
    param_errors = np.sqrt(np.diag(covariance))
    chi2 = float(np.sum(((counts - gflat(centres, *params)) / np.sqrt(np.maximum(counts, 1))) ** 2))
    ndf = len(centres) - 4
    curve_x = np.linspace(centres[0], centres[-1], 600)
    ax.plot(curve_x, gflat(curve_x, *params), "r-", lw=1.2, label="Gaussian + flat")
    signal_pairs = params[0] * abs(params[2]) * np.sqrt(2 * np.pi) / 2.0
    text = (f"\nmean = {params[1]:.2f} ± {param_errors[1]:.2f} ns"
            f"\nsigma = {abs(params[2]):.2f} ± {param_errors[2]:.2f} ns"
            f"\nflat = {params[3]:.2f} ± {param_errors[3]:.2f} /2 ns"
            f"\nsignal ≈ {signal_pairs:.0f} pairs"
            f"\nchi2/ndf = {chi2:.0f}/{ndf}")
    summary["peak_fit"][sector_name] = dict(mean=float(params[1]), mean_err=float(param_errors[1]),
                                            sigma=float(abs(params[2])), sigma_err=float(param_errors[2]),
                                            flat=float(params[3]), signal_pairs=float(signal_pairs),
                                            chi2=chi2, ndf=ndf)
    return text


def plot_peak_per_sector(dt, sector, peak, runs, summary, outdir):
    """c2: dt within +-150 ns of the peak, one panel per CRV sector type, fitted when there are > 30 pairs."""
    fig, axs = plt.subplots(2, 2, figsize=(13, 8.5))
    summary["peak_fit"] = {}
    for ax, sector_type in zip(axs.flat, (1, 2, 3, 4)):
        in_sector = sector == sector_type
        bin_edges = np.arange(peak - PEAK_HALF_RANGE, peak + PEAK_HALF_RANGE + PEAK_BIN, PEAK_BIN)
        counts, edges = np.histogram(dt[in_sector], bin_edges)
        centres = 0.5 * (edges[1:] + edges[:-1])
        ax.errorbar(centres, counts, yerr=np.sqrt(np.maximum(counts, 1)), fmt=".", ms=3, lw=0.6, color="k",
                    label="data")
        text = f"pairs in window: {int(counts.sum())}"
        if counts.sum() > MIN_PAIRS_FOR_FIT:
            try:
                text += fit_sector_peak(ax, centres, counts, SECTOR[sector_type], summary)
            except Exception as ex:
                text += f"\nfit failed: {ex}"
        ax.text(0.02, 0.97, text, transform=ax.transAxes, va="top", fontsize=8.5, bbox=dict(fc="w", ec="0.7"))
        ax.set_title(f"CRV sector {SECTOR[sector_type]}", fontsize=10)
        ax.set_xlabel("t_calo - t_CRV [ns]")
        ax.set_ylabel("pairs / 2 ns")
        ax.grid(alpha=0.3)
        ax.legend(loc="upper right", fontsize=8)
    fig.suptitle(f"CRV-calorimeter coincidence peak per CRV sector (calo hits > {EMIN:g} MeV), runs {', '.join(runs)}",
                 fontsize=11)
    fig.tight_layout()
    fig.savefig(outdir / "c2_dt_peak.png", dpi=110)
    plt.close(fig)


def main():
    """Pairs of all prefixes, wide plot, per-sector peak fits, summary."""
    args = parse_args()
    outdir = pathlib.Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    runs = [pathlib.Path(p).name for p in args.prefixes]

    crv_tables, calo_tables = zip(*[load(p) for p in args.prefixes])
    dt_per_run, sector_per_run = [], []
    for crv, calo in zip(crv_tables, calo_tables):
        dt_run, sector_run, _ = pairs(crv, calo)
        dt_per_run.append(dt_run)
        sector_per_run.append(sector_run)
    dt = np.concatenate(dt_per_run)
    sector = np.concatenate(sector_per_run)

    summary = dict(runs=runs, emin_MeV=EMIN, n_crv_clusters=int(sum(len(c) for c in crv_tables)),
                   n_pairs=int(len(dt)))
    peak = plot_wide(dt, runs, summary, outdir)
    plot_peak_per_sector(dt, sector, peak, runs, summary, outdir)

    json.dump(summary, open(outdir / "summary.json", "w"), indent=1)
    print(json.dumps(summary, indent=1))


if __name__ == "__main__":
    main()
