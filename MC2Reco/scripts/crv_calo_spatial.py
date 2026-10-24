"""Spatial check of the CRV-calorimeter coincidences: does the track seen in the calorimeter point to the CRV cluster?

Selection
  Coincidences are selected in time only: CRV cluster + calorimeter hits (eDep > EMIN) with |dt - mean| < NSIG sigma of
  the peak of crv_calo_coinc.py (its summary.json, per CRV sector). One hit per crystal is kept; at least MINCRY
  crystals in one disk and a vertical extent of at least 30 mm are required.
Calorimeter track
  A straight line through the centres of the hit crystals (x, y of crystalPos, z = crystal front face + 100 mm, the
  middle of the crystal), fitted as x(y) and z(y) by least squares (energy-weighted), and extrapolated to the height
  of the CRV cluster.
  - z: needs both disks hit (two z values); compared with the CRV z of EX and T2 (bars stacked along z).
  - x: one disk is enough; compared with the CRV x of T2 when the position along the bar was computed from the two
    sides (posCalc).
Coordinates: Mu2e detector frame (mm), origin at the tracker centre: crystal front face of disk 0 at z = 1744.9.
Each run is passed with its own coincidence summary (the CRV-calo time offset can change from run to run).

Inputs
  for each run, PREFIX=COINCDIR: the tables of crv_calo_flatten.C (<PREFIX>_crv.csv, <PREFIX>_calo.csv) and the
  directory with the summary.json of crv_calo_coinc.py for that run.

Outputs (in OUTDIR)
  s1_dz.png            z residual z_CRV - z_predicted for EX and T2, two-disk tracks, with a Gaussian + flat fit
  s2_z_scatter.png     z_CRV against z predicted, EX and T2
  s3_dx_T2.png         x residual for T2 with the position along the bar, with fit, and x_CRV against x predicted
  summary_spatial.json runs, cuts, track counts, fit results, the time-peak fits used; printed without the peak fits

Usage
  python3 crv_calo_spatial.py OUTDIR PREFIX=COINCDIR [PREFIX=COINCDIR ...]
  e.g. python3 Offline/MC2Reco/scripts/crv_calo_spatial.py out/spatial flat/r125400=out/coinc_125400
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
from crv_calo_coinc import curve_fit, gflat, EMIN, SECTOR  # noqa: E402

NSIG = 3.0                   # time window around the coincidence peak, in sigma of the peak fit
HALF_LEN = 100.0             # mm, half crystal length: track point at the middle of the crystal
MINCRY = 2                   # crystals per disk
MIN_VERTICAL_EXTENT = 30     # mm, minimum spread in y of the hit crystals (to fit x(y), z(y))
SECTOR_EX = 1                # sectorType of EX
SECTOR_T2 = 3                # sectorType of T2

TRACK_DTYPE = [("sector", int), ("xc", float), ("zc", float), ("posCalc", int), ("xp", float), ("zp", float),
               ("ncry", int), ("two", int), ("dxdy", float)]


def parse_args():
    """Command line: output directory and one PREFIX=COINCDIR per run."""
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("outdir", help="output directory (created if missing)")
    parser.add_argument("runs", nargs="*", metavar="PREFIX=COINCDIR",
                        help="table prefix of crv_calo_flatten.C and directory of the crv_calo_coinc.py summary")
    return parser.parse_args()


def weighted_line_fit(design, weight_matrix, values):
    """Weighted least squares of values = a + b y (normal equations); returns [a, b]."""
    return np.linalg.lstsq(design.T @ weight_matrix @ design, design.T @ weight_matrix @ values, rcond=None)[0]


def tracks(crv, calo, peak_fit):
    """Calorimeter track of every CRV cluster with in-time hits, extrapolated to the CRV height.

    One row per accepted CRV cluster: sector, CRV x and z, posCalc, predicted x and z (z is NaN without two disks),
    number of crystals, two-disk flag, slope dx/dy.
    """
    calo = calo[calo["eDep"] > EMIN]
    calo = calo[np.argsort(calo["ev"], kind="stable")]
    calo_event = calo["ev"]
    # hits of the event of each CRV cluster: calo[first[i]:last[i]]
    first = np.searchsorted(calo_event, crv["ev"], "left")
    last = np.searchsorted(calo_event, crv["ev"], "right")
    rows = []
    for i in range(len(crv)):
        sector_name = SECTOR.get(int(crv["sector"][i]))
        if sector_name not in peak_fit or last[i] <= first[i]:
            continue
        fit = peak_fit[sector_name]
        hits = calo[first[i]:last[i]]
        dt = hits["time"] - crv["time"][i]
        hits = hits[np.abs(dt - fit["mean"]) < NSIG * fit["sigma"]]
        if len(hits) < MINCRY:
            continue
        _, first_of_crystal = np.unique(hits["crystal"], return_index=True)
        hits = hits[first_of_crystal]                                       # one entry per crystal
        z = hits["z"] + HALF_LEN
        disks = np.unique(np.round(hits["z"]))
        crystals_per_disk = [np.sum(np.round(hits["z"]) == d) for d in disks]
        if max(crystals_per_disk) < MINCRY:
            continue
        y, energy = hits["y"], hits["eDep"]
        if np.ptp(y) < MIN_VERTICAL_EXTENT:                                 # needs some vertical extent
            continue
        design = np.c_[np.ones(len(y)), y]
        weight_matrix = np.diag(energy)
        x_line = weighted_line_fit(design, weight_matrix, hits["x"])
        y_crv = crv["y"][i]
        x_predicted = x_line[0] + x_line[1] * y_crv
        two_disks = len(disks) == 2 and min(crystals_per_disk) >= MINCRY
        z_predicted = np.nan
        if two_disks:
            z_line = weighted_line_fit(design, weight_matrix, z)
            z_predicted = z_line[0] + z_line[1] * y_crv
        rows.append((int(crv["sector"][i]), crv["x"][i], crv["z"][i], int(crv["posCalc"][i]), x_predicted,
                     z_predicted, len(hits), int(two_disks), x_line[1]))
    return np.array(rows, dtype=TRACK_DTYPE)


def fitplot(ax, residuals, lo, hi, bin_width, label, summary, key):
    """Histogram of residuals in [lo, hi] with a Gaussian + flat fit (> 30 entries); stores the fit as summary[key].

    Starting values: amplitude = max, mean at the highest bin, sigma 100 mm, flat = median (at least 0.5).
    """
    bin_edges = np.arange(lo, hi + bin_width, bin_width)
    counts, edges = np.histogram(residuals, bin_edges)
    centres = 0.5 * (edges[1:] + edges[:-1])
    ax.errorbar(centres, counts, yerr=np.sqrt(np.maximum(counts, 1)), fmt=".", ms=3, lw=0.6, color="k", label="data")
    in_range = residuals[(residuals > lo) & (residuals < hi)]
    text = (f"entries: {len(residuals)}\nmedian = {np.median(residuals):.0f} mm"
            f"\nrms (in range) = {np.std(in_range):.0f} mm")
    if counts.sum() > 30:
        try:
            p0 = [counts.max(), centres[np.argmax(counts)], 100, max(np.median(counts), 0.5)]
            params, covariance = curve_fit(gflat, centres, counts, p0=p0, sigma=np.sqrt(np.maximum(counts, 1)))
            param_errors = np.sqrt(np.diag(covariance))
            chi2 = float(np.sum(((counts - gflat(centres, *params)) / np.sqrt(np.maximum(counts, 1))) ** 2))
            ndf = len(centres) - 4
            curve_x = np.linspace(lo, hi, 600)
            ax.plot(curve_x, gflat(curve_x, *params), "r-", lw=1.2, label="Gaussian + flat")
            text += (f"\nfit mean = {params[1]:.0f} ± {param_errors[1]:.0f} mm"
                     f"\nfit sigma = {abs(params[2]):.0f} ± {param_errors[2]:.0f} mm"
                     f"\nflat = {params[3]:.1f} ± {param_errors[3]:.1f}"
                     f"\nchi2/ndf = {chi2:.0f}/{ndf}")
            summary[key] = dict(n=int(len(residuals)), median=float(np.median(residuals)), mean=float(params[1]),
                                mean_err=float(param_errors[1]), sigma=float(abs(params[2])),
                                sigma_err=float(param_errors[2]), chi2=chi2, ndf=ndf)
        except Exception as ex:
            text += f"\nfit failed: {ex}"
    ax.text(0.02, 0.97, text, transform=ax.transAxes, va="top", fontsize=8.5, bbox=dict(fc="w", ec="0.7"))
    ax.set_xlabel(label)
    ax.set_ylabel(f"coincidences / {bin_width:g} mm")
    ax.grid(alpha=0.3)
    ax.legend(loc="upper right", fontsize=8)


def load_runs(run_specs):
    """Tracks of every run (each with its own time-peak fit); returns (run names, all tracks, fits used per run)."""
    specs = [spec.split("=") for spec in run_specs]
    runs = [pathlib.Path(prefix).name for prefix, _ in specs]
    track_tables = []
    fits = {}
    for prefix, coinc_dir in specs:
        peak_fit = json.load(open(pathlib.Path(coinc_dir) / "summary.json"))["peak_fit"]
        fits[pathlib.Path(prefix).name] = peak_fit
        crv = np.atleast_1d(np.genfromtxt(f"{prefix}_crv.csv", delimiter=",", names=True))
        calo = np.atleast_1d(np.genfromtxt(f"{prefix}_calo.csv", delimiter=",", names=True))
        track_tables.append(tracks(crv, calo, peak_fit))
    return runs, np.concatenate(track_tables), fits


def plot_z_residuals(trk, title, summary, outdir):
    """s1: z residual of the two-disk tracks for EX and T2."""
    fig, axs = plt.subplots(1, 2, figsize=(14, 5.5))
    for ax, sector_type in zip(axs, (SECTOR_EX, SECTOR_T2)):
        selected = (trk["sector"] == sector_type) & (trk["two"] == 1)
        residuals = trk["zc"][selected] - trk["zp"][selected]
        fitplot(ax, residuals, -2000, 2000, 50, "z_CRV - z_predicted [mm]", summary, f"dz_{SECTOR[sector_type]}")
        ax.set_title(f"CRV {SECTOR[sector_type]}: z residual, calorimeter tracks through both disks", fontsize=10)
    fig.suptitle(title, fontsize=10)
    fig.tight_layout()
    fig.savefig(outdir / "s1_dz.png", dpi=110)
    plt.close(fig)


def plot_z_scatter(trk, title, outdir):
    """s2: z of the CRV cluster against the z predicted by the two-disk tracks, EX and T2."""
    fig, axs = plt.subplots(1, 2, figsize=(13, 6))
    for ax, sector_type in zip(axs, (SECTOR_EX, SECTOR_T2)):
        selected = (trk["sector"] == sector_type) & (trk["two"] == 1)
        ax.plot(trk["zp"][selected], trk["zc"][selected], ".", ms=2, alpha=0.4, color="k")
        limits = (-4000, 2000)
        ax.plot(limits, limits, "r--", lw=1, label="z_CRV = z_predicted")
        ax.set_xlim(*limits)
        ax.set_ylim(*limits)
        ax.set_aspect("equal")
        ax.grid(alpha=0.3)
        ax.legend(fontsize=8)
        ax.set_xlabel("z predicted at the CRV height from the calorimeter track [mm]")
        ax.set_ylabel(f"z of the CRV {SECTOR[sector_type]} cluster [mm]")
        ax.set_title(f"CRV {SECTOR[sector_type]} ({selected.sum()} two-disk tracks)", fontsize=10)
    fig.suptitle(title + "\ndetector frame: crystal front face of disk 0 at z = 1745 mm", fontsize=10)
    fig.tight_layout()
    fig.savefig(outdir / "s2_z_scatter.png", dpi=110)
    plt.close(fig)


def plot_x_t2(trk, title, summary, outdir):
    """s3: x residual for T2 with the position along the bar from the two sides (one disk is enough)."""
    fig, axs = plt.subplots(1, 2, figsize=(14, 5.5))
    selected = (trk["sector"] == SECTOR_T2) & (trk["posCalc"] == 1)
    residuals = trk["xc"][selected] - trk["xp"][selected]
    fitplot(axs[0], residuals, -2000, 2000, 50, "x_CRV - x_predicted [mm]", summary, "dx_T2")
    axs[0].set_title("CRV T2 (position along the bar from the two sides): x residual", fontsize=10)
    axs[1].plot(trk["xp"][selected], trk["xc"][selected], ".", ms=2, alpha=0.4, color="k")
    limits = (-2500, 2500)
    axs[1].plot(limits, limits, "r--", lw=1, label="x_CRV = x_predicted")
    axs[1].set_xlim(*limits)
    axs[1].set_ylim(*limits)
    axs[1].set_aspect("equal")
    axs[1].grid(alpha=0.3)
    axs[1].legend(fontsize=8)
    axs[1].set_xlabel("x predicted at the CRV height [mm]")
    axs[1].set_ylabel("x of the CRV T2 cluster [mm]")
    fig.suptitle(title, fontsize=10)
    fig.tight_layout()
    fig.savefig(outdir / "s3_dx_T2.png", dpi=110)
    plt.close(fig)


def main():
    """Tracks of all runs, z and x comparisons, summary."""
    args = parse_args()
    outdir = pathlib.Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    runs, trk, fits = load_runs(args.runs)
    summary = dict(runs=runs, nsig=NSIG, emin_MeV=EMIN, n_tracks=int(len(trk)), n_two_disk=int(trk["two"].sum()),
                   peak_fit_used=fits)
    title = f"runs {', '.join(runs)}; time coincidence |dt-mean| < {NSIG:g} sigma, calo hits > {EMIN:g} MeV"

    plot_z_residuals(trk, title, summary, outdir)
    plot_z_scatter(trk, title, outdir)
    plot_x_t2(trk, title, summary, outdir)

    json.dump(summary, open(outdir / "summary_spatial.json", "w"), indent=1)
    print(json.dumps({k: v for k, v in summary.items() if k != "peak_fit_used"}, indent=1))


if __name__ == "__main__":
    main()
