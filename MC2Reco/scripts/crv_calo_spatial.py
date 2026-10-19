"""Spatial check of the CRV-calorimeter coincidences: does the track seen in the calorimeter point to the CRV cluster?

Coincidences are selected in time only: CRV cluster + calorimeter hits (eDep > EMIN) with |dt - mean| < NSIG sigma of the peak of
crv_calo_coinc.py (summary.json, per sector). The calorimeter track is a straight line through the centres of the hit crystals
(x, y of crystalPos, z = crystal front face + 100 mm, the middle of the crystal), fitted as x(y) and z(y) by least squares
(energy-weighted), and extrapolated to the height of the CRV cluster.
  - z: needs both disks hit (two z values); compared with the CRV z of EX and T2 (bars stacked along z).
  - x: one disk is enough; compared with the CRV x of T2 when the position along the bar was computed from the two sides (posCalc).
Coordinates: Mu2e detector frame (mm), origin at the tracker centre: crystal front face of disk 0 at z = 1744.9.
Figures: s1 residuals in z, s2 z_CRV vs z predicted, s3 residuals in x (T2), with fits (parameters on the plots). summary_spatial.json.
Each run is passed with its own coincidence summary (the CRV-calo time offset can change from run to run).
usage: python3 crv_calo_spatial.py OUTDIR PREFIX=COINCDIR [PREFIX=COINCDIR ...]
"""
import sys, json, pathlib, numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
sys.path.insert(0, str(pathlib.Path(__file__).parent)); from crv_calo_coinc import curve_fit, gflat, EMIN, SECTOR

NSIG = 3.0
HALF_LEN = 100.0     # mm, half crystal length: track point at the middle of the crystal
MINCRY = 2           # crystals per disk

def tracks(c, h, fit):
    h = h[h["eDep"] > EMIN]; h = h[np.argsort(h["ev"], kind="stable")]; ev = h["ev"]
    lo = np.searchsorted(ev, c["ev"], "left"); hi = np.searchsorted(ev, c["ev"], "right"); rows = []
    for i in range(len(c)):
        name = SECTOR.get(int(c["sector"][i]))
        if name not in fit or hi[i] <= lo[i]: continue
        f = fit[name]; s = h[lo[i]:hi[i]]; dt = s["time"] - c["time"][i]; s = s[np.abs(dt - f["mean"]) < NSIG * f["sigma"]]
        if len(s) < MINCRY: continue
        _, k = np.unique(s["crystal"], return_index=True); s = s[k]          # one entry per crystal
        z = s["z"] + HALF_LEN; disks = np.unique(np.round(s["z"]))
        nd = [np.sum(np.round(s["z"]) == d) for d in disks]
        if max(nd) < MINCRY: continue
        y, w = s["y"], s["eDep"]
        if np.ptp(y) < 30: continue                                             # needs some vertical extent
        A = np.c_[np.ones(len(y)), y]; W = np.diag(w)
        px = np.linalg.lstsq(A.T @ W @ A, A.T @ W @ s["x"], rcond=None)[0]
        yc = c["y"][i]; xp = px[0] + px[1] * yc
        two = len(disks) == 2 and min(nd) >= MINCRY
        zp = np.nan
        if two:
            pz = np.linalg.lstsq(A.T @ W @ A, A.T @ W @ z, rcond=None)[0]; zp = pz[0] + pz[1] * yc
        rows.append((int(c["sector"][i]), c["x"][i], c["z"][i], int(c["posCalc"][i]), xp, zp, len(s), int(two), px[1]))
    return np.array(rows, dtype=[("sector", int), ("xc", float), ("zc", float), ("posCalc", int), ("xp", float), ("zp", float), ("ncry", int), ("two", int), ("dxdy", float)])

def fitplot(ax, r, lo, hi, bw, label, S, key):
    b = np.arange(lo, hi + bw, bw); n, e = np.histogram(r, b); x = 0.5 * (e[1:] + e[:-1])
    ax.errorbar(x, n, yerr=np.sqrt(np.maximum(n, 1)), fmt=".", ms=3, lw=0.6, color="k", label="data")
    txt = f"entries: {len(r)}\nmedian = {np.median(r):.0f} mm\nrms (in range) = {np.std(r[(r > lo) & (r < hi)]):.0f} mm"
    if n.sum() > 30:
        try:
            p, cv = curve_fit(gflat, x, n, p0=[n.max(), x[np.argmax(n)], 100, max(np.median(n), 0.5)], sigma=np.sqrt(np.maximum(n, 1)))
            pe = np.sqrt(np.diag(cv)); chi2 = float(np.sum(((n - gflat(x, *p)) / np.sqrt(np.maximum(n, 1))) ** 2)); ndf = len(x) - 4
            xx = np.linspace(lo, hi, 600); ax.plot(xx, gflat(xx, *p), "r-", lw=1.2, label="Gaussian + flat")
            txt += f"\nfit mean = {p[1]:.0f} ± {pe[1]:.0f} mm\nfit sigma = {abs(p[2]):.0f} ± {pe[2]:.0f} mm\nflat = {p[3]:.1f} ± {pe[3]:.1f}\nchi2/ndf = {chi2:.0f}/{ndf}"
            S[key] = dict(n=int(len(r)), median=float(np.median(r)), mean=float(p[1]), mean_err=float(pe[1]), sigma=float(abs(p[2])), sigma_err=float(pe[2]), chi2=chi2, ndf=ndf)
        except Exception as ex: txt += f"\nfit failed: {ex}"
    ax.text(0.02, 0.97, txt, transform=ax.transAxes, va="top", fontsize=8.5, bbox=dict(fc="w", ec="0.7"))
    ax.set_xlabel(label); ax.set_ylabel(f"coincidences / {bw:g} mm"); ax.grid(alpha=0.3); ax.legend(loc="upper right", fontsize=8)

def main():
    out = pathlib.Path(sys.argv[1]); out.mkdir(parents=True, exist_ok=True)
    args = [a.split("=") for a in sys.argv[2:]]; runs = [pathlib.Path(p).name for p, _ in args]; T = []; fits = {}
    for p, cd in args:
        fit = json.load(open(pathlib.Path(cd) / "summary.json"))["peak_fit"]; fits[pathlib.Path(p).name] = fit
        c = np.atleast_1d(np.genfromtxt(f"{p}_crv.csv", delimiter=",", names=True)); h = np.atleast_1d(np.genfromtxt(f"{p}_calo.csv", delimiter=",", names=True))
        T.append(tracks(c, h, fit))
    t = np.concatenate(T); S = dict(runs=runs, nsig=NSIG, emin_MeV=EMIN, n_tracks=int(len(t)), n_two_disk=int(t["two"].sum()), peak_fit_used=fits)
    tt = f"runs {', '.join(runs)}; time coincidence |dt-mean| < {NSIG:g} sigma, calo hits > {EMIN:g} MeV"
    # s1, s2: z (two disks)
    fig, axs = plt.subplots(1, 2, figsize=(14, 5.5))
    for a, st in zip(axs, (1, 3)):
        m = (t["sector"] == st) & (t["two"] == 1); r = t["zc"][m] - t["zp"][m]
        fitplot(a, r, -2000, 2000, 50, "z_CRV - z_predicted [mm]", S, f"dz_{SECTOR[st]}"); a.set_title(f"CRV {SECTOR[st]}: z residual, calorimeter tracks through both disks", fontsize=10)
    fig.suptitle(tt, fontsize=10); fig.tight_layout(); fig.savefig(out / "s1_dz.png", dpi=110); plt.close(fig)
    fig, axs = plt.subplots(1, 2, figsize=(13, 6))
    for a, st in zip(axs, (1, 3)):
        m = (t["sector"] == st) & (t["two"] == 1)
        a.plot(t["zp"][m], t["zc"][m], ".", ms=2, alpha=0.4, color="k"); lim = (-4000, 2000); a.plot(lim, lim, "r--", lw=1, label="z_CRV = z_predicted")
        a.set_xlim(*lim); a.set_ylim(*lim); a.set_aspect("equal"); a.grid(alpha=0.3); a.legend(fontsize=8)
        a.set_xlabel("z predicted at the CRV height from the calorimeter track [mm]"); a.set_ylabel(f"z of the CRV {SECTOR[st]} cluster [mm]"); a.set_title(f"CRV {SECTOR[st]} ({m.sum()} two-disk tracks)", fontsize=10)
    fig.suptitle(tt + "\ndetector frame: crystal front face of disk 0 at z = 1745 mm", fontsize=10); fig.tight_layout(); fig.savefig(out / "s2_z_scatter.png", dpi=110); plt.close(fig)
    # s3: x (T2 with position along the bar), one disk enough
    fig, axs = plt.subplots(1, 2, figsize=(14, 5.5))
    m = (t["sector"] == 3) & (t["posCalc"] == 1); r = t["xc"][m] - t["xp"][m]
    fitplot(axs[0], r, -2000, 2000, 50, "x_CRV - x_predicted [mm]", S, "dx_T2"); axs[0].set_title("CRV T2 (position along the bar from the two sides): x residual", fontsize=10)
    axs[1].plot(t["xp"][m], t["xc"][m], ".", ms=2, alpha=0.4, color="k"); lim = (-2500, 2500); axs[1].plot(lim, lim, "r--", lw=1, label="x_CRV = x_predicted")
    axs[1].set_xlim(*lim); axs[1].set_ylim(*lim); axs[1].set_aspect("equal"); axs[1].grid(alpha=0.3); axs[1].legend(fontsize=8)
    axs[1].set_xlabel("x predicted at the CRV height [mm]"); axs[1].set_ylabel("x of the CRV T2 cluster [mm]")
    fig.suptitle(tt, fontsize=10); fig.tight_layout(); fig.savefig(out / "s3_dx_T2.png", dpi=110); plt.close(fig)
    json.dump(S, open(out / "summary_spatial.json", "w"), indent=1)
    print(json.dumps({k: v for k, v in S.items() if k != "peak_fit_used"}, indent=1))

if __name__ == "__main__":
    main()
