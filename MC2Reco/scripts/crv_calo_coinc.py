"""Time coincidences between CRV coincidence clusters and calorimeter hits (runs with the CRV, KPP extracted position).

Input: the tables of crv_calo_flatten.C (<prefix>_crv.csv, <prefix>_calo.csv), one or more prefixes.
For every event with CRV coincidences: every pair (CRV cluster, calorimeter hit with eDep > EMIN) gives dt = t_calo - t_CRV.
Figures: c1 wide dt (find the peak), c2 dt around the peak per CRV sector type with a Gaussian + flat fit (parameters on the plot).
Output: <outdir>/c1_dt_wide.png, c2_dt_peak.png, summary.json
usage: python3 crv_calo_coinc.py OUTDIR PREFIX [PREFIX ...]
"""
import sys, json, pathlib, numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

def curve_fit(f, x, y, p0, sigma, absolute_sigma=True, nit=200):
    """weighted least squares, Levenberg-Marquardt with numerical derivatives (numpy only); returns (p, covariance)"""
    p = np.array(p0, float); w = 1.0 / np.asarray(sigma, float); lam = 1e-3
    r = lambda q: (y - f(x, *q)) * w
    def jac(q):
        J = np.empty((len(x), len(q)))
        for j in range(len(q)):
            h = 1e-6 * max(abs(q[j]), 1e-3); d = q.copy(); d[j] += h; J[:, j] = (f(x, *d) - f(x, *q)) * w / h
        return J
    c0 = np.sum(r(p) ** 2)
    for _ in range(nit):
        J = jac(p); A = J.T @ J; g = J.T @ r(p)
        step = np.linalg.solve(A + lam * np.diag(np.diag(A)), g); c1 = np.sum(r(p + step) ** 2)
        if c1 < c0:
            p, lam, done = p + step, lam / 10, abs(c0 - c1) < 1e-8 * max(c0, 1); c0 = c1
            if done: break
        else: lam *= 10
    J = jac(p); return p, np.linalg.inv(J.T @ J)

EMIN = 5.0          # MeV, calorimeter hit threshold (Offline energy scale)
SECTOR = {1: "EX", 2: "T1", 3: "T2", 4: "M1-M8"}   # sectorType as in crv_counters_extracted_v04.txt

def load(prefix):
    c = np.genfromtxt(f"{prefix}_crv.csv", delimiter=",", names=True)
    h = np.genfromtxt(f"{prefix}_calo.csv", delimiter=",", names=True)
    return np.atleast_1d(c), np.atleast_1d(h)

def pairs(c, h):
    h = h[h["eDep"] > EMIN]
    o = np.argsort(h["ev"], kind="stable"); h = h[o]; ev = h["ev"]
    lo = np.searchsorted(ev, c["ev"], "left"); hi = np.searchsorted(ev, c["ev"], "right")
    dt, sec, idx = [], [], []
    for i in range(len(c)):
        if hi[i] > lo[i]:
            t = h["time"][lo[i]:hi[i]] - c["time"][i]; dt.append(t); sec.append(np.full(len(t), c["sector"][i])); idx.append(np.full(len(t), i))
    return (np.concatenate(dt), np.concatenate(sec).astype(int), np.concatenate(idx)) if dt else (np.array([]),) * 3

def gflat(x, a, mu, s, b): return a * np.exp(-0.5 * ((x - mu) / s) ** 2) + b

def main():
    out = pathlib.Path(sys.argv[1]); out.mkdir(parents=True, exist_ok=True); runs = [pathlib.Path(p).name for p in sys.argv[2:]]
    C, H = zip(*[load(p) for p in sys.argv[2:]])
    DT, SEC = [], []
    for c, h in zip(C, H):
        d, s, _ = pairs(c, h); DT.append(d); SEC.append(s)
    dt, sec = np.concatenate(DT), np.concatenate(SEC)
    S = dict(runs=runs, emin_MeV=EMIN, n_crv_clusters=int(sum(len(c) for c in C)), n_pairs=int(len(dt)))
    # c1: wide
    fig, ax = plt.subplots(figsize=(11, 5))
    b = np.arange(-5000, 5000 + 20, 20); n, e = np.histogram(dt, b); x = 0.5 * (e[1:] + e[:-1])
    ax.errorbar(x, n, yerr=np.sqrt(n), fmt=".", ms=3, lw=0.6, color="k")
    k = np.argmax(n); peak = x[k]; S["wide_peak_ns"] = float(peak)
    ax.set_xlabel("t_calo - t_CRV [ns]  (calorimeter hit time - CRV coincidence time, raw, no offsets)"); ax.set_ylabel("pairs / 20 ns")
    ax.set_title(f"CRV-calorimeter time difference, all pairs (calo hits > {EMIN:g} MeV), runs {', '.join(runs)}", fontsize=10)
    ax.text(0.02, 0.92, f"pairs: {len(dt)}\nhighest bin at {peak:.0f} ns", transform=ax.transAxes, va="top", fontsize=9, bbox=dict(fc="w", ec="0.7"))
    ax.set_yscale("log"); ax.grid(alpha=0.3); fig.tight_layout(); fig.savefig(out / "c1_dt_wide.png", dpi=110); plt.close(fig)
    # c2: around the peak, per sector
    fig, axs = plt.subplots(2, 2, figsize=(13, 8.5)); S["peak_fit"] = {}
    for a, st in zip(axs.flat, (1, 2, 3, 4)):
        m = sec == st; b = np.arange(peak - 150, peak + 150 + 2, 2); n, e = np.histogram(dt[m], b); x = 0.5 * (e[1:] + e[:-1])
        a.errorbar(x, n, yerr=np.sqrt(np.maximum(n, 1)), fmt=".", ms=3, lw=0.6, color="k", label="data")
        txt = f"pairs in window: {int(n.sum())}"
        if n.sum() > 30:
            try:
                p0 = [n.max() - np.median(n), x[np.argmax(n)], 5, max(np.median(n), 0.1)]
                p, cv = curve_fit(gflat, x, n, p0=p0, sigma=np.sqrt(np.maximum(n, 1)), absolute_sigma=True)
                pe = np.sqrt(np.diag(cv)); chi2 = float(np.sum(((n - gflat(x, *p)) / np.sqrt(np.maximum(n, 1))) ** 2)); ndf = len(x) - 4
                xx = np.linspace(x[0], x[-1], 600); a.plot(xx, gflat(xx, *p), "r-", lw=1.2, label="Gaussian + flat")
                nsig = p[0] * abs(p[2]) * np.sqrt(2 * np.pi) / 2.0
                txt += f"\nmean = {p[1]:.2f} ± {pe[1]:.2f} ns\nsigma = {abs(p[2]):.2f} ± {pe[2]:.2f} ns\nflat = {p[3]:.2f} ± {pe[3]:.2f} /2 ns\nsignal ≈ {nsig:.0f} pairs\nchi2/ndf = {chi2:.0f}/{ndf}"
                S["peak_fit"][SECTOR[st]] = dict(mean=float(p[1]), mean_err=float(pe[1]), sigma=float(abs(p[2])), sigma_err=float(pe[2]), flat=float(p[3]), signal_pairs=float(nsig), chi2=chi2, ndf=ndf)
            except Exception as ex: txt += f"\nfit failed: {ex}"
        a.text(0.02, 0.97, txt, transform=a.transAxes, va="top", fontsize=8.5, bbox=dict(fc="w", ec="0.7"))
        a.set_title(f"CRV sector {SECTOR[st]}", fontsize=10); a.set_xlabel("t_calo - t_CRV [ns]"); a.set_ylabel("pairs / 2 ns"); a.grid(alpha=0.3); a.legend(loc="upper right", fontsize=8)
    fig.suptitle(f"CRV-calorimeter coincidence peak per CRV sector (calo hits > {EMIN:g} MeV), runs {', '.join(runs)}", fontsize=11)
    fig.tight_layout(); fig.savefig(out / "c2_dt_peak.png", dpi=110); plt.close(fig)
    json.dump(S, open(out / "summary.json", "w"), indent=1); print(json.dumps(S, indent=1))

if __name__ == "__main__":
    main()
