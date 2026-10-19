"""mu- capture with the delayed deposits in the crystals around the stopping crystal (tree "delayed" of MC2MuCapture).

One entry per neighbouring crystal deposit (SiPM 0 of the crystal only, so that L and R are not counted twice), ring 1 = within 80 mm of
the candidate crystal (not on the track), ring 2 = 250-400 mm (accidentals). SiPMs without a time offset from the tracks (calibration
stats file) are left out: their deposits can be displaced by board-level offsets. Samples: all STOP, STOP with the mu+ veto, THROUGH;
rates per candidate per ns. The ring-1 spectrum of all STOP is fitted with A exp(-t/tauC) + B exp(-t/2197 ns) + C (capture, mu+ decay,
flat) from several start times (the stability of tauC with the start tells a real fast component from the tail of the prompt deposits).
Output: OUTDIR/{n1_dt.png, n2_fit.png, summary.json}
usage: python3 neighbour_delayed_analysis.py DELAYED.csv CAPTURE.csv CALIB_STATS.csv OUTDIR [EMIN_MeV=0]
"""
import sys, json, pathlib, numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
sys.path.insert(0, str(pathlib.Path(__file__).parent)); from crv_calo_coinc import curve_fit
D = np.genfromtxt(sys.argv[1], delimiter=",", names=True); C = np.genfromtxt(sys.argv[2], delimiter=",", names=True)
st = np.genfromtxt(sys.argv[3], delimiter=",", names=True); OUT = pathlib.Path(sys.argv[4]); OUT.mkdir(parents=True, exist_ok=True)
EMIN = float(sys.argv[5]) if len(sys.argv) > 5 else 0.
good = set(st["roid"][st["timeFromTracks"] == 1].astype(int))
D = D[(D["sipm"] % 2 == 0) & (D["e"] > EMIN) & np.isin(D["sipm"].astype(int), np.array(sorted(good)))]
edges = np.r_[np.arange(-1000, 0, 50), np.arange(0, 500, 10), np.arange(500, 2000, 100), np.arange(2000, 10001, 500)]; xc = 0.5 * (edges[1:] + edges[:-1]); bw = np.diff(edges)
samples = {"STOP (all)": (C["kind"] == 0, D["kind"] == 0), "STOP with a mu+ decay tag (control)": ((C["kind"] == 0) & (C["nMichel"] > 0), (D["kind"] == 0) & (D["nMichel"] > 0)), "THROUGH": (C["kind"] == 1, D["kind"] == 1)}
N = dict(emin=EMIN, sipms_used=len(good)); R = {}
fig, axs = plt.subplots(1, 2, figsize=(15, 5.5))
for (name, (mc, md)), col in zip(samples.items(), ("k", "C3", "0.5")):
    nc = int(mc.sum())
    for ring, a in ((1, axs[0]), (2, axs[1])):
        d = D[md & (D["ring"] == ring)]["dt"]; h, _ = np.histogram(d, edges); r = h / bw / max(nc, 1); e = np.sqrt(np.maximum(h, 1)) / bw / max(nc, 1)
        flat = ((d >= -1000) & (d < -50)).sum() / 950. / max(nc, 1); R[(name, ring)] = (r, e, flat, h, nc)
        k = h > 0; a.errorbar(xc[k], r[k], e[k], fmt=".", ms=3, color=col, label=f"{name}: {nc} candidates, before the muon {flat:.1e} /ns"); a.axhline(flat, color=col, ls=":", lw=0.8)
    N[name] = dict(candidates=nc)
for a, t in zip(axs, ("ring 1 (within 80 mm, not on the track)", "ring 2 (250-400 mm, accidentals)")):
    a.set_xscale("symlog", linthresh=100); a.set_yscale("log"); a.set_ylim(1e-9, 0.2); a.set_xlabel("dt = t(deposit) - t(candidate crystal) [ns]"); a.set_ylabel("deposits per candidate per ns"); a.grid(alpha=0.3); a.legend(fontsize=8); a.set_title(t, fontsize=10)
fig.suptitle(f"Delayed deposits around the candidate crystals (SiPM 0 of each crystal, calibrated SiPMs only, E > {EMIN:g} MeV)", fontsize=10.5)
fig.tight_layout(); fig.savefig(OUT / "n1_dt.png", dpi=110); plt.close(fig)
r, e, flat, h, nc = R[("STOP (all)", 1)]
model = lambda t, A, tauC, B, Cc: A * np.exp(-t / tauC) + B * np.exp(-t / 2197.) + Cc
m2 = lambda t, B, Cc: B * np.exp(-t / 2197.) + Cc
fig, ax = plt.subplots(figsize=(10.5, 5.8)); m0 = (xc >= 20) & (h > 0)
ax.errorbar(xc[m0], r[m0], e[m0], fmt="o", ms=3, color="k", label="STOP (all), ring 1")
rT, eT, _, hT, _ = R[("THROUGH", 1)]; ax.errorbar(xc[m0], rT[m0], eT[m0], fmt=".", ms=3, color="0.6", label="THROUGH, ring 1")
N["fits"] = {}; lines = []
for TLO, col in ((20, "r"), (50, "m"), (100, "b")):
    m = (xc >= TLO) & (h > 0)
    try:
        p, cv = curve_fit(model, xc[m], r[m], p0=[3e-5, 80., 1e-5, 1e-8], sigma=e[m]); pe = np.sqrt(np.diag(cv)); c1 = float(np.sum(((r[m] - model(xc[m], *p)) / e[m]) ** 2))
        q, _ = curve_fit(m2, xc[m], r[m], p0=[1e-5, 1e-8], sigma=e[m]); c0 = float(np.sum(((r[m] - m2(xc[m], *q)) / e[m]) ** 2))
        tt = np.linspace(TLO, 10000, 2000); ax.plot(tt, model(tt, *p), "-", color=col, lw=1, label=f"fit from {TLO} ns")
        N["fits"][str(TLO)] = dict(tauC=float(p[1]), tauC_err=float(pe[1]), A=float(p[0]), A_err=float(pe[0]), B=float(p[2]), B_err=float(pe[2]), C=float(p[3]), chi2=c1, ndf=int(m.sum() - 4), chi2_no_capture=c0)
        lines.append(f"from {TLO:3d} ns: tauC = {p[1]:.0f} ± {pe[1]:.0f} ns, A = {p[0]:.1e} ± {pe[0]:.0e}, chi2/ndf {c1:.0f}/{m.sum()-4} (no capture {c0:.0f})")
    except Exception as exn: lines.append(f"from {TLO} ns: fit failed ({exn})")
ax.text(0.02, 0.03, "A exp(-t/tauC) + B exp(-t/2197 ns) + C\n" + "\n".join(lines), transform=ax.transAxes, fontsize=7.5, va="bottom", bbox=dict(fc="w", ec="0.7"))
ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlabel("dt [ns]"); ax.set_ylabel("deposits per candidate per ns"); ax.grid(alpha=0.3); ax.legend(fontsize=8, loc="upper right")
fig.tight_layout(); fig.savefig(OUT / "n2_fit.png", dpi=110); plt.close(fig)
json.dump(N, open(OUT / "summary.json", "w"), indent=1); print(json.dumps(N, indent=1))
