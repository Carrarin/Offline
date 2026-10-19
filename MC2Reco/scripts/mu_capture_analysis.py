"""mu- capture from the output of MC2MuCapture (capture_to_csv.C): second pulses in the STOP crystals over the THROUGH controls,
corrected for the detection efficiency measured by injection, and the decay time.

Selection: second pulse found if maxS > SCUT and A2 > A2MIN ADC; delay = tau. STOP sample with the mu+ veto (nMichel == 0) and without.
Controls: THROUGH reweighted to the first-pulse amplitude distribution of the STOP sample (the fake rate depends on the amplitude).
Excess per STOP crystal: f_stop(t) - f_through(t) (fractions per delay bin); efficiency eff(t) = fraction of the injected pulses found
with |tau - injDelay| < 15 ns (pulses with A2 above A2MIN); excess / eff fitted with N exp(-t/tauC) + C on [TLO, THI] ns.
Output: OUTDIR/{c1_maxS.png, c2_delay.png, c3_eff.png, c4_fit.png, summary.json}
usage: python3 mu_capture_analysis.py PREFIX OUTDIR [SCUT=5] [A2MIN=30] [TLO=20] [THI=600]
"""
import sys, json, pathlib, numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
sys.path.insert(0, str(pathlib.Path(__file__).parent)); from crv_calo_coinc import curve_fit
P, OUT = sys.argv[1], pathlib.Path(sys.argv[2]); OUT.mkdir(parents=True, exist_ok=True)
SCUT = float(sys.argv[3]) if len(sys.argv) > 3 else 5.; A2MIN = float(sys.argv[4]) if len(sys.argv) > 4 else 30.
TLO = float(sys.argv[5]) if len(sys.argv) > 5 else 20.; THI = float(sys.argv[6]) if len(sys.argv) > 6 else 600.; BW = 20.
C = np.genfromtxt(f"{P}_capture.csv", delimiter=",", names=True); J = np.genfromtxt(f"{P}_inject.csv", delimiter=",", names=True)
S = C[C["kind"] == 0]; T = C[C["kind"] == 1]; Sv = S[S["nMichel"] == 0]
N = dict(stop=int(len(S)), stop_veto=int(len(Sv)), through=int(len(T)), inject=int(len(J)), scut=SCUT, a2min=A2MIN, tlo=TLO, thi=THI)
# amplitude reweighting of the controls
ab = np.quantile(np.r_[S["amp"], T["amp"]], np.linspace(0, 1, 21)); ab[-1] += 1
def weights(ref, ctl):
    hr, _ = np.histogram(ref["amp"], ab); hc, _ = np.histogram(ctl["amp"], ab); r = np.where(hc > 0, hr / np.maximum(hc, 1), 0) * len(ctl) / max(len(ref), 1)
    return r[np.clip(np.digitize(ctl["amp"], ab) - 1, 0, len(r) - 1)]
found = lambda X: (X["maxS"] > SCUT) & (X["A2"] > A2MIN) & (X["tau"] >= 0)
# c1: maxS
fig, ax = plt.subplots(figsize=(10, 5)); b = np.linspace(-5, 40, 91)
for X, lab, w in ((S, "STOP", None), (Sv, "STOP, mu+ veto", None), (T, "THROUGH (reweighted)", weights(S, T))):
    ax.hist(X["maxS"], b, weights=(w if w is not None else np.ones(len(X))) / max(len(X), 1), histtype="step", lw=1.3, label=f"{lab}: {len(X)}")
ax.axvline(SCUT, color="r", ls="--", label=f"cut maxS > {SCUT:g}"); ax.set_yscale("log"); ax.set_xlabel("max of the matched filter on the residual, S"); ax.set_ylabel("fraction of crystals / bin"); ax.grid(alpha=0.3); ax.legend(fontsize=8)
fig.tight_layout(); fig.savefig(OUT / "c1_maxS.png", dpi=110); plt.close(fig)
# c2: delays
bins = np.arange(0, 1500 + BW, BW); x = 0.5 * (bins[1:] + bins[:-1])
def frac(X, w=None):
    m = found(X); ww = (w if w is not None else np.ones(len(X))); h, _ = np.histogram(X["tau"][m], bins, weights=ww[m]); h2, _ = np.histogram(X["tau"][m], bins, weights=ww[m] ** 2)
    n = ww.sum(); return h / n, np.sqrt(h2) / n
fs, es = frac(S); fsv, esv = frac(Sv); ft, et = frac(T, weights(S, T)); ftv, etv = frac(T, weights(Sv, T))
N["frac_found"] = dict(stop=float(found(S).mean()), stop_veto=float(found(Sv).mean()), through=float(np.average(found(T), weights=weights(S, T))))
fig, axs = plt.subplots(1, 2, figsize=(14, 5))
for a, (f1, e1, f2, e2, lab) in zip(axs, ((fs, es, ft, et, "STOP (all)"), (fsv, esv, ftv, etv, "STOP with the mu+ veto"))):
    a.errorbar(x, f1, e1, fmt=".", color="k", ms=3, label=lab); a.errorbar(x, f2, e2, fmt=".", color="0.6", ms=3, label="THROUGH (reweighted)")
    a.set_yscale("log"); a.set_xlabel("delay of the second pulse [ns]"); a.set_ylabel(f"fraction of crystals / {BW:g} ns"); a.grid(alpha=0.3); a.legend(fontsize=8)
fig.tight_layout(); fig.savefig(OUT / "c2_delay.png", dpi=110); plt.close(fig)
# c3: efficiency from injection
Jm = J[J["injA2"] > A2MIN]; ok = (Jm["maxS"] > SCUT) & (Jm["A2"] > A2MIN) & (np.abs(Jm["tau"] - Jm["injDelay"]) < 15)
he, _ = np.histogram(Jm["injDelay"], bins); hk, _ = np.histogram(Jm["injDelay"][ok], bins); eff = np.where(he > 0, hk / np.maximum(he, 1), np.nan); eerr = np.sqrt(np.maximum(eff * (1 - eff), 0) / np.maximum(he, 1))
fig, ax = plt.subplots(figsize=(10, 4.5)); ax.errorbar(x, eff, eerr, fmt="o", ms=3, color="k"); ax.set_ylim(0, 1.05); ax.set_xlim(0, 650); ax.grid(alpha=0.3)
ax.set_xlabel("injected delay [ns]"); ax.set_ylabel("efficiency"); ax.set_title(f"injected pulses A2 {A2MIN:g}-400 ADC in THROUGH crystals: found with maxS > {SCUT:g}, A2 > {A2MIN:g}, |tau - delay| < 15 ns", fontsize=9)
fig.tight_layout(); fig.savefig(OUT / "c3_eff.png", dpi=110); plt.close(fig)
# c4: excess / efficiency and fit
fig, axs = plt.subplots(1, 2, figsize=(14, 5)); N["fit"] = {}
for a, (f1, e1, f2, e2, lab) in zip(axs, ((fs, es, ft, et, "all STOP"), (fsv, esv, ftv, etv, "STOP with the mu+ veto"))):
    ex = (f1 - f2) / eff; ee = np.hypot(e1, e2) / eff; m = (x >= TLO) & (x <= THI) & np.isfinite(ex) & (ee > 0)
    a.errorbar(x[m], ex[m], ee[m], fmt="o", ms=3, color="k", label="(STOP - THROUGH) / efficiency")
    model = lambda t, n0, tau, c: n0 * np.exp(-t / tau) + c
    try:
        p, cv = curve_fit(model, x[m], ex[m], p0=[ex[m][0] * 1.5, 80., max(ex[m][-5:].mean(), 1e-7)], sigma=ee[m]); pe = np.sqrt(np.diag(cv))
        chi2 = float(np.sum(((ex[m] - model(x[m], *p)) / ee[m]) ** 2)); ndf = int(m.sum() - 3)
        tt = np.linspace(TLO, THI, 300); a.plot(tt, model(tt, *p), "r-", label="N exp(-t/tau) + C")
        a.text(0.35, 0.95, f"tau = {p[1]:.1f} ± {pe[1]:.1f} ns\nN = {p[0]:.2e} ± {pe[0]:.1e}\nC = {p[2]:.2e} ± {pe[2]:.1e}\nchi2/ndf = {chi2:.1f}/{ndf}\nfit {TLO:g}-{THI:g} ns", transform=a.transAxes, va="top", fontsize=9, bbox=dict(fc="w", ec="0.7"))
        N["fit"][lab] = dict(tau=float(p[1]), tau_err=float(pe[1]), N=float(p[0]), C=float(p[2]), C_err=float(pe[2]), chi2=chi2, ndf=ndf)
    except Exception as exn: a.text(0.4, 0.9, f"fit failed: {exn}", transform=a.transAxes)
    a.set_title(lab, fontsize=10); a.set_xlabel("delay of the second pulse [ns]"); a.set_ylabel(f"excess per STOP crystal / {BW:g} ns (efficiency corrected)"); a.grid(alpha=0.3); a.legend(fontsize=8, loc="upper right")
fig.tight_layout(); fig.savefig(OUT / "c4_fit.png", dpi=110); plt.close(fig)
json.dump(N, open(OUT / "summary.json", "w"), indent=1); print(json.dumps(N, indent=1))
