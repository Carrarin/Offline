#!/usr/bin/env python3
"""Checks on the Monte Carlo truth of the stopped muons (MC2StopTruth trees of the hall production): where they stop, what their
decay / capture products leave in the crystals and when, and whether the counts close against the decays counted by the
reconstruction of the same production (sim_vs_data.py).

usage: PYTHONPATH=$ROOTSYS/lib python3 Offline/MC2Reco/scripts/stop_truth_check.py <truth dir> <livetime s> <reco-counted decays> <out dir>
       (plain ROOT 6.36 Python, not inside the muse environment)
"""
import glob, json, math, pathlib, sys
import numpy as np
import ROOT
ROOT.gROOT.SetBatch(True)
DECAY, DIO, NCAP = 14, 114, 115           # ProcessCode ids (Offline MCDataProducts/inc/ProcessCode.hh)
TD, LIVE, NRECO, OUT = pathlib.Path(sys.argv[1]), float(sys.argv[2]), float(sys.argv[3]), pathlib.Path(sys.argv[4]); OUT.mkdir(parents=True, exist_ok=True)
ch = ROOT.TChain("truth/stops"); [ch.Add(f) for f in sorted(glob.glob(str(TD / "truth_*.root")))]
cols = {k: [] for k in ("q", "inCry", "endX", "endY", "endZ", "eTot", "eStop", "eNear", "eFar", "ePrompt", "tdec", "tcap", "nCryMu")}
win, swin = [], []
for e in ch:
    dc, dt = list(e.dCode), list(e.dT)
    for k in ("q", "inCry", "endX", "endY", "endZ", "eTot", "eStop", "eNear", "eFar", "ePrompt", "nCryMu"): cols[k].append(getattr(e, k))
    cols["tdec"].append(min([t for c, t in zip(dc, dt) if c in (DECAY, DIO)], default=-1.0))
    cols["tcap"].append(min([t for c, t in zip(dc, dt) if c == NCAP], default=-1.0))
    win.append(list(e.eWin)); swin.append(list(e.eStopWin))
d = {k: np.array(v) for k, v in cols.items()}; win, swin = np.array(win), np.array(swin)
q, inC, eT = d["q"], d["inCry"].astype(bool), d["eTot"]; R = {}
R["stops"] = {f"{nm}_{lab}": int(((q == qq) & m).sum()) for qq, nm in ((1, "mu+"), (-1, "mu-")) for lab, m in (("inCry", inC), ("outside", ~inC))}
R["stops_visible5MeV"] = {f"{nm}_{lab}": int(((q == qq) & m & (eT > 5)).sum()) for qq, nm in ((1, "mu+"), (-1, "mu-")) for lab, m in (("inCry", inC), ("outside", ~inC))}
R["rates_Hz"] = {k: v / LIVE for k, v in R["stops"].items()}
m = (q > 0) & inC
R["mu+_inCry_products"] = {"median_MeV": float(np.median(eT[m])), "frac_gt20": float(np.mean(eT[m] > 20)), "mean_stop_near_far": [float(d[k][m].mean()) for k in ("eStop", "eNear", "eFar")]}
late = (d["tdec"] >= 1000) & (d["tdec"] < 20000)
R["decays_1_20us_products_gt5MeV"] = {lab: int((late & mm & (eT > 5)).sum()) for lab, mm in (("inCry", inC), ("all", np.ones_like(inC)))}
R["predicted_observed_x_window0785"] = {k: 0.785 * v for k, v in R["decays_1_20us_products_gt5MeV"].items()}
R["reco_counted"] = NRECO
tp = d["tdec"][(q > 0) & (d["tdec"] >= 0)]
cap = (q < 0) & inC & (d["tcap"] >= 0)
tm = np.where(d["tcap"] >= 0, d["tcap"], d["tdec"])[(q < 0) & inC & ((d["tcap"] >= 0) | (d["tdec"] >= 0))]
def tau_ml(t, lo, hi):
    t = t[(t >= lo) & (t < hi)]; g = np.linspace(max(1e-3, 0.002 * hi), 5 * hi, 20000)
    v = [-len(t) * math.log(x * (math.exp(-lo / x) - math.exp(-hi / x))) - t.sum() / x for x in g]; tau = g[int(np.argmax(v))]
    return float(tau), float(tau / math.sqrt(len(t))), int(len(t))
R["tau_mu+_ns"] = tau_ml(tp, 0, 20000); R["tau_mu-_inCry_ns"] = tau_ml(tm, 0, 2000)
R["mu-_inCry_capture"] = {"n": int(cap.sum()), "of_mu-_inCry": int(((q < 0) & inC).sum()), "dio": int(((q < 0) & inC & (d["tdec"] >= 0)).sum()),
    "E_frac_gt1_gt5_gt10": [float(np.mean(eT[cap] > x)) for x in (1, 5, 10)], "E_median_84_99": [float(np.percentile(eT[cap], p)) for p in (50, 84, 99)],
    "E_by_window_mean": win[cap].mean(0).tolist(), "E_stop_crystal_by_window_mean": swin[cap].mean(0).tolist(),
    "E_stop_near_far_mean": [float(d[k][cap].mean()) for k in ("eStop", "eNear", "eFar")],
    "frac_gt5MeV_stop_crystal_within_500ns": float(np.mean(swin[cap][:, 0] + swin[cap][:, 1] > 5))}
(OUT / "stop_truth_check.json").write_text(json.dumps(R, indent=1)); print(json.dumps(R, indent=1))
np.savez(OUT / "stop_truth.npz", **d, win=win, swin=swin)

import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
fig, ax = plt.subplots(1, 3, figsize=(15, 4.4))
ax[0].hist(tm, bins=np.linspace(0, 600, 61), histtype="step", color="k", label=f"μ⁻ stopped in a crystal: τ = {R['tau_mu-_inCry_ns'][0]:.0f} ± {R['tau_mu-_inCry_ns'][1]:.0f} ns")
ax[0].set_yscale("log"); ax[0].set_xlabel("capture or DIO time after the stop [ns]"); ax[0].set_ylabel("μ⁻ / 10 ns"); ax[0].legend(fontsize=7.5); ax[0].set_title("Geant4 (ShieldingM): μ⁻ lifetime in CsI", fontsize=10)
bins = np.linspace(0, 80, 81)
ax[1].hist(eT[cap], bins=bins, histtype="step", color="C0", label="μ⁻ nuclear capture"); ax[1].hist(eT[(q > 0) & inC & (d["tdec"] >= 0)], bins=bins, histtype="step", color="C3", label="μ⁺ decay")
ax[1].set_yscale("log"); ax[1].set_xlabel("energy of the products in the crystals [MeV, Geant4]"); ax[1].set_ylabel("stopped muons / MeV"); ax[1].legend(fontsize=7.5); ax[1].set_title("what the products leave (no detector response)", fontsize=10)
lab = ["0–50 ns", "50–500 ns", "0.5–20 µs", "> 20 µs"]; x = np.arange(4)
ax[2].bar(x - 0.2, win[cap].mean(0), 0.4, color="C0", label="all crystals"); ax[2].bar(x + 0.2, swin[cap].mean(0), 0.4, color="C2", label="stop crystal")
ax[2].set_xticks(x); ax[2].set_xticklabels(lab); ax[2].set_ylabel("mean energy per capture [MeV]"); ax[2].legend(fontsize=7.5); ax[2].set_title("μ⁻ capture: when the energy arrives", fontsize=10)
fig.suptitle("Monte Carlo truth of the stopped muons, hall production", fontsize=11); fig.tight_layout(); fig.savefig(OUT / "stop_truth.png", dpi=125)
