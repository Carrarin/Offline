#!/usr/bin/env python3
"""Signature of mu- captures and mu- decays in orbit (DIO) in the hall Monte Carlo, with the selection used on the data (MC2MuCapture)
and the truth of the same events (MC2StopTruth), from the files of mc2_digi_reco_capture.fcl.

Candidate (tree capture/capture): kind 0 STOP (last crystal of a one-disk track, next crystal empty), kind 1 THROUGH (middle crystal).
Truth class of a candidate: a stopped muon of the same event whose stop crystal (or, if it stopped outside the crystals, the last
crystal it crossed) is within 50 mm of the candidate crystal in the same disk: mu- capture, mu- DIO, mu+ decay; otherwise "no stop".
Observables (as on the data):
  second pulse in the candidate crystal: found = maxS > 5 and A2 > 30 ADC (MC2MuCapture example thresholds); A2 in MeV = A2 / 22.7;
  delayed cluster (tree capture/cluster, crystals within 80 mm + second pulse): dt in [30, 300] ns with E > 5 MeV / > 15 MeV;
  Michel: nMichel >= 1 (CaloHit within 80 mm, 1-20 us later).
Scaled to the data: per STOP candidate, and to the 79709 STOP / 157770 THROUGH of run 125367 (raw files 901-931, notebook).
usage: PYTHONPATH=$ROOTSYS/lib python3 Offline/MC2Reco/scripts/capture_signature.py <digicap dir> <DMAP.txt> <out dir>
"""
import glob, json, math, pathlib, sys
from collections import Counter, defaultdict
import numpy as np
import ROOT
ROOT.gROOT.SetBatch(True)
D, DMAP, OUT = pathlib.Path(sys.argv[1]), pathlib.Path(sys.argv[2]), pathlib.Path(sys.argv[3]); OUT.mkdir(parents=True, exist_ok=True)
DECAY, DIO, NCAP = 14, 114, 115
pos = {}
for ln in DMAP.read_text().splitlines()[1:]:
    p = ln.split()
    if len(p) > 17 and p[17] == "CAL" and p[2] in ("0", "1"): pos[int(p[14])] = (int(p[2]), float(p[12]), float(p[13]))
def near(c1, c2, r=50.):
    a, b = pos.get(c1), pos.get(c2)
    return a is not None and b is not None and a[0] == b[0] and math.hypot(a[1] - b[1], a[2] - b[2]) <= r
CL = ["mu- capture", "mu- DIO", "mu+ decay", "no stop"]
rows = []
for fn in sorted(glob.glob(str(D / "digi_*.root"))):
    f = ROOT.TFile.Open(fn)
    tr = defaultdict(list)
    for e in f.Get("truth/stops"):
        dc, dt = list(e.dCode), list(e.dT)
        if e.q < 0:
            tcap = min([t for c, t in zip(dc, dt) if c == NCAP], default=-1.0); tdio = min([t for c, t in zip(dc, dt) if c == DIO], default=-1.0)
            cls = "mu- DIO" if tdio >= 0 else ("mu- capture" if tcap >= 0 else None)
        else: cls = "mu+ decay"
        if cls: tr[e.event].append((cls, e.cry, bool(e.inCry), e.eTot, list(e.eWin), list(e.eStopWin)))
    clu = defaultdict(list)
    for c in f.Get("capture/cluster"): clu[(c.event, c.cry, c.kind)].append((c.dt, c.E, c.Estop))
    for c in f.Get("capture/capture"):
        lab, inc, eT, ew, esw = "no stop", False, 0., [0] * 4, [0] * 4
        for cls, cry, ic, et, w, sw in tr.get(c.event, []):
            if cry == c.cry or near(cry, c.cry): lab, inc, eT, ew, esw = cls, ic, et, w, sw; break
        cs = clu.get((c.event, c.cry, c.kind), [])
        rows.append(dict(kind=c.kind, cls=lab, inCry=inc, eTrue=eT, eTrueWin=ew, eTrueStop=esw, maxS=c.maxS, A2=c.A2, tau=c.tau, Q=c.Q, amp=c.amp,
                         nMichel=c.nMichel, cl5=any(30 <= x[0] <= 300 and x[1] > 5 for x in cs), cl15=any(30 <= x[0] <= 300 and x[1] > 15 for x in cs),
                         clE=max([x[1] for x in cs if 30 <= x[0] <= 300], default=0.)))
    f.Close()
R = {}
for kind, kn in ((0, "STOP"), (1, "THROUGH")):
    rk = [r for r in rows if r["kind"] == kind]; n = len(rk); cnt = Counter(r["cls"] for r in rk)
    out = {"candidates": n, "classes": {c: cnt.get(c, 0) for c in CL}, "fraction": {c: round(cnt.get(c, 0) / n, 4) for c in CL}, "signature": {}}
    for c in CL:
        rc = [r for r in rk if r["cls"] == c]
        if not rc: continue
        sp = [r for r in rc if r["maxS"] > 5 and r["A2"] > 30]
        out["signature"][c] = {"n": len(rc), "in_crystal": round(np.mean([r["inCry"] for r in rc]), 3),
            "second_pulse": round(len(sp) / len(rc), 3), "second_pulse_A2_MeV_median": round(float(np.median([r["A2"] / 22.7 for r in sp])), 1) if sp else None,
            "second_pulse_tau_ns_median": round(float(np.median([r["tau"] for r in sp])), 0) if sp else None,
            "cluster_30_300ns_gt5MeV": round(np.mean([r["cl5"] for r in rc]), 3), "cluster_30_300ns_gt15MeV": round(np.mean([r["cl15"] for r in rc]), 3),
            "michel_1_20us": round(np.mean([r["nMichel"] >= 1 for r in rc]), 3),
            "true_E_products_median_MeV": round(float(np.median([r["eTrue"] for r in rc])), 1),
            "true_E_0_500ns_mean_MeV": round(float(np.mean([r["eTrueWin"][0] + r["eTrueWin"][1] for r in rc])), 2)}
    R[kn] = out
nS, nT = R["STOP"]["candidates"], R["THROUGH"]["candidates"]
R["scaled_to_125367"] = {"STOP": {c: round(79709 * R["STOP"]["classes"][c] / nS) for c in CL}, "THROUGH": {c: round(157770 * R["THROUGH"]["classes"][c] / nT) for c in CL}}
exp = {}
for c in CL:
    s = R["STOP"]["signature"].get(c)
    if s: exp[c] = {k: round(R["scaled_to_125367"]["STOP"][c] * s[k]) for k in ("second_pulse", "cluster_30_300ns_gt5MeV", "cluster_30_300ns_gt15MeV", "michel_1_20us")}
R["expected_in_125367_STOP"] = exp
(OUT / "capture_signature.json").write_text(json.dumps(R, indent=1)); print(json.dumps(R, indent=1))

import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
fig, ax = plt.subplots(1, 3, figsize=(15, 4.4)); col = {"mu- capture": "C0", "mu- DIO": "C1", "mu+ decay": "C3", "no stop": "0.45"}
st = [r for r in rows if r["kind"] == 0]
for c in CL:
    rc = [r for r in st if r["cls"] == c]
    if not rc: continue
    ax[0].hist([r["A2"] / 22.7 for r in rc if r["maxS"] > 5], bins=np.linspace(0, 40, 41), histtype="step", color=col[c], density=True, label=f"{c} ({len(rc)})")
    ax[1].hist([r["clE"] for r in rc if r["clE"] > 0], bins=np.linspace(0, 60, 61), histtype="step", color=col[c], density=True, label=c)
    ax[2].hist([r["tau"] for r in rc if r["maxS"] > 5 and r["A2"] > 30], bins=np.linspace(0, 450, 46), histtype="step", color=col[c], density=True, label=c)
ax[0].set_xlabel("second pulse A2 [MeV] (maxS > 5)"); ax[0].set_ylabel("normalised"); ax[0].set_yscale("log"); ax[0].legend(fontsize=7.5); ax[0].set_title("second pulse in the stop crystal", fontsize=10)
ax[1].set_xlabel("delayed cluster energy, 30–300 ns [MeV]"); ax[1].set_yscale("log"); ax[1].legend(fontsize=7.5); ax[1].set_title("cluster in the crystals around (and the stop crystal)", fontsize=10)
ax[2].set_xlabel("delay of the second pulse [ns]"); ax[2].legend(fontsize=7.5); ax[2].set_title("when the second pulse comes", fontsize=10)
fig.suptitle("Hall Monte Carlo, STOP candidates of MC2MuCapture by truth class", fontsize=11); fig.tight_layout(); fig.savefig(OUT / "capture_signature.png", dpi=125)
