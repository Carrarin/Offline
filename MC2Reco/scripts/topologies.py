#!/usr/bin/env python3
"""Topologies of the stopped muons in the hall Monte Carlo, in simple terms, and two ways of finding the mu- capture:
  A. stop crystal: second pulse of MC2MuCapture (maxS > 5, A2 > 30 ADC)
  B. neighbours: a crystal within 80 mm of the stop crystal (same disk), not on the muon track (no hit within +-30 ns of the stop),
     with a CaloHit of both SiPMs 30-500 ns after the stop hit; energy = the largest such crystal and the sum.
Truth class of each STOP candidate as capture_signature.py. Event displays: one example per topology, crystals around the stop crystal
in three time windows (prompt -30..+30 ns, 30-500 ns, 1-20 us).
usage: PYTHONPATH=$ROOTSYS/lib python3 Offline/MC2Reco/scripts/topologies.py <digicap dir> <DMAP.txt> <out dir>"""
import glob, json, math, pathlib, sys
from collections import Counter, defaultdict
import numpy as np
import ROOT
D, DMAP, OUT = pathlib.Path(sys.argv[1]), pathlib.Path(sys.argv[2]), pathlib.Path(sys.argv[3]); OUT.mkdir(parents=True, exist_ok=True)
DECAY, DIO, NCAP = 14, 114, 115
pos = {}
for ln in DMAP.read_text().splitlines()[1:]:
    p = ln.split()
    if len(p) > 17 and p[17] == "CAL" and p[2] in ("0", "1"): pos[int(p[14])] = (int(p[2]), float(p[12]), float(p[13]))
def dist(a, b):
    A, B = pos.get(a), pos.get(b)
    return 1e9 if (A is None or B is None or A[0] != B[0]) else math.hypot(A[1] - B[1], A[2] - B[2])
CL = ["mu- capture", "mu- DIO", "mu+ decay", "no stop"]
rows, ex = [], defaultdict(list)
for fn in sorted(glob.glob(str(D / "digi_*.root"))):
    f = ROOT.TFile.Open(fn); tr = defaultdict(list); hits = {}
    for e in f.Get("truth/stops"):
        dc = list(e.dCode)
        cls = ("mu- DIO" if DIO in dc else ("mu- capture" if NCAP in dc else None)) if e.q < 0 else "mu+ decay"
        if cls: tr[e.event].append((cls, e.cry, bool(e.inCry), list(e.eStopWin), list(e.eNearWin), list(e.eFarWin)))
    for h in f.Get("hits/hits"): hits[h.event] = (list(h.cry), list(h.t), list(h.e), list(h.nSiPM))
    for c in f.Get("capture/capture"):
        if c.kind != 0 or c.event not in hits: continue
        cls, inc, tw = "no stop", False, None
        for cl, cry, ic, sw, nw, fw in tr.get(c.event, []):
            if cry == c.cry or dist(cry, c.cry) <= 50: cls, inc, tw = cl, ic, (sw, nw, fw); break
        hc, ht, he, hn = hits[c.event]
        t0 = [t for cr, t, en in zip(hc, ht, he) if cr == c.cry]
        if not t0: continue
        t0 = min(t0, key=lambda t: abs(t - c.t1)) if c.t1 > 0 else min(t0)
        prompt = {cr for cr, t in zip(hc, ht) if abs(t - t0) < 30}
        nb = defaultdict(float)
        for cr, t, en, ns in zip(hc, ht, he, hn):
            if cr == c.cry or cr in prompt or ns < 2 or dist(cr, c.cry) > 80: continue
            if 30 <= t - t0 <= 500: nb[cr] += en
        late = any(1000 <= t - t0 <= 20000 and dist(cr, c.cry) <= 80 and en > 5 for cr, t, en in zip(hc, ht, he))
        r = dict(cls=cls, inCry=inc, A=(c.maxS > 5 and c.A2 > 30), A2=c.A2 / 22.7, nbMax=max(nb.values(), default=0.), nbSum=sum(nb.values()), late=late, truthWin=tw)
        rows.append(r)
        if len(ex[cls]) < 300: ex[cls].append((fn, c.event, c.cry, t0, r))
    f.Close()
res = {}
for cl in CL:
    rc = [r for r in rows if r["cls"] == cl]; n = len(rc)
    if not n: continue
    res[cl] = {"n": n, "A_stop_crystal_second_pulse": round(np.mean([r["A"] for r in rc]), 4),
               "B_neighbour_ge3MeV": round(np.mean([r["nbMax"] >= 3 for r in rc]), 4), "B_neighbour_ge5MeV": round(np.mean([r["nbMax"] >= 5 for r in rc]), 4),
               "A_or_B3": round(np.mean([r["A"] or r["nbMax"] >= 3 for r in rc]), 4),
               "B3_and_no_late_electron": round(np.mean([(r["nbMax"] >= 3) and not r["late"] for r in rc]), 4),
               "late_electron_1_20us": round(np.mean([r["late"] for r in rc]), 4)}
    tw = [r["truthWin"] for r in rc if r["truthWin"]]
    if tw:
        a = np.array(tw)   # [n, 3 (stop, near, far), 4 windows]
        res[cl]["truth_mean_MeV_0_500ns"] = {"stop": round(float((a[:, 0, 0] + a[:, 0, 1]).mean()), 2), "near": round(float((a[:, 1, 0] + a[:, 1, 1]).mean()), 2), "far": round(float((a[:, 2, 0] + a[:, 2, 1]).mean()), 2)}
CLS = [c for c in CL if c in res]; nS = sum(res[c]["n"] for c in CLS)
for key in ("A_stop_crystal_second_pulse", "B_neighbour_ge3MeV", "B3_and_no_late_electron"):
    tot = sum(res[c]["n"] * res[c][key] for c in CLS); cap = res["mu- capture"]["n"] * res["mu- capture"][key]
    res[f"purity_{key}"] = round(cap / tot, 3) if tot else None
    res[f"per_79709_STOP_{key}"] = {c: round(79709 / nS * res[c]["n"] * res[c][key]) for c in CLS}
(OUT / "topologies.json").write_text(json.dumps(res, indent=1)); print(json.dumps(res, indent=1))

# ---- event displays: one example per topology ----
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
pick = {}
def first(cl, cond):
    for x in ex[cl]:
        if cond(x[4]): return x
pick["μ⁻ capture (stopped in the crystal)"] = first("mu- capture", lambda r: r["inCry"] and r["nbMax"] >= 3)
pick["μ⁻ decay in orbit (DIO)"] = first("mu- DIO", lambda r: r["nbMax"] >= 5) or first("mu- DIO", lambda r: True)
pick["μ⁺ decay (Michel electron after µs)"] = first("mu+ decay", lambda r: r["inCry"] and r["late"])
pick["muon stopped outside the crystals"] = first("mu+ decay", lambda r: (not r["inCry"]) and r["late"])
pick["no stop: the muon leaves the disk"] = first("no stop", lambda r: True)
wins = [("prompt −30…+30 ns: the muon", -30, 30), ("30–500 ns: capture / DIO", 30, 500), ("1–20 µs: Michel electron", 1000, 20000)]
fig, ax = plt.subplots(len(pick), 3, figsize=(13, 4.2 * len(pick)))
for i, (title, x) in enumerate(pick.items()):
    if x is None:
        for j in range(3): ax[i, j].axis("off")
        continue
    fn, ev, c0, t0, r = x
    f = ROOT.TFile.Open(fn); h = next(hh for hh in f.Get("hits/hits") if hh.event == ev); hc, ht, he = list(h.cry), list(h.t), list(h.e); f.Close()
    d0, x0, y0 = pos[c0]
    for j, (wl, lo, hi) in enumerate(wins):
        a = ax[i, j]; E = defaultdict(float)
        for cr, t, en in zip(hc, ht, he):
            if lo <= t - t0 < hi and cr in pos and pos[cr][0] == d0: E[cr] += en
        for cr, (dd, xx, yy) in pos.items():
            if dd != d0 or abs(xx - x0) > 260 or abs(yy - y0) > 260: continue
            en = E.get(cr, 0.)
            a.add_patch(Rectangle((xx - 17, yy - 17), 34, 34, facecolor=plt.cm.viridis(min(en / 40, 1)) if en > 0.5 else "white", edgecolor="0.8", lw=0.5))
            if en > 0.5: a.text(xx, yy, f"{en:.0f}", ha="center", va="center", fontsize=6, color="w" if en < 25 else "k")
        a.add_patch(Rectangle((x0 - 17, y0 - 17), 34, 34, fill=False, edgecolor="C3", lw=1.8))
        a.set_xlim(x0 - 260, x0 + 260); a.set_ylim(y0 - 260, y0 + 260); a.set_aspect("equal"); a.set_xticks([]); a.set_yticks([])
        a.set_title(wl if i == 0 else "", fontsize=10)
        if j == 0: a.set_ylabel(title, fontsize=10)
fig.suptitle("Hall Monte Carlo: one event per topology (red square = crystal where the track ends; numbers = MeV)", fontsize=11)
fig.tight_layout(rect=(0, 0, 1, 0.975)); fig.savefig(OUT / "topologies.png", dpi=110)
