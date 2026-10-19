#!/usr/bin/env python3
"""Hall MC (mc2_digi_reco_capture.fcl) with the cuts of the 125367 data tables of the notebook, STOP candidates split by truth class.
Second pulse: maxS > 5, A2 / 23.4 ADC/MeV in bins; clusters (capture/cluster) with E >= 10 MeV in 30-300 ns and 1-10 us.
Scaled to 79709 STOP and 157770 THROUGH. usage: PYTHONPATH=$ROOTSYS/lib python3 ... <digicap dir> <DMAP.txt> <out json>"""
import glob, json, math, sys
from collections import Counter, defaultdict
import ROOT
D, DMAP, OUT = sys.argv[1:4]
DIO, NCAP = 114, 115
pos = {}
for ln in open(DMAP).read().splitlines()[1:]:
    p = ln.split()
    if len(p) > 17 and p[17] == "CAL" and p[2] in ("0", "1"): pos[int(p[14])] = (int(p[2]), float(p[12]), float(p[13]))
def near(a, b):
    A, B = pos.get(a), pos.get(b); return bool(A and B and A[0] == B[0] and math.hypot(A[1] - B[1], A[2] - B[2]) <= 50)
CL = ["mu- capture", "mu- DIO", "mu+ decay", "no stop"]
cnt = {0: Counter(), 1: Counter()}; ncand = Counter()
for fn in sorted(glob.glob(D + "/digi_*.root")):
    f = ROOT.TFile.Open(fn); tr = defaultdict(list)
    for e in f.Get("truth/stops"):
        dc = list(e.dCode)
        cls = ("mu- DIO" if DIO in dc else ("mu- capture" if NCAP in dc else None)) if e.q < 0 else "mu+ decay"
        if cls: tr[e.event].append((cls, e.cry))
    lab = {}
    for c in f.Get("capture/capture"):
        l = next((cls for cls, cry in tr.get(c.event, []) if cry == c.cry or near(cry, c.cry)), "no stop")
        lab[(c.event, c.cry, c.kind)] = l; ncand[c.kind] += 1
        if c.maxS > 5:
            m = c.A2 / 23.4
            b = "5-10" if 5 <= m < 10 else ("10-20" if 10 <= m < 20 else ("20-40" if 20 <= m < 40 else (">40" if m >= 40 else None)))
            if b: cnt[c.kind][(b, l)] += 1
    seen = set()
    for c in f.Get("capture/cluster"):
        if c.E < 10: continue
        w = "cl>=10 30-300ns" if 30 <= c.dt <= 300 else ("cl>=10 1-10us" if 1000 <= c.dt <= 10000 else None)
        k = (c.event, c.cry, c.kind, w)
        if w is None or k in seen: continue
        seen.add(k); cnt[c.kind][(w, lab.get((c.event, c.cry, c.kind), "no stop"))] += 1
    f.Close()
S, T = 79709 / ncand[0], 157770 / ncand[1]
data = {"5-10": [277, 54], "10-20": [94, 2], "20-40": [12, 0], ">40": [0, 0], "cl>=10 30-300ns": [85, 11], "cl>=10 1-10us": [623, 69]}
res = {}
print("%-16s %8s %7s %6s %6s %7s %7s | %s" % ("observable", "MC STOP", "capture", "DIO", "mu+", "no stop", "MC THR", "data STOP, THROUGH"))
for b in data:
    parts = [cnt[0][(b, l)] * S for l in CL]; thr = sum(cnt[1][(b, l)] for l in CL) * T
    res[b] = {"mc_stop": round(sum(parts)), "by_class": dict(zip(CL, [round(x) for x in parts])), "mc_through": round(thr), "data_stop_through": data[b]}
    print("%-16s %8.0f %7.0f %6.0f %6.0f %7.0f %7.0f | %s" % (b, sum(parts), *parts, thr, data[b]))
json.dump(res, open(OUT, "w"), indent=1)
