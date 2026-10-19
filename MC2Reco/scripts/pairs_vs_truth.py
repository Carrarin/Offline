#!/usr/bin/env python3
"""Each reconstructed near pair (MC2StoppedMuons, cuts of book chapter 3: cls 1, nD <= 6, aD <= 3000 ADC, dt 1-20 us) of the hall
Monte Carlo matched to the truth of the same event (MC2StopTruth): which stopped muon, if any, has a decay / capture at that delay.
usage: PYTHONPATH=$ROOTSYS/lib python3 Offline/MC2Reco/scripts/pairs_vs_truth.py <digi dir> <truth dir> <out json>"""
import glob, json, pathlib, sys
from collections import Counter, defaultdict
import numpy as np
import ROOT
DIGI, TRU, OUT = pathlib.Path(sys.argv[1]), pathlib.Path(sys.argv[2]), pathlib.Path(sys.argv[3])
DECAY, DIO, NCAP = 14, 114, 115
cls = Counter(); per_event = Counter(); dts = defaultdict(list)
for fd in sorted(glob.glob(str(DIGI / "digi_*.root"))):
    n = fd.split("_")[-1].split(".")[0]; ft = TRU / f"truth_{n}.root"
    tr = defaultdict(list); f = ROOT.TFile.Open(str(ft))
    for e in f.Get("truth/stops"):
        dc, dt = list(e.dCode), list(e.dT)
        tdec = min([t for c, t in zip(dc, dt) if c in (DECAY, DIO)], default=-1.0); tcap = min([t for c, t in zip(dc, dt) if c == NCAP], default=-1.0)
        tr[e.event].append((e.q, bool(e.inCry), tdec, tcap, e.eTot))
    f.Close()
    f = ROOT.TFile.Open(fd)
    for p in f.Get("stopped/pairs"):
        if not (p.cls == 1 and p.nD <= 6 and p.aD <= 3000 and 1 <= p.dt < 20): continue
        ev = int(p.entry); per_event[(n, ev)] += 1; best = "no stopped muon in the event"
        cands = tr.get(ev, [])
        for q, inc, tdec, tcap, eT in cands:
            if tdec >= 0 and abs(tdec / 1000 - p.dt) < 0.2: best = ("mu+ decay" if q > 0 else "mu- DIO") + (" (in crystal)" if inc else " (outside crystals)"); break
        else:
            if cands: best = "stopped muon(s), no decay at that delay" + (" [mu- capture in event]" if any(c[0] < 0 and c[3] >= 0 for c in cands) else "")
        cls[best] += 1; dts[best].append(p.dt)
    f.Close()
tot = sum(cls.values()); mult = Counter(min(v, 4) for v in per_event.values())
R = {"pairs": tot, "classes": {k: [v, round(v / tot, 3), round(float(np.mean(dts[k])), 2)] for k, v in cls.most_common()},
     "pairs_per_event": {str(k): v for k, v in sorted(mult.items())}}
OUT.write_text(json.dumps(R, indent=1)); print(json.dumps(R, indent=1))
