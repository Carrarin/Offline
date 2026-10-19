#!/usr/bin/env python3
"""Hall Monte Carlo (Mu2eG4 + CRY, digitized and reconstructed with mc2_digi_reco.fcl) against the cosmic data of run 125329:
track rates and observed decays of stopped muons, with no free parameter (normalisation = CRY livetime).

MC   : digi_*.root of run_digi_reco.sh: dump/tracks (MC2CaloCosmicTrackDump; type 1, 2 one-disk, 3 two-disk) and stopped/pairs
       (MC2StoppedMuons). Livetime: sum of "Cosmic livetime" of the Geant4 job logs.
Decays observed: near pairs (cls 1) with nD <= 6 and aD <= 3000 ADC, as book chapter 3 / 12; binned Poisson likelihood fit of
       N0 exp(-t/tau) + B on dt in [1, 20] us, 0.25 us bins; observed decays = integral of the exponential on [1, 20] us.
Data : 125329 from mc2-ana (output/book/ch05, ch03 results.json, output/cosmic_tracks_run125329/summary.txt): tracks of the mc2tools
       tagger and pairs of cosmic_delayed_coincidence.C (same selection; the data macro also asks |tL - tR| < 10 ns per crystal).
usage: PYTHONPATH=$ROOTSYS/lib python3 Offline/MC2Reco/scripts/sim_vs_data.py <production dir> <digi dir> <mc2-ana dir> <out dir>
"""
import glob, json, math, pathlib, re, sys
import numpy as np
import ROOT
ROOT.gROOT.SetBatch(True); ROOT.gErrorIgnoreLevel = ROOT.kError
PROD, DIGI, ANA, OUT = map(pathlib.Path, sys.argv[1:5]); OUT.mkdir(parents=True, exist_ok=True)

live = sum(float(m.group(1)) for f in PROD.glob("job_*.log") for m in re.finditer(r"Cosmic livetime: ([0-9.eE+-]+) s", f.read_text()))
files = sorted(glob.glob(str(DIGI / "digi_*.root")))
ntype = {1: 0, 2: 0, 3: 0}; dts = []
for fn in files:
    f = ROOT.TFile.Open(fn)
    t = f.Get("dump/tracks")
    for e in t:
        ntype[e.type] = ntype.get(e.type, 0) + 1
    p = f.Get("stopped/pairs")
    for e in p:
        if e.cls == 1 and e.nD <= 6 and e.aD <= 3000: dts.append(e.dt)
    f.Close()
dts = np.array(dts)

TMIN, TMAX, TBIN = 1.0, 20.0, 0.25
h = ROOT.TH1D("h", "", int((TMAX - TMIN) / TBIN), TMIN, TMAX); h.SetDirectory(0)
for x in dts[(dts >= TMIN) & (dts < TMAX)]: h.Fill(x)
fx = ROOT.TF1("fx", "[0]*exp(-x/[1])+[2]", TMIN, TMAX); fx.SetParameters(max(1, h.GetBinContent(1)) * 1.5, 2.2, 0.1); fx.SetParLimits(1, 0.3, 20); fx.SetParLimits(2, 0, 1e6)
r = h.Fit(fx, "LSQN0R")
N0, tau, B = (fx.GetParameter(i) for i in range(3)); eN0, etau = fx.GetParError(0), fx.GetParError(1)
sig = N0 * tau / TBIN * (math.exp(-TMIN / tau) - math.exp(-TMAX / tau))
cov = r.GetCovarianceMatrix()
# error on sig from N0 and tau (numerical derivatives)
def S(n0, ta): return n0 * ta / TBIN * (math.exp(-TMIN / ta) - math.exp(-TMAX / ta))
dN = S(N0 * 1.001, tau) / sig - 1; dT = (S(N0, tau * 1.001) - sig) / (sig * 0.001 * tau)
g = np.array([sig / N0, dT * sig, 0.0]); C = np.array([[cov(i, j) for j in range(3)] for i in range(3)]); esig = float(np.sqrt(g @ C @ g))

# data 125329
c5 = json.loads((ANA / "output/book/ch05/results.json").read_text()); c3 = json.loads((ANA / "output/book/ch03/results.json").read_text())
ts = (ANA / "output/cosmic_tracks_run125329/summary.txt").read_text()
m = re.search(r"disk 0 only (\d+), disk 1 only (\d+), both disks (\d+)", ts); T = c5["ch05.run_s"]; eps_laser = 1 - c5["ch05.dead_laser_frac"]
data = {"one_disk_Hz": (int(m.group(1)) + int(m.group(2))) / T, "two_disk_Hz": int(m.group(3)) / T,
        "decays_Hz": c3["ch03.main.signal"] / c3["ch03.run_s"], "decays_err_Hz": c3["ch03.main.signal"] ** 0.5 / c3["ch03.run_s"], "run_s": T}
n1 = ntype[1] + ntype[2]
mc = {"livetime_s": live, "files": len(files), "one_disk_Hz": n1 / live * eps_laser, "one_disk_err_Hz": n1 ** 0.5 / live,
      "two_disk_Hz": ntype[3] / live * eps_laser, "two_disk_err_Hz": ntype[3] ** 0.5 / live, "tracks": ntype,
      "near_pairs_1_20us": int(h.Integral()), "tau_us": tau, "tau_err_us": etau, "B_per_bin": B,
      "decays_Hz": sig / live * eps_laser, "decays_err_Hz": esig / live, "decays_counted": sig}
res = {"mc": mc, "data": data, "ratio_data_mc": {k: data[k] / mc[k] for k in ("one_disk_Hz", "two_disk_Hz", "decays_Hz")},
       "decays_per_1000_onedisk": {"mc": 1000 * mc["decays_Hz"] / mc["one_disk_Hz"], "data": 1000 * data["decays_Hz"] / data["one_disk_Hz"]}}
(OUT / "sim_vs_data.json").write_text(json.dumps(res, indent=1)); print(json.dumps(res, indent=1))

import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
fig, ax = plt.subplots(1, 2, figsize=(15, 4.6))
c = np.array([h.GetBinCenter(i) for i in range(1, h.GetNbinsX() + 1)]); y = np.array([h.GetBinContent(i) for i in range(1, h.GetNbinsX() + 1)])
ax[0].errorbar(c, y, np.sqrt(np.clip(y, 1, None)), fmt="o", ms=4, color="k", label=f"hall MC, near pairs ({live:.0f} s live)")
ax[0].plot(c, [fx.Eval(x) for x in c], color="C0", lw=0.9, label=f"fit: τ = {tau:.2f} ± {etau:.2f} µs")
ax[0].set_yscale("log"); ax[0].set_xlabel("delay of the delayed deposit [µs]"); ax[0].set_ylabel("pairs / 0.25 µs"); ax[0].legend(fontsize=7.5)
ax[0].set_title("stopped muons in the hall Monte Carlo", fontsize=10)
lab = ["one-disk tracks", "two-disk tracks", "decays 1–20 µs"]; km = ["one_disk_Hz", "two_disk_Hz", "decays_Hz"]
mv = [mc[k] for k in km]; me = [mc["one_disk_err_Hz"], mc["two_disk_err_Hz"], mc["decays_err_Hz"]]; dv = [data[k] for k in km]
xx = np.arange(3)
ax[1].errorbar(xx - 0.1, np.array(dv) / np.array(mv), np.array(me) / np.array(mv) * np.array(dv) / np.array(mv), fmt="o", ms=5, color="C3", label="data 125329 / hall MC")
ax[1].axhline(1, color="0.45", lw=0.9); ax[1].set_xticks(xx); ax[1].set_xticklabels(lab); ax[1].set_ylabel("data / MC"); ax[1].legend(fontsize=7.5)
ax[1].set_ylim(0, 1.6); ax[1].set_title("no free parameter: normalisation = CRY livetime", fontsize=10)
fig.suptitle("Hall Monte Carlo (Mu2eG4 + CRY) against run 125329", fontsize=11); fig.tight_layout(); fig.savefig(OUT / "sim_vs_data.png", dpi=125)
