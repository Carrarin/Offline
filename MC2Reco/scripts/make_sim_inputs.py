#!/usr/bin/env python3
"""Detector response of the MC2 calorimeter for the Offline digitization (CaloShowerROMaker + Offline CaloDigiMaker with the threshold per SiPM), from measurements
that do not use the Michel electrons (19/10/2026):

  CalSimCrystals_<run>_laser.txt  TABLE CalSimCrystals (crid,LRU,pePerMeV0,pePerMeV1,ADCPerMeV0,ADCPerMeV1):
      ADCPerMeV of each SiPM = 1 / adc2mev of the MIP calibration of the run (mc2-ana energy_calibration, fitted amplitudes,
      one-disk tracks); pePerMeV = ADCPerMeV / p0 of the crystal from the laser intensity scan (laser_lr_photostat, run 125328,
      book chapter 9; median for crystals without a fit). SiPMs without a MIP fit: median adc2mev. LRU = 0 (the measured 5.8 %
      is not linear and its sign in the Offline z is not checked).
  CalSimCrystals_<run>_laser_fcl.txt  the same, in the format of CalSimParams.fileName (useDb false).
  CalEnergyCalib_<run>_sim.txt    TABLE CalEnergyCalib, ADC2MeV = 1/ADCPerMeV of the simulation: the calibration of the simulated data.
  sim_thresholds_<run>.csv        roid,thrRaw_adc: threshold of each SiPM above the pedestal, 0.5 % percentile of the hit
      amplitude (largest sample - pedestal), as book chapter 10, on the first nThr hits of cosmic_hit_features_<run>.
  pulse_cosmic_r125328_v2.root    h_waveform: the cosmic template v2 (gtempl of splines_cosmic_r125328_v2.root, firmware 1.8)
      on 0.05 ns bins from where it rises above 0.1 % of the peak; CaloPulseUtil normalises it to its maximum.

usage: PYTHONPATH=$ROOTSYS/lib python3 Offline/MC2Reco/scripts/make_sim_inputs.py <mc2-ana dir> <run, e.g. r125329> <out dir> [nThr]
"""
import csv, json, pathlib, sys
import numpy as np
import ROOT
ROOT.gROOT.SetBatch(True)

ANA, RUN, OUT = pathlib.Path(sys.argv[1]), sys.argv[2], pathlib.Path(sys.argv[3]); NTHR = int(sys.argv[4]) if len(sys.argv) > 4 else 30_000_000
OUT.mkdir(parents=True, exist_ok=True); here = pathlib.Path(__file__).resolve().parents[1]
mip = list(csv.DictReader(open(next((ANA / "output").glob(f"energy_calibration_cal_{RUN}__*/sipm_mip.csv")))))
key2roid = {int(r["board"]) * 20 + int(r["chan"]): int(r["sipmid"]) for r in mip}
a2m = {int(r["sipmid"]): float(r["adc2mev"]) for r in mip if int(r["fit_ok"]) == 1 and float(r["adc2mev"]) > 0}
a2m_med = float(np.median(list(a2m.values())))
cry_of = {int(r["sipmid"]): int(r["cry"]) for r in mip}

# p0 per crystal from the laser scan
lp = list(csv.DictReader(open(ANA / "output" / "laser_lr_photostat_run125328" / "per_crystal.csv")))
p0 = {}
for r in lp:
    k = int(r["keyL"])
    if k in key2roid and float(r["p0"]) > 0: p0[cry_of[key2roid[k]]] = float(r["p0"])
p0_med = float(np.median(list(p0.values())))
NCRY = 1348
lines = ["TABLE CalSimCrystals",
         f"#crid,LRU,pePerMeV0,pePerMeV1,ADCPerMeV0,ADCPerMeV1  make_sim_inputs.py: ADCPerMeV = 1/adc2mev (MIP, {RUN}; median {1/a2m_med:.2f} where not fitted); "
         f"pePerMeV = ADCPerMeV / p0 (laser 125328, median p0 {p0_med:.3f} ADC/pe where not fitted); LRU = 0"]
npe = []
for c in range(NCRY):
    adc = [1.0 / a2m.get(2 * c + s, a2m_med) for s in (0, 1)]; pp = p0.get(c, p0_med); pe = [x / pp for x in adc]; npe += pe
    lines.append(f"{c},0.0000,{pe[0]:.3f},{pe[1]:.3f},{adc[0]:.4f},{adc[1]:.4f}")
(OUT / f"CalSimCrystals_{RUN}_laser.txt").write_text("\n".join(lines) + "\n")
# the same table in the CalSimParams fcl format (CalSimParams.fileName, useDb false: "crid LRU npe0 npe1 adc0 adc1", no header)
(OUT / f"CalSimCrystals_{RUN}_laser_fcl.txt").write_text("".join(" ".join(l.split(",")) + "\n" for l in lines[2:]))
# CalEnergyCalib for the reconstruction of the simulation: ADC2MeV = 1/ADCPerMeV of the simulation, SiPM by SiPM
e = ["TABLE CalEnergyCalib", f"#roid,ADC2MeV  make_sim_inputs.py: 1/ADCPerMeV of CalSimCrystals_{RUN}_laser.txt"]
for l in lines[2:]:
    c, _, _, _, a0, a1 = l.split(","); e += [f"{2 * int(c)},{1 / float(a0):.6f}", f"{2 * int(c) + 1},{1 / float(a1):.6f}"]
e += [f"{r},{a2m_med:.6f}" for r in range(2 * NCRY, 2720)]   # the 24 channels beyond the crystals (as the machine tables)
(OUT / f"CalEnergyCalib_{RUN}_sim.txt").write_text("\n".join(e) + "\n")

# thresholds per SiPM (book chapter 10)
f = str(next((ANA / "output").glob(f"cosmic_hit_features_{RUN}/features.root")))
d = ROOT.RDataFrame("hits", f).Range(NTHR).Filter("type==0").AsNumpy(["key", "amp"])
key, amp = d["key"].astype(int), d["amp"].astype(float)
order = np.argsort(key, kind="stable"); ks = key[order]; st = np.r_[0, np.flatnonzero(np.diff(ks)) + 1]; en = np.r_[st[1:], len(ks)]
thr = {}
for s, e in zip(st, en):
    if e - s < 200 or ks[s] not in key2roid: continue
    thr[key2roid[ks[s]]] = float(np.percentile(amp[order[s:e]], 0.5))
with open(OUT / f"sim_thresholds_{RUN}.csv", "w") as fo:
    fo.write("roid,thrRaw_adc\n"); [fo.write(f"{r},{thr[r]:.1f}\n") for r in sorted(thr)]

# pulse shape
tf = ROOT.TFile.Open(str(here / "data" / "templates" / "splines_cosmic_r125328_v2.root")); g = tf.Get("gtempl")
x = np.array([g.GetX()[i] for i in range(g.GetN())]); y = np.array([g.GetY()[i] for i in range(g.GetN())])
ipk = int(np.argmax(y)); i0 = ipk
while i0 > 0 and y[i0 - 1] > 1e-3 * y[ipk]: i0 -= 1
x0, x1 = x[i0], x[-1]; nb = int(round((x1 - x0) / 0.05))
h = ROOT.TH1F("h_waveform", "cosmic template v2 (r125328), t from the start of the rise;t [ns];shape", nb, 0.0, x1 - x0)
for i in range(1, nb + 1): h.SetBinContent(i, max(0.0, g.Eval(x0 + h.GetBinCenter(i))))
fo = ROOT.TFile(str(OUT / "pulse_cosmic_r125328_v2.root"), "RECREATE"); h.Write(); fo.Close()

s = {"run": RUN, "sipm_mip_fitted": len(a2m), "ADCPerMeV_median": 1 / a2m_med, "crystals_p0_laser": len(p0), "p0_median": p0_med,
     "pePerMeV_median": float(np.median(npe)), "pePerMeV_16_84": [float(np.percentile(npe, 16)), float(np.percentile(npe, 84))],
     "thresholds": len(thr), "thr_median_adc": float(np.median(list(thr.values()))), "thr_16_84": [float(np.percentile(list(thr.values()), q)) for q in (16, 84)],
     "hits_for_thresholds": int(len(key)), "pulse": {"start_ns_from_peak": float(x0), "end_ns": float(x1), "peak_ns_from_start": float(-x0), "bins": nb}}
(OUT / f"sim_inputs_{RUN}.json").write_text(json.dumps(s, indent=1)); print(json.dumps(s, indent=1))
