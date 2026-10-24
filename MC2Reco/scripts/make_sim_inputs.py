#!/usr/bin/env python3
"""Detector response of the MC2 calorimeter for the Offline digitization (CaloShowerROMaker + Offline CaloDigiMaker
with the threshold per SiPM), from measurements that do not use the Michel electrons (19/10/2026).

Inputs (from the mc2-ana area <mc2-ana dir>/output):
  energy_calibration_cal_<run>__*/sipm_mip.csv   MIP calibration of the run (board, chan, sipmid, cry, adc2mev, fit_ok)
  laser_lr_photostat_run125328/per_crystal.csv   photostatistics p0 [ADC/pe] per crystal from the laser intensity scan
  cosmic_hit_features_<run>/features.root        tree "hits" (key, amp, type) of the cosmic hits of the run
  and, from this package, data/templates/splines_cosmic_r125328_v2.root (graph gtempl, the cosmic template v2).

Outputs (in <out dir>):
  CalSimCrystals_<run>_laser.txt  TABLE CalSimCrystals (crid,LRU,pePerMeV0,pePerMeV1,ADCPerMeV0,ADCPerMeV1):
      ADCPerMeV of each SiPM = 1 / adc2mev of the MIP calibration of the run (mc2-ana energy_calibration, fitted
      amplitudes, one-disk tracks); pePerMeV = ADCPerMeV / p0 of the crystal from the laser intensity scan
      (laser_lr_photostat, run 125328, book chapter 9; median for crystals without a fit). SiPMs without a MIP fit:
      median adc2mev. LRU = 0 (the measured 5.8 % is not linear and its sign in the Offline z is not checked).
  CalSimCrystals_<run>_laser_fcl.txt  the same, in the format of CalSimParams.fileName (useDb false).
  CalEnergyCalib_<run>_sim.txt    TABLE CalEnergyCalib, ADC2MeV = 1/ADCPerMeV of the simulation: the calibration of
      the simulated data.
  sim_thresholds_<run>.csv        roid,thrRaw_adc: threshold of each SiPM above the pedestal, 0.5 % percentile of the
      hit amplitude (largest sample - pedestal), as book chapter 10, on the first nThr hits of
      cosmic_hit_features_<run>.
  pulse_cosmic_r125328_v2.root    h_waveform: the cosmic template v2 (gtempl of splines_cosmic_r125328_v2.root,
      firmware 1.8) on 0.05 ns bins from where it rises above 0.1 % of the peak; CaloPulseUtil normalises it to its
      maximum.
  sim_inputs_<run>.json           summary (also printed).

usage: PYTHONPATH=$ROOTSYS/lib python3 Offline/MC2Reco/scripts/make_sim_inputs.py <mc2-ana dir> <run, e.g. r125329>
           <out dir> [nThr]
example: PYTHONPATH=$ROOTSYS/lib python3 Offline/MC2Reco/scripts/make_sim_inputs.py ~/laser_analysis/mc2-ana r125329 \\
           /data1/mu2e/mc2/work/sim_inputs
"""
import argparse
import csv
import json
import pathlib

import numpy as np
import ROOT

NCRY = 1348                    # crystals of the two disks
NROID = 2720                   # readout ids of the machine tables (2 * NCRY SiPMs + 24 channels beyond the crystals)
DEFAULT_NTHR = 30_000_000      # hits of cosmic_hit_features used for the thresholds
THR_PERCENTILE = 0.5           # threshold = this percentile [%] of the hit amplitude of the SiPM
THR_MIN_HITS = 200             # SiPMs with fewer hits get no threshold
PULSE_BIN_NS = 0.05            # bin width of h_waveform [ns]
PULSE_START_FRACTION = 1e-3    # the pulse starts where the template rises above this fraction of the peak


def parse_args():
    """Command line: mc2-ana area, run, output directory, number of hits for the thresholds."""
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("ana", type=pathlib.Path, help="mc2-ana area (with output/)")
    parser.add_argument("run", help="run label as in the mc2-ana outputs, e.g. r125329")
    parser.add_argument("out", type=pathlib.Path, help="output directory")
    parser.add_argument("nthr", nargs="?", type=int, default=DEFAULT_NTHR,
                        help="hits of cosmic_hit_features used for the thresholds (default %(default)d)")
    return parser.parse_args()


def read_mip_calibration(ana, run):
    """MIP calibration of the run: (board*20+chan -> roid, roid -> adc2mev of the good fits, median adc2mev,
    roid -> crystal)."""
    with open(next((ana / "output").glob(f"energy_calibration_cal_{run}__*/sipm_mip.csv"))) as csv_file:
        mip = list(csv.DictReader(csv_file))
    key_to_roid = {int(r["board"]) * 20 + int(r["chan"]): int(r["sipmid"]) for r in mip}
    adc2mev = {int(r["sipmid"]): float(r["adc2mev"]) for r in mip if int(r["fit_ok"]) == 1 and float(r["adc2mev"]) > 0}
    adc2mev_median = float(np.median(list(adc2mev.values())))
    crystal_of_roid = {int(r["sipmid"]): int(r["cry"]) for r in mip}
    return key_to_roid, adc2mev, adc2mev_median, crystal_of_roid


def read_laser_p0(ana, key_to_roid, crystal_of_roid):
    """p0 [ADC/pe] per crystal from the laser scan of run 125328 (keyed by the left SiPM), and its median."""
    with open(ana / "output" / "laser_lr_photostat_run125328" / "per_crystal.csv") as csv_file:
        laser_rows = list(csv.DictReader(csv_file))
    p0 = {}
    for row in laser_rows:
        key = int(row["keyL"])
        if key in key_to_roid and float(row["p0"]) > 0:
            p0[crystal_of_roid[key_to_roid[key]]] = float(row["p0"])
    p0_median = float(np.median(list(p0.values())))
    return p0, p0_median


def sim_crystals_lines(run, adc2mev, adc2mev_median, p0, p0_median):
    """Lines of TABLE CalSimCrystals (two header lines, then one per crystal) and the pePerMeV of every SiPM."""
    lines = ["TABLE CalSimCrystals",
             f"#crid,LRU,pePerMeV0,pePerMeV1,ADCPerMeV0,ADCPerMeV1  make_sim_inputs.py: ADCPerMeV = 1/adc2mev "
             f"(MIP, {run}; median {1/adc2mev_median:.2f} where not fitted); "
             f"pePerMeV = ADCPerMeV / p0 (laser 125328, median p0 {p0_median:.3f} ADC/pe where not fitted); LRU = 0"]
    pe_per_mev_all = []
    for crystal in range(NCRY):
        # roid = 2 * crystal + SiPM index
        adc_per_mev = [1.0 / adc2mev.get(2 * crystal + sipm, adc2mev_median) for sipm in (0, 1)]
        crystal_p0 = p0.get(crystal, p0_median)
        pe_per_mev = [x / crystal_p0 for x in adc_per_mev]
        pe_per_mev_all += pe_per_mev
        lines.append(f"{crystal},0.0000,{pe_per_mev[0]:.3f},{pe_per_mev[1]:.3f},"
                     f"{adc_per_mev[0]:.4f},{adc_per_mev[1]:.4f}")
    return lines, pe_per_mev_all


def write_sim_crystals(out, run, lines):
    """CalSimCrystals as DbService table and in the CalSimParams fcl format."""
    (out / f"CalSimCrystals_{run}_laser.txt").write_text("\n".join(lines) + "\n")
    # the same table in the CalSimParams fcl format (CalSimParams.fileName, useDb false:
    # "crid LRU npe0 npe1 adc0 adc1", no header)
    (out / f"CalSimCrystals_{run}_laser_fcl.txt").write_text("".join(" ".join(line.split(",")) + "\n"
                                                                     for line in lines[2:]))


def write_sim_energy_calib(out, run, lines, adc2mev_median):
    """CalEnergyCalib for the reconstruction of the simulation: ADC2MeV = 1/ADCPerMeV of the simulation, SiPM by
    SiPM."""
    table = ["TABLE CalEnergyCalib",
             f"#roid,ADC2MeV  make_sim_inputs.py: 1/ADCPerMeV of CalSimCrystals_{run}_laser.txt"]
    for line in lines[2:]:
        crystal, _, _, _, adc_per_mev0, adc_per_mev1 = line.split(",")
        table += [f"{2 * int(crystal)},{1 / float(adc_per_mev0):.6f}",
                  f"{2 * int(crystal) + 1},{1 / float(adc_per_mev1):.6f}"]
    # the 24 channels beyond the crystals (as the machine tables)
    table += [f"{roid},{adc2mev_median:.6f}" for roid in range(2 * NCRY, NROID)]
    (out / f"CalEnergyCalib_{run}_sim.txt").write_text("\n".join(table) + "\n")


def sipm_thresholds(ana, run, nthr, key_to_roid):
    """Threshold per SiPM (book chapter 10): 0.5 % percentile of the amplitude of the cosmic hits (type 0) among the
    first nthr hits. Returns (roid -> threshold [ADC], number of hits read)."""
    features_file = str(next((ana / "output").glob(f"cosmic_hit_features_{run}/features.root")))
    columns = ROOT.RDataFrame("hits", features_file).Range(nthr).Filter("type==0").AsNumpy(["key", "amp"])
    key, amp = columns["key"].astype(int), columns["amp"].astype(float)
    # group the hits by key: sort, then [start, end) of each run of equal keys
    order = np.argsort(key, kind="stable")
    sorted_keys = key[order]
    starts = np.r_[0, np.flatnonzero(np.diff(sorted_keys)) + 1]
    ends = np.r_[starts[1:], len(sorted_keys)]
    thresholds = {}
    for start, end in zip(starts, ends):
        if end - start < THR_MIN_HITS or sorted_keys[start] not in key_to_roid:
            continue
        thresholds[key_to_roid[sorted_keys[start]]] = float(np.percentile(amp[order[start:end]], THR_PERCENTILE))
    return thresholds, len(key)


def write_thresholds(out, run, thresholds):
    """sim_thresholds_<run>.csv: roid,thrRaw_adc."""
    with open(out / f"sim_thresholds_{run}.csv", "w") as csv_out:
        csv_out.write("roid,thrRaw_adc\n")
        for roid in sorted(thresholds):
            csv_out.write(f"{roid},{thresholds[roid]:.1f}\n")


def write_pulse_shape(out, package_dir):
    """h_waveform from the cosmic template v2: from where it rises above 0.1 % of the peak to its end, 0.05 ns bins.
    Returns (start of the pulse w.r.t. the template time [ns], end [ns], number of bins)."""
    template_file = ROOT.TFile.Open(str(package_dir / "data" / "templates" / "splines_cosmic_r125328_v2.root"))
    template = template_file.Get("gtempl")
    x = np.array([template.GetX()[i] for i in range(template.GetN())])
    y = np.array([template.GetY()[i] for i in range(template.GetN())])
    i_peak = int(np.argmax(y))
    i_start = i_peak
    while i_start > 0 and y[i_start - 1] > PULSE_START_FRACTION * y[i_peak]:
        i_start -= 1
    x_start, x_end = x[i_start], x[-1]
    nbins = int(round((x_end - x_start) / PULSE_BIN_NS))
    # the histogram is booked while the template file is the current directory; it is written to the output file
    hist = ROOT.TH1F("h_waveform", "cosmic template v2 (r125328), t from the start of the rise;t [ns];shape",
                     nbins, 0.0, x_end - x_start)
    for i in range(1, nbins + 1):
        hist.SetBinContent(i, max(0.0, template.Eval(x_start + hist.GetBinCenter(i))))
    pulse_file = ROOT.TFile(str(out / "pulse_cosmic_r125328_v2.root"), "RECREATE")
    hist.Write()
    pulse_file.Close()
    return x_start, x_end, nbins


def main():
    args = parse_args()
    ROOT.gROOT.SetBatch(True)
    args.out.mkdir(parents=True, exist_ok=True)
    package_dir = pathlib.Path(__file__).resolve().parents[1]   # MC2Reco

    key_to_roid, adc2mev, adc2mev_median, crystal_of_roid = read_mip_calibration(args.ana, args.run)
    p0, p0_median = read_laser_p0(args.ana, key_to_roid, crystal_of_roid)

    lines, pe_per_mev_all = sim_crystals_lines(args.run, adc2mev, adc2mev_median, p0, p0_median)
    write_sim_crystals(args.out, args.run, lines)
    write_sim_energy_calib(args.out, args.run, lines, adc2mev_median)

    thresholds, n_hits = sipm_thresholds(args.ana, args.run, args.nthr, key_to_roid)
    write_thresholds(args.out, args.run, thresholds)

    x_start, x_end, nbins = write_pulse_shape(args.out, package_dir)

    summary = {"run": args.run, "sipm_mip_fitted": len(adc2mev), "ADCPerMeV_median": 1 / adc2mev_median,
               "crystals_p0_laser": len(p0), "p0_median": p0_median,
               "pePerMeV_median": float(np.median(pe_per_mev_all)),
               "pePerMeV_16_84": [float(np.percentile(pe_per_mev_all, 16)), float(np.percentile(pe_per_mev_all, 84))],
               "thresholds": len(thresholds), "thr_median_adc": float(np.median(list(thresholds.values()))),
               "thr_16_84": [float(np.percentile(list(thresholds.values()), q)) for q in (16, 84)],
               "hits_for_thresholds": int(n_hits),
               "pulse": {"start_ns_from_peak": float(x_start), "end_ns": float(x_end),
                         "peak_ns_from_start": float(-x_start), "bins": nbins}}
    (args.out / f"sim_inputs_{args.run}.json").write_text(json.dumps(summary, indent=1))
    print(json.dumps(summary, indent=1))


if __name__ == "__main__":
    main()
