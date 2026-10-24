#!/usr/bin/env python3
"""Report of the raw-data check of a run (module 1 of the run validation): reads the histograms of
MC2RawCheck_module.cc and writes one HTML page with, at the top, every check with its criterion and its result
(OK / ATTENZIONE / NO), then the plots.

Input:  the ROOT file of mc2_raw_check.fcl (directory rawCheck/MC2RawCheck); crystal positions from the Offline file
        CalorimeterGeom/data/crystalPos.txt (cvmfs, CRYPOS below). PyROOT needed.
Output: <output dir>/raw_check_<run>.html (in Italian: checks table, rates, list of dead/noisy channels, plots),
        <output dir>/img/*.png, <output dir>/summary.json (parameters, criteria, checks).
        On stdout: "written <html>" and one line per check "<result> <check>: <value>".

Checks (criteria in CRIT, decided before looking at the run; they are printed in the page):
   1 events and rate per time slice          6 ErrorFlags of the hits
   2 event-window tag increasing             7 dead / noisy channels, coincidence with the other SiPM
   3 ROC blocks present                      8 waveform length and saturation, Baseline field
   4 tags of the ROC blocks and of the hits  9 pedestal stability and noise
   5 hits dropped by the Offline decoder    10 board jumps with the laser

usage: python3 Offline/MC2Reco/scripts/raw_check_report.py <raw_check_<run>.root> <run> <output dir>
example: python3 Offline/MC2Reco/scripts/raw_check_report.py out/raw_check_125367.root 125367 out/raw_check_125367
"""
import argparse
import json
import pathlib

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.colors import LogNorm  # noqa: E402
import ROOT  # noqa: E402

CRYPOS = "/cvmfs/mu2e.opensciencegrid.org/Musings/Offline/v13_40_00/Offline/CalorimeterGeom/data/crystalPos.txt"
NB = 160          # boards
NC = 20           # channels per board
NK = 3200         # channel keys, key = board * NC + channel
CRYSTALS_PER_DISK = 674

# ---- criteria (decided before looking at the run; change them here, they are printed in the page) ----
CRIT = dict(bad_hit_ok=1e-3, bad_hit_no=1e-2, flag_ok=1e-3, flag_no=1e-2, tag_shift_ok=1e-4, block_missing_ok=1e-3,
            rate_dev_ok=0.2, dead_min_hits=1, noisy_factor=5.0, partner_min=0.5, sat_frac_max=0.01, ped_shift_ok=2.0,
            ped_shift_no=5.0, ped_sigmas=5.0, jump_samples=1)
# decoder result codes of the Offline (isHitGood), bin index of hErrTot / hErr / hErrSlice
ERRNAME = {0: "Good", 1: "BeginMarker", 2: "LastSampleMarker", 3: "WaveformSize", 4: "NumberOfSamples",
           5: "MaxSampleIndex", 6: "BoardID", 7: "ChannelID", 8: "ErrorFlags", 9: "InvalidChannel"}
RESULT_COLOR = {"OK": "#1d7a4c", "ATTENZIONE": "#9a6a00", "NO": "#c9531f"}
CSS = ("body{margin:0;background:#fbfbf9;color:#16181d;font:15px/1.5 Charter,Georgia,serif}"
       "main{max-width:80rem;margin:0 auto;padding:2rem 1rem}h1,h2{font-family:Helvetica,Arial,sans-serif}"
       "table{border-collapse:collapse;width:100%;font-size:13px}"
       "td,th{border-bottom:1px solid #ddd;padding:4px 6px;vertical-align:top;text-align:left}img{max-width:100%}"
       "figcaption{color:#555;font-size:13px;margin-bottom:1.5rem}"
       "@media (prefers-color-scheme:dark){body{background:#15171b;color:#e8e9ec}td,th{border-color:#333}}")


def channel_name(key):
    """Short name of a channel key: b<board> c<channel>."""
    return f"b{key//NC} c{key%NC}"


def median_bin(counts):
    """Index of the bin where the cumulative of counts reaches half of the total."""
    return np.searchsorted(np.cumsum(counts), counts.sum() / 2)


class RawCheckReport:
    """The histograms of one run, the checks computed from them, the figures and the page."""

    # ---------------- input ----------------

    def __init__(self, root_path, run, out_dir):
        """Open the file, read the parameters of the module, the channel map and the crystal positions."""
        self.root_path = root_path
        self.run = run
        self.out_dir = out_dir
        self.img_dir = out_dir / "img"
        self.img_dir.mkdir(parents=True, exist_ok=True)
        self.root_file = ROOT.TFile.Open(root_path)
        self.directory = self.root_file.Get("rawCheck/MC2RawCheck")
        params = self.hist1d("hParams")
        # event-window length [us] (assumed), time-slice length [min], number of slices
        self.window_us, self.slice_minutes, self.n_slices = params[0], params[1], int(params[2])
        self.events_per_slice = self.hist1d("hEvSlice")
        used = np.nonzero(self.events_per_slice)[0]
        self.n_slices_used = int(used.max()) + 1 if len(used) else 0
        self.slice_centre_min = (np.arange(self.n_slices) + 0.5) * self.slice_minutes
        # offline SiPM id of every channel key (-1: not in the map)
        self.key_offline = self.hist1d("hKeyOffline").astype(int) - 1
        crystal_positions = np.loadtxt(CRYPOS)
        self.crystal_x, self.crystal_y = crystal_positions[:, 2], crystal_positions[:, 3]
        self.checks = []
        self.figures = []

    def hist1d(self, name):
        """Bin contents of a TH1 (no under/overflow)."""
        h = self.directory.Get(name)
        return np.array([h.GetBinContent(i) for i in range(1, h.GetNbinsX() + 1)])

    def hist2d(self, name):
        """Bin contents of a TH2 as [x bin, y bin]."""
        h = self.directory.Get(name)
        return np.array([[h.GetBinContent(i, j) for j in range(1, h.GetNbinsY() + 1)]
                         for i in range(1, h.GetNbinsX() + 1)])

    def hist3d(self, name):
        """Bin contents of a TH3 as [x bin, y bin, z bin]."""
        h = self.directory.Get(name)
        return np.array([[[h.GetBinContent(i, j, k) for k in range(1, h.GetNbinsZ() + 1)]
                          for j in range(1, h.GetNbinsY() + 1)]
                         for i in range(1, h.GetNbinsX() + 1)])

    def has_hist(self, name):
        """True if the histogram is in the file (older versions of the module do not book all of them)."""
        return bool(self.directory.Get(name))

    def crystal_of_key(self, key):
        """Offline crystal of a channel key (-1 if the key is not in the map)."""
        offline_id = self.key_offline[key]
        return (offline_id // 2) if offline_id >= 0 else -1

    def add_check(self, name, value, ok, criterion, warn=None, detail=""):
        """One row of the checks table: OK if ok, otherwise ATTENZIONE if warn, otherwise NO."""
        result = "OK" if ok else ("ATTENZIONE" if warn else "NO")
        self.checks.append(dict(check=name, value=value, criterion=criterion, result=result, detail=detail))

    # ---------------- checks ----------------

    def run_checks(self):
        """All the checks, in the order of the table."""
        self.check_events_and_rate()
        self.check_tag_increasing()
        self.check_roc_blocks()
        self.check_tags()
        self.check_decoder_errors()
        self.check_error_flags()
        self.check_occupancy()
        self.check_waveforms()
        self.check_pedestals()
        self.check_board_jumps()

    def check_events_and_rate(self):
        """1 events and rate: the rate of every time slice (the last one excluded, it is partial) close to the
        median."""
        self.laser_per_slice = self.hist1d("hLaserSlice")
        self.n_events = int(self.events_per_slice.sum())
        tag_range = self.hist1d("hTagRange") if self.has_hist("hTagRange") else np.zeros(2)
        # duration [min] from the first and last event-window tag, otherwise from the slices
        if tag_range[1] > tag_range[0]:
            self.duration_min = (tag_range[1] - tag_range[0]) * self.window_us * 1e-6 / 60
        else:
            self.duration_min = self.n_slices_used * self.slice_minutes
        self.rate_per_slice = self.events_per_slice[:self.n_slices_used] / (self.slice_minutes * 60)
        full = self.rate_per_slice[:-1] if self.n_slices_used > 1 else self.rate_per_slice   # the last slice is partial
        deviation = float(np.max(np.abs(full / np.median(full) - 1))) if len(full) else 0
        self.add_check("eventi e tasso nel tempo",
                       f"{self.n_events} eventi in {self.duration_min:.1f} min (dai tag), "
                       f"{np.median(full):.1f} eventi/s (mediana delle fette)",
                       deviation <= CRIT["rate_dev_ok"],
                       f"tasso di ogni fetta entro ±{100*CRIT['rate_dev_ok']:.0f} % della mediana "
                       f"(ultima fetta esclusa)",
                       warn=True, detail=f"scarto massimo {100*deviation:.0f} %")

    def check_tag_increasing(self):
        """2 event-window tag increasing from one event to the next."""
        tag_step = self.hist1d("hTagStep")
        bad_steps = int(tag_step[:21].sum())   # bins -20..0
        self.add_check("tag dell'event window crescente", f"{bad_steps} passi ≤ 0 su {int(tag_step.sum())}",
                       bad_steps == 0, "nessun evento con tag uguale o minore del precedente")

    def check_roc_blocks(self):
        """3 every ROC block (DTC, link) seen in the run present in (almost) every event."""
        blocks = self.hist2d("hBlocks")
        present = blocks > 0
        missing = [(dtc, link, int(self.n_events - blocks[dtc, link]))
                   for dtc in range(blocks.shape[0]) for link in range(blocks.shape[1])
                   if present[dtc, link] and blocks[dtc, link] < self.n_events * (1 - CRIT["block_missing_ok"])]
        n_blocks = int(present.sum())
        self.add_check("blocchi ROC presenti",
                       f"{n_blocks} blocchi (DTC, link) visti; {len(missing)} mancano in più dello "
                       f"{100*CRIT['block_missing_ok']:.1f} % degli eventi",
                       len(missing) == 0,
                       f"ogni blocco presente in almeno il {100*(1-CRIT['block_missing_ok']):.1f} % degli eventi",
                       detail="; ".join(f"DTC {a} link {b}: manca in {c}" for a, b, c in missing[:10]))

    def check_tags(self):
        """4 event-window tag of the ROC blocks and of the hits equal to the tag of the event (shift bins -4..+4)."""
        block_tag_shift = self.hist2d("hBlockTagShift")
        block_shift_fraction = float(1 - block_tag_shift[:, 4].sum() / max(block_tag_shift.sum(), 1))
        hit_tag_shift = self.hist2d("hHitTagShift")
        hits_per_key = hit_tag_shift.sum(1)
        # the tag offset of the run (most hits): 0 from firmware 1.7, +1 before
        run_offset = int(np.argmax(hit_tag_shift.sum(0))) - 4
        shifted_lanes = [k for k in range(NK)
                         if hits_per_key[k] > 0
                         and 1 - hit_tag_shift[k, run_offset + 4] / hits_per_key[k] > CRIT["tag_shift_ok"]]
        self.add_check("tag dei blocchi ROC = tag dell'evento", f"frazione diversa {block_shift_fraction:.2e}",
                       block_shift_fraction <= CRIT["tag_shift_ok"], f"≤ {CRIT['tag_shift_ok']:.0e}")
        self.add_check("tag degli hit (corsie sfasate)",
                       f"scarto del run: tag dell'hit − tag dell'evento = {run_offset:+d} (dalla maggioranza degli "
                       f"hit); {len(shifted_lanes)} canali con hit a uno scarto diverso",
                       len(shifted_lanes) == 0,
                       f"nessun canale con più di {CRIT['tag_shift_ok']:.0e} degli hit a uno scarto diverso da "
                       f"quello del run",
                       detail=", ".join(f"scheda {k//NC} ch {k%NC}" for k in shifted_lanes[:20]))

    def check_decoder_errors(self):
        """5 hits dropped by the decoder of the Offline (isHitGood), by reason."""
        self.errors_total = self.hist1d("hErrTot")
        n_hits = self.errors_total.sum()
        bad_fraction = float(self.errors_total[1:].sum() / max(n_hits, 1))
        self.errors_per_key = self.hist2d("hErr")
        by_name = {ERRNAME[i]: int(v) for i, v in enumerate(self.errors_total) if i > 0 and v > 0}
        self.add_check("hit scartati dal decoder dell'Offline (isHitGood)",
                       f"{int(self.errors_total[1:].sum())} su {int(n_hits)} = {bad_fraction:.2e}; {by_name}",
                       bad_fraction <= CRIT["bad_hit_ok"],
                       f"OK ≤ {CRIT['bad_hit_ok']:.0e}, NO > {CRIT['bad_hit_no']:.0e}",
                       warn=bad_fraction <= CRIT["bad_hit_no"],
                       detail=f"eventi con almeno un hit scartato: {int(self.hist1d('hEvFailSlice').sum())} "
                              f"({100*self.hist1d('hEvFailSlice').sum()/max(self.n_events,1):.3f} %)")

    def check_error_flags(self):
        """6 ErrorFlags field of the hits (4 bits), for the hits kept and dropped by the decoder."""
        self.flags_total = self.hist1d("hFlagsTot")
        flag_fraction = float(self.flags_total[1:].sum() / max(self.flags_total.sum(), 1))
        kept, dropped = self.hist2d("hFlagBitKept"), self.hist2d("hFlagBitDropped")
        self.add_check("flag di errore degli hit (ErrorFlags, 4 bit)",
                       f"hit con almeno un bit: {int(self.flags_total[1:].sum())} = {flag_fraction:.2e}; "
                       f"per bit (accettati / scartati): "
                       + ", ".join(f"bit {b}: {int(kept[:, b].sum())}/{int(dropped[:, b].sum())}" for b in range(4)),
                       flag_fraction <= CRIT["flag_ok"],
                       f"OK ≤ {CRIT['flag_ok']:.0e}, NO > {CRIT['flag_no']:.0e}",
                       warn=flag_fraction <= CRIT["flag_no"])

    def check_occupancy(self):
        """7 occupancy: dead channels, noisy channels, coincidence with the other SiPM of the crystal."""
        self.occupancy, partner = self.hist1d("hOcc"), self.hist1d("hPartner")
        occupancy = self.occupancy
        self.mapped = self.key_offline >= 0
        mapped = self.mapped
        self.is_lyso = self.hist2d("hDet")[:, 1] > 0   # channels of the LYSO (CAPHRI) crystals
        self.dead = [k for k in range(NK) if mapped[k] and occupancy[k] < CRIT["dead_min_hits"]]
        self.median_occupancy = np.median(occupancy[mapped & (occupancy > 0)])
        median_occupancy = self.median_occupancy
        self.noisy = [k for k in range(NK) if mapped[k] and occupancy[k] > CRIT["noisy_factor"] * median_occupancy]
        # fraction of the hits of a channel with a hit of the other SiPM of the crystal in the same event
        self.partner_fraction = np.where(occupancy > 0, partner / np.maximum(occupancy, 1), 0)
        partner_fraction = self.partner_fraction
        low_partner = [k for k in range(NK)
                       if mapped[k] and occupancy[k] > 100 and partner_fraction[k] < CRIT["partner_min"]]
        self.add_check("canali morti", f"{len(self.dead)} canali della mappa senza hit", len(self.dead) == 0,
                       "nessun canale della mappa senza hit (nei run corti può essere statistica)", warn=True,
                       detail=", ".join(channel_name(k) for k in self.dead[:30]))
        self.add_check("canali rumorosi",
                       f"{len(self.noisy)} canali con più di {CRIT['noisy_factor']:.0f} volte la mediana degli hit",
                       len(self.noisy) == 0, "nessuno", warn=True,
                       detail=", ".join(f"{channel_name(k)} ({occupancy[k]/median_occupancy:.0f}×)"
                                        for k in self.noisy[:20]))
        self.add_check("coincidenza con l'altro SiPM del cristallo",
                       f"mediana {np.median(partner_fraction[mapped & (occupancy > 100)]):.2f}; "
                       f"{len(low_partner)} canali sotto {CRIT['partner_min']}",
                       len(low_partner) == 0,
                       f"ogni canale con > 100 hit ha un hit dell'altro SiPM nello stesso evento in almeno il "
                       f"{100*CRIT['partner_min']:.0f} % dei casi",
                       warn=True,
                       detail=", ".join(f"{channel_name(k)} ({partner_fraction[k]:.2f})" for k in low_partner[:20]))

    def check_waveforms(self):
        """8 waveform length and saturated hits; the Baseline field of the packet (informative)."""
        samples_per_key = self.hist2d("hNSamp")
        samples_edges = np.linspace(0, 1024, 129)
        self.samples_centres = 0.5 * (samples_edges[1:] + samples_edges[:-1])
        self.samples_all = samples_per_key.sum(0)
        samples_median = float(self.samples_centres[median_bin(self.samples_all)])
        saturated = self.hist1d("hSat")
        good_hits = self.errors_per_key[:, 0]
        saturated_fraction = np.where(good_hits > 0, saturated / np.maximum(good_hits, 1), 0)
        high_saturation = [k for k in range(NK)
                           if good_hits[k] > 100 and saturated_fraction[k] > CRIT["sat_frac_max"]]
        self.add_check("forma d'onda: lunghezza e saturazione",
                       f"campioni per hit, mediana ~{samples_median:.0f}; "
                       f"hit saturati {saturated.sum()/max(good_hits.sum(),1):.2e}; {len(high_saturation)} canali "
                       f"con più dell'{100*CRIT['sat_frac_max']:.0f} % saturati",
                       len(high_saturation) == 0,
                       f"nessun canale con più dell'{100*CRIT['sat_frac_max']:.0f} % di hit saturati "
                       f"(con il laser è atteso)",
                       warn=True)
        baseline = self.hist1d("hBaseline")
        self.add_check("campo Baseline del pacchetto",
                       f"valore medio {np.sum(np.arange(len(baseline))*baseline)/max(baseline.sum(),1):.1f} ADC "
                       f"su {int(baseline.sum())} hit",
                       True, "informativo (se è sempre 0 il firmware non lo riempie)")

    def check_pedestals(self):
        """9 pedestals: stability in time (per slice against the run mean) and hit-to-hit noise."""
        ped_sum, ped2_sum, ped_n = self.hist2d("hPedSum"), self.hist2d("hPed2Sum"), self.hist2d("hPedN")
        self.ped_ok = ped_n.sum(1) > 20
        ped_ok = self.ped_ok
        self.ped_run = np.where(ped_ok, ped_sum.sum(1) / np.maximum(ped_n.sum(1), 1), np.nan)
        ped_run = self.ped_run
        self.hit_to_hit = np.where(ped_ok, np.sqrt(np.maximum(ped2_sum.sum(1) / np.maximum(ped_n.sum(1), 1)
                                                              - ped_run ** 2, 0)), np.nan)
        self.ped_slice = np.where(ped_n > 10, ped_sum / np.maximum(ped_n, 1), np.nan)
        # a shift counts only if it is both larger than the criterion and statistically significant: |shift| > 5 sigma,
        # sigma = hit-to-hit rms / sqrt(hits in the slice)
        sigma_slice = self.hit_to_hit[:, None] / np.sqrt(np.maximum(ped_n, 1))
        deviation_slice = np.abs(self.ped_slice - ped_run[:, None])
        significant = np.where(np.isfinite(deviation_slice) & (deviation_slice > CRIT["ped_sigmas"] * sigma_slice),
                               deviation_slice, 0.0)
        self.ped_shift = np.nanmax(significant, axis=1)
        shift = self.ped_shift
        self.ped_moved = [k for k in range(NK) if ped_ok[k] and shift[k] > CRIT["ped_shift_ok"]]
        big_moves = [k for k in self.ped_moved if shift[k] > CRIT["ped_shift_no"]]
        self.add_check("piedistalli stabili nel tempo",
                       f"{len(self.ped_moved)} canali si spostano più di {CRIT['ped_shift_ok']} ADC (e più di "
                       f"{CRIT['ped_sigmas']:.0f}σ) in una fetta di {self.slice_minutes:.0f} min ({len(big_moves)} "
                       f"più di {CRIT['ped_shift_no']} ADC)",
                       len(self.ped_moved) == 0,
                       f"OK: nessuno spostamento oltre ±{CRIT['ped_shift_ok']} ADC che sia anche oltre "
                       f"{CRIT['ped_sigmas']:.0f}σ (σ = rumore / √hit della fetta); NO: qualcuno oltre "
                       f"±{CRIT['ped_shift_no']} ADC",
                       warn=len(big_moves) == 0,
                       detail=", ".join(f"{channel_name(k)} ({shift[k]:.1f})"
                                        for k in sorted(self.ped_moved, key=lambda k: -shift[k])[:20]))
        self.add_check("piedistallo e rumore",
                       f"piedistallo mediano {np.nanmedian(ped_run):.1f} ADC (dispersione fra canali "
                       f"{np.nanstd(ped_run):.1f}); rumore del piedistallo da un hit all'altro, mediana "
                       f"{np.nanmedian(self.hit_to_hit):.2f} ADC",
                       True, "informativo")

    def check_board_jumps(self):
        """10 board jumps with the laser: median laser time of each board per slice against its run median, in
        samples (hBoardJump: [board, slice, value -20..+20])."""
        board_jump = self.hist3d("hBoardJump")
        n_laser = int(self.hist1d("hLaserSlice").sum())
        if n_laser > 0:
            values = np.arange(-20, 21)
            jumps = []
            for board in range(NB):
                total = board_jump[board].sum(1)
                if total.sum() == 0:
                    continue
                median_per_slice = [values[median_bin(board_jump[board, j])] if board_jump[board, j].sum() > 0
                                    else None for j in range(self.n_slices)]
                median_all = values[median_bin(board_jump[board].sum(0))]
                for j, median in enumerate(median_per_slice):
                    if median is not None and abs(median - median_all) >= CRIT["jump_samples"]:
                        jumps.append((board, j, int(median - median_all)))
            self.add_check("salti di scheda con il laser",
                           f"{n_laser} impulsi laser; {len(jumps)} (scheda, fetta) con il tempo mediano spostato di "
                           f"≥ {CRIT['jump_samples']} campione",
                           len(jumps) == 0, "nessuna scheda che salta di un campione intero dentro il run",
                           detail=", ".join(f"scheda {b} a {self.slice_centre_min[j]:.0f} min: {s:+d}"
                                            for b, j, s in jumps[:20]))
        else:
            self.add_check("salti di scheda con il laser", "nessun impulso laser nel run", True,
                           "non applicabile: senza laser i salti si controllano dopo la ricostruzione (modulo 2)")

    # ---------------- figures ----------------

    def save(self, fig, name, caption):
        """Write a figure to img/ and remember it (with its caption) for the page."""
        fig.tight_layout()
        fig.savefig(self.img_dir / name, dpi=110)
        plt.close(fig)
        self.figures.append((name, caption))

    def make_figures(self):
        """All the figures: rates and occupancy first (the most important), then the others."""
        self.prepare_channel_quantities()
        self.plot_rates()
        self.plot_occupancy_maps()
        self.plot_rate_distribution()
        self.make_rate_texts()
        self.plot_time()
        self.plot_errors()
        self.plot_maps()
        self.plot_pedestals()
        self.plot_moved_pedestals()
        self.plot_waveforms()

    def prepare_channel_quantities(self):
        """Duration [s], rate per channel, and crystal, disk and SiPM of every channel key."""
        self.duration_s = max(self.duration_min * 60.0, 1e-9)
        self.rate_channel = self.occupancy / self.duration_s   # Hz per channel, good hits of non-laser events
        self.sensor = np.where(self.key_offline >= 0, self.key_offline % 2, -1)
        self.crystal = np.array([self.crystal_of_key(k) for k in range(NK)])
        self.disk = np.where(self.crystal >= 0, (self.crystal >= CRYSTALS_PER_DISK).astype(int), -1)

    def plot_rates(self):
        """o1: event rate, hit rate per disk, LYSO and dropped hits in time (hRateFine, 10 s bins regrouped to ~50
        points)."""
        fig, ax = plt.subplots(1, 3, figsize=(18, 4.3))
        if self.has_hist("hRateFine"):
            # columns of hRateFine: 0 all events, 1 laser events, 2-3 good hits disk 0-1, 4 LYSO, 5 dropped hits
            rate_fine = self.hist2d("hRateFine")
            duration_s = self.duration_s
            n_fine = int(np.ceil(duration_s / 10.0)) if duration_s > 0 else 1
            n_fine = max(1, min(n_fine, rate_fine.shape[0]))
            group = max(1, int(np.ceil(n_fine / 50.0)))
            n_bins = int(np.ceil(n_fine / group))
            edges = np.arange(n_bins + 1) * group * 10.0
            bin_length = np.minimum(edges[1:], duration_s) - edges[:-1]
            bin_length = np.maximum(bin_length, 1e-3)
            bin_centre_min = (edges[:-1] + 0.5 * bin_length) / 60.0

            def rebin(column):
                """Counts of one column of hRateFine summed in groups of `group` 10-s bins."""
                values = rate_fine[:n_bins * group, column]
                values = np.pad(values, (0, n_bins * group - len(values)))
                return values.reshape(n_bins, group).sum(1)

            def plot_rate(a, column, label, color, scale=1.0):
                """Rate (counts / bin length * scale) with its statistical error; returns the rates."""
                counts = rebin(column)
                rates = counts / bin_length * scale
                errors = np.sqrt(counts) / bin_length * scale
                a.errorbar(bin_centre_min, rates, errors, fmt="o-", ms=3, lw=1, color=color, label=label)
                return rates

            rate_events = plot_rate(ax[0], 0, "all events", "#1f62b8")
            if rate_fine[:, 1].sum() > 0:
                plot_rate(ax[0], 1, "laser events", "#9a6a00")
            ax[0].set_ylim(0, 1.25 * max(np.max(rate_events), 1e-6))
            ax[0].set_ylabel("events / s")
            ax[0].set_title(f"event rate ({group*10:.0f} s bins)", fontsize=10)
            rate_disk0 = plot_rate(ax[1], 2, "disk 0", "#1f62b8", 1e-3)
            rate_disk1 = plot_rate(ax[1], 3, "disk 1", "#c9531f", 1e-3)
            ax[1].set_ylim(0, 1.25 * max(np.max(rate_disk0), np.max(rate_disk1), 1e-6))
            ax[1].set_ylabel("good hits / s [kHz]")
            ax[1].set_title("hit rate of the non-laser events, per disk", fontsize=10)
            rate_lyso = plot_rate(ax[2], 4, "LYSO (good hits, non-laser)", "#1d7a4c")
            rate_dropped = plot_rate(ax[2], 5, "dropped hits (all events)", "#c9531f")
            ax[2].set_ylim(0, 1.25 * max(np.max(rate_lyso), np.max(rate_dropped), 1e-6))
            ax[2].set_ylabel("hits / s")
            ax[2].set_title("LYSO hits and hits dropped by the decoder", fontsize=10)
            for a in ax:
                a.set_xlabel("time from the first event [min]")
                a.grid(alpha=.3)
                a.legend(fontsize=8)
        self.save(fig, "o1_rates.png",
                  "Tassi nel tempo, con l'errore statistico e la scala da 0: eventi al secondo (e impulsi laser), hit "
                  "buoni per disco negli eventi senza laser, hit dei LYSO e hit scartati dal decoder. Gli intervalli "
                  "si adattano alla durata del run (~50 punti).")

    def plot_occupancy_maps(self):
        """o2: hit rate per channel on the two disks, for SiPM 0 and SiPM 1, with dead and noisy channels marked."""
        fig, ax = plt.subplots(2, 2, figsize=(14, 11.5))
        rate_channel = self.rate_channel
        positive = rate_channel[self.mapped & (rate_channel > 0)]
        if len(positive):
            vmin, vmax = (np.percentile(positive, 2), np.percentile(positive, 98))
        else:
            vmin, vmax = (1e-3, 1)
        cx, cy, crystal = self.crystal_x, self.crystal_y, self.crystal
        for i, disk in enumerate((0, 1)):
            for j, sipm in enumerate((0, 1)):
                a = ax[i, j]
                keys = [k for k in range(NK)
                        if crystal[k] >= 0 and crystal[k] < len(cx) and self.disk[k] == disk and self.sensor[k] == sipm]
                xs, ys, vs = cx[crystal[keys]], cy[crystal[keys]], np.maximum(rate_channel[keys], vmin * 0.5)
                sc = a.scatter(xs, ys, c=vs, s=22, marker="s", cmap="viridis", norm=LogNorm(vmin=vmin, vmax=vmax))
                dead = [k for k in keys if self.occupancy[k] < CRIT["dead_min_hits"]]
                a.scatter(cx[crystal[dead]], cy[crystal[dead]], s=40, marker="x", color="red",
                          label=f"no hits ({len(dead)})")
                noisy = [k for k in keys if k in self.noisy]
                a.scatter(cx[crystal[noisy]], cy[crystal[noisy]], s=60, facecolors="none", edgecolors="magenta",
                          label=f"> {CRIT['noisy_factor']:.0f}× median ({len(noisy)})")
                a.set_aspect("equal")
                a.set_title(f"disk {disk}, SiPM {sipm}: hit rate per channel", fontsize=10)
                a.set_xlabel("x [mm]")
                a.set_ylabel("y [mm]")
                a.legend(fontsize=8, loc="lower right")
                plt.colorbar(sc, ax=a, fraction=0.046, label="good hits / s (non-laser events)")
        self.save(fig, "o2_occupancy_maps.png",
                  "Occupazione: tasso di hit buoni per canale [Hz], negli eventi senza laser, sui due dischi e per i "
                  "due SiPM di ogni cristallo (posizioni dall'Offline); croci rosse: canali senza hit; cerchi: canali "
                  "rumorosi.")

    def plot_rate_distribution(self):
        """o3: distribution of the channel rates and rate per board."""
        fig, ax = plt.subplots(1, 2, figsize=(15, 4.3))
        median_rate = self.median_occupancy / self.duration_s
        log_rate = np.log10(np.maximum(self.rate_channel[self.mapped], 1e-4))
        ax[0].hist(log_rate, 100, color="#1f62b8")
        ax[0].axvline(np.log10(max(self.median_occupancy / self.duration_s, 1e-4)), color="k", ls="--",
                      label=f"median {median_rate:.3g} Hz")
        ax[0].axvline(np.log10(max(CRIT['noisy_factor'] * self.median_occupancy / self.duration_s, 1e-4)),
                      color="magenta", ls=":", label=f"noisy > {CRIT['noisy_factor']:.0f}× median")
        ax[0].set_xlabel("log10(hit rate per channel / Hz)")
        ax[0].set_ylabel("channels")
        ax[0].set_title("distribution of the channel rates (channels with 0 hits at the left edge)", fontsize=10)
        ax[0].legend(fontsize=8)
        ax[0].grid(alpha=.3)
        rate_board = self.rate_channel.reshape(NB, NC).sum(1)
        mapped_per_board = self.mapped.reshape(NB, NC).sum(1)
        ax[1].bar(np.arange(NB)[mapped_per_board > 0], rate_board[mapped_per_board > 0], color="#1d7a4c")
        ax[1].set_xlabel("board")
        ax[1].set_ylabel("hit rate [Hz] (sum of its channels)")
        ax[1].set_title("hit rate per board", fontsize=10)
        ax[1].grid(alpha=.3)
        self.save(fig, "o3_rate_distribution.png", "Distribuzione dei tassi dei canali e tasso per scheda.")

    def channel_row(self, key):
        """One row of the table of dead / noisy channels."""
        c = self.crystal[key]
        relative = self.occupancy[key] / self.median_occupancy if self.median_occupancy else 0
        return (f"<tr><td>{key//NC}</td><td>{key%NC}</td><td>{self.key_offline[key]}</td><td>{c}</td>"
                f"<td>{self.disk[key]}</td><td>{self.sensor[key]}</td><td>{'LYSO' if self.is_lyso[key] else 'CsI'}"
                f"</td><td class=\"num\">{self.rate_channel[key]:.3g}</td><td class=\"num\">{relative:.2f}</td></tr>")

    def make_rate_texts(self):
        """Table of the dead and noisy channels (sorted by rate) and the sentence summarising the rates."""
        self.channel_table = ("<table><thead><tr><th>scheda</th><th>canale</th><th>SiPM offline</th><th>cristallo</th>"
                              "<th>disco</th><th>SiPM</th><th>tipo</th><th>tasso [Hz]</th><th>/ mediana</th></tr>"
                              "</thead><tbody>"
                              + "".join(self.channel_row(k)
                                        for k in sorted(set(self.dead) | set(self.noisy),
                                                        key=lambda k: -self.rate_channel[k]))
                              + "</tbody></table>")
        csi_hits = self.mapped & ~self.is_lyso
        self.rate_summary = (f"Durata {self.duration_min:.1f} min (dai tag dell'event window). Eventi: "
                             f"{self.n_events/self.duration_s:.1f} /s ({int(self.laser_per_slice.sum())} laser). "
                             f"Hit buoni negli eventi senza laser: "
                             f"{self.rate_channel[csi_hits & (self.disk == 0)].sum()/1e3:.2f} kHz sul disco 0, "
                             f"{self.rate_channel[csi_hits & (self.disk == 1)].sum()/1e3:.2f} kHz sul disco 1; "
                             f"tasso mediano di un canale {self.median_occupancy/self.duration_s:.3g} Hz.")

    def plot_time(self):
        """r1: event (and laser) rate per slice, dropped hits by reason and ErrorFlags bits in time."""
        n_used = self.n_slices_used
        t_slice = self.slice_centre_min[:n_used]
        fig, ax = plt.subplots(1, 3, figsize=(16, 4))
        ax[0].step(t_slice, self.rate_per_slice, where="mid")
        ax[0].set_xlabel(f"time [min] (slices of {self.slice_minutes:.0f} min)")
        ax[0].set_ylabel("events / s")
        ax[0].set_title("event rate", fontsize=10)
        ax[0].step(t_slice, self.laser_per_slice[:n_used] / (self.slice_minutes * 60), where="mid",
                   label="laser events / s")
        ax[0].legend(fontsize=8)
        errors_slice = self.hist2d("hErrSlice")
        for code in range(1, 10):
            if errors_slice[:, code].sum() > 0:
                ax[1].step(t_slice, errors_slice[:n_used, code] / np.maximum(errors_slice[:n_used].sum(1), 1),
                           where="mid", label=ERRNAME[code])
        ax[1].set_yscale("log")
        ax[1].set_xlabel("time [min]")
        ax[1].set_ylabel("fraction of the hits")
        ax[1].set_title("hits dropped by the Offline decoder, by reason", fontsize=10)
        if errors_slice[:, 1:].sum() > 0:
            ax[1].legend(fontsize=8)
        flag_bits_slice = self.hist2d("hFlagBitSlice")
        for bit in range(4):
            ax[2].step(t_slice, flag_bits_slice[:n_used, bit], where="mid", label=f"bit {bit}")
        ax[2].set_xlabel("time [min]")
        ax[2].set_ylabel("hits with the bit set")
        ax[2].set_title("ErrorFlags of the hits, bit by bit", fontsize=10)
        ax[2].legend(fontsize=8)
        for a in ax:
            a.grid(alpha=.3)
        self.save(fig, "r1_time.png",
                  "Nel tempo: tasso degli eventi (e degli impulsi laser), hit scartati dal decoder per motivo, flag di "
                  "errore degli hit bit per bit.")

    def plot_errors(self):
        """r2: dropped hits by reason and per board, ErrorFlags values."""
        fig, ax = plt.subplots(1, 3, figsize=(16, 4.4))
        ax[0].bar(range(1, 10), self.errors_total[1:], color="#c9531f")
        ax[0].set_xticks(range(1, 10))
        ax[0].set_xticklabels([ERRNAME[i] for i in range(1, 10)], rotation=45, ha="right", fontsize=8)
        ax[0].set_ylabel("hits")
        ax[0].set_title("dropped hits by reason", fontsize=10)
        dropped_per_board = self.errors_per_key[:, 1:].sum(1).reshape(NB, NC).sum(1)
        ax[1].bar(range(NB), dropped_per_board, color="#c9531f")
        ax[1].set_xlabel("board")
        ax[1].set_ylabel("dropped hits")
        ax[1].set_title("dropped hits per board", fontsize=10)
        ax[2].bar(range(16), self.flags_total, color="#1f62b8")
        ax[2].set_yscale("log")
        ax[2].set_xlabel("ErrorFlags value")
        ax[2].set_ylabel("hits")
        ax[2].set_title("ErrorFlags values (all hits)", fontsize=10)
        for a in ax:
            a.grid(alpha=.3)
        self.save(fig, "r2_errors.png",
                  "Hit scartati dal decoder dell'Offline per motivo e per scheda; valori del campo ErrorFlags.")

    def disk_map(self, ax, values, title, cmap="viridis", vmin=None, vmax=None, label=""):
        """Map of a per-key quantity on the two disks (ax[0], ax[1]), one point per crystal (SiPM 0)."""
        for disk in (0, 1):
            a = ax[disk]
            xs, ys, vs = [], [], []
            for k in range(NK):
                c = self.crystal_of_key(k)
                # one point per crystal (SiPM 0); offline ids beyond the crystals (e.g. PIN) skipped
                if (c < 0 or c >= len(self.crystal_x) or (c >= CRYSTALS_PER_DISK) != bool(disk)
                        or self.key_offline[k] % 2):
                    continue
                xs.append(self.crystal_x[c])
                ys.append(self.crystal_y[c])
                vs.append(values[k])
            sc = a.scatter(xs, ys, c=vs, s=14, marker="s", cmap=cmap, vmin=vmin, vmax=vmax)
            a.set_aspect("equal")
            a.set_title(f"{title}, disk {disk}", fontsize=10)
            a.set_xlabel("x [mm]")
            a.set_ylabel("y [mm]")
            plt.colorbar(sc, ax=a, fraction=0.046, label=label)

    def plot_maps(self):
        """r3: occupancy relative to the median and partner coincidence on the two disks."""
        fig, ax = plt.subplots(1, 4, figsize=(20, 4.6))
        relative = np.where(self.median_occupancy > 0, self.occupancy / self.median_occupancy, 0)
        self.disk_map(ax[:2], relative, "hits / median (SiPM 0)", vmin=0, vmax=2, label="occupancy / median")
        self.disk_map(ax[2:], self.partner_fraction, "partner coincidence (SiPM 0)", vmin=0, vmax=1, label="fraction")
        self.save(fig, "r3_maps.png",
                  "Mappe dei due dischi (posizioni dei cristalli dall'Offline): numero di hit rispetto alla mediana, e "
                  "frazione di hit con un hit dell'altro SiPM del cristallo nello stesso evento.")

    def plot_pedestals(self):
        """r4: pedestal and its noise per channel, median pedestal shift of each board in time."""
        ped_ok = self.ped_ok
        fig, ax = plt.subplots(1, 3, figsize=(16, 4.2))
        ax[0].hist(self.ped_run[ped_ok], 80, color="#1f62b8")
        ax[0].set_xlabel("pedestal [ADC]")
        ax[0].set_ylabel("channels")
        ax[0].set_title("pedestal per channel", fontsize=10)
        ax[1].hist(self.hit_to_hit[ped_ok], 80, color="#1d7a4c")
        ax[1].set_xlabel("pedestal rms from hit to hit [ADC]")
        ax[1].set_ylabel("channels")
        ax[1].set_title("noise of the pedestal", fontsize=10)
        board_shift = np.full((NB, self.n_slices), np.nan)
        for board in range(NB):
            keys = [k for k in range(board * NC, (board + 1) * NC) if ped_ok[k]]
            if keys:
                board_shift[board] = np.nanmedian(self.ped_slice[keys] - self.ped_run[keys][:, None], axis=0)
        im = ax[2].imshow(board_shift[:, :self.n_slices_used], aspect="auto", origin="lower", cmap="RdBu_r",
                          vmin=-2, vmax=2, extent=[0, self.n_slices_used * self.slice_minutes, -0.5, NB - 0.5])
        plt.colorbar(im, ax=ax[2], label="median pedestal shift [ADC]")
        ax[2].set_xlabel("time [min]")
        ax[2].set_ylabel("board")
        ax[2].set_title("pedestal stability: median shift per board", fontsize=10)
        for a in ax[:2]:
            a.grid(alpha=.3)
        self.save(fig, "r4_pedestals.png",
                  "Piedistallo di ogni canale, suo rumore da un hit all'altro, e stabilità nel tempo (spostamento "
                  "mediano dei canali di ogni scheda rispetto alla media del run).")

    def plot_moved_pedestals(self):
        """r5 (only if some pedestal moves): pedestal in time of the 8 channels that move most."""
        if not self.ped_moved:
            return
        n_used = self.n_slices_used
        fig, ax = plt.subplots(figsize=(9, 4))
        for k in sorted(self.ped_moved, key=lambda k: -self.ped_shift[k])[:8]:
            ax.plot(self.slice_centre_min[:n_used], self.ped_slice[k, :n_used] - self.ped_run[k], "o-", ms=3,
                    label=f"board {k//NC} ch {k%NC}")
        ax.axhline(0, color="k", lw=.8)
        ax.set_xlabel("time [min]")
        ax.set_ylabel("pedestal - run mean [ADC]")
        ax.grid(alpha=.3)
        ax.legend(fontsize=8)
        ax.set_title("channels whose pedestal moves most", fontsize=10)
        self.save(fig, "r5_ped_moved.png", "I canali il cui piedistallo si sposta di più nel tempo.")

    def plot_waveforms(self):
        """r6: waveform length, position of the largest sample, hit tag relative to the event tag."""
        fig, ax = plt.subplots(1, 3, figsize=(16, 4))
        ax[0].step(self.samples_centres, self.samples_all, where="mid")
        ax[0].set_yscale("log")
        ax[0].set_xlabel("samples in the waveform")
        ax[0].set_ylabel("hits")
        ax[0].set_title("waveform length", fontsize=10)
        peak_index = self.hist1d("hPeakIdx")
        ax[1].step(range(len(peak_index)), peak_index, where="mid")
        ax[1].set_yscale("log")
        ax[1].set_xlabel("index of the largest sample")
        ax[1].set_title("position of the maximum", fontsize=10)
        hit_tag_shift = self.hist2d("hHitTagShift").sum(0)
        ax[2].bar(range(-4, 5), hit_tag_shift, color="#9a6a00")
        ax[2].set_yscale("log")
        ax[2].set_xlabel("hit tag - event tag (clipped at ±4)")
        ax[2].set_title("event-window tag of the hits", fontsize=10)
        for a in ax:
            a.grid(alpha=.3)
        self.save(fig, "r6_waveforms.png",
                  "Lunghezza delle forme d'onda, posizione del massimo, tag degli hit rispetto all'evento.")

    # ---------------- outputs ----------------

    def write_summary(self):
        """summary.json: run, file, parameters, criteria, checks, number of events."""
        summary = dict(run=self.run, file=self.root_path,
                       params=dict(windowMicroseconds=self.window_us, sliceMinutes=self.slice_minutes),
                       criteria=CRIT, checks=self.checks, n_events=self.n_events)
        (self.out_dir / "summary.json").write_text(json.dumps(summary, indent=1, default=float))

    def figures_html(self, figures):
        """<figure> elements of a list of (file name, caption)."""
        return ''.join(f'<figure><img src="img/{n}" alt=""><figcaption>{c}</figcaption></figure>' for n, c in figures)

    def write_html(self):
        """raw_check_<run>.html: checks table, rates and occupancy, dead / noisy channels, other plots."""
        rows = "".join(f'<tr><td>{c["check"]}</td>'
                       f'<td style="color:{RESULT_COLOR[c["result"]]};font-weight:bold">{c["result"]}</td>'
                       f'<td>{c["value"]}</td><td>{c["criterion"]}</td><td>{c["detail"]}</td></tr>'
                       for c in self.checks)
        run = self.run
        html = ('<!doctype html><html lang="it"><head><meta charset="utf-8">'
                '<meta name="viewport" content="width=device-width, initial-scale=1">'
                f'<title>Controllo dati grezzi {run}</title><style>{CSS}</style></head><body><main>\n'
                f'<h1>Controllo dei dati grezzi, run {run}</h1>\n'
                "<p>Modulo art <code>MC2RawCheck</code> (mc2reco, Offline v13_40_00) sul file raw, prima della "
                "ricostruzione; legge i frammenti DTC come il decoder dell'Offline e tiene anche gli hit che il "
                f"decoder scarta. File: <code>{self.root_path}</code>. Asse dei tempi: {self.window_us:.0f} µs per "
                f"event window (ipotesi), fette di {self.slice_minutes:.0f} min. Criteri decisi prima di guardare il "
                "run (in testa allo script).</p>\n"
                '<h2>Controlli</h2><table><thead><tr><th>controllo</th><th>esito</th><th>valore</th>'
                f'<th>criterio</th><th>dettaglio</th></tr></thead><tbody>{rows}</tbody></table>\n'
                f'<h2>Tassi e occupazione</h2><p>{self.rate_summary}</p>{self.figures_html(self.figures[:3])}\n'
                f'<h3>Canali senza hit o rumorosi</h3>{self.channel_table}\n'
                f'<h2>Altri plot</h2>{self.figures_html(self.figures[3:])}\n'
                '</main></body></html>')
        html_path = self.out_dir / f"raw_check_{run}.html"
        html_path.write_text(html, encoding="utf-8")
        print("written", html_path)

    def print_checks(self):
        """One line per check on stdout."""
        for c in self.checks:
            print(f"{c['result']:10s} {c['check']}: {c['value']}")


def parse_args():
    """Command line: ROOT file of MC2RawCheck, run, output directory."""
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("root_file", help="raw_check_<run>.root (histograms of MC2RawCheck)")
    parser.add_argument("run", help="run number (used in the page title and file name)")
    parser.add_argument("out", type=pathlib.Path, help="output directory")
    return parser.parse_args()


def main():
    args = parse_args()
    report = RawCheckReport(args.root_file, args.run, args.out)
    report.run_checks()
    report.make_figures()
    report.write_summary()
    report.write_html()
    report.print_checks()


if __name__ == "__main__":
    main()
