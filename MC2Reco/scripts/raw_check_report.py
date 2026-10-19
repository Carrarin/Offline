#!/usr/bin/env python3
"""Report of the raw-data check of a run (module 1 of the run validation): reads the histograms of MC2RawCheck_module.cc and writes
one HTML page with, at the top, every check with its criterion and its result (OK / ATTENZIONE / NO), then the plots.

usage: python3 Offline/MC2Reco/scripts/raw_check_report.py <raw_check_<run>.root> <run> <output dir>
       (PyROOT needed; crystal positions from the Offline file CalorimeterGeom/data/crystalPos.txt)
Output: <output dir>/raw_check_<run>.html, img/*.png, summary.json
"""
import json, pathlib, sys
import numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
import ROOT

ROOTF, RUN, OUT = sys.argv[1], sys.argv[2], pathlib.Path(sys.argv[3]); IMG = OUT / "img"; IMG.mkdir(parents=True, exist_ok=True)
CRYPOS = "/cvmfs/mu2e.opensciencegrid.org/Musings/Offline/v13_40_00/Offline/CalorimeterGeom/data/crystalPos.txt"
NB, NC, NK = 160, 20, 3200
# ---- criteria (decided before looking at the run; change them here, they are printed in the page) ----
CRIT = dict(bad_hit_ok=1e-3, bad_hit_no=1e-2, flag_ok=1e-3, flag_no=1e-2, tag_shift_ok=1e-4, block_missing_ok=1e-3, rate_dev_ok=0.2,
            dead_min_hits=1, noisy_factor=5.0, partner_min=0.5, sat_frac_max=0.01, ped_shift_ok=2.0, ped_shift_no=5.0, ped_sigmas=5.0, jump_samples=1)
ERRNAME = {0: "Good", 1: "BeginMarker", 2: "LastSampleMarker", 3: "WaveformSize", 4: "NumberOfSamples", 5: "MaxSampleIndex", 6: "BoardID", 7: "ChannelID", 8: "ErrorFlags", 9: "InvalidChannel"}

f = ROOT.TFile.Open(ROOTF); D = f.Get("rawCheck/MC2RawCheck")
def h1(n): h = D.Get(n); return np.array([h.GetBinContent(i) for i in range(1, h.GetNbinsX() + 1)])
def h2(n): h = D.Get(n); return np.array([[h.GetBinContent(i, j) for j in range(1, h.GetNbinsY() + 1)] for i in range(1, h.GetNbinsX() + 1)])
def h3(n):
    h = D.Get(n); return np.array([[[h.GetBinContent(i, j, k) for k in range(1, h.GetNbinsZ() + 1)] for j in range(1, h.GetNbinsY() + 1)] for i in range(1, h.GetNbinsX() + 1)])
P = h1("hParams"); WIN, SLM, NSL = P[0], P[1], int(P[2])
ev = h1("hEvSlice"); used = np.nonzero(ev)[0]; nsl_used = int(used.max()) + 1 if len(used) else 0; tsl = (np.arange(NSL) + 0.5) * SLM
keyoff = h1("hKeyOffline").astype(int) - 1
cp = np.loadtxt(CRYPOS); cx, cy = cp[:, 2], cp[:, 3]
def cry_of(key): o = keyoff[key]; return (o // 2) if o >= 0 else -1
checks = []
def check(name, value, ok, crit, warn=None, detail=""):
    res = "OK" if ok else ("ATTENZIONE" if warn else "NO"); checks.append(dict(check=name, value=value, criterion=crit, result=res, detail=detail))

# 1 events and rate
la_ = h1("hLaserSlice"); nev = int(ev.sum()); tr = h1("hTagRange") if D.Get("hTagRange") else np.zeros(2); dur = (tr[1] - tr[0]) * WIN * 1e-6 / 60 if tr[1] > tr[0] else nsl_used * SLM; rate = ev[:nsl_used] / (SLM * 60)
full = rate[:-1] if nsl_used > 1 else rate   # the last slice is partial
dev = float(np.max(np.abs(full / np.median(full) - 1))) if len(full) else 0
check("eventi e tasso nel tempo", f"{nev} eventi in {dur:.1f} min (dai tag), {np.median(full):.1f} eventi/s (mediana delle fette)", dev <= CRIT["rate_dev_ok"], f"tasso di ogni fetta entro ±{100*CRIT['rate_dev_ok']:.0f} % della mediana (ultima fetta esclusa)", warn=True, detail=f"scarto massimo {100*dev:.0f} %")
# 2 tag monotonic
ts = h1("hTagStep"); bad_steps = int(ts[:21].sum())   # bins -20..0
check("tag dell'event window crescente", f"{bad_steps} passi ≤ 0 su {int(ts.sum())}", bad_steps == 0, "nessun evento con tag uguale o minore del precedente")
# 3 ROC blocks
bl = h2("hBlocks"); present = bl > 0; miss = [(d_, l_, int(nev - bl[d_, l_])) for d_ in range(bl.shape[0]) for l_ in range(bl.shape[1]) if present[d_, l_] and bl[d_, l_] < nev * (1 - CRIT["block_missing_ok"])]
nbl = int(present.sum())
check("blocchi ROC presenti", f"{nbl} blocchi (DTC, link) visti; {len(miss)} mancano in più dello {100*CRIT['block_missing_ok']:.1f} % degli eventi", len(miss) == 0, f"ogni blocco presente in almeno il {100*(1-CRIT['block_missing_ok']):.1f} % degli eventi", detail="; ".join(f"DTC {a} link {b}: manca in {c}" for a, b, c in miss[:10]))
# 4 tags
bts = h2("hBlockTagShift"); bfrac = float(1 - bts[:, 4].sum() / max(bts.sum(), 1))
hts = h2("hHitTagShift"); hn = hts.sum(1); run_off = int(np.argmax(hts.sum(0))) - 4   # the tag offset of the run (most hits): 0 from firmware 1.7, +1 before
lanes = [k for k in range(NK) if hn[k] > 0 and 1 - hts[k, run_off + 4] / hn[k] > CRIT["tag_shift_ok"]]
check("tag dei blocchi ROC = tag dell'evento", f"frazione diversa {bfrac:.2e}", bfrac <= CRIT["tag_shift_ok"], f"≤ {CRIT['tag_shift_ok']:.0e}")
check("tag degli hit (corsie sfasate)", f"scarto del run: tag dell'hit − tag dell'evento = {run_off:+d} (dalla maggioranza degli hit); {len(lanes)} canali con hit a uno scarto diverso", len(lanes) == 0, f"nessun canale con più di {CRIT['tag_shift_ok']:.0e} degli hit a uno scarto diverso da quello del run",
      detail=", ".join(f"scheda {k//NC} ch {k%NC}" for k in lanes[:20]))
# 5 decoder errors
err = h1("hErrTot"); ntot = err.sum(); fbad = float(err[1:].sum() / max(ntot, 1)); eperk = h2("hErr")
byname = {ERRNAME[i]: int(v) for i, v in enumerate(err) if i > 0 and v > 0}
check("hit scartati dal decoder dell'Offline (isHitGood)", f"{int(err[1:].sum())} su {int(ntot)} = {fbad:.2e}; {byname}", fbad <= CRIT["bad_hit_ok"], f"OK ≤ {CRIT['bad_hit_ok']:.0e}, NO > {CRIT['bad_hit_no']:.0e}", warn=fbad <= CRIT["bad_hit_no"],
      detail=f"eventi con almeno un hit scartato: {int(h1('hEvFailSlice').sum())} ({100*h1('hEvFailSlice').sum()/max(nev,1):.3f} %)")
# 6 ErrorFlags
fl = h1("hFlagsTot"); fflag = float(fl[1:].sum() / max(fl.sum(), 1)); kept, drop = h2("hFlagBitKept"), h2("hFlagBitDropped")
check("flag di errore degli hit (ErrorFlags, 4 bit)", f"hit con almeno un bit: {int(fl[1:].sum())} = {fflag:.2e}; per bit (accettati / scartati): " + ", ".join(f"bit {b}: {int(kept[:, b].sum())}/{int(drop[:, b].sum())}" for b in range(4)),
      fflag <= CRIT["flag_ok"], f"OK ≤ {CRIT['flag_ok']:.0e}, NO > {CRIT['flag_no']:.0e}", warn=fflag <= CRIT["flag_no"])
# 7 occupancy
occ, part = h1("hOcc"), h1("hPartner"); mapped = keyoff >= 0; caphri = h2("hDet")[:, 1] > 0
dead = [k for k in range(NK) if mapped[k] and occ[k] < CRIT["dead_min_hits"]]
medocc = np.median(occ[mapped & (occ > 0)]); noisy = [k for k in range(NK) if mapped[k] and occ[k] > CRIT["noisy_factor"] * medocc]
pf = np.where(occ > 0, part / np.maximum(occ, 1), 0); lowp = [k for k in range(NK) if mapped[k] and occ[k] > 100 and pf[k] < CRIT["partner_min"]]
check("canali morti", f"{len(dead)} canali della mappa senza hit", len(dead) == 0, "nessun canale della mappa senza hit (nei run corti può essere statistica)", warn=True, detail=", ".join(f"b{k//NC} c{k%NC}" for k in dead[:30]))
check("canali rumorosi", f"{len(noisy)} canali con più di {CRIT['noisy_factor']:.0f} volte la mediana degli hit", len(noisy) == 0, "nessuno", warn=True, detail=", ".join(f"b{k//NC} c{k%NC} ({occ[k]/medocc:.0f}×)" for k in noisy[:20]))
check("coincidenza con l'altro SiPM del cristallo", f"mediana {np.median(pf[mapped & (occ > 100)]):.2f}; {len(lowp)} canali sotto {CRIT['partner_min']}", len(lowp) == 0, f"ogni canale con > 100 hit ha un hit dell'altro SiPM nello stesso evento in almeno il {100*CRIT['partner_min']:.0f} % dei casi", warn=True,
      detail=", ".join(f"b{k//NC} c{k%NC} ({pf[k]:.2f})" for k in lowp[:20]))
# 8 waveform
ns = h2("hNSamp"); nsx = np.linspace(0, 1024, 129); nsc = 0.5 * (nsx[1:] + nsx[:-1]); nsall = ns.sum(0); nsmed = float(nsc[np.searchsorted(np.cumsum(nsall), nsall.sum() / 2)])
sat = h1("hSat"); goodk = eperk[:, 0]; satf = np.where(goodk > 0, sat / np.maximum(goodk, 1), 0); hisat = [k for k in range(NK) if goodk[k] > 100 and satf[k] > CRIT["sat_frac_max"]]
check("forma d'onda: lunghezza e saturazione", f"campioni per hit, mediana ~{nsmed:.0f}; hit saturati {sat.sum()/max(goodk.sum(),1):.2e}; {len(hisat)} canali con più dell'{100*CRIT['sat_frac_max']:.0f} % saturati", len(hisat) == 0, f"nessun canale con più dell'{100*CRIT['sat_frac_max']:.0f} % di hit saturati (con il laser è atteso)", warn=True)
base = h1("hBaseline"); check("campo Baseline del pacchetto", f"valore medio {np.sum(np.arange(len(base))*base)/max(base.sum(),1):.1f} ADC su {int(base.sum())} hit", True, "informativo (se è sempre 0 il firmware non lo riempie)")
# 9 pedestals
ps, p2, pn, nz = h2("hPedSum"), h2("hPed2Sum"), h2("hPedN"), h2("hNoiseSum")
ok = pn.sum(1) > 20; pedall = np.where(ok, ps.sum(1) / np.maximum(pn.sum(1), 1), np.nan)
hit2hit = np.where(ok, np.sqrt(np.maximum(p2.sum(1) / np.maximum(pn.sum(1), 1) - pedall ** 2, 0)), np.nan)
pedsl = np.where(pn > 10, ps / np.maximum(pn, 1), np.nan)
# a shift counts only if it is both larger than the criterion and statistically significant: |shift| > 5 sigma, sigma = hit-to-hit rms / sqrt(hits in the slice)
sig_sl = hit2hit[:, None] / np.sqrt(np.maximum(pn, 1)); dev_sl = np.abs(pedsl - pedall[:, None]); signif = np.where(np.isfinite(dev_sl) & (dev_sl > CRIT["ped_sigmas"] * sig_sl), dev_sl, 0.0)
shift = np.nanmax(signif, axis=1)
moved = [k for k in range(NK) if ok[k] and shift[k] > CRIT["ped_shift_ok"]]; bigmov = [k for k in moved if shift[k] > CRIT["ped_shift_no"]]
check("piedistalli stabili nel tempo", f"{len(moved)} canali si spostano più di {CRIT['ped_shift_ok']} ADC (e più di {CRIT['ped_sigmas']:.0f}σ) in una fetta di {SLM:.0f} min ({len(bigmov)} più di {CRIT['ped_shift_no']} ADC)", len(moved) == 0,
      f"OK: nessuno spostamento oltre ±{CRIT['ped_shift_ok']} ADC che sia anche oltre {CRIT['ped_sigmas']:.0f}σ (σ = rumore / √hit della fetta); NO: qualcuno oltre ±{CRIT['ped_shift_no']} ADC", warn=len(bigmov) == 0, detail=", ".join(f"b{k//NC} c{k%NC} ({shift[k]:.1f})" for k in sorted(moved, key=lambda k: -shift[k])[:20]))
check("piedistallo e rumore", f"piedistallo mediano {np.nanmedian(pedall):.1f} ADC (dispersione fra canali {np.nanstd(pedall):.1f}); rumore del piedistallo da un hit all'altro, mediana {np.nanmedian(hit2hit):.2f} ADC", True, "informativo")
# 10 board jumps with the laser
bj = h3("hBoardJump"); nlas = int(h1("hLaserSlice").sum())
if nlas > 0:
    vals = np.arange(-20, 21); jumps = []
    for b in range(NB):
        tot = bj[b].sum(1)
        if tot.sum() == 0: continue
        med_sl = [vals[np.searchsorted(np.cumsum(bj[b, j]), bj[b, j].sum() / 2)] if bj[b, j].sum() > 0 else None for j in range(NSL)]
        allm = vals[np.searchsorted(np.cumsum(bj[b].sum(0)), bj[b].sum() / 2)]
        for j, m in enumerate(med_sl):
            if m is not None and abs(m - allm) >= CRIT["jump_samples"]: jumps.append((b, j, int(m - allm)))
    check("salti di scheda con il laser", f"{nlas} impulsi laser; {len(jumps)} (scheda, fetta) con il tempo mediano spostato di ≥ {CRIT['jump_samples']} campione", len(jumps) == 0, "nessuna scheda che salta di un campione intero dentro il run",
          detail=", ".join(f"scheda {b} a {tsl[j]:.0f} min: {s:+d}" for b, j, s in jumps[:20]))
else:
    check("salti di scheda con il laser", "nessun impulso laser nel run", True, "non applicabile: senza laser i salti si controllano dopo la ricostruzione (modulo 2)")

# ---------------- figures ----------------
figs = []
def save(fig, name, cap): fig.tight_layout(); fig.savefig(IMG / name, dpi=110); plt.close(fig); figs.append((name, cap))
# ---- rates and occupancy (first: the most important) ----
dur_s = max(dur * 60.0, 1e-9)
slen = np.full(NSL, SLM * 60.0)
if nsl_used > 0: slen[nsl_used - 1] = max(dur_s - (nsl_used - 1) * SLM * 60.0, 1.0)   # the last slice is partial
rate_ch = occ / dur_s                                            # Hz per channel, good hits of non-laser events
osl = h2("hOccSlice") if D.Get("hOccSlice") else None
sensor = np.where(keyoff >= 0, keyoff % 2, -1); crys = np.array([cry_of(k) for k in range(NK)]); dsk_k = np.where(crys >= 0, (crys >= 674).astype(int), -1)
fig, ax = plt.subplots(1, 3, figsize=(18, 4.3))
if D.Get("hRateFine"):
    rf = h2("hRateFine"); nfine = int(np.ceil(dur_s / 10.0)) if dur_s > 0 else 1; nfine = max(1, min(nfine, rf.shape[0]))
    grp = max(1, int(np.ceil(nfine / 50.0))); nb_ = int(np.ceil(nfine / grp)); edges = np.arange(nb_ + 1) * grp * 10.0
    def rebin(col): v = rf[:nb_ * grp, col]; v = np.pad(v, (0, nb_ * grp - len(v))); return v.reshape(nb_, grp).sum(1)
    blen = np.minimum(edges[1:], dur_s) - edges[:-1]; blen = np.maximum(blen, 1e-3); tc_ = 0.5 * (edges[1:] + np.minimum(edges[1:], dur_s)) / 60.0
    tc_ = (edges[:-1] + 0.5 * blen) / 60.0
    def plot_rate(a, col, label, color, scale=1.0, unit=""):
        n_ = rebin(col); r_ = n_ / blen * scale; e_ = np.sqrt(n_) / blen * scale
        a.errorbar(tc_, r_, e_, fmt="o-", ms=3, lw=1, color=color, label=label); return r_
    r_ev = plot_rate(ax[0], 0, "all events", "#1f62b8")
    if rf[:, 1].sum() > 0: plot_rate(ax[0], 1, "laser events", "#9a6a00")
    ax[0].set_ylim(0, 1.25 * max(np.max(r_ev), 1e-6)); ax[0].set_ylabel("events / s"); ax[0].set_title(f"event rate ({grp*10:.0f} s bins)", fontsize=10)
    r0 = plot_rate(ax[1], 2, "disk 0", "#1f62b8", 1e-3); r1 = plot_rate(ax[1], 3, "disk 1", "#c9531f", 1e-3)
    ax[1].set_ylim(0, 1.25 * max(np.max(r0), np.max(r1), 1e-6)); ax[1].set_ylabel("good hits / s [kHz]"); ax[1].set_title("hit rate of the non-laser events, per disk", fontsize=10)
    rl = plot_rate(ax[2], 4, "LYSO (good hits, non-laser)", "#1d7a4c"); rb_ = plot_rate(ax[2], 5, "dropped hits (all events)", "#c9531f")
    ax[2].set_ylim(0, 1.25 * max(np.max(rl), np.max(rb_), 1e-6)); ax[2].set_ylabel("hits / s"); ax[2].set_title("LYSO hits and hits dropped by the decoder", fontsize=10)
    for a in ax: a.set_xlabel("time from the first event [min]"); a.grid(alpha=.3); a.legend(fontsize=8)
save(fig, "o1_rates.png", "Tassi nel tempo, con l'errore statistico e la scala da 0: eventi al secondo (e impulsi laser), hit buoni per disco negli eventi senza laser, hit dei LYSO e hit scartati dal decoder. Gli intervalli si adattano alla durata del run (~50 punti).")
fig, ax = plt.subplots(2, 2, figsize=(14, 11.5))
pos = rate_ch[mapped & (rate_ch > 0)]; vmin, vmax = (np.percentile(pos, 2), np.percentile(pos, 98)) if len(pos) else (1e-3, 1)
from matplotlib.colors import LogNorm
for i, d_ in enumerate((0, 1)):
    for j, sn in enumerate((0, 1)):
        a = ax[i, j]; m_ = [k for k in range(NK) if crys[k] >= 0 and crys[k] < len(cx) and dsk_k[k] == d_ and sensor[k] == sn]
        xs, ys, vs = cx[crys[m_]], cy[crys[m_]], np.maximum(rate_ch[m_], vmin * 0.5)
        sc = a.scatter(xs, ys, c=vs, s=22, marker="s", cmap="viridis", norm=LogNorm(vmin=vmin, vmax=vmax))
        dead_ = [k for k in m_ if occ[k] < CRIT["dead_min_hits"]]; a.scatter(cx[crys[dead_]], cy[crys[dead_]], s=40, marker="x", color="red", label=f"no hits ({len(dead_)})")
        noisy_ = [k for k in m_ if k in noisy]; a.scatter(cx[crys[noisy_]], cy[crys[noisy_]], s=60, facecolors="none", edgecolors="magenta", label=f"> {CRIT['noisy_factor']:.0f}× median ({len(noisy_)})")
        a.set_aspect("equal"); a.set_title(f"disk {d_}, SiPM {sn}: hit rate per channel", fontsize=10); a.set_xlabel("x [mm]"); a.set_ylabel("y [mm]"); a.legend(fontsize=8, loc="lower right")
        plt.colorbar(sc, ax=a, fraction=0.046, label="good hits / s (non-laser events)")
save(fig, "o2_occupancy_maps.png", "Occupazione: tasso di hit buoni per canale [Hz], negli eventi senza laser, sui due dischi e per i due SiPM di ogni cristallo (posizioni dall'Offline); croci rosse: canali senza hit; cerchi: canali rumorosi.")
fig, ax = plt.subplots(1, 2, figsize=(15, 4.3))
lr = np.log10(np.maximum(rate_ch[mapped], 1e-4)); ax[0].hist(lr, 100, color="#1f62b8")
ax[0].axvline(np.log10(max(medocc / dur_s, 1e-4)), color="k", ls="--", label=f"median {medocc/dur_s:.3g} Hz"); ax[0].axvline(np.log10(max(CRIT['noisy_factor'] * medocc / dur_s, 1e-4)), color="magenta", ls=":", label=f"noisy > {CRIT['noisy_factor']:.0f}× median")
ax[0].set_xlabel("log10(hit rate per channel / Hz)"); ax[0].set_ylabel("channels"); ax[0].set_title("distribution of the channel rates (channels with 0 hits at the left edge)", fontsize=10); ax[0].legend(fontsize=8); ax[0].grid(alpha=.3)
rb = rate_ch.reshape(NB, NC).sum(1); nb_ch = mapped.reshape(NB, NC).sum(1); ax[1].bar(np.arange(NB)[nb_ch > 0], rb[nb_ch > 0], color="#1d7a4c")
ax[1].set_xlabel("board"); ax[1].set_ylabel("hit rate [Hz] (sum of its channels)"); ax[1].set_title("hit rate per board", fontsize=10); ax[1].grid(alpha=.3)
save(fig, "o3_rate_distribution.png", "Distribuzione dei tassi dei canali e tasso per scheda.")
def chan_row(k): c = crys[k]; return f"<tr><td>{k//NC}</td><td>{k%NC}</td><td>{keyoff[k]}</td><td>{c}</td><td>{dsk_k[k]}</td><td>{sensor[k]}</td><td>{'LYSO' if caphri[k] else 'CsI'}</td><td class=\"num\">{rate_ch[k]:.3g}</td><td class=\"num\">{(occ[k]/medocc if medocc else 0):.2f}</td></tr>"
chan_table = ("<table><thead><tr><th>scheda</th><th>canale</th><th>SiPM offline</th><th>cristallo</th><th>disco</th><th>SiPM</th><th>tipo</th><th>tasso [Hz]</th><th>/ mediana</th></tr></thead><tbody>"
              + "".join(chan_row(k) for k in sorted(set(dead) | set(noisy), key=lambda k: -rate_ch[k])) + "</tbody></table>")
rate_summary = (f"Durata {dur:.1f} min (dai tag dell'event window). Eventi: {nev/dur_s:.1f} /s ({int(la_.sum())} laser). Hit buoni negli eventi senza laser: "
                f"{rate_ch[mapped & ~caphri & (dsk_k == 0)].sum()/1e3:.2f} kHz sul disco 0, {rate_ch[mapped & ~caphri & (dsk_k == 1)].sum()/1e3:.2f} kHz sul disco 1; tasso mediano di un canale {medocc/dur_s:.3g} Hz.")
fig, ax = plt.subplots(1, 3, figsize=(16, 4))
ax[0].step(tsl[:nsl_used], rate, where="mid"); ax[0].set_xlabel(f"time [min] (slices of {SLM:.0f} min)"); ax[0].set_ylabel("events / s"); ax[0].set_title("event rate", fontsize=10)
la = la_; ax[0].step(tsl[:nsl_used], la[:nsl_used] / (SLM * 60), where="mid", label="laser events / s"); ax[0].legend(fontsize=8)
es = h2("hErrSlice"); [ax[1].step(tsl[:nsl_used], es[:nsl_used, c] / np.maximum(es[:nsl_used].sum(1), 1), where="mid", label=ERRNAME[c]) for c in range(1, 10) if es[:, c].sum() > 0]
ax[1].set_yscale("log"); ax[1].set_xlabel("time [min]"); ax[1].set_ylabel("fraction of the hits"); ax[1].set_title("hits dropped by the Offline decoder, by reason", fontsize=10)
if es[:, 1:].sum() > 0: ax[1].legend(fontsize=8)
fb = h2("hFlagBitSlice"); [ax[2].step(tsl[:nsl_used], fb[:nsl_used, b], where="mid", label=f"bit {b}") for b in range(4)]
ax[2].set_xlabel("time [min]"); ax[2].set_ylabel("hits with the bit set"); ax[2].set_title("ErrorFlags of the hits, bit by bit", fontsize=10); ax[2].legend(fontsize=8)
for a in ax: a.grid(alpha=.3)
save(fig, "r1_time.png", "Nel tempo: tasso degli eventi (e degli impulsi laser), hit scartati dal decoder per motivo, flag di errore degli hit bit per bit.")
fig, ax = plt.subplots(1, 3, figsize=(16, 4.4))
ax[0].bar(range(1, 10), err[1:], color="#c9531f"); ax[0].set_xticks(range(1, 10)); ax[0].set_xticklabels([ERRNAME[i] for i in range(1, 10)], rotation=45, ha="right", fontsize=8); ax[0].set_ylabel("hits"); ax[0].set_title("dropped hits by reason", fontsize=10)
bb = eperk[:, 1:].sum(1).reshape(NB, NC).sum(1); ax[1].bar(range(NB), bb, color="#c9531f"); ax[1].set_xlabel("board"); ax[1].set_ylabel("dropped hits"); ax[1].set_title("dropped hits per board", fontsize=10)
ax[2].bar(range(16), fl, color="#1f62b8"); ax[2].set_yscale("log"); ax[2].set_xlabel("ErrorFlags value"); ax[2].set_ylabel("hits"); ax[2].set_title("ErrorFlags values (all hits)", fontsize=10)
for a in ax: a.grid(alpha=.3)
save(fig, "r2_errors.png", "Hit scartati dal decoder dell'Offline per motivo e per scheda; valori del campo ErrorFlags.")
def disk_map(ax, val, title, cmap="viridis", vmin=None, vmax=None, lab=""):
    for dsk in (0, 1):
        a = ax[dsk]; xs, ys, vs = [], [], []
        for k in range(NK):
            c = cry_of(k)
            if c < 0 or c >= len(cx) or (c >= 674) != bool(dsk) or keyoff[k] % 2: continue     # one point per crystal (SiPM 0); offline ids beyond the crystals (e.g. PIN) skipped
            xs.append(cx[c]); ys.append(cy[c]); vs.append(val[k])
        sc = a.scatter(xs, ys, c=vs, s=14, marker="s", cmap=cmap, vmin=vmin, vmax=vmax); a.set_aspect("equal"); a.set_title(f"{title}, disk {dsk}", fontsize=10); a.set_xlabel("x [mm]"); a.set_ylabel("y [mm]")
        plt.colorbar(sc, ax=a, fraction=0.046, label=lab)
fig, ax = plt.subplots(1, 4, figsize=(20, 4.6)); rel = np.where(medocc > 0, occ / medocc, 0)
disk_map(ax[:2], rel, "hits / median (SiPM 0)", vmin=0, vmax=2, lab="occupancy / median"); disk_map(ax[2:], pf, "partner coincidence (SiPM 0)", vmin=0, vmax=1, lab="fraction")
save(fig, "r3_maps.png", "Mappe dei due dischi (posizioni dei cristalli dall'Offline): numero di hit rispetto alla mediana, e frazione di hit con un hit dell'altro SiPM del cristallo nello stesso evento.")
fig, ax = plt.subplots(1, 3, figsize=(16, 4.2))
ax[0].hist(pedall[ok], 80, color="#1f62b8"); ax[0].set_xlabel("pedestal [ADC]"); ax[0].set_ylabel("channels"); ax[0].set_title("pedestal per channel", fontsize=10)
ax[1].hist(hit2hit[ok], 80, color="#1d7a4c"); ax[1].set_xlabel("pedestal rms from hit to hit [ADC]"); ax[1].set_ylabel("channels"); ax[1].set_title("noise of the pedestal", fontsize=10)
bdev = np.full((NB, NSL), np.nan)
for b in range(NB):
    ks = [k for k in range(b * NC, (b + 1) * NC) if ok[k]]
    if ks: bdev[b] = np.nanmedian(pedsl[ks] - pedall[ks][:, None], axis=0)
im = ax[2].imshow(bdev[:, :nsl_used], aspect="auto", origin="lower", cmap="RdBu_r", vmin=-2, vmax=2, extent=[0, nsl_used * SLM, -0.5, NB - 0.5]); plt.colorbar(im, ax=ax[2], label="median pedestal shift [ADC]")
ax[2].set_xlabel("time [min]"); ax[2].set_ylabel("board"); ax[2].set_title("pedestal stability: median shift per board", fontsize=10)
for a in ax[:2]: a.grid(alpha=.3)
save(fig, "r4_pedestals.png", "Piedistallo di ogni canale, suo rumore da un hit all'altro, e stabilità nel tempo (spostamento mediano dei canali di ogni scheda rispetto alla media del run).")
if moved:
    fig, ax = plt.subplots(figsize=(9, 4))
    for k in sorted(moved, key=lambda k: -shift[k])[:8]: ax.plot(tsl[:nsl_used], pedsl[k, :nsl_used] - pedall[k], "o-", ms=3, label=f"board {k//NC} ch {k%NC}")
    ax.axhline(0, color="k", lw=.8); ax.set_xlabel("time [min]"); ax.set_ylabel("pedestal - run mean [ADC]"); ax.grid(alpha=.3); ax.legend(fontsize=8); ax.set_title("channels whose pedestal moves most", fontsize=10)
    save(fig, "r5_ped_moved.png", "I canali il cui piedistallo si sposta di più nel tempo.")
fig, ax = plt.subplots(1, 3, figsize=(16, 4))
ax[0].step(nsc, nsall, where="mid"); ax[0].set_yscale("log"); ax[0].set_xlabel("samples in the waveform"); ax[0].set_ylabel("hits"); ax[0].set_title("waveform length", fontsize=10)
pk = h1("hPeakIdx"); ax[1].step(range(len(pk)), pk, where="mid"); ax[1].set_yscale("log"); ax[1].set_xlabel("index of the largest sample"); ax[1].set_title("position of the maximum", fontsize=10)
hb = h2("hHitTagShift").sum(0); ax[2].bar(range(-4, 5), hb, color="#9a6a00"); ax[2].set_yscale("log"); ax[2].set_xlabel("hit tag - event tag (clipped at ±4)"); ax[2].set_title("event-window tag of the hits", fontsize=10)
for a in ax: a.grid(alpha=.3)
save(fig, "r6_waveforms.png", "Lunghezza delle forme d'onda, posizione del massimo, tag degli hit rispetto all'evento.")

S = dict(run=RUN, file=ROOTF, params=dict(windowMicroseconds=WIN, sliceMinutes=SLM), criteria=CRIT, checks=checks, n_events=nev)
(OUT / "summary.json").write_text(json.dumps(S, indent=1, default=float))
col = {"OK": "#1d7a4c", "ATTENZIONE": "#9a6a00", "NO": "#c9531f"}
rows = "".join(f'<tr><td>{c["check"]}</td><td style="color:{col[c["result"]]};font-weight:bold">{c["result"]}</td><td>{c["value"]}</td><td>{c["criterion"]}</td><td>{c["detail"]}</td></tr>' for c in checks)
css = "body{margin:0;background:#fbfbf9;color:#16181d;font:15px/1.5 Charter,Georgia,serif}main{max-width:80rem;margin:0 auto;padding:2rem 1rem}h1,h2{font-family:Helvetica,Arial,sans-serif}table{border-collapse:collapse;width:100%;font-size:13px}td,th{border-bottom:1px solid #ddd;padding:4px 6px;vertical-align:top;text-align:left}img{max-width:100%}figcaption{color:#555;font-size:13px;margin-bottom:1.5rem}@media (prefers-color-scheme:dark){body{background:#15171b;color:#e8e9ec}td,th{border-color:#333}}"
html = f"""<!doctype html><html lang="it"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1"><title>Controllo dati grezzi {RUN}</title><style>{css}</style></head><body><main>
<h1>Controllo dei dati grezzi, run {RUN}</h1>
<p>Modulo art <code>MC2RawCheck</code> (mc2reco, Offline v13_40_00) sul file raw, prima della ricostruzione; legge i frammenti DTC come il decoder dell'Offline e tiene anche gli hit che il decoder scarta. File: <code>{ROOTF}</code>. Asse dei tempi: {WIN:.0f} µs per event window (ipotesi), fette di {SLM:.0f} min. Criteri decisi prima di guardare il run (in testa allo script).</p>
<h2>Controlli</h2><table><thead><tr><th>controllo</th><th>esito</th><th>valore</th><th>criterio</th><th>dettaglio</th></tr></thead><tbody>{rows}</tbody></table>
<h2>Tassi e occupazione</h2><p>{rate_summary}</p>{''.join(f'<figure><img src="img/{n}" alt=""><figcaption>{c}</figcaption></figure>' for n, c in figs[:3])}
<h3>Canali senza hit o rumorosi</h3>{chan_table}
<h2>Altri plot</h2>{''.join(f'<figure><img src="img/{n}" alt=""><figcaption>{c}</figcaption></figure>' for n, c in figs[3:])}
</main></body></html>"""
(OUT / f"raw_check_{RUN}.html").write_text(html, encoding="utf-8"); print("written", OUT / f"raw_check_{RUN}.html")
for c in checks: print(f"{c['result']:10s} {c['check']}: {c['value']}")
