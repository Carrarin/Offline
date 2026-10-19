#!/usr/bin/env python3
"""Writes the fcl of the reconstruction comparison job (MC2RecoCompare) for one run.

One job, the same CaloDigis for every method:
  reference : MC2CaloRecoDigiMaker (cosmic v2, [-80,-31] ns) -> CaloHitMaker -> MC2CaloCosmicTrackFinder (the tracks of all methods)
  methods   : Offline CaloRecoDigiMaker as in the pass1, the same with refitLeadingEdge, and our fit (MC2CaloRecoDigiMaker) with every
              template of TEMPLATES and every window end of ENDS (start -80 ns); one family per template.
Time and energy tables of the run as Offline conditions (DbService.textFile, calCalib.useDb): every method subtracts the same
per-SiPM offsets; MC2RecoCompare then calibrates the residual offsets of each method on the even events.

Reference options: --ref-merge NS (CaloHitMaker time4Merge of the reference, default 6 ns as Offline; 50 ns keeps crystals whose
L and R the reference fit disagrees on) and --ref-raw-amp (track finder cut on the raw waveform amplitude instead of the fit one).
--fix FAMILY=END (repeatable): that family is run only with the window end END (e.g. chosen by an earlier job on the even events)
instead of the whole scan ENDS.
--only FAM1,FAM2,...: run only these template families (the two Offline methods are always there).
--walk: walk correction w(A) of the SiPM raw amplitude, fitted with the offsets on the even events and subtracted on both parities.
--extra FAMILY:START:END (repeatable): one more method of that family with the window [START, END] (e.g. a later start, to test
what the samples before the edge give); named FAMILY_START_END.
Families with perChannel = True use one template per SiPM (histograms tmpl_bBBB_cCC of the file), the common one where missing.

Usage (from the muse area):
  Offline/MC2Reco/scripts/make_reco_compare_fcl.py [--ref-merge 50] [--ref-raw-amp] [--fix v3_r125367=-30 ...] <time table> <energy table> <out dir> > <out dir>/compare.fcl
"""
import argparse
import os

TEMPLATES = [   # family, template file (Offline/MC2Reco/data/templates), one template per SiPM
    ("v3_r125367", "Offline/MC2Reco/data/templates/splines_cosmic_r125367_v3.root", False),
    ("v2_r125328", "Offline/MC2Reco/data/templates/splines_cosmic_r125328_v2.root", False),
    ("offline2020", "Offline/MC2Reco/data/templates/splines_offline_CsI2020.root", False),
    ("run39", "Offline/MC2Reco/data/templates/splines4_run39.root", False),
    ("stefano_r125163", "Offline/MC2Reco/data/templates/splines_stefano_r125163.root", False),
    ("stefano_r125163_perch", "Offline/MC2Reco/data/templates/splines_stefano_r125163.root", True),
    ("v3_r125367_perch", "Offline/MC2Reco/data/templates/templates_cosmic_r125367_v3_perchannel.root", True),
    ("v3high_r125367", "Offline/MC2Reco/data/templates/splines_cosmic_r125367_v3high.root", False),   # v3 recipe, hits 1200-3000 ADC only
    # amplitude model T0 + u T1 + u^2 T2, u = ln(A / 700 ADC) (mc2-ana scripts/template_perf/template_amp_model.py): shapeMode 1 = shape
    # tied to the fitted amplitude, 2 = free shape coefficient of T1
    ("ampmodel_r125367", "Offline/MC2Reco/data/templates/ampmodel_cosmic_r125367.root", False, 1),
    ("ampfree_r125367", "Offline/MC2Reco/data/templates/ampmodel_cosmic_r125367.root", False, 2),
    # shapeMode 3 = shape tied to A plus a deviation D T1 of each hit, penalised with width shapeSigma (units of ln A)
    ("ampcons02_r125367", "Offline/MC2Reco/data/templates/ampmodel_cosmic_r125367.root", False, 3, 0.2),
    ("ampcons05_r125367", "Offline/MC2Reco/data/templates/ampmodel_cosmic_r125367.root", False, 3, 0.5),
    ("ampcons10_r125367", "Offline/MC2Reco/data/templates/ampmodel_cosmic_r125367.root", False, 3, 1.0),
]
ENDS = [-40, -36, -34, -32, -31, -30, -28, -26, -24, -20]
FIT_LO = -80.


def shape_opts(sm):   # extra fcl of an amplitude-model family: shapeMode [, shapeSigma]
    return (" shapeMode : %d" % sm[0] if sm else "") + (" shapeSigma : %.2f" % sm[1] if len(sm) > 1 else "")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--ref-merge", type=float, default=6., help="time4Merge of the reference CaloHitMaker [ns]")
    ap.add_argument("--ref-raw-amp", action="store_true", help="track finder cut on the raw amplitude")
    ap.add_argument("--fix", action="append", default=[], metavar="FAMILY=END", help="run FAMILY only at window end END [ns]")
    ap.add_argument("--only", default="", help="comma-separated template families to run (default: all)")
    ap.add_argument("--extra", action="append", default=[], metavar="FAMILY:START:END", help="one more method with window [START, END] [ns]")
    ap.add_argument("--walk", action="store_true", help="MC2RecoCompare fits a walk w(A) with the offsets on the even events and subtracts it")
    ap.add_argument("ttab"); ap.add_argument("etab"); ap.add_argument("out")
    a = ap.parse_args()
    ttab, etab, out = (os.path.relpath(os.path.abspath(p)) for p in (a.ttab, a.etab, a.out))
    fixed = {}
    for f in a.fix:
        fam, end = f.split("=")
        if fam not in [t[0] for t in TEMPLATES]:
            ap.error("unknown family " + fam)
        fixed[fam] = [int(end)]
    extra = []
    for x in a.extra:
        fam, lo, hi = x.split(":")
        if fam not in [t[0] for t in TEMPLATES]:
            ap.error("unknown family " + fam)
        extra.append((fam, int(lo), int(hi)))
    prod, path, methods = [], [], []
    prod.append('      CaloDigis             : { @table::DAQ.producers.CaloDigis }')
    prod.append('      RefRecoDigi           : { @table::MC2Reco.MC2CaloRecoDigiMaker useCalCalibTimeOffset : true }')
    prod.append('      RefCaloHit            : { @table::CaloReco.CaloHitMaker caloRecoDigiCollection : RefRecoDigi time4Merge : %.1f }' % a.ref_merge)
    prod.append('      RefTracks             : { @table::MC2Reco.MC2CaloCosmicTrackFinder caloHitCollection : RefCaloHit rawAmplitude : %s }'
                % ("true" if a.ref_raw_amp else "false"))
    prod.append('      OfflineRecoDigi       : { @table::CaloReco.CaloRecoDigiMaker caloDigiCollection : CaloDigis UseProtonBunchTime : false }')
    prod.append('      OfflineEdgeRecoDigi   : { @table::MC2Reco.CaloRecoDigiMakerEdge UseProtonBunchTime : false }')
    path += ["CaloDigis", "RefRecoDigi", "RefCaloHit", "RefTracks", "OfflineRecoDigi", "OfflineEdgeRecoDigi"]
    methods.append(("offline_pass1", "offline", "OfflineRecoDigi"))
    methods.append(("offline_edge", "offline_edge", "OfflineEdgeRecoDigi"))
    only = [f for f in a.only.split(",") if f]
    for fam, tf, perch, *sm in TEMPLATES:
        if only and fam not in only:
            continue
        shape = shape_opts(sm)
        for e in fixed.get(fam, ENDS):
            lab = "MC2%s%d" % (fam.replace("_", ""), -e)
            prod.append('      %-21s : { @table::MC2Reco.MC2CaloRecoDigiMaker templateFile : "%s" perChannel : %s fitLo : %.1f fitHi : %.1f noise : 0.7 useCalCalibTimeOffset : true%s }'
                        % (lab, tf, "true" if perch else "false", FIT_LO, e, shape))
            path.append(lab)
            methods.append(("%s_%d" % (fam, e), fam, lab))
    for fam, lo, hi in extra:
        tf, perch, *sm = [t[1:] for t in TEMPLATES if t[0] == fam][0]
        lab = "MC2%sS%dE%d" % (fam.replace("_", ""), -lo, -hi)
        prod.append('      %-21s : { @table::MC2Reco.MC2CaloRecoDigiMaker templateFile : "%s" perChannel : %s fitLo : %.1f fitHi : %.1f noise : 0.7 useCalCalibTimeOffset : true%s }'
                    % (lab, tf, "true" if perch else "false", lo, hi, shape_opts(sm)))
        path.append(lab)
        methods.append(("%s_start%d_%d" % (fam, lo, hi), fam + "_start", lab))
    mlist = ",\n".join('            { name : "%s" family : "%s" caloRecoDigiCollection : "%s" }' % m for m in methods)
    print(f"""# Reconstruction comparison job, written by Offline/MC2Reco/scripts/make_reco_compare_fcl.py (do not edit: change the script)
# Tables: {ttab}, {etab}; fixed window ends: {fixed if fixed else "none"}; reference: time4Merge {a.ref_merge} ns, raw-amplitude cut {a.ref_raw_amp}
#include "Offline/fcl/minimalMessageService.fcl"
#include "Offline/fcl/standardServices.fcl"
#include "Offline/DAQ/fcl/prolog.fcl"
#include "Offline/CaloReco/fcl/prolog.fcl"
#include "Offline/MC2Reco/fcl/prolog.fcl"
process_name : MC2RecoCompare
source   : {{ module_type : RootInput }}
services : @local::Services.Reco
physics : {{
   producers : {{
{chr(10).join(prod)}
   }}
   analyzers : {{
      compare : {{
         module_type         : MC2RecoCompare
         caloTrackCollection : RefTracks
         summaryFile         : "{out}/compare_summary.txt"
         writeTree           : true
         walkCorrection      : {"true" if a.walk else "false"}
         methods : [
{mlist}
         ]
      }}
   }}
   reco : [ {", ".join(path)} ]
   ana  : [ compare ]
   trigger_paths : [ reco ]
   end_paths     : [ ana ]
}}
services.GeometryService.inputFile : "Offline/Mu2eG4/geom/geom_common_extracted.txt"
services.ProditionsService.calCalib.useDb : true
services.DbService.textFile : [ "{etab}", "{ttab}" ]
services.scheduler.wantSummary : true""")


if __name__ == "__main__":
    main()
