#!/usr/bin/env python3
"""Writes the fcl of the reconstruction comparison job (MC2RecoCompare) for one run.

One job, the same CaloDigis for every method:
  reference : MC2CaloRecoDigiMaker (cosmic v2, [-80,-31] ns) -> CaloHitMaker -> MC2CaloCosmicTrackFinder
              (the tracks of all methods)
  methods   : Offline CaloRecoDigiMaker as in the pass1, the same with refitLeadingEdge, and our fit
              (MC2CaloRecoDigiMaker) with every template of TEMPLATES and every window end of ENDS (start -80 ns);
              one family per template.
Time and energy tables of the run as Offline conditions (DbService.textFile, calCalib.useDb): every method subtracts
the same per-SiPM offsets; MC2RecoCompare then calibrates the residual offsets of each method on the even events.

Inputs:
  <time table>    CalTimeCalib text table of the run (DbService format)
  <energy table>  CalEnergyCalib text table of the run (DbService format)
  <out dir>       directory of the job output: the fcl writes <out dir>/compare_summary.txt there
  All three paths are written into the fcl relative to the current directory (the muse area).

Output: the fcl, on stdout.

Reference options: --ref-merge NS (CaloHitMaker time4Merge of the reference, default 6 ns as Offline; 50 ns keeps
crystals whose L and R the reference fit disagrees on) and --ref-raw-amp (track finder cut on the raw waveform
amplitude instead of the fit one).
--fix FAMILY=END (repeatable): that family is run only with the window end END (e.g. chosen by an earlier job on the
even events) instead of the whole scan ENDS.
--only FAM1,FAM2,...: run only these template families (the two Offline methods are always there).
--walk: walk correction w(A) of the SiPM raw amplitude, fitted with the offsets on the even events and subtracted on
both parities.
--extra FAMILY:START:END (repeatable): one more method of that family with the window [START, END] (e.g. a later
start, to test what the samples before the edge give); named FAMILY_START_END.
Families with perChannel = True use one template per SiPM (histograms tmpl_bBBB_cCC of the file), the common one
where missing.

Usage (from the muse area):
  Offline/MC2Reco/scripts/make_reco_compare_fcl.py [--ref-merge 50] [--ref-raw-amp] [--fix v3_r125367=-30 ...]
      <time table> <energy table> <out dir> > <out dir>/compare.fcl
Example:
  Offline/MC2Reco/scripts/make_reco_compare_fcl.py --fix v3_r125367=-30 out/125367/CalTimeCalib_125367.txt \\
      out/125367/CalEnergyCalib_125367.txt out/125367/compare > out/125367/compare/compare.fcl
"""
import argparse
import os

# Template families: (family, template file, one template per SiPM [, shapeMode [, shapeSigma]]).
# Files in Offline/MC2Reco/data/templates.
TEMPLATES = [
    ("v3_r125367", "Offline/MC2Reco/data/templates/splines_cosmic_r125367_v3.root", False),
    ("v2_r125328", "Offline/MC2Reco/data/templates/splines_cosmic_r125328_v2.root", False),
    ("offline2020", "Offline/MC2Reco/data/templates/splines_offline_CsI2020.root", False),
    ("run39", "Offline/MC2Reco/data/templates/splines4_run39.root", False),
    ("stefano_r125163", "Offline/MC2Reco/data/templates/splines_stefano_r125163.root", False),
    ("stefano_r125163_perch", "Offline/MC2Reco/data/templates/splines_stefano_r125163.root", True),
    ("v3_r125367_perch", "Offline/MC2Reco/data/templates/templates_cosmic_r125367_v3_perchannel.root", True),
    # v3 recipe, hits 1200-3000 ADC only
    ("v3high_r125367", "Offline/MC2Reco/data/templates/splines_cosmic_r125367_v3high.root", False),
    # amplitude model T0 + u T1 + u^2 T2, u = ln(A / 700 ADC) (mc2-ana scripts/template_perf/template_amp_model.py):
    # shapeMode 1 = shape tied to the fitted amplitude, 2 = free shape coefficient of T1
    ("ampmodel_r125367", "Offline/MC2Reco/data/templates/ampmodel_cosmic_r125367.root", False, 1),
    ("ampfree_r125367", "Offline/MC2Reco/data/templates/ampmodel_cosmic_r125367.root", False, 2),
    # shapeMode 3 = shape tied to A plus a deviation D T1 of each hit, penalised with width shapeSigma (units of ln A)
    ("ampcons02_r125367", "Offline/MC2Reco/data/templates/ampmodel_cosmic_r125367.root", False, 3, 0.2),
    ("ampcons05_r125367", "Offline/MC2Reco/data/templates/ampmodel_cosmic_r125367.root", False, 3, 0.5),
    ("ampcons10_r125367", "Offline/MC2Reco/data/templates/ampmodel_cosmic_r125367.root", False, 3, 1.0),
]
# Window ends of the scan [ns] and the common window start [ns] of our fit
ENDS = [-40, -36, -34, -32, -31, -30, -28, -26, -24, -20]
FIT_LO = -80.

# Producer line of one of our fit methods (MC2CaloRecoDigiMaker with a template file and a fit window)
MC2_PRODUCER_FORMAT = ('      %-21s : { @table::MC2Reco.MC2CaloRecoDigiMaker templateFile : "%s" perChannel : %s '
                       'fitLo : %.1f fitHi : %.1f noise : 0.7 useCalCalibTimeOffset : true%s }')


def shape_opts(shape_params):
    """Extra fcl of an amplitude-model family: shapeMode [, shapeSigma] (empty for a plain template family)."""
    mode_text = " shapeMode : %d" % shape_params[0] if shape_params else ""
    sigma_text = " shapeSigma : %.2f" % shape_params[1] if len(shape_params) > 1 else ""
    return mode_text + sigma_text


def parse_args():
    """Command line."""
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--ref-merge", type=float, default=6., help="time4Merge of the reference CaloHitMaker [ns]")
    parser.add_argument("--ref-raw-amp", action="store_true", help="track finder cut on the raw amplitude")
    parser.add_argument("--fix", action="append", default=[], metavar="FAMILY=END",
                        help="run FAMILY only at window end END [ns]")
    parser.add_argument("--only", default="", help="comma-separated template families to run (default: all)")
    parser.add_argument("--extra", action="append", default=[], metavar="FAMILY:START:END",
                        help="one more method with window [START, END] [ns]")
    parser.add_argument("--walk", action="store_true",
                        help="MC2RecoCompare fits a walk w(A) with the offsets on the even events and subtracts it")
    parser.add_argument("ttab")
    parser.add_argument("etab")
    parser.add_argument("out")
    return parser, parser.parse_args()


def fixed_window_ends(parser, fix_options):
    """--fix FAMILY=END options -> {family: [END]}; unknown families are a command-line error."""
    families = [t[0] for t in TEMPLATES]
    fixed = {}
    for option in fix_options:
        family, end = option.split("=")
        if family not in families:
            parser.error("unknown family " + family)
        fixed[family] = [int(end)]
    return fixed


def extra_windows(parser, extra_options):
    """--extra FAMILY:START:END options -> [(family, start, end)]; unknown families are a command-line error."""
    families = [t[0] for t in TEMPLATES]
    extra = []
    for option in extra_options:
        family, start, end = option.split(":")
        if family not in families:
            parser.error("unknown family " + family)
        extra.append((family, int(start), int(end)))
    return extra


def reference_and_offline_producers(args, producers, path, methods):
    """The digis, the reference chain (reco digis, hits, tracks) and the two Offline methods."""
    producers.append('      CaloDigis             : { @table::DAQ.producers.CaloDigis }')
    producers.append('      RefRecoDigi           : { @table::MC2Reco.MC2CaloRecoDigiMaker '
                     'useCalCalibTimeOffset : true }')
    producers.append('      RefCaloHit            : { @table::CaloReco.CaloHitMaker '
                     'caloRecoDigiCollection : RefRecoDigi time4Merge : %.1f }' % args.ref_merge)
    producers.append('      RefTracks             : { @table::MC2Reco.MC2CaloCosmicTrackFinder caloHitCollection : '
                     'RefCaloHit rawAmplitude : %s }'
                     % ("true" if args.ref_raw_amp else "false"))
    producers.append('      OfflineRecoDigi       : { @table::CaloReco.CaloRecoDigiMaker '
                     'caloDigiCollection : CaloDigis UseProtonBunchTime : false }')
    producers.append('      OfflineEdgeRecoDigi   : { @table::MC2Reco.CaloRecoDigiMakerEdge '
                     'UseProtonBunchTime : false }')
    path += ["CaloDigis", "RefRecoDigi", "RefCaloHit", "RefTracks", "OfflineRecoDigi", "OfflineEdgeRecoDigi"]
    methods.append(("offline_pass1", "offline", "OfflineRecoDigi"))
    methods.append(("offline_edge", "offline_edge", "OfflineEdgeRecoDigi"))


def template_scan_producers(only, fixed, producers, path, methods):
    """Our fit with every selected template family and every window end (or the fixed end of the family)."""
    for family, template_file, per_channel, *shape_params in TEMPLATES:
        if only and family not in only:
            continue
        shape = shape_opts(shape_params)
        for end in fixed.get(family, ENDS):
            label = "MC2%s%d" % (family.replace("_", ""), -end)
            producers.append(MC2_PRODUCER_FORMAT
                             % (label, template_file, "true" if per_channel else "false", FIT_LO, end, shape))
            path.append(label)
            methods.append(("%s_%d" % (family, end), family, label))


def extra_window_producers(extra, producers, path, methods):
    """One more method for each --extra FAMILY:START:END, with the window [START, END]."""
    for family, start, end in extra:
        template_file, per_channel, *shape_params = [t[1:] for t in TEMPLATES if t[0] == family][0]
        label = "MC2%sS%dE%d" % (family.replace("_", ""), -start, -end)
        producers.append(MC2_PRODUCER_FORMAT
                         % (label, template_file, "true" if per_channel else "false", start, end,
                            shape_opts(shape_params)))
        path.append(label)
        methods.append(("%s_start%d_%d" % (family, start, end), family + "_start", label))


def fcl_text(args, time_table, energy_table, out_dir, fixed, producers, path, methods):
    """The whole fcl of the job."""
    method_list = ",\n".join('            { name : "%s" family : "%s" caloRecoDigiCollection : "%s" }' % m
                             for m in methods)
    return f"""# Reconstruction comparison job, written by Offline/MC2Reco/scripts/make_reco_compare_fcl.py \
(do not edit: change the script)
# Tables: {time_table}, {energy_table}; fixed window ends: {fixed if fixed else "none"}; \
reference: time4Merge {args.ref_merge} ns, raw-amplitude cut {args.ref_raw_amp}
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
{chr(10).join(producers)}
   }}
   analyzers : {{
      compare : {{
         module_type         : MC2RecoCompare
         caloTrackCollection : RefTracks
         summaryFile         : "{out_dir}/compare_summary.txt"
         writeTree           : true
         walkCorrection      : {"true" if args.walk else "false"}
         methods : [
{method_list}
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
services.DbService.textFile : [ "{energy_table}", "{time_table}" ]
services.scheduler.wantSummary : true"""


def main():
    parser, args = parse_args()
    # DbService finds its tables in the search path: paths relative to the current directory (the muse area)
    time_table, energy_table, out_dir = (os.path.relpath(os.path.abspath(p)) for p in (args.ttab, args.etab, args.out))
    fixed = fixed_window_ends(parser, args.fix)
    extra = extra_windows(parser, args.extra)

    producers, path, methods = [], [], []
    reference_and_offline_producers(args, producers, path, methods)
    only = [family for family in args.only.split(",") if family]
    template_scan_producers(only, fixed, producers, path, methods)
    extra_window_producers(extra, producers, path, methods)

    print(fcl_text(args, time_table, energy_table, out_dir, fixed, producers, path, methods))


if __name__ == "__main__":
    main()
