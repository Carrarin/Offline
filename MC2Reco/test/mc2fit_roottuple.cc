// Standalone check of MC2TemplateFit (no art): fit the hits of an MC2 roottuple and write the sidet-like tree read by
// mc2-ana scripts/sidet_track_crystals.C, to compare with the analysis fit on the cosmic tracks.
// Usage: mc2fit_roottuple <roottuple> <template file> <output> <first> <last> [minAmp=55] [fitLo=-80] [fitHi=-31]
#include "Offline/MC2Reco/inc/MC2TemplateFit.hh"
#include <TChain.h>
#include <TFile.h>
#include <TTree.h>
#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <vector>
int main(int argc, char** argv) {
  if (argc < 6) { fprintf(stderr, "usage: %s roottuple template output first last [minAmp fitLo fitHi]\n", argv[0]); return 1; }
  mu2e::MC2TemplateFit::Config cfg; cfg.templateFile = argv[2];
  const double minAmp = argc > 6 ? atof(argv[6]) : 55.; if (argc > 7) cfg.fitLo = atof(argv[7]); if (argc > 8) cfg.fitHi = atof(argv[8]);
  mu2e::MC2TemplateFit fitter(cfg);
  TChain c("caloAnalyzer/tree"); c.Add(argv[1]);
  const int M = 5000, MA = 1500000; static int nh, b[M], ch[M], tt[M], nofs[M], fs[M], nsamp, adc[MA]; static Long64_t eh[M]; Long64_t ew;
  c.SetBranchStatus("*", 0);
  for (auto s : {"currentDTCEventWindow", "nhits", "boardID", "chanID", "ewhit", "timetot", "nofsamples", "firstsample", "nsamples", "ADC"}) c.SetBranchStatus(s, 1);
  c.SetBranchAddress("currentDTCEventWindow", &ew); c.SetBranchAddress("nhits", &nh); c.SetBranchAddress("boardID", b); c.SetBranchAddress("chanID", ch);
  c.SetBranchAddress("ewhit", eh); c.SetBranchAddress("timetot", tt); c.SetBranchAddress("nofsamples", nofs); c.SetBranchAddress("firstsample", fs);
  c.SetBranchAddress("nsamples", &nsamp); c.SetBranchAddress("ADC", adc);
  TFile fo(argv[3], "RECREATE"); TTree t("sidet", "MC2TemplateFit on MC2 roottuples");
  static int oN, oB[M], oC[M], oNs[M]; static Long64_t oEW, oEh[M]; static double oTv[M], oTt[M], oVm[M], oChi[M], oFit[M][3];
  t.Branch("dtcEW", &oEW, "dtcEW/L"); t.Branch("nHits", &oN, "nHits/I"); t.Branch("iBoa", oB, "iBoa[nHits]/I"); t.Branch("iCha", oC, "iCha[nHits]/I");
  t.Branch("EWch", oEh, "EWch[nHits]/L"); t.Branch("Tval", oTv, "Tval[nHits]/D"); t.Branch("templTime", oTt, "templTime[nHits]/D");
  t.Branch("Vmax", oVm, "Vmax[nHits]/D"); t.Branch("templChi2", oChi, "templChi2[nHits]/D"); t.Branch("templFit", oFit, "templFit[nHits][3]/D");
  t.Branch("nSamples", oNs, "nSamples[nHits]/I");
  long nFit = 0, nFail = 0; std::vector<int> wf;
  const Long64_t first = atoll(argv[4]), last = std::min<Long64_t>(atoll(argv[5]), c.GetEntries() - 1);
  for (Long64_t i = first; i <= last; ++i) {
    c.GetEntry(i); if (nh > M || nsamp > MA) continue; oEW = ew; oN = 0;
    for (int j = 0; j < nh; ++j) {
      const int n = nofs[j], o = fs[j]; if (n < 8 || o < 0 || o + n > nsamp) continue;
      double ped = 0; for (int q = 0; q < 4; ++q) ped += adc[o + q]; ped /= 4.; int mx = adc[o]; for (int q = 1; q < n; ++q) mx = std::max(mx, adc[o + q]);
      if (mx - ped < minAmp) continue;
      wf.assign(adc + o, adc + o + n); const auto r = fitter.fit(wf); ++nFit;
      const int k = oN++; oB[k] = b[j]; oC[k] = ch[j]; oEh[k] = eh[j]; oTv[k] = 5. * tt[j]; oVm[k] = mx - ped; oNs[k] = n;
      if (!r.ok) { ++nFail; oChi[k] = -1; oTt[k] = 0; oFit[k][0] = oFit[k][1] = oFit[k][2] = 0; continue; }
      oChi[k] = std::max(r.chi2ndf, 1e-6); oTt[k] = r.time; oFit[k][0] = r.amplitude; oFit[k][1] = r.time; oFit[k][2] = r.baseline;
    }
    t.Fill();
  }
  t.Write(); fo.Close(); printf("entries %lld-%lld: hits fitted %ld, failed %ld\n", first, last, nFit, nFail); return 0;
}
