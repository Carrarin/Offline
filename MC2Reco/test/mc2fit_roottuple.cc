//
// mc2fit_roottuple: standalone check of MC2TemplateFit (no art) on an MC2 roottuple.
//
// It fits the hits of the roottuple (tree caloAnalyzer/tree) with MC2TemplateFit and writes the sidet-like tree read by
// mc2-ana scripts/sidet_track_crystals.C, to compare with the analysis fit on the cosmic tracks.
//
// Usage: mc2fit_roottuple <roottuple> <template file> <output> <first> <last> [minAmp=55] [fitLo=-80] [fitHi=-31]
//
// Steps, for every entry between first and last:
//   1. entries with more hits or samples than the buffers are skipped (not written)
//   2. for every hit with at least 8 samples inside the waveform buffer:
//        pedestal = mean of the first 4 samples, Vmax = largest sample - pedestal; hits with Vmax < minAmp are skipped
//   3. template fit of the raw waveform (MC2TemplateFit, pedestal and saturation handled by the fit)
//
// Output tree "sidet", one entry per roottuple entry:
//   dtcEW                   DTC event window
//   nHits                   hits written
//   iBoa, iCha, EWch        board, channel, event window of the hit
//   Tval                    5 ns * timetot (start of the hit)
//   templTime               fitted peak time, ns from the first sample (0 if the fit failed)
//   Vmax                    largest sample - pedestal (ADC)
//   templChi2               chi2/ndf of the fit, at least 1e-6 (-1 if the fit failed)
//   templFit[3]             amplitude, time, residual baseline (0 if the fit failed)
//   nSamples                samples of the hit
// At the end it prints the number of fitted and failed hits.
//
#include "Offline/MC2Reco/inc/MC2TemplateFit.hh"
#include <TChain.h>
#include <TFile.h>
#include <TTree.h>
#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <vector>

namespace {

  constexpr int kMaxHits = 5000;          // hits per entry
  constexpr int kMaxSamples = 1500000;    // ADC samples per entry

  // Branches read from the roottuple.
  struct InputEntry {
    Long64_t eventWindow;
    int nHits;
    int board[kMaxHits];
    int channel[kMaxHits];
    Long64_t hitEventWindow[kMaxHits];
    int timeTot[kMaxHits];
    int nSamplesOfHit[kMaxHits];
    int firstSample[kMaxHits];
    int nSamples;
    int adc[kMaxSamples];
  };

  // Branches of the output tree.
  struct OutputEntry {
    Long64_t eventWindow;
    int nHits;
    int board[kMaxHits];
    int channel[kMaxHits];
    Long64_t hitEventWindow[kMaxHits];
    double startTime[kMaxHits];
    double fitTime[kMaxHits];
    double vMax[kMaxHits];
    double chi2[kMaxHits];
    double fit[kMaxHits][3];
    int nSamples[kMaxHits];
  };

  void connectInput(TChain& chain, InputEntry& in) {
    chain.SetBranchStatus("*", 0);
    for (auto name : {"currentDTCEventWindow", "nhits", "boardID", "chanID", "ewhit", "timetot", "nofsamples",
                      "firstsample", "nsamples", "ADC"}) {
      chain.SetBranchStatus(name, 1);
    }
    chain.SetBranchAddress("currentDTCEventWindow", &in.eventWindow);
    chain.SetBranchAddress("nhits", &in.nHits);
    chain.SetBranchAddress("boardID", in.board);
    chain.SetBranchAddress("chanID", in.channel);
    chain.SetBranchAddress("ewhit", in.hitEventWindow);
    chain.SetBranchAddress("timetot", in.timeTot);
    chain.SetBranchAddress("nofsamples", in.nSamplesOfHit);
    chain.SetBranchAddress("firstsample", in.firstSample);
    chain.SetBranchAddress("nsamples", &in.nSamples);
    chain.SetBranchAddress("ADC", in.adc);
  }

  void bookOutput(TTree& tree, OutputEntry& out) {
    tree.Branch("dtcEW", &out.eventWindow, "dtcEW/L");
    tree.Branch("nHits", &out.nHits, "nHits/I");
    tree.Branch("iBoa", out.board, "iBoa[nHits]/I");
    tree.Branch("iCha", out.channel, "iCha[nHits]/I");
    tree.Branch("EWch", out.hitEventWindow, "EWch[nHits]/L");
    tree.Branch("Tval", out.startTime, "Tval[nHits]/D");
    tree.Branch("templTime", out.fitTime, "templTime[nHits]/D");
    tree.Branch("Vmax", out.vMax, "Vmax[nHits]/D");
    tree.Branch("templChi2", out.chi2, "templChi2[nHits]/D");
    tree.Branch("templFit", out.fit, "templFit[nHits][3]/D");
    tree.Branch("nSamples", out.nSamples, "nSamples[nHits]/I");
  }

}

int main(int argc, char** argv) {
  if (argc < 6) {
    fprintf(stderr, "usage: %s roottuple template output first last [minAmp fitLo fitHi]\n", argv[0]);
    return 1;
  }

  // Fit configuration from the command line
  mu2e::MC2TemplateFit::Config config;
  config.templateFile = argv[2];
  const double minAmp = argc > 6 ? atof(argv[6]) : 55.;
  if (argc > 7) {
    config.fitLo = atof(argv[7]);
  }
  if (argc > 8) {
    config.fitHi = atof(argv[8]);
  }
  mu2e::MC2TemplateFit fitter(config);

  // Input roottuple and output tree (static: the buffers are large)
  static InputEntry in;
  static OutputEntry out;
  TChain chain("caloAnalyzer/tree");
  chain.Add(argv[1]);
  connectInput(chain, in);
  TFile outputFile(argv[3], "RECREATE");
  TTree tree("sidet", "MC2TemplateFit on MC2 roottuples");
  bookOutput(tree, out);

  long nFitted = 0;
  long nFailed = 0;
  std::vector<int> waveform;
  const Long64_t first = atoll(argv[4]);
  const Long64_t last = std::min<Long64_t>(atoll(argv[5]), chain.GetEntries() - 1);
  for (Long64_t entry = first; entry <= last; ++entry) {
    chain.GetEntry(entry);
    // 1. entries larger than the buffers
    if (in.nHits > kMaxHits || in.nSamples > kMaxSamples) {
      continue;
    }
    out.eventWindow = in.eventWindow;
    out.nHits = 0;

    for (int j = 0; j < in.nHits; ++j) {
      // 2. hit inside the buffer, pedestal and largest sample, amplitude cut
      const int nSamples = in.nSamplesOfHit[j];
      const int offset = in.firstSample[j];
      if (nSamples < 8 || offset < 0 || offset + nSamples > in.nSamples) {
        continue;
      }
      double pedestal = 0;
      for (int q = 0; q < 4; ++q) {
        pedestal += in.adc[offset + q];
      }
      pedestal /= 4.;
      int maxSample = in.adc[offset];
      for (int q = 1; q < nSamples; ++q) {
        maxSample = std::max(maxSample, in.adc[offset + q]);
      }
      if (maxSample - pedestal < minAmp) {
        continue;
      }

      // 3. template fit of the raw waveform
      waveform.assign(in.adc + offset, in.adc + offset + nSamples);
      const auto result = fitter.fit(waveform);
      ++nFitted;
      const int k = out.nHits++;
      out.board[k] = in.board[j];
      out.channel[k] = in.channel[j];
      out.hitEventWindow[k] = in.hitEventWindow[j];
      out.startTime[k] = 5. * in.timeTot[j];
      out.vMax[k] = maxSample - pedestal;
      out.nSamples[k] = nSamples;
      if (!result.ok) {
        ++nFailed;
        out.chi2[k] = -1;
        out.fitTime[k] = 0;
        out.fit[k][0] = out.fit[k][1] = out.fit[k][2] = 0;
        continue;
      }
      out.chi2[k] = std::max(result.chi2ndf, 1e-6);
      out.fitTime[k] = result.time;
      out.fit[k][0] = result.amplitude;
      out.fit[k][1] = result.time;
      out.fit[k][2] = result.baseline;
    }
    tree.Fill();
  }

  tree.Write();
  outputFile.Close();
  printf("entries %lld-%lld: hits fitted %ld, failed %ld\n", first, last, nFitted, nFailed);
  return 0;
}
