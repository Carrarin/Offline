// capture_examples_to_json.C: tree "examples" and the pulse template of MC2MuCapture -> json, for
// capture_examples_page.py (no PyROOT needed).
//
// Usage (muse environment, or any ROOT):
//
//   root -l -b -q 'Offline/MC2Reco/scripts/capture_examples_to_json.C("pass.root","out.json")'
//
//   e.g. root -l -b -q 'Offline/MC2Reco/scripts/capture_examples_to_json.C("out/physics.root","out/examples.json")'
//
// Arguments: in = TFileService output of a job with the MC2MuCapture analyzer "capture" (tree capture/examples and
//            histogram capture/htempl); out = json file to write.
//
// Output: <out>, one json object
//   {"template": {"x0": <low edge>, "bw": <bin width>, "v": [bin contents]},
//    "examples": [ one object per entry, one per line, with the scalar branches
//                  cls, kind, cry, disk, run, subrun, event, nMichel, amp, maxS, tau, A2, Q, t0a, t0b, A1a, A1b, t1a,
//                  t1b, tref and the waveforms y0, y1 (arrays) ]}
// Prints the number of examples.

#include <TFile.h>
#include <TH1D.h>
#include <TTree.h>

#include <cstdio>
#include <vector>

// The template histogram: "x0" (low edge of the axis), "bw" (bin width) and the bin contents "v"
void writeTemplate(FILE* json, TH1D* templateHist) {
  fprintf(json, "{\"template\": {\"x0\": %g, \"bw\": %g, \"v\": [", templateHist->GetXaxis()->GetXmin(),
          templateHist->GetBinWidth(1));
  for (int i = 1; i <= templateHist->GetNbinsX(); ++i) {
    fprintf(json, "%s%.5f", i > 1 ? "," : "", templateHist->GetBinContent(i));
  }
  fprintf(json, "]}, \"examples\": [\n");
}

// A waveform as a json array of numbers with one decimal
void writeWaveform(FILE* json, const std::vector<float>& samples) {
  for (size_t i = 0; i < samples.size(); ++i) {
    fprintf(json, "%s%.1f", i ? "," : "", samples[i]);
  }
}

void capture_examples_to_json(const char* in, const char* out) {
  TFile inputFile(in);
  TTree* tree = (TTree*)inputFile.Get("capture/examples");
  TH1D* templateHist = (TH1D*)inputFile.Get("capture/htempl");
  if (!tree || !templateHist) {
    printf("ERROR: no examples or template in %s\n", in);
    return;
  }

  // branches of the examples tree
  int cls, kind, cry, disk, run, subrun, event, nMichel;
  float amp, maxS, tau, A2, Q, t0a, t0b, A1a, A1b, t1a, t1b, tref;
  std::vector<float>* y0 = nullptr;
  std::vector<float>* y1 = nullptr;
  tree->SetBranchAddress("cls", &cls);
  tree->SetBranchAddress("kind", &kind);
  tree->SetBranchAddress("cry", &cry);
  tree->SetBranchAddress("disk", &disk);
  tree->SetBranchAddress("run", &run);
  tree->SetBranchAddress("subrun", &subrun);
  tree->SetBranchAddress("event", &event);
  tree->SetBranchAddress("nMichel", &nMichel);
  tree->SetBranchAddress("amp", &amp);
  tree->SetBranchAddress("maxS", &maxS);
  tree->SetBranchAddress("tau", &tau);
  tree->SetBranchAddress("A2", &A2);
  tree->SetBranchAddress("Q", &Q);
  tree->SetBranchAddress("y0", &y0);
  tree->SetBranchAddress("y1", &y1);
  tree->SetBranchAddress("t0a", &t0a);
  tree->SetBranchAddress("t0b", &t0b);
  tree->SetBranchAddress("A1a", &A1a);
  tree->SetBranchAddress("A1b", &A1b);
  tree->SetBranchAddress("t1a", &t1a);
  tree->SetBranchAddress("t1b", &t1b);
  tree->SetBranchAddress("tref", &tref);

  FILE* json = fopen(out, "w");
  writeTemplate(json, templateHist);

  // one example per line, separated by commas
  for (Long64_t i = 0; i < tree->GetEntries(); ++i) {
    tree->GetEntry(i);
    fprintf(json,
            "%s{\"cls\":%d,\"kind\":%d,\"cry\":%d,\"disk\":%d,\"run\":%d,\"subrun\":%d,\"event\":%d,\"nMichel\":%d,"
            "\"amp\":%.1f,\"maxS\":%.2f,\"tau\":%.1f,\"A2\":%.1f,\"Q\":%.4f,"
            "\"t0a\":%.1f,\"t0b\":%.1f,\"A1a\":%.1f,\"A1b\":%.1f,\"t1a\":%.2f,\"t1b\":%.2f,\"tref\":%.2f,\"y0\":[",
            i ? "," : "", cls, kind, cry, disk, run, subrun, event, nMichel, amp, maxS, tau, A2, Q, t0a, t0b, A1a, A1b,
            t1a, t1b, tref);
    writeWaveform(json, *y0);
    fprintf(json, "],\"y1\":[");
    writeWaveform(json, *y1);
    fprintf(json, "]}\n");
  }

  fprintf(json, "]}\n");
  fclose(json);
  printf("written %s: %lld examples\n", out, tree->GetEntries());
}
