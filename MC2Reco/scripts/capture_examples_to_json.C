// capture_examples_to_json.C -- tree "examples" and the template of MC2MuCapture -> json, for capture_examples_page.py.
// Usage: root -l -b -q 'Offline/MC2Reco/scripts/capture_examples_to_json.C("pass.root","out.json")'
#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <cstdio>
#include <vector>
void capture_examples_to_json(const char* in, const char* out) {
  TFile f(in); TTree* t = (TTree*)f.Get("capture/examples"); TH1D* h = (TH1D*)f.Get("capture/htempl");
  if (!t || !h) { printf("ERROR: no examples or template in %s\n", in); return; }
  int cls, kind, cry, disk, run, sub, evt, nMi; float amp, maxS, tau, A2, Q, t0a, t0b, A1a, A1b, t1a, t1b, tref; std::vector<float> *y0 = nullptr, *y1 = nullptr;
  t->SetBranchAddress("cls", &cls); t->SetBranchAddress("kind", &kind); t->SetBranchAddress("cry", &cry); t->SetBranchAddress("disk", &disk); t->SetBranchAddress("run", &run);
  t->SetBranchAddress("subrun", &sub); t->SetBranchAddress("event", &evt); t->SetBranchAddress("nMichel", &nMi); t->SetBranchAddress("amp", &amp); t->SetBranchAddress("maxS", &maxS);
  t->SetBranchAddress("tau", &tau); t->SetBranchAddress("A2", &A2); t->SetBranchAddress("Q", &Q); t->SetBranchAddress("y0", &y0); t->SetBranchAddress("y1", &y1);
  t->SetBranchAddress("t0a", &t0a); t->SetBranchAddress("t0b", &t0b); t->SetBranchAddress("A1a", &A1a); t->SetBranchAddress("A1b", &A1b); t->SetBranchAddress("t1a", &t1a); t->SetBranchAddress("t1b", &t1b); t->SetBranchAddress("tref", &tref);
  FILE* o = fopen(out, "w"); fprintf(o, "{\"template\": {\"x0\": %g, \"bw\": %g, \"v\": [", h->GetXaxis()->GetXmin(), h->GetBinWidth(1));
  for (int i = 1; i <= h->GetNbinsX(); ++i) fprintf(o, "%s%.5f", i > 1 ? "," : "", h->GetBinContent(i)); fprintf(o, "]}, \"examples\": [\n");
  for (Long64_t i = 0; i < t->GetEntries(); ++i) { t->GetEntry(i);
    fprintf(o, "%s{\"cls\":%d,\"kind\":%d,\"cry\":%d,\"disk\":%d,\"run\":%d,\"subrun\":%d,\"event\":%d,\"nMichel\":%d,\"amp\":%.1f,\"maxS\":%.2f,\"tau\":%.1f,\"A2\":%.1f,\"Q\":%.4f,"
            "\"t0a\":%.1f,\"t0b\":%.1f,\"A1a\":%.1f,\"A1b\":%.1f,\"t1a\":%.2f,\"t1b\":%.2f,\"tref\":%.2f,\"y0\":[", i ? "," : "", cls, kind, cry, disk, run, sub, evt, nMi, amp, maxS, tau, A2, Q, t0a, t0b, A1a, A1b, t1a, t1b, tref);
    for (size_t q = 0; q < y0->size(); ++q) fprintf(o, "%s%.1f", q ? "," : "", (*y0)[q]); fprintf(o, "],\"y1\":[");
    for (size_t q = 0; q < y1->size(); ++q) fprintf(o, "%s%.1f", q ? "," : "", (*y1)[q]); fprintf(o, "]}\n"); }
  fprintf(o, "]}\n"); fclose(o); printf("written %s: %lld examples\n", out, t->GetEntries());
}
