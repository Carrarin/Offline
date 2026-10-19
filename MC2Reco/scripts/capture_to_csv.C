// capture_to_csv.C -- trees "capture" and "inject" of MC2MuCapture -> csv, for mu_capture_analysis.py (no PyROOT needed).
// Usage: root -l -b -q 'Offline/MC2Reco/scripts/capture_to_csv.C("pass.root","out/prefix")'
#include <TFile.h>
#include <TTree.h>
#include <cstdio>
void dumpTree(TTree* t, const char* fn, bool inj) {
  int kind, cry, disk, ns, nMi = 0; float amp, t1, maxS, tau, A2, Q, miDt = -1, miE = 0, injD = 0, injA2 = 0;
  t->SetBranchAddress("kind", &kind); t->SetBranchAddress("cry", &cry); t->SetBranchAddress("disk", &disk); t->SetBranchAddress("ns", &ns);
  t->SetBranchAddress("amp", &amp); t->SetBranchAddress("t1", &t1); t->SetBranchAddress("maxS", &maxS); t->SetBranchAddress("tau", &tau); t->SetBranchAddress("A2", &A2); t->SetBranchAddress("Q", &Q);
  t->SetBranchAddress("nMichel", &nMi); t->SetBranchAddress("michelDt", &miDt); t->SetBranchAddress("michelE", &miE);
  if (inj) { t->SetBranchAddress("injDelay", &injD); t->SetBranchAddress("injA2", &injA2); }
  FILE* f = fopen(fn, "w"); fprintf(f, "kind,cry,disk,ns,amp,t1,maxS,tau,A2,Q,nMichel,michelDt,michelE,injDelay,injA2\n");
  for (Long64_t i = 0; i < t->GetEntries(); ++i) { t->GetEntry(i);
    fprintf(f, "%d,%d,%d,%d,%.2f,%.2f,%.3f,%.2f,%.2f,%.5f,%d,%.1f,%.2f,%.2f,%.2f\n", kind, cry, disk, ns, amp, t1, maxS, tau, A2, Q, nMi, miDt, miE, injD, injA2); }
  fclose(f); printf("%s: %lld rows\n", fn, t->GetEntries());
}
void capture_to_csv(const char* in, const char* out) {
  TFile f(in); TTree* c = (TTree*)f.Get("capture/capture"); TTree* j = (TTree*)f.Get("capture/inject");
  if (!c || !j) { printf("ERROR: no capture trees in %s\n", in); return; }
  dumpTree(c, Form("%s_capture.csv", out), false); dumpTree(j, Form("%s_inject.csv", out), true);
}
