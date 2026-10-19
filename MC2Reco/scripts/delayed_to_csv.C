// delayed_to_csv.C -- tree "delayed" of MC2MuCapture -> csv (no PyROOT needed).
// Usage: root -l -b -q 'Offline/MC2Reco/scripts/delayed_to_csv.C("pass.root","out.csv")'
#include <TFile.h>
#include <TTree.h>
#include <cstdio>
void delayed_to_csv(const char* in, const char* out) {
  TFile f(in); TTree* t = (TTree*)f.Get("capture/delayed"); if (!t) { printf("ERROR: no capture/delayed in %s\n", in); return; }
  int kind, cry, disk, nMi, ring, nb, sipm; float dt, e, dist, stopE;
  t->SetBranchAddress("kind", &kind); t->SetBranchAddress("cry", &cry); t->SetBranchAddress("disk", &disk); t->SetBranchAddress("nMichel", &nMi); t->SetBranchAddress("ring", &ring);
  t->SetBranchAddress("nbCry", &nb); t->SetBranchAddress("sipm", &sipm); t->SetBranchAddress("dt", &dt); t->SetBranchAddress("e", &e); t->SetBranchAddress("dist", &dist); t->SetBranchAddress("stopE", &stopE);
  FILE* o = fopen(out, "w"); fprintf(o, "kind,cry,disk,nMichel,ring,nbCry,sipm,dt,e,dist,stopE\n");
  for (Long64_t i = 0; i < t->GetEntries(); ++i) { t->GetEntry(i); fprintf(o, "%d,%d,%d,%d,%d,%d,%d,%.2f,%.3f,%.1f,%.2f\n", kind, cry, disk, nMi, ring, nb, sipm, dt, e, dist, stopE); }
  fclose(o); printf("written %s: %lld rows\n", out, t->GetEntries());
}
