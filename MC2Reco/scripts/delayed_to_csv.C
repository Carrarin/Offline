// delayed_to_csv.C: tree "delayed" of MC2MuCapture -> csv, for neighbour_delayed_analysis.py (no PyROOT needed).
//
// Usage (muse environment, or any ROOT):
//
//   root -l -b -q 'Offline/MC2Reco/scripts/delayed_to_csv.C("pass.root","out.csv")'
//
//   e.g. root -l -b -q 'Offline/MC2Reco/scripts/delayed_to_csv.C("out/physics.root","out/delayed.csv")'
//
// Arguments: in = TFileService output of a job with the MC2MuCapture analyzer "capture" (tree capture/delayed);
//            out = csv file to write.
//
// Output: <out>, one row per entry of the tree, with the header
//   kind,cry,disk,nMichel,ring,nbCry,sipm,dt,e,dist,stopE
// (same names as the branches; dt, e, dist, stopE printed with 2, 3, 1, 2 decimals). Prints the number of rows.
//
// Called by run_machine.sh (step 6).

#include <TFile.h>
#include <TTree.h>

#include <cstdio>

void delayed_to_csv(const char* in, const char* out) {
  TFile inputFile(in);
  TTree* tree = (TTree*)inputFile.Get("capture/delayed");
  if (!tree) {
    printf("ERROR: no capture/delayed in %s\n", in);
    return;
  }

  // branches of the tree
  int kind, cry, disk, nMichel, ring, nbCry, sipm;
  float dt, e, dist, stopE;
  tree->SetBranchAddress("kind", &kind);
  tree->SetBranchAddress("cry", &cry);
  tree->SetBranchAddress("disk", &disk);
  tree->SetBranchAddress("nMichel", &nMichel);
  tree->SetBranchAddress("ring", &ring);
  tree->SetBranchAddress("nbCry", &nbCry);
  tree->SetBranchAddress("sipm", &sipm);
  tree->SetBranchAddress("dt", &dt);
  tree->SetBranchAddress("e", &e);
  tree->SetBranchAddress("dist", &dist);
  tree->SetBranchAddress("stopE", &stopE);

  // one csv row per entry
  FILE* csv = fopen(out, "w");
  fprintf(csv, "kind,cry,disk,nMichel,ring,nbCry,sipm,dt,e,dist,stopE\n");
  for (Long64_t i = 0; i < tree->GetEntries(); ++i) {
    tree->GetEntry(i);
    fprintf(csv, "%d,%d,%d,%d,%d,%d,%d,%.2f,%.3f,%.1f,%.2f\n", kind, cry, disk, nMichel, ring, nbCry, sipm, dt, e, dist,
            stopE);
  }
  fclose(csv);
  printf("written %s: %lld rows\n", out, tree->GetEntries());
}
