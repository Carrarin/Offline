// capture_to_csv.C: trees "capture" and "inject" of MC2MuCapture -> csv, for mu_capture_analysis.py (no PyROOT needed).
//
// Usage (muse environment, or any ROOT):
//
//   root -l -b -q 'Offline/MC2Reco/scripts/capture_to_csv.C("pass.root","out/prefix")'
//
//   e.g. root -l -b -q 'Offline/MC2Reco/scripts/capture_to_csv.C("out/physics.root","out/physics")'
//
// Arguments: in = TFileService output of a job with the MC2MuCapture analyzer "capture" (trees capture/capture and
//            capture/inject); out = prefix of the csv files.
//
// Output: <out>_capture.csv and <out>_inject.csv, one row per entry, with the same header
//   kind,cry,disk,ns,amp,t1,maxS,tau,A2,Q,nMichel,michelDt,michelE,injDelay,injA2
// injDelay and injA2 exist only in the inject tree: in <out>_capture.csv they are 0. A branch missing in the file keeps
// its default (nMichel 0, michelDt -1, michelE 0). Prints the number of rows of each file.
//
// Called by run_machine.sh (step 6).

#include <TFile.h>
#include <TTree.h>

#include <cstdio>

// One tree -> one csv file; withInjection = read also the branches of the injected pulses (tree "inject")
void dumpTree(TTree* tree, const char* csvName, bool withInjection) {
  int kind, cry, disk, ns;
  int nMichel = 0;
  float amp, t1, maxS, tau, A2, Q;
  float michelDt = -1;
  float michelE = 0;
  float injDelay = 0;
  float injA2 = 0;
  tree->SetBranchAddress("kind", &kind);
  tree->SetBranchAddress("cry", &cry);
  tree->SetBranchAddress("disk", &disk);
  tree->SetBranchAddress("ns", &ns);
  tree->SetBranchAddress("amp", &amp);
  tree->SetBranchAddress("t1", &t1);
  tree->SetBranchAddress("maxS", &maxS);
  tree->SetBranchAddress("tau", &tau);
  tree->SetBranchAddress("A2", &A2);
  tree->SetBranchAddress("Q", &Q);
  tree->SetBranchAddress("nMichel", &nMichel);
  tree->SetBranchAddress("michelDt", &michelDt);
  tree->SetBranchAddress("michelE", &michelE);
  if (withInjection) {
    tree->SetBranchAddress("injDelay", &injDelay);
    tree->SetBranchAddress("injA2", &injA2);
  }

  FILE* csv = fopen(csvName, "w");
  fprintf(csv, "kind,cry,disk,ns,amp,t1,maxS,tau,A2,Q,nMichel,michelDt,michelE,injDelay,injA2\n");
  for (Long64_t i = 0; i < tree->GetEntries(); ++i) {
    tree->GetEntry(i);
    fprintf(csv, "%d,%d,%d,%d,%.2f,%.2f,%.3f,%.2f,%.2f,%.5f,%d,%.1f,%.2f,%.2f,%.2f\n", kind, cry, disk, ns, amp, t1,
            maxS, tau, A2, Q, nMichel, michelDt, michelE, injDelay, injA2);
  }
  fclose(csv);
  printf("%s: %lld rows\n", csvName, tree->GetEntries());
}

void capture_to_csv(const char* in, const char* out) {
  TFile inputFile(in);
  TTree* captureTree = (TTree*)inputFile.Get("capture/capture");
  TTree* injectTree = (TTree*)inputFile.Get("capture/inject");
  if (!captureTree || !injectTree) {
    printf("ERROR: no capture trees in %s\n", in);
    return;
  }
  dumpTree(captureTree, Form("%s_capture.csv", out), false);
  dumpTree(injectTree, Form("%s_inject.csv", out), true);
}
