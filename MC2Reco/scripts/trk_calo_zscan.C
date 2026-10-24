// trk_calo_zscan.C: where are the calorimeter disks along z, according to the tracker lines?
//
// Input: the "proj" tree of MC2TrkCaloProjection: per track and disk, the tracker line at the disk plane of the job
// (x, y, z, t, direction ux, uy, uz, speed) and all the calorimeter hits of the event (hitDisk, hitX, hitY, hitT).
//
// Method: for each disk and each trial shift z of its plane, the line is moved by z along itself (path length
// s = z / uz) and every calorimeter hit of that disk is compared with it. Among the in-time pairs
// (|t_hit - t_line - dtOffset| < dtWindow), the fraction with the hit crystal within rMatch of the line is maximal at
// the effective z of the disk. The best shift comes from a parabola fitted around the maximum (+-4 trial points); it
// is the shift with respect to the plane used in the job, i.e. the diskZOffset to add in the fcl.
//
// Usage (any ROOT; compiled with + is faster):
//
//   root -l -b -q 'Offline/MC2Reco/scripts/trk_calo_zscan.C("trk_calo_<run>.root","out/zscan_<run>")'
//
// Arguments: in = TFileService output with trkCalo/proj; out = prefix of the output files;
//            dtOffset (2 ns), dtWindow (20 ns) = in-time window; rMatch (100 mm) = distance hit crystal - line;
//            zLo, zHi, zStep (-1500, 1500, 25 mm) = trial shifts.
//
// Output: <out>.png (fraction vs shift, both disks) and <out>.txt (the result lines, also printed):
//   disk <d>: plane of the job z = ... mm; best shift ... mm (maximum ... at ...); pairs ...
//   diskZOffset to add to the one of the job: [ ..., ... ]; disk separation changes by ... mm

#include <TAxis.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TFile.h>
#include <TGraph.h>
#include <TLegend.h>
#include <TLine.h>
#include <TTree.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <vector>

void trk_calo_zscan(const char* in, const char* out, double dtOffset = 2., double dtWindow = 20., double rMatch = 100.,
                    double zLo = -1500., double zHi = 1500., double zStep = 25.) {
  TFile inputFile(in);
  TTree* tree = (TTree*)inputFile.Get("trkCalo/proj");
  if (!tree) {
    printf("ERROR: no trkCalo/proj in %s\n", in);
    return;
  }

  // branches: the line at the disk plane of the job, and the calorimeter hits of the event
  int disk;
  float x, y, z, lineT, ux, uy, uz, speed;
  std::vector<int>* hitDisk = nullptr;
  std::vector<float>* hitX = nullptr;
  std::vector<float>* hitY = nullptr;
  std::vector<float>* hitT = nullptr;
  tree->SetBranchAddress("disk", &disk);
  tree->SetBranchAddress("x", &x);
  tree->SetBranchAddress("y", &y);
  tree->SetBranchAddress("z", &z);
  tree->SetBranchAddress("t", &lineT);
  tree->SetBranchAddress("ux", &ux);
  tree->SetBranchAddress("uy", &uy);
  tree->SetBranchAddress("uz", &uz);
  tree->SetBranchAddress("speed", &speed);
  tree->SetBranchAddress("hitDisk", &hitDisk);
  tree->SetBranchAddress("hitX", &hitX);
  tree->SetBranchAddress("hitY", &hitY);
  tree->SetBranchAddress("hitT", &hitT);

  // trial shifts of the disk plane
  const int nz = int((zHi - zLo) / zStep) + 1;
  std::vector<double> zShifts(nz);
  for (int i = 0; i < nz; ++i) {
    zShifts[i] = zLo + i * zStep;
  }

  // per disk and trial shift: in-time pairs, and in-time pairs within rMatch
  std::vector<double> nInTime[2] = {std::vector<double>(nz, 0), std::vector<double>(nz, 0)};
  std::vector<double> nMatched[2] = {std::vector<double>(nz, 0), std::vector<double>(nz, 0)};
  double planeZ[2] = {NAN, NAN};   // plane of the job (z of the first line of each disk)
  long nPairs[2] = {0, 0};         // (line, hit of the same disk) pairs, all times

  for (Long64_t entry = 0; entry < tree->GetEntries(); ++entry) {
    tree->GetEntry(entry);
    if (disk < 0 || disk > 1 || std::abs(uz) < 1e-6) {
      continue;
    }
    if (std::isnan(planeZ[disk])) {
      planeZ[disk] = z;
    }
    for (size_t j = 0; j < hitDisk->size(); ++j) {
      if ((*hitDisk)[j] != disk) {
        continue;
      }
      for (int i = 0; i < nz; ++i) {
        const double s = zShifts[i] / uz;   // s: path length to the trial plane
        const double xl = x + ux * s;
        const double yl = y + uy * s;
        const double tl = lineT + s / speed;
        if (std::abs((*hitT)[j] - tl - dtOffset) > dtWindow) {
          continue;
        }
        nInTime[disk][i] += 1;
        if (std::hypot((*hitX)[j] - xl, (*hitY)[j] - yl) < rMatch) {
          nMatched[disk][i] += 1;
        }
      }
      ++nPairs[disk];
    }
  }

  // per disk: fraction vs shift, parabola around the maximum, result lines on screen and in <out>.txt
  FILE* resultFile = fopen(Form("%s.txt", out), "w");
  TCanvas canvas("c", "", 1100, 500);
  TGraph* graphs[2];
  double bestShift[2] = {NAN, NAN};
  for (int d = 0; d < 2; ++d) {
    std::vector<double> fraction(nz);
    for (int i = 0; i < nz; ++i) {
      fraction[i] = nInTime[d][i] > 0 ? nMatched[d][i] / nInTime[d][i] : 0;
    }
    const int iMax = std::max_element(fraction.begin(), fraction.end()) - fraction.begin();

    graphs[d] = new TGraph(nz, zShifts.data(), fraction.data());
    graphs[d]->SetLineColor(d == 0 ? 4 : 2);
    graphs[d]->SetMarkerColor(d == 0 ? 4 : 2);
    graphs[d]->SetMarkerStyle(20);
    graphs[d]->SetMarkerSize(0.5);

    TF1 parabola(Form("pf%d", d), "pol2", zShifts[std::max(iMax - 4, 0)], zShifts[std::min(iMax + 4, nz - 1)]);
    graphs[d]->Fit(&parabola, "QRN0");
    bestShift[d] = -parabola.GetParameter(1) / (2 * parabola.GetParameter(2));

    for (FILE* output : {stdout, resultFile}) {
      if (output) {
        fprintf(output,
                "disk %d: plane of the job z = %.1f mm; best shift %+.0f mm (maximum %.3f at %+.0f); pairs %ld\n", d,
                planeZ[d], bestShift[d], fraction[iMax], zShifts[iMax], nPairs[d]);
      }
    }
  }
  for (FILE* output : {stdout, resultFile}) {
    if (output) {
      fprintf(output,
              "diskZOffset to add to the one of the job: [ %.0f, %.0f ]; disk separation changes by %+.0f mm\n",
              bestShift[0], bestShift[1], bestShift[1] - bestShift[0]);
    }
  }
  if (resultFile) {
    fclose(resultFile);
  }

  // figure: both disks, y range up to 1.1 x the highest fraction
  graphs[0]->SetTitle(Form("tracker lines vs calo hits;shift of the disk plane with respect to the job [mm];"
                           "fraction of in-time pairs within %.0f mm",
                           rMatch));
  const double maxFraction0 = *std::max_element(graphs[0]->GetY(), graphs[0]->GetY() + nz);
  const double maxFraction1 = *std::max_element(graphs[1]->GetY(), graphs[1]->GetY() + nz);
  graphs[0]->GetYaxis()->SetRangeUser(0, 1.1 * std::max(maxFraction0, maxFraction1));
  graphs[0]->Draw("APL");
  graphs[1]->Draw("PL same");
  TLegend legend(0.13, 0.75, 0.6, 0.88);
  legend.AddEntry(graphs[0], Form("disk 0: best shift %+.0f mm", bestShift[0]), "pl");
  legend.AddEntry(graphs[1], Form("disk 1: best shift %+.0f mm", bestShift[1]), "pl");
  legend.Draw();
  canvas.SaveAs(Form("%s.png", out));
}
