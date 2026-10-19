// trk_calo_zscan.C -- where are the calorimeter disks along z, according to the tracker lines?
// Input: the "proj" tree of MC2TrkCaloProjection (the line at the plane: x, y, z, t, ux, uy, uz, speed; all calo hits of the event).
// For each disk and each trial z, the line is moved to that z and every calo hit of that disk is compared with it: the fraction of
// in-time pairs (|t_hit - t_line - dtOffset| < dtWindow) with the hit crystal within rMatch of the line is maximal at the effective z
// of the disk. Best z from a parabola around the maximum; result = shift with respect to the plane used in the job (diskZOffset
// to add in the fcl).
// Output: <out>.png (fraction vs z, both disks), <out>.txt (results). Usage:
//   root -l -b -q 'Offline/MC2Reco/scripts/trk_calo_zscan.C("trk_calo_<run>.root","out/zscan_<run>")'
#include <TFile.h>
#include <TTree.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLine.h>
#include <TF1.h>
#include <TAxis.h>
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <vector>
void trk_calo_zscan(const char* in, const char* out, double dtOffset = 2., double dtWindow = 20., double rMatch = 100.,
                    double zLo = -1500., double zHi = 1500., double zStep = 25.) {
  TFile f(in); TTree* t = (TTree*)f.Get("trkCalo/proj");
  if (!t) { printf("ERROR: no trkCalo/proj in %s\n", in); return; }
  int disk; float x, y, z, tt, ux, uy, uz, speed; std::vector<int>* hd = nullptr; std::vector<float> *hx = nullptr, *hy = nullptr, *ht = nullptr;
  t->SetBranchAddress("disk", &disk); t->SetBranchAddress("x", &x); t->SetBranchAddress("y", &y); t->SetBranchAddress("z", &z); t->SetBranchAddress("t", &tt);
  t->SetBranchAddress("ux", &ux); t->SetBranchAddress("uy", &uy); t->SetBranchAddress("uz", &uz); t->SetBranchAddress("speed", &speed);
  t->SetBranchAddress("hitDisk", &hd); t->SetBranchAddress("hitX", &hx); t->SetBranchAddress("hitY", &hy); t->SetBranchAddress("hitT", &ht);
  const int nz = int((zHi - zLo) / zStep) + 1; std::vector<double> zs(nz); for (int i = 0; i < nz; ++i) zs[i] = zLo + i * zStep;
  std::vector<double> nin[2] = {std::vector<double>(nz, 0), std::vector<double>(nz, 0)}, nok[2] = {std::vector<double>(nz, 0), std::vector<double>(nz, 0)};
  double z0[2] = {NAN, NAN};
  long npair[2] = {0, 0};
  for (Long64_t e = 0; e < t->GetEntries(); ++e) {
    t->GetEntry(e); if (disk < 0 || disk > 1 || std::abs(uz) < 1e-6) continue; if (std::isnan(z0[disk])) z0[disk] = z;
    for (size_t j = 0; j < hd->size(); ++j) {
      if ((*hd)[j] != disk) continue;
      for (int i = 0; i < nz; ++i) {
        const double s = zs[i] / uz; const double xl = x + ux * s, yl = y + uy * s, tl = tt + s / speed;  // s: path length to the trial plane
        if (std::abs((*ht)[j] - tl - dtOffset) > dtWindow) continue;
        nin[disk][i] += 1; if (std::hypot((*hx)[j] - xl, (*hy)[j] - yl) < rMatch) nok[disk][i] += 1;
      }
      ++npair[disk];
    }
  }
  FILE* fo = fopen(Form("%s.txt", out), "w");
  TCanvas c("c", "", 1100, 500); TGraph* g[2]; double best[2] = {NAN, NAN};
  for (int d = 0; d < 2; ++d) {
    std::vector<double> fr(nz); for (int i = 0; i < nz; ++i) fr[i] = nin[d][i] > 0 ? nok[d][i] / nin[d][i] : 0;
    const int im = std::max_element(fr.begin(), fr.end()) - fr.begin();
    g[d] = new TGraph(nz, zs.data(), fr.data()); g[d]->SetLineColor(d == 0 ? 4 : 2); g[d]->SetMarkerColor(d == 0 ? 4 : 2); g[d]->SetMarkerStyle(20); g[d]->SetMarkerSize(0.5);
    TF1 pf(Form("pf%d", d), "pol2", zs[std::max(im - 4, 0)], zs[std::min(im + 4, nz - 1)]); g[d]->Fit(&pf, "QRN0");
    best[d] = -pf.GetParameter(1) / (2 * pf.GetParameter(2));
    for (FILE* o : {stdout, fo}) if (o) fprintf(o, "disk %d: plane of the job z = %.1f mm; best shift %+.0f mm (maximum %.3f at %+.0f); pairs %ld\n", d, z0[d], best[d], fr[im], zs[im], npair[d]);
  }
  for (FILE* o : {stdout, fo}) if (o) fprintf(o, "diskZOffset to add to the one of the job: [ %.0f, %.0f ]; disk separation changes by %+.0f mm\n", best[0], best[1], best[1] - best[0]);
  if (fo) fclose(fo);
  g[0]->SetTitle(Form("tracker lines vs calo hits;shift of the disk plane with respect to the job [mm];fraction of in-time pairs within %.0f mm", rMatch));
  g[0]->GetYaxis()->SetRangeUser(0, 1.1 * std::max(*std::max_element(g[0]->GetY(), g[0]->GetY() + nz), *std::max_element(g[1]->GetY(), g[1]->GetY() + nz)));
  g[0]->Draw("APL"); g[1]->Draw("PL same");
  TLegend l(0.13, 0.75, 0.6, 0.88); l.AddEntry(g[0], Form("disk 0: best shift %+.0f mm", best[0]), "pl"); l.AddEntry(g[1], Form("disk 1: best shift %+.0f mm", best[1]), "pl"); l.Draw();
  c.SaveAs(Form("%s.png", out));
}
