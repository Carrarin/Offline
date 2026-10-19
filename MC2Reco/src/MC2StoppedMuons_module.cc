// MC2StoppedMuons: muons that stop in the calorimeter and their delayed deposits (decay electrons of mu+, Michel spectrum),
// on calibrated CaloHits. Port of the mc2tools macro cosmic_delayed_coincidence.C (same selection, same tree), so that
// michel_scale.py reads the output as it is; in addition the Michel spectra per disk and the lifetime fit at the end of the job.
//
// Crystal : CaloHit with both SiPMs (nSiPMs == 2) and the amplitude of each SiPM >= minAmp ADC counts (energy / CalCalib
//           ADC2MeV); time = hit time (our template, after the time offsets of the conditions); amplitude = mean of the two SiPM
//           amplitudes; energy = CaloHit energy (MeV, mean of the two SiPMs, from the energy table of the conditions);
//           position = crystal x, y (tracker frame of the Offline geometry). More than one hit per crystal per event is allowed.
// Cluster : crystals in time order, a new cluster after a gap > clusGap ns.
// Pair    : prompt cluster P and later cluster D of the same event window, dt = t_D - t_P in [dtMin, dtMax] us, P starting after
//           tEdge us from the window start and dtMax before its end (window length windowNs): full acceptance in dt;
//           near (cls 1) = some crystal of D within dNear mm (same disk) of some crystal of P; far (cls 2) = all crystals of D farther
//           than dFar mm (accidental control). Same branches as mc2tools pairs.root (entry = art event number).
// Michel  : per disk of the delayed cluster, energy eD of the near pairs with dt in [michelDtMin, dtMax] us (mu+ only: the mu- are
//           captured in < 0.1 us in CsI), and of the far pairs; accidentals subtracted with k = near/far in dt > dtNorm us (flat).
// Lifetime: binned Poisson likelihood fit N0 exp(-t/tau) + B of the near pairs, dt in [fitMin, fitMax] us, all and per disk;
//           results printed and written to summaryFile.
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art_root_io/TFileService.h"
#include "fhiclcpp/types/Atom.h"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/CalorimeterGeom/inc/Calorimeter.hh"
#include "Offline/CaloConditions/inc/CalCalib.hh"
#include "Offline/ProditionsService/inc/ProditionsHandle.hh"
#include "Offline/RecoDataProducts/inc/CaloHit.hh"
#include "Offline/RecoDataProducts/inc/CaloRecoDigi.hh"
#include "TF1.h"
#include "TH1F.h"
#include "TTree.h"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <map>
#include <vector>

namespace mu2e {
  class MC2StoppedMuons : public art::EDAnalyzer {
    public:
      struct Config {
        using Name = fhicl::Name; using Comment = fhicl::Comment;
        fhicl::Atom<art::InputTag> hitTag     { Name("caloHitCollection"), Comment("calibrated calorimeter hits"), "MC2CaloHitMaker" };
        fhicl::Atom<double> minAmp      { Name("minAmp"), Comment("min amplitude of each SiPM [ADC]"), 100. };
        fhicl::Atom<double> clusGap     { Name("clusGap"), Comment("time gap that starts a new cluster [ns]"), 30. };
        fhicl::Atom<double> dtMin       { Name("dtMin"), Comment("min delay [us]"), 0.1 };
        fhicl::Atom<double> dtMax       { Name("dtMax"), Comment("max delay [us]"), 20. };
        fhicl::Atom<double> tEdge       { Name("tEdge"), Comment("min prompt time from the window start [us]"), 1. };
        fhicl::Atom<double> windowNs    { Name("windowNs"), Comment("event window length [ns]"), 97495. };
        fhicl::Atom<double> dNear       { Name("dNear"), Comment("near: max distance prompt-delayed crystal [mm]"), 80. };
        fhicl::Atom<double> dFar        { Name("dFar"), Comment("far: min distance [mm]"), 300. };
        fhicl::Atom<double> michelDtMin { Name("michelDtMin"), Comment("min delay for the Michel spectra [us]"), 1. };
        fhicl::Atom<double> dtNorm      { Name("dtNorm"), Comment("near/far normalisation for dt > dtNorm [us]"), 15. };
        fhicl::Atom<double> fitMin      { Name("fitMin"), Comment("lifetime fit from [us]"), 1. };
        fhicl::Atom<double> fitMax      { Name("fitMax"), Comment("lifetime fit to [us]"), 20. };
        fhicl::Atom<std::string> summaryFile { Name("summaryFile"), Comment("text summary"), "stopped_muons_summary.txt" };
      };
      using Parameters = art::EDAnalyzer::Table<Config>;
      explicit MC2StoppedMuons(const Parameters& conf);
      void beginJob() override;
      void analyze(const art::Event& event) override;
      void endJob() override;
    private:
      struct Cry { int disk; float x, y; double t, a, e; int id; };
      art::InputTag tag_; double minAmp_, gap_, dtMin_, dtMax_, tEdge_, win_, dNear_, dFar_, mDt_, dtNorm_, fitMin_, fitMax_; std::string sumFile_;
      ProditionsHandle<CalCalib> calibH_;
      TTree* t_ = nullptr;
      long long entry_; float dt_, dmin_, tP_, aP_, aD_, aDmax_, xD_, yD_, xLo_, yLo_, xHi_, yHi_, xP_, yP_, eP_, eD_, eDmax_; int nP_, nD_, disk_, same_, cls_, nCl_, diskP_;
      TH1F *hDt_[3][2] = {{nullptr}}, *hE_[2][2] = {{nullptr}}, *hESub_[2] = {nullptr, nullptr}, *hDtNorm_[2] = {nullptr, nullptr};
      long nEv_ = 0, nCry_ = 0, nClus_ = 0, nPrompt_ = 0;
  };

  MC2StoppedMuons::MC2StoppedMuons(const Parameters& c) : art::EDAnalyzer(c), tag_(c().hitTag()), minAmp_(c().minAmp()), gap_(c().clusGap()),
    dtMin_(c().dtMin()), dtMax_(c().dtMax()), tEdge_(c().tEdge()), win_(c().windowNs()), dNear_(c().dNear()), dFar_(c().dFar()), mDt_(c().michelDtMin()),
    dtNorm_(c().dtNorm()), fitMin_(c().fitMin()), fitMax_(c().fitMax()), sumFile_(c().summaryFile()) { consumes<CaloHitCollection>(tag_); }

  void MC2StoppedMuons::beginJob() {
    art::ServiceHandle<art::TFileService> tfs;
    t_ = tfs->make<TTree>("pairs", "delayed coincidences (branches of mc2tools cosmic_delayed_coincidence.C)");
    t_->Branch("entry", &entry_); t_->Branch("dt", &dt_); t_->Branch("dmin", &dmin_); t_->Branch("tP", &tP_); t_->Branch("aP", &aP_);
    t_->Branch("aD", &aD_); t_->Branch("aDmax", &aDmax_); t_->Branch("nP", &nP_); t_->Branch("nD", &nD_); t_->Branch("disk", &disk_); t_->Branch("sameCry", &same_);
    t_->Branch("cls", &cls_); t_->Branch("nCl", &nCl_); t_->Branch("xD", &xD_); t_->Branch("yD", &yD_); t_->Branch("xLo", &xLo_); t_->Branch("yLo", &yLo_);
    t_->Branch("xHi", &xHi_); t_->Branch("yHi", &yHi_); t_->Branch("xP", &xP_); t_->Branch("yP", &yP_); t_->Branch("diskP", &diskP_);
    t_->Branch("eP", &eP_); t_->Branch("eD", &eD_); t_->Branch("eDmax", &eDmax_);
    const char* cn[3] = {"all", "disk0", "disk1"};
    for (int k = 0; k < 3; ++k) {
      hDt_[k][0] = tfs->make<TH1F>(Form("hDtNear_%s", cn[k]), Form("near pairs, %s;delay [#mus];pairs / 0.2 #mus", cn[k]), 100, 0, 20);
      hDt_[k][1] = tfs->make<TH1F>(Form("hDtFar_%s", cn[k]), Form("far pairs, %s;delay [#mus];pairs / 0.2 #mus", cn[k]), 100, 0, 20);
    }
    for (int d = 0; d < 2; ++d) {
      hE_[d][0] = tfs->make<TH1F>(Form("hMichelNear_disk%d", d), Form("disk %d, near pairs, delay > %.1f #mus;delayed energy [MeV];pairs / 2 MeV", d, mDt_), 50, 0, 100);
      hE_[d][1] = tfs->make<TH1F>(Form("hMichelFar_disk%d", d), Form("disk %d, far pairs, delay > %.1f #mus;delayed energy [MeV];pairs / 2 MeV", d, mDt_), 50, 0, 100);
      hESub_[d] = tfs->make<TH1F>(Form("hMichel_disk%d", d), Form("disk %d, near - k far (accidentals subtracted);delayed energy [MeV];pairs / 2 MeV", d), 50, 0, 100);
      hDtNorm_[d] = tfs->make<TH1F>(Form("hNorm_disk%d", d), Form("disk %d: near and far pairs with delay > %.0f #mus;0 near, 1 far;pairs", d, dtNorm_), 2, -0.5, 1.5);
    }
  }

  void MC2StoppedMuons::analyze(const art::Event& event) {
    ++nEv_;
    const Calorimeter& cal = *GeomHandle<Calorimeter>(); const CalCalib& calib = calibH_.get(event.id());
    const auto& hits = *event.getValidHandle<CaloHitCollection>(tag_);
    std::vector<Cry> cr;
    for (const auto& h : hits) {
      if (h.nSiPMs() < 2) continue;
      double a = 0; bool ok = true; int n = 0;
      for (const auto& rd : h.recoCaloDigis()) { const double a2m = calib.ADC2MeV(rd->SiPMID()); const double am = a2m > 0 ? rd->energyDep() / a2m : 0;
        if (am < minAmp_) { ok = false; break; }
        a += am; ++n; }
      if (!ok || n == 0) continue;
      const auto& c = cal.crystal(h.crystalID()); const auto p = cal.mu2eToTracker(c.position());
      cr.push_back({int(c.diskID()), float(p.x()), float(p.y()), h.time(), a / n, h.energyDep(), h.crystalID()});
    }
    nCry_ += cr.size(); if (cr.size() < 2) return;
    std::sort(cr.begin(), cr.end(), [](const Cry& p, const Cry& q) { return p.t < q.t; });
    std::vector<std::pair<int, int>> cl; int s0 = 0;
    for (size_t h = 1; h <= cr.size(); ++h) if (h == cr.size() || cr[h].t - cr[h - 1].t > gap_) { cl.push_back({s0, (int)h}); s0 = h; }
    nClus_ += cl.size(); if (cl.size() < 2) return;
    for (size_t p = 0; p < cl.size(); ++p) {
      const double tP = cr[cl[p].first].t; if (tP < tEdge_ * 1e3 || tP > win_ - dtMax_ * 1e3) continue; ++nPrompt_;
      double aP = 0, eP = 0; for (int h = cl[p].first; h < cl[p].second; ++h) { aP += cr[h].a; eP += cr[h].e; }
      for (size_t d = p + 1; d < cl.size(); ++d) {
        const double dt = (cr[cl[d].first].t - tP) * 1e-3; if (dt < dtMin_) continue; if (dt > dtMax_) break;
        double dmin = 1e9, aD = 0, aDm = 0, eD = 0, eDm = 0; int same = 0; const int disk = cr[cl[d].first].disk;
        for (int h = cl[d].first; h < cl[d].second; ++h) { aD += cr[h].a; aDm = std::max(aDm, cr[h].a); eD += cr[h].e; eDm = std::max(eDm, cr[h].e);
          for (int g = cl[p].first; g < cl[p].second; ++g) { if (cr[h].disk != cr[g].disk) continue; dmin = std::min(dmin, double(std::hypot(cr[h].x - cr[g].x, cr[h].y - cr[g].y))); if (cr[h].id == cr[g].id) same = 1; } }
        const int cls = dmin <= dNear_ ? 1 : (dmin > dFar_ ? 2 : 0); if (!cls) continue;
        double sx = 0, sy = 0, sw = 0; for (int h = cl[d].first; h < cl[d].second; ++h) { sx += cr[h].a * cr[h].x; sy += cr[h].a * cr[h].y; sw += cr[h].a; }
        int gLo = cl[p].first, gHi = cl[p].first; double px = 0, py = 0, pw = 0; std::map<int, int> nDisk;
        for (int g = cl[p].first; g < cl[p].second; ++g) { if (cr[g].y < cr[gLo].y) gLo = g; if (cr[g].y > cr[gHi].y) gHi = g; px += cr[g].a * cr[g].x; py += cr[g].a * cr[g].y; pw += cr[g].a; nDisk[cr[g].disk]++; }
        entry_ = event.event(); dt_ = dt; dmin_ = dmin; tP_ = tP * 1e-3; aP_ = aP; aD_ = aD; aDmax_ = aDm; nP_ = cl[p].second - cl[p].first; nD_ = cl[d].second - cl[d].first;
        disk_ = disk; same_ = same; cls_ = cls; nCl_ = cl.size(); xD_ = sx / sw; yD_ = sy / sw; xLo_ = cr[gLo].x; yLo_ = cr[gLo].y; xHi_ = cr[gHi].x; yHi_ = cr[gHi].y;
        xP_ = px / pw; yP_ = py / pw; diskP_ = nDisk.size() > 1 ? 2 : nDisk.begin()->first; eP_ = eP; eD_ = eD; eDmax_ = eDm; t_->Fill();
        const int ic = cls == 1 ? 0 : 1; hDt_[0][ic]->Fill(dt); if (disk == 0 || disk == 1) hDt_[1 + disk][ic]->Fill(dt);
        if ((disk == 0 || disk == 1) && dt >= mDt_) hE_[disk][ic]->Fill(eD);
        if ((disk == 0 || disk == 1) && dt >= dtNorm_) hDtNorm_[disk]->Fill(ic);
      }
    }
  }

  void MC2StoppedMuons::endJob() {
    FILE* fs = fopen(sumFile_.c_str(), "w");
    auto out = [&](const char* fmt, auto... a) { printf(fmt, a...); if (fs) fprintf(fs, fmt, a...); };
    out("MC2StoppedMuons: events %ld, crystals %ld, clusters %ld, prompt clusters in acceptance %ld\n", nEv_, nCry_, nClus_, nPrompt_);
    out("cuts: minAmp %.0f ADC, clusGap %.0f ns, dt %.2f-%.1f us, tEdge %.1f us, near <= %.0f mm, far > %.0f mm\n", minAmp_, gap_, dtMin_, dtMax_, tEdge_, dNear_, dFar_);
    const char* cn[3] = {"all", "disk 0", "disk 1"};
    for (int k = 0; k < 3; ++k) {
      TH1F* h = hDt_[k][0]; if (h->Integral() < 20) { out("lifetime %s: too few near pairs (%.0f)\n", cn[k], h->Integral()); continue; }
      TF1 f(Form("fexp%d", k), "[0]*exp(-x/[1])+[2]", fitMin_, fitMax_); f.SetParameters(std::max(1., h->GetBinContent(h->FindBin(fitMin_)) * 1.5), 2.2, std::max(0.1, h->GetBinContent(h->GetNbinsX())));
      f.SetParLimits(1, 0.1, 20.); h->Fit(&f, "LQR0");
      out("lifetime %s: tau = %.3f +- %.3f us, N0 = %.1f +- %.1f per bin, B = %.2f +- %.2f per bin, chi2/ndf = %.1f/%d, near pairs %.0f, far pairs %.0f\n", cn[k],
          f.GetParameter(1), f.GetParError(1), f.GetParameter(0), f.GetParError(0), f.GetParameter(2), f.GetParError(2), f.GetChisquare(), f.GetNDF(), h->Integral(), hDt_[k][1]->Integral());
    }
    for (int d = 0; d < 2; ++d) {
      const double nn = hDtNorm_[d]->GetBinContent(1), nf = hDtNorm_[d]->GetBinContent(2), k = nf > 0 ? nn / nf : 0;
      hESub_[d]->Add(hE_[d][0], hE_[d][1], 1., -k);
      out("Michel disk %d: near %.0f, far %.0f (delay > %.1f us), accidental scale k = %.3f (near/far for delay > %.0f us: %.0f/%.0f); after subtraction %.0f pairs, mean %.1f MeV\n",
          d, hE_[d][0]->Integral(), hE_[d][1]->Integral(), mDt_, k, dtNorm_, nn, nf, hESub_[d]->Integral(), hESub_[d]->GetMean());
    }
    if (fs) fclose(fs);
  }
}
DEFINE_ART_MODULE(mu2e::MC2StoppedMuons)
