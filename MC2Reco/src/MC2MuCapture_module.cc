// MC2MuCapture: search for the nuclear capture of stopped mu- in CsI: a second pulse in the crystal where a cosmic muon stops.
// A detector, not a time fit: the first pulse is removed and what is left is searched for "something with the shape of a pulse".
//
// Samples (as mc2tools capture_double_pulse.C), from the calorimeter tracks of the job (MC2CaloCosmicTrack):
//   STOP    (kind 0): one-disk tracks with >= minCrystals crystals; the last crystal along the downward direction of the track in the disk
//                     plane, if the crystal one pitch (34.3 mm) further along that direction exists and has no CaloDigi on either SiPM;
//   THROUGH (kind 1): tracks with >= minThrough crystals, the crystal in the middle (the muon went through: no capture there).
// Waveforms: the CaloDigis of the two SiPMs of the crystal (through the CaloHit of the track); pedestal = mean of the first 4 samples;
//   rejected if shorter than minSamples, saturated (raw >= satADC) or with amplitude < minAmp ADC.
// Template: built at the end of the job from the THROUGH waveforms of the run (amplitude 300-1800 ADC), aligned at the 50 % crossing,
//   iterated nIterTemplate times (as mc2tools); 0.5 ns bins over [-60, xMax] ns.
// Detection, per crystal:
//   1. first pulse of each SiPM: amplitude and time from a template fit on [-20, +40] ns around the 50 % crossing (leading edge and peak);
//   2. residual r = y - A1 T(t - t1) of each SiPM, normalised to its A1, the two SiPMs averaged on the common sample grid (t0 of the
//      CaloDigis): the light of a capture reaches both SiPMs, the electronic noise does not;
//   3. matched filter S(tau) = sum r T(t - t1 - tau) / (sigma_r sqrt(sum T^2)), tau scanned from dtMin to the end of the waveform - 30 ns
//      in 2.5 ns steps; maxS and its tau (delay of the second pulse), A2 = amplitude of the second pulse at tau [ADC, mean of the SiPMs];
//      sigma_r from the noise per sample (noiseADC, white noise assumed in this version: maxS is a ranking, the threshold comes from the
//      THROUGH controls);
//   4. tail charge excess Q = sum of r over t - t1 >= tailStart ns, in units of the first pulse (independent of the filter and of tau).
// Efficiency (injection): every THROUGH crystal is processed a second time with an injected pulse a T(t - t1 - d) added to both SiPMs,
//   d uniform in [dtMin, injMaxDelay] ns, A2 uniform in [injMinA2, injMaxA2] ADC (fixed seed); tree "inject" with the true and found
//   values -> efficiency(delay, A2) for a given cut on maxS.
// mu+ veto: for every candidate, the CaloHits of the event in the same disk within vetoDist mm of the crystal and later than its hit by
//   [vetoMin, vetoMax] ns (decay electron of a mu+, tau = 2.2 us; below vetoMin the capture itself lives): nMichel, the delay of the
//   first one and their summed energy are stored; the candidates are not removed (compare with and without the veto).
// Delayed deposits in the other crystals (tree "delayed"): for every candidate, the CaloRecoDigis of the event (single SiPMs, time with
//   the conditions' offsets, energy in MeV) in the same disk, not on the track and not the candidate crystal, with dt = t_digi - t_hit of
//   the candidate in [nbDtMin, nbDtMax] ns; ring 1 = within nbNear mm of the candidate crystal, ring 2 = between nbFarMin and nbFarMax
//   mm (accidental control). No first pulse to remove, and no limit from the waveform length: capture (~85 ns), mu+ decays (2.2 us) and
//   accidentals (flat) in one distribution.
// Delayed cluster (tree "cluster"; decay electrons of the mu- in the muonic atom (DIO) and of early mu+, in the tail of the stop):
//   members = the second pulse of the candidate crystal (found when maxS > secS; energy = A2 x mean ADC2MeV of its two SiPMs,
//   time = t_hit + tau) and the crystals within clusR mm in the same disk, not on the track, with a pulse on both SiPMs
//   (CaloRecoDigis of the two SiPMs within pairDt ns; energy = mean of the two, MeV from the conditions; time = mean), dt from the
//   hit of the candidate in [dtMin, clusDtMax] ns; members in time order, a new cluster after a gap > clusGap ns; one entry per
//   cluster with E >= clusEmin MeV: dt (energy weighted), E, the part in the candidate crystal, number of crystals, inWf (the cluster
//   time is inside the candidate waveform, i.e. its second pulse could have been seen). The threshold on E for the physics is
//   eCut (counts printed at the end); THROUGH candidates (kind 1) give the same clusters without a stop: tail and accidentals.
// Output (TFileService): tree "capture" (kind, crystal, disk, amp, t1, maxS, tau, A2, Q, ns, nMichel, michelDt, michelE, run/event),
//   tree "inject", the template.
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art_root_io/TFileService.h"
#include "fhiclcpp/types/Atom.h"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/CalorimeterGeom/inc/Calorimeter.hh"
#include "Offline/CaloConditions/inc/CalCalib.hh"
#include "Offline/ProditionsService/inc/ProditionsHandle.hh"
#include "Offline/DataProducts/inc/CaloConst.hh"
#include "Offline/RecoDataProducts/inc/CaloDigi.hh"
#include "Offline/RecoDataProducts/inc/CaloRecoDigi.hh"
#include "Offline/RecoDataProducts/inc/CaloHit.hh"
#include "Offline/MC2Reco/inc/MC2CaloCosmicTrack.hh"
#include "TH1D.h"
#include "TRandom3.h"
#include "TTree.h"
#include <algorithm>
#include <cmath>
#include <map>
#include <set>
#include <vector>

namespace mu2e {
  namespace {
    const double kNs = 5.0, XLO = -60., BW = 0.5;
    struct Tmpl {
      std::vector<double> v; int nbin = 0;
      double at(double x) const { const double u = (x - XLO) / BW - 0.5; if (u <= 0) return v.empty() ? 0 : v.front(); const int i = int(u); if (i >= nbin - 1) return 0.; const double f = u - i; return v[i] * (1 - f) + v[i + 1] * f; }
    };
    struct Wf { std::vector<float> y; double t0; float amp; double t50; };   // y normalised to the amplitude, t0 = time of sample 0 [ns]
    struct Dep { float dt, e; int cry; };   // a crystal deposit, dt from the hit of the candidate [ns], energy [MeV]
    struct Cand { int kind, cry, disk, run, sub, evt; Wf w[2]; int nMichel = 0; float michelDt = -1, michelE = 0; double a2m[2] = {0, 0}; double tEnd = 0, stopE = 0; std::vector<Dep> nb; };
    double t50of(const std::vector<float>& y) {
      const int im = int(std::max_element(y.begin(), y.end()) - y.begin());
      for (int q = im; q > 0; --q) if (y[q - 1] < 0.5 && y[q] >= 0.5) return (q - 1 + (0.5 - y[q - 1]) / (y[q] - y[q - 1])) * kNs;
      return im * kNs; }
    Tmpl build(const std::vector<const Wf*>& w, const std::vector<double>& ta, int nbin) {
      std::vector<double> s(nbin, 0), n(nbin, 0);
      for (size_t i = 0; i < w.size(); ++i) for (size_t q = 0; q < w[i]->y.size(); ++q) { const double x = q * kNs - ta[i]; const int b = int((x - XLO) / BW); if (b < 0 || b >= nbin) continue; s[b] += w[i]->y[q]; n[b] += 1; }
      std::vector<double> raw(nbin, 0); for (int b = 0; b < nbin; ++b) raw[b] = n[b] > 0 ? s[b] / n[b] : 0;
      Tmpl t; t.nbin = nbin; t.v.assign(nbin, 0);
      for (int b = 0; b < nbin; ++b) { double a = 0; int c = 0; for (int d = -2; d <= 2; ++d) if (b + d >= 0 && b + d < nbin && n[b + d] > 0) { a += raw[b + d]; ++c; } t.v[b] = c ? a / c : 0; }
      const double pk = *std::max_element(t.v.begin(), t.v.end()); if (pk > 0) for (auto& x : t.v) x /= pk;
      const int im = int(std::max_element(t.v.begin(), t.v.end()) - t.v.begin()); double x50 = 0;
      for (int b = im; b > 0; --b) if (t.v[b - 1] < 0.5 && t.v[b] >= 0.5) { x50 = XLO + (b - 0.5 + (0.5 - t.v[b - 1]) / (t.v[b] - t.v[b - 1])) * BW; break; }
      Tmpl u; u.nbin = nbin; u.v.resize(nbin); for (int b = 0; b < nbin; ++b) { const double x = XLO + (b + 0.5) * BW + x50; u.v[b] = (x - XLO) / BW - 0.5 >= nbin - 1 ? 0 : t.at(x); }
      return u; }
    // first pulse: template fit on [lo, hi] around the time, scan +-5 ns around the 50 % crossing in 0.1 ns steps
    void fit1(const Wf& w, const Tmpl& T, double lo, double hi, double& tb, double& Ab) {
      auto chi = [&](double t, double& A) { double sty = 0, stt = 0; for (size_t q = 0; q < w.y.size(); ++q) { const double x = q * kNs - t; if (x < lo || x > hi) continue; const double m = T.at(x); sty += m * w.y[q]; stt += m * m; }
        A = stt > 0 ? sty / stt : 0; double c = 0; for (size_t q = 0; q < w.y.size(); ++q) { const double x = q * kNs - t; if (x < lo || x > hi) continue; const double r = w.y[q] - A * T.at(x); c += r * r; } return c; };
      double best = 1e30, A; tb = w.t50; for (double t = w.t50 - 5; t <= w.t50 + 5 + 1e-9; t += 0.1) { const double c = chi(t, A); if (c < best) { best = c; tb = t; } }
      chi(tb, Ab); }
  }

  class MC2MuCapture : public art::EDAnalyzer {
    public:
      struct Config {
        using Name = fhicl::Name; using Comment = fhicl::Comment;
        fhicl::Atom<art::InputTag> trkTag      { Name("caloTrackCollection"), Comment("MC2CaloCosmicTracks"), "MC2CaloCosmicTrackFinder" };
        fhicl::Atom<art::InputTag> digiTag     { Name("caloDigiCollection"), Comment("CaloDigis (all waveforms of the event)"), "CaloDigis" };
        fhicl::Atom<art::InputTag> hitTag      { Name("caloHitCollection"), Comment("CaloHits (mu+ veto)"), "MC2CaloHitMaker" };
        fhicl::Atom<art::InputTag> rdTag       { Name("caloRecoDigiCollection"), Comment("reco digis (delayed deposits in the other crystals)"), "MC2CaloRecoDigiMaker" };
        fhicl::Atom<double> nbNear        { Name("nbNear"), Comment("ring 1: max distance [mm]"), 80. };
        fhicl::Atom<double> nbFarMin      { Name("nbFarMin"), Comment("ring 2: min distance [mm]"), 250. };
        fhicl::Atom<double> nbFarMax      { Name("nbFarMax"), Comment("ring 2: max distance [mm]"), 400. };
        fhicl::Atom<double> nbDtMin       { Name("nbDtMin"), Comment("min dt [ns]"), -1000. };
        fhicl::Atom<double> nbDtMax       { Name("nbDtMax"), Comment("max dt [ns]"), 10000. };
        fhicl::Atom<double> vetoDist      { Name("vetoDist"), Comment("mu+ veto: max distance from the crystal [mm]"), 80. };
        fhicl::Atom<double> vetoMin       { Name("vetoMin"), Comment("mu+ veto: min delay [ns]"), 1000. };
        fhicl::Atom<double> vetoMax       { Name("vetoMax"), Comment("mu+ veto: max delay [ns]"), 20000. };
        fhicl::Atom<int>    minCrystals   { Name("minCrystals"), Comment("crystals of a one-disk track for a STOP"), 4 };
        fhicl::Atom<int>    minThrough    { Name("minThrough"), Comment("crystals of a track for a THROUGH"), 6 };
        fhicl::Atom<int>    maxThrough    { Name("maxThrough"), Comment("max THROUGH crystals kept"), 400000 };
        fhicl::Atom<int>    minSamples    { Name("minSamples"), Comment("min waveform length [samples]"), 16 };
        fhicl::Atom<double> minAmp        { Name("minAmp"), Comment("min first-pulse amplitude [ADC]"), 100. };
        fhicl::Atom<int>    satADC        { Name("satADC"), Comment("saturation [raw ADC]"), 4090 };
        fhicl::Atom<double> xMax          { Name("xMax"), Comment("template range [ns]"), 1500. };
        fhicl::Atom<int>    nIterTemplate { Name("nIterTemplate"), Comment("template iterations"), 2 };
        fhicl::Atom<double> dtMin         { Name("dtMin"), Comment("min delay of the second pulse [ns]"), 10. };
        fhicl::Atom<double> noiseADC      { Name("noiseADC"), Comment("noise per sample [ADC]"), 3.75 };
        fhicl::Atom<double> tailStart     { Name("tailStart"), Comment("start of the tail charge [ns after t1]"), 60. };
        fhicl::Atom<double> injMaxDelay   { Name("injMaxDelay"), Comment("injection: max delay [ns]"), 600. };
        fhicl::Atom<double> injMinA2      { Name("injMinA2"), Comment("injection: min A2 [ADC]"), 20. };
        fhicl::Atom<double> injMaxA2      { Name("injMaxA2"), Comment("injection: max A2 [ADC]"), 400. };
        fhicl::Atom<int>    seed          { Name("seed"), Comment("injection seed"), 12345 };
        fhicl::Atom<int>    nExamples     { Name("nExamples"), Comment("waveforms saved per class (tree examples): STOP found x5, THROUGH found, STOP not found"), 24 };
        fhicl::Atom<double> exampleS      { Name("exampleS"), Comment("found = maxS > exampleS and A2 > exampleA2"), 5. };
        fhicl::Atom<double> exampleA2     { Name("exampleA2"), Comment("[ADC]"), 30. };
        fhicl::Atom<double> secS          { Name("secS"), Comment("cluster: second pulse of the candidate crystal if maxS > secS"), 5. };
        fhicl::Atom<double> clusR         { Name("clusR"), Comment("cluster: max distance of the other crystals from the candidate [mm]"), 80. };
        fhicl::Atom<double> pairDt        { Name("pairDt"), Comment("cluster: max |tL - tR| of the two SiPMs of a crystal [ns]"), 10. };
        fhicl::Atom<double> clusGap       { Name("clusGap"), Comment("cluster: new cluster after a gap > clusGap [ns]"), 10. };
        fhicl::Atom<double> clusDtMax     { Name("clusDtMax"), Comment("cluster: max dt from the stop [ns]"), 20000. };
        fhicl::Atom<double> clusEmin      { Name("clusEmin"), Comment("cluster: min energy written to the tree [MeV]"), 1. };
        fhicl::Atom<double> eCut          { Name("eCut"), Comment("cluster: energy threshold of the physics sample [MeV] (summary counts)"), 15. };
      };
      using Parameters = art::EDAnalyzer::Table<Config>;
      explicit MC2MuCapture(const Parameters& c);
      void beginJob() override;
      void analyze(const art::Event& event) override;
      void endJob() override;
    private:
      bool waveform(const CaloDigi& d, Wf& w) const;
      void detect(const Cand& c, const Tmpl& T, double& amp, double& t1, double& maxS, double& tau, double& A2, double& Q, int& ns, double* t1s = nullptr, double* A1s = nullptr) const;
      void clusters(const Cand& c, double maxS, double tau, double A2);
      int cKind_, cCry_, cDisk_, cRun_, cSub_, cEvt_, cN_, cInWf_, cNMi_; float cDt_, cE_, cEstop_, cStopE_, cMaxS_; long nCut_[2] = {0, 0};
      art::InputTag trkTag_, digiTag_, hitTag_, rdTag_; double vetoD_, vetoLo_, vetoHi_, nbNear_, nbFarMin_, nbFarMax_, nbDtMin_, nbDtMax_; int minCry_, minThr_, maxThr_, minSamp_; double minAmp_; int sat_; double xMax_; int nIt_; double dtMin_, noise_, tail_, injMax_, injA2lo_, injA2hi_; int seed_; int nEx_; double exS_, exA2_;
      ProditionsHandle<CalCalib> calibH_; double secS_, clusR_, pairDt_, clusGap_, clusDtMax_, clusEmin_, eCut_; TTree* tc_ = nullptr;
      std::vector<Cand> cand_; long nThr_ = 0, nStopTried_ = 0, nStopNext_ = 0;
      TTree *t_ = nullptr, *ti_ = nullptr, *te_ = nullptr, *td_ = nullptr; TH1D* hT_ = nullptr;
      int dKind_, dCry_, dDisk_, dNMi_, dRing_, dNbCry_, dSiPM_; float dDt_, dE_, dDist_, dStopE_;
      int exCls_; int nExSaved_[3] = {0, 0, 0}; std::vector<float> exY0_, exY1_; float exT0a_, exT0b_, exA1a_, exA1b_, exT1a_, exT1b_, exTref_;
      int kind_, cry_, disk_, run_, sub_, evt_, ns_, nMi_; float amp_, t1_, maxS_, tau_, A2_, Q_, injD_, injA2_, miDt_, miE_;
  };

  MC2MuCapture::MC2MuCapture(const Parameters& c) : art::EDAnalyzer(c), trkTag_(c().trkTag()), digiTag_(c().digiTag()), hitTag_(c().hitTag()), rdTag_(c().rdTag()),
    vetoD_(c().vetoDist()), vetoLo_(c().vetoMin()), vetoHi_(c().vetoMax()), nbNear_(c().nbNear()), nbFarMin_(c().nbFarMin()), nbFarMax_(c().nbFarMax()), nbDtMin_(c().nbDtMin()), nbDtMax_(c().nbDtMax()), minCry_(c().minCrystals()),
    minThr_(c().minThrough()), maxThr_(c().maxThrough()), minSamp_(c().minSamples()), minAmp_(c().minAmp()), sat_(c().satADC()), xMax_(c().xMax()), nIt_(c().nIterTemplate()),
    dtMin_(c().dtMin()), noise_(c().noiseADC()), tail_(c().tailStart()), injMax_(c().injMaxDelay()), injA2lo_(c().injMinA2()), injA2hi_(c().injMaxA2()), seed_(c().seed()), nEx_(c().nExamples()), exS_(c().exampleS()), exA2_(c().exampleA2()),
    secS_(c().secS()), clusR_(c().clusR()), pairDt_(c().pairDt()), clusGap_(c().clusGap()), clusDtMax_(c().clusDtMax()), clusEmin_(c().clusEmin()), eCut_(c().eCut()) {
    consumes<MC2CaloCosmicTrackCollection>(trkTag_); consumes<CaloDigiCollection>(digiTag_); consumes<CaloHitCollection>(hitTag_); consumes<CaloRecoDigiCollection>(rdTag_);
  }

  void MC2MuCapture::beginJob() {
    art::ServiceHandle<art::TFileService> tfs;
    t_ = tfs->make<TTree>("capture", "second-pulse search in STOP (kind 0) and THROUGH (kind 1) crystals");
    ti_ = tfs->make<TTree>("inject", "THROUGH crystals with an injected second pulse");
    for (TTree* t : {t_, ti_}) { t->Branch("kind", &kind_); t->Branch("cry", &cry_); t->Branch("disk", &disk_); t->Branch("run", &run_); t->Branch("subrun", &sub_); t->Branch("event", &evt_);
      t->Branch("amp", &amp_); t->Branch("t1", &t1_); t->Branch("maxS", &maxS_); t->Branch("tau", &tau_); t->Branch("A2", &A2_); t->Branch("Q", &Q_); t->Branch("ns", &ns_);
      t->Branch("nMichel", &nMi_); t->Branch("michelDt", &miDt_); t->Branch("michelE", &miE_); }
    ti_->Branch("injDelay", &injD_); ti_->Branch("injA2", &injA2_);
    td_ = tfs->make<TTree>("delayed", "single-SiPM deposits in the other crystals of the disk, dt from the candidate (ring 1 near, 2 far)");
    td_->Branch("kind", &dKind_); td_->Branch("cry", &dCry_); td_->Branch("disk", &dDisk_); td_->Branch("nMichel", &dNMi_); td_->Branch("ring", &dRing_);
    td_->Branch("nbCry", &dNbCry_); td_->Branch("sipm", &dSiPM_); td_->Branch("dt", &dDt_); td_->Branch("e", &dE_); td_->Branch("dist", &dDist_); td_->Branch("stopE", &dStopE_);
    te_ = tfs->make<TTree>("examples", "example waveforms: cls 0 STOP found, 1 THROUGH found, 2 STOP not found");
    te_->Branch("cls", &exCls_); te_->Branch("kind", &kind_); te_->Branch("cry", &cry_); te_->Branch("disk", &disk_); te_->Branch("run", &run_); te_->Branch("subrun", &sub_); te_->Branch("event", &evt_);
    te_->Branch("amp", &amp_); te_->Branch("maxS", &maxS_); te_->Branch("tau", &tau_); te_->Branch("A2", &A2_); te_->Branch("Q", &Q_); te_->Branch("nMichel", &nMi_);
    te_->Branch("y0", &exY0_); te_->Branch("y1", &exY1_); te_->Branch("t0a", &exT0a_); te_->Branch("t0b", &exT0b_);   // raw - pedestal [ADC]; t0 of the digis [ns]
    te_->Branch("A1a", &exA1a_); te_->Branch("A1b", &exA1b_); te_->Branch("t1a", &exT1a_); te_->Branch("t1b", &exT1b_); te_->Branch("tref", &exTref_);
    tc_ = tfs->make<TTree>("cluster", "delayed clusters around STOP (kind 0) and THROUGH (kind 1) crystals: second pulse of the crystal + neighbours");
    tc_->Branch("kind", &cKind_); tc_->Branch("cry", &cCry_); tc_->Branch("disk", &cDisk_); tc_->Branch("run", &cRun_); tc_->Branch("subrun", &cSub_); tc_->Branch("event", &cEvt_);
    tc_->Branch("dt", &cDt_); tc_->Branch("E", &cE_); tc_->Branch("Estop", &cEstop_); tc_->Branch("n", &cN_); tc_->Branch("inWf", &cInWf_);
    tc_->Branch("stopE", &cStopE_); tc_->Branch("maxS", &cMaxS_); tc_->Branch("nMichel", &cNMi_);
    hT_ = tfs->make<TH1D>("htempl", "template from the THROUGH crystals of the run;time - 50% crossing [ns];normalised", int((xMax_ - XLO) / BW), XLO, xMax_);
  }

  bool MC2MuCapture::waveform(const CaloDigi& d, Wf& w) const {
    const auto& a = d.waveform(); if ((int)a.size() < minSamp_) return false;
    double p = 0; for (int q = 0; q < 4; ++q) p += a[q]; p /= 4.; int im = 0;
    for (size_t q = 0; q < a.size(); ++q) { if (a[q] >= sat_) return false; if (a[q] > a[im]) im = q; }
    const double amp = a[im] - p; if (amp < minAmp_) return false;
    w.y.resize(a.size()); for (size_t q = 0; q < a.size(); ++q) w.y[q] = float((a[q] - p) / amp);
    w.amp = amp; w.t0 = d.t0() * kNs; w.t50 = t50of(w.y); return true;
  }

  void MC2MuCapture::analyze(const art::Event& event) {
    const Calorimeter& cal = *GeomHandle<Calorimeter>();
    const auto& trks = *event.getValidHandle<MC2CaloCosmicTrackCollection>(trkTag_);
    if (trks.empty()) return;
    const auto& digis = *event.getValidHandle<CaloDigiCollection>(digiTag_);
    const auto& allHits = *event.getValidHandle<CaloHitCollection>(hitTag_);
    const auto& rdigis = *event.getValidHandle<CaloRecoDigiCollection>(rdTag_);
    const CalCalib& calib = calibH_.get(event.id());
    std::set<int> crysWithDigi; for (const auto& d : digis) crysWithDigi.insert(d.SiPMID() / CaloConst::_nSiPMPerCrystal);
    for (const auto& tk : trks) {
      if ((int)tk.hits.size() < std::min(minCry_, minThr_)) continue;
      double ux = tk.ux, uy = tk.uy; const double nn = std::hypot(ux, uy); if (nn <= 0) continue; ux /= nn; uy /= nn; if (uy > 0) { ux = -ux; uy = -uy; }   // downward
      std::vector<std::pair<double, size_t>> sk;
      for (size_t i = 0; i < tk.hits.size(); ++i) { const auto p = cal.mu2eToTracker(cal.crystal(tk.hits[i]->crystalID()).position()); sk.push_back({p.x() * ux + p.y() * uy, i}); }
      std::sort(sk.begin(), sk.end());
      auto addCand = [&](int kind, size_t ih) {
        const auto& h = tk.hits[ih]; Cand c; c.kind = kind; c.cry = h->crystalID(); c.disk = cal.crystal(c.cry).diskID(); c.run = event.run(); c.sub = event.subRun(); c.evt = event.event();
        int got = 0;
        for (const auto& rd : h->recoCaloDigis()) { const int s = rd->SiPMID() % CaloConst::_nSiPMPerCrystal; if (waveform(*rd->caloDigiPtr(), c.w[s])) ++got; }
        if (got != 2) return;
        for (const auto& rd : h->recoCaloDigis()) c.a2m[rd->SiPMID() % CaloConst::_nSiPMPerCrystal] = calib.ADC2MeV(rd->SiPMID());
        c.tEnd = std::min(c.w[0].t0 + kNs * (c.w[0].y.size() - 1), c.w[1].t0 + kNs * (c.w[1].y.size() - 1)) - h->time();   // approximate: digi t0 vs calibrated hit time
        // mu+ veto: later deposits near the crystal
        const auto pc = cal.mu2eToTracker(cal.crystal(c.cry).position());
        for (const auto& ah : allHits) {
          const double dt = ah.time() - h->time(); if (dt < vetoLo_ || dt > vetoHi_) continue;
          const auto& ac = cal.crystal(ah.crystalID()); if ((int)ac.diskID() != c.disk) continue;
          const auto pa = cal.mu2eToTracker(ac.position()); if (std::hypot(pa.x() - pc.x(), pa.y() - pc.y()) > vetoD_) continue;
          ++c.nMichel; c.michelE += ah.energyDep(); if (c.michelDt < 0 || dt < c.michelDt) c.michelDt = dt; }
        // delayed deposits in the other crystals (written here: no waveform needed)
        std::set<int> onTrk; for (const auto& hp : tk.hits) onTrk.insert(hp->crystalID());
        for (const auto& rd : rdigis) {
          const int nc = rd.SiPMID() / CaloConst::_nSiPMPerCrystal; if (nc == c.cry || onTrk.count(nc)) continue;
          const auto& ac = cal.crystal(nc); if ((int)ac.diskID() != c.disk) continue;
          const double dt = rd.time() - h->time(); if (dt < nbDtMin_ || dt > nbDtMax_) continue;
          const auto pa = cal.mu2eToTracker(ac.position()); const double dd = std::hypot(pa.x() - pc.x(), pa.y() - pc.y());
          const int ring = dd <= nbNear_ ? 1 : ((dd > nbFarMin_ && dd <= nbFarMax_) ? 2 : 0); if (!ring) continue;
          dKind_ = c.kind; dCry_ = c.cry; dDisk_ = c.disk; dNMi_ = c.nMichel; dRing_ = ring; dNbCry_ = nc; dSiPM_ = rd.SiPMID(); dDt_ = dt; dE_ = rd.energyDep(); dDist_ = dd; dStopE_ = h->energyDep();
          td_->Fill(); }
        {   // crystal deposits for the delayed cluster: both SiPMs within pairDt
          std::map<int, std::vector<const CaloRecoDigi*>> bySiPM[2];
          for (const auto& rd : rdigis) {
            const int nc = rd.SiPMID() / CaloConst::_nSiPMPerCrystal; if (nc == c.cry || onTrk.count(nc)) continue;
            const auto& ac = cal.crystal(nc); if ((int)ac.diskID() != c.disk) continue;
            const auto pa = cal.mu2eToTracker(ac.position()); if (std::hypot(pa.x() - pc.x(), pa.y() - pc.y()) > clusR_) continue;
            bySiPM[rd.SiPMID() % CaloConst::_nSiPMPerCrystal][nc].push_back(&rd); }
          for (const auto& [nc, v0] : bySiPM[0]) { auto it = bySiPM[1].find(nc); if (it == bySiPM[1].end()) continue;
            for (const auto* a : v0) { const CaloRecoDigi* b = nullptr; for (const auto* x : it->second) if (std::abs(x->time() - a->time()) < pairDt_ && (!b || std::abs(x->time() - a->time()) < std::abs(b->time() - a->time()))) b = x;
              if (!b) continue;
              const double dt = 0.5 * (a->time() + b->time()) - h->time();
              if (dt < dtMin_ || dt > clusDtMax_) continue;
              c.nb.push_back({float(dt), float(0.5 * (a->energyDep() + b->energyDep())), nc}); } } }
        c.stopE = h->energyDep();
        cand_.push_back(std::move(c)); };
      if ((tk.type == 1 || tk.type == 2) && (int)sk.size() >= minCry_) {
        ++nStopTried_;
        const size_t il = sk.back().second; const int cl = tk.hits[il]->crystalID(); const int dsk = cal.crystal(cl).diskID();
        const auto pl = cal.mu2eToTracker(cal.crystal(cl).position()); const double px = pl.x() + 34.3 * ux, py = pl.y() + 34.3 * uy;
        int cnext = -1; double db = 20.;
        for (unsigned k = 0; k < cal.nCrystals(); ++k) { if ((int)cal.crystal(k).diskID() != dsk) continue; const auto q = cal.mu2eToTracker(cal.crystal(k).position()); const double d = std::hypot(q.x() - px, q.y() - py); if (d < db) { db = d; cnext = k; } }
        bool onTrack = false; for (const auto& hp : tk.hits) if (hp->crystalID() == cnext) onTrack = true;
        if (cnext >= 0 && !onTrack) { if (crysWithDigi.count(cnext)) ++nStopNext_; else addCand(0, il); }
      }
      if ((int)sk.size() >= minThr_ && nThr_ < maxThr_) { addCand(1, sk[sk.size() / 2].second); ++nThr_; }
    }
  }

  void MC2MuCapture::detect(const Cand& c, const Tmpl& T, double& amp, double& t1o, double& maxS, double& tau, double& A2, double& Q, int& ns, double* t1s, double* A1s) const {
    // first pulse of each SiPM
    double t1[2], A1[2];
    for (int s = 0; s < 2; ++s) fit1(c.w[s], T, -20, 40, t1[s], A1[s]);
    if (t1s) for (int s = 0; s < 2; ++s) { t1s[s] = t1[s]; A1s[s] = A1[s] * c.w[s].amp; }
    // common grid: absolute times (t0 of the digi + sample), residual normalised to the first pulse, averaged over the SiPMs
    const double ta = c.w[0].t0 + t1[0], tb = c.w[1].t0 + t1[1], tref = 0.5 * (ta + tb);
    const double lo = std::max(c.w[0].t0, c.w[1].t0), hi = std::min(c.w[0].t0 + kNs * (c.w[0].y.size() - 1), c.w[1].t0 + kNs * (c.w[1].y.size() - 1));
    std::vector<double> tt, r; for (double t = lo; t <= hi + 1e-6; t += kNs) {
      double v = 0; for (int s = 0; s < 2; ++s) { const double qf = (t - c.w[s].t0) / kNs; const int q = int(std::lround(qf)); if (q < 0 || q >= (int)c.w[s].y.size()) { v = NAN; break; }
        const double tl = s == 0 ? ta : tb; v += (c.w[s].y[q] - A1[s] * T.at(t - tl)) / (A1[s] > 0 ? A1[s] : 1.); }
      if (std::isnan(v)) continue;
      tt.push_back(t - tref); r.push_back(0.5 * v); }
    amp = 0.5 * (A1[0] * c.w[0].amp + A1[1] * c.w[1].amp); t1o = tref; ns = r.size();
    const double sig = 0.5 * std::hypot(noise_ / std::max(1e-3, A1[0] * c.w[0].amp), noise_ / std::max(1e-3, A1[1] * c.w[1].amp));   // per sample, normalised
    maxS = -1e9; tau = -1; A2 = 0; Q = 0;
    for (size_t i = 0; i < tt.size(); ++i) if (tt[i] >= tail_) Q += r[i];
    if (tt.empty()) return;
    const double tEnd = tt.back() - 30.;
    for (double d = dtMin_; d <= tEnd; d += 2.5) {
      double sr = 0, s2 = 0; for (size_t i = 0; i < tt.size(); ++i) { const double m = T.at(tt[i] - d); if (m == 0) continue; sr += r[i] * m; s2 += m * m; }
      if (s2 <= 0) continue;
      const double S = sr / (sig * std::sqrt(s2)); if (S > maxS) { maxS = S; tau = d; A2 = sr / s2 * amp; } }
  }

  void MC2MuCapture::clusters(const Cand& c, double maxS, double tau, double A2) {
    std::vector<Dep> m = c.nb;
    double eStop2 = 0; if (maxS > secS_ && A2 > 0) { eStop2 = A2 * 0.5 * (c.a2m[0] + c.a2m[1]); m.push_back({float(tau), float(eStop2), c.cry}); }
    std::sort(m.begin(), m.end(), [](const Dep& a, const Dep& b) { return a.dt < b.dt; });
    for (size_t i = 0; i < m.size();) {
      size_t j = i + 1; while (j < m.size() && m[j].dt - m[j - 1].dt <= clusGap_) ++j;
      double e = 0, et = 0, es = 0; std::set<int> cr;
      for (size_t k = i; k < j; ++k) { e += m[k].e; et += m[k].e * m[k].dt; cr.insert(m[k].cry); if (m[k].cry == c.cry) es += m[k].e; }
      if (e >= clusEmin_ && e > 0) {
        cKind_ = c.kind; cCry_ = c.cry; cDisk_ = c.disk; cRun_ = c.run; cSub_ = c.sub; cEvt_ = c.evt; cDt_ = et / e; cE_ = e; cEstop_ = es; cN_ = cr.size();
        cInWf_ = cDt_ < c.tEnd - 30.; cStopE_ = c.stopE; cMaxS_ = maxS; cNMi_ = c.nMichel; tc_->Fill();
        if (e >= eCut_ && c.kind <= 1) ++nCut_[c.kind]; }
      i = j; }
  }

  void MC2MuCapture::endJob() {
    const int nbin = int((xMax_ - XLO) / BW);
    std::vector<const Wf*> wt; for (auto& c : cand_) if (c.kind == 1) for (int s = 0; s < 2; ++s) if (c.w[s].amp >= 300 && c.w[s].amp <= 1800) wt.push_back(&c.w[s]);
    printf("MC2MuCapture: candidates %zu (STOP tried %ld, rejected because the next crystal has a digi %ld; THROUGH %ld); template from %zu waveforms\n",
           cand_.size(), nStopTried_, nStopNext_, nThr_, wt.size());
    if (wt.size() < 200) { printf("MC2MuCapture: too few THROUGH waveforms for the template, nothing done\n"); return; }
    std::vector<double> ph; for (auto* w : wt) ph.push_back(std::fmod(w->t50, kNs)); std::vector<double> srt = ph; std::sort(srt.begin(), srt.end());
    std::vector<double> ta(wt.size()); for (size_t i = 0; i < wt.size(); ++i) ta[i] = wt[i]->t50 - ph[i] + kNs * double(std::lower_bound(srt.begin(), srt.end(), ph[i]) - srt.begin()) / srt.size();
    Tmpl T = build(wt, ta, nbin);
    for (int it = 0; it < nIt_; ++it) { for (size_t i = 0; i < wt.size(); ++i) { double tb, A; fit1(*wt[i], T, -20, 40, tb, A); ta[i] = tb; } T = build(wt, ta, nbin); }
    for (int b = 0; b < nbin; ++b) hT_->SetBinContent(b + 1, T.v[b]);
    TRandom3 rnd(seed_); long nStop = 0, nThr = 0;
    for (const auto& c : cand_) {
      double amp, t1, mS, tau, A2, Q; int ns; double t1s[2], A1s[2]; detect(c, T, amp, t1, mS, tau, A2, Q, ns, t1s, A1s);
      kind_ = c.kind; cry_ = c.cry; disk_ = c.disk; run_ = c.run; sub_ = c.sub; evt_ = c.evt; amp_ = amp; t1_ = t1; maxS_ = mS; tau_ = tau; A2_ = A2; Q_ = Q; ns_ = ns;
      nMi_ = c.nMichel; miDt_ = c.michelDt; miE_ = c.michelE; t_->Fill();
      clusters(c, mS, tau, A2);
      {   // examples
        const bool fnd = mS > exS_ && A2 > exA2_; const int cls = c.kind == 0 ? (fnd ? 0 : 2) : (fnd ? 1 : -1);
        const int lim = cls == 0 ? 5 * nEx_ : nEx_;
        if (cls >= 0 && nExSaved_[cls] < lim) { ++nExSaved_[cls]; exCls_ = cls;
          exY0_.assign(c.w[0].y.begin(), c.w[0].y.end()); for (auto& v : exY0_) v *= c.w[0].amp; exY1_.assign(c.w[1].y.begin(), c.w[1].y.end()); for (auto& v : exY1_) v *= c.w[1].amp;
          exT0a_ = c.w[0].t0; exT0b_ = c.w[1].t0; exA1a_ = A1s[0]; exA1b_ = A1s[1]; exT1a_ = t1s[0]; exT1b_ = t1s[1]; exTref_ = t1; te_->Fill(); } }
      (c.kind == 0 ? nStop : nThr)++;
      if (c.kind == 1) {   // injection
        Cand ci = c; const double d = rnd.Uniform(dtMin_, injMax_), a2 = rnd.Uniform(injA2lo_, injA2hi_);
        for (int s = 0; s < 2; ++s) { double tb, A; fit1(c.w[s], T, -20, 40, tb, A);
          for (size_t q = 0; q < ci.w[s].y.size(); ++q) ci.w[s].y[q] += float(a2 / c.w[s].amp * T.at(q * kNs - tb - d)); }
        detect(ci, T, amp, t1, mS, tau, A2, Q, ns);
        kind_ = 2; amp_ = amp; t1_ = t1; maxS_ = mS; tau_ = tau; A2_ = A2; Q_ = Q; ns_ = ns; injD_ = d; injA2_ = a2; ti_->Fill();
      }
    }
    printf("MC2MuCapture: STOP crystals analysed %ld, THROUGH %ld (each also with an injected pulse)\n", nStop, nThr);
    printf("MC2MuCapture: delayed clusters with E >= %.1f MeV: STOP %ld, THROUGH %ld\n", eCut_, nCut_[0], nCut_[1]);
  }
}
DEFINE_ART_MODULE(mu2e::MC2MuCapture)
