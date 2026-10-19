// MC2RecoCompare: performance of several reconstructions of the same calorimeter digis, measured on the same cosmic tracks.
//
// Every reconstruction (method) is a CaloRecoDigi collection made in the same job from the same CaloDigis (Offline
// CaloRecoDigiMaker, MC2CaloRecoDigiMaker with different templates or fit windows, ...). The tracks come from one reference
// reconstruction (MC2CaloCosmicTrackFinder): every method is measured on the same crystals, none chooses its own tracks.
// A method's time of a SiPM is the CaloRecoDigi of that method made from the same CaloDigi (largest energy if more than one).
//
// Common variables (from the raw waveform of the CaloDigi, the same for all methods):
//   amplitude = largest sample - pedestal (mean of the first nPed samples) [ADC]; saturated = largest sample >= satADC;
//   phase     = 50 % crossing of the rising edge (linear interpolation) modulo 5 ns: where the pulse falls in the sample.
// Crystals of a track are used if the reference has both SiPMs and neither is saturated.
//
// Protocol (the same for every method):
//   track line (as the mc2tools track check, docs lr_vs_tracks): two-disk tracks t = T - s/c (3D direction known, only T free);
//     one-disk tracks t = T + k s with T and k free (s = position along the direction in the disk plane: the inclination in z, hence
//     the slope, is unknown);
//   even events (track parity = event number % 2 == 0): per-SiPM time offsets of the method, global least squares on the tracks
//     t_kj = line_k(s_kj) + o_j (alternating least squares: line of every track, then offsets; mean offset 0, progressive outlier cut);
//     walkCorrection: t_kj = line_k(s_kj) + o_j + w(A_kj), w a table in ln of the SiPM raw amplitude A (walkBins logarithmic bins
//     from 150 to 2100 ADC, linear interpolation between bin centres, constant outside, mean 0 over the SiPMs), estimated in the same
//     alternating least squares; then subtracted from every SiPM time on the odd events too, before the line fit.
//   odd events: figures of merit with those offsets. The even events give the same figures, used only to choose the fit window
//     within a family (the method of a family with the smallest track residual on the even events is marked).
// Figures of merit:
//   track residual  crystal time T = (tL + tR)/2 (offsets subtracted) minus the line fitted to the crystals of the track (crystals
//                   farther than maxDevOthers from the line are left out of the fit, iteratively; at least minCrystals in the fit),
//                   studentised: divided by sqrt(1 - h) for a crystal in the fit, sqrt(1 + h) for one left out, h = its leverage.
//                   Its width is then the resolution of one crystal against the others, without the error of the line. It is the
//                   error that matters at Mu2e and it does not cancel what is common to L and R.
//                   Reported: half of the 16-84 % interval (main figure), sigma of a gaussian fit of the core, fraction > 1 ns.
//   L-R             ((tL - oL) - (tR - oR)) / sqrt2: resolution of one SiPM, blind to errors common to L and R.
//   phase bias      rms over 10 phase bins of the median track residual (and of L-R): a correct template gives a flat profile.
//   walk            median track residual in amplitude bins: its spread is an amplitude-dependent bias of the method.
//   efficiency      reference crystals in which the method has both SiPMs; chi2/ndf median on the track SiPMs (method's own);
//                   hChi2Amp: chi2/ndf versus the SiPM raw amplitude. chi2/ndf in log bins (100 per decade, 1e-2 to 1e3), amplitude
//                   in 20 ADC bins: millions of SiPMs, so the binning must not set the resolution of the figures.
// Output: summaryFile (text tables), histograms per method in the TFileService file, optional tree of the odd crystals.
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art_root_io/TFileService.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"
#include "fhiclcpp/types/Table.h"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/CalorimeterGeom/inc/Calorimeter.hh"
#include "Offline/CaloConditions/inc/CalDAQMap.hh"
#include "Offline/ProditionsService/inc/ProditionsHandle.hh"
#include "Offline/DataProducts/inc/CaloConst.hh"
#include "Offline/RecoDataProducts/inc/CaloRecoDigi.hh"
#include "Offline/RecoDataProducts/inc/CaloDigi.hh"
#include "Offline/MC2Reco/inc/MC2CaloCosmicTrack.hh"
#include "TF1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TTree.h"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>
#include <map>
#include <unordered_map>
#include <vector>

namespace mu2e {
  namespace {
    const double kC = 299.792458;   // mm/ns
    const float kNaN = std::numeric_limits<float>::quiet_NaN();
    const std::vector<double> kAmpBins = {150, 300, 500, 800, 1200, 1600, 4100};
    struct Cry { int track; int crystal; int sipm[2]; float s; float amp[2]; float phase[2]; bool sat; };
    struct Trk { int first, n; int parity; int type; };
    // line t = T + k s through points (s, t) with weights use (0/1); fixedK: k = -1/c (two-disk tracks), else k free.
    // Returns false if it cannot be fitted. h(s) = leverage of a point at s: 1/n + (s - sbar)^2 / Sxx (free slope) or 1/n.
    struct Line { double T = 0, k = 0, sbar = 0, sxx = 0; int n = 0; bool fixedK = false;
      double at(double s) const { return T + k * s; }
      double lev(double s) const { return n > 0 ? 1. / n + (fixedK || sxx <= 0 ? 0. : (s - sbar) * (s - sbar) / sxx) : 1.; } };
    bool fitLine(const std::vector<double>& s, const std::vector<double>& t, const std::vector<char>& use, bool fixedK, Line& L) {
      L = Line(); L.fixedK = fixedK; double ss = 0, st = 0; int n = 0;
      for (size_t i = 0; i < s.size(); ++i) if (use[i]) { ss += s[i]; st += t[i] + (fixedK ? s[i] / kC : 0.); ++n; }
      if (n < (fixedK ? 1 : 3)) return false;
      L.n = n; L.sbar = ss / n;
      if (fixedK) { L.k = -1. / kC; L.T = st / n; return true; }
      const double tbar = st / n; double sxx = 0, sxt = 0;
      for (size_t i = 0; i < s.size(); ++i) if (use[i]) { sxx += (s[i] - L.sbar) * (s[i] - L.sbar); sxt += (s[i] - L.sbar) * (t[i] - tbar); }
      if (sxx < 100.) return false;   // crystals spread over less than ~10 mm: no slope
      L.sxx = sxx; L.k = sxt / sxx; L.T = tbar - L.k * L.sbar; return true; }
    double quantile(std::vector<double> v, double q) {
      if (v.empty()) return kNaN;
      const size_t k = std::min(v.size() - 1, size_t(q * (v.size() - 1) + 0.5));
      std::nth_element(v.begin(), v.begin() + k, v.end()); return v[k]; }
    double half1684(const std::vector<double>& v) { return v.size() < 20 ? kNaN : 0.5 * (quantile(v, 0.84) - quantile(v, 0.16)); }
    double median(const std::vector<double>& v) { return quantile(v, 0.5); }
    // sigma of a gaussian fit of the core (+-1.5 half-width of 16-84 % around the median)
    double coreSigma(const std::vector<double>& v) {
      if (v.size() < 50) return kNaN;
      const double m = median(v), w = half1684(v); if (!(w > 0)) return kNaN;
      TH1F h("hcore", "", 200, m - 4 * w, m + 4 * w); h.SetDirectory(nullptr); for (double x : v) h.Fill(x);
      TF1 f("gcore", "gaus", m - 1.5 * w, m + 1.5 * w); f.SetParameters(h.GetMaximum(), m, w);
      if (h.Fit(&f, "QRN0") != 0) return kNaN;
      return std::abs(f.GetParameter(2)); }
    // rms over phase bins of the median of y (bins with >= 30 entries)
    double phaseBias(const std::vector<double>& y, const std::vector<double>& ph) {
      const int nb = 10; std::vector<std::vector<double>> b(nb);
      for (size_t i = 0; i < y.size(); ++i) { int k = int(ph[i] / 5. * nb); if (k >= 0 && k < nb) b[k].push_back(y[i]); }
      std::vector<double> m; for (auto& v : b) if (v.size() >= 30) m.push_back(median(v));
      if (m.size() < 3) return kNaN;
      double a = 0; for (double x : m) a += x; a /= m.size();
      double s = 0; for (double x : m) s += (x - a) * (x - a); return std::sqrt(s / m.size()); }
  }

  class MC2RecoCompare : public art::EDAnalyzer {
    public:
      struct MethodConfig {
        using Name = fhicl::Name; using Comment = fhicl::Comment;
        fhicl::Atom<std::string>   name   { Name("name"), Comment("name of the method in the tables") };
        fhicl::Atom<std::string>   family { Name("family"), Comment("methods of one family differ only in the fit window: the best on the even events is marked") };
        fhicl::Atom<art::InputTag> tag    { Name("caloRecoDigiCollection"), Comment("CaloRecoDigi collection of the method") };
      };
      struct Config {
        using Name = fhicl::Name; using Comment = fhicl::Comment;
        fhicl::Atom<art::InputTag> trkTag        { Name("caloTrackCollection"), Comment("reference MC2CaloCosmicTracks"), "MC2CaloCosmicTrackFinder" };
        fhicl::Sequence<fhicl::Table<MethodConfig>> methods { Name("methods"), Comment("reconstructions to compare") };
        fhicl::Atom<double>        depthInCrystal{ Name("depthInCrystal"), Comment("crystal centre behind the front face [mm]"), 100. };
        fhicl::Sequence<double>    diskZOffset   { Name("diskZOffset"), Comment("added to the Offline z of disk 0 and 1 [mm]"), std::vector<double>{0., 0.} };
        fhicl::Atom<int>           nPed          { Name("nPed"), Comment("samples for the pedestal of the raw amplitude"), 4 };
        fhicl::Atom<int>           satADC        { Name("satADC"), Comment("largest raw sample >= satADC: saturated"), 4090 };
        fhicl::Atom<int>           nIterTime     { Name("nIterTime"), Comment("iterations of the offset least squares"), 25 };
        fhicl::Atom<double>        outlierCut    { Name("outlierCut"), Comment("final max |residual| in the least squares [ns]"), 5. };
        fhicl::Atom<bool>          walkCorr      { Name("walkCorrection"), Comment("fit a walk w(A) of the SiPM raw amplitude with the offsets and subtract it"), false };
        fhicl::Atom<int>           walkBins      { Name("walkBins"), Comment("logarithmic amplitude bins of the walk table, 150-2100 ADC"), 24 };
        fhicl::Atom<int>           minEntriesTime{ Name("minEntriesTime"), Comment("min entries per SiPM for an offset"), 10 };
        fhicl::Atom<int>           minCrystals   { Name("minCrystals"), Comment("min usable crystals of a track (tested one included)"), 4 };
        fhicl::Atom<double>        maxDevOthers  { Name("maxDevOthers"), Comment("other crystals farther than this from the track median are not used [ns]"), 5. };
        fhicl::Atom<std::string>   summaryFile   { Name("summaryFile"), Comment("text tables"), "mc2_reco_compare.txt" };
        fhicl::Atom<bool>          writeTree     { Name("writeTree"), Comment("tree of the odd crystals, every method"), false };
      };
      using Parameters = art::EDAnalyzer::Table<Config>;
      explicit MC2RecoCompare(const Parameters& conf);
      void beginJob() override;
      void analyze(const art::Event& event) override;
      void endJob() override;

    private:
      struct Method { std::string name, family; art::InputTag tag; };
      struct Fom { double eff = kNaN, chi2 = kNaN, resHalf = kNaN, resCore = kNaN, resTail = kNaN, lrHalf = kNaN, lrCore = kNaN,
                   phRes = kNaN, phLR = kNaN, walk = kNaN; long n = 0; std::vector<double> walkBins, resBins; };
      art::InputTag trkTag_; std::vector<Method> m_; double depth_; std::vector<double> dzOff_; int nPed_, sat_, nIt_; double cut_;
      int minT_, minCry_; double maxDev_; std::string sumFile_; bool writeTree_; bool walkCorr_; int nW_;
      std::vector<double> w_;   // walk table of the current method (empty: none)
      double wAt(double a) const;   // walk at the SiPM raw amplitude a, from w_
      ProditionsHandle<CalDAQMap> mapH_;
      std::vector<Cry> cry_; std::vector<Trk> trk_; std::vector<float> t_;   // t_[(c * M + m) * 2 + side]: time of method m, NaN if none
      std::vector<TH1F*> hChi2_; std::vector<TH2F*> hChi2Amp_;
      long nEvents_ = 0;
      void fom(int m, const std::vector<double>& o, const std::vector<char>& cal, int parity, Fom& f, bool fill);
      std::vector<TH1F*> hRes_, hLR_; std::vector<TH2F*> hResPhase_, hResAmp_, hLRPhase_;
      TTree* tree_ = nullptr; int trM_ = 0, trCry_ = 0, trType_ = 0, trN_ = 0; float trRes_ = 0, trLR_ = 0, trAmp_ = 0, trPh_ = 0;
  };

  MC2RecoCompare::MC2RecoCompare(const Parameters& conf) : art::EDAnalyzer(conf),
    trkTag_(conf().trkTag()), depth_(conf().depthInCrystal()), dzOff_(conf().diskZOffset()), nPed_(conf().nPed()), sat_(conf().satADC()),
    nIt_(conf().nIterTime()), cut_(conf().outlierCut()), minT_(conf().minEntriesTime()), minCry_(conf().minCrystals()),
    maxDev_(conf().maxDevOthers()), sumFile_(conf().summaryFile()), writeTree_(conf().writeTree()),
    walkCorr_(conf().walkCorr()), nW_(conf().walkBins()) {
    if (dzOff_.size() != 2) throw cet::exception("MC2RECOCOMPARE") << "diskZOffset needs two values\n";
    for (const auto& mc : conf().methods()) m_.push_back({mc.name(), mc.family(), mc.tag()});
    if (m_.empty()) throw cet::exception("MC2RECOCOMPARE") << "no methods\n";
    consumes<MC2CaloCosmicTrackCollection>(trkTag_);
    for (const auto& m : m_) consumes<CaloRecoDigiCollection>(m.tag);
  }

  void MC2RecoCompare::beginJob() {
    art::ServiceHandle<art::TFileService> tfs;
    std::vector<double> lc(501); for (int i = 0; i <= 500; ++i) lc[i] = std::pow(10., -2. + i / 100.);   // chi2/ndf log bins
    for (const auto& m : m_) {
      art::TFileDirectory d = tfs->mkdir(m.name);
      hChi2_.push_back(d.make<TH1F>("hChi2", (m.name + ": track SiPMs;#chi^{2}/ndf;SiPMs").c_str(), 500, lc.data()));
      hChi2Amp_.push_back(d.make<TH2F>("hChi2Amp", (m.name + ": track SiPMs;SiPM raw amplitude [ADC];#chi^{2}/ndf").c_str(), 205, 0, 4100, 500, lc.data()));
      hRes_.push_back(d.make<TH1F>("hRes", (m.name + ": odd events;t_{crystal} - t_{other crystals} [ns];crystals").c_str(), 800, -4, 4));
      hLR_.push_back(d.make<TH1F>("hLR", (m.name + ": odd events;(t_{L} - t_{R})/#sqrt{2} [ns];crystals").c_str(), 800, -4, 4));
      hResPhase_.push_back(d.make<TH2F>("hResPhase", (m.name + ": odd events;phase of L [ns];track residual [ns]").c_str(), 50, 0, 5, 200, -2, 2));
      hLRPhase_.push_back(d.make<TH2F>("hLRPhase", (m.name + ": odd events;phase of L [ns];(t_{L} - t_{R})/#sqrt{2} [ns]").c_str(), 50, 0, 5, 200, -2, 2));
      hResAmp_.push_back(d.make<TH2F>("hResAmp", (m.name + ": odd events;crystal amplitude (mean of L and R, raw) [ADC];track residual [ns]").c_str(), 82, 0, 4100, 200, -2, 2));
    }
    if (writeTree_) {
      tree_ = tfs->make<TTree>("crystals", "odd-event track crystals, one entry per method");
      tree_->Branch("method", &trM_); tree_->Branch("crystal", &trCry_); tree_->Branch("res", &trRes_); tree_->Branch("lr", &trLR_);
      tree_->Branch("amp", &trAmp_); tree_->Branch("phase", &trPh_); tree_->Branch("trackType", &trType_); tree_->Branch("nFit", &trN_);
    }
  }

  void MC2RecoCompare::analyze(const art::Event& event) {
    ++nEvents_;
    const Calorimeter& cal = *GeomHandle<Calorimeter>();
    const int M = m_.size();
    // per method: CaloDigi key -> CaloRecoDigi (largest energy)
    std::vector<std::unordered_map<size_t, const CaloRecoDigi*>> byDigi(M);
    for (int m = 0; m < M; ++m)
      for (const auto& rd : *event.getValidHandle<CaloRecoDigiCollection>(m_[m].tag)) {
        auto& p = byDigi[m][rd.caloDigiPtr().key()]; if (!p || rd.energyDep() > p->energyDep()) p = &rd; }
    const int parity = event.id().event() % 2;
    for (const auto& tk : *event.getValidHandle<MC2CaloCosmicTrackCollection>(trkTag_)) {
      Trk tr{int(cry_.size()), 0, parity, tk.type}; const double inPlane = std::hypot(tk.ux, tk.uy);
      for (const auto& hp : tk.hits) {
        if (hp->recoCaloDigis().size() != 2) continue;
        const auto& c = cal.crystal(hp->crystalID()); const int d = c.diskID(); const auto p = cal.mu2eToTracker(c.position());
        const double cz = p.z() + depth_ + dzOff_[d];
        double s;   // position along the line, up
        if (tk.type == 3) s = (p.x() - tk.cx) * tk.ux + (p.y() - tk.cy) * tk.uy + (cz - tk.cz) * tk.uz;
        else s = ((p.x() - tk.cx) * tk.ux + (p.y() - tk.cy) * tk.uy) / (inPlane > 0 ? inPlane : 1.);
        Cry cr{int(trk_.size()), hp->crystalID(), {-1, -1}, float(s), {0, 0}, {0, 0}, false};
        std::vector<float> tt(2 * M, kNaN); bool okSides = true;
        for (const auto& rp : hp->recoCaloDigis()) {
          const auto& dg = rp->caloDigiPtr(); const int side = rp->SiPMID() % 2; if (cr.sipm[side] >= 0) { okSides = false; break; }
          cr.sipm[side] = rp->SiPMID();
          const auto& w = dg->waveform(); if ((int)w.size() <= nPed_ + 1) { okSides = false; break; }
          double ped = 0; for (int q = 0; q < nPed_; ++q) ped += w[q]; ped /= nPed_;
          const int im = std::max_element(w.begin(), w.end()) - w.begin(); const double a = w[im] - ped;
          if (w[im] >= sat_) cr.sat = true;
          cr.amp[side] = a;
          double ph = kNaN; for (int q = im; q > 0; --q) if (w[q - 1] - ped < 0.5 * a) {
            const double y0 = w[q - 1] - ped, y1 = w[q] - ped; const double t50 = 5. * (q - 1 + (0.5 * a - y0) / (y1 - y0)); ph = std::fmod(t50, 5.); break; }
          cr.phase[side] = ph;
          for (int m = 0; m < M; ++m) { auto it = byDigi[m].find(dg.key()); if (it == byDigi[m].end()) continue;
            tt[2 * m + side] = it->second->time();
            if (it->second->ndf() > 0) { const double c2 = it->second->chi2() / it->second->ndf(); hChi2_[m]->Fill(c2); hChi2Amp_[m]->Fill(a, c2); } }
        }
        if (!okSides || cr.sipm[0] < 0 || cr.sipm[1] < 0) continue;
        cry_.push_back(cr); t_.insert(t_.end(), tt.begin(), tt.end()); ++tr.n;
      }
      if (tr.n > 0) trk_.push_back(tr);
    }
  }

  double MC2RecoCompare::wAt(double a) const {
    if (w_.empty() || !(a > 0)) return 0.;
    const double x = std::log(a / 150.) / std::log(2100. / 150.) * nW_ - 0.5;   // in units of bins, 0 = centre of the first
    if (x <= 0) return w_.front();
    if (x >= nW_ - 1) return w_.back();
    const int i = int(x); return w_[i] + (x - i) * (w_[i + 1] - w_[i]);
  }

  // figures of merit of method m on the tracks of one parity, with offsets o (cal[j]: SiPM j calibrated)
  void MC2RecoCompare::fom(int m, const std::vector<double>& o, const std::vector<char>& cal, int parity, Fom& f, bool fill) {
    const int M = m_.size(); const double sq2 = std::sqrt(2.);
    std::vector<double> res, lr, ph, amp, phRes; long nRef = 0, nOk = 0;
    std::vector<std::vector<double>> resB(kAmpBins.size() - 1);
    std::vector<double> tc, sc; std::vector<int> ic;
    for (const auto& tr : trk_) {
      if (tr.parity != parity) continue;
      tc.clear(); sc.clear(); ic.clear();
      for (int i = tr.first; i < tr.first + tr.n; ++i) {
        const Cry& c = cry_[i]; if (c.sat) continue; ++nRef;
        const float tL = t_[(size_t(i) * M + m) * 2], tR = t_[(size_t(i) * M + m) * 2 + 1];
        if (std::isnan(tL) || std::isnan(tR)) continue;
        ++nOk;
        if (!cal[c.sipm[0]] || !cal[c.sipm[1]]) continue;
        const double a = tL - o[c.sipm[0]] - wAt(c.amp[0]), b = tR - o[c.sipm[1]] - wAt(c.amp[1]);
        const double d = (a - b) / sq2; lr.push_back(d); ph.push_back(c.phase[0]);
        if (fill) { hLR_[m]->Fill(d); hLRPhase_[m]->Fill(c.phase[0], d); }
        tc.push_back(0.5 * (a + b)); sc.push_back(c.s); ic.push_back(i);
      }
      if ((int)tc.size() < minCry_) continue;
      // line through the crystals, outliers (> maxDevOthers) left out iteratively
      std::vector<char> use(tc.size(), 1); Line L; bool ok = fitLine(sc, tc, use, tr.type == 3, L);
      for (int it = 0; ok && it < 3; ++it) { bool ch = false;
        for (size_t q = 0; q < tc.size(); ++q) { const char u = std::abs(tc[q] - L.at(sc[q])) < maxDev_; if (u != use[q]) { use[q] = u; ch = true; } }
        if (!ch) break;
        ok = fitLine(sc, tc, use, tr.type == 3, L); }
      if (!ok || L.n < minCry_) continue;
      for (size_t k = 0; k < tc.size(); ++k) {
        const double h = L.lev(sc[k]); const double den = use[k] ? 1. - h : 1. + h; if (den <= 0.05) continue;
        const double r = (tc[k] - L.at(sc[k])) / std::sqrt(den); const Cry& c = cry_[ic[k]]; const double am = 0.5 * (c.amp[0] + c.amp[1]);
        res.push_back(r); amp.push_back(am); phRes.push_back(c.phase[0]);
        for (size_t b = 0; b + 1 < kAmpBins.size(); ++b) if (am >= kAmpBins[b] && am < kAmpBins[b + 1]) resB[b].push_back(r);
        if (fill) { hRes_[m]->Fill(r); hResPhase_[m]->Fill(c.phase[0], r); hResAmp_[m]->Fill(am, r);
          if (tree_) { trM_ = m; trCry_ = c.crystal; trRes_ = r; trAmp_ = am; trPh_ = c.phase[0]; trType_ = tr.type; trN_ = L.n;
            trLR_ = ((t_[(size_t(ic[k]) * M + m) * 2] - o[c.sipm[0]] - wAt(c.amp[0])) - (t_[(size_t(ic[k]) * M + m) * 2 + 1] - o[c.sipm[1]] - wAt(c.amp[1]))) / sq2; tree_->Fill(); } }
      }
    }
    f.n = res.size(); f.eff = nRef ? double(nOk) / nRef : kNaN;
    f.resHalf = half1684(res); f.resCore = coreSigma(res);
    long nt = 0; for (double r : res) if (std::abs(r) > 1.) ++nt; f.resTail = res.empty() ? kNaN : double(nt) / res.size();
    f.lrHalf = half1684(lr); f.lrCore = coreSigma(lr); f.phLR = phaseBias(lr, ph);
    f.walkBins.clear(); f.resBins.clear(); double wmin = 1e9, wmax = -1e9;
    for (auto& v : resB) { const bool ok = v.size() >= 50; const double md = ok ? median(v) : kNaN; f.walkBins.push_back(md); f.resBins.push_back(ok ? half1684(v) : kNaN);
      if (ok) { wmin = std::min(wmin, md); wmax = std::max(wmax, md); } }
    f.walk = wmax > wmin ? wmax - wmin : kNaN;
    f.phRes = phaseBias(res, phRes);
    const TH1F* h = hChi2_[m]; if (h->GetEntries() > 0) { double q = 0.5; double x; const_cast<TH1F*>(h)->GetQuantiles(1, &x, &q); f.chi2 = x; }
  }

  void MC2RecoCompare::endJob() {
    const int M = m_.size(), NS = CaloConst::_nChannel;
    std::vector<Fom> even(M), odd(M); std::vector<int> nCal(M, 0); std::vector<std::vector<double>> walkTab(M);
    for (int m = 0; m < M; ++m) {
      // ---- offsets of method m: alternating least squares on the even tracks (SiPM times; both SiPMs present) ----
      struct E { int track; int sipm; double s, t, a; bool use; };
      std::vector<E> e;
      for (size_t k = 0; k < trk_.size(); ++k) { const Trk& tr = trk_[k]; if (tr.parity != 0) continue;
        for (int i = tr.first; i < tr.first + tr.n; ++i) { const Cry& c = cry_[i]; if (c.sat) continue;
          for (int sd = 0; sd < 2; ++sd) { const float t = t_[(size_t(i) * M + m) * 2 + sd]; if (!std::isnan(t)) e.push_back({int(k), c.sipm[sd], c.s, t, c.amp[sd], true}); } } }
      std::vector<std::vector<size_t>> byTrk(trk_.size()); for (size_t i = 0; i < e.size(); ++i) byTrk[e[i].track].push_back(i);
      std::vector<double> o(NS, 0.), sO(NS); std::vector<int> nO(NS); std::vector<Line> line(trk_.size()); std::vector<char> lok(trk_.size(), 0);
      std::vector<double> ls, lt; std::vector<char> lu; w_.clear();
      auto trackLines = [&]() {   // line of every track from its SiPM times minus the current offsets
        for (size_t k = 0; k < trk_.size(); ++k) { if (byTrk[k].empty()) { lok[k] = 0; continue; }
          ls.clear(); lt.clear(); lu.clear();
          for (size_t i : byTrk[k]) { ls.push_back(e[i].s); lt.push_back(e[i].t - o[e[i].sipm] - wAt(e[i].a)); lu.push_back(e[i].use); }
          lok[k] = fitLine(ls, lt, lu, trk_[k].type == 3, line[k]); } };
      trackLines();
      for (int it = 0; it < nIt_; ++it) {
        if (it >= 3 && it % 3 == 0) {   // progressive outlier cut: 50, 17, 5.6 ns, ..., then outlierCut
          const double cnow = std::max(cut_, 50. / std::pow(3., (it - 3) / 3));
          for (auto& x : e) x.use = lok[x.track] && std::abs(x.t - line[x.track].at(x.s) - o[x.sipm] - wAt(x.a)) < cnow; }
        std::fill(sO.begin(), sO.end(), 0.); std::fill(nO.begin(), nO.end(), 0);
        for (auto& x : e) if (x.use && lok[x.track] && line[x.track].n > 1) { sO[x.sipm] += x.t - line[x.track].at(x.s) - wAt(x.a); ++nO[x.sipm]; }
        double mm = 0; int nm = 0; for (int j = 0; j < NS; ++j) if (nO[j] >= minT_) { o[j] = sO[j] / nO[j]; mm += o[j]; ++nm; } else o[j] = 0;
        if (nm) { mm /= nm; for (int j = 0; j < NS; ++j) if (nO[j] >= minT_) o[j] -= mm; }
        if (walkCorr_) {   // walk: mean of t - line - o per amplitude bin (bins with >= 50 entries, else the previous value), mean 0
          std::vector<double> sw(nW_, 0.); std::vector<long> nw(nW_, 0); double tot = 0; long ntot = 0;
          for (auto& x : e) { if (!(x.use && lok[x.track] && line[x.track].n > 1) || !(x.a > 0)) continue;
            const int b = std::clamp(int(std::log(x.a / 150.) / std::log(2100. / 150.) * nW_), 0, nW_ - 1);
            sw[b] += x.t - line[x.track].at(x.s) - o[x.sipm]; ++nw[b]; }
          std::vector<double> w = w_.empty() ? std::vector<double>(nW_, 0.) : w_;
          for (int b = 0; b < nW_; ++b) if (nw[b] >= 50) { w[b] = sw[b] / nw[b]; tot += w[b] * nw[b]; ntot += nw[b]; }
          if (ntot) for (auto& v : w) v -= tot / ntot;
          w_ = w;
        }
        trackLines();
      }
      std::vector<char> cal(NS, 0); for (int j = 0; j < NS; ++j) if (nO[j] >= minT_) { cal[j] = 1; ++nCal[m]; }
      fom(m, o, cal, 0, even[m], false);
      fom(m, o, cal, 1, odd[m], true);
      walkTab[m] = w_;
    }
    // ---- best of each family on the even events ----
    std::map<std::string, int> best;
    for (int m = 0; m < M; ++m) { auto it = best.find(m_[m].family);
      if (!std::isnan(even[m].resHalf) && (it == best.end() || even[m].resHalf < even[it->second].resHalf)) best[m_[m].family] = m; }
    // ---- tables ----
    FILE* f = fopen(sumFile_.c_str(), "w"); if (!f) throw cet::exception("MC2RECOCOMPARE") << "cannot write " << sumFile_ << "\n";
    long nTrk[2] = {0, 0}, nCry[2] = {0, 0}; for (const auto& tr : trk_) { ++nTrk[tr.parity]; nCry[tr.parity] += tr.n; }
    fprintf(f, "MC2RecoCompare: %ld events; reference tracks %ld even / %ld odd, crystals with both SiPMs %ld / %ld\n", nEvents_, nTrk[0], nTrk[1], nCry[0], nCry[1]);
    fprintf(f, "offsets per method from the even tracks (global least squares), figures of merit on the odd tracks; * = best of its family on the even tracks\n");
    fprintf(f, "res = crystal time - line of the track (T - s/c two disks, T + k s one disk), studentised [ns]: half of 16-84 %% (main figure), gaussian core, fraction > 1 ns\n");
    fprintf(f, "LR = (tL - tR)/sqrt2 [ns]; phase = rms over the sampling phase of the median residual [ns]; walk = spread of the median residual over the amplitude bins [ns]\n\n");
    fprintf(f, "%-28s %-14s | %5s %6s %6s | %8s | %8s %8s %6s | %7s %7s | %7s %7s | %6s | %8s\n", "method", "family", "eff", "chi2", "nCalib",
            "res even", "res", "core", ">1ns", "LR", "LR core", "ph res", "ph LR", "walk", "crystals");
    for (int m = 0; m < M; ++m) { const Fom& a = odd[m]; const bool b = best.count(m_[m].family) && best[m_[m].family] == m;
      fprintf(f, "%-28s %-14s | %5.3f %6.2f %6d | %8.4f | %7.4f%s %8.4f %6.3f | %7.4f %7.4f | %7.4f %7.4f | %6.3f | %8ld\n", m_[m].name.c_str(), m_[m].family.c_str(),
              a.eff, a.chi2, nCal[m], even[m].resHalf, a.resHalf, b ? "*" : " ", a.resCore, a.resTail, a.lrHalf, a.lrCore, a.phRes, a.phLR, a.walk, a.n); }
    fprintf(f, "\nper amplitude bin (crystal amplitude = mean of L and R raw amplitudes [ADC]), odd tracks: half 16-84 %% of the residual / median residual [ns]\n");
    fprintf(f, "%-28s", "method"); for (size_t b = 0; b + 1 < kAmpBins.size(); ++b) fprintf(f, " | %4.0f-%-4.0f      ", kAmpBins[b], kAmpBins[b + 1]); fprintf(f, "\n");
    for (int m = 0; m < M; ++m) { fprintf(f, "%-28s", m_[m].name.c_str());
      for (size_t b = 0; b < odd[m].resBins.size(); ++b) fprintf(f, " | %6.3f / %+6.3f", odd[m].resBins[b], odd[m].walkBins[b]);
      fprintf(f, "\n"); }
    if (walkCorr_) {   // walk tables from the even tracks, subtracted from the SiPM times
      fprintf(f, "\nwalk correction w(A) [ns] from the even tracks, %d logarithmic bins of the SiPM raw amplitude from 150 to 2100 ADC (bin centres):\n%-28s", nW_, "method");
      for (int b = 0; b < nW_; ++b) fprintf(f, " %6.0f", 150. * std::pow(2100. / 150., (b + 0.5) / nW_));
      fprintf(f, "\n");
      for (int m = 0; m < M; ++m) { fprintf(f, "%-28s", m_[m].name.c_str()); for (double v : walkTab[m]) fprintf(f, " %+6.3f", v); fprintf(f, "\n"); }
    }
    fclose(f);
    printf("MC2RecoCompare: %zu methods, %ld + %ld tracks -> %s\n", m_.size(), nTrk[0], nTrk[1], sumFile_.c_str());
  }
}
DEFINE_ART_MODULE(mu2e::MC2RecoCompare)
