// MC2CaloCalibFromTracks: our calibration machine. Per-SiPM time offsets and energy scale of the calorimeter from the cosmic tracks
// produced in the same job (MC2CaloCosmicTrack, MC2CaloCosmicTrackFinder), written as Offline condition tables (CalTimeCalib,
// CalEnergyCalib), the format of the calorimeter group's tables, so that the two can be compared line by line.
//
// Starting point: the CalCalib conditions the reconstruction used (e.g. the tables of the calorimeter group loaded with
// DbService.textFile); the hit times are already corrected with those offsets (MC2CaloRecoDigiMaker useCalCalibTimeOffset) and the
// SiPM amplitude is A = CaloRecoDigi energy / CalCalib ADC2MeV.
//
// Time (global least squares, as mc2tools time_calibration lsq): for every track k and every SiPM j of its crystals,
//   t_kj = T_k - s_kj / c + o_j,  s_kj = position of the crystal centre along the line (up), c = speed of light;
//   T_k and o_j by alternating least squares (nIterTime), with the mean of the o_j fixed to 0 (the absolute alignment of the
//   starting tables is kept); from iteration nIterBeforeCut, every 5 iterations, entries with |residual| > 200, 67, 22, 7.4 ns, then
//   outlierCut, are left out (progressive: offsets of hundreds of ns, as in a pass without tables, can converge first).
//   Before it, a board step robust to jumps (see endJob). One-disk tracks: s along the line in the disk plane (3D direction unknown).
//   New offset = starting CalCalib offset + o_j (SiPMs with fewer than minEntriesTime entries: starting offset, flagged in stats).
// Energy (as mc2tools energy_calibration.py): chord of the line in the 34 x 34 mm cross-section of the crystal (x-y); two-disk
//   tracks: chord / sqrt(1 - uz^2); chords < minChordFrac * 34 mm left out; A34 = A * dp(34 mm) / dp(path), dp = Landau most
//   probable loss in CsI (PDG); per SiPM, Landau (x) Gauss (ROOT tutorial langaus) on the one-disk tracks, histogram 120 bins over
//   [0, 3 median], fit range [0.5, 2.5] x peak, >= minEntriesEnergy entries, accepted if error < 10 % of the MPV;
//   ADC2MeV = dp(34 mm) / MPV (SiPMs without a good fit: median ADC2MeV of the fitted SiPMs, flagged).
// Optional (kalSeedCollection not empty): tracker lines vs the calo hits of the tracks, t_SiPM(corrected) - t_line, histogram only
//   (calorimeter-tracker alignment check).
// Detector response for the simulation (all pulses of the run, not only the tracks):
//   baseline, noise: mean and rms of the first nPed samples of every CaloDigi waveform;
//   threshold: per SiPM, spectrum of the pulse amplitude (A of our fit, and raw = largest sample - baseline), rising edge fitted
//     with 0.5 (1 + erf((x - T) / (sqrt2 w))) N exp(-(x - T) / lambda); T = threshold [ADC], w = width of the edge (noise on the
//     threshold); in MeV with the new ADC2MeV of the SiPM. Fallback: half-maximum crossing of the rising edge (w = -1).
//   effective photoelectrons: per crystal, from the crystals of the tracks with both SiPMs, amplitudes brought to the mean gain of
//     the two, p0 = (Var(L - R) - noise variance) / <L + R> [ADC per effective photoelectron] (robust variance); pePerMeV of each
//     SiPM = ADCPerMeV / p0. Crystals without enough entries: median p0 (flagged). Var(L - R) also contains the left/right sharing
//     with the transverse position: p0 is an upper estimate.
// Output: timeTableFile, energyTableFile (TABLE CalTimeCalib / CalEnergyCalib, roid,value), simTableFile (TABLE CalSimCrystals:
//   crid, LRU (defaultLRU, not measured here), pePerMeV0, pePerMeV1, ADCPerMeV0, ADCPerMeV1), baselineTableFile (TABLE
//   CalBaselines: roid, baseline, threshold [raw ADC]), statsFile and responseFile (csv), TFileService histograms.
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art_root_io/TFileService.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/CalorimeterGeom/inc/Calorimeter.hh"
#include "Offline/CaloConditions/inc/CalCalib.hh"
#include "Offline/CaloConditions/inc/CalDAQMap.hh"
#include "Offline/ProditionsService/inc/ProditionsHandle.hh"
#include "Offline/DataProducts/inc/CaloConst.hh"
#include "Offline/RecoDataProducts/inc/KalSeed.hh"
#include "Offline/RecoDataProducts/inc/CaloRecoDigi.hh"
#include "Offline/RecoDataProducts/inc/CaloDigi.hh"
#include "KinKal/Trajectory/KinematicLine.hh"
#include "Offline/MC2Reco/inc/MC2CaloCosmicTrack.hh"
#include "TF1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <vector>

namespace mu2e {
  namespace {
    const double kC = 299.792458;   // mm/ns
    // Landau most probable energy loss in CsI for a path x (mm), beta = 1 (PDG; same constants as mc2tools energy_calibration.py)
    double dp(double x) { const double K = 0.307075, ZA = 0.41569, RHO = 4.510, HWP = 39.46e-6, ME = 0.51099895, J = 0.200;
      const double xi = 0.5 * K * ZA * RHO * x / 10.0; return xi * (std::log(2 * ME * xi / (HWP * HWP)) + J); }
    double langaus(double* x, double* par) {   // ROOT tutorial langaus: par = width (Landau), MPV, area, sigma (Gauss)
      const double invsq2pi = 0.3989422804014, mpshift = -0.22278298, np = 100.0, sc = 5.0;
      const double mpc = par[1] - mpshift * par[0], xlow = x[0] - sc * par[3], xupp = x[0] + sc * par[3], step = (xupp - xlow) / np; double sum = 0;
      for (double i = 1.0; i <= np / 2; i++) {
        double xx = xlow + (i - .5) * step; double fl = TMath::Landau(xx, mpc, par[0]) / par[0]; sum += fl * TMath::Gaus(x[0], xx, par[3]);
        xx = xupp - (i - .5) * step; fl = TMath::Landau(xx, mpc, par[0]) / par[0]; sum += fl * TMath::Gaus(x[0], xx, par[3]); }
      return par[2] * step * sum * invsq2pi / par[3];
    }
    // chord of the 2D line (point p, unit direction d) in the square of side a centred at (cx, cy); 0 if it misses
    double chord(double px, double py, double dx, double dy, double cx, double cy, double a) {
      double t0 = -1e30, t1 = 1e30; const double lo[2] = {cx - a / 2, cy - a / 2}, hi[2] = {cx + a / 2, cy + a / 2}, p[2] = {px, py}, d[2] = {dx, dy};
      for (int q = 0; q < 2; ++q) {
        if (std::abs(d[q]) < 1e-12) { if (p[q] < lo[q] || p[q] > hi[q]) return 0; continue; }
        double a1 = (lo[q] - p[q]) / d[q], a2 = (hi[q] - p[q]) / d[q]; if (a1 > a2) std::swap(a1, a2); t0 = std::max(t0, a1); t1 = std::min(t1, a2); }
      return t1 > t0 ? t1 - t0 : 0;
    }
    struct TE { int track, sipm, board; double s, t; bool use; };
  }

  class MC2CaloCalibFromTracks : public art::EDAnalyzer {
    public:
      struct Config {
        using Name = fhicl::Name; using Comment = fhicl::Comment;
        fhicl::Atom<art::InputTag> trkTag          { Name("caloTrackCollection"), Comment("MC2CaloCosmicTracks"), "MC2CaloCosmicTrackFinder" };
        fhicl::Atom<std::string>   kalSeedTag      { Name("kalSeedCollection"), Comment("tracker KalSeeds for the alignment check; empty = none"), "" };
        fhicl::Atom<double>        depthInCrystal  { Name("depthInCrystal"), Comment("crystal centre behind the front face [mm]"), 100. };
        fhicl::Sequence<double>    diskZOffset     { Name("diskZOffset"), Comment("added to the Offline z of disk 0 and 1 [mm]"), std::vector<double>{0., 0.} };
        fhicl::Atom<int>           nIterTime       { Name("nIterTime"), Comment("iterations of the time least squares"), 40 };
        fhicl::Atom<int>           nIterBeforeCut  { Name("nIterBeforeCut"), Comment("iterations before the outlier cut"), 10 };
        fhicl::Atom<double>        outlierCut      { Name("outlierCut"), Comment("max |residual| after the first iterations [ns]"), 5. };
        fhicl::Atom<int>           minEntriesTime  { Name("minEntriesTime"), Comment("min entries per SiPM for a time offset"), 10 };
        fhicl::Atom<double>        crystalSide     { Name("crystalSide"), Comment("crystal side [mm]"), 34. };
        fhicl::Atom<double>        minChordFrac    { Name("minChordFrac"), Comment("min chord / side"), 0.6 };
        fhicl::Atom<int>           minEntriesEnergy{ Name("minEntriesEnergy"), Comment("min entries per SiPM for the MIP fit"), 80 };
        fhicl::Atom<double>        rMatchTracker   { Name("rMatchTracker"), Comment("crystal centre - tracker line [mm]"), 60. };
        fhicl::Atom<std::string>   timeTableFile   { Name("timeTableFile"), Comment("output CalTimeCalib table"), "CalTimeCalib_mc2.txt" };
        fhicl::Atom<std::string>   energyTableFile { Name("energyTableFile"), Comment("output CalEnergyCalib table"), "CalEnergyCalib_mc2.txt" };
        fhicl::Atom<std::string>   statsFile       { Name("statsFile"), Comment("per-SiPM statistics (csv)"), "calib_mc2_stats.csv" };
        fhicl::Atom<art::InputTag> recoDigiTag     { Name("caloRecoDigiCollection"), Comment("all reco digis of the event (thresholds, baselines)"), "MC2CaloRecoDigiMaker" };
        fhicl::Atom<int>           nPed            { Name("nPed"), Comment("samples for the baseline"), 4 };
        fhicl::Atom<double>        thrMaxADC       { Name("thrMaxADC"), Comment("range of the threshold spectra [ADC]"), 400. };
        fhicl::Atom<int>           minPulsesThr    { Name("minPulsesThr"), Comment("min pulses per SiPM for the threshold"), 200 };
        fhicl::Atom<int>           minEntriesPE    { Name("minEntriesPE"), Comment("min track crystals for p0"), 30 };
        fhicl::Atom<double>        defaultLRU      { Name("defaultLRU"), Comment("LRU written in CalSimCrystals (not measured here)"), 0.05 };
        fhicl::Atom<std::string>   simTableFile    { Name("simTableFile"), Comment("output CalSimCrystals table"), "CalSimCrystals_mc2.txt" };
        fhicl::Atom<std::string>   baselineTableFile { Name("baselineTableFile"), Comment("output CalBaselines table"), "CalBaselines_mc2.txt" };
        fhicl::Atom<std::string>   responseFile    { Name("responseFile"), Comment("per-SiPM response for the simulation (csv)"), "response_mc2.csv" };
      };
      using Parameters = art::EDAnalyzer::Table<Config>;
      explicit MC2CaloCalibFromTracks(const Parameters& conf);
      void beginJob() override;
      void analyze(const art::Event& event) override;
      void endJob() override;
    private:
      art::InputTag trkTag_; std::string kalTag_; double depth_; std::vector<double> dzOff_; int nIt_, nItCut_; double cut_; int minT_;
      double side_, minChord_; int minE_; double rMatch_; std::string tFile_, eFile_, sFile_;
      art::InputTag rdTag_; int nPed_; double thrMax_; int minThr_, minPE_; double lru_; std::string simFile_, baseFile_, respFile_;
      std::vector<double> bSum_, b2Sum_; std::vector<long> bN_; struct LR { int crid; float a0, a1; }; std::vector<LR> lr_;
      TH2F *hAmpFit_ = nullptr, *hAmpRaw_ = nullptr; TH1F *hP0_ = nullptr, *hThr_ = nullptr, *hThrMeV_ = nullptr;
      ProditionsHandle<CalCalib> calibH_; ProditionsHandle<CalDAQMap> mapH_; std::vector<int> board_;
      std::vector<TE> te_; int ntrk_ = 0; std::vector<std::vector<float>> a34_;
      std::vector<double> prevT_, prevE_; bool havePrev_ = false;
      TH1F *hTres_ = nullptr, *hTresBefore_ = nullptr, *hOff_ = nullptr, *hMPV_ = nullptr, *hChord_ = nullptr, *hTrkDt_ = nullptr; TH2F* hTrkDtSiPM_ = nullptr;
  };

  MC2CaloCalibFromTracks::MC2CaloCalibFromTracks(const Parameters& conf) : art::EDAnalyzer(conf),
    trkTag_(conf().trkTag()), kalTag_(conf().kalSeedTag()), depth_(conf().depthInCrystal()), dzOff_(conf().diskZOffset()), nIt_(conf().nIterTime()),
    nItCut_(conf().nIterBeforeCut()), cut_(conf().outlierCut()), minT_(conf().minEntriesTime()), side_(conf().crystalSide()), minChord_(conf().minChordFrac()),
    minE_(conf().minEntriesEnergy()), rMatch_(conf().rMatchTracker()), tFile_(conf().timeTableFile()), eFile_(conf().energyTableFile()), sFile_(conf().statsFile()),
    rdTag_(conf().recoDigiTag()), nPed_(conf().nPed()), thrMax_(conf().thrMaxADC()), minThr_(conf().minPulsesThr()), minPE_(conf().minEntriesPE()),
    lru_(conf().defaultLRU()), simFile_(conf().simTableFile()), baseFile_(conf().baselineTableFile()), respFile_(conf().responseFile()),
    bSum_(CaloConst::_nChannel, 0.), b2Sum_(CaloConst::_nChannel, 0.), bN_(CaloConst::_nChannel, 0), a34_(CaloConst::_nChannel) {
    if (dzOff_.size() != 2) throw cet::exception("MC2CALOCALIB") << "diskZOffset needs two values\n";
    consumes<MC2CaloCosmicTrackCollection>(trkTag_); consumes<CaloRecoDigiCollection>(rdTag_); if (!kalTag_.empty()) consumes<KalSeedCollection>(art::InputTag(kalTag_));
  }

  void MC2CaloCalibFromTracks::beginJob() {
    art::ServiceHandle<art::TFileService> tfs;
    hTresBefore_ = tfs->make<TH1F>("hTresBefore", "before the least squares (starting offsets);t_{SiPM} - track time [ns];SiPM entries", 400, -20, 20);
    hTres_ = tfs->make<TH1F>("hTres", "after the least squares;t_{SiPM} - track time [ns];SiPM entries", 400, -20, 20);
    hOff_ = tfs->make<TH1F>("hOff", ";correction to the starting time offset [ns];SiPMs", 400, -20, 20);
    hMPV_ = tfs->make<TH1F>("hMPV", ";MIP MPV per SiPM [ADC per 34 mm];SiPMs", 200, 0, 2000);
    hChord_ = tfs->make<TH1F>("hChord", ";path in the crystal [mm];crystals", 120, 0, 120);
    hAmpFit_ = tfs->make<TH2F>("hAmpFit", "all pulses;offline SiPM id;amplitude of the fit [ADC]", CaloConst::_nChannel, 0, CaloConst::_nChannel, int(thrMax_ / 2), 0, thrMax_);
    hAmpRaw_ = tfs->make<TH2F>("hAmpRaw", "all pulses;offline SiPM id;largest sample - baseline [ADC]", CaloConst::_nChannel, 0, CaloConst::_nChannel, int(thrMax_ / 2), 0, thrMax_);
    hThr_ = tfs->make<TH1F>("hThr", ";threshold (raw) [ADC];SiPMs", 200, 0, thrMax_);
    hThrMeV_ = tfs->make<TH1F>("hThrMeV", ";threshold [MeV];SiPMs", 200, 0, 20);
    hP0_ = tfs->make<TH1F>("hP0", ";p_{0} [ADC per effective photoelectron];crystals", 200, 0, 5);
    if (!kalTag_.empty()) {
      hTrkDt_ = tfs->make<TH1F>("hTrkDt", "calo hits of the calo tracks vs the tracker line;t_{SiPM} - t_{line} [ns];SiPM entries", 1000, -3000, 3000);
      hTrkDtSiPM_ = tfs->make<TH2F>("hTrkDtSiPM", ";offline SiPM id;t_{SiPM} - t_{line} [ns]", CaloConst::_nChannel, 0, CaloConst::_nChannel, 600, -3000, 3000);
    }
  }

  void MC2CaloCalibFromTracks::analyze(const art::Event& event) {
    const Calorimeter& cal = *GeomHandle<Calorimeter>(); const CalCalib& calib = calibH_.get(event.id());
    if (board_.empty()) { const CalDAQMap& dm = mapH_.get(event.id()); board_.resize(CaloConst::_nChannel); for (int i = 0; i < (int)CaloConst::_nChannel; ++i) board_[i] = dm.rawId(CaloSiPMId(i)).dirac(); }
    if (!havePrev_) { prevT_.resize(CaloConst::_nChannel); prevE_.resize(CaloConst::_nChannel);
      for (int i = 0; i < (int)CaloConst::_nChannel; ++i) { prevT_[i] = calib.timeoffset(i); prevE_[i] = calib.ADC2MeV(i); } havePrev_ = true; }
    for (const auto& rd : *event.getValidHandle<CaloRecoDigiCollection>(rdTag_)) {
      const int id = rd.SiPMID(); if (id < 0 || id >= (int)CaloConst::_nChannel) continue;
      const auto& w = rd.caloDigiPtr()->waveform(); if ((int)w.size() <= nPed_) continue;
      double b = 0; for (int q = 0; q < nPed_; ++q) b += w[q]; b /= nPed_; double v = 0; for (int q = 0; q < nPed_; ++q) v += (w[q] - b) * (w[q] - b);
      bSum_[id] += b; b2Sum_[id] += v / nPed_; ++bN_[id];
      const double mx = *std::max_element(w.begin(), w.end());
      hAmpRaw_->Fill(id + 0.5, mx - b); if (prevE_[id] > 0) hAmpFit_->Fill(id + 0.5, rd.energyDep() / prevE_[id]);
    }
    const auto& trks = *event.getValidHandle<MC2CaloCosmicTrackCollection>(trkTag_);
    const KalSeedCollection* seeds = nullptr;
    if (!kalTag_.empty()) { art::Handle<KalSeedCollection> h; if (event.getByLabel(art::InputTag(kalTag_), h)) seeds = h.product(); }
    const double dp34 = dp(side_);
    for (const auto& tk : trks) {
      const int k = ntrk_++; const double inPlane = std::hypot(tk.ux, tk.uy);
      for (const auto& hp : tk.hits) {
        const auto& c = cal.crystal(hp->crystalID()); const int d = c.diskID(); const auto p = cal.mu2eToTracker(c.position());
        const double cz = p.z() + depth_ + dzOff_[d];
        double s;   // position along the line, up
        if (tk.type == 3) s = (p.x() - tk.cx) * tk.ux + (p.y() - tk.cy) * tk.uy + (cz - tk.cz) * tk.uz;
        else s = ((p.x() - tk.cx) * tk.ux + (p.y() - tk.cy) * tk.uy) / (inPlane > 0 ? inPlane : 1.);
        // path in the crystal
        double path = 0;
        if (inPlane > 0) { const double ch = chord(tk.cx, tk.cy, tk.ux / inPlane, tk.uy / inPlane, p.x(), p.y(), side_);
          path = (tk.type == 3) ? ch / std::sqrt(std::max(1e-6, 1. - double(tk.uz) * tk.uz)) : ch; if (ch < minChord_ * side_) path = 0; else hChord_->Fill(path); }
        if (hp->recoCaloDigis().size() == 2) {   // photostatistics: both SiPMs of the crystal
          float a[2] = {-1, -1};
          for (const auto& rd : hp->recoCaloDigis()) { const int id = rd->SiPMID(); if (id >= 0 && id < (int)CaloConst::_nChannel && prevE_[id] > 0) a[id % 2] = rd->energyDep() / prevE_[id]; }
          if (a[0] > 0 && a[1] > 0) lr_.push_back({hp->crystalID(), a[0], a[1]});
        }
        for (const auto& rd : hp->recoCaloDigis()) {
          const int id = rd->SiPMID(); if (id < 0 || id >= (int)CaloConst::_nChannel) continue;
          te_.push_back({k, id, board_[id], s, rd->time(), true});
          if (path > 0 && tk.type < 3 && prevE_[id] > 0) a34_[id].push_back(rd->energyDep() / prevE_[id] * dp34 / dp(path));
        }
        // tracker check
        if (seeds) for (const auto& ks : *seeds) {
          if (!ks.kinematicLineFit() || ks.segments().empty()) continue;
          const KalSegment* best = nullptr; double zb = -1e30;
          for (const auto& sg : ks.segments()) { const double z = sg.position3().Z(); if (z > zb) { zb = z; best = &sg; } }
          const auto line = best->kinematicLine(); const double tref = 0.5 * (best->tmin() + best->tmax()); const auto p0 = line.position3(tref); const auto v = line.velocity(tref);
          const double vm2 = v.Mag2(); if (vm2 <= 0) continue;
          const double dx = p.x() - p0.X(), dy = p.y() - p0.Y(), dz = cz - p0.Z(); const double u = (dx * v.X() + dy * v.Y() + dz * v.Z()) / vm2;
          const double ex = dx - v.X() * u, ey = dy - v.Y() * u, ez = dz - v.Z() * u; if (std::sqrt(ex * ex + ey * ey + ez * ez) > rMatch_) continue;
          for (const auto& rd : hp->recoCaloDigis()) { const double r = rd->time() - (tref + u); hTrkDt_->Fill(r); hTrkDtSiPM_->Fill(rd->SiPMID() + 0.5, r); }
        }
      }
    }
  }

  void MC2CaloCalibFromTracks::endJob() {
    const int NS = CaloConst::_nChannel;
    // ---- time, step 1: board offsets robust to jumps (e.g. after a re-initialisation): track time = median of its SiPMs (insensitive
    // to a minority of jumped boards), per board the mode of the residuals over +-3 us (10 ns sliding window) then the mean within
    // +-10 ns; 3 iterations; the mean board correction is fixed to 0
    const int NB = CaloConst::_nDIRAC; std::vector<double> ob(NB, 0.); std::vector<int> jumpN(NB, 0);
    {
      std::vector<std::vector<size_t>> byTrk(ntrk_); for (size_t i = 0; i < te_.size(); ++i) byTrk[te_[i].track].push_back(i);
      std::vector<double> Tm(ntrk_, 0.);
      for (int itb = 0; itb < 3; ++itb) {
        for (int k = 0; k < ntrk_; ++k) { std::vector<double> v; for (size_t i : byTrk[k]) v.push_back(te_[i].t - ob[te_[i].board] + te_[i].s / kC);
          if (v.empty()) continue;
          std::nth_element(v.begin(), v.begin() + v.size() / 2, v.end()); Tm[k] = v[v.size() / 2]; }
        std::vector<std::vector<int>> h(NB, std::vector<int>(3000, 0));   // 2 ns bins over [-3000, 3000] ns
        for (auto& e : te_) { if (byTrk[e.track].size() < 4) continue; const double r = e.t - ob[e.board] - (Tm[e.track] - e.s / kC); const int b = int((r + 3000.) / 2.); if (b >= 0 && b < 3000) h[e.board][b]++; }
        double msum = 0; int mn = 0;
        for (int bd = 0; bd < NB; ++bd) {
          int best = -1, bw = 0; for (int b = 0; b + 4 < 3000; ++b) { int w = 0; for (int q = 0; q < 5; ++q) w += h[bd][b + q]; if (w > bw) { bw = w; best = b + 2; } }
          if (bw < 10) continue;
          double c = -3000. + 2. * best + 1., sw = 0, sx = 0;
          for (int b = std::max(0, best - 5); b <= std::min(2999, best + 5); ++b) { const double x = -3000. + 2. * b + 1.; sw += h[bd][b]; sx += h[bd][b] * x; }
          if (sw > 0) c = sx / sw;
          ob[bd] += c; if (std::abs(c) > 20.) ++jumpN[bd];
          msum += ob[bd]; ++mn; }
        if (mn) for (int bd = 0; bd < NB; ++bd) ob[bd] -= msum / mn;
      }
      int nj = 0; for (int bd = 0; bd < NB; ++bd) if (std::abs(ob[bd]) > 20.) { ++nj; printf("MC2CaloCalibFromTracks: board %d jumped by %.1f ns with respect to the start tables\n", bd, ob[bd]); }
      printf("MC2CaloCalibFromTracks: board step: %d boards with a correction > 20 ns\n", nj);
    }
    // ---- time, step 2: alternating least squares per SiPM, starting from the board offsets ----
    std::vector<double> o(NS, 0.), T(ntrk_, 0.); std::vector<int> nT(ntrk_), nO(NS);
    for (int j = 0; j < NS; ++j) o[j] = ob[board_.empty() ? 0 : board_[j]];
    std::vector<double> sT(ntrk_), sO(NS);
    auto trackTimes = [&]() { std::fill(sT.begin(), sT.end(), 0.); std::fill(nT.begin(), nT.end(), 0);
      for (auto& e : te_) if (e.use) { sT[e.track] += e.t - o[e.sipm] + e.s / kC; ++nT[e.track]; }
      for (int k = 0; k < ntrk_; ++k) T[k] = nT[k] ? sT[k] / nT[k] : 0; };
    trackTimes(); for (auto& e : te_) hTresBefore_->Fill(e.t - (T[e.track] - e.s / kC));
    for (int it = 0; it < nIt_; ++it) {
      if (it >= nItCut_ && (it - nItCut_) % 5 == 0) {   // progressive outlier cut: 200, 67, 22, 7.4 ns, then outlierCut
        const double cnow = std::max(cut_, 200. / std::pow(3., (it - nItCut_) / 5));
        for (auto& e : te_) if (std::abs(e.t - (T[e.track] - e.s / kC) - o[e.sipm]) > cnow) e.use = false; }
      std::fill(sO.begin(), sO.end(), 0.); std::fill(nO.begin(), nO.end(), 0);
      for (auto& e : te_) if (e.use && nT[e.track] > 1) { sO[e.sipm] += e.t - (T[e.track] - e.s / kC); ++nO[e.sipm]; }
      double m = 0; int nm = 0; for (int j = 0; j < NS; ++j) if (nO[j] >= minT_) { o[j] = sO[j] / nO[j]; m += o[j]; ++nm; } else o[j] = ob[board_[j]];   // few entries: board offset
      if (nm) { m /= nm; for (int j = 0; j < NS; ++j) if (nO[j] >= minT_) o[j] -= m; }
      trackTimes();
    }
    for (auto& e : te_) if (e.use) hTres_->Fill(e.t - (T[e.track] - e.s / kC) - o[e.sipm]);
    // ---- energy: MIP fit per SiPM ----
    const double dp34 = dp(side_); std::vector<double> mpv(NS, -1), mpvErr(NS, -1), chi2(NS, -1); std::vector<int> nE(NS, 0);
    for (int j = 0; j < NS; ++j) {
      auto& v = a34_[j]; nE[j] = v.size(); if ((int)v.size() < minE_) continue;
      std::vector<float> w(v); std::nth_element(w.begin(), w.begin() + w.size() / 2, w.end()); const double med = w[w.size() / 2]; if (med <= 0) continue;
      TH1F h("hmip", "", 120, 0, 3 * med); h.SetDirectory(nullptr); for (float x : v) h.Fill(x);
      const double pk = h.GetBinCenter(h.GetMaximumBin());
      TF1 f("mc2lg", langaus, 0.5 * pk, 2.5 * pk, 4); f.SetParameters(0.08 * pk, pk, h.Integral("width"), 0.08 * pk);
      f.SetParLimits(0, 0.01 * pk, 0.5 * pk); f.SetParLimits(1, 0.5 * pk, 1.8 * pk); f.SetParLimits(3, 0.005 * pk, 0.6 * pk);
      h.Fit(&f, "QRN0");
      if (f.GetParameter(1) > 0 && f.GetParError(1) > 0 && f.GetParError(1) < 0.1 * f.GetParameter(1)) {
        mpv[j] = f.GetParameter(1); mpvErr[j] = f.GetParError(1); chi2[j] = f.GetNDF() > 0 ? f.GetChisquare() / f.GetNDF() : -1; hMPV_->Fill(mpv[j]); }
    }
    // new ADC2MeV per SiPM (used below for MeV and pePerMeV)
    // SiPMs without a good MIP fit: median ADC2MeV of the fitted SiPMs (as mc2tools energy_calibration.py)
    std::vector<double> fitted; for (int j = 0; j < NS; ++j) if (mpv[j] > 0) fitted.push_back(dp34 / mpv[j]);
    double medE = 0.; if (!fitted.empty()) { std::nth_element(fitted.begin(), fitted.begin() + fitted.size() / 2, fitted.end()); medE = fitted[fitted.size() / 2]; }
    std::vector<double> newE(NS); for (int j = 0; j < NS; ++j) newE[j] = mpv[j] > 0 ? dp34 / mpv[j] : (medE > 0 ? medE : (havePrev_ ? prevE_[j] : 0.));
    // ---- thresholds ----
    auto edge = [&](TH2F* h2, int j, double& T, double& w) {
      T = -1; w = -1; TH1D* h = h2->ProjectionY("pj", j + 1, j + 1); if (h->Integral() < minThr_) { delete h; return; }
      const int im = h->GetMaximumBin(); const double m = h->GetBinContent(im); int k = 1; while (k < im && h->GetBinContent(k) < 0.5 * m) ++k;
      const double c1 = h->GetBinContent(std::max(1, k - 1)), c2 = h->GetBinContent(k), x1 = h->GetBinCenter(std::max(1, k - 1)), x2 = h->GetBinCenter(k);
      const double th = (c2 > c1) ? x1 + (0.5 * m - c1) / (c2 - c1) * (x2 - x1) : x2; T = th;
      TF1 f("edge", "[0]*0.5*(1+TMath::Erf((x-[1])/(sqrt(2.)*[2])))*exp(-(x-[1])/[3])", std::max(0., th - 30), std::min(thrMax_, th + 80));
      f.SetParameters(m, th, 3., 40.); f.SetParLimits(2, 0.2, 30.); f.SetParLimits(3, 2., 1e4);
      if (h->Fit(&f, "QRN0") == 0 && std::abs(f.GetParameter(1) - th) < 20) { T = f.GetParameter(1); w = f.GetParameter(2); }
      delete h; };
    std::vector<double> thrRaw(NS, -1), thrRawW(NS, -1), thrFit(NS, -1), thrFitW(NS, -1);
    for (int j = 0; j < NS; ++j) { edge(hAmpRaw_, j, thrRaw[j], thrRawW[j]); edge(hAmpFit_, j, thrFit[j], thrFitW[j]);
      if (thrRaw[j] > 0) hThr_->Fill(thrRaw[j]);
      if (thrFit[j] > 0 && newE[j] > 0) hThrMeV_->Fill(thrFit[j] * newE[j]); }
    // ---- effective photoelectrons per crystal ----
    const int NC = CaloConst::_nCrystal; std::vector<std::vector<std::pair<float, float>>> byC(NC);
    for (auto& e : lr_) if (e.crid >= 0 && e.crid < NC) byC[e.crid].push_back({e.a0, e.a1});
    std::vector<double> p0(NC, -1), p0e(NC, -1); std::vector<double> allp;
    for (int c = 0; c < NC; ++c) {
      auto& v = byC[c]; if ((int)v.size() < minPE_) continue; const int j0 = 2 * c, j1 = 2 * c + 1; if (newE[j0] <= 0 || newE[j1] <= 0) continue;
      const double g0 = 1. / newE[j0], g1 = 1. / newE[j1], g = 0.5 * (g0 + g1);
      std::vector<double> d, sm; for (auto& x : v) { const double a0 = x.first * g / g0, a1 = x.second * g / g1; d.push_back(a0 - a1); sm.push_back(a0 + a1); }
      std::vector<double> dd(d); std::nth_element(dd.begin(), dd.begin() + dd.size() / 2, dd.end()); const double md = dd[dd.size() / 2];
      std::vector<double> ad; for (double x : d) ad.push_back(std::abs(x - md)); std::nth_element(ad.begin(), ad.begin() + ad.size() / 2, ad.end());
      const double var = std::pow(1.4826 * ad[ad.size() / 2], 2); double ms = 0; for (double x : sm) ms += x; ms /= sm.size();
      const double vn = (bN_[j0] ? b2Sum_[j0] / bN_[j0] : 0) + (bN_[j1] ? b2Sum_[j1] / bN_[j1] : 0);
      if (ms > 0 && var > vn) { p0[c] = (var - vn) / ms; p0e[c] = p0[c] * std::sqrt(2. / v.size()); allp.push_back(p0[c]); hP0_->Fill(p0[c]); }
    }
    double p0med = -1; if (!allp.empty()) { std::nth_element(allp.begin(), allp.begin() + allp.size() / 2, allp.end()); p0med = allp[allp.size() / 2]; }
    // ---- tables ----
    FILE* ft = fopen(tFile_.c_str(), "w"); FILE* fe = fopen(eFile_.c_str(), "w"); FILE* fs = fopen(sFile_.c_str(), "w");
    if (!ft || !fe || !fs) throw cet::exception("MC2CALOCALIB") << "cannot write the output tables\n";
    fprintf(ft, "TABLE CalTimeCalib\n#roid,timeOffset  MC2CaloCalibFromTracks: starting tables + least squares on %d calorimeter cosmic tracks\n", ntrk_);
    fprintf(fe, "TABLE CalEnergyCalib\n#roid,ADC2MeV  MC2CaloCalibFromTracks: MIP MPV per 34 mm (Landau x Gauss), ADC2MeV = %.4f MeV / MPV\n", dp34);
    fprintf(fs, "roid,nTime,timeCorrection_ns,startTimeOffset_ns,newTimeOffset_ns,timeFromTracks,nEnergy,mpv_adc,mpv_err,chi2ndf,startADC2MeV,newADC2MeV,energyFromTracks\n");
    int nTok = 0, nEok = 0;
    for (int j = 0; j < NS; ++j) {
      const bool tok = nO[j] >= minT_, eok = mpv[j] > 0;
      const double corr = tok ? o[j] : ob[board_[j]];   // SiPMs without enough entries: correction of their board (jumps)
      const double nt = havePrev_ ? prevT_[j] + corr : corr, pe = havePrev_ ? prevE_[j] : 0, ne = newE[j];
      if (tok) { hOff_->Fill(o[j]); ++nTok; } if (eok) ++nEok;
      fprintf(ft, "%d,%.4f\n", j, nt); fprintf(fe, "%d,%.6f\n", j, ne);
      fprintf(fs, "%d,%d,%.4f,%.4f,%.4f,%d,%d,%.2f,%.2f,%.3f,%.6f,%.6f,%d\n", j, nO[j], corr, havePrev_ ? prevT_[j] : 0., nt, (int)tok, nE[j], mpv[j], mpvErr[j], chi2[j], pe, ne, (int)eok);
    }
    for (int j = NS; j < (int)CaloConst::_nChannelDB; ++j) { fprintf(ft, "%d,0.0000\n", j); fprintf(fe, "%d,%.6f\n", j, medE); }   // spare channels of the DB tables
    fclose(ft); fclose(fe); fclose(fs);
    FILE* fsim = fopen(simFile_.c_str(), "w"); FILE* fb = fopen(baseFile_.c_str(), "w"); FILE* fr = fopen(respFile_.c_str(), "w");
    if (!fsim || !fb || !fr) throw cet::exception("MC2CALOCALIB") << "cannot write the response tables\n";
    fprintf(fsim, "TABLE CalSimCrystals\n#crid,LRU,pePerMeV0,pePerMeV1,ADCPerMeV0,ADCPerMeV1  MC2CaloCalibFromTracks: pePerMeV = ADCPerMeV / p0 (L-R photostatistics; median p0 %.3f ADC/pe where not measured), LRU = %.3f not measured\n", p0med, lru_);
    for (int c = 0; c < NC; ++c) { const double pp = p0[c] > 0 ? p0[c] : p0med; const double ad0 = newE[2 * c] > 0 ? 1. / newE[2 * c] : 0, ad1 = newE[2 * c + 1] > 0 ? 1. / newE[2 * c + 1] : 0;
      fprintf(fsim, "%d,%.4f,%.3f,%.3f,%.4f,%.4f\n", c, lru_, pp > 0 ? ad0 / pp : 0., pp > 0 ? ad1 / pp : 0., ad0, ad1); }
    fprintf(fb, "TABLE CalBaselines\n#roid,baseline,threshold  MC2CaloCalibFromTracks: baseline = mean of the first %d samples [ADC], threshold = rising edge of the raw amplitude (largest sample - baseline) [ADC]\n", nPed_);
    fprintf(fr, "roid,crid,nPulses,baseline_adc,noise_adc,thrRaw_adc,thrRawWidth_adc,thrFit_adc,thrFitWidth_adc,ADC2MeV,thr_MeV,thrWidth_MeV,p0_adc_per_pe,p0_err,p0_measured,pePerMeV\n");
    int nThr = 0;
    for (int j = 0; j < NS; ++j) { const double base = bN_[j] ? bSum_[j] / bN_[j] : 0, noise = bN_[j] ? std::sqrt(b2Sum_[j] / bN_[j]) : 0; const int c = j / 2;
      fprintf(fb, "%d,%.3f,%.3f\n", j, base, thrRaw[j] > 0 ? thrRaw[j] : 0.); if (thrFit[j] > 0) ++nThr;
      const double pp = (c < NC && p0[c] > 0) ? p0[c] : p0med;
      fprintf(fr, "%d,%d,%ld,%.3f,%.3f,%.2f,%.2f,%.2f,%.2f,%.6f,%.3f,%.3f,%.4f,%.4f,%d,%.3f\n", j, c, bN_[j], base, noise, thrRaw[j], thrRawW[j], thrFit[j], thrFitW[j], newE[j],
              thrFit[j] > 0 ? thrFit[j] * newE[j] : -1., thrFitW[j] > 0 ? thrFitW[j] * newE[j] : -1., pp, (c < NC) ? p0e[c] : -1., (int)(c < NC && p0[c] > 0), (pp > 0 && newE[j] > 0) ? 1. / newE[j] / pp : 0.); }
    for (int j = NS; j < (int)CaloConst::_nChannelDB; ++j) fprintf(fb, "%d,0.000,0.000\n", j);
    fclose(fsim); fclose(fb); fclose(fr);
    printf("MC2CaloCalibFromTracks: thresholds for %d SiPMs (median %.1f ADC raw); p0 measured for %zu crystals, median %.3f ADC per effective photoelectron -> %s, %s, %s\n",
           nThr, hThr_->GetMean(), allp.size(), p0med, simFile_.c_str(), baseFile_.c_str(), respFile_.c_str());
    printf("MC2CaloCalibFromTracks: %d tracks, %zu SiPM time entries; time offsets for %d SiPMs, MIP fits for %d SiPMs -> %s, %s, %s\n",
           ntrk_, te_.size(), nTok, nEok, tFile_.c_str(), eFile_.c_str(), sFile_.c_str());
    printf("MC2CaloCalibFromTracks: residual rms before %.3f ns, after %.3f ns (within +-20 ns)\n", hTresBefore_->GetRMS(), hTres_->GetRMS());
  }
}
DEFINE_ART_MODULE(mu2e::MC2CaloCalibFromTracks)
