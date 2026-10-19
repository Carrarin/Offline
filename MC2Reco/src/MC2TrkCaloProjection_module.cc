// MC2TrkCaloProjection: projects the straight-line tracker tracks (KalSeeds of KinKal KinematicLineFit, "KKLine") onto the
// calorimeter disks and looks for the calorimeter hits of the same cosmic, in time and in space.
//
// For each KalSeed with the KKLine flag:
//   line   = the KinematicLine of the segment closest to the calorimeter (largest z); straight line, speed from the fit
//            (no field: the momentum is the seed value, so the speed is ~c);
//   plane  = z of the crystal front face of each disk (Offline Calorimeter geometry, tracker frame) + depthInCrystal
//            + diskZOffset[disk] (fcl: to follow a set-up that differs from the Offline geometry);
//   at the plane: predicted x, y, t; inRing = radius between the inner and outer crystal radius of the disk (fcl);
//   matching: calo hits of that disk with |t_hit - t_pred - dtOffset| < dtWindow; their energy-weighted centroid (crystal
//            x, y), the distance of the nearest hit crystal, the summed energy.
// Output (TFileService):
//   tree "proj": one entry per (track, disk), with the line itself (point and direction at the plane) so that a different
//                plane can be tried offline, and all the calo hits of the event (vectors), for z scans and checks;
//   histograms: dt of all (track, hit) pairs in the disk (wide and zoom), dx and dy of the matched centroid.
// Frames: tracker (detector) frame, mm, ns. Crystal positions as in EventNtuple: cal.mu2eToTracker(crystal.position()).
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art_root_io/TFileService.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/CalorimeterGeom/inc/Calorimeter.hh"
#include "Offline/RecoDataProducts/inc/KalSeed.hh"
#include "Offline/RecoDataProducts/inc/CaloHit.hh"
#include "KinKal/Trajectory/KinematicLine.hh"
#include "TTree.h"
#include "TH1F.h"
#include <cmath>
#include <vector>

namespace mu2e {
  class MC2TrkCaloProjection : public art::EDAnalyzer {
    public:
      struct Config {
        using Name = fhicl::Name; using Comment = fhicl::Comment;
        fhicl::Atom<art::InputTag>  kalSeedTag     { Name("kalSeedCollection"), Comment("KalSeeds of the straight-line fit"), "KKLine" };
        fhicl::Atom<art::InputTag>  caloHitTag     { Name("caloHitCollection"), Comment("calorimeter hits"), "CaloHitMaker" };
        fhicl::Atom<double>         depthInCrystal { Name("depthInCrystal"), Comment("plane depth behind the crystal front face [mm]"), 100. };
        fhicl::Sequence<double>     diskZOffset    { Name("diskZOffset"), Comment("added to the Offline z of each disk [mm]"), std::vector<double>{0., 0.} };
        fhicl::Atom<double>         rIn            { Name("rIn"),  Comment("inner radius of the crystal ring [mm]"), 374. };
        fhicl::Atom<double>         rOut           { Name("rOut"), Comment("outer radius of the crystal ring [mm]"), 660. };
        fhicl::Atom<double>         eMin           { Name("eMin"), Comment("calo hit threshold [MeV]"), 5. };
        fhicl::Atom<double>         dtOffset       { Name("dtOffset"), Comment("t_calo - t_track offset [ns]"), 0. };
        fhicl::Atom<double>         dtWindow       { Name("dtWindow"), Comment("half width of the time window [ns]"), 25. };
        fhicl::Atom<int>            diagLevel      { Name("diagLevel"), Comment("diagnostics"), 0 };
      };
      using Parameters = art::EDAnalyzer::Table<Config>;
      explicit MC2TrkCaloProjection(const Parameters& conf);
      void beginJob() override;
      void analyze(const art::Event& event) override;
    private:
      art::InputTag kalTag_, hitTag_;
      double depth_, rIn_, rOut_, eMin_, dtOff_, dtWin_; std::vector<double> dzOff_; int diag_;
      TTree* t_ = nullptr;
      TH1F *hDtWide_ = nullptr, *hDt_ = nullptr, *hDx_[2] = {nullptr, nullptr}, *hDy_[2] = {nullptr, nullptr};
      // tree variables
      int run_, sub_, evt_, itrk_, ntrk_, disk_, nact_, ndof_, inRing_, nMatch_, idNear_;
      float chi2_, fitcon_, xp_, yp_, zp_, tp_, ux_, uy_, uz_, speed_, eSum_, xc_, yc_, dNear_, tMatch_;
      std::vector<int> hId_, hDisk_; std::vector<float> hX_, hY_, hZ_, hT_, hE_;
  };

  MC2TrkCaloProjection::MC2TrkCaloProjection(const Parameters& conf) : art::EDAnalyzer(conf),
    kalTag_(conf().kalSeedTag()), hitTag_(conf().caloHitTag()), depth_(conf().depthInCrystal()), rIn_(conf().rIn()), rOut_(conf().rOut()),
    eMin_(conf().eMin()), dtOff_(conf().dtOffset()), dtWin_(conf().dtWindow()), dzOff_(conf().diskZOffset()), diag_(conf().diagLevel()) {
    if (dzOff_.size() != 2) throw cet::exception("MC2TRKCALO") << "diskZOffset needs two values\n";
    consumes<KalSeedCollection>(kalTag_); consumes<CaloHitCollection>(hitTag_);
  }

  void MC2TrkCaloProjection::beginJob() {
    art::ServiceHandle<art::TFileService> tfs;
    t_ = tfs->make<TTree>("proj", "tracker straight lines projected onto the calorimeter disks (one entry per track and disk)");
    t_->Branch("run", &run_); t_->Branch("subrun", &sub_); t_->Branch("event", &evt_); t_->Branch("itrk", &itrk_); t_->Branch("ntrk", &ntrk_);
    t_->Branch("disk", &disk_); t_->Branch("nactive", &nact_); t_->Branch("ndof", &ndof_); t_->Branch("chi2", &chi2_); t_->Branch("fitcon", &fitcon_);
    t_->Branch("x", &xp_); t_->Branch("y", &yp_); t_->Branch("z", &zp_); t_->Branch("t", &tp_);
    t_->Branch("ux", &ux_); t_->Branch("uy", &uy_); t_->Branch("uz", &uz_); t_->Branch("speed", &speed_);
    t_->Branch("inRing", &inRing_); t_->Branch("nMatch", &nMatch_); t_->Branch("eSum", &eSum_); t_->Branch("xc", &xc_); t_->Branch("yc", &yc_);
    t_->Branch("tMatch", &tMatch_); t_->Branch("dNear", &dNear_); t_->Branch("idNear", &idNear_);
    t_->Branch("hitId", &hId_); t_->Branch("hitDisk", &hDisk_); t_->Branch("hitX", &hX_); t_->Branch("hitY", &hY_); t_->Branch("hitZ", &hZ_);
    t_->Branch("hitT", &hT_); t_->Branch("hitE", &hE_);
    hDtWide_ = tfs->make<TH1F>("hDtWide", "all (track, calo hit) pairs in the disk;t_{calo} - t_{track} [ns];pairs", 500, -250000., 250000.);
    hDt_ = tfs->make<TH1F>("hDt", "all (track, calo hit) pairs in the disk;t_{calo} - t_{track} [ns];pairs / 1 ns", 200, -100., 100.);
    for (int d = 0; d < 2; ++d) {
      hDx_[d] = tfs->make<TH1F>(Form("hDx_disk%d", d), Form("disk %d, track in the ring, matched;x_{calo} - x_{track} [mm];tracks", d), 120, -600., 600.);
      hDy_[d] = tfs->make<TH1F>(Form("hDy_disk%d", d), Form("disk %d, track in the ring, matched;y_{calo} - y_{track} [mm];tracks", d), 120, -600., 600.);
    }
  }

  void MC2TrkCaloProjection::analyze(const art::Event& event) {
    const Calorimeter& cal = *GeomHandle<Calorimeter>();
    const auto& seeds = *event.getValidHandle<KalSeedCollection>(kalTag_);
    const auto& hits = *event.getValidHandle<CaloHitCollection>(hitTag_);
    run_ = event.run(); sub_ = event.subRun(); evt_ = event.event(); ntrk_ = seeds.size();
    // calo hits of the event, tracker frame; and the front-face z of each disk
    hId_.clear(); hDisk_.clear(); hX_.clear(); hY_.clear(); hZ_.clear(); hT_.clear(); hE_.clear();
    double zFace[2] = {NAN, NAN};
    for (unsigned i = 0; i < cal.nCrystals() && (std::isnan(zFace[0]) || std::isnan(zFace[1])); ++i) {
      const auto& c = cal.crystal(i); const int d = c.diskID(); if (d < 0 || d > 1 || !std::isnan(zFace[d])) continue;
      zFace[d] = cal.mu2eToTracker(c.position()).z();
    }
    for (const auto& h : hits) {
      if (h.energyDep() < eMin_) continue;
      const auto& c = cal.crystal(h.crystalID()); const auto p = cal.mu2eToTracker(c.position());
      hId_.push_back(h.crystalID()); hDisk_.push_back(c.diskID()); hX_.push_back(p.x()); hY_.push_back(p.y()); hZ_.push_back(p.z());
      hT_.push_back(h.time()); hE_.push_back(h.energyDep());
    }
    itrk_ = -1;
    for (const auto& ks : seeds) {
      ++itrk_;
      if (!ks.kinematicLineFit() || ks.segments().empty()) continue;
      // segment closest to the calorimeter (largest z)
      const KalSegment* best = nullptr; double zbest = -1e30;
      for (const auto& s : ks.segments()) { const double z = s.position3().Z(); if (z > zbest) { zbest = z; best = &s; } }
      const KinKal::KinematicLine line = best->kinematicLine();
      const double tref = 0.5 * (best->tmin() + best->tmax()); const auto p0 = line.position3(tref); const auto v = line.velocity(tref);
      const double vmag = std::sqrt(v.Mag2()); if (vmag <= 0 || std::abs(v.Z()) < 1e-9) continue;
      nact_ = ks.nHits(true); ndof_ = ks.nDOF(); chi2_ = ks.chisquared(); fitcon_ = ks.fitConsistency();
      ux_ = v.X() / vmag; uy_ = v.Y() / vmag; uz_ = v.Z() / vmag; speed_ = vmag;
      for (int d = 0; d < 2; ++d) {
        disk_ = d; zp_ = zFace[d] + depth_ + dzOff_[d];
        const double tt = tref + (zp_ - p0.Z()) / v.Z(); xp_ = p0.X() + v.X() * (tt - tref); yp_ = p0.Y() + v.Y() * (tt - tref); tp_ = tt;
        const double r = std::hypot(xp_, yp_); inRing_ = (r > rIn_ && r < rOut_);
        nMatch_ = 0; eSum_ = 0; double sx = 0, sy = 0, st = 0; dNear_ = 1e9; idNear_ = -1;
        for (size_t j = 0; j < hId_.size(); ++j) {
          if (hDisk_[j] != d) continue;
          const double dt = hT_[j] - tp_; hDtWide_->Fill(dt); hDt_->Fill(dt);
          if (std::abs(dt - dtOff_) > dtWin_) continue;
          ++nMatch_; eSum_ += hE_[j]; sx += hE_[j] * hX_[j]; sy += hE_[j] * hY_[j]; st += hE_[j] * hT_[j];
          const double dd = std::hypot(hX_[j] - xp_, hY_[j] - yp_); if (dd < dNear_) { dNear_ = dd; idNear_ = hId_[j]; }
        }
        if (nMatch_ > 0) { xc_ = sx / eSum_; yc_ = sy / eSum_; tMatch_ = st / eSum_; if (inRing_) { hDx_[d]->Fill(xc_ - xp_); hDy_[d]->Fill(yc_ - yp_); } }
        else { xc_ = yc_ = tMatch_ = NAN; dNear_ = NAN; }
        t_->Fill();
      }
      if (diag_ > 0) printf("MC2TrkCaloProjection: event %d track %d nactive %d\n", evt_, itrk_, nact_);
    }
  }
}
DEFINE_ART_MODULE(mu2e::MC2TrkCaloProjection)
