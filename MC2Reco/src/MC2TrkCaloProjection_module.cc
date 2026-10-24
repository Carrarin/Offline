//
// MC2TrkCaloProjection: projection of the straight-line tracker tracks onto the calorimeter disks, and search for the
// calorimeter hits of the same cosmic, in time and in space.
//
// Inputs:
//   KalSeedCollection  (kalSeedCollection, default "KKLine"): KalSeeds of the KinKal KinematicLineFit
//   CaloHitCollection  (caloHitCollection, default "CaloHitMaker")
//
// Steps, for every event:
//
//   1. z of the crystal front face of each disk (Offline Calorimeter geometry, tracker frame)
//   2. calorimeter hits of the event with E >= eMin, with the crystal position in the tracker frame
//   3. for each KalSeed with the KKLine flag:
//        line  = the KinematicLine of the segment closest to the calorimeter (largest z); straight line, speed from the
//                fit (no field: the momentum is the seed value, so the speed is ~c)
//      and for each disk:
//        plane = front-face z + depthInCrystal + diskZOffset[disk] (diskZOffset: to follow a set-up that differs from
//                the Offline geometry)
//        at the plane: predicted x, y, t; inRing = radius between the inner and outer crystal radius of the disk (fcl)
//        matching: the calo hits of that disk with |t_hit - t_pred - dtOffset| < dtWindow; their energy-weighted
//                centroid (crystal x, y) and time, the distance of the nearest hit crystal, the summed energy
//
// Output (TFileService):
//   tree "proj", one entry per (track, disk):
//     run, subrun, event, itrk (index of the KalSeed), ntrk (KalSeeds in the event), disk
//     nactive, ndof, chi2, fitcon         fit quality of the KalSeed
//     x, y, z, t                          point of the line at the plane
//     ux, uy, uz, speed                   direction and speed of the line, so that a different plane can be tried
//                                         offline
//     inRing                              the point is inside the crystal ring
//     nMatch, eSum, xc, yc, tMatch        matched hits: number, summed energy, energy-weighted x, y, t (NAN if none)
//     dNear, idNear                       distance in x-y and crystal id of the nearest matched hit (NAN, -1 if none)
//     hitId, hitDisk, hitX, hitY, hitZ, hitT, hitE   all the calo hits of the event (vectors), for z scans and checks
//   histograms:
//     hDtWide, hDt          t_calo - t_track of all (track, hit) pairs in the same disk (wide and zoom)
//     hDx_disk<d>, hDy_disk<d>   x and y of the matched centroid minus the track, for tracks in the ring
//
// Main parameters: depthInCrystal, diskZOffset, rIn, rOut, eMin, dtOffset, dtWindow.
//
// Frames: tracker (detector) frame, mm, ns. Crystal positions as in EventNtuple: cal.mu2eToTracker(crystal.position()),
// where Crystal::position() is the centre of the front face.
//
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art_root_io/TFileService.h"
#include "cetlib_except/exception.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"

#include "Offline/CalorimeterGeom/inc/Calorimeter.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/RecoDataProducts/inc/CaloHit.hh"
#include "Offline/RecoDataProducts/inc/KalSeed.hh"
#include "KinKal/Trajectory/KinematicLine.hh"
#include "TH1F.h"
#include "TString.h"
#include "TTree.h"

#include <cmath>
#include <cstdio>
#include <vector>

namespace mu2e {

  class MC2TrkCaloProjection : public art::EDAnalyzer {
    public:
      struct Config {
        using Name = fhicl::Name;
        using Comment = fhicl::Comment;

        // Inputs
        fhicl::Atom<art::InputTag> kalSeedTag{Name("kalSeedCollection"),
            Comment("KalSeeds of the straight-line fit"), "KKLine"};
        fhicl::Atom<art::InputTag> caloHitTag{Name("caloHitCollection"), Comment("calorimeter hits"), "CaloHitMaker"};

        // Projection plane
        fhicl::Atom<double> depthInCrystal{Name("depthInCrystal"),
            Comment("plane depth behind the crystal front face [mm]"), 100.};
        fhicl::Sequence<double> diskZOffset{Name("diskZOffset"),
            Comment("added to the Offline z of each disk [mm]"), std::vector<double>{0., 0.}};

        // Crystal ring
        fhicl::Atom<double> rIn{Name("rIn"), Comment("inner radius of the crystal ring [mm]"), 374.};
        fhicl::Atom<double> rOut{Name("rOut"), Comment("outer radius of the crystal ring [mm]"), 660.};

        // Matching of the calo hits
        fhicl::Atom<double> eMin{Name("eMin"), Comment("calo hit threshold [MeV]"), 5.};
        fhicl::Atom<double> dtOffset{Name("dtOffset"), Comment("t_calo - t_track offset [ns]"), 0.};
        fhicl::Atom<double> dtWindow{Name("dtWindow"), Comment("half width of the time window [ns]"), 25.};

        fhicl::Atom<int> diagLevel{Name("diagLevel"), Comment("diagnostics"), 0};
      };
      using Parameters = art::EDAnalyzer::Table<Config>;

      explicit MC2TrkCaloProjection(const Parameters& config);
      void beginJob() override;
      void analyze(const art::Event& event) override;

    private:
      void bookTree(art::TFileService& tfs);
      void bookHistograms(art::TFileService& tfs);
      static void findDiskFrontFaces(const Calorimeter& cal, double zFace[2]);
      void collectCaloHits(const CaloHitCollection& hits, const Calorimeter& cal);
      static const KalSegment* segmentClosestToCalo(const KalSeed& seed);
      void projectOntoPlane(double zPlane, double tRef, const KinKal::VEC3& point, const KinKal::VEC3& velocity);
      void matchCaloHits(int disk);

      // Inputs
      art::InputTag kalSeedTag_;
      art::InputTag caloHitTag_;

      // Parameters
      double depthInCrystal_;
      double rIn_;
      double rOut_;
      double eMin_;
      double dtOffset_;
      double dtWindow_;
      std::vector<double> diskZOffset_;
      int diagLevel_;

      // Output: tree and histograms
      TTree* tree_ = nullptr;
      TH1F* hDtWide_ = nullptr;
      TH1F* hDt_ = nullptr;
      TH1F* hDx_[2] = {nullptr, nullptr};
      TH1F* hDy_[2] = {nullptr, nullptr};

      // Tree variables: event and track
      int run_;
      int subRun_;
      int event_;
      int trackIndex_;
      int nTracks_;
      int disk_;
      int nActive_;
      int nDof_;
      float chi2_;
      float fitConsistency_;
      // Tree variables: line at the plane
      float xPlane_;
      float yPlane_;
      float zPlane_;
      float tPlane_;
      float ux_;
      float uy_;
      float uz_;
      float speed_;
      int inRing_;
      // Tree variables: matched calo hits
      int nMatch_;
      float eSum_;
      float xCentroid_;
      float yCentroid_;
      float tMatch_;
      float dNearest_;
      int idNearest_;
      // Tree variables: all the calo hits of the event (E >= eMin), tracker frame
      std::vector<int> hitId_;
      std::vector<int> hitDisk_;
      std::vector<float> hitX_;
      std::vector<float> hitY_;
      std::vector<float> hitZ_;
      std::vector<float> hitT_;
      std::vector<float> hitE_;
  };

  MC2TrkCaloProjection::MC2TrkCaloProjection(const Parameters& config) :
    art::EDAnalyzer(config),
    kalSeedTag_(config().kalSeedTag()),
    caloHitTag_(config().caloHitTag()),
    depthInCrystal_(config().depthInCrystal()),
    rIn_(config().rIn()),
    rOut_(config().rOut()),
    eMin_(config().eMin()),
    dtOffset_(config().dtOffset()),
    dtWindow_(config().dtWindow()),
    diskZOffset_(config().diskZOffset()),
    diagLevel_(config().diagLevel())
  {
    if (diskZOffset_.size() != 2) {
      throw cet::exception("MC2TRKCALO") << "diskZOffset needs two values\n";
    }
    consumes<KalSeedCollection>(kalSeedTag_);
    consumes<CaloHitCollection>(caloHitTag_);
  }

  void MC2TrkCaloProjection::beginJob() {
    art::ServiceHandle<art::TFileService> tfs;
    bookTree(*tfs);
    bookHistograms(*tfs);
  }

  void MC2TrkCaloProjection::bookTree(art::TFileService& tfs) {
    tree_ = tfs.make<TTree>("proj",
        "tracker straight lines projected onto the calorimeter disks (one entry per track and disk)");
    tree_->Branch("run", &run_);
    tree_->Branch("subrun", &subRun_);
    tree_->Branch("event", &event_);
    tree_->Branch("itrk", &trackIndex_);
    tree_->Branch("ntrk", &nTracks_);
    tree_->Branch("disk", &disk_);
    tree_->Branch("nactive", &nActive_);
    tree_->Branch("ndof", &nDof_);
    tree_->Branch("chi2", &chi2_);
    tree_->Branch("fitcon", &fitConsistency_);
    tree_->Branch("x", &xPlane_);
    tree_->Branch("y", &yPlane_);
    tree_->Branch("z", &zPlane_);
    tree_->Branch("t", &tPlane_);
    tree_->Branch("ux", &ux_);
    tree_->Branch("uy", &uy_);
    tree_->Branch("uz", &uz_);
    tree_->Branch("speed", &speed_);
    tree_->Branch("inRing", &inRing_);
    tree_->Branch("nMatch", &nMatch_);
    tree_->Branch("eSum", &eSum_);
    tree_->Branch("xc", &xCentroid_);
    tree_->Branch("yc", &yCentroid_);
    tree_->Branch("tMatch", &tMatch_);
    tree_->Branch("dNear", &dNearest_);
    tree_->Branch("idNear", &idNearest_);
    tree_->Branch("hitId", &hitId_);
    tree_->Branch("hitDisk", &hitDisk_);
    tree_->Branch("hitX", &hitX_);
    tree_->Branch("hitY", &hitY_);
    tree_->Branch("hitZ", &hitZ_);
    tree_->Branch("hitT", &hitT_);
    tree_->Branch("hitE", &hitE_);
  }

  void MC2TrkCaloProjection::bookHistograms(art::TFileService& tfs) {
    hDtWide_ = tfs.make<TH1F>("hDtWide", "all (track, calo hit) pairs in the disk;t_{calo} - t_{track} [ns];pairs",
                              500, -250000., 250000.);
    hDt_ = tfs.make<TH1F>("hDt", "all (track, calo hit) pairs in the disk;t_{calo} - t_{track} [ns];pairs / 1 ns",
                          200, -100., 100.);
    for (int disk = 0; disk < 2; ++disk) {
      hDx_[disk] = tfs.make<TH1F>(Form("hDx_disk%d", disk),
          Form("disk %d, track in the ring, matched;x_{calo} - x_{track} [mm];tracks", disk), 120, -600., 600.);
      hDy_[disk] = tfs.make<TH1F>(Form("hDy_disk%d", disk),
          Form("disk %d, track in the ring, matched;y_{calo} - y_{track} [mm];tracks", disk), 120, -600., 600.);
    }
  }

  // Step 1: z (tracker frame) of the crystal front face of each disk, from the first crystal found in that disk.
  void MC2TrkCaloProjection::findDiskFrontFaces(const Calorimeter& cal, double zFace[2]) {
    zFace[0] = NAN;
    zFace[1] = NAN;
    for (unsigned i = 0; i < cal.nCrystals() && (std::isnan(zFace[0]) || std::isnan(zFace[1])); ++i) {
      const auto& crystal = cal.crystal(i);
      const int disk = crystal.diskID();
      if (disk < 0 || disk > 1 || !std::isnan(zFace[disk])) {
        continue;
      }
      zFace[disk] = cal.mu2eToTracker(crystal.position()).z();
    }
  }

  // Step 2: calo hits above eMin, with the front-face centre of their crystal in the tracker frame.
  void MC2TrkCaloProjection::collectCaloHits(const CaloHitCollection& hits, const Calorimeter& cal) {
    hitId_.clear();
    hitDisk_.clear();
    hitX_.clear();
    hitY_.clear();
    hitZ_.clear();
    hitT_.clear();
    hitE_.clear();
    for (const auto& hit : hits) {
      if (hit.energyDep() < eMin_) {
        continue;
      }
      const auto& crystal = cal.crystal(hit.crystalID());
      const auto position = cal.mu2eToTracker(crystal.position());
      hitId_.push_back(hit.crystalID());
      hitDisk_.push_back(crystal.diskID());
      hitX_.push_back(position.x());
      hitY_.push_back(position.y());
      hitZ_.push_back(position.z());
      hitT_.push_back(hit.time());
      hitE_.push_back(hit.energyDep());
    }
  }

  // Segment of the KalSeed closest to the calorimeter: the one with the largest z.
  const KalSegment* MC2TrkCaloProjection::segmentClosestToCalo(const KalSeed& seed) {
    const KalSegment* closest = nullptr;
    double zClosest = -1e30;
    for (const auto& segment : seed.segments()) {
      const double z = segment.position3().Z();
      if (z > zClosest) {
        zClosest = z;
        closest = &segment;
      }
    }
    return closest;
  }

  // Point of the straight line at z = zPlane (the tree variables are filled in float and used as such afterwards).
  void MC2TrkCaloProjection::projectOntoPlane(double zPlane, double tRef, const KinKal::VEC3& point,
                                              const KinKal::VEC3& velocity) {
    zPlane_ = zPlane;
    const double tAtPlane = tRef + (zPlane_ - point.Z()) / velocity.Z();
    xPlane_ = point.X() + velocity.X() * (tAtPlane - tRef);
    yPlane_ = point.Y() + velocity.Y() * (tAtPlane - tRef);
    tPlane_ = tAtPlane;
    const double radius = std::hypot(xPlane_, yPlane_);
    inRing_ = (radius > rIn_ && radius < rOut_);
  }

  // Calo hits of this disk in the time window around the predicted time: energy-weighted centroid and time, nearest
  // hit crystal. The dt histograms take all the hits of the disk, before the time window.
  void MC2TrkCaloProjection::matchCaloHits(int disk) {
    nMatch_ = 0;
    eSum_ = 0;
    double sumEx = 0;
    double sumEy = 0;
    double sumEt = 0;
    dNearest_ = 1e9;
    idNearest_ = -1;
    for (size_t j = 0; j < hitId_.size(); ++j) {
      if (hitDisk_[j] != disk) {
        continue;
      }
      const double dt = hitT_[j] - tPlane_;
      hDtWide_->Fill(dt);
      hDt_->Fill(dt);
      if (std::abs(dt - dtOffset_) > dtWindow_) {
        continue;
      }
      ++nMatch_;
      eSum_ += hitE_[j];
      sumEx += hitE_[j] * hitX_[j];
      sumEy += hitE_[j] * hitY_[j];
      sumEt += hitE_[j] * hitT_[j];
      const double distance = std::hypot(hitX_[j] - xPlane_, hitY_[j] - yPlane_);
      if (distance < dNearest_) {
        dNearest_ = distance;
        idNearest_ = hitId_[j];
      }
    }

    if (nMatch_ > 0) {
      xCentroid_ = sumEx / eSum_;
      yCentroid_ = sumEy / eSum_;
      tMatch_ = sumEt / eSum_;
      if (inRing_) {
        hDx_[disk]->Fill(xCentroid_ - xPlane_);
        hDy_[disk]->Fill(yCentroid_ - yPlane_);
      }
    } else {
      xCentroid_ = yCentroid_ = tMatch_ = NAN;
      dNearest_ = NAN;
    }
  }

  void MC2TrkCaloProjection::analyze(const art::Event& event) {
    const Calorimeter& cal = *GeomHandle<Calorimeter>();
    const auto& seeds = *event.getValidHandle<KalSeedCollection>(kalSeedTag_);
    const auto& hits = *event.getValidHandle<CaloHitCollection>(caloHitTag_);
    run_ = event.run();
    subRun_ = event.subRun();
    event_ = event.event();
    nTracks_ = seeds.size();

    // 1-2. front face of the disks and calo hits of the event
    double zFace[2];
    findDiskFrontFaces(cal, zFace);
    collectCaloHits(hits, cal);

    // 3. each straight-line track, projected onto each disk
    trackIndex_ = -1;
    for (const auto& seed : seeds) {
      ++trackIndex_;
      if (!seed.kinematicLineFit() || seed.segments().empty()) {
        continue;
      }
      const KalSegment* segment = segmentClosestToCalo(seed);
      const KinKal::KinematicLine line = segment->kinematicLine();
      const double tRef = 0.5 * (segment->tmin() + segment->tmax());
      const auto point = line.position3(tRef);
      const auto velocity = line.velocity(tRef);
      const double speed = std::sqrt(velocity.Mag2());
      if (speed <= 0 || std::abs(velocity.Z()) < 1e-9) {
        continue;
      }

      nActive_ = seed.nHits(true);
      nDof_ = seed.nDOF();
      chi2_ = seed.chisquared();
      fitConsistency_ = seed.fitConsistency();
      ux_ = velocity.X() / speed;
      uy_ = velocity.Y() / speed;
      uz_ = velocity.Z() / speed;
      speed_ = speed;

      for (int disk = 0; disk < 2; ++disk) {
        disk_ = disk;
        projectOntoPlane(zFace[disk] + depthInCrystal_ + diskZOffset_[disk], tRef, point, velocity);
        matchCaloHits(disk);
        tree_->Fill();
      }
      if (diagLevel_ > 0) {
        printf("MC2TrkCaloProjection: event %d track %d nactive %d\n", event_, trackIndex_, nActive_);
      }
    }
  }

}

DEFINE_ART_MODULE(mu2e::MC2TrkCaloProjection)
