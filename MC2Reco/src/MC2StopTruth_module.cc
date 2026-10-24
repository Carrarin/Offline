//
// MC2StopTruth: Monte Carlo truth of the muons that stop in the calorimeter crystals (hall production, Mu2eG4): what
// the decay or the capture leaves, and when. One tree entry per stopped mu+ or mu-.
//
// Input: SimParticles and CaloShowerSteps as saved by the Primary job (label compressDetStepMCs).
//
// Definitions:
//   Stopped muon : SimParticle with |pdg| = 13, end kinetic energy <= maxEndKE, with CaloShowerSteps of its own.
//                  If its end point is inside a crystal box (34 x 34 x 200 mm around the crystal centre; the centre is
//                  the front face Crystal::position() + half the crystal length along the crystal axis), inCry = 1 and
//                  that is the stop crystal; otherwise (stopped in the material around the crystals) inCry = 0 and the
//                  reference crystal is the last one the muon crossed. Its own deposit (steps up to its end + 1 ns) is
//                  the prompt deposit. (Most muons with calo steps that reach zero energy stop outside the crystals:
//                  floor, supports.)
//   Daughters    : first generation (parent = the muon): pdg, creation code, kinetic energy at the start, start time -
//                  muon end time. In Geant4 11 a mu- is first captured in an atom (muonic atom), which then decays
//                  (DIO electron) or is captured by the nucleus: the delay of the daughters measures that lifetime.
//   Deposits     : CaloShowerSteps of every descendant (parent chain up to the muon) AND the muon's own steps later
//                  than its end + 1 ns (CaloShowerStepMaker gives the energy of the particles born inside the
//                  calorimeter to the particle that entered it), by crystal: stop crystal, neighbours (front-face
//                  centre within dNear mm, same disk), far; energy (G4, MeV), first time and energy-weighted time after
//                  the stop; also in time windows after the stop: 0-50 ns, 50-500 ns, 0.5-20 us, > 20 us.
//
// Output: TTree "stops", one entry per stopped muon:
//   run, subrun, event
//   q                        muon charge (-1 mu-, +1 mu+)
//   cry, disk                stop (or reference) crystal and its disk
//   inCry                    1 if the end point is inside the box of crystal cry
//   x, y                     front-face centre of crystal cry (tracker frame) [mm]
//   endX, endY, endZ         muon end point (tracker frame) [mm]
//   endKE, endCode           muon end kinetic energy [MeV] and stopping process code
//   tStop                    muon end time [ns]
//   dtEnd                    muon end time - time of its last calo step [ns]
//   eMu, nCryMu              energy of all the muon's own steps [MeV], number of crystals they touch
//   ePrompt                  energy of the muon's own steps up to tStop + 1 ns [MeV]
//   nDau, tDauMin            number of first-generation daughters, earliest daughter delay [ns] (1e9 if none)
//   dPdg, dCode, dKE, dT     per daughter: pdg, creation code, start kinetic energy [MeV], delay [ns]
//   eTot, eStop, eNear, eFar delayed deposits [MeV]: all crystals, stop crystal, neighbours, others
//   tDep1                    delay of the first delayed step with more than 0.01 MeV [ns] (1e9 if none)
//   tDepMean                 energy-weighted delay of the delayed steps [ns] (-1 if no energy)
//   eWin, eStopWin, eNearWin, eFarWin [4]   the same energies in the four time windows
// At the end of the job: number of stopped mu- and mu+ and the creation processes of their daughters (printed).
//
// Parameters: maxEndKE, dNear (stopDt is accepted but not used, for old configurations).
//
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art_root_io/TFileService.h"
#include "fhiclcpp/types/Atom.h"

#include "Offline/CalorimeterGeom/inc/Calorimeter.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/MCDataProducts/inc/CaloShowerStep.hh"
#include "Offline/MCDataProducts/inc/ProcessCode.hh"
#include "Offline/MCDataProducts/inc/SimParticle.hh"
#include "TTree.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <map>
#include <set>
#include <vector>

namespace mu2e {

  class MC2StopTruth : public art::EDAnalyzer {
    public:
      struct Config {
        using Name = fhicl::Name;
        using Comment = fhicl::Comment;

        // Input
        fhicl::Atom<art::InputTag> simTag{Name("simParticles"), Comment("SimParticle collection"),
            "compressDetStepMCs"};
        fhicl::Atom<art::InputTag> stepTag{Name("caloShowerSteps"), Comment("CaloShowerStep collection"),
            "compressDetStepMCs"};

        // Selection and deposits
        fhicl::Atom<double> maxEndKE{Name("maxEndKE"), Comment("max end kinetic energy of a stopped muon [MeV]"), 1.};
        fhicl::Atom<double> stopDt{Name("stopDt"), Comment("not used (kept for old configurations)"), 5.};
        fhicl::Atom<double> dNear{Name("dNear"), Comment("neighbour crystal: centre within [mm] of the stop crystal"),
            80.};
      };
      using Parameters = art::EDAnalyzer::Table<Config>;

      explicit MC2StopTruth(const Parameters& config);
      void beginJob() override;
      void analyze(const art::Event& event) override;
      void endJob() override;

    private:
      using StepsByParticle = std::map<unsigned long, std::vector<const CaloShowerStep*>>;

      static StepsByParticle stepsByParticle(const CaloShowerStepCollection& steps);
      int crystalContaining(const Calorimeter& calorimeter, const CLHEP::Hep3Vector& pointInTracker) const;
      void fillMuon(const art::Event& event, const SimParticle& muon, const CLHEP::Hep3Vector& endInTracker,
                    const CaloShowerStep* lastStep, int crystal, int disk, const CLHEP::Hep3Vector& crystalFront,
                    double muonEnergy, size_t nMuonCrystals);
      void fillDaughters(const SimParticleCollection& sims, const cet::map_vector_key& muonKey);
      void resetDeposits();
      void addDeposit(const Calorimeter& calorimeter, const CaloShowerStep& step, int stopCrystal, int stopDisk,
                      const CLHEP::Hep3Vector& stopCrystalFront, double& energyTimeSum);
      static bool isDescendant(const SimParticle& particle, const cet::map_vector_key& muonKey);

      // Crystal box around the crystal centre [mm]
      static constexpr double halfXY_ = 17.;
      static constexpr double halfZ_ = 100.;

      // Input
      art::InputTag simTag_;
      art::InputTag stepTag_;

      // Parameters
      double maxEndKE_;
      double dNear_;

      // Output tree "stops" (one member per branch)
      TTree* tree_ = nullptr;
      int run_, subrun_, event_;
      int q_, cry_, disk_, inCry_;
      float ePrompt_, dtEnd_;
      int endCode_;
      float endX_, endY_, endZ_, endKE_;
      float x_, y_, tStop_, eMu_;
      int nCryMu_;
      int nDau_;
      float tDauMin_;
      std::vector<int> dPdg_, dCode_;
      std::vector<float> dKE_, dT_;
      float eTot_, eStop_, eNear_, eFar_, tDep1_, tDepMean_;
      float eWin_[4], eStopWin_[4], eNearWin_[4], eFarWin_[4];   // windows after the stop: 0-50 ns, 50-500 ns,
                                                                  // 0.5-20 us, > 20 us

      // Counters, [0] mu-, [1] mu+
      std::map<int, long> daughterCodeCount_[2];   // creation code -> number of daughters
      long nStopped_[2] = {0, 0};
  };

  MC2StopTruth::MC2StopTruth(const Parameters& config) :
    art::EDAnalyzer(config),
    simTag_(config().simTag()),
    stepTag_(config().stepTag()),
    maxEndKE_(config().maxEndKE()),
    dNear_(config().dNear())
  {
    consumes<SimParticleCollection>(simTag_);
    consumes<CaloShowerStepCollection>(stepTag_);
  }

  void MC2StopTruth::beginJob() {
    art::ServiceHandle<art::TFileService> tfs;
    tree_ = tfs->make<TTree>("stops", "muons stopped in the crystals: truth of the decay / capture products");
    tree_->Branch("run", &run_);
    tree_->Branch("subrun", &subrun_);
    tree_->Branch("event", &event_);
    tree_->Branch("q", &q_);
    tree_->Branch("cry", &cry_);
    tree_->Branch("disk", &disk_);
    tree_->Branch("inCry", &inCry_);
    tree_->Branch("ePrompt", &ePrompt_);
    tree_->Branch("dtEnd", &dtEnd_);
    tree_->Branch("endCode", &endCode_);
    tree_->Branch("endX", &endX_);
    tree_->Branch("endY", &endY_);
    tree_->Branch("endZ", &endZ_);
    tree_->Branch("endKE", &endKE_);
    tree_->Branch("x", &x_);
    tree_->Branch("y", &y_);
    tree_->Branch("tStop", &tStop_);
    tree_->Branch("eMu", &eMu_);
    tree_->Branch("nCryMu", &nCryMu_);
    tree_->Branch("nDau", &nDau_);
    tree_->Branch("tDauMin", &tDauMin_);
    tree_->Branch("dPdg", &dPdg_);
    tree_->Branch("dCode", &dCode_);
    tree_->Branch("dKE", &dKE_);
    tree_->Branch("dT", &dT_);
    tree_->Branch("eTot", &eTot_);
    tree_->Branch("eStop", &eStop_);
    tree_->Branch("eNear", &eNear_);
    tree_->Branch("eFar", &eFar_);
    tree_->Branch("tDep1", &tDep1_);
    tree_->Branch("tDepMean", &tDepMean_);
    tree_->Branch("eWin", eWin_, "eWin[4]/F");
    tree_->Branch("eStopWin", eStopWin_, "eStopWin[4]/F");
    tree_->Branch("eNearWin", eNearWin_, "eNearWin[4]/F");
    tree_->Branch("eFarWin", eFarWin_, "eFarWin[4]/F");
  }

  void MC2StopTruth::analyze(const art::Event& event) {
    const Calorimeter& calorimeter = *GeomHandle<Calorimeter>();
    const auto simHandle = event.getValidHandle<SimParticleCollection>(simTag_);
    const auto& sims = *simHandle;
    const auto& steps = *event.getValidHandle<CaloShowerStepCollection>(stepTag_);
    const StepsByParticle stepsOf = stepsByParticle(steps);

    for (const auto& [muonKey, muon] : sims) {
      // 1. stopped muon with calo steps of its own
      if (std::abs(muon.pdgId()) != 13 || muon.endKineticEnergy() > maxEndKE_) {
        continue;
      }
      const auto muonSteps = stepsOf.find(muonKey.asUint());
      if (muonSteps == stepsOf.end()) {
        continue;
      }

      // 2. the muon's own steps: energy, crystals, last step in time
      double muonEnergy = 0;
      std::set<int> muonCrystals;
      const CaloShowerStep* lastStep = nullptr;
      for (auto* step : muonSteps->second) {
        muonEnergy += step->energyDepG4();
        muonCrystals.insert(step->volumeG4ID());
        if (!lastStep || step->time() > lastStep->time()) {
          lastStep = step;
        }
      }

      // 3. stop crystal: the crystal box that contains the end point; otherwise the last crystal the muon crossed
      const auto endInTracker = calorimeter.mu2eToTracker(muon.endPosition());
      int stopCrystal = crystalContaining(calorimeter, endInTracker);
      inCry_ = stopCrystal >= 0;
      if (stopCrystal < 0) {
        stopCrystal = lastStep->volumeG4ID();
      }
      if (stopCrystal < 0 || stopCrystal >= (int)calorimeter.nCrystals()) {
        continue;
      }
      const auto stopCrystalFront = calorimeter.mu2eToTracker(calorimeter.crystal(stopCrystal).position());
      const int stopDisk = int(calorimeter.crystal(stopCrystal).diskID());

      fillMuon(event, muon, endInTracker, lastStep, stopCrystal, stopDisk, stopCrystalFront, muonEnergy,
               muonCrystals.size());

      // 4. first-generation daughters
      fillDaughters(sims, muonKey);

      // 5. delayed deposits: the muon's own steps after its end + 1 ns (the others are prompt), then all descendants
      resetDeposits();
      double energyTimeSum = 0;
      for (auto* step : muonSteps->second) {
        if (step->time() > tStop_ + 1.) {
          addDeposit(calorimeter, *step, stopCrystal, stopDisk, stopCrystalFront, energyTimeSum);
        } else {
          ePrompt_ += step->energyDepG4();
        }
      }
      for (const auto& [key, particle] : sims) {
        if (key == muonKey) {
          continue;
        }
        if (!isDescendant(particle, muonKey)) {
          continue;
        }
        const auto particleSteps = stepsOf.find(key.asUint());
        if (particleSteps == stepsOf.end()) {
          continue;
        }
        for (auto* step : particleSteps->second) {
          addDeposit(calorimeter, *step, stopCrystal, stopDisk, stopCrystalFront, energyTimeSum);
        }
      }
      tDepMean_ = eTot_ > 0 ? energyTimeSum / eTot_ : -1;

      tree_->Fill();
    }
  }

  // CaloShowerSteps of each SimParticle (by key).
  MC2StopTruth::StepsByParticle MC2StopTruth::stepsByParticle(const CaloShowerStepCollection& steps) {
    StepsByParticle stepsOf;
    for (const auto& step : steps) {
      if (step.simParticle().isNonnull()) {
        stepsOf[step.simParticle().key()].push_back(&step);
      }
    }
    return stepsOf;
  }

  // First crystal whose box (34 x 34 x 200 mm around its centre) contains the point (tracker frame); -1 if none.
  // Crystal::position() is the centre of the front face: the centre is half a crystal length further along the axis.
  int MC2StopTruth::crystalContaining(const Calorimeter& calorimeter, const CLHEP::Hep3Vector& pointInTracker) const {
    for (int c = 0; c < (int)calorimeter.nCrystals(); ++c) {
      const auto& diskInfo = calorimeter.disk(calorimeter.crystal(c).diskID()).diskInfo();
      const auto centre = calorimeter.mu2eToTracker(calorimeter.crystal(c).position()
                                                    + 0.5 * diskInfo.crystalZLength() * diskInfo.crystalDirection());
      if (std::fabs(pointInTracker.x() - centre.x()) < halfXY_ && std::fabs(pointInTracker.y() - centre.y()) < halfXY_
          && std::fabs(pointInTracker.z() - centre.z()) < halfZ_) {
        return c;
      }
    }
    return -1;
  }

  // Event, muon and stop-crystal branches.
  void MC2StopTruth::fillMuon(const art::Event& event, const SimParticle& muon, const CLHEP::Hep3Vector& endInTracker,
                              const CaloShowerStep* lastStep, int crystal, int disk,
                              const CLHEP::Hep3Vector& crystalFront, double muonEnergy, size_t nMuonCrystals) {
    dtEnd_ = muon.endGlobalTime() - lastStep->time();
    endX_ = endInTracker.x();
    endY_ = endInTracker.y();
    endZ_ = endInTracker.z();
    endKE_ = muon.endKineticEnergy();
    endCode_ = muon.stoppingCode().id();
    run_ = event.run();
    subrun_ = event.subRun();
    event_ = event.event();
    q_ = muon.pdgId() > 0 ? -1 : +1;   // pdg 13 = mu-
    cry_ = crystal;
    disk_ = disk;
    x_ = crystalFront.x();
    y_ = crystalFront.y();
    tStop_ = muon.endGlobalTime();
    eMu_ = muonEnergy;
    nCryMu_ = nMuonCrystals;
  }

  // First-generation daughters of the muon: pdg, creation code, start kinetic energy, delay from the stop (tStop_).
  void MC2StopTruth::fillDaughters(const SimParticleCollection& sims, const cet::map_vector_key& muonKey) {
    dPdg_.clear();
    dCode_.clear();
    dKE_.clear();
    dT_.clear();
    tDauMin_ = 1e9;
    for (const auto& [key, daughter] : sims) {
      if (daughter.parent().isNull() || daughter.parent().key() != muonKey.asUint()) {
        continue;
      }
      dPdg_.push_back(daughter.pdgId());
      dCode_.push_back(daughter.creationCode().id());
      const double mass = daughter.startMomentum().m();
      dKE_.push_back(daughter.startMomentum().e() - mass);
      const float delay = daughter.startGlobalTime() - tStop_;
      dT_.push_back(delay);
      tDauMin_ = std::min(tDauMin_, delay);
      daughterCodeCount_[q_ < 0 ? 0 : 1][daughter.creationCode().id()]++;
    }
    nDau_ = dPdg_.size();
    ++nStopped_[q_ < 0 ? 0 : 1];
  }

  void MC2StopTruth::resetDeposits() {
    eTot_ = eStop_ = eNear_ = eFar_ = 0;
    tDep1_ = 1e9;
    for (int w = 0; w < 4; ++w) {
      eWin_[w] = eStopWin_[w] = eNearWin_[w] = eFarWin_[w] = 0;
    }
    ePrompt_ = 0;
  }

  // One delayed step: total, time window, and stop crystal / neighbour (same disk, front-face centres within dNear) /
  // far. energyTimeSum accumulates energy * delay for tDepMean.
  void MC2StopTruth::addDeposit(const Calorimeter& calorimeter, const CaloShowerStep& step, int stopCrystal,
                                int stopDisk, const CLHEP::Hep3Vector& stopCrystalFront, double& energyTimeSum) {
    const double energy = step.energyDepG4(), delay = step.time() - tStop_;
    const int crystal = step.volumeG4ID();
    if (crystal < 0 || crystal >= (int)calorimeter.nCrystals()) {
      return;
    }
    eTot_ += energy;
    energyTimeSum += energy * delay;
    if (delay < tDep1_ && energy > 0.01) {
      tDep1_ = delay;
    }
    const int w = delay < 50 ? 0 : (delay < 500 ? 1 : (delay < 20000 ? 2 : 3));
    eWin_[w] += energy;
    if (crystal == stopCrystal) {
      eStop_ += energy;
      eStopWin_[w] += energy;
      return;
    }
    const auto crystalFront = calorimeter.mu2eToTracker(calorimeter.crystal(crystal).position());
    if (int(calorimeter.crystal(crystal).diskID()) == stopDisk
        && std::hypot(crystalFront.x() - stopCrystalFront.x(), crystalFront.y() - stopCrystalFront.y()) <= dNear_) {
      eNear_ += energy;
      eNearWin_[w] += energy;
    } else {
      eFar_ += energy;
      eFarWin_[w] += energy;
    }
  }

  // True if the muon is in the parent chain of the particle (at most 40 generations up).
  bool MC2StopTruth::isDescendant(const SimParticle& particle, const cet::map_vector_key& muonKey) {
    art::Ptr<SimParticle> parent = particle.parent();
    for (int depth = 0; depth < 40 && parent.isNonnull(); ++depth) {
      if (parent.key() == muonKey.asUint()) {
        return true;
      }
      parent = parent->parent();
    }
    return false;
  }

  void MC2StopTruth::endJob() {
    printf("MC2StopTruth: stopped mu- %ld, mu+ %ld\n", nStopped_[0], nStopped_[1]);
    for (int k = 0; k < 2; ++k) {
      for (auto& [code, count] : daughterCodeCount_[k]) {
        printf("  %s daughters created by %-28s %ld\n", k == 0 ? "mu-" : "mu+",
               ProcessCode(static_cast<ProcessCode::enum_type>(code)).name().c_str(), count);
      }
    }
  }

}

DEFINE_ART_MODULE(mu2e::MC2StopTruth)
