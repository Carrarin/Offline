//
// MC2MuXrayAna: muonic X-rays of the mu- that stop in the calorimeter, seen by the LYSO crystals (CAPHRI), and the
// prompt energy around the stop point in the CsI. Made for the Geant4 gun simulation of mu- at rest.
//
// Inputs (labels configurable: gun job g4run / CaloShowerStepMaker / MC2CaloHitMaker, hall production
// compressDetStepMCs):
//   SimParticleCollection     (simParticles)
//   CaloShowerStepCollection  (caloShowerSteps)
//   CaloHitCollection         (caloHits)
//
// Steps, for every stopped mu- (pdg 13, end kinetic energy <= maxEndKE, stopped by muMinusCaptureAtRest):
//
//   1. stop point: end point of the muon (tracker frame); stop crystal = the crystal box (34 x 34 x 200 mm around its
//      centre) that contains it (-1 if none: stopped in the material around the crystals); nearest LYSO and its
//      distance in x-y. Crystal::position() is the centre of the FRONT face: the centre of the crystal is half a
//      crystal length further along the crystal direction.
//   2. first-generation photons of the muon: EMCascade (atomic cascade of the muonic atom: the X-rays) and
//      NuclearCapture, with energy and start time - stop time
//   3. prompt Geant4 energy (|t - stop| < promptWin) in the stop crystal and in its CsI neighbours (centre within
//      nbDist mm in x-y, same disk): sum and the largest single crystal (the muon's own steps included)
//   4. Geant4 energy in each LYSO from all the CaloShowerSteps, prompt and total
//   5. reco: CaloHits of the LYSO crystals
//
// Output (TFileService), tree "xray", one entry per stopped mu-:
//   run, subrun, event
//   endX, endY, endZ, tStop     stop point (tracker frame, mm) and time (ns)
//   stopCry                     stop crystal (-1: none)
//   nearLyso, dLyso             nearest LYSO crystal and its distance in x-y (mm)
//   xE, xT                      EMCascade photons: energy (MeV), start time - tStop (ns)
//   nE, nT                      NuclearCapture photons: energy, start time - tStop
//   eLyPrompt[4], eLyAll[4]     Geant4 energy in each LYSO (order of lysoCrystals): prompt and total
//   eStopPrompt                 prompt Geant4 energy in the stop crystal
//   eNbPrompt, eNbMax, nbMaxCry prompt Geant4 energy in the CsI neighbours: sum, largest crystal and its id
//   hCry, hE, hT                CaloHits of the LYSO: crystal, energy, time - tStop
//
// Main parameters: lysoCrystals, maxEndKE, promptWin, nbDist.
//
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art_root_io/TFileService.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"

#include "Offline/CalorimeterGeom/inc/Calorimeter.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/MCDataProducts/inc/CaloShowerStep.hh"
#include "Offline/MCDataProducts/inc/ProcessCode.hh"
#include "Offline/MCDataProducts/inc/SimParticle.hh"
#include "Offline/RecoDataProducts/inc/CaloHit.hh"
#include "CLHEP/Vector/ThreeVector.h"
#include "TTree.h"

#include <cmath>
#include <vector>

namespace mu2e {

  class MC2MuXrayAna : public art::EDAnalyzer {
    public:
      struct Config {
        using Name = fhicl::Name;
        using Comment = fhicl::Comment;

        // Inputs
        fhicl::Atom<art::InputTag> simTag{Name("simParticles"), Comment("SimParticle collection"),
            "compressDetStepMCs"};
        fhicl::Atom<art::InputTag> stepTag{Name("caloShowerSteps"), Comment("CaloShowerStep collection"),
            "compressDetStepMCs"};
        fhicl::Atom<art::InputTag> hitTag{Name("caloHits"), Comment("CaloHit collection"), "MC2CaloHitMaker"};

        // Crystals
        fhicl::Sequence<int> lyso{Name("lysoCrystals"), Comment("crystal ids of the LYSO (CAPHRI)"),
            std::vector<int>{582, 609, 610, 637}};

        // Selection of the stopped muons and energy sums
        fhicl::Atom<double> maxEndKE{Name("maxEndKE"), Comment("max end kinetic energy of a stopped muon [MeV]"), 1.};
        fhicl::Atom<double> promptWin{Name("promptWin"),
            Comment("prompt window |t - stop| [ns] for the Geant4 energy in the LYSO"), 10.};
        fhicl::Atom<double> nbDist{Name("nbDist"),
            Comment("CsI neighbours of the stop crystal: centre within [mm] in x-y"), 50.};
      };
      using Parameters = art::EDAnalyzer::Table<Config>;

      explicit MC2MuXrayAna(const Parameters& config);
      void beginJob() override;
      void analyze(const art::Event& event) override;

    private:
      bool isStoppedMuMinus(const SimParticle& particle) const;
      int lysoIndex(int crystalId) const;
      void findStopCrystal(const Calorimeter& cal, const CLHEP::Hep3Vector& endPosition);
      void collectPhotons(const SimParticleCollection& particles, const cet::map_vector_key& muonKey);
      void sumPromptEnergyAroundStop(const CaloShowerStepCollection& steps, const Calorimeter& cal);
      void sumLysoEnergy(const CaloShowerStepCollection& steps);
      void collectLysoHits(const CaloHitCollection& hits);

      // Inputs
      art::InputTag simTag_;
      art::InputTag stepTag_;
      art::InputTag hitTag_;

      // Parameters
      std::vector<int> lysoCrystals_;
      double maxEndKE_;
      double promptWindow_;
      double neighbourDistance_;
      // Half sizes of the crystal box around its centre (mm)
      static constexpr double halfXY_ = 17.;
      static constexpr double halfZ_ = 100.;

      // Output tree and its variables
      TTree* tree_ = nullptr;
      int run_;
      int subRun_;
      int event_;
      // stop point
      float endX_;
      float endY_;
      float endZ_;
      float tStop_;
      int stopCrystal_;
      int nearestLyso_;
      float dLyso_;
      // first-generation photons
      std::vector<float> xrayE_;
      std::vector<float> xrayT_;
      std::vector<float> captureE_;
      std::vector<float> captureT_;
      // Geant4 energies
      float eLysoPrompt_[4];
      float eLysoAll_[4];
      float eStopPrompt_;
      float eNeighboursPrompt_;
      float eNeighbourMax_;
      int neighbourMaxCrystal_;
      // reco hits of the LYSO
      std::vector<int> hitCrystal_;
      std::vector<float> hitE_;
      std::vector<float> hitT_;
  };

  MC2MuXrayAna::MC2MuXrayAna(const Parameters& config) :
    art::EDAnalyzer(config),
    simTag_(config().simTag()),
    stepTag_(config().stepTag()),
    hitTag_(config().hitTag()),
    lysoCrystals_(config().lyso()),
    maxEndKE_(config().maxEndKE()),
    promptWindow_(config().promptWin()),
    neighbourDistance_(config().nbDist())
  {
    consumes<SimParticleCollection>(simTag_);
    consumes<CaloShowerStepCollection>(stepTag_);
    consumes<CaloHitCollection>(hitTag_);
  }

  void MC2MuXrayAna::beginJob() {
    art::ServiceHandle<art::TFileService> tfs;
    tree_ = tfs->make<TTree>("xray", "stopped mu-: muonic X-rays (truth) and LYSO energies (truth and reco)");
    tree_->Branch("run", &run_);
    tree_->Branch("subrun", &subRun_);
    tree_->Branch("event", &event_);
    tree_->Branch("endX", &endX_);
    tree_->Branch("endY", &endY_);
    tree_->Branch("endZ", &endZ_);
    tree_->Branch("tStop", &tStop_);
    tree_->Branch("stopCry", &stopCrystal_);
    tree_->Branch("nearLyso", &nearestLyso_);
    tree_->Branch("dLyso", &dLyso_);
    tree_->Branch("xE", &xrayE_);
    tree_->Branch("xT", &xrayT_);
    tree_->Branch("nE", &captureE_);
    tree_->Branch("nT", &captureT_);
    tree_->Branch("eLyPrompt", eLysoPrompt_, "eLyPrompt[4]/F");
    tree_->Branch("eLyAll", eLysoAll_, "eLyAll[4]/F");
    tree_->Branch("eStopPrompt", &eStopPrompt_);
    tree_->Branch("eNbPrompt", &eNeighboursPrompt_);
    tree_->Branch("eNbMax", &eNeighbourMax_);
    tree_->Branch("nbMaxCry", &neighbourMaxCrystal_);
    tree_->Branch("hCry", &hitCrystal_);
    tree_->Branch("hE", &hitE_);
    tree_->Branch("hT", &hitT_);
  }

  bool MC2MuXrayAna::isStoppedMuMinus(const SimParticle& particle) const {
    const bool rejected = particle.pdgId() != 13 || particle.endKineticEnergy() > maxEndKE_
                       || particle.stoppingCode().id() != ProcessCode::muMinusCaptureAtRest;
    return !rejected;
  }

  // Position of the crystal in lysoCrystals, -1 if it is not a LYSO.
  int MC2MuXrayAna::lysoIndex(int crystalId) const {
    for (size_t i = 0; i < lysoCrystals_.size(); ++i) {
      if (lysoCrystals_[i] == crystalId) {
        return int(i);
      }
    }
    return -1;
  }

  // Step 1: crystal whose box contains the stop point (the last one found, -1 if none) and nearest LYSO in x-y
  // (with the stop point in double precision, not the float tree variables).
  void MC2MuXrayAna::findStopCrystal(const Calorimeter& cal, const CLHEP::Hep3Vector& endPosition) {
    stopCrystal_ = -1;
    nearestLyso_ = -1;
    dLyso_ = 1e9;
    for (int c = 0; c < (int)cal.nCrystals(); ++c) {
      const auto& crystal = cal.crystal(c);
      const auto& diskInfo = cal.disk(crystal.diskID()).diskInfo();
      // centre of the crystal: position() is the centre of the front face
      const auto centre = cal.mu2eToTracker(crystal.position()
                                            + 0.5 * diskInfo.crystalZLength() * diskInfo.crystalDirection());
      if (std::fabs(endPosition.x() - centre.x()) < halfXY_ && std::fabs(endPosition.y() - centre.y()) < halfXY_
          && std::fabs(endPosition.z() - centre.z()) < halfZ_) {
        stopCrystal_ = c;
      }
      if (lysoIndex(c) >= 0) {
        const float distance = std::hypot(endPosition.x() - centre.x(), endPosition.y() - centre.y());
        if (distance < dLyso_) {
          dLyso_ = distance;
          nearestLyso_ = c;
        }
      }
    }
  }

  // Step 2: photons whose parent is the muon, by creation process.
  void MC2MuXrayAna::collectPhotons(const SimParticleCollection& particles, const cet::map_vector_key& muonKey) {
    xrayE_.clear();
    xrayT_.clear();
    captureE_.clear();
    captureT_.clear();
    for (const auto& [key, particle] : particles) {
      if (particle.pdgId() != 22 || particle.parent().isNull() || particle.parent().key() != muonKey.asUint()) {
        continue;
      }
      const float energy = particle.startMomentum().e();
      const float dt = particle.startGlobalTime() - tStop_;
      if (particle.creationCode().id() == ProcessCode::EMCascade) {
        xrayE_.push_back(energy);
        xrayT_.push_back(dt);
      }
      if (particle.creationCode().id() == ProcessCode::NuclearCapture) {
        captureE_.push_back(energy);
        captureT_.push_back(dt);
      }
    }
  }

  // Step 3: prompt Geant4 energy in the stop crystal and in the CsI crystals of the same disk whose front-face centre
  // is within neighbourDistance in x-y of that of the stop crystal (LYSO left out); the largest neighbour.
  void MC2MuXrayAna::sumPromptEnergyAroundStop(const CaloShowerStepCollection& steps, const Calorimeter& cal) {
    eStopPrompt_ = eNeighboursPrompt_ = eNeighbourMax_ = 0;
    neighbourMaxCrystal_ = -1;
    if (stopCrystal_ < 0) {
      return;
    }

    std::vector<float> eNeighbour(cal.nCrystals(), 0.f);
    for (const auto& step : steps) {
      const int c = step.volumeG4ID();
      if (c < 0 || c >= (int)cal.nCrystals() || lysoIndex(c) >= 0 || std::fabs(step.time() - tStop_) >= promptWindow_) {
        continue;
      }
      if (c == stopCrystal_) {
        eStopPrompt_ += step.energyDepG4();
        continue;
      }
      if (cal.crystal(c).diskID() != cal.crystal(stopCrystal_).diskID()) {
        continue;
      }
      const auto position = cal.mu2eToTracker(cal.crystal(c).position());
      const auto stopPosition = cal.mu2eToTracker(cal.crystal(stopCrystal_).position());
      if (std::hypot(position.x() - stopPosition.x(), position.y() - stopPosition.y()) <= neighbourDistance_) {
        eNeighboursPrompt_ += step.energyDepG4();
        eNeighbour[c] += step.energyDepG4();
      }
    }

    for (int c = 0; c < (int)cal.nCrystals(); ++c) {
      if (eNeighbour[c] > eNeighbourMax_) {
        eNeighbourMax_ = eNeighbour[c];
        neighbourMaxCrystal_ = c;
      }
    }
  }

  // Step 4: Geant4 energy in each LYSO, all the steps and the prompt ones.
  void MC2MuXrayAna::sumLysoEnergy(const CaloShowerStepCollection& steps) {
    for (int i = 0; i < 4; ++i) {
      eLysoPrompt_[i] = eLysoAll_[i] = 0;
    }
    for (const auto& step : steps) {
      const int iLyso = lysoIndex(step.volumeG4ID());
      if (iLyso < 0 || iLyso > 3) {
        continue;
      }
      eLysoAll_[iLyso] += step.energyDepG4();
      if (std::fabs(step.time() - tStop_) < promptWindow_) {
        eLysoPrompt_[iLyso] += step.energyDepG4();
      }
    }
  }

  // Step 5: reco hits of the LYSO crystals, time relative to the stop.
  void MC2MuXrayAna::collectLysoHits(const CaloHitCollection& hits) {
    hitCrystal_.clear();
    hitE_.clear();
    hitT_.clear();
    for (const auto& hit : hits) {
      if (lysoIndex(hit.crystalID()) >= 0) {
        hitCrystal_.push_back(hit.crystalID());
        hitE_.push_back(hit.energyDep());
        hitT_.push_back(hit.time() - tStop_);
      }
    }
  }

  void MC2MuXrayAna::analyze(const art::Event& event) {
    const Calorimeter& cal = *GeomHandle<Calorimeter>();
    const auto& particles = *event.getValidHandle<SimParticleCollection>(simTag_);
    const auto& steps = *event.getValidHandle<CaloShowerStepCollection>(stepTag_);
    const auto& hits = *event.getValidHandle<CaloHitCollection>(hitTag_);

    for (const auto& [key, muon] : particles) {
      if (!isStoppedMuMinus(muon)) {
        continue;
      }
      const auto endPosition = cal.mu2eToTracker(muon.endPosition());
      run_ = event.run();
      subRun_ = event.subRun();
      event_ = event.event();
      endX_ = endPosition.x();
      endY_ = endPosition.y();
      endZ_ = endPosition.z();
      tStop_ = muon.endGlobalTime();

      findStopCrystal(cal, endPosition);
      collectPhotons(particles, key);
      sumPromptEnergyAroundStop(steps, cal);
      sumLysoEnergy(steps);
      collectLysoHits(hits);
      tree_->Fill();
    }
  }

}

DEFINE_ART_MODULE(mu2e::MC2MuXrayAna)
