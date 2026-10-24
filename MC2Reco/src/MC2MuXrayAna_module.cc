// MC2MuXrayAna: muonic X-rays of the mu- that stop in the calorimeter, seen by the LYSO crystals (CAPHRI). One tree entry per
// stopped mu- (pdg 13, end kinetic energy < maxEndKE, stopped by muMinusCaptureAtRest).
//
// Stop        : end point of the muon (tracker frame); stop crystal = the crystal box (34 x 34 x 200 mm around its centre = front face + half length) that contains
//               it (-1 if none: stopped in the material around the crystals); nearest LYSO and its distance in x-y.
// Truth       : first-generation photons of the muon: EMCascade (atomic cascade of the muonic atom: X-rays) and NuclearCapture, with
//               energy and start time - stop time; Geant4 energy in each LYSO from all CaloShowerSteps, prompt (|t - stop| < promptWin)
//               and total; prompt Geant4 energy in the stop crystal and in the CsI neighbours (centre within nbDist mm, same disk):
//               sum and the largest single crystal (the muon's own steps included).
// Reco        : CaloHits of the LYSO crystals (crystal, energy, time - stop time).
// Inputs      : SimParticles, CaloShowerSteps, CaloHits (labels configurable: gun job g4run / CaloShowerStepMaker / MC2CaloHitMaker,
//               hall production compressDetStepMCs).
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art_root_io/TFileService.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/CalorimeterGeom/inc/Calorimeter.hh"
#include "Offline/MCDataProducts/inc/SimParticle.hh"
#include "Offline/MCDataProducts/inc/CaloShowerStep.hh"
#include "Offline/MCDataProducts/inc/ProcessCode.hh"
#include "Offline/RecoDataProducts/inc/CaloHit.hh"
#include "TTree.h"
#include <cmath>
#include <vector>

namespace mu2e {
  class MC2MuXrayAna : public art::EDAnalyzer {
    public:
      struct Config {
        using Name = fhicl::Name; using Comment = fhicl::Comment;
        fhicl::Atom<art::InputTag> simTag  { Name("simParticles"), Comment("SimParticle collection"), "compressDetStepMCs" };
        fhicl::Atom<art::InputTag> stepTag { Name("caloShowerSteps"), Comment("CaloShowerStep collection"), "compressDetStepMCs" };
        fhicl::Atom<art::InputTag> hitTag  { Name("caloHits"), Comment("CaloHit collection"), "MC2CaloHitMaker" };
        fhicl::Sequence<int> lyso          { Name("lysoCrystals"), Comment("crystal ids of the LYSO (CAPHRI)"), std::vector<int>{582, 609, 610, 637} };
        fhicl::Atom<double> maxEndKE       { Name("maxEndKE"), Comment("max end kinetic energy of a stopped muon [MeV]"), 1. };
        fhicl::Atom<double> promptWin      { Name("promptWin"), Comment("prompt window |t - stop| [ns] for the Geant4 energy in the LYSO"), 10. };
        fhicl::Atom<double> nbDist         { Name("nbDist"), Comment("CsI neighbours of the stop crystal: centre within [mm] in x-y"), 50. };
      };
      using Parameters = art::EDAnalyzer::Table<Config>;
      explicit MC2MuXrayAna(const Parameters& c);
      void beginJob() override;
      void analyze(const art::Event& e) override;
    private:
      art::InputTag simTag_, stepTag_, hitTag_; std::vector<int> lyso_; double maxKE_, prompt_, nbDist_, halfXY_ = 17., halfZ_ = 100.;
      float eStopPr_, eNbPr_, eNbMax_; int nbMaxCry_;
      TTree* t_ = nullptr;
      int run_, sub_, evt_, stopCry_, nearLyso_; float endX_, endY_, endZ_, tStop_, dLyso_;
      float eLyPrompt_[4], eLyAll_[4];
      std::vector<float> xE_, xT_, nE_, nT_, hE_, hT_; std::vector<int> hCry_;
  };

  MC2MuXrayAna::MC2MuXrayAna(const Parameters& c) : art::EDAnalyzer(c), simTag_(c().simTag()), stepTag_(c().stepTag()), hitTag_(c().hitTag()),
    lyso_(c().lyso()), maxKE_(c().maxEndKE()), prompt_(c().promptWin()), nbDist_(c().nbDist()) {
    consumes<SimParticleCollection>(simTag_); consumes<CaloShowerStepCollection>(stepTag_); consumes<CaloHitCollection>(hitTag_);
  }

  void MC2MuXrayAna::beginJob() {
    art::ServiceHandle<art::TFileService> tfs;
    t_ = tfs->make<TTree>("xray", "stopped mu-: muonic X-rays (truth) and LYSO energies (truth and reco)");
    t_->Branch("run", &run_); t_->Branch("subrun", &sub_); t_->Branch("event", &evt_);
    t_->Branch("endX", &endX_); t_->Branch("endY", &endY_); t_->Branch("endZ", &endZ_); t_->Branch("tStop", &tStop_);
    t_->Branch("stopCry", &stopCry_); t_->Branch("nearLyso", &nearLyso_); t_->Branch("dLyso", &dLyso_);
    t_->Branch("xE", &xE_); t_->Branch("xT", &xT_); t_->Branch("nE", &nE_); t_->Branch("nT", &nT_);
    t_->Branch("eLyPrompt", eLyPrompt_, "eLyPrompt[4]/F"); t_->Branch("eLyAll", eLyAll_, "eLyAll[4]/F");
    t_->Branch("eStopPrompt", &eStopPr_); t_->Branch("eNbPrompt", &eNbPr_); t_->Branch("eNbMax", &eNbMax_); t_->Branch("nbMaxCry", &nbMaxCry_);
    t_->Branch("hCry", &hCry_); t_->Branch("hE", &hE_); t_->Branch("hT", &hT_);
  }

  void MC2MuXrayAna::analyze(const art::Event& e) {
    const Calorimeter& cal = *GeomHandle<Calorimeter>();
    const auto& sims = *e.getValidHandle<SimParticleCollection>(simTag_);
    const auto& steps = *e.getValidHandle<CaloShowerStepCollection>(stepTag_);
    const auto& hits = *e.getValidHandle<CaloHitCollection>(hitTag_);
    auto lysoIndex = [&](int c) { for (size_t i = 0; i < lyso_.size(); ++i) if (lyso_[i] == c) return int(i); return -1; };
    for (const auto& [key, sp] : sims) {
      if (sp.pdgId() != 13 || sp.endKineticEnergy() > maxKE_ || sp.stoppingCode().id() != ProcessCode::muMinusCaptureAtRest) continue;
      const auto pe = cal.mu2eToTracker(sp.endPosition());
      run_ = e.run(); sub_ = e.subRun(); evt_ = e.event(); endX_ = pe.x(); endY_ = pe.y(); endZ_ = pe.z(); tStop_ = sp.endGlobalTime();
      stopCry_ = -1; nearLyso_ = -1; dLyso_ = 1e9;
      for (int c = 0; c < (int)cal.nCrystals(); ++c) {
        const auto pc = cal.mu2eToTracker(cal.crystal(c).position() + 0.5 * cal.disk(cal.crystal(c).diskID()).diskInfo().crystalZLength() * cal.disk(cal.crystal(c).diskID()).diskInfo().crystalDirection());   // centre: position() is the front face
        if (std::fabs(pe.x() - pc.x()) < halfXY_ && std::fabs(pe.y() - pc.y()) < halfXY_ && std::fabs(pe.z() - pc.z()) < halfZ_) stopCry_ = c;
        if (lysoIndex(c) >= 0) { const float d = std::hypot(pe.x() - pc.x(), pe.y() - pc.y()); if (d < dLyso_) { dLyso_ = d; nearLyso_ = c; } }
      }
      xE_.clear(); xT_.clear(); nE_.clear(); nT_.clear();
      for (const auto& [k2, s2] : sims) {
        if (s2.pdgId() != 22 || s2.parent().isNull() || s2.parent().key() != key.asUint()) continue;
        const float E = s2.startMomentum().e(), dt = s2.startGlobalTime() - tStop_;
        if (s2.creationCode().id() == ProcessCode::EMCascade) { xE_.push_back(E); xT_.push_back(dt); }
        if (s2.creationCode().id() == ProcessCode::NuclearCapture) { nE_.push_back(E); nT_.push_back(dt); }
      }
      for (int i = 0; i < 4; ++i) eLyPrompt_[i] = eLyAll_[i] = 0;
      eStopPr_ = eNbPr_ = eNbMax_ = 0; nbMaxCry_ = -1; std::vector<float> eNb(cal.nCrystals(), 0.f);
      for (const auto& s : steps) {
        const int c = s.volumeG4ID(); if (c < 0 || c >= (int)cal.nCrystals() || lysoIndex(c) >= 0 || stopCry_ < 0 || std::fabs(s.time() - tStop_) >= prompt_) continue;
        if (c == stopCry_) { eStopPr_ += s.energyDepG4(); continue; }
        if (cal.crystal(c).diskID() != cal.crystal(stopCry_).diskID()) continue;
        const auto p1 = cal.mu2eToTracker(cal.crystal(c).position()), p0 = cal.mu2eToTracker(cal.crystal(stopCry_).position());
        if (std::hypot(p1.x() - p0.x(), p1.y() - p0.y()) <= nbDist_) { eNbPr_ += s.energyDepG4(); eNb[c] += s.energyDepG4(); } }
      for (int c = 0; c < (int)cal.nCrystals(); ++c) if (eNb[c] > eNbMax_) { eNbMax_ = eNb[c]; nbMaxCry_ = c; }
      for (const auto& s : steps) {
        const int li = lysoIndex(s.volumeG4ID()); if (li < 0 || li > 3) continue;
        eLyAll_[li] += s.energyDepG4(); if (std::fabs(s.time() - tStop_) < prompt_) eLyPrompt_[li] += s.energyDepG4();
      }
      hCry_.clear(); hE_.clear(); hT_.clear();
      for (const auto& h : hits) if (lysoIndex(h.crystalID()) >= 0) { hCry_.push_back(h.crystalID()); hE_.push_back(h.energyDep()); hT_.push_back(h.time() - tStop_); }
      t_->Fill();
    }
  }
}
DEFINE_ART_MODULE(mu2e::MC2MuXrayAna)
