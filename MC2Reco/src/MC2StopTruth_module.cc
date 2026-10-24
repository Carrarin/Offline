// MC2StopTruth: Monte Carlo truth of the muons that stop in the calorimeter crystals (hall production, Mu2eG4): what the decay or
// the capture leaves, and when. One tree entry per stopped mu+ or mu-.
//
// Stopped muon: SimParticle with |pdg| = 13, end kinetic energy < maxEndKE, with CaloShowerSteps of its own, whose end point is inside
//               a crystal box (34 x 34 x 200 mm around the crystal centre = front face (Crystal::position) + half length): inCry = 1 and that is the stop crystal; otherwise (stopped
//               in the material around the crystals) inCry = 0 and the reference crystal is the last one the muon crossed. Its own deposit (all its steps) =
//               prompt. (Most muons with calo steps that reach zero energy stop outside the crystals: floor, supports.)
// Daughters   : first generation (parent = the muon): pdg, creation code, kinetic energy at start, start time - muon end time.
//               mu- in Geant4 11 is first captured in an atom (muonic atom), which then decays (DIO electron) or is captured by the
//               nucleus: the delay of the daughters measures that lifetime.
// Deposits    : CaloShowerSteps of every descendant (parent chain up to the muon) AND the muon's own steps later than its end + 1 ns
//               (CaloShowerStepMaker gives the energy of particles born inside the calorimeter to the particle that entered it), by crystal: stop crystal, neighbours (centre within
//               dNear mm, same disk), far; energy (G4, MeV), first time and energy-weighted time after the stop; also in time windows.
// Inputs      : SimParticles and CaloShowerSteps as saved by the Primary job (label compressDetStepMCs).
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art_root_io/TFileService.h"
#include "fhiclcpp/types/Atom.h"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/CalorimeterGeom/inc/Calorimeter.hh"
#include "Offline/MCDataProducts/inc/SimParticle.hh"
#include "Offline/MCDataProducts/inc/CaloShowerStep.hh"
#include "Offline/MCDataProducts/inc/ProcessCode.hh"
#include "TTree.h"
#include <cmath>
#include <cstdio>
#include <map>
#include <set>
#include <vector>

namespace mu2e {
  class MC2StopTruth : public art::EDAnalyzer {
    public:
      struct Config {
        using Name = fhicl::Name; using Comment = fhicl::Comment;
        fhicl::Atom<art::InputTag> simTag   { Name("simParticles"), Comment("SimParticle collection"), "compressDetStepMCs" };
        fhicl::Atom<art::InputTag> stepTag  { Name("caloShowerSteps"), Comment("CaloShowerStep collection"), "compressDetStepMCs" };
        fhicl::Atom<double> maxEndKE { Name("maxEndKE"), Comment("max end kinetic energy of a stopped muon [MeV]"), 1. };
        fhicl::Atom<double> stopDt   { Name("stopDt"), Comment("not used (kept for old configurations)"), 5. };
        fhicl::Atom<double> dNear    { Name("dNear"), Comment("neighbour crystal: centre within [mm] of the stop crystal"), 80. };
      };
      using Parameters = art::EDAnalyzer::Table<Config>;
      explicit MC2StopTruth(const Parameters& c);
      void beginJob() override;
      void analyze(const art::Event& e) override;
      void endJob() override;
    private:
      art::InputTag simTag_, stepTag_; double maxKE_, stopDt_, dNear_, halfXY_ = 17., halfZ_ = 100.;
      TTree* t_ = nullptr;
      int run_, sub_, evt_, q_, cry_, disk_, nCryMu_, nDau_, endCode_, inCry_;
      float eNearWin_[4], eFarWin_[4], x_, y_, ePrompt_, dtEnd_, endX_, endY_, endZ_, endKE_, tStop_, eMu_, tDauMin_, eTot_, eStop_, eNear_, eFar_, tDep1_, tDepMean_, eWin_[4], eStopWin_[4];
      std::vector<int> dPdg_, dCode_; std::vector<float> dKE_, dT_;
      std::map<int, long> codeCount_[2]; long nStop_[2] = {0, 0};
  };

  MC2StopTruth::MC2StopTruth(const Parameters& c) : art::EDAnalyzer(c), simTag_(c().simTag()), stepTag_(c().stepTag()), maxKE_(c().maxEndKE()),
    stopDt_(c().stopDt()), dNear_(c().dNear()) {
    consumes<SimParticleCollection>(simTag_); consumes<CaloShowerStepCollection>(stepTag_);
  }

  void MC2StopTruth::beginJob() {
    art::ServiceHandle<art::TFileService> tfs;
    t_ = tfs->make<TTree>("stops", "muons stopped in the crystals: truth of the decay / capture products");
    t_->Branch("run", &run_); t_->Branch("subrun", &sub_); t_->Branch("event", &evt_); t_->Branch("q", &q_); t_->Branch("cry", &cry_); t_->Branch("disk", &disk_);
    t_->Branch("inCry", &inCry_); t_->Branch("ePrompt", &ePrompt_); t_->Branch("dtEnd", &dtEnd_); t_->Branch("endCode", &endCode_); t_->Branch("endX", &endX_); t_->Branch("endY", &endY_); t_->Branch("endZ", &endZ_); t_->Branch("endKE", &endKE_); t_->Branch("x", &x_); t_->Branch("y", &y_); t_->Branch("tStop", &tStop_); t_->Branch("eMu", &eMu_); t_->Branch("nCryMu", &nCryMu_);
    t_->Branch("nDau", &nDau_); t_->Branch("tDauMin", &tDauMin_); t_->Branch("dPdg", &dPdg_); t_->Branch("dCode", &dCode_); t_->Branch("dKE", &dKE_); t_->Branch("dT", &dT_);
    t_->Branch("eTot", &eTot_); t_->Branch("eStop", &eStop_); t_->Branch("eNear", &eNear_); t_->Branch("eFar", &eFar_); t_->Branch("tDep1", &tDep1_); t_->Branch("tDepMean", &tDepMean_);
    t_->Branch("eWin", eWin_, "eWin[4]/F"); t_->Branch("eStopWin", eStopWin_, "eStopWin[4]/F"); t_->Branch("eNearWin", eNearWin_, "eNearWin[4]/F"); t_->Branch("eFarWin", eFarWin_, "eFarWin[4]/F");   // windows after the stop: 0-50 ns, 50-500 ns, 0.5-20 us, > 20 us
  }

  void MC2StopTruth::analyze(const art::Event& e) {
    const Calorimeter& cal = *GeomHandle<Calorimeter>();
    const auto simH = e.getValidHandle<SimParticleCollection>(simTag_); const auto& sims = *simH;
    const auto& steps = *e.getValidHandle<CaloShowerStepCollection>(stepTag_);
    // steps of each SimParticle
    std::map<unsigned long, std::vector<const CaloShowerStep*>> byPart;
    for (const auto& s : steps) if (s.simParticle().isNonnull()) byPart[s.simParticle().key()].push_back(&s);
    for (const auto& [key, sp] : sims) {
      if (std::abs(sp.pdgId()) != 13 || sp.endKineticEnergy() > maxKE_) continue;
      const auto it = byPart.find(key.asUint()); if (it == byPart.end()) continue;
      double emu = 0; std::set<int> cmu; const CaloShowerStep* last = nullptr;
      for (auto* s : it->second) { emu += s->energyDepG4(); cmu.insert(s->volumeG4ID()); if (!last || s->time() > last->time()) last = s; }
      // stop crystal: the crystal box (34 x 34 x 200 mm around its centre) that contains the end point
      const auto pe = cal.mu2eToTracker(sp.endPosition()); int c0 = -1;
      for (int c = 0; c < (int)cal.nCrystals(); ++c) {
        const auto pc = cal.mu2eToTracker(cal.crystal(c).position() + 0.5 * cal.disk(cal.crystal(c).diskID()).diskInfo().crystalZLength() * cal.disk(cal.crystal(c).diskID()).diskInfo().crystalDirection());   // centre: position() is the front face
        if (std::fabs(pe.x() - pc.x()) < halfXY_ && std::fabs(pe.y() - pc.y()) < halfXY_ && std::fabs(pe.z() - pc.z()) < halfZ_) { c0 = c; break; }
      }
      inCry_ = c0 >= 0;
      if (c0 < 0) c0 = last->volumeG4ID();   // stopped outside the crystals: reference = the last crystal the muon crossed
      if (c0 < 0 || c0 >= (int)cal.nCrystals()) continue;
      dtEnd_ = sp.endGlobalTime() - last->time(); endX_ = pe.x(); endY_ = pe.y(); endZ_ = pe.z(); endKE_ = sp.endKineticEnergy(); endCode_ = sp.stoppingCode().id();
      const auto p0 = cal.mu2eToTracker(cal.crystal(c0).position());
      const int d0 = int(cal.crystal(c0).diskID());
      run_ = e.run(); sub_ = e.subRun(); evt_ = e.event(); q_ = sp.pdgId() > 0 ? -1 : +1;   // pdg 13 = mu-
      cry_ = c0; disk_ = d0; x_ = p0.x(); y_ = p0.y(); tStop_ = sp.endGlobalTime(); eMu_ = emu; nCryMu_ = cmu.size();
      // first-generation daughters
      dPdg_.clear(); dCode_.clear(); dKE_.clear(); dT_.clear(); tDauMin_ = 1e9;
      for (const auto& [k2, s2] : sims) {
        if (s2.parent().isNull() || s2.parent().key() != key.asUint()) continue;
        dPdg_.push_back(s2.pdgId()); dCode_.push_back(s2.creationCode().id());
        const double m = s2.startMomentum().m(); dKE_.push_back(s2.startMomentum().e() - m);
        const float dt = s2.startGlobalTime() - tStop_; dT_.push_back(dt); tDauMin_ = std::min(tDauMin_, dt);
        codeCount_[q_ < 0 ? 0 : 1][s2.creationCode().id()]++;
      }
      nDau_ = dPdg_.size(); ++nStop_[q_ < 0 ? 0 : 1];
      // deposits of all descendants
      eTot_ = eStop_ = eNear_ = eFar_ = 0; tDep1_ = 1e9; double tw = 0; for (int w = 0; w < 4; ++w) eWin_[w] = eStopWin_[w] = eNearWin_[w] = eFarWin_[w] = 0;
      ePrompt_ = 0;
      auto addDep = [&](const CaloShowerStep* s) {
          const double ed = s->energyDepG4(), dt = s->time() - tStop_; const int c = s->volumeG4ID();
          if (c < 0 || c >= (int)cal.nCrystals()) return;
          eTot_ += ed; tw += ed * dt; if (dt < tDep1_ && ed > 0.01) tDep1_ = dt;
          const int w = dt < 50 ? 0 : (dt < 500 ? 1 : (dt < 20000 ? 2 : 3)); eWin_[w] += ed;
          if (c == c0) { eStop_ += ed; eStopWin_[w] += ed; }
          else { const auto pc = cal.mu2eToTracker(cal.crystal(c).position());
            if (int(cal.crystal(c).diskID()) == d0 && std::hypot(pc.x() - p0.x(), pc.y() - p0.y()) <= dNear_) { eNear_ += ed; eNearWin_[w] += ed; } else { eFar_ += ed; eFarWin_[w] += ed; } }
      };
      for (auto* s : it->second) { if (s->time() > tStop_ + 1.) addDep(s); else ePrompt_ += s->energyDepG4(); }
      for (const auto& [k2, s2] : sims) {
        if (k2 == key) continue;
        bool desc = false; art::Ptr<SimParticle> p = s2.parent();
        for (int depth = 0; depth < 40 && p.isNonnull(); ++depth) { if (p.key() == key.asUint()) { desc = true; break; } p = p->parent(); }
        if (!desc) continue;
        const auto jt = byPart.find(k2.asUint()); if (jt == byPart.end()) continue;
        for (auto* s : jt->second) addDep(s);
      }
      tDepMean_ = eTot_ > 0 ? tw / eTot_ : -1;
      t_->Fill();
    }
  }

  void MC2StopTruth::endJob() {
    printf("MC2StopTruth: stopped mu- %ld, mu+ %ld\n", nStop_[0], nStop_[1]);
    for (int k = 0; k < 2; ++k) for (auto& [c, n] : codeCount_[k])
      printf("  %s daughters created by %-28s %ld\n", k == 0 ? "mu-" : "mu+", ProcessCode(static_cast<ProcessCode::enum_type>(c)).name().c_str(), n);
  }
}
DEFINE_ART_MODULE(mu2e::MC2StopTruth)
