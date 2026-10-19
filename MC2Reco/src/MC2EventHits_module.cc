// MC2EventHits: all calorimeter hits of selected events, for event displays (topologies of stopped muons, data and Monte Carlo).
// Events: eventList = text file with "run subrun event" per line (empty: every event, up to maxEvents).
// Tree "hits", one entry per event: CaloHits (crystal, disk, x, y tracker frame, time [ns], energy [MeV], number of SiPMs) and the
// CaloRecoDigis of single SiPMs (SiPM id, crystal, time, energy) - the delayed deposits of MC2MuCapture use those.
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art_root_io/TFileService.h"
#include "fhiclcpp/types/Atom.h"
#include "Offline/ConfigTools/inc/ConfigFileLookupPolicy.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/CalorimeterGeom/inc/Calorimeter.hh"
#include "Offline/RecoDataProducts/inc/CaloHit.hh"
#include "Offline/RecoDataProducts/inc/CaloRecoDigi.hh"
#include "TTree.h"
#include <fstream>
#include <set>
#include <tuple>
#include <vector>

namespace mu2e {
  class MC2EventHits : public art::EDAnalyzer {
    public:
      struct Config {
        using Name = fhicl::Name; using Comment = fhicl::Comment;
        fhicl::Atom<art::InputTag> hitTag { Name("caloHitCollection"), Comment("CaloHits"), "MC2CaloHitMaker" };
        fhicl::Atom<art::InputTag> rdTag  { Name("caloRecoDigiCollection"), Comment("CaloRecoDigis"), "MC2CaloRecoDigiMaker" };
        fhicl::Atom<std::string> eventList { Name("eventList"), Comment("file with run subrun event per line; empty = all events"), "" };
        fhicl::Atom<int> maxEvents { Name("maxEvents"), Comment("max events written when eventList is empty"), 1000 };
      };
      using Parameters = art::EDAnalyzer::Table<Config>;
      explicit MC2EventHits(const Parameters& c);
      void beginJob() override;
      void analyze(const art::Event& e) override;
    private:
      art::InputTag hitTag_, rdTag_; std::set<std::tuple<int, int, int>> list_; bool all_; int maxEv_, nOut_ = 0;
      TTree* t_ = nullptr; int run_, sub_, evt_;
      std::vector<int> hc_, hd_, hn_, rs_, rc_; std::vector<float> hx_, hy_, ht_, he_, rt_, re_;
  };

  MC2EventHits::MC2EventHits(const Parameters& c) : art::EDAnalyzer(c), hitTag_(c().hitTag()), rdTag_(c().rdTag()), all_(c().eventList().empty()), maxEv_(c().maxEvents()) {
    consumes<CaloHitCollection>(hitTag_); consumes<CaloRecoDigiCollection>(rdTag_);
    if (!all_) { ConfigFileLookupPolicy resolve; std::ifstream in(resolve(c().eventList())); int r, s, ev; while (in >> r >> s >> ev) list_.insert({r, s, ev}); }
  }
  void MC2EventHits::beginJob() {
    art::ServiceHandle<art::TFileService> tfs; t_ = tfs->make<TTree>("hits", "calorimeter hits of selected events");
    t_->Branch("run", &run_); t_->Branch("subrun", &sub_); t_->Branch("event", &evt_);
    t_->Branch("cry", &hc_); t_->Branch("disk", &hd_); t_->Branch("nSiPM", &hn_); t_->Branch("x", &hx_); t_->Branch("y", &hy_); t_->Branch("t", &ht_); t_->Branch("e", &he_);
    t_->Branch("rdSiPM", &rs_); t_->Branch("rdCry", &rc_); t_->Branch("rdT", &rt_); t_->Branch("rdE", &re_);
  }
  void MC2EventHits::analyze(const art::Event& e) {
    run_ = e.run(); sub_ = e.subRun(); evt_ = e.event();
    if (all_ ? nOut_ >= maxEv_ : !list_.count({run_, sub_, evt_})) return;
    const Calorimeter& cal = *GeomHandle<Calorimeter>();
    hc_.clear(); hd_.clear(); hn_.clear(); hx_.clear(); hy_.clear(); ht_.clear(); he_.clear(); rs_.clear(); rc_.clear(); rt_.clear(); re_.clear();
    for (const auto& h : *e.getValidHandle<CaloHitCollection>(hitTag_)) {
      const auto& c = cal.crystal(h.crystalID()); const auto p = cal.mu2eToTracker(c.position());
      hc_.push_back(h.crystalID()); hd_.push_back(c.diskID()); hn_.push_back(h.nSiPMs()); hx_.push_back(p.x()); hy_.push_back(p.y()); ht_.push_back(h.time()); he_.push_back(h.energyDep());
    }
    for (const auto& r : *e.getValidHandle<CaloRecoDigiCollection>(rdTag_)) {
      const int id = r.SiPMID(); rs_.push_back(id); rc_.push_back(id / 2); rt_.push_back(r.time()); re_.push_back(r.energyDep());
    }
    t_->Fill(); ++nOut_;
  }
}
DEFINE_ART_MODULE(mu2e::MC2EventHits)
