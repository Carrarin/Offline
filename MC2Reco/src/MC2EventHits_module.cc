//
// MC2EventHits: all calorimeter hits of selected events, for event displays (topologies of stopped muons, data and
// Monte Carlo).
//
// Inputs: CaloHitCollection (caloHitCollection) and CaloRecoDigiCollection (caloRecoDigiCollection).
//
// Event selection: eventList = text file with "run subrun event" per line; with an empty eventList, every event up to
// maxEvents written.
//
// Output: TTree "hits" in the TFileService file, one entry per selected event:
//   run, subrun, event
//   CaloHits, one element per hit:      cry (crystal id), disk, nSiPM (number of SiPMs in the hit),
//                                       x, y (crystal position in the tracker frame [mm]), t [ns], e [MeV]
//   CaloRecoDigis of single SiPMs:      rdSiPM (SiPM id), rdCry (crystal = SiPM id / 2), rdT [ns], rdE [MeV]
//   (the delayed deposits of MC2MuCapture use the single-SiPM reco digis).
//
// Steps per event: selection, CaloHits, CaloRecoDigis, fill.
//
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
#include <string>
#include <tuple>
#include <vector>

namespace mu2e {

  class MC2EventHits : public art::EDAnalyzer {
    public:
      struct Config {
        using Name = fhicl::Name;
        using Comment = fhicl::Comment;

        // Inputs
        fhicl::Atom<art::InputTag> caloHitTag{Name("caloHitCollection"), Comment("CaloHits"), "MC2CaloHitMaker"};
        fhicl::Atom<art::InputTag> recoDigiTag{Name("caloRecoDigiCollection"), Comment("CaloRecoDigis"),
            "MC2CaloRecoDigiMaker"};

        // Event selection
        fhicl::Atom<std::string> eventList{Name("eventList"),
            Comment("file with run subrun event per line; empty = all events"), ""};
        fhicl::Atom<int> maxEvents{Name("maxEvents"), Comment("max events written when eventList is empty"), 1000};
      };
      using Parameters = art::EDAnalyzer::Table<Config>;

      explicit MC2EventHits(const Parameters& config);
      void beginJob() override;
      void analyze(const art::Event& event) override;

    private:
      void readEventList(const std::string& fileName);
      bool isSelected() const;
      void clearHits();
      void fillCaloHits(const art::Event& event);
      void fillRecoDigis(const art::Event& event);

      // Inputs
      art::InputTag caloHitTag_;
      art::InputTag recoDigiTag_;

      // Event selection
      std::set<std::tuple<int, int, int>> selectedEvents_;   // (run, subrun, event) from eventList
      bool allEvents_;                                       // no eventList: every event up to maxEvents
      int maxEvents_;
      int nWritten_ = 0;

      // Output tree and its branches
      TTree* tree_ = nullptr;
      int run_;
      int subRun_;
      int event_;
      std::vector<int> hitCrystal_;
      std::vector<int> hitDisk_;
      std::vector<int> hitNSiPMs_;
      std::vector<float> hitX_;
      std::vector<float> hitY_;
      std::vector<float> hitTime_;
      std::vector<float> hitEnergy_;
      std::vector<int> digiSiPM_;
      std::vector<int> digiCrystal_;
      std::vector<float> digiTime_;
      std::vector<float> digiEnergy_;
  };

  MC2EventHits::MC2EventHits(const Parameters& config) :
    art::EDAnalyzer(config),
    caloHitTag_(config().caloHitTag()),
    recoDigiTag_(config().recoDigiTag()),
    allEvents_(config().eventList().empty()),
    maxEvents_(config().maxEvents())
  {
    consumes<CaloHitCollection>(caloHitTag_);
    consumes<CaloRecoDigiCollection>(recoDigiTag_);
    if (!allEvents_) {
      readEventList(config().eventList());
    }
  }

  // Event list: "run subrun event" per line, read until the first line that does not parse.
  void MC2EventHits::readEventList(const std::string& fileName) {
    ConfigFileLookupPolicy resolve;
    std::ifstream in(resolve(fileName));
    int run;
    int subRun;
    int event;
    while (in >> run >> subRun >> event) {
      selectedEvents_.insert({run, subRun, event});
    }
  }

  void MC2EventHits::beginJob() {
    art::ServiceHandle<art::TFileService> tfs;
    tree_ = tfs->make<TTree>("hits", "calorimeter hits of selected events");
    tree_->Branch("run", &run_);
    tree_->Branch("subrun", &subRun_);
    tree_->Branch("event", &event_);
    // CaloHits
    tree_->Branch("cry", &hitCrystal_);
    tree_->Branch("disk", &hitDisk_);
    tree_->Branch("nSiPM", &hitNSiPMs_);
    tree_->Branch("x", &hitX_);
    tree_->Branch("y", &hitY_);
    tree_->Branch("t", &hitTime_);
    tree_->Branch("e", &hitEnergy_);
    // CaloRecoDigis
    tree_->Branch("rdSiPM", &digiSiPM_);
    tree_->Branch("rdCry", &digiCrystal_);
    tree_->Branch("rdT", &digiTime_);
    tree_->Branch("rdE", &digiEnergy_);
  }

  void MC2EventHits::analyze(const art::Event& event) {
    run_ = event.run();
    subRun_ = event.subRun();
    event_ = event.event();
    if (!isSelected()) {
      return;
    }
    clearHits();
    fillCaloHits(event);
    fillRecoDigis(event);
    tree_->Fill();
    ++nWritten_;
  }

  // Without eventList: the first maxEvents events; with eventList: the events in the list (run_, subRun_, event_ set).
  bool MC2EventHits::isSelected() const {
    if (allEvents_) {
      return nWritten_ < maxEvents_;
    }
    return selectedEvents_.count({run_, subRun_, event_}) > 0;
  }

  void MC2EventHits::clearHits() {
    hitCrystal_.clear();
    hitDisk_.clear();
    hitNSiPMs_.clear();
    hitX_.clear();
    hitY_.clear();
    hitTime_.clear();
    hitEnergy_.clear();
    digiSiPM_.clear();
    digiCrystal_.clear();
    digiTime_.clear();
    digiEnergy_.clear();
  }

  // CaloHits, with the position of the crystal (Crystal::position) in the tracker frame.
  void MC2EventHits::fillCaloHits(const art::Event& event) {
    const Calorimeter& calorimeter = *GeomHandle<Calorimeter>();
    for (const auto& hit : *event.getValidHandle<CaloHitCollection>(caloHitTag_)) {
      const auto& crystal = calorimeter.crystal(hit.crystalID());
      const auto position = calorimeter.mu2eToTracker(crystal.position());
      hitCrystal_.push_back(hit.crystalID());
      hitDisk_.push_back(crystal.diskID());
      hitNSiPMs_.push_back(hit.nSiPMs());
      hitX_.push_back(position.x());
      hitY_.push_back(position.y());
      hitTime_.push_back(hit.time());
      hitEnergy_.push_back(hit.energyDep());
    }
  }

  // CaloRecoDigis of single SiPMs; crystal id = SiPM id / 2.
  void MC2EventHits::fillRecoDigis(const art::Event& event) {
    for (const auto& recoDigi : *event.getValidHandle<CaloRecoDigiCollection>(recoDigiTag_)) {
      const int sipmId = recoDigi.SiPMID();
      digiSiPM_.push_back(sipmId);
      digiCrystal_.push_back(sipmId / 2);
      digiTime_.push_back(recoDigi.time());
      digiEnergy_.push_back(recoDigi.energyDep());
    }
  }

}  // namespace mu2e

DEFINE_ART_MODULE(mu2e::MC2EventHits)
