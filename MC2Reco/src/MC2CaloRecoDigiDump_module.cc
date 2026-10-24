//
// MC2CaloRecoDigiDump: writes the CaloRecoDigi of one collection to a flat tree with the branches read by the MC2
// analysis (mc2-ana scripts/sidet_track_crystals.C), so that the cosmic-track timing check runs on the art
// reconstruction.
//
// Input: CaloRecoDigiCollection (caloRecoDigiCollection); CalCalib conditions (ADC2MeV); the raw -> offline SiPM map
// file (rawToOfflineMap, "raw offline" per line), used in reverse to go from the offline SiPM id back to board and
// channel.
//
// Output: TTree treeName (default "sidet") in the TFileService file, one entry per event:
//   dtcEW     = event number (= DTC event-window tag for the MC2 raw data)
//   nHits     = number of reco digis written (at most kMaxDigis; digis of SiPMs not in the map are skipped)
//   per reco digi:
//     iBoa, iCha  DIRAC board index and channel (raw id = board * 20 + channel)
//     EWch        = dtcEW
//     Tval        = 0
//     templTime   = CaloRecoDigi time [ns]
//     Vmax        = eDep / ADC2MeV [ADC] (0 if ADC2MeV <= 0)
//     templChi2   = chi2/ndf (at least 1e-6; 1e-6 if ndf = 0)
//     templFit    = (amplitude [ADC], time [ns], 0)
//     nSamples    = waveform length of the CaloDigi
//     SiPMID, pileUp
//
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/types/Atom.h"
#include "cetlib_except/exception.h"

#include "Offline/CaloConditions/inc/CalCalib.hh"
#include "Offline/ConfigTools/inc/ConfigFileLookupPolicy.hh"
#include "Offline/ProditionsService/inc/ProditionsHandle.hh"
#include "Offline/RecoDataProducts/inc/CaloRecoDigi.hh"

#include "TTree.h"

#include <algorithm>
#include <fstream>
#include <string>
#include <unordered_map>

namespace mu2e {

  class MC2CaloRecoDigiDump : public art::EDAnalyzer {
    public:
      struct Config {
        using Name = fhicl::Name;
        using Comment = fhicl::Comment;

        fhicl::Atom<art::InputTag> recoDigiTag{Name("caloRecoDigiCollection"), Comment("CaloRecoDigi collection")};
        fhicl::Atom<std::string> mapFile{Name("rawToOfflineMap"), Comment("raw SiPM id -> offline SiPM id"),
            "Offline/CaloConditions/data/nominal.txt"};
        fhicl::Atom<std::string> treeName{Name("treeName"), Comment("output tree name"), "sidet"};
      };

      explicit MC2CaloRecoDigiDump(const art::EDAnalyzer::Table<Config>& config);
      void analyze(const art::Event& event) override;

    private:
      static constexpr int kMaxDigis = 5000;   // size of the per-digi arrays of the tree

      void readRawToOfflineMap(const std::string& fileName);
      void bookTree(const std::string& treeName);

      // Inputs
      art::InputTag recoDigiTag_;
      std::unordered_map<int, int> rawIdOfOffline_;   // offline SiPMID -> raw id = board * 20 + channel
      ProditionsHandle<CalCalib> calCalibHandle_;

      // Output tree and its branches
      TTree* tree_ = nullptr;
      Long64_t eventWindow_ = 0;
      int nHits_ = 0;
      int board_[kMaxDigis];
      int channel_[kMaxDigis];
      int nSamples_[kMaxDigis];
      int sipmId_[kMaxDigis];
      int pileUp_[kMaxDigis];
      Long64_t eventWindowPerHit_[kMaxDigis];
      double tval_[kMaxDigis];
      double templTime_[kMaxDigis];
      double vmax_[kMaxDigis];
      double chi2ndf_[kMaxDigis];
      double templFit_[kMaxDigis][3];
  };

  MC2CaloRecoDigiDump::MC2CaloRecoDigiDump(const art::EDAnalyzer::Table<Config>& config) :
    art::EDAnalyzer{config},
    recoDigiTag_(config().recoDigiTag())
  {
    readRawToOfflineMap(config().mapFile());
    bookTree(config().treeName());
  }

  // Map file: "raw offline" per line, read until the first line that does not parse; stored as offline -> raw.
  void MC2CaloRecoDigiDump::readRawToOfflineMap(const std::string& fileName) {
    ConfigFileLookupPolicy resolve;
    const std::string path = resolve(fileName);
    std::ifstream in(path);
    if (!in) {
      throw cet::exception("MC2DUMP") << "cannot open " << path << "\n";
    }
    int rawId;
    int offlineId;
    while (in >> rawId >> offlineId) {
      rawIdOfOffline_[offlineId] = rawId;
    }
  }

  void MC2CaloRecoDigiDump::bookTree(const std::string& treeName) {
    art::ServiceHandle<art::TFileService> tfs;
    tree_ = tfs->make<TTree>(treeName.c_str(), "CaloRecoDigi for the MC2 timing checks");
    tree_->Branch("dtcEW", &eventWindow_, "dtcEW/L");
    tree_->Branch("nHits", &nHits_, "nHits/I");
    tree_->Branch("iBoa", board_, "iBoa[nHits]/I");
    tree_->Branch("iCha", channel_, "iCha[nHits]/I");
    tree_->Branch("EWch", eventWindowPerHit_, "EWch[nHits]/L");
    tree_->Branch("Tval", tval_, "Tval[nHits]/D");
    tree_->Branch("templTime", templTime_, "templTime[nHits]/D");
    tree_->Branch("Vmax", vmax_, "Vmax[nHits]/D");
    tree_->Branch("templChi2", chi2ndf_, "templChi2[nHits]/D");
    tree_->Branch("templFit", templFit_, "templFit[nHits][3]/D");
    tree_->Branch("nSamples", nSamples_, "nSamples[nHits]/I");
    tree_->Branch("SiPMID", sipmId_, "SiPMID[nHits]/I");
    tree_->Branch("pileUp", pileUp_, "pileUp[nHits]/I");
  }

  void MC2CaloRecoDigiDump::analyze(const art::Event& event) {
    const auto& recoDigis = *event.getValidHandle<CaloRecoDigiCollection>(recoDigiTag_);
    const CalCalib& calCalib = calCalibHandle_.get(event.id());
    eventWindow_ = event.event();
    nHits_ = 0;

    for (const auto& recoDigi : recoDigis) {
      if (nHits_ >= kMaxDigis) {
        break;
      }
      const int sipmId = recoDigi.SiPMID();
      const auto raw = rawIdOfOffline_.find(sipmId);
      if (raw == rawIdOfOffline_.end()) {
        continue;
      }
      // amplitude in ADC back from the calibrated energy
      const double adc2MeV = calCalib.ADC2MeV(sipmId);
      const double amplitude = adc2MeV > 0 ? recoDigi.energyDep() / adc2MeV : 0.;

      const int i = nHits_++;
      board_[i] = raw->second / 20;
      channel_[i] = raw->second % 20;
      eventWindowPerHit_[i] = eventWindow_;
      sipmId_[i] = sipmId;
      pileUp_[i] = recoDigi.pileUp();
      tval_[i] = 0;
      templTime_[i] = recoDigi.time();
      vmax_[i] = amplitude;
      chi2ndf_[i] = recoDigi.ndf() > 0 ? std::max(1e-6, double(recoDigi.chi2()) / recoDigi.ndf()) : 1e-6;
      templFit_[i][0] = amplitude;
      templFit_[i][1] = recoDigi.time();
      templFit_[i][2] = 0;
      nSamples_[i] = recoDigi.caloDigiPtr()->waveform().size();
    }
    tree_->Fill();
  }

}  // namespace mu2e

DEFINE_ART_MODULE(mu2e::MC2CaloRecoDigiDump)
