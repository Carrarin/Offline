//
// MC2CaloRecoDigiDump: writes the CaloRecoDigi of one collection to a flat tree "sidet" with the branches read by the MC2
// analysis (mc2-ana scripts/sidet_track_crystals.C), so that the cosmic-track timing check runs on the art reconstruction:
//   dtcEW = event number (= DTC event window tag for the MC2 raw data), nHits, and per reco digi: iBoa, iCha (DIRAC board
//   index and channel, from the offline SiPMID with the inverse of the raw -> offline map file), EWch (= dtcEW),
//   Tval = 0, templTime = CaloRecoDigi time [ns], Vmax = eDep / ADC2MeV [ADC], templChi2 = chi2/ndf, templFit[0..2] =
//   (amplitude ADC, time, 0), nSamples = waveform length, SiPMID, pileUp.
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
#include <fstream>
#include <unordered_map>

namespace mu2e {

  class MC2CaloRecoDigiDump : public art::EDAnalyzer {
    public:
      struct Config {
        using Name = fhicl::Name; using Comment = fhicl::Comment;
        fhicl::Atom<art::InputTag> recoDigiTag { Name("caloRecoDigiCollection"), Comment("CaloRecoDigi collection") };
        fhicl::Atom<std::string>   mapFile     { Name("rawToOfflineMap"), Comment("raw SiPM id -> offline SiPM id"), "Offline/CaloConditions/data/nominal.txt" };
        fhicl::Atom<std::string>   treeName    { Name("treeName"), Comment("output tree name"), "sidet" };
      };
      explicit MC2CaloRecoDigiDump(const art::EDAnalyzer::Table<Config>& config);
      void analyze(const art::Event& event) override;

    private:
      static constexpr int M = 5000;
      art::InputTag tag_;
      std::unordered_map<int, int> rawOf_;   // offline SiPMID -> raw id = board * 20 + channel
      ProditionsHandle<CalCalib> calCalibHandle_;
      TTree* t_ = nullptr;
      Long64_t ew_ = 0; int n_ = 0;
      int boa_[M], cha_[M], ns_[M], id_[M], pu_[M]; Long64_t ewch_[M]; double tv_[M], tt_[M], vm_[M], chi_[M], fit_[M][3];
  };

  MC2CaloRecoDigiDump::MC2CaloRecoDigiDump(const art::EDAnalyzer::Table<Config>& config) :
    art::EDAnalyzer{config}, tag_(config().recoDigiTag()) {
    ConfigFileLookupPolicy resolve; const std::string fn = resolve(config().mapFile());
    std::ifstream in(fn); if (!in) throw cet::exception("MC2DUMP") << "cannot open " << fn << "\n";
    int raw, off; while (in >> raw >> off) rawOf_[off] = raw;
    art::ServiceHandle<art::TFileService> tfs;
    t_ = tfs->make<TTree>(config().treeName().c_str(), "CaloRecoDigi for the MC2 timing checks");
    t_->Branch("dtcEW", &ew_, "dtcEW/L"); t_->Branch("nHits", &n_, "nHits/I");
    t_->Branch("iBoa", boa_, "iBoa[nHits]/I"); t_->Branch("iCha", cha_, "iCha[nHits]/I"); t_->Branch("EWch", ewch_, "EWch[nHits]/L");
    t_->Branch("Tval", tv_, "Tval[nHits]/D"); t_->Branch("templTime", tt_, "templTime[nHits]/D"); t_->Branch("Vmax", vm_, "Vmax[nHits]/D");
    t_->Branch("templChi2", chi_, "templChi2[nHits]/D"); t_->Branch("templFit", fit_, "templFit[nHits][3]/D");
    t_->Branch("nSamples", ns_, "nSamples[nHits]/I"); t_->Branch("SiPMID", id_, "SiPMID[nHits]/I"); t_->Branch("pileUp", pu_, "pileUp[nHits]/I");
  }

  void MC2CaloRecoDigiDump::analyze(const art::Event& event) {
    const auto& digis = *event.getValidHandle<CaloRecoDigiCollection>(tag_);
    const CalCalib& calib = calCalibHandle_.get(event.id());
    ew_ = event.event(); n_ = 0;
    for (const auto& d : digis) {
      if (n_ >= M) break;
      const int id = d.SiPMID(); const auto it = rawOf_.find(id); if (it == rawOf_.end()) continue;
      const double a2m = calib.ADC2MeV(id); const double amp = a2m > 0 ? d.energyDep() / a2m : 0.;
      const int k = n_++;
      boa_[k] = it->second / 20; cha_[k] = it->second % 20; ewch_[k] = ew_; id_[k] = id; pu_[k] = d.pileUp();
      tv_[k] = 0; tt_[k] = d.time(); vm_[k] = amp; chi_[k] = d.ndf() > 0 ? std::max(1e-6, double(d.chi2()) / d.ndf()) : 1e-6;
      fit_[k][0] = amp; fit_[k][1] = d.time(); fit_[k][2] = 0; ns_[k] = d.caloDigiPtr()->waveform().size();
    }
    t_->Fill();
  }

}

DEFINE_ART_MODULE(mu2e::MC2CaloRecoDigiDump)
