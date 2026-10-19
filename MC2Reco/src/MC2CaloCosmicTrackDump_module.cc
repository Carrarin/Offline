// MC2CaloCosmicTrackDump: flat tree "tracks" of the MC2CaloCosmicTracks, with the branches and units of the mc2tools tagger
// tracks.root (macros/cosmic_track_tagger.C), so that the mc2tools stages after `tracks` (track_time_check.C, ...) read it as it is:
//   run, ew (= art event number = DTC event window tag, the roottuple currentDTCEventWindow), type, n0, n1, nTimeClusters,
//   time, timeRms [samples of 5 ns], cx, cy, cz [mm, tagger frame: z = 0 at the crystal front face of disk 0 of the Offline
//   geometry, so disk 1 at 700 mm], dx, dy, dz (direction, up), lineRms, length [mm], ampSum [ADC counts],
//   hitKey (raw id of SiPM 0 of the crystal = 20 * board + ROC channel, the tagger key), hitTime [samples], hitAmp [ADC counts,
//   mean of the two SiPMs].
// Plus, in the Offline tracker frame and ns: zTrk (cz in the tracker frame), timeNs, dtDisk, eSum [MeV], hitCrystal, hitTimeNs, hitE.
// Histograms: crystal time - track time, crystals per track, angle from the vertical (one-disk), time between disks.
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art_root_io/TFileService.h"
#include "fhiclcpp/types/Atom.h"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/CalorimeterGeom/inc/Calorimeter.hh"
#include "Offline/CaloConditions/inc/CalCalib.hh"
#include "Offline/CaloConditions/inc/CalDAQMap.hh"
#include "Offline/ProditionsService/inc/ProditionsHandle.hh"
#include "Offline/RecoDataProducts/inc/CaloRecoDigi.hh"
#include "Offline/MC2Reco/inc/MC2CaloCosmicTrack.hh"
#include "TTree.h"
#include "TH1F.h"
#include <cmath>
#include <vector>

namespace mu2e {
  class MC2CaloCosmicTrackDump : public art::EDAnalyzer {
    public:
      struct Config {
        using Name = fhicl::Name; using Comment = fhicl::Comment;
        fhicl::Atom<art::InputTag> trkTag { Name("trackCollection"), Comment("MC2CaloCosmicTracks"), "MC2CaloCosmicTrackFinder" };
        fhicl::Atom<double>        nsPerSample { Name("nsPerSample"), Comment("for the tagger units"), 5. };
      };
      using Parameters = art::EDAnalyzer::Table<Config>;
      explicit MC2CaloCosmicTrackDump(const Parameters& conf);
      void beginJob() override;
      void analyze(const art::Event& event) override;
    private:
      art::InputTag tag_; double nsps_; ProditionsHandle<CalCalib> calibH_; ProditionsHandle<CalDAQMap> mapH_;
      TTree* t_ = nullptr; TH1F *hTres_ = nullptr, *hNc_[2] = {nullptr, nullptr}, *hAng_ = nullptr, *hDtDisk_ = nullptr;
      int run_, type_, n0_, n1_, ncl_; Long64_t ew_; double time_, trms_, cx_, cy_, cz_, dx_, dy_, dz_, rms_, len_, amp_, zTrk_, tns_, dtd_, esum_;
      std::vector<int> key_, cry_; std::vector<double> ht_, ha_, htns_, he_;
  };

  MC2CaloCosmicTrackDump::MC2CaloCosmicTrackDump(const Parameters& conf) : art::EDAnalyzer(conf), tag_(conf().trkTag()), nsps_(conf().nsPerSample()) {
    consumes<MC2CaloCosmicTrackCollection>(tag_);
  }

  void MC2CaloCosmicTrackDump::beginJob() {
    art::ServiceHandle<art::TFileService> tfs;
    t_ = tfs->make<TTree>("tracks", "calorimeter cosmic tracks (MC2CaloCosmicTrackFinder), branches of the mc2tools tagger");
    t_->Branch("run", &run_); t_->Branch("ew", &ew_); t_->Branch("type", &type_); t_->Branch("n0", &n0_); t_->Branch("n1", &n1_); t_->Branch("nTimeClusters", &ncl_);
    t_->Branch("time", &time_); t_->Branch("timeRms", &trms_); t_->Branch("cx", &cx_); t_->Branch("cy", &cy_); t_->Branch("cz", &cz_);
    t_->Branch("dx", &dx_); t_->Branch("dy", &dy_); t_->Branch("dz", &dz_); t_->Branch("lineRms", &rms_); t_->Branch("length", &len_); t_->Branch("ampSum", &amp_);
    t_->Branch("hitKey", &key_); t_->Branch("hitTime", &ht_); t_->Branch("hitAmp", &ha_);
    t_->Branch("zTrk", &zTrk_); t_->Branch("timeNs", &tns_); t_->Branch("dtDisk", &dtd_); t_->Branch("eSum", &esum_);
    t_->Branch("hitCrystal", &cry_); t_->Branch("hitTimeNs", &htns_); t_->Branch("hitE", &he_);
    hTres_ = tfs->make<TH1F>("hTres", ";crystal time - track time [ns];crystals", 200, -10, 10);
    hNc_[0] = tfs->make<TH1F>("hNc1", "one-disk tracks;crystals on the track;tracks", 40, 0, 40);
    hNc_[1] = tfs->make<TH1F>("hNc3", "two-disk tracks;crystals on the track;tracks", 40, 0, 40);
    hAng_ = tfs->make<TH1F>("hAng", ";angle from the vertical in the x-y plane [deg];one-disk tracks", 90, 0, 90);
    hDtDisk_ = tfs->make<TH1F>("hDtDisk", ";mean time of disk 1 - disk 0 [ns];two-disk tracks", 120, -15, 15);
  }

  void MC2CaloCosmicTrackDump::analyze(const art::Event& event) {
    const Calorimeter& cal = *GeomHandle<Calorimeter>();
    const CalCalib& calib = calibH_.get(event.id()); const CalDAQMap& dmap = mapH_.get(event.id());
    double zFace0 = NAN;   // crystal front face of disk 0, tracker frame, Offline geometry (origin of the tagger z)
    for (unsigned i = 0; i < cal.nCrystals(); ++i) if (cal.crystal(i).diskID() == 0) { zFace0 = cal.mu2eToTracker(cal.crystal(i).position()).z(); break; }
    const auto& trks = *event.getValidHandle<MC2CaloCosmicTrackCollection>(tag_);
    for (const auto& tk : trks) {
      run_ = event.run(); ew_ = event.event(); type_ = tk.type; n0_ = tk.n0; n1_ = tk.n1; ncl_ = tk.nTimeClusters;
      time_ = tk.time / nsps_; trms_ = tk.timeRms / nsps_; tns_ = tk.time; dtd_ = tk.dtDisk; esum_ = tk.eSum;
      cx_ = tk.cx; cy_ = tk.cy; zTrk_ = tk.cz; cz_ = tk.cz - zFace0; dx_ = tk.ux; dy_ = tk.uy; dz_ = tk.uz; rms_ = tk.lineRms; len_ = tk.length;
      key_.clear(); cry_.clear(); ht_.clear(); ha_.clear(); htns_.clear(); he_.clear(); amp_ = 0;
      for (const auto& hp : tk.hits) {
        const int cid = hp->crystalID(); double a = 0; int na = 0;
        for (const auto& rd : hp->recoCaloDigis()) { const double a2m = calib.ADC2MeV(rd->SiPMID()); if (a2m > 0) { a += rd->energyDep() / a2m; ++na; } }
        if (na) a /= na;
        const CaloRawSiPMId raw = dmap.rawId(CaloSiPMId(cid * CaloConst::_nSiPMPerCrystal));   // SiPM 0 of the crystal
        key_.push_back(raw.dirac() * CaloConst::_nChPerDIRAC + raw.ROCchannel()); cry_.push_back(cid);
        ht_.push_back(hp->time() / nsps_); htns_.push_back(hp->time()); ha_.push_back(a); he_.push_back(hp->energyDep()); amp_ += a;
        hTres_->Fill(hp->time() - tk.time);
      }
      hNc_[tk.type == 3 ? 1 : 0]->Fill(tk.hits.size());
      if (tk.type < 3) hAng_->Fill(std::atan2(std::fabs(tk.ux), std::fabs(tk.uy)) * 180. / M_PI); else hDtDisk_->Fill(tk.dtDisk);
      t_->Fill();
    }
  }
}
DEFINE_ART_MODULE(mu2e::MC2CaloCosmicTrackDump)
