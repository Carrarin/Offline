//
// MC2CaloCosmicTrackDump: flat tree of the MC2CaloCosmicTracks in the format of the mc2tools tagger.
//
// The tree "tracks" has the branches and units of the tagger tracks.root (macros/cosmic_track_tagger.C), so that the
// mc2tools stages after `tracks` (track_time_check.C, ...) read it as it is.
//
// Input:  MC2CaloCosmicTrackCollection (trackCollection); CalCalib (ADC2MeV, for the amplitudes in ADC counts) and
//         CalDAQMap (raw ids) conditions; Calorimeter geometry.
// Output: TTree "tracks", one entry per track:
//   tagger branches
//     run, ew (= art event number = DTC event window tag, the roottuple currentDTCEventWindow), type, n0, n1,
//     nTimeClusters, time, timeRms [samples of nsPerSample = 5 ns],
//     cx, cy, cz [mm, tagger frame: z = 0 at the crystal front face of disk 0 of the Offline geometry, so disk 1 at
//     700 mm], dx, dy, dz (direction, pointing up), lineRms, length [mm], ampSum [ADC counts],
//     hitKey (raw id of SiPM 0 of the crystal = 20 * board + ROC channel, the tagger key), hitTime [samples],
//     hitAmp [ADC counts, mean of the two SiPMs]
//   plus, in the Offline tracker frame and ns
//     zTrk (cz in the tracker frame), timeNs, dtDisk, eSum [MeV], hitCrystal, hitTimeNs, hitE [MeV]
// Histograms: crystal time - track time, crystals per track (one-disk, two-disk), angle from the vertical (one-disk),
//   time between disks (two-disk).
//
// Steps, for every track: track branches, crystal branches (and crystal time residuals), histograms, tree entry.
//
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
        using Name = fhicl::Name;
        using Comment = fhicl::Comment;

        fhicl::Atom<art::InputTag> trkTag{Name("trackCollection"), Comment("MC2CaloCosmicTracks"),
            "MC2CaloCosmicTrackFinder"};
        fhicl::Atom<double> nsPerSample{Name("nsPerSample"), Comment("for the tagger units"), 5.};
      };
      using Parameters = art::EDAnalyzer::Table<Config>;

      explicit MC2CaloCosmicTrackDump(const Parameters& conf);
      void beginJob() override;
      void analyze(const art::Event& event) override;

    private:
      static double frontFaceZDisk0(const Calorimeter& calorimeter);
      static double meanSiPMAmplitude(const CaloHit& hit, const CalCalib& calCalib);
      void fillTrackBranches(const art::Event& event, const MC2CaloCosmicTrack& track, double zFace0);
      void fillCrystalBranches(const MC2CaloCosmicTrack& track, const CalCalib& calCalib, const CalDAQMap& daqMap);
      void fillTrackHistograms(const MC2CaloCosmicTrack& track);

      // Parameters
      art::InputTag trackTag_;
      double nsPerSample_;

      // Conditions
      ProditionsHandle<CalCalib> calCalibHandle_;
      ProditionsHandle<CalDAQMap> daqMapHandle_;

      // Histograms
      TH1F* hTimeResidual_ = nullptr;                        // crystal time - track time
      TH1F* hCrystalsPerTrack_[2] = {nullptr, nullptr};      // [0] one-disk, [1] two-disk tracks
      TH1F* hAngle_ = nullptr;                               // angle from the vertical, one-disk tracks
      TH1F* hDtDisk_ = nullptr;                              // time disk 1 - disk 0, two-disk tracks

      // Tree "tracks" and its branches (branch name in the comment when different)
      TTree* tree_ = nullptr;
      int run_;
      Long64_t eventWindow_;      // ew
      int type_;
      int n0_;
      int n1_;
      int nTimeClusters_;
      double time_;               // [samples]
      double timeRms_;            // timeRms [samples]
      double cx_;
      double cy_;
      double cz_;                 // tagger frame
      double dx_;
      double dy_;
      double dz_;
      double lineRms_;
      double length_;
      double ampSum_;
      double zTracker_;           // zTrk
      double timeNs_;
      double dtDisk_;
      double eSum_;
      std::vector<int> hitKey_;
      std::vector<int> hitCrystal_;
      std::vector<double> hitTime_;     // [samples]
      std::vector<double> hitAmp_;
      std::vector<double> hitTimeNs_;
      std::vector<double> hitE_;
  };

  MC2CaloCosmicTrackDump::MC2CaloCosmicTrackDump(const Parameters& conf) :
    art::EDAnalyzer(conf),
    trackTag_(conf().trkTag()),
    nsPerSample_(conf().nsPerSample())
  {
    consumes<MC2CaloCosmicTrackCollection>(trackTag_);
  }

  void MC2CaloCosmicTrackDump::beginJob() {
    art::ServiceHandle<art::TFileService> tfs;
    tree_ = tfs->make<TTree>("tracks",
        "calorimeter cosmic tracks (MC2CaloCosmicTrackFinder), branches of the mc2tools tagger");

    // tagger branches
    tree_->Branch("run", &run_);
    tree_->Branch("ew", &eventWindow_);
    tree_->Branch("type", &type_);
    tree_->Branch("n0", &n0_);
    tree_->Branch("n1", &n1_);
    tree_->Branch("nTimeClusters", &nTimeClusters_);
    tree_->Branch("time", &time_);
    tree_->Branch("timeRms", &timeRms_);
    tree_->Branch("cx", &cx_);
    tree_->Branch("cy", &cy_);
    tree_->Branch("cz", &cz_);
    tree_->Branch("dx", &dx_);
    tree_->Branch("dy", &dy_);
    tree_->Branch("dz", &dz_);
    tree_->Branch("lineRms", &lineRms_);
    tree_->Branch("length", &length_);
    tree_->Branch("ampSum", &ampSum_);
    tree_->Branch("hitKey", &hitKey_);
    tree_->Branch("hitTime", &hitTime_);
    tree_->Branch("hitAmp", &hitAmp_);

    // Offline frame and units
    tree_->Branch("zTrk", &zTracker_);
    tree_->Branch("timeNs", &timeNs_);
    tree_->Branch("dtDisk", &dtDisk_);
    tree_->Branch("eSum", &eSum_);
    tree_->Branch("hitCrystal", &hitCrystal_);
    tree_->Branch("hitTimeNs", &hitTimeNs_);
    tree_->Branch("hitE", &hitE_);

    hTimeResidual_ = tfs->make<TH1F>("hTres", ";crystal time - track time [ns];crystals", 200, -10, 10);
    hCrystalsPerTrack_[0] = tfs->make<TH1F>("hNc1", "one-disk tracks;crystals on the track;tracks", 40, 0, 40);
    hCrystalsPerTrack_[1] = tfs->make<TH1F>("hNc3", "two-disk tracks;crystals on the track;tracks", 40, 0, 40);
    hAngle_ = tfs->make<TH1F>("hAng", ";angle from the vertical in the x-y plane [deg];one-disk tracks", 90, 0, 90);
    hDtDisk_ = tfs->make<TH1F>("hDtDisk", ";mean time of disk 1 - disk 0 [ns];two-disk tracks", 120, -15, 15);
  }

  void MC2CaloCosmicTrackDump::analyze(const art::Event& event) {
    const Calorimeter& calorimeter = *GeomHandle<Calorimeter>();
    const CalCalib& calCalib = calCalibHandle_.get(event.id());
    const CalDAQMap& daqMap = daqMapHandle_.get(event.id());
    const double zFace0 = frontFaceZDisk0(calorimeter);

    const auto& tracks = *event.getValidHandle<MC2CaloCosmicTrackCollection>(trackTag_);
    for (const auto& track : tracks) {
      fillTrackBranches(event, track, zFace0);
      fillCrystalBranches(track, calCalib, daqMap);
      fillTrackHistograms(track);
      tree_->Fill();
    }
  }

  // Crystal front face of disk 0 in the tracker frame, Offline geometry: the origin of the tagger z.
  // (Crystal::position() is the centre of the front face.) NAN if disk 0 has no crystals.
  double MC2CaloCosmicTrackDump::frontFaceZDisk0(const Calorimeter& calorimeter) {
    double zFace0 = NAN;
    for (unsigned i = 0; i < calorimeter.nCrystals(); ++i) {
      if (calorimeter.crystal(i).diskID() == 0) {
        zFace0 = calorimeter.mu2eToTracker(calorimeter.crystal(i).position()).z();
        break;
      }
    }
    return zFace0;
  }

  // Mean amplitude of the SiPMs of a crystal [ADC counts]: energy / ADC2MeV, averaged over the SiPMs with a valid
  // calibration (0 if none). This is the fit amplitude, also when the finder cut on the raw amplitude.
  double MC2CaloCosmicTrackDump::meanSiPMAmplitude(const CaloHit& hit, const CalCalib& calCalib) {
    double amplitude = 0;
    int nAmplitudes = 0;
    for (const auto& recoDigi : hit.recoCaloDigis()) {
      const double adcToMeV = calCalib.ADC2MeV(recoDigi->SiPMID());
      if (adcToMeV > 0) {
        amplitude += recoDigi->energyDep() / adcToMeV;
        ++nAmplitudes;
      }
    }
    if (nAmplitudes) {
      amplitude /= nAmplitudes;
    }
    return amplitude;
  }

  // Track-level branches; times converted to samples and cz to the tagger frame.
  void MC2CaloCosmicTrackDump::fillTrackBranches(const art::Event& event, const MC2CaloCosmicTrack& track,
                                                 double zFace0) {
    run_ = event.run();
    eventWindow_ = event.event();
    type_ = track.type;
    n0_ = track.n0;
    n1_ = track.n1;
    nTimeClusters_ = track.nTimeClusters;
    time_ = track.time / nsPerSample_;
    timeRms_ = track.timeRms / nsPerSample_;
    timeNs_ = track.time;
    dtDisk_ = track.dtDisk;
    eSum_ = track.eSum;
    cx_ = track.cx;
    cy_ = track.cy;
    zTracker_ = track.cz;
    cz_ = track.cz - zFace0;
    dx_ = track.ux;
    dy_ = track.uy;
    dz_ = track.uz;
    lineRms_ = track.lineRms;
    length_ = track.length;
  }

  // Crystal branches (vectors, one element per crystal of the track), their amplitude sum ampSum, and the
  // crystal time residuals.
  void MC2CaloCosmicTrackDump::fillCrystalBranches(const MC2CaloCosmicTrack& track, const CalCalib& calCalib,
                                                   const CalDAQMap& daqMap) {
    hitKey_.clear();
    hitCrystal_.clear();
    hitTime_.clear();
    hitAmp_.clear();
    hitTimeNs_.clear();
    hitE_.clear();
    ampSum_ = 0;
    for (const auto& hit : track.hits) {
      const int crystalId = hit->crystalID();
      const double amplitude = meanSiPMAmplitude(*hit, calCalib);
      // tagger key: raw id of SiPM 0 of the crystal
      const CaloRawSiPMId rawId = daqMap.rawId(CaloSiPMId(crystalId * CaloConst::_nSiPMPerCrystal));
      hitKey_.push_back(rawId.dirac() * CaloConst::_nChPerDIRAC + rawId.ROCchannel());
      hitCrystal_.push_back(crystalId);
      hitTime_.push_back(hit->time() / nsPerSample_);
      hitTimeNs_.push_back(hit->time());
      hitAmp_.push_back(amplitude);
      hitE_.push_back(hit->energyDep());
      ampSum_ += amplitude;
      hTimeResidual_->Fill(hit->time() - track.time);
    }
  }

  void MC2CaloCosmicTrackDump::fillTrackHistograms(const MC2CaloCosmicTrack& track) {
    hCrystalsPerTrack_[track.type == 3 ? 1 : 0]->Fill(track.hits.size());
    if (track.type < 3) {
      hAngle_->Fill(std::atan2(std::fabs(track.ux), std::fabs(track.uy)) * 180. / M_PI);
    } else {
      hDtDisk_->Fill(track.dtDisk);
    }
  }

}

DEFINE_ART_MODULE(mu2e::MC2CaloCosmicTrackDump)
