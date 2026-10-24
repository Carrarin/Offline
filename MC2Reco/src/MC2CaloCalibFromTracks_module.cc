//
// MC2CaloCalibFromTracks: our calibration machine. Per-SiPM time offsets and energy scale of the calorimeter from the
// cosmic tracks produced in the same job (MC2CaloCosmicTrack, MC2CaloCosmicTrackFinder), written as Offline condition
// tables (CalTimeCalib, CalEnergyCalib) in the format of the calorimeter group's tables, so that the two can be
// compared line by line. It also measures the detector response for the simulation (baselines, noise, thresholds,
// effective photoelectrons).
//
// Starting point: the CalCalib conditions the reconstruction used (e.g. the tables of the calorimeter group loaded with
// DbService.textFile). The hit times are already corrected with those offsets (MC2CaloRecoDigiMaker
// useCalCalibTimeOffset) and the SiPM amplitude is A = CaloRecoDigi energy / CalCalib ADC2MeV.
//
// Inputs (products):
//   caloTrackCollection     MC2CaloCosmicTrackCollection, the calorimeter cosmic tracks
//   caloRecoDigiCollection  CaloRecoDigiCollection, all reco digis of the event (baselines, noise, thresholds)
//   kalSeedCollection       optional KalSeedCollection, tracker lines for the calorimeter-tracker alignment check
//
// Per event (analyze):
//   1. at the first event: board (DIRAC) of every SiPM from CalDAQMap, starting time offsets and ADC2MeV from CalCalib;
//   2. all pulses: baseline and noise (mean and rms of the first nPed samples of every CaloDigi waveform), amplitude
//      spectra per SiPM (A of our fit and raw = largest sample - baseline) for the thresholds;
//   3. every crystal of every track: position s of the crystal centre along the line, path in the crystal, one time
//      entry per SiPM (time calibration), amplitude per 34 mm (energy calibration), left/right amplitudes when both
//      SiPMs are there (photostatistics), and the time difference to the tracker lines (alignment check).
//
// At the end of the job (endJob):
//   Time, step 1 (board step, robust to jumps, e.g. after a re-initialisation of a board): track time = median of its
//     SiPMs (insensitive to a minority of jumped boards); per board the mode of the residuals over +-3 us (10 ns
//     sliding window), then the mean within +-10 ns of it; 3 iterations; the mean board correction is fixed to 0.
//   Time, step 2 (global least squares, as mc2tools time_calibration lsq): for every track k and every SiPM j of its
//     crystals,
//       t_kj = T_k - s_kj / c + o_j,  s_kj = position of the crystal centre along the line (up), c = speed of light;
//     T_k and o_j by alternating least squares (nIterTime), starting from the board offsets, with the mean of the o_j
//     fixed to 0 (the absolute alignment of the starting tables is kept). From iteration nIterBeforeCut, every 5
//     iterations, entries with |residual| > 200, 67, 22, 7.4 ns, then outlierCut, are left out (progressive: offsets of
//     hundreds of ns, as in a pass without tables, can converge first). One-disk tracks: s along the line in the disk
//     plane (the 3D direction is unknown).
//     New offset = starting CalCalib offset + o_j. SiPMs with fewer than minEntriesTime entries: starting offset + the
//     correction of their board (flagged in the stats file).
//   Energy (as mc2tools energy_calibration.py): chord of the line in the 34 x 34 mm cross-section of the crystal (x-y);
//     two-disk tracks: chord / sqrt(1 - uz^2); chords < minChordFrac * 34 mm are left out;
//     A34 = A * dp(34 mm) / dp(path), dp = Landau most probable loss in CsI (PDG). Per SiPM, Landau (x) Gauss (ROOT
//     tutorial langaus) on the one-disk tracks only, histogram of 120 bins over [0, 3 median], fit range
//     [0.5, 2.5] x peak, >= minEntriesEnergy entries, accepted if the error is < 10 % of the MPV.
//     ADC2MeV = dp(34 mm) / MPV (SiPMs without a good fit: median ADC2MeV of the fitted SiPMs, flagged).
//   Thresholds: per SiPM, the rising edge of the amplitude spectrum (raw and fit) fitted with
//       0.5 (1 + erf((x - T) / (sqrt2 w))) N exp(-(x - T) / lambda);
//     T = threshold [ADC], w = width of the edge (noise on the threshold); in MeV with the new ADC2MeV of the SiPM.
//     Fallback: half-maximum crossing of the rising edge (w = -1).
//   Effective photoelectrons: per crystal, from the crystals of the tracks with both SiPMs, amplitudes brought to the
//     mean gain of the two, p0 = (Var(L - R) - noise variance) / <L + R> [ADC per effective photoelectron] (robust
//     variance); pePerMeV of each SiPM = ADCPerMeV / p0. Crystals without enough entries: median p0 (flagged).
//     Var(L - R) also contains the left/right sharing with the transverse position: p0 is an upper estimate.
//   Optional (kalSeedCollection not empty): tracker lines vs the calo hits of the tracks, t_SiPM(corrected) - t_line,
//     histograms only (calorimeter-tracker alignment check).
//
// Outputs (text files, read by DbService and by scripts: formats are fixed):
//   timeTableFile      TABLE CalTimeCalib     roid,timeOffset
//   energyTableFile    TABLE CalEnergyCalib   roid,ADC2MeV
//   statsFile          csv, per SiPM: entries, corrections, MIP fit, flags
//   simTableFile       TABLE CalSimCrystals   crid, LRU (defaultLRU, not measured here), pePerMeV0, pePerMeV1,
//                                             ADCPerMeV0, ADCPerMeV1
//   baselineTableFile  TABLE CalBaselines     roid, baseline, threshold [raw ADC]
//   responseFile       csv, per SiPM: pulses, baseline, noise, thresholds, ADC2MeV, p0, pePerMeV
// TFileService histograms: hTresBefore, hTres (time residuals before/after the least squares), hOff (corrections to
//   the starting offsets), hMPV, hChord, hAmpFit, hAmpRaw (amplitude spectra per SiPM), hThr, hThrMeV, hP0, and with
//   kalSeedCollection hTrkDt, hTrkDtSiPM.
//
// Main parameters: depthInCrystal and diskZOffset (crystal centres), nIterTime, nIterBeforeCut, outlierCut,
// minEntriesTime (time), crystalSide, minChordFrac, minEntriesEnergy (energy), nPed, thrMaxADC, minPulsesThr,
// minEntriesPE (response).
//
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art_root_io/TFileService.h"
#include "cetlib_except/exception.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/CalorimeterGeom/inc/Calorimeter.hh"
#include "Offline/CaloConditions/inc/CalCalib.hh"
#include "Offline/CaloConditions/inc/CalDAQMap.hh"
#include "Offline/ProditionsService/inc/ProditionsHandle.hh"
#include "Offline/DataProducts/inc/CaloConst.hh"
#include "Offline/RecoDataProducts/inc/KalSeed.hh"
#include "Offline/RecoDataProducts/inc/CaloHit.hh"
#include "Offline/RecoDataProducts/inc/CaloRecoDigi.hh"
#include "Offline/RecoDataProducts/inc/CaloDigi.hh"
#include "KinKal/Trajectory/KinematicLine.hh"
#include "Offline/MC2Reco/inc/MC2CaloCosmicTrack.hh"
#include "TF1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <vector>

namespace mu2e {
  namespace {
    const double kC = 299.792458;   // speed of light [mm/ns]

    // Landau most probable energy loss in CsI for a path x (mm), beta = 1
    // (PDG; same constants as mc2tools energy_calibration.py).
    double mostProbableLoss(double x) {
      const double K = 0.307075;
      const double ZA = 0.41569;
      const double RHO = 4.510;
      const double HWP = 39.46e-6;
      const double ME = 0.51099895;
      const double J = 0.200;
      const double xi = 0.5 * K * ZA * RHO * x / 10.0;
      return xi * (std::log(2 * ME * xi / (HWP * HWP)) + J);
    }

    // Landau (x) Gauss, from the ROOT tutorial langaus.
    // par = width (Landau), MPV, area, sigma (Gauss).
    double langaus(double* x, double* par) {
      const double invsq2pi = 0.3989422804014;
      const double mpshift = -0.22278298;   // Landau maximum is at mpshift in the ROOT parametrisation
      const double np = 100.0;              // number of convolution steps
      const double sc = 5.0;                // convolution extends to +-sc Gaussian sigmas
      const double mpc = par[1] - mpshift * par[0];
      const double xlow = x[0] - sc * par[3];
      const double xupp = x[0] + sc * par[3];
      const double step = (xupp - xlow) / np;
      double sum = 0;
      for (double i = 1.0; i <= np / 2; i++) {
        double xx = xlow + (i - .5) * step;
        double landau = TMath::Landau(xx, mpc, par[0]) / par[0];
        sum += landau * TMath::Gaus(x[0], xx, par[3]);
        xx = xupp - (i - .5) * step;
        landau = TMath::Landau(xx, mpc, par[0]) / par[0];
        sum += landau * TMath::Gaus(x[0], xx, par[3]);
      }
      return par[2] * step * sum * invsq2pi / par[3];
    }

    // Chord of the 2D line (point p, unit direction d) in the square of side a centred at (cx, cy); 0 if it misses.
    // Slab method: intersection of the parameter intervals inside the square along x and along y.
    double chord(double px, double py, double dx, double dy, double cx, double cy, double a) {
      double tIn = -1e30;
      double tOut = 1e30;
      const double lo[2] = {cx - a / 2, cy - a / 2};
      const double hi[2] = {cx + a / 2, cy + a / 2};
      const double p[2] = {px, py};
      const double d[2] = {dx, dy};
      for (int axis = 0; axis < 2; ++axis) {
        if (std::abs(d[axis]) < 1e-12) {   // parallel to this side
          if (p[axis] < lo[axis] || p[axis] > hi[axis]) {
            return 0;
          }
          continue;
        }
        double t1 = (lo[axis] - p[axis]) / d[axis];
        double t2 = (hi[axis] - p[axis]) / d[axis];
        if (t1 > t2) {
          std::swap(t1, t2);
        }
        tIn = std::max(tIn, t1);
        tOut = std::min(tOut, t2);
      }
      return tOut > tIn ? tOut - tIn : 0;
    }

    // One SiPM time on one track, for the time calibration.
    struct TimeEntry {
      int track;       // index of the track in the job
      int sipm;        // offline SiPM id
      int board;       // DIRAC board of the SiPM
      double s;        // position of the crystal centre along the line (up) [mm]
      double t;        // SiPM time (corrected with the starting offsets) [ns]
      bool use;        // false once removed by the outlier cut
    };

    // Amplitudes of the two SiPMs of one crystal of a track (photostatistics).
    struct LeftRightAmplitudes {
      int crystalId;
      float a0;        // SiPM 0 [ADC]
      float a1;        // SiPM 1 [ADC]
    };

    // Per-SiPM result of the MIP fits.
    struct MipFits {
      std::vector<double> mpv;        // MPV [ADC per 34 mm], -1 without a good fit
      std::vector<double> mpvErr;
      std::vector<double> chi2ndf;
      std::vector<int> nEntries;      // amplitude entries of the SiPM
    };

    // Per-SiPM rising edges of the amplitude spectra.
    struct Thresholds {
      std::vector<double> raw;        // raw amplitude (largest sample - baseline) [ADC], -1 if not measured
      std::vector<double> rawWidth;
      std::vector<double> fit;        // amplitude of the fit [ADC]
      std::vector<double> fitWidth;
    };

    // Per-crystal effective photoelectron scale.
    struct Photostatistics {
      std::vector<double> p0;         // ADC per effective photoelectron, -1 if not measured
      std::vector<double> p0Err;
      std::vector<double> measured;   // the measured p0 values
      double median = -1;             // median of the measured p0, used where not measured
    };
  }

  class MC2CaloCalibFromTracks : public art::EDAnalyzer {
    public:
      struct Config {
        using Name = fhicl::Name;
        using Comment = fhicl::Comment;

        // Inputs
        fhicl::Atom<art::InputTag> trkTag{Name("caloTrackCollection"), Comment("MC2CaloCosmicTracks"),
            "MC2CaloCosmicTrackFinder"};
        fhicl::Atom<art::InputTag> recoDigiTag{Name("caloRecoDigiCollection"),
            Comment("all reco digis of the event (thresholds, baselines)"), "MC2CaloRecoDigiMaker"};
        fhicl::Atom<std::string> kalSeedTag{Name("kalSeedCollection"),
            Comment("tracker KalSeeds for the alignment check; empty = none"), ""};

        // Geometry
        fhicl::Atom<double> depthInCrystal{Name("depthInCrystal"), Comment("crystal centre behind the front face [mm]"),
            100.};
        fhicl::Sequence<double> diskZOffset{Name("diskZOffset"), Comment("added to the Offline z of disk 0 and 1 [mm]"),
            std::vector<double>{0., 0.}};

        // Time calibration
        fhicl::Atom<int> nIterTime{Name("nIterTime"), Comment("iterations of the time least squares"), 40};
        fhicl::Atom<int> nIterBeforeCut{Name("nIterBeforeCut"), Comment("iterations before the outlier cut"), 10};
        fhicl::Atom<double> outlierCut{Name("outlierCut"), Comment("max |residual| after the first iterations [ns]"),
            5.};
        fhicl::Atom<int> minEntriesTime{Name("minEntriesTime"), Comment("min entries per SiPM for a time offset"), 10};

        // Energy calibration
        fhicl::Atom<double> crystalSide{Name("crystalSide"), Comment("crystal side [mm]"), 34.};
        fhicl::Atom<double> minChordFrac{Name("minChordFrac"), Comment("min chord / side"), 0.6};
        fhicl::Atom<int> minEntriesEnergy{Name("minEntriesEnergy"), Comment("min entries per SiPM for the MIP fit"),
            80};

        // Calorimeter-tracker alignment check
        fhicl::Atom<double> rMatchTracker{Name("rMatchTracker"), Comment("crystal centre - tracker line [mm]"), 60.};

        // Response for the simulation
        fhicl::Atom<int> nPed{Name("nPed"), Comment("samples for the baseline"), 4};
        fhicl::Atom<double> thrMaxADC{Name("thrMaxADC"), Comment("range of the threshold spectra [ADC]"), 400.};
        fhicl::Atom<int> minPulsesThr{Name("minPulsesThr"), Comment("min pulses per SiPM for the threshold"), 200};
        fhicl::Atom<int> minEntriesPE{Name("minEntriesPE"), Comment("min track crystals for p0"), 30};
        fhicl::Atom<double> defaultLRU{Name("defaultLRU"), Comment("LRU written in CalSimCrystals (not measured here)"),
            0.05};

        // Output files
        fhicl::Atom<std::string> timeTableFile{Name("timeTableFile"), Comment("output CalTimeCalib table"),
            "CalTimeCalib_mc2.txt"};
        fhicl::Atom<std::string> energyTableFile{Name("energyTableFile"), Comment("output CalEnergyCalib table"),
            "CalEnergyCalib_mc2.txt"};
        fhicl::Atom<std::string> statsFile{Name("statsFile"), Comment("per-SiPM statistics (csv)"),
            "calib_mc2_stats.csv"};
        fhicl::Atom<std::string> simTableFile{Name("simTableFile"), Comment("output CalSimCrystals table"),
            "CalSimCrystals_mc2.txt"};
        fhicl::Atom<std::string> baselineTableFile{Name("baselineTableFile"), Comment("output CalBaselines table"),
            "CalBaselines_mc2.txt"};
        fhicl::Atom<std::string> responseFile{Name("responseFile"),
            Comment("per-SiPM response for the simulation (csv)"), "response_mc2.csv"};
      };
      using Parameters = art::EDAnalyzer::Table<Config>;

      explicit MC2CaloCalibFromTracks(const Parameters& config);
      void beginJob() override;
      void analyze(const art::Event& event) override;
      void endJob() override;

    private:
      // analyze steps
      void readBoardMap(const art::Event& event);
      void readStartTables(const CalCalib& calCalib);
      void accumulatePulses(const art::Event& event);
      void accumulateTracks(const art::Event& event, const Calorimeter& calorimeter);
      void accumulateTrackCrystal(int trackIndex, const MC2CaloCosmicTrack& track, double inPlane, const CaloHit& hit,
                                  const Calorimeter& calorimeter, const KalSeedCollection* seeds);
      double positionAlongTrack(const MC2CaloCosmicTrack& track, double inPlane, double x, double y, double z) const;
      double pathInCrystal(const MC2CaloCosmicTrack& track, double inPlane, double x, double y);
      void storeLeftRightAmplitudes(const CaloHit& hit);
      void fillTrackerTimeDifferences(const CaloHit& hit, double x, double y, double z,
                                      const KalSeedCollection& seeds);

      // endJob steps
      std::vector<double> fitBoardOffsets() const;
      void fitSiPMOffsets(const std::vector<double>& boardOffset, std::vector<double>& sipmOffset,
                          std::vector<int>& nTimeEntries);
      MipFits fitMipPeaks();
      std::vector<double> newADC2MeV(const MipFits& mip, double& medianADC2MeV) const;
      void fitThresholdEdge(TH2F* spectra, int sipmId, double& threshold, double& width) const;
      Thresholds measureThresholds(const std::vector<double>& adc2MeV);
      Photostatistics measurePhotostatistics(const std::vector<double>& adc2MeV);
      void writeCalibrationTables(const std::vector<double>& boardOffset, const std::vector<double>& sipmOffset,
                                  const std::vector<int>& nTimeEntries, const MipFits& mip,
                                  const std::vector<double>& adc2MeV, double medianADC2MeV,
                                  int& nTimeOk, int& nEnergyOk);
      int writeResponseTables(const std::vector<double>& adc2MeV, const Thresholds& thresholds,
                              const Photostatistics& photostatistics);

      // Inputs
      art::InputTag trkTag_;
      art::InputTag recoDigiTag_;
      std::string kalSeedTag_;                        // empty: no tracker check

      // Parameters
      double depthInCrystal_;
      std::vector<double> diskZOffset_;
      int nIterTime_;
      int nIterBeforeCut_;
      double outlierCut_;
      int minEntriesTime_;
      double crystalSide_;
      double minChordFrac_;
      int minEntriesEnergy_;
      double rMatchTracker_;
      int nPed_;
      double thrMaxADC_;
      int minPulsesThr_;
      int minEntriesPE_;
      double defaultLRU_;
      std::string timeTableFile_;
      std::string energyTableFile_;
      std::string statsFile_;
      std::string simTableFile_;
      std::string baselineTableFile_;
      std::string responseFile_;

      // Conditions
      ProditionsHandle<CalCalib> calCalibHandle_;
      ProditionsHandle<CalDAQMap> daqMapHandle_;
      std::vector<int> boardOfSiPM_;                  // DIRAC board by offline SiPM id (filled at the first event)
      std::vector<double> startTimeOffset_;           // CalCalib time offset at the first event [ns]
      std::vector<double> startADC2MeV_;              // CalCalib ADC2MeV at the first event
      bool haveStartTables_ = false;

      // Accumulated over the job: all pulses
      std::vector<double> baselineSum_;               // sum of the per-pulse baselines, by SiPM
      std::vector<double> noiseVarianceSum_;          // sum of the per-pulse baseline variances, by SiPM
      std::vector<long> nPulses_;

      // Accumulated over the job: tracks
      int nTracks_ = 0;
      std::vector<TimeEntry> timeEntries_;
      std::vector<std::vector<float>> amplitudePer34mm_;   // A34 by SiPM (one-disk tracks)
      std::vector<LeftRightAmplitudes> leftRightAmplitudes_;

      // Histograms
      TH1F* hTresBefore_ = nullptr;
      TH1F* hTres_ = nullptr;
      TH1F* hOff_ = nullptr;
      TH1F* hMPV_ = nullptr;
      TH1F* hChord_ = nullptr;
      TH2F* hAmpFit_ = nullptr;
      TH2F* hAmpRaw_ = nullptr;
      TH1F* hThr_ = nullptr;
      TH1F* hThrMeV_ = nullptr;
      TH1F* hP0_ = nullptr;
      TH1F* hTrkDt_ = nullptr;
      TH2F* hTrkDtSiPM_ = nullptr;
  };

  MC2CaloCalibFromTracks::MC2CaloCalibFromTracks(const Parameters& config) :
    art::EDAnalyzer(config),
    trkTag_(config().trkTag()),
    recoDigiTag_(config().recoDigiTag()),
    kalSeedTag_(config().kalSeedTag()),
    depthInCrystal_(config().depthInCrystal()),
    diskZOffset_(config().diskZOffset()),
    nIterTime_(config().nIterTime()),
    nIterBeforeCut_(config().nIterBeforeCut()),
    outlierCut_(config().outlierCut()),
    minEntriesTime_(config().minEntriesTime()),
    crystalSide_(config().crystalSide()),
    minChordFrac_(config().minChordFrac()),
    minEntriesEnergy_(config().minEntriesEnergy()),
    rMatchTracker_(config().rMatchTracker()),
    nPed_(config().nPed()),
    thrMaxADC_(config().thrMaxADC()),
    minPulsesThr_(config().minPulsesThr()),
    minEntriesPE_(config().minEntriesPE()),
    defaultLRU_(config().defaultLRU()),
    timeTableFile_(config().timeTableFile()),
    energyTableFile_(config().energyTableFile()),
    statsFile_(config().statsFile()),
    simTableFile_(config().simTableFile()),
    baselineTableFile_(config().baselineTableFile()),
    responseFile_(config().responseFile()),
    baselineSum_(CaloConst::_nChannel, 0.),
    noiseVarianceSum_(CaloConst::_nChannel, 0.),
    nPulses_(CaloConst::_nChannel, 0),
    amplitudePer34mm_(CaloConst::_nChannel)
  {
    if (diskZOffset_.size() != 2) {
      throw cet::exception("MC2CALOCALIB") << "diskZOffset needs two values\n";
    }
    consumes<MC2CaloCosmicTrackCollection>(trkTag_);
    consumes<CaloRecoDigiCollection>(recoDigiTag_);
    if (!kalSeedTag_.empty()) {
      consumes<KalSeedCollection>(art::InputTag(kalSeedTag_));
    }
  }

  void MC2CaloCalibFromTracks::beginJob() {
    art::ServiceHandle<art::TFileService> tfs;
    const int nSiPM = CaloConst::_nChannel;
    hTresBefore_ = tfs->make<TH1F>("hTresBefore",
        "before the least squares (starting offsets);t_{SiPM} - track time [ns];SiPM entries", 400, -20, 20);
    hTres_ = tfs->make<TH1F>("hTres", "after the least squares;t_{SiPM} - track time [ns];SiPM entries",
        400, -20, 20);
    hOff_ = tfs->make<TH1F>("hOff", ";correction to the starting time offset [ns];SiPMs", 400, -20, 20);
    hMPV_ = tfs->make<TH1F>("hMPV", ";MIP MPV per SiPM [ADC per 34 mm];SiPMs", 200, 0, 2000);
    hChord_ = tfs->make<TH1F>("hChord", ";path in the crystal [mm];crystals", 120, 0, 120);
    hAmpFit_ = tfs->make<TH2F>("hAmpFit", "all pulses;offline SiPM id;amplitude of the fit [ADC]",
        nSiPM, 0, nSiPM, int(thrMaxADC_ / 2), 0, thrMaxADC_);
    hAmpRaw_ = tfs->make<TH2F>("hAmpRaw", "all pulses;offline SiPM id;largest sample - baseline [ADC]",
        nSiPM, 0, nSiPM, int(thrMaxADC_ / 2), 0, thrMaxADC_);
    hThr_ = tfs->make<TH1F>("hThr", ";threshold (raw) [ADC];SiPMs", 200, 0, thrMaxADC_);
    hThrMeV_ = tfs->make<TH1F>("hThrMeV", ";threshold [MeV];SiPMs", 200, 0, 20);
    hP0_ = tfs->make<TH1F>("hP0", ";p_{0} [ADC per effective photoelectron];crystals", 200, 0, 5);
    if (!kalSeedTag_.empty()) {
      hTrkDt_ = tfs->make<TH1F>("hTrkDt",
          "calo hits of the calo tracks vs the tracker line;t_{SiPM} - t_{line} [ns];SiPM entries",
          1000, -3000, 3000);
      hTrkDtSiPM_ = tfs->make<TH2F>("hTrkDtSiPM", ";offline SiPM id;t_{SiPM} - t_{line} [ns]",
          nSiPM, 0, nSiPM, 600, -3000, 3000);
    }
  }

  //--------------------------------------------------------------------------------------------------------------------
  // Per event
  //--------------------------------------------------------------------------------------------------------------------

  void MC2CaloCalibFromTracks::analyze(const art::Event& event) {
    const Calorimeter& calorimeter = *GeomHandle<Calorimeter>();
    const CalCalib& calCalib = calCalibHandle_.get(event.id());
    if (boardOfSiPM_.empty()) {
      readBoardMap(event);
    }
    if (!haveStartTables_) {
      readStartTables(calCalib);
    }
    accumulatePulses(event);
    accumulateTracks(event, calorimeter);
  }

  // DIRAC board of every SiPM (for the board step of the time calibration).
  void MC2CaloCalibFromTracks::readBoardMap(const art::Event& event) {
    const CalDAQMap& daqMap = daqMapHandle_.get(event.id());
    boardOfSiPM_.resize(CaloConst::_nChannel);
    for (int i = 0; i < (int)CaloConst::_nChannel; ++i) {
      boardOfSiPM_[i] = daqMap.rawId(CaloSiPMId(i)).dirac();
    }
  }

  // Starting time offsets and ADC2MeV: the CalCalib conditions of the first event.
  void MC2CaloCalibFromTracks::readStartTables(const CalCalib& calCalib) {
    startTimeOffset_.resize(CaloConst::_nChannel);
    startADC2MeV_.resize(CaloConst::_nChannel);
    for (int i = 0; i < (int)CaloConst::_nChannel; ++i) {
      startTimeOffset_[i] = calCalib.timeoffset(i);
      startADC2MeV_[i] = calCalib.ADC2MeV(i);
    }
    haveStartTables_ = true;
  }

  // All pulses of the event: baseline and noise from the first nPed samples, amplitude spectra for the thresholds.
  void MC2CaloCalibFromTracks::accumulatePulses(const art::Event& event) {
    for (const auto& recoDigi : *event.getValidHandle<CaloRecoDigiCollection>(recoDigiTag_)) {
      const int sipmId = recoDigi.SiPMID();
      if (sipmId < 0 || sipmId >= (int)CaloConst::_nChannel) {
        continue;
      }
      const auto& waveform = recoDigi.caloDigiPtr()->waveform();
      if ((int)waveform.size() <= nPed_) {
        continue;
      }
      double baseline = 0;
      for (int i = 0; i < nPed_; ++i) {
        baseline += waveform[i];
      }
      baseline /= nPed_;
      double sumSquares = 0;
      for (int i = 0; i < nPed_; ++i) {
        sumSquares += (waveform[i] - baseline) * (waveform[i] - baseline);
      }
      baselineSum_[sipmId] += baseline;
      noiseVarianceSum_[sipmId] += sumSquares / nPed_;
      ++nPulses_[sipmId];

      const double largestSample = *std::max_element(waveform.begin(), waveform.end());
      hAmpRaw_->Fill(sipmId + 0.5, largestSample - baseline);
      if (startADC2MeV_[sipmId] > 0) {
        hAmpFit_->Fill(sipmId + 0.5, recoDigi.energyDep() / startADC2MeV_[sipmId]);
      }
    }
  }

  // Every crystal of every calorimeter track (and the tracker lines, if requested).
  void MC2CaloCalibFromTracks::accumulateTracks(const art::Event& event, const Calorimeter& calorimeter) {
    const auto& tracks = *event.getValidHandle<MC2CaloCosmicTrackCollection>(trkTag_);
    const KalSeedCollection* seeds = nullptr;
    if (!kalSeedTag_.empty()) {
      art::Handle<KalSeedCollection> seedsHandle;
      if (event.getByLabel(art::InputTag(kalSeedTag_), seedsHandle)) {
        seeds = seedsHandle.product();
      }
    }
    for (const auto& track : tracks) {
      const int trackIndex = nTracks_++;
      const double inPlane = std::hypot(track.ux, track.uy);   // length of the x-y projection of the direction
      for (const auto& hit : track.hits) {
        accumulateTrackCrystal(trackIndex, track, inPlane, *hit, calorimeter, seeds);
      }
    }
  }

  void MC2CaloCalibFromTracks::accumulateTrackCrystal(int trackIndex, const MC2CaloCosmicTrack& track, double inPlane,
                                                      const CaloHit& hit, const Calorimeter& calorimeter,
                                                      const KalSeedCollection* seeds) {
    // Crystal centre in the tracker frame (position() is the front face) with the disk z offsets
    const auto& crystal = calorimeter.crystal(hit.crystalID());
    const int disk = crystal.diskID();
    const auto position = calorimeter.mu2eToTracker(crystal.position());
    const double crystalZ = position.z() + depthInCrystal_ + diskZOffset_[disk];

    const double s = positionAlongTrack(track, inPlane, position.x(), position.y(), crystalZ);
    const double path = pathInCrystal(track, inPlane, position.x(), position.y());

    if (hit.recoCaloDigis().size() == 2) {
      storeLeftRightAmplitudes(hit);
    }

    // Time entry of every SiPM; amplitude per 34 mm for the MIP fit (one-disk tracks with a good chord)
    const double lossSide = mostProbableLoss(crystalSide_);
    for (const auto& recoDigi : hit.recoCaloDigis()) {
      const int sipmId = recoDigi->SiPMID();
      if (sipmId < 0 || sipmId >= (int)CaloConst::_nChannel) {
        continue;
      }
      timeEntries_.push_back({trackIndex, sipmId, boardOfSiPM_[sipmId], s, recoDigi->time(), true});
      if (path > 0 && track.type < 3 && startADC2MeV_[sipmId] > 0) {
        amplitudePer34mm_[sipmId].push_back(
            recoDigi->energyDep() / startADC2MeV_[sipmId] * lossSide / mostProbableLoss(path));
      }
    }

    if (seeds) {
      fillTrackerTimeDifferences(hit, position.x(), position.y(), crystalZ, *seeds);
    }
  }

  // Position of the crystal centre along the line, positive up. Two-disk tracks: 3D projection; one-disk tracks:
  // projection in the disk plane, normalised to the x-y length of the direction.
  double MC2CaloCalibFromTracks::positionAlongTrack(const MC2CaloCosmicTrack& track, double inPlane,
                                                    double x, double y, double z) const {
    if (track.type == 3) {
      return (x - track.cx) * track.ux + (y - track.cy) * track.uy + (z - track.cz) * track.uz;
    }
    return ((x - track.cx) * track.ux + (y - track.cy) * track.uy) / (inPlane > 0 ? inPlane : 1.);
  }

  // Path of the line in the crystal: chord in the x-y cross-section, divided by sqrt(1 - uz^2) for two-disk tracks.
  // 0 for chords shorter than minChordFrac * side; the accepted paths are filled in hChord.
  double MC2CaloCalibFromTracks::pathInCrystal(const MC2CaloCosmicTrack& track, double inPlane, double x, double y) {
    double path = 0;
    if (inPlane > 0) {
      const double chordXY = chord(track.cx, track.cy, track.ux / inPlane, track.uy / inPlane, x, y, crystalSide_);
      path = (track.type == 3) ? chordXY / std::sqrt(std::max(1e-6, 1. - double(track.uz) * track.uz)) : chordXY;
      if (chordXY < minChordFrac_ * crystalSide_) {
        path = 0;
      } else {
        hChord_->Fill(path);
      }
    }
    return path;
  }

  // Photostatistics: amplitudes of the two SiPMs of the crystal, kept if both are positive.
  void MC2CaloCalibFromTracks::storeLeftRightAmplitudes(const CaloHit& hit) {
    float amplitude[2] = {-1, -1};
    for (const auto& recoDigi : hit.recoCaloDigis()) {
      const int sipmId = recoDigi->SiPMID();
      if (sipmId >= 0 && sipmId < (int)CaloConst::_nChannel && startADC2MeV_[sipmId] > 0) {
        amplitude[sipmId % 2] = recoDigi->energyDep() / startADC2MeV_[sipmId];
      }
    }
    if (amplitude[0] > 0 && amplitude[1] > 0) {
      leftRightAmplitudes_.push_back({hit.crystalID(), amplitude[0], amplitude[1]});
    }
  }

  // Calorimeter-tracker alignment check: for every tracker line passing within rMatchTracker of the crystal centre,
  // SiPM time - time of the line at the point closest to the crystal centre.
  void MC2CaloCalibFromTracks::fillTrackerTimeDifferences(const CaloHit& hit, double x, double y, double z,
                                                          const KalSeedCollection& seeds) {
    for (const auto& seed : seeds) {
      if (!seed.kinematicLineFit() || seed.segments().empty()) {
        continue;
      }
      // the segment with the largest z
      const KalSegment* segment = nullptr;
      double segmentZ = -1e30;
      for (const auto& candidate : seed.segments()) {
        const double candidateZ = candidate.position3().Z();
        if (candidateZ > segmentZ) {
          segmentZ = candidateZ;
          segment = &candidate;
        }
      }
      const auto line = segment->kinematicLine();
      const double tRef = 0.5 * (segment->tmin() + segment->tmax());
      const auto lineRef = line.position3(tRef);
      const auto velocity = line.velocity(tRef);
      const double velocity2 = velocity.Mag2();
      if (velocity2 <= 0) {
        continue;
      }
      // u = time from tRef to the point of the line closest to the crystal centre; e = distance vector
      const double dx = x - lineRef.X();
      const double dy = y - lineRef.Y();
      const double dz = z - lineRef.Z();
      const double u = (dx * velocity.X() + dy * velocity.Y() + dz * velocity.Z()) / velocity2;
      const double ex = dx - velocity.X() * u;
      const double ey = dy - velocity.Y() * u;
      const double ez = dz - velocity.Z() * u;
      if (std::sqrt(ex * ex + ey * ey + ez * ez) > rMatchTracker_) {
        continue;
      }
      for (const auto& recoDigi : hit.recoCaloDigis()) {
        const double dt = recoDigi->time() - (tRef + u);
        hTrkDt_->Fill(dt);
        hTrkDtSiPM_->Fill(recoDigi->SiPMID() + 0.5, dt);
      }
    }
  }

  //--------------------------------------------------------------------------------------------------------------------
  // End of the job: fits and tables
  //--------------------------------------------------------------------------------------------------------------------

  void MC2CaloCalibFromTracks::endJob() {
    // Time
    const std::vector<double> boardOffset = fitBoardOffsets();
    std::vector<double> sipmOffset;
    std::vector<int> nTimeEntries;
    fitSiPMOffsets(boardOffset, sipmOffset, nTimeEntries);

    // Energy
    const MipFits mip = fitMipPeaks();
    double medianADC2MeV = 0.;
    const std::vector<double> adc2MeV = newADC2MeV(mip, medianADC2MeV);

    // Response for the simulation
    const Thresholds thresholds = measureThresholds(adc2MeV);
    const Photostatistics photostatistics = measurePhotostatistics(adc2MeV);

    // Tables
    int nTimeOk = 0;
    int nEnergyOk = 0;
    writeCalibrationTables(boardOffset, sipmOffset, nTimeEntries, mip, adc2MeV, medianADC2MeV, nTimeOk, nEnergyOk);
    const int nThresholds = writeResponseTables(adc2MeV, thresholds, photostatistics);

    printf("MC2CaloCalibFromTracks: thresholds for %d SiPMs (median %.1f ADC raw); p0 measured for %zu crystals, "
           "median %.3f ADC per effective photoelectron -> %s, %s, %s\n",
           nThresholds, hThr_->GetMean(), photostatistics.measured.size(), photostatistics.median,
           simTableFile_.c_str(), baselineTableFile_.c_str(), responseFile_.c_str());
    printf("MC2CaloCalibFromTracks: %d tracks, %zu SiPM time entries; time offsets for %d SiPMs, "
           "MIP fits for %d SiPMs -> %s, %s, %s\n",
           nTracks_, timeEntries_.size(), nTimeOk, nEnergyOk,
           timeTableFile_.c_str(), energyTableFile_.c_str(), statsFile_.c_str());
    printf("MC2CaloCalibFromTracks: residual rms before %.3f ns, after %.3f ns (within +-20 ns)\n",
           hTresBefore_->GetRMS(), hTres_->GetRMS());
  }

  // Time, step 1: one offset per board, robust to jumps (see the header). Returns the offsets by board, mean 0.
  std::vector<double> MC2CaloCalibFromTracks::fitBoardOffsets() const {
    const int nBoards = CaloConst::_nDIRAC;
    const int nBins = 3000;   // residual histogram: 2 ns bins over [-3000, 3000] ns
    std::vector<double> boardOffset(nBoards, 0.);

    std::vector<std::vector<size_t>> entriesOfTrack(nTracks_);
    for (size_t i = 0; i < timeEntries_.size(); ++i) {
      entriesOfTrack[timeEntries_[i].track].push_back(i);
    }
    std::vector<double> trackTime(nTracks_, 0.);

    for (int iteration = 0; iteration < 3; ++iteration) {
      // track time: median of its SiPMs (upper median), with the current board offsets
      for (int k = 0; k < nTracks_; ++k) {
        std::vector<double> times;
        for (size_t i : entriesOfTrack[k]) {
          times.push_back(timeEntries_[i].t - boardOffset[timeEntries_[i].board] + timeEntries_[i].s / kC);
        }
        if (times.empty()) {
          continue;
        }
        std::nth_element(times.begin(), times.begin() + times.size() / 2, times.end());
        trackTime[k] = times[times.size() / 2];
      }

      // residuals per board, tracks with at least 4 SiPMs
      std::vector<std::vector<int>> residuals(nBoards, std::vector<int>(nBins, 0));
      for (auto& entry : timeEntries_) {
        if (entriesOfTrack[entry.track].size() < 4) {
          continue;
        }
        const double residual = entry.t - boardOffset[entry.board] - (trackTime[entry.track] - entry.s / kC);
        const int bin = int((residual + 3000.) / 2.);
        if (bin >= 0 && bin < nBins) {
          residuals[entry.board][bin]++;
        }
      }

      // per board: mode (most populated 5-bin window), then mean within +-5 bins of it
      double offsetSum = 0;
      int nFitted = 0;
      for (int board = 0; board < nBoards; ++board) {
        const std::vector<int>& histogram = residuals[board];
        int modeBin = -1;
        int modeContent = 0;
        for (int b = 0; b + 4 < nBins; ++b) {
          int windowContent = 0;
          for (int q = 0; q < 5; ++q) {
            windowContent += histogram[b + q];
          }
          if (windowContent > modeContent) {
            modeContent = windowContent;
            modeBin = b + 2;
          }
        }
        if (modeContent < 10) {
          continue;
        }
        double correction = -3000. + 2. * modeBin + 1.;
        double sumWeights = 0;
        double sumX = 0;
        for (int b = std::max(0, modeBin - 5); b <= std::min(2999, modeBin + 5); ++b) {
          const double x = -3000. + 2. * b + 1.;
          sumWeights += histogram[b];
          sumX += histogram[b] * x;
        }
        if (sumWeights > 0) {
          correction = sumX / sumWeights;
        }
        boardOffset[board] += correction;
        offsetSum += boardOffset[board];
        ++nFitted;
      }

      // mean board offset (over the fitted boards) fixed to 0
      if (nFitted) {
        for (int board = 0; board < nBoards; ++board) {
          boardOffset[board] -= offsetSum / nFitted;
        }
      }
    }

    int nJumped = 0;
    for (int board = 0; board < nBoards; ++board) {
      if (std::abs(boardOffset[board]) > 20.) {
        ++nJumped;
        printf("MC2CaloCalibFromTracks: board %d jumped by %.1f ns with respect to the start tables\n",
               board, boardOffset[board]);
      }
    }
    printf("MC2CaloCalibFromTracks: board step: %d boards with a correction > 20 ns\n", nJumped);
    return boardOffset;
  }

  // Time, step 2: alternating least squares for the track times T_k and the SiPM offsets o_j, starting from the board
  // offsets, with the progressive outlier cut. Fills hTresBefore and hTres.
  void MC2CaloCalibFromTracks::fitSiPMOffsets(const std::vector<double>& boardOffset, std::vector<double>& sipmOffset,
                                              std::vector<int>& nTimeEntries) {
    const int nSiPM = CaloConst::_nChannel;
    std::vector<double>& o = sipmOffset;
    std::vector<int>& nO = nTimeEntries;
    o.assign(nSiPM, 0.);
    nO.assign(nSiPM, 0);
    std::vector<double> T(nTracks_, 0.);
    std::vector<int> nT(nTracks_);
    std::vector<double> sumT(nTracks_);
    std::vector<double> sumO(nSiPM);

    for (int j = 0; j < nSiPM; ++j) {
      o[j] = boardOffset[boardOfSiPM_.empty() ? 0 : boardOfSiPM_[j]];
    }

    // T_k = mean over the used entries of t + s / c - o_j
    auto updateTrackTimes = [&]() {
      std::fill(sumT.begin(), sumT.end(), 0.);
      std::fill(nT.begin(), nT.end(), 0);
      for (auto& entry : timeEntries_) {
        if (entry.use) {
          sumT[entry.track] += entry.t - o[entry.sipm] + entry.s / kC;
          ++nT[entry.track];
        }
      }
      for (int k = 0; k < nTracks_; ++k) {
        T[k] = nT[k] ? sumT[k] / nT[k] : 0;
      }
    };

    updateTrackTimes();
    for (auto& entry : timeEntries_) {
      hTresBefore_->Fill(entry.t - (T[entry.track] - entry.s / kC));
    }

    for (int iteration = 0; iteration < nIterTime_; ++iteration) {
      // progressive outlier cut: 200, 67, 22, 7.4 ns, then outlierCut
      if (iteration >= nIterBeforeCut_ && (iteration - nIterBeforeCut_) % 5 == 0) {
        const double cutNow = std::max(outlierCut_, 200. / std::pow(3., (iteration - nIterBeforeCut_) / 5));
        for (auto& entry : timeEntries_) {
          if (std::abs(entry.t - (T[entry.track] - entry.s / kC) - o[entry.sipm]) > cutNow) {
            entry.use = false;
          }
        }
      }

      // o_j = mean over the used entries (of tracks with more than one entry) of t - (T_k - s / c)
      std::fill(sumO.begin(), sumO.end(), 0.);
      std::fill(nO.begin(), nO.end(), 0);
      for (auto& entry : timeEntries_) {
        if (entry.use && nT[entry.track] > 1) {
          sumO[entry.sipm] += entry.t - (T[entry.track] - entry.s / kC);
          ++nO[entry.sipm];
        }
      }
      double meanOffset = 0;
      int nMean = 0;
      for (int j = 0; j < nSiPM; ++j) {
        if (nO[j] >= minEntriesTime_) {
          o[j] = sumO[j] / nO[j];
          meanOffset += o[j];
          ++nMean;
        } else {
          o[j] = boardOffset[boardOfSiPM_[j]];   // few entries: board offset
        }
      }
      // mean of the fitted o_j fixed to 0
      if (nMean) {
        meanOffset /= nMean;
        for (int j = 0; j < nSiPM; ++j) {
          if (nO[j] >= minEntriesTime_) {
            o[j] -= meanOffset;
          }
        }
      }
      updateTrackTimes();
    }

    for (auto& entry : timeEntries_) {
      if (entry.use) {
        hTres_->Fill(entry.t - (T[entry.track] - entry.s / kC) - o[entry.sipm]);
      }
    }
  }

  // Energy: Landau (x) Gauss fit of the amplitude per 34 mm of every SiPM. Fills hMPV.
  MipFits MC2CaloCalibFromTracks::fitMipPeaks() {
    const int nSiPM = CaloConst::_nChannel;
    MipFits mip;
    mip.mpv.assign(nSiPM, -1);
    mip.mpvErr.assign(nSiPM, -1);
    mip.chi2ndf.assign(nSiPM, -1);
    mip.nEntries.assign(nSiPM, 0);
    for (int j = 0; j < nSiPM; ++j) {
      auto& amplitudes = amplitudePer34mm_[j];
      mip.nEntries[j] = amplitudes.size();
      if ((int)amplitudes.size() < minEntriesEnergy_) {
        continue;
      }
      std::vector<float> sorted(amplitudes);
      std::nth_element(sorted.begin(), sorted.begin() + sorted.size() / 2, sorted.end());
      const double median = sorted[sorted.size() / 2];
      if (median <= 0) {
        continue;
      }

      TH1F histogram("hmip", "", 120, 0, 3 * median);
      histogram.SetDirectory(nullptr);
      for (float x : amplitudes) {
        histogram.Fill(x);
      }
      const double peak = histogram.GetBinCenter(histogram.GetMaximumBin());

      TF1 fit("mc2lg", langaus, 0.5 * peak, 2.5 * peak, 4);
      fit.SetParameters(0.08 * peak, peak, histogram.Integral("width"), 0.08 * peak);
      fit.SetParLimits(0, 0.01 * peak, 0.5 * peak);
      fit.SetParLimits(1, 0.5 * peak, 1.8 * peak);
      fit.SetParLimits(3, 0.005 * peak, 0.6 * peak);
      histogram.Fit(&fit, "QRN0");

      // accepted if the MPV error is below 10 %
      if (fit.GetParameter(1) > 0 && fit.GetParError(1) > 0 && fit.GetParError(1) < 0.1 * fit.GetParameter(1)) {
        mip.mpv[j] = fit.GetParameter(1);
        mip.mpvErr[j] = fit.GetParError(1);
        mip.chi2ndf[j] = fit.GetNDF() > 0 ? fit.GetChisquare() / fit.GetNDF() : -1;
        hMPV_->Fill(mip.mpv[j]);
      }
    }
    return mip;
  }

  // New ADC2MeV per SiPM = dp(34 mm) / MPV (used below for MeV and pePerMeV). SiPMs without a good MIP fit: median
  // ADC2MeV of the fitted SiPMs (as mc2tools energy_calibration.py), else the starting value.
  std::vector<double> MC2CaloCalibFromTracks::newADC2MeV(const MipFits& mip, double& medianADC2MeV) const {
    const int nSiPM = CaloConst::_nChannel;
    const double lossSide = mostProbableLoss(crystalSide_);
    std::vector<double> fitted;
    for (int j = 0; j < nSiPM; ++j) {
      if (mip.mpv[j] > 0) {
        fitted.push_back(lossSide / mip.mpv[j]);
      }
    }
    medianADC2MeV = 0.;
    if (!fitted.empty()) {
      std::nth_element(fitted.begin(), fitted.begin() + fitted.size() / 2, fitted.end());
      medianADC2MeV = fitted[fitted.size() / 2];
    }
    std::vector<double> adc2MeV(nSiPM);
    for (int j = 0; j < nSiPM; ++j) {
      adc2MeV[j] = mip.mpv[j] > 0 ? lossSide / mip.mpv[j]
                                  : (medianADC2MeV > 0 ? medianADC2MeV : (haveStartTables_ ? startADC2MeV_[j] : 0.));
    }
    return adc2MeV;
  }

  // Rising edge of the amplitude spectrum of one SiPM: half-maximum crossing, then the erf x exponential fit around
  // it. threshold = -1 with fewer than minPulsesThr pulses; width = -1 if the fit is not used.
  void MC2CaloCalibFromTracks::fitThresholdEdge(TH2F* spectra, int sipmId, double& threshold, double& width) const {
    threshold = -1;
    width = -1;
    TH1D* spectrum = spectra->ProjectionY("pj", sipmId + 1, sipmId + 1);
    if (spectrum->Integral() < minPulsesThr_) {
      delete spectrum;
      return;
    }
    // half-maximum crossing, linear interpolation between the two bins around it
    const int maxBin = spectrum->GetMaximumBin();
    const double maxContent = spectrum->GetBinContent(maxBin);
    int k = 1;
    while (k < maxBin && spectrum->GetBinContent(k) < 0.5 * maxContent) {
      ++k;
    }
    const double c1 = spectrum->GetBinContent(std::max(1, k - 1));
    const double c2 = spectrum->GetBinContent(k);
    const double x1 = spectrum->GetBinCenter(std::max(1, k - 1));
    const double x2 = spectrum->GetBinCenter(k);
    const double halfMax = (c2 > c1) ? x1 + (0.5 * maxContent - c1) / (c2 - c1) * (x2 - x1) : x2;
    threshold = halfMax;

    // edge fit: used if it converges within 20 ADC of the half-maximum crossing
    TF1 edge("edge", "[0]*0.5*(1+TMath::Erf((x-[1])/(sqrt(2.)*[2])))*exp(-(x-[1])/[3])",
             std::max(0., halfMax - 30), std::min(thrMaxADC_, halfMax + 80));
    edge.SetParameters(maxContent, halfMax, 3., 40.);
    edge.SetParLimits(2, 0.2, 30.);
    edge.SetParLimits(3, 2., 1e4);
    if (spectrum->Fit(&edge, "QRN0") == 0 && std::abs(edge.GetParameter(1) - halfMax) < 20) {
      threshold = edge.GetParameter(1);
      width = edge.GetParameter(2);
    }
    delete spectrum;
  }

  // Thresholds of all SiPMs, raw and fit amplitude. Fills hThr (raw) and hThrMeV (fit, with the new ADC2MeV).
  Thresholds MC2CaloCalibFromTracks::measureThresholds(const std::vector<double>& adc2MeV) {
    const int nSiPM = CaloConst::_nChannel;
    Thresholds thresholds;
    thresholds.raw.assign(nSiPM, -1);
    thresholds.rawWidth.assign(nSiPM, -1);
    thresholds.fit.assign(nSiPM, -1);
    thresholds.fitWidth.assign(nSiPM, -1);
    for (int j = 0; j < nSiPM; ++j) {
      fitThresholdEdge(hAmpRaw_, j, thresholds.raw[j], thresholds.rawWidth[j]);
      fitThresholdEdge(hAmpFit_, j, thresholds.fit[j], thresholds.fitWidth[j]);
      if (thresholds.raw[j] > 0) {
        hThr_->Fill(thresholds.raw[j]);
      }
      if (thresholds.fit[j] > 0 && adc2MeV[j] > 0) {
        hThrMeV_->Fill(thresholds.fit[j] * adc2MeV[j]);
      }
    }
    return thresholds;
  }

  // Effective photoelectrons per crystal from the left/right amplitudes of the track crystals (see the header).
  // Fills hP0.
  Photostatistics MC2CaloCalibFromTracks::measurePhotostatistics(const std::vector<double>& adc2MeV) {
    const int nCrystals = CaloConst::_nCrystal;
    std::vector<std::vector<std::pair<float, float>>> amplitudesOfCrystal(nCrystals);
    for (auto& entry : leftRightAmplitudes_) {
      if (entry.crystalId >= 0 && entry.crystalId < nCrystals) {
        amplitudesOfCrystal[entry.crystalId].push_back({entry.a0, entry.a1});
      }
    }

    Photostatistics result;
    result.p0.assign(nCrystals, -1);
    result.p0Err.assign(nCrystals, -1);
    for (int c = 0; c < nCrystals; ++c) {
      auto& amplitudes = amplitudesOfCrystal[c];
      if ((int)amplitudes.size() < minEntriesPE_) {
        continue;
      }
      const int j0 = 2 * c;
      const int j1 = 2 * c + 1;
      if (adc2MeV[j0] <= 0 || adc2MeV[j1] <= 0) {
        continue;
      }
      // amplitudes brought to the mean gain g of the two SiPMs
      const double g0 = 1. / adc2MeV[j0];
      const double g1 = 1. / adc2MeV[j1];
      const double g = 0.5 * (g0 + g1);
      std::vector<double> difference;
      std::vector<double> sum;
      for (auto& pair : amplitudes) {
        const double a0 = pair.first * g / g0;
        const double a1 = pair.second * g / g1;
        difference.push_back(a0 - a1);
        sum.push_back(a0 + a1);
      }

      // robust variance of L - R: (1.4826 x median absolute deviation)^2
      std::vector<double> sorted(difference);
      std::nth_element(sorted.begin(), sorted.begin() + sorted.size() / 2, sorted.end());
      const double medianDifference = sorted[sorted.size() / 2];
      std::vector<double> absDeviation;
      for (double x : difference) {
        absDeviation.push_back(std::abs(x - medianDifference));
      }
      std::nth_element(absDeviation.begin(), absDeviation.begin() + absDeviation.size() / 2, absDeviation.end());
      const double variance = std::pow(1.4826 * absDeviation[absDeviation.size() / 2], 2);

      double meanSum = 0;
      for (double x : sum) {
        meanSum += x;
      }
      meanSum /= sum.size();

      // noise variance of the two SiPMs
      const double noiseVariance = (nPulses_[j0] ? noiseVarianceSum_[j0] / nPulses_[j0] : 0)
                                 + (nPulses_[j1] ? noiseVarianceSum_[j1] / nPulses_[j1] : 0);
      if (meanSum > 0 && variance > noiseVariance) {
        result.p0[c] = (variance - noiseVariance) / meanSum;
        result.p0Err[c] = result.p0[c] * std::sqrt(2. / amplitudes.size());
        result.measured.push_back(result.p0[c]);
        hP0_->Fill(result.p0[c]);
      }
    }

    result.median = -1;
    if (!result.measured.empty()) {
      std::vector<double>& measured = result.measured;
      std::nth_element(measured.begin(), measured.begin() + measured.size() / 2, measured.end());
      result.median = measured[measured.size() / 2];
    }
    return result;
  }

  // CalTimeCalib, CalEnergyCalib and the stats csv. Fills hOff; counts the SiPMs with a time offset and a MIP fit.
  void MC2CaloCalibFromTracks::writeCalibrationTables(const std::vector<double>& boardOffset,
                                                      const std::vector<double>& sipmOffset,
                                                      const std::vector<int>& nTimeEntries, const MipFits& mip,
                                                      const std::vector<double>& adc2MeV, double medianADC2MeV,
                                                      int& nTimeOk, int& nEnergyOk) {
    const int nSiPM = CaloConst::_nChannel;
    const double lossSide = mostProbableLoss(crystalSide_);
    FILE* timeTable = fopen(timeTableFile_.c_str(), "w");
    FILE* energyTable = fopen(energyTableFile_.c_str(), "w");
    FILE* stats = fopen(statsFile_.c_str(), "w");
    if (!timeTable || !energyTable || !stats) {
      throw cet::exception("MC2CALOCALIB") << "cannot write the output tables\n";
    }
    fprintf(timeTable, "TABLE CalTimeCalib\n#roid,timeOffset  MC2CaloCalibFromTracks: starting tables + least squares "
            "on %d calorimeter cosmic tracks\n", nTracks_);
    fprintf(energyTable, "TABLE CalEnergyCalib\n#roid,ADC2MeV  MC2CaloCalibFromTracks: MIP MPV per 34 mm "
            "(Landau x Gauss), ADC2MeV = %.4f MeV / MPV\n", lossSide);
    fprintf(stats, "roid,nTime,timeCorrection_ns,startTimeOffset_ns,newTimeOffset_ns,timeFromTracks,nEnergy,mpv_adc,"
            "mpv_err,chi2ndf,startADC2MeV,newADC2MeV,energyFromTracks\n");

    nTimeOk = 0;
    nEnergyOk = 0;
    for (int j = 0; j < nSiPM; ++j) {
      const bool timeOk = nTimeEntries[j] >= minEntriesTime_;
      const bool energyOk = mip.mpv[j] > 0;
      // SiPMs without enough entries: correction of their board (jumps)
      const double correction = timeOk ? sipmOffset[j] : boardOffset[boardOfSiPM_[j]];
      const double newTime = haveStartTables_ ? startTimeOffset_[j] + correction : correction;
      const double startEnergy = haveStartTables_ ? startADC2MeV_[j] : 0;
      const double newEnergy = adc2MeV[j];
      if (timeOk) {
        hOff_->Fill(sipmOffset[j]);
        ++nTimeOk;
      }
      if (energyOk) {
        ++nEnergyOk;
      }
      fprintf(timeTable, "%d,%.4f\n", j, newTime);
      fprintf(energyTable, "%d,%.6f\n", j, newEnergy);
      fprintf(stats, "%d,%d,%.4f,%.4f,%.4f,%d,%d,%.2f,%.2f,%.3f,%.6f,%.6f,%d\n",
              j, nTimeEntries[j], correction, haveStartTables_ ? startTimeOffset_[j] : 0., newTime, (int)timeOk,
              mip.nEntries[j], mip.mpv[j], mip.mpvErr[j], mip.chi2ndf[j], startEnergy, newEnergy, (int)energyOk);
    }
    // spare channels of the DB tables
    for (int j = nSiPM; j < (int)CaloConst::_nChannelDB; ++j) {
      fprintf(timeTable, "%d,0.0000\n", j);
      fprintf(energyTable, "%d,%.6f\n", j, medianADC2MeV);
    }
    fclose(timeTable);
    fclose(energyTable);
    fclose(stats);
  }

  // CalSimCrystals, CalBaselines and the response csv. Returns the number of SiPMs with a fit-amplitude threshold.
  int MC2CaloCalibFromTracks::writeResponseTables(const std::vector<double>& adc2MeV, const Thresholds& thresholds,
                                                  const Photostatistics& photostatistics) {
    const int nSiPM = CaloConst::_nChannel;
    const int nCrystals = CaloConst::_nCrystal;
    const std::vector<double>& p0 = photostatistics.p0;
    const double p0Median = photostatistics.median;
    FILE* simTable = fopen(simTableFile_.c_str(), "w");
    FILE* baselineTable = fopen(baselineTableFile_.c_str(), "w");
    FILE* response = fopen(responseFile_.c_str(), "w");
    if (!simTable || !baselineTable || !response) {
      throw cet::exception("MC2CALOCALIB") << "cannot write the response tables\n";
    }

    // CalSimCrystals, per crystal: ADCPerMeV = 1 / ADC2MeV, pePerMeV = ADCPerMeV / p0
    fprintf(simTable, "TABLE CalSimCrystals\n#crid,LRU,pePerMeV0,pePerMeV1,ADCPerMeV0,ADCPerMeV1  "
            "MC2CaloCalibFromTracks: pePerMeV = ADCPerMeV / p0 (L-R photostatistics; median p0 %.3f ADC/pe where not "
            "measured), LRU = %.3f not measured\n", p0Median, defaultLRU_);
    for (int c = 0; c < nCrystals; ++c) {
      const double p0Used = p0[c] > 0 ? p0[c] : p0Median;
      const double adcPerMeV0 = adc2MeV[2 * c] > 0 ? 1. / adc2MeV[2 * c] : 0;
      const double adcPerMeV1 = adc2MeV[2 * c + 1] > 0 ? 1. / adc2MeV[2 * c + 1] : 0;
      fprintf(simTable, "%d,%.4f,%.3f,%.3f,%.4f,%.4f\n", c, defaultLRU_,
              p0Used > 0 ? adcPerMeV0 / p0Used : 0., p0Used > 0 ? adcPerMeV1 / p0Used : 0., adcPerMeV0, adcPerMeV1);
    }

    // CalBaselines and response csv, per SiPM
    fprintf(baselineTable, "TABLE CalBaselines\n#roid,baseline,threshold  MC2CaloCalibFromTracks: baseline = mean of "
            "the first %d samples [ADC], threshold = rising edge of the raw amplitude (largest sample - baseline) "
            "[ADC]\n", nPed_);
    fprintf(response, "roid,crid,nPulses,baseline_adc,noise_adc,thrRaw_adc,thrRawWidth_adc,thrFit_adc,thrFitWidth_adc,"
            "ADC2MeV,thr_MeV,thrWidth_MeV,p0_adc_per_pe,p0_err,p0_measured,pePerMeV\n");
    int nThresholds = 0;
    for (int j = 0; j < nSiPM; ++j) {
      const double baseline = nPulses_[j] ? baselineSum_[j] / nPulses_[j] : 0;
      const double noise = nPulses_[j] ? std::sqrt(noiseVarianceSum_[j] / nPulses_[j]) : 0;
      const int c = j / 2;
      fprintf(baselineTable, "%d,%.3f,%.3f\n", j, baseline, thresholds.raw[j] > 0 ? thresholds.raw[j] : 0.);
      if (thresholds.fit[j] > 0) {
        ++nThresholds;
      }
      const double p0Used = (c < nCrystals && p0[c] > 0) ? p0[c] : p0Median;
      fprintf(response, "%d,%d,%ld,%.3f,%.3f,%.2f,%.2f,%.2f,%.2f,%.6f,%.3f,%.3f,%.4f,%.4f,%d,%.3f\n",
              j, c, nPulses_[j], baseline, noise,
              thresholds.raw[j], thresholds.rawWidth[j], thresholds.fit[j], thresholds.fitWidth[j], adc2MeV[j],
              thresholds.fit[j] > 0 ? thresholds.fit[j] * adc2MeV[j] : -1.,
              thresholds.fitWidth[j] > 0 ? thresholds.fitWidth[j] * adc2MeV[j] : -1.,
              p0Used, (c < nCrystals) ? photostatistics.p0Err[c] : -1., (int)(c < nCrystals && p0[c] > 0),
              (p0Used > 0 && adc2MeV[j] > 0) ? 1. / adc2MeV[j] / p0Used : 0.);
    }
    // spare channels of the DB table
    for (int j = nSiPM; j < (int)CaloConst::_nChannelDB; ++j) {
      fprintf(baselineTable, "%d,0.000,0.000\n", j);
    }
    fclose(simTable);
    fclose(baselineTable);
    fclose(response);
    return nThresholds;
  }
}

DEFINE_ART_MODULE(mu2e::MC2CaloCalibFromTracks)
