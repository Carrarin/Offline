//
// MC2RecoCompare: performance of several reconstructions of the same calorimeter digis, measured on the same cosmic
// tracks.
//
// Every reconstruction (method) is a CaloRecoDigi collection made in the same job from the same CaloDigis (Offline
// CaloRecoDigiMaker, MC2CaloRecoDigiMaker with different templates or fit windows, ...). The tracks come from one
// reference reconstruction (MC2CaloCosmicTrackFinder): every method is measured on the same crystals, none chooses its
// own tracks. A method's time of a SiPM is the CaloRecoDigi of that method made from the same CaloDigi (the one with
// the largest energy if more than one). The job is written by MC2Reco/scripts/make_reco_compare_fcl.py.
//
// Inputs:
//   caloTrackCollection     MC2CaloCosmicTrackCollection of the reference (its hits give the crystals and CaloDigis)
//   methods[i].caloRecoDigiCollection   CaloRecoDigiCollection of each method
//
// Common variables (from the raw waveform of the CaloDigi, the same for all methods):
//   amplitude = largest sample - pedestal (mean of the first nPed samples) [ADC]; saturated = largest sample >= satADC;
//   phase     = 50 % crossing of the rising edge (linear interpolation) modulo 5 ns: where the pulse falls in the
//               sample.
// Crystals of a track are used if the reference has both SiPMs and neither is saturated.
//
// Protocol (the same for every method):
//   track line (as the mc2tools track check, docs lr_vs_tracks):
//     two-disk tracks t = T - s/c (3D direction known, only T free);
//     one-disk tracks t = T + k s with T and k free (s = position along the direction in the disk plane: the
//     inclination in z, hence the slope, is unknown).
//   even events (track parity = event number % 2 == 0): per-SiPM time offsets of the method, global least squares on
//     the tracks t_kj = line_k(s_kj) + o_j (alternating least squares: line of every track, then offsets; mean offset
//     0, progressive outlier cut);
//     walkCorrection: t_kj = line_k(s_kj) + o_j + w(A_kj), w a table in ln of the SiPM raw amplitude A (walkBins
//     logarithmic bins from 150 to 2100 ADC, linear interpolation between bin centres, constant outside, mean 0 over
//     the SiPMs), estimated in the same alternating least squares; then subtracted from every SiPM time on the odd
//     events too, before the line fit.
//   odd events: figures of merit with those offsets. The even events give the same figures, used only to choose the
//     fit window within a family (the method of a family with the smallest track residual on the even events is
//     marked).
//
// Figures of merit:
//   track residual  crystal time T = (tL + tR)/2 (offsets subtracted) minus the line fitted to the crystals of the
//                   track (crystals farther than maxDevOthers from the line are left out of the fit, iteratively; at
//                   least minCrystals in the fit), studentised: divided by sqrt(1 - h) for a crystal in the fit,
//                   sqrt(1 + h) for one left out, h = its leverage. Its width is then the resolution of one crystal
//                   against the others, without the error of the line. It is the error that matters at Mu2e and it
//                   does not cancel what is common to L and R.
//                   Reported: half of the 16-84 % interval (main figure), sigma of a gaussian fit of the core,
//                   fraction > 1 ns.
//   L-R             ((tL - oL) - (tR - oR)) / sqrt2: resolution of one SiPM, blind to errors common to L and R.
//   phase bias      rms over 10 phase bins of the median track residual (and of L-R): a correct template gives a flat
//                   profile.
//   walk            median track residual in amplitude bins: its spread is an amplitude-dependent bias of the method.
//   efficiency      reference crystals in which the method has both SiPMs; chi2/ndf median on the track SiPMs (the
//                   method's own chi2).
//
// Outputs:
//   summaryFile     text tables: one line per method (figures of merit on the odd events, track residual on the even
//                   events, * = best of its family), the residual per crystal amplitude bin, the walk tables.
//   TFileService, one directory per method:
//     hChi2       chi2/ndf of the method on the track SiPMs (log bins, 100 per decade, 1e-2 to 1e3)
//     hChi2Amp    chi2/ndf versus the SiPM raw amplitude (20 ADC bins): millions of SiPMs, so the binning must not
//                 set the resolution of the figures
//     hRes, hLR   track residual and (tL - tR)/sqrt2, odd events
//     hResPhase, hLRPhase, hResAmp   the same versus the phase of L and the crystal amplitude
//   tree "crystals" (writeTree): one entry per odd-event crystal and method; branches method (index in methods),
//     crystal (id), res (track residual), lr ((tL - tR)/sqrt2), amp (mean raw amplitude of L and R), phase (of L),
//     trackType (1, 2 one disk, 3 two disks), nFit (crystals in the line fit).
//
// Steps:
//   analyze  for every reference track, its usable crystals: position along the track, raw amplitude and phase of
//            the two SiPMs, time of every method (all kept in memory until the end of the job).
//   endJob   for every method: offsets (and walk) on the even tracks, figures of merit on the even and odd tracks;
//            then the best method of each family and the summary tables.
//
// Main parameters: methods (name, family, caloRecoDigiCollection), nIterTime and outlierCut (offset least squares),
// walkCorrection and walkBins, minEntriesTime, minCrystals, maxDevOthers, depthInCrystal and diskZOffset (crystal
// positions, as given to the track finder).
//
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art_root_io/TFileService.h"
#include "cetlib_except/exception.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"
#include "fhiclcpp/types/Table.h"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/CalorimeterGeom/inc/Calorimeter.hh"
#include "Offline/CaloConditions/inc/CalDAQMap.hh"
#include "Offline/ProditionsService/inc/ProditionsHandle.hh"
#include "Offline/DataProducts/inc/CaloConst.hh"
#include "Offline/RecoDataProducts/inc/CaloRecoDigi.hh"
#include "Offline/RecoDataProducts/inc/CaloDigi.hh"
#include "Offline/MC2Reco/inc/MC2CaloCosmicTrack.hh"
#include "TF1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TTree.h"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>
#include <map>
#include <unordered_map>
#include <vector>

namespace mu2e {
  namespace {
    const double kC = 299.792458;   // mm/ns
    const float kNaN = std::numeric_limits<float>::quiet_NaN();
    // crystal amplitude bins of the residual table (mean of the L and R raw amplitudes) [ADC]
    const std::vector<double> kAmpBins = {150, 300, 500, 800, 1200, 1600, 4100};

    // A usable crystal of a reference track (both SiPMs), with the raw-waveform variables common to all methods.
    struct TrackCrystal {
      int track;          // index in tracks_
      int crystal;        // crystal id
      int sipm[2];        // SiPM id of side 0 (L) and side 1 (R)
      float s;            // position along the track [mm]
      float amp[2];       // raw amplitude of each SiPM [ADC]
      float phase[2];     // phase of each SiPM in the 5 ns sample [ns]
      bool saturated;     // one of the SiPMs saturated
    };

    // A reference track: its crystals are crystals_[first .. first + n - 1].
    struct RefTrack {
      int first;
      int n;
      int parity;   // event number % 2
      int type;     // 1, 2 one disk, 3 two disks
    };

    // Line t = intercept + slope * s through points (s, t) with weights use (0/1).
    // fixedSlope: slope = -1/c (two-disk tracks), otherwise the slope is free.
    // leverage(s) of a point at s: 1/n + (s - sMean)^2 / sxx (free slope) or 1/n.
    struct Line {
      double intercept = 0;
      double slope = 0;
      double sMean = 0;
      double sxx = 0;
      int n = 0;
      bool fixedSlope = false;

      double at(double s) const { return intercept + slope * s; }

      double leverage(double s) const {
        return n > 0 ? 1. / n + (fixedSlope || sxx <= 0 ? 0. : (s - sMean) * (s - sMean) / sxx) : 1.;
      }
    };

    // Fits the line; returns false if it cannot be fitted (too few points, or no lever arm for a free slope).
    bool fitLine(const std::vector<double>& s, const std::vector<double>& t, const std::vector<char>& use,
                 bool fixedSlope, Line& line) {
      line = Line();
      line.fixedSlope = fixedSlope;
      double sumS = 0;
      double sumT = 0;
      int n = 0;
      for (size_t i = 0; i < s.size(); ++i) {
        if (use[i]) {
          sumS += s[i];
          sumT += t[i] + (fixedSlope ? s[i] / kC : 0.);
          ++n;
        }
      }
      if (n < (fixedSlope ? 1 : 3)) {
        return false;
      }
      line.n = n;
      line.sMean = sumS / n;
      if (fixedSlope) {
        line.slope = -1. / kC;
        line.intercept = sumT / n;
        return true;
      }
      const double tMean = sumT / n;
      double sxx = 0;
      double sxt = 0;
      for (size_t i = 0; i < s.size(); ++i) {
        if (use[i]) {
          sxx += (s[i] - line.sMean) * (s[i] - line.sMean);
          sxt += (s[i] - line.sMean) * (t[i] - tMean);
        }
      }
      if (sxx < 100.) {   // crystals spread over less than ~10 mm: no slope
        return false;
      }
      line.sxx = sxx;
      line.slope = sxt / sxx;
      line.intercept = tMean - line.slope * line.sMean;
      return true;
    }

    // Value at fraction q of the sorted values (nearest rank).
    double quantile(std::vector<double> values, double q) {
      if (values.empty()) {
        return kNaN;
      }
      const size_t k = std::min(values.size() - 1, size_t(q * (values.size() - 1) + 0.5));
      std::nth_element(values.begin(), values.begin() + k, values.end());
      return values[k];
    }

    // Half of the 16-84 % interval (NaN with fewer than 20 values).
    double halfWidth1684(const std::vector<double>& values) {
      return values.size() < 20 ? kNaN : 0.5 * (quantile(values, 0.84) - quantile(values, 0.16));
    }

    double median(const std::vector<double>& values) { return quantile(values, 0.5); }

    // Sigma of a gaussian fit of the core (+-1.5 half-width of 16-84 % around the median).
    double coreSigma(const std::vector<double>& values) {
      if (values.size() < 50) {
        return kNaN;
      }
      const double center = median(values);
      const double halfWidth = halfWidth1684(values);
      if (!(halfWidth > 0)) {
        return kNaN;
      }
      TH1F histogram("hcore", "", 200, center - 4 * halfWidth, center + 4 * halfWidth);
      histogram.SetDirectory(nullptr);
      for (double x : values) {
        histogram.Fill(x);
      }
      TF1 gaussian("gcore", "gaus", center - 1.5 * halfWidth, center + 1.5 * halfWidth);
      gaussian.SetParameters(histogram.GetMaximum(), center, halfWidth);
      if (histogram.Fit(&gaussian, "QRN0") != 0) {
        return kNaN;
      }
      return std::abs(gaussian.GetParameter(2));
    }

    // rms over 10 phase bins (0-5 ns) of the median of y (bins with >= 30 entries; NaN with fewer than 3 bins).
    double phaseBias(const std::vector<double>& y, const std::vector<double>& phase) {
      const int nBins = 10;
      std::vector<std::vector<double>> yByBin(nBins);
      for (size_t i = 0; i < y.size(); ++i) {
        int bin = int(phase[i] / 5. * nBins);
        if (bin >= 0 && bin < nBins) {
          yByBin[bin].push_back(y[i]);
        }
      }
      std::vector<double> medians;
      for (auto& values : yByBin) {
        if (values.size() >= 30) {
          medians.push_back(median(values));
        }
      }
      if (medians.size() < 3) {
        return kNaN;
      }
      double mean = 0;
      for (double x : medians) {
        mean += x;
      }
      mean /= medians.size();
      double sumSquares = 0;
      for (double x : medians) {
        sumSquares += (x - mean) * (x - mean);
      }
      return std::sqrt(sumSquares / medians.size());
    }
  }

  class MC2RecoCompare : public art::EDAnalyzer {
    public:
      struct MethodConfig {
        using Name = fhicl::Name;
        using Comment = fhicl::Comment;
        fhicl::Atom<std::string> name{Name("name"), Comment("name of the method in the tables")};
        fhicl::Atom<std::string> family{Name("family"),
            Comment("methods of one family differ only in the fit window: the best on the even events is marked")};
        fhicl::Atom<art::InputTag> tag{Name("caloRecoDigiCollection"),
            Comment("CaloRecoDigi collection of the method")};
      };

      struct Config {
        using Name = fhicl::Name;
        using Comment = fhicl::Comment;

        // Inputs
        fhicl::Atom<art::InputTag> trackTag{Name("caloTrackCollection"), Comment("reference MC2CaloCosmicTracks"),
            "MC2CaloCosmicTrackFinder"};
        fhicl::Sequence<fhicl::Table<MethodConfig>> methods{Name("methods"), Comment("reconstructions to compare")};

        // Crystal positions (as given to the track finder)
        fhicl::Atom<double> depthInCrystal{Name("depthInCrystal"), Comment("crystal centre behind the front face [mm]"),
            100.};
        fhicl::Sequence<double> diskZOffset{Name("diskZOffset"), Comment("added to the Offline z of disk 0 and 1 [mm]"),
            std::vector<double>{0., 0.}};

        // Raw waveform
        fhicl::Atom<int> nPed{Name("nPed"), Comment("samples for the pedestal of the raw amplitude"), 4};
        fhicl::Atom<int> satADC{Name("satADC"), Comment("largest raw sample >= satADC: saturated"), 4090};

        // Offsets (and walk) on the even events
        fhicl::Atom<int> nIterTime{Name("nIterTime"), Comment("iterations of the offset least squares"), 25};
        fhicl::Atom<double> outlierCut{Name("outlierCut"), Comment("final max |residual| in the least squares [ns]"),
            5.};
        fhicl::Atom<bool> walkCorr{Name("walkCorrection"),
            Comment("fit a walk w(A) of the SiPM raw amplitude with the offsets and subtract it"), false};
        fhicl::Atom<int> walkBins{Name("walkBins"),
            Comment("logarithmic amplitude bins of the walk table, 150-2100 ADC"), 24};
        fhicl::Atom<int> minEntriesTime{Name("minEntriesTime"), Comment("min entries per SiPM for an offset"), 10};

        // Track residual
        fhicl::Atom<int> minCrystals{Name("minCrystals"),
            Comment("min usable crystals of a track (tested one included)"), 4};
        fhicl::Atom<double> maxDevOthers{Name("maxDevOthers"),
            Comment("other crystals farther than this from the track median are not used [ns]"), 5.};

        // Output
        fhicl::Atom<std::string> summaryFile{Name("summaryFile"), Comment("text tables"), "mc2_reco_compare.txt"};
        fhicl::Atom<bool> writeTree{Name("writeTree"), Comment("tree of the odd crystals, every method"), false};
      };

      using Parameters = art::EDAnalyzer::Table<Config>;
      explicit MC2RecoCompare(const Parameters& config);
      void beginJob() override;
      void analyze(const art::Event& event) override;
      void endJob() override;

    private:
      struct Method {
        std::string name;
        std::string family;
        art::InputTag tag;
      };

      // Figures of merit of one method on the tracks of one parity.
      struct FiguresOfMerit {
        double efficiency = kNaN;      // reference crystals in which the method has both SiPMs
        double chi2 = kNaN;            // median chi2/ndf (hChi2)
        double resHalf = kNaN;         // track residual: half 16-84 %
        double resCore = kNaN;         // track residual: gaussian core sigma
        double resTail = kNaN;         // track residual: fraction > 1 ns
        double lrHalf = kNaN;          // (tL - tR)/sqrt2: half 16-84 %
        double lrCore = kNaN;          // (tL - tR)/sqrt2: gaussian core sigma
        double phaseBiasRes = kNaN;    // phase bias of the track residual
        double phaseBiasLR = kNaN;     // phase bias of L-R
        double walk = kNaN;            // spread of the median residual over the amplitude bins
        long n = 0;                    // track residuals
        std::vector<double> medianByAmpBin;   // median residual per bin of kAmpBins
        std::vector<double> widthByAmpBin;    // half 16-84 % of the residual per bin of kAmpBins
      };

      // One SiPM time of an even track, an entry of the offset least squares.
      struct SiPMTime {
        int track;
        int sipm;
        double s;           // position along the track
        double t;           // time of the method
        double amplitude;   // raw amplitude
        bool use;           // not an outlier
      };

      // Per method: CaloDigi key -> CaloRecoDigi of the method (largest energy).
      using RecoDigiIndex = std::vector<std::unordered_map<size_t, const CaloRecoDigi*>>;

      // analyze steps
      RecoDigiIndex indexRecoDigisByDigi(const art::Event& event) const;
      void storeTrack(const MC2CaloCosmicTrack& track, int parity, const RecoDigiIndex& recoDigis,
                      const Calorimeter& calorimeter);
      double positionAlongTrack(const MC2CaloCosmicTrack& track, const Crystal& crystal,
                                const Calorimeter& calorimeter) const;

      // endJob steps
      std::vector<SiPMTime> evenTrackSiPMTimes(int method) const;
      void fitTrackLines(const std::vector<SiPMTime>& entries, const std::vector<std::vector<size_t>>& entriesByTrack,
                         const std::vector<double>& offsets, std::vector<Line>& lines,
                         std::vector<char>& lineOk) const;
      std::vector<int> fitOffsets(int method, std::vector<double>& offsets);
      void updateWalkTable(const std::vector<SiPMTime>& entries, const std::vector<Line>& lines,
                           const std::vector<char>& lineOk, const std::vector<double>& offsets);
      void figuresOfMerit(int method, const std::vector<double>& offsets, const std::vector<char>& calibrated,
                          int parity, FiguresOfMerit& fom, bool fill);
      std::map<std::string, int> bestOfEachFamily(const std::vector<FiguresOfMerit>& even) const;
      void writeSummary(const std::vector<FiguresOfMerit>& even, const std::vector<FiguresOfMerit>& odd,
                        const std::vector<int>& nCalibrated, const std::vector<std::vector<double>>& walkTables,
                        const std::map<std::string, int>& best) const;

      // helpers
      float methodTime(int crystalIndex, int method, int side) const;
      double calibratedTime(int crystalIndex, int method, int side, const std::vector<double>& offsets) const;
      double walkAt(double amplitude) const;

      // Parameters
      art::InputTag trackTag_;
      std::vector<Method> methods_;
      double depthInCrystal_;
      std::vector<double> diskZOffset_;
      int nPed_;
      int satADC_;
      int nIterTime_;
      double outlierCut_;
      int minEntriesTime_;
      int minCrystals_;
      double maxDevOthers_;
      std::string summaryFile_;
      bool writeTree_;
      bool walkCorrection_;
      int nWalkBins_;

      ProditionsHandle<CalDAQMap> daqMapHandle_;

      // State: the usable crystals of all reference tracks, kept until endJob
      std::vector<TrackCrystal> crystals_;
      std::vector<RefTrack> tracks_;
      std::vector<float> times_;   // times_[(crystal * nMethods + method) * 2 + side]: time of the method, NaN if none
      long nEvents_ = 0;
      std::vector<double> walkTable_;   // walk table of the method being processed in endJob (empty: none)

      // Output: histograms per method
      std::vector<TH1F*> hChi2_;
      std::vector<TH2F*> hChi2Amp_;
      std::vector<TH1F*> hRes_;
      std::vector<TH1F*> hLR_;
      std::vector<TH2F*> hResPhase_;
      std::vector<TH2F*> hResAmp_;
      std::vector<TH2F*> hLRPhase_;

      // Output: tree of the odd crystals
      TTree* tree_ = nullptr;
      int treeMethod_ = 0;
      int treeCrystal_ = 0;
      int treeTrackType_ = 0;
      int treeNFit_ = 0;
      float treeRes_ = 0;
      float treeLR_ = 0;
      float treeAmp_ = 0;
      float treePhase_ = 0;
  };

  MC2RecoCompare::MC2RecoCompare(const Parameters& config) :
    art::EDAnalyzer(config),
    trackTag_(config().trackTag()),
    depthInCrystal_(config().depthInCrystal()),
    diskZOffset_(config().diskZOffset()),
    nPed_(config().nPed()),
    satADC_(config().satADC()),
    nIterTime_(config().nIterTime()),
    outlierCut_(config().outlierCut()),
    minEntriesTime_(config().minEntriesTime()),
    minCrystals_(config().minCrystals()),
    maxDevOthers_(config().maxDevOthers()),
    summaryFile_(config().summaryFile()),
    writeTree_(config().writeTree()),
    walkCorrection_(config().walkCorr()),
    nWalkBins_(config().walkBins())
  {
    if (diskZOffset_.size() != 2) {
      throw cet::exception("MC2RECOCOMPARE") << "diskZOffset needs two values\n";
    }
    for (const auto& methodConfig : config().methods()) {
      methods_.push_back({methodConfig.name(), methodConfig.family(), methodConfig.tag()});
    }
    if (methods_.empty()) {
      throw cet::exception("MC2RECOCOMPARE") << "no methods\n";
    }
    consumes<MC2CaloCosmicTrackCollection>(trackTag_);
    for (const auto& method : methods_) {
      consumes<CaloRecoDigiCollection>(method.tag);
    }
  }

  void MC2RecoCompare::beginJob() {
    art::ServiceHandle<art::TFileService> tfs;
    // chi2/ndf log bins: 100 per decade from 1e-2 to 1e3
    std::vector<double> chi2Edges(501);
    for (int i = 0; i <= 500; ++i) {
      chi2Edges[i] = std::pow(10., -2. + i / 100.);
    }
    for (const auto& method : methods_) {
      art::TFileDirectory dir = tfs->mkdir(method.name);
      hChi2_.push_back(dir.make<TH1F>("hChi2", (method.name + ": track SiPMs;#chi^{2}/ndf;SiPMs").c_str(),
                                      500, chi2Edges.data()));
      hChi2Amp_.push_back(dir.make<TH2F>("hChi2Amp",
          (method.name + ": track SiPMs;SiPM raw amplitude [ADC];#chi^{2}/ndf").c_str(),
          205, 0, 4100, 500, chi2Edges.data()));
      hRes_.push_back(dir.make<TH1F>("hRes",
          (method.name + ": odd events;t_{crystal} - t_{other crystals} [ns];crystals").c_str(), 800, -4, 4));
      hLR_.push_back(dir.make<TH1F>("hLR",
          (method.name + ": odd events;(t_{L} - t_{R})/#sqrt{2} [ns];crystals").c_str(), 800, -4, 4));
      hResPhase_.push_back(dir.make<TH2F>("hResPhase",
          (method.name + ": odd events;phase of L [ns];track residual [ns]").c_str(), 50, 0, 5, 200, -2, 2));
      hLRPhase_.push_back(dir.make<TH2F>("hLRPhase",
          (method.name + ": odd events;phase of L [ns];(t_{L} - t_{R})/#sqrt{2} [ns]").c_str(),
          50, 0, 5, 200, -2, 2));
      hResAmp_.push_back(dir.make<TH2F>("hResAmp",
          (method.name + ": odd events;crystal amplitude (mean of L and R, raw) [ADC];track residual [ns]").c_str(),
          82, 0, 4100, 200, -2, 2));
    }
    if (writeTree_) {
      tree_ = tfs->make<TTree>("crystals", "odd-event track crystals, one entry per method");
      tree_->Branch("method", &treeMethod_);
      tree_->Branch("crystal", &treeCrystal_);
      tree_->Branch("res", &treeRes_);
      tree_->Branch("lr", &treeLR_);
      tree_->Branch("amp", &treeAmp_);
      tree_->Branch("phase", &treePhase_);
      tree_->Branch("trackType", &treeTrackType_);
      tree_->Branch("nFit", &treeNFit_);
    }
  }

  // ------------------------------------------------------------------------------------------------------------------
  // analyze: store the usable crystals of every reference track with the times of every method
  // ------------------------------------------------------------------------------------------------------------------

  void MC2RecoCompare::analyze(const art::Event& event) {
    ++nEvents_;
    const Calorimeter& calorimeter = *GeomHandle<Calorimeter>();
    const RecoDigiIndex recoDigis = indexRecoDigisByDigi(event);
    const int parity = event.id().event() % 2;
    for (const auto& track : *event.getValidHandle<MC2CaloCosmicTrackCollection>(trackTag_)) {
      storeTrack(track, parity, recoDigis, calorimeter);
    }
  }

  // Per method: CaloDigi key -> CaloRecoDigi made from it (the one with the largest energy).
  MC2RecoCompare::RecoDigiIndex MC2RecoCompare::indexRecoDigisByDigi(const art::Event& event) const {
    const int nMethods = methods_.size();
    RecoDigiIndex index(nMethods);
    for (int m = 0; m < nMethods; ++m) {
      for (const auto& recoDigi : *event.getValidHandle<CaloRecoDigiCollection>(methods_[m].tag)) {
        auto& stored = index[m][recoDigi.caloDigiPtr().key()];
        if (!stored || recoDigi.energyDep() > stored->energyDep()) {
          stored = &recoDigi;
        }
      }
    }
    return index;
  }

  // Position of the crystal centre along the track, pointing up [mm]: along the 3D direction for two-disk tracks,
  // along the direction in the disk plane for one-disk tracks.
  double MC2RecoCompare::positionAlongTrack(const MC2CaloCosmicTrack& track, const Crystal& crystal,
                                            const Calorimeter& calorimeter) const {
    const double inPlane = std::hypot(track.ux, track.uy);
    const int disk = crystal.diskID();
    const auto position = calorimeter.mu2eToTracker(crystal.position());
    const double z = position.z() + depthInCrystal_ + diskZOffset_[disk];
    double s;
    if (track.type == 3) {
      s = (position.x() - track.cx) * track.ux + (position.y() - track.cy) * track.uy + (z - track.cz) * track.uz;
    } else {
      s = ((position.x() - track.cx) * track.ux + (position.y() - track.cy) * track.uy) / (inPlane > 0 ? inPlane : 1.);
    }
    return s;
  }

  // The crystals of one track with both SiPMs (one each side, raw waveform longer than the pedestal): raw amplitude
  // and phase of each SiPM, time of every method. The chi2/ndf histograms are filled here for every SiPM read.
  void MC2RecoCompare::storeTrack(const MC2CaloCosmicTrack& track, int parity, const RecoDigiIndex& recoDigis,
                                  const Calorimeter& calorimeter) {
    const int nMethods = methods_.size();
    RefTrack stored{int(crystals_.size()), 0, parity, track.type};
    for (const auto& hit : track.hits) {
      if (hit->recoCaloDigis().size() != 2) {
        continue;
      }
      const double s = positionAlongTrack(track, calorimeter.crystal(hit->crystalID()), calorimeter);
      TrackCrystal crystal{int(tracks_.size()), hit->crystalID(), {-1, -1}, float(s), {0, 0}, {0, 0}, false};
      std::vector<float> times(2 * nMethods, kNaN);
      bool sidesOk = true;
      for (const auto& refRecoDigi : hit->recoCaloDigis()) {
        const auto& digiPtr = refRecoDigi->caloDigiPtr();
        const int side = refRecoDigi->SiPMID() % 2;
        if (crystal.sipm[side] >= 0) {
          sidesOk = false;
          break;
        }
        crystal.sipm[side] = refRecoDigi->SiPMID();

        // raw amplitude, saturation and phase from the waveform
        const auto& waveform = digiPtr->waveform();
        if ((int)waveform.size() <= nPed_ + 1) {
          sidesOk = false;
          break;
        }
        double pedestal = 0;
        for (int q = 0; q < nPed_; ++q) {
          pedestal += waveform[q];
        }
        pedestal /= nPed_;
        const int iMax = std::max_element(waveform.begin(), waveform.end()) - waveform.begin();
        const double amplitude = waveform[iMax] - pedestal;
        if (waveform[iMax] >= satADC_) {
          crystal.saturated = true;
        }
        crystal.amp[side] = amplitude;
        // 50 % crossing of the rising edge: first sample below half the amplitude going back from the peak
        double phase = kNaN;
        for (int q = iMax; q > 0; --q) {
          if (waveform[q - 1] - pedestal < 0.5 * amplitude) {
            const double y0 = waveform[q - 1] - pedestal;
            const double y1 = waveform[q] - pedestal;
            const double t50 = 5. * (q - 1 + (0.5 * amplitude - y0) / (y1 - y0));
            phase = std::fmod(t50, 5.);
            break;
          }
        }
        crystal.phase[side] = phase;

        // time of every method from the same CaloDigi, and its chi2/ndf
        for (int m = 0; m < nMethods; ++m) {
          const auto found = recoDigis[m].find(digiPtr.key());
          if (found == recoDigis[m].end()) {
            continue;
          }
          const CaloRecoDigi* recoDigi = found->second;
          times[2 * m + side] = recoDigi->time();
          if (recoDigi->ndf() > 0) {
            const double chi2ndf = recoDigi->chi2() / recoDigi->ndf();
            hChi2_[m]->Fill(chi2ndf);
            hChi2Amp_[m]->Fill(amplitude, chi2ndf);
          }
        }
      }
      if (!sidesOk || crystal.sipm[0] < 0 || crystal.sipm[1] < 0) {
        continue;
      }
      crystals_.push_back(crystal);
      times_.insert(times_.end(), times.begin(), times.end());
      ++stored.n;
    }
    if (stored.n > 0) {
      tracks_.push_back(stored);
    }
  }

  // ------------------------------------------------------------------------------------------------------------------
  // Helpers on the stored crystals
  // ------------------------------------------------------------------------------------------------------------------

  // Time of a method for one SiPM of a stored crystal (NaN if the method has none).
  float MC2RecoCompare::methodTime(int crystalIndex, int method, int side) const {
    const int nMethods = methods_.size();
    return times_[(size_t(crystalIndex) * nMethods + method) * 2 + side];
  }

  // Time of a method for one SiPM minus its offset and the walk at its raw amplitude.
  double MC2RecoCompare::calibratedTime(int crystalIndex, int method, int side,
                                        const std::vector<double>& offsets) const {
    const TrackCrystal& crystal = crystals_[crystalIndex];
    return methodTime(crystalIndex, method, side) - offsets[crystal.sipm[side]] - walkAt(crystal.amp[side]);
  }

  // Walk at the SiPM raw amplitude, from walkTable_ (0 if none): linear interpolation between the bin centres,
  // constant outside.
  double MC2RecoCompare::walkAt(double amplitude) const {
    if (walkTable_.empty() || !(amplitude > 0)) {
      return 0.;
    }
    // in units of bins, 0 = centre of the first
    const double x = std::log(amplitude / 150.) / std::log(2100. / 150.) * nWalkBins_ - 0.5;
    if (x <= 0) {
      return walkTable_.front();
    }
    if (x >= nWalkBins_ - 1) {
      return walkTable_.back();
    }
    const int i = int(x);
    return walkTable_[i] + (x - i) * (walkTable_[i + 1] - walkTable_[i]);
  }

  // ------------------------------------------------------------------------------------------------------------------
  // endJob: offsets and figures of merit of every method, summary tables
  // ------------------------------------------------------------------------------------------------------------------

  void MC2RecoCompare::endJob() {
    const int nMethods = methods_.size();
    const int nSiPMs = CaloConst::_nChannel;
    std::vector<FiguresOfMerit> even(nMethods);
    std::vector<FiguresOfMerit> odd(nMethods);
    std::vector<int> nCalibrated(nMethods, 0);
    std::vector<std::vector<double>> walkTables(nMethods);
    for (int m = 0; m < nMethods; ++m) {
      // offsets (and walk table) of the method on the even tracks
      std::vector<double> offsets(nSiPMs, 0.);
      const std::vector<int> nEntries = fitOffsets(m, offsets);
      std::vector<char> calibrated(nSiPMs, 0);
      for (int j = 0; j < nSiPMs; ++j) {
        if (nEntries[j] >= minEntriesTime_) {
          calibrated[j] = 1;
          ++nCalibrated[m];
        }
      }
      // figures of merit: even tracks (choice of the window within a family), odd tracks (results, histograms)
      figuresOfMerit(m, offsets, calibrated, 0, even[m], false);
      figuresOfMerit(m, offsets, calibrated, 1, odd[m], true);
      walkTables[m] = walkTable_;
    }
    writeSummary(even, odd, nCalibrated, walkTables, bestOfEachFamily(even));
  }

  // SiPM times of a method on the even tracks (crystals not saturated, the method has the time).
  std::vector<MC2RecoCompare::SiPMTime> MC2RecoCompare::evenTrackSiPMTimes(int method) const {
    std::vector<SiPMTime> entries;
    for (size_t k = 0; k < tracks_.size(); ++k) {
      const RefTrack& track = tracks_[k];
      if (track.parity != 0) {
        continue;
      }
      for (int i = track.first; i < track.first + track.n; ++i) {
        const TrackCrystal& crystal = crystals_[i];
        if (crystal.saturated) {
          continue;
        }
        for (int side = 0; side < 2; ++side) {
          const float t = methodTime(i, method, side);
          if (!std::isnan(t)) {
            entries.push_back({int(k), crystal.sipm[side], crystal.s, t, crystal.amp[side], true});
          }
        }
      }
    }
    return entries;
  }

  // Line of every track from its SiPM times minus the current offsets and walk.
  void MC2RecoCompare::fitTrackLines(const std::vector<SiPMTime>& entries,
                                     const std::vector<std::vector<size_t>>& entriesByTrack,
                                     const std::vector<double>& offsets, std::vector<Line>& lines,
                                     std::vector<char>& lineOk) const {
    std::vector<double> s;
    std::vector<double> t;
    std::vector<char> use;
    for (size_t k = 0; k < tracks_.size(); ++k) {
      if (entriesByTrack[k].empty()) {
        lineOk[k] = 0;
        continue;
      }
      s.clear();
      t.clear();
      use.clear();
      for (size_t i : entriesByTrack[k]) {
        const SiPMTime& x = entries[i];
        s.push_back(x.s);
        t.push_back(x.t - offsets[x.sipm] - walkAt(x.amplitude));
        use.push_back(x.use);
      }
      lineOk[k] = fitLine(s, t, use, tracks_[k].type == 3, lines[k]);
    }
  }

  // Per-SiPM offsets of a method on the even tracks: alternating least squares (line of every track, then the offsets
  // as the mean residual of each SiPM, mean offset 0 over the SiPMs with >= minEntriesTime entries), with a
  // progressive outlier cut and, with walkCorrection, the walk table (walkTable_). Returns the entries per SiPM of the
  // last iteration.
  std::vector<int> MC2RecoCompare::fitOffsets(int method, std::vector<double>& offsets) {
    const int nSiPMs = CaloConst::_nChannel;
    std::vector<SiPMTime> entries = evenTrackSiPMTimes(method);
    std::vector<std::vector<size_t>> entriesByTrack(tracks_.size());
    for (size_t i = 0; i < entries.size(); ++i) {
      entriesByTrack[entries[i].track].push_back(i);
    }
    std::vector<double> sumResiduals(nSiPMs);
    std::vector<int> nEntries(nSiPMs);
    std::vector<Line> lines(tracks_.size());
    std::vector<char> lineOk(tracks_.size(), 0);
    walkTable_.clear();

    fitTrackLines(entries, entriesByTrack, offsets, lines, lineOk);
    for (int iter = 0; iter < nIterTime_; ++iter) {
      // progressive outlier cut every third iteration: 50, 17, 5.6 ns, ..., then outlierCut
      if (iter >= 3 && iter % 3 == 0) {
        const double cut = std::max(outlierCut_, 50. / std::pow(3., (iter - 3) / 3));
        for (auto& x : entries) {
          x.use = lineOk[x.track]
                  && std::abs(x.t - lines[x.track].at(x.s) - offsets[x.sipm] - walkAt(x.amplitude)) < cut;
        }
      }

      // offsets: mean residual from the line of each SiPM, then mean offset 0
      std::fill(sumResiduals.begin(), sumResiduals.end(), 0.);
      std::fill(nEntries.begin(), nEntries.end(), 0);
      for (auto& x : entries) {
        if (x.use && lineOk[x.track] && lines[x.track].n > 1) {
          sumResiduals[x.sipm] += x.t - lines[x.track].at(x.s) - walkAt(x.amplitude);
          ++nEntries[x.sipm];
        }
      }
      double meanOffset = 0;
      int nOffsets = 0;
      for (int j = 0; j < nSiPMs; ++j) {
        if (nEntries[j] >= minEntriesTime_) {
          offsets[j] = sumResiduals[j] / nEntries[j];
          meanOffset += offsets[j];
          ++nOffsets;
        } else {
          offsets[j] = 0;
        }
      }
      if (nOffsets) {
        meanOffset /= nOffsets;
        for (int j = 0; j < nSiPMs; ++j) {
          if (nEntries[j] >= minEntriesTime_) {
            offsets[j] -= meanOffset;
          }
        }
      }

      if (walkCorrection_) {
        updateWalkTable(entries, lines, lineOk, offsets);
      }
      fitTrackLines(entries, entriesByTrack, offsets, lines, lineOk);
    }
    return nEntries;
  }

  // Walk table: mean of t - line - offset per amplitude bin (bins with >= 50 entries, otherwise the previous value),
  // then mean 0 over the entries.
  void MC2RecoCompare::updateWalkTable(const std::vector<SiPMTime>& entries, const std::vector<Line>& lines,
                                       const std::vector<char>& lineOk, const std::vector<double>& offsets) {
    std::vector<double> sumByBin(nWalkBins_, 0.);
    std::vector<long> nByBin(nWalkBins_, 0);
    double total = 0;
    long nTotal = 0;
    for (auto& x : entries) {
      if (!(x.use && lineOk[x.track] && lines[x.track].n > 1) || !(x.amplitude > 0)) {
        continue;
      }
      const int bin = std::clamp(int(std::log(x.amplitude / 150.) / std::log(2100. / 150.) * nWalkBins_),
                                 0, nWalkBins_ - 1);
      sumByBin[bin] += x.t - lines[x.track].at(x.s) - offsets[x.sipm];
      ++nByBin[bin];
    }
    std::vector<double> walk = walkTable_.empty() ? std::vector<double>(nWalkBins_, 0.) : walkTable_;
    for (int bin = 0; bin < nWalkBins_; ++bin) {
      if (nByBin[bin] >= 50) {
        walk[bin] = sumByBin[bin] / nByBin[bin];
        total += walk[bin] * nByBin[bin];
        nTotal += nByBin[bin];
      }
    }
    if (nTotal) {
      for (auto& value : walk) {
        value -= total / nTotal;
      }
    }
    walkTable_ = walk;
  }

  // Figures of merit of a method on the tracks of one parity, with the given offsets (calibrated[j]: SiPM j has an
  // offset). fill: fill the histograms and the tree (odd events).
  void MC2RecoCompare::figuresOfMerit(int method, const std::vector<double>& offsets,
                                      const std::vector<char>& calibrated, int parity, FiguresOfMerit& fom,
                                      bool fill) {
    const int m = method;
    const double sqrt2 = std::sqrt(2.);
    std::vector<double> residuals;
    std::vector<double> lrValues;
    std::vector<double> lrPhases;
    std::vector<double> residualPhases;
    std::vector<std::vector<double>> residualsByAmpBin(kAmpBins.size() - 1);
    long nReference = 0;
    long nBothSiPMs = 0;
    std::vector<double> crystalTimes;
    std::vector<double> crystalS;
    std::vector<int> crystalIndices;
    for (const auto& track : tracks_) {
      if (track.parity != parity) {
        continue;
      }

      // crystal times of the track (offsets and walk subtracted) and L-R
      crystalTimes.clear();
      crystalS.clear();
      crystalIndices.clear();
      for (int i = track.first; i < track.first + track.n; ++i) {
        const TrackCrystal& crystal = crystals_[i];
        if (crystal.saturated) {
          continue;
        }
        ++nReference;
        const float tL = methodTime(i, m, 0);
        const float tR = methodTime(i, m, 1);
        if (std::isnan(tL) || std::isnan(tR)) {
          continue;
        }
        ++nBothSiPMs;
        if (!calibrated[crystal.sipm[0]] || !calibrated[crystal.sipm[1]]) {
          continue;
        }
        const double left = calibratedTime(i, m, 0, offsets);
        const double right = calibratedTime(i, m, 1, offsets);
        const double lr = (left - right) / sqrt2;
        lrValues.push_back(lr);
        lrPhases.push_back(crystal.phase[0]);
        if (fill) {
          hLR_[m]->Fill(lr);
          hLRPhase_[m]->Fill(crystal.phase[0], lr);
        }
        crystalTimes.push_back(0.5 * (left + right));
        crystalS.push_back(crystal.s);
        crystalIndices.push_back(i);
      }
      if ((int)crystalTimes.size() < minCrystals_) {
        continue;
      }

      // line through the crystals, outliers (> maxDevOthers) left out iteratively
      std::vector<char> use(crystalTimes.size(), 1);
      Line line;
      const bool fixedSlope = track.type == 3;
      bool ok = fitLine(crystalS, crystalTimes, use, fixedSlope, line);
      for (int iter = 0; ok && iter < 3; ++iter) {
        bool changed = false;
        for (size_t q = 0; q < crystalTimes.size(); ++q) {
          const char inFit = std::abs(crystalTimes[q] - line.at(crystalS[q])) < maxDevOthers_;
          if (inFit != use[q]) {
            use[q] = inFit;
            changed = true;
          }
        }
        if (!changed) {
          break;
        }
        ok = fitLine(crystalS, crystalTimes, use, fixedSlope, line);
      }
      if (!ok || line.n < minCrystals_) {
        continue;
      }

      // studentised residual of every crystal
      for (size_t k = 0; k < crystalTimes.size(); ++k) {
        const double h = line.leverage(crystalS[k]);
        const double denominator = use[k] ? 1. - h : 1. + h;
        if (denominator <= 0.05) {
          continue;
        }
        const double residual = (crystalTimes[k] - line.at(crystalS[k])) / std::sqrt(denominator);
        const int i = crystalIndices[k];
        const TrackCrystal& crystal = crystals_[i];
        const double amplitude = 0.5 * (crystal.amp[0] + crystal.amp[1]);
        residuals.push_back(residual);
        residualPhases.push_back(crystal.phase[0]);
        for (size_t b = 0; b + 1 < kAmpBins.size(); ++b) {
          if (amplitude >= kAmpBins[b] && amplitude < kAmpBins[b + 1]) {
            residualsByAmpBin[b].push_back(residual);
          }
        }
        if (fill) {
          hRes_[m]->Fill(residual);
          hResPhase_[m]->Fill(crystal.phase[0], residual);
          hResAmp_[m]->Fill(amplitude, residual);
          if (tree_) {
            treeMethod_ = m;
            treeCrystal_ = crystal.crystal;
            treeRes_ = residual;
            treeAmp_ = amplitude;
            treePhase_ = crystal.phase[0];
            treeTrackType_ = track.type;
            treeNFit_ = line.n;
            treeLR_ = (calibratedTime(i, m, 0, offsets) - calibratedTime(i, m, 1, offsets)) / sqrt2;
            tree_->Fill();
          }
        }
      }
    }

    // summary figures
    fom.n = residuals.size();
    fom.efficiency = nReference ? double(nBothSiPMs) / nReference : kNaN;
    fom.resHalf = halfWidth1684(residuals);
    fom.resCore = coreSigma(residuals);
    long nTail = 0;
    for (double residual : residuals) {
      if (std::abs(residual) > 1.) {
        ++nTail;
      }
    }
    fom.resTail = residuals.empty() ? kNaN : double(nTail) / residuals.size();
    fom.lrHalf = halfWidth1684(lrValues);
    fom.lrCore = coreSigma(lrValues);
    fom.phaseBiasLR = phaseBias(lrValues, lrPhases);

    // walk: spread of the median residual over the amplitude bins with >= 50 entries
    fom.medianByAmpBin.clear();
    fom.widthByAmpBin.clear();
    double walkMin = 1e9;
    double walkMax = -1e9;
    for (auto& values : residualsByAmpBin) {
      const bool enough = values.size() >= 50;
      const double binMedian = enough ? median(values) : kNaN;
      fom.medianByAmpBin.push_back(binMedian);
      fom.widthByAmpBin.push_back(enough ? halfWidth1684(values) : kNaN);
      if (enough) {
        walkMin = std::min(walkMin, binMedian);
        walkMax = std::max(walkMax, binMedian);
      }
    }
    fom.walk = walkMax > walkMin ? walkMax - walkMin : kNaN;
    fom.phaseBiasRes = phaseBias(residuals, residualPhases);

    // median chi2/ndf of the method (histogram of the whole job, both parities)
    TH1F* hChi2 = hChi2_[m];
    if (hChi2->GetEntries() > 0) {
      double q = 0.5;
      double x;
      hChi2->GetQuantiles(1, &x, &q);
      fom.chi2 = x;
    }
  }

  // Per family: the method with the smallest track residual (resHalf) on the even events.
  std::map<std::string, int> MC2RecoCompare::bestOfEachFamily(const std::vector<FiguresOfMerit>& even) const {
    std::map<std::string, int> best;
    for (int m = 0; m < (int)methods_.size(); ++m) {
      const auto found = best.find(methods_[m].family);
      if (!std::isnan(even[m].resHalf) && (found == best.end() || even[m].resHalf < even[found->second].resHalf)) {
        best[methods_[m].family] = m;
      }
    }
    return best;
  }

  void MC2RecoCompare::writeSummary(const std::vector<FiguresOfMerit>& even, const std::vector<FiguresOfMerit>& odd,
                                    const std::vector<int>& nCalibrated,
                                    const std::vector<std::vector<double>>& walkTables,
                                    const std::map<std::string, int>& best) const {
    const int nMethods = methods_.size();
    FILE* out = fopen(summaryFile_.c_str(), "w");
    if (!out) {
      throw cet::exception("MC2RECOCOMPARE") << "cannot write " << summaryFile_ << "\n";
    }
    long nTracks[2] = {0, 0};
    long nCrystals[2] = {0, 0};
    for (const auto& track : tracks_) {
      ++nTracks[track.parity];
      nCrystals[track.parity] += track.n;
    }

    // header
    fprintf(out, "MC2RecoCompare: %ld events; reference tracks %ld even / %ld odd, "
                 "crystals with both SiPMs %ld / %ld\n",
            nEvents_, nTracks[0], nTracks[1], nCrystals[0], nCrystals[1]);
    fprintf(out, "offsets per method from the even tracks (global least squares), figures of merit on the odd "
                 "tracks; * = best of its family on the even tracks\n");
    fprintf(out, "res = crystal time - line of the track (T - s/c two disks, T + k s one disk), studentised [ns]: "
                 "half of 16-84 %% (main figure), gaussian core, fraction > 1 ns\n");
    fprintf(out, "LR = (tL - tR)/sqrt2 [ns]; phase = rms over the sampling phase of the median residual [ns]; "
                 "walk = spread of the median residual over the amplitude bins [ns]\n\n");

    // figures of merit, one line per method
    fprintf(out, "%-28s %-14s | %5s %6s %6s | %8s | %8s %8s %6s | %7s %7s | %7s %7s | %6s | %8s\n",
            "method", "family", "eff", "chi2", "nCalib", "res even", "res", "core", ">1ns", "LR", "LR core",
            "ph res", "ph LR", "walk", "crystals");
    for (int m = 0; m < nMethods; ++m) {
      const FiguresOfMerit& fom = odd[m];
      const std::string& family = methods_[m].family;
      const bool isBest = best.count(family) && best.at(family) == m;
      fprintf(out, "%-28s %-14s | %5.3f %6.2f %6d | %8.4f | %7.4f%s %8.4f %6.3f | %7.4f %7.4f | %7.4f %7.4f | "
                   "%6.3f | %8ld\n",
              methods_[m].name.c_str(), family.c_str(), fom.efficiency, fom.chi2, nCalibrated[m], even[m].resHalf,
              fom.resHalf, isBest ? "*" : " ", fom.resCore, fom.resTail, fom.lrHalf, fom.lrCore, fom.phaseBiasRes,
              fom.phaseBiasLR, fom.walk, fom.n);
    }

    // residual per crystal amplitude bin
    fprintf(out, "\nper amplitude bin (crystal amplitude = mean of L and R raw amplitudes [ADC]), odd tracks: "
                 "half 16-84 %% of the residual / median residual [ns]\n");
    fprintf(out, "%-28s", "method");
    for (size_t b = 0; b + 1 < kAmpBins.size(); ++b) {
      fprintf(out, " | %4.0f-%-4.0f      ", kAmpBins[b], kAmpBins[b + 1]);
    }
    fprintf(out, "\n");
    for (int m = 0; m < nMethods; ++m) {
      fprintf(out, "%-28s", methods_[m].name.c_str());
      for (size_t b = 0; b < odd[m].widthByAmpBin.size(); ++b) {
        fprintf(out, " | %6.3f / %+6.3f", odd[m].widthByAmpBin[b], odd[m].medianByAmpBin[b]);
      }
      fprintf(out, "\n");
    }

    // walk tables from the even tracks, subtracted from the SiPM times
    if (walkCorrection_) {
      fprintf(out, "\nwalk correction w(A) [ns] from the even tracks, %d logarithmic bins of the SiPM raw amplitude "
                   "from 150 to 2100 ADC (bin centres):\n%-28s",
              nWalkBins_, "method");
      for (int b = 0; b < nWalkBins_; ++b) {
        fprintf(out, " %6.0f", 150. * std::pow(2100. / 150., (b + 0.5) / nWalkBins_));
      }
      fprintf(out, "\n");
      for (int m = 0; m < nMethods; ++m) {
        fprintf(out, "%-28s", methods_[m].name.c_str());
        for (double value : walkTables[m]) {
          fprintf(out, " %+6.3f", value);
        }
        fprintf(out, "\n");
      }
    }
    fclose(out);
    printf("MC2RecoCompare: %zu methods, %ld + %ld tracks -> %s\n", methods_.size(), nTracks[0], nTracks[1],
           summaryFile_.c_str());
  }
}

DEFINE_ART_MODULE(mu2e::MC2RecoCompare)
