//
// MC2CaloRecoDigiMaker: CaloDigi -> CaloRecoDigi with the MC2 template fit (MC2TemplateFit) and per-SiPM time offsets.
//
// It replaces mu2e::CaloRecoDigiMaker (Offline/CaloReco) for the MC2 cosmic-ray runs and writes the same data
// product, so that CaloHitMaker and everything downstream run unchanged: set CaloHitMaker.caloRecoDigiCollection to
// this module label.
// Both reconstructions can run in the same job (two collections) and be compared on the same events.
//
// For every CaloDigi:
//
//   1. template fit of the waveform (MC2TemplateFit): amplitude A, time of the peak, chi2/ndf
//      (with perChannel, the template of that SiPM if there is one, otherwise the common template)
//   2. the digi is dropped if the fit failed, if A < minAmplitude or if chi2/ndf > maxChi2Cut
//   3. calibrated values:
//        time    = 5 ns * t0(CaloDigi) + fitted peak time - [proton bunch time - calo-tracker offset] - tcorr(SiPM)
//                  (the bracket only with UseProtonBunchTime, exactly as in CaloRecoDigiMaker)
//        eDep    = A * ADC2MeV(SiPM) of the CalCalib conditions, as in CaloRecoDigiMaker
//        timeErr = sqrt((timeErrA / A)^2 + timeErrC^2), from the MC2 single-SiPM resolution on cosmics
//
// Per-SiPM time offset tcorr, one of the two (not both, exception):
//   useCalCalibTimeOffset : the CalCalib time offset (CalTimeCalib table, e.g. from DbService.textFile),
//                           exactly as CaloRecoDigiMaker
//   timeOffsetFile        : text file with "SiPMID tcorr[ns]" per line (# for comments); 0 for the SiPMs not listed
//
// Per-channel templates (perChannel): histograms perChannelFormat (board, channel of the raw id, e.g. tmpl_b012_c07) of
// templateFile, as written by mc2-ana laser_template_stream.C; they are read at the first event (the DAQ map is a
// condition).
//
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/types/Atom.h"
#include "cetlib_except/exception.h"

#include "Offline/CaloConditions/inc/CalCalib.hh"
#include "Offline/CaloConditions/inc/CalDAQMap.hh"
#include "Offline/ConfigTools/inc/ConfigFileLookupPolicy.hh"
#include "Offline/DAQConditions/inc/EventTiming.hh"
#include "Offline/DataProducts/inc/CaloConst.hh"
#include "Offline/ProditionsService/inc/ProditionsHandle.hh"
#include "Offline/RecoDataProducts/inc/CaloDigi.hh"
#include "Offline/RecoDataProducts/inc/CaloRecoDigi.hh"
#include "Offline/RecoDataProducts/inc/ProtonBunchTime.hh"
#include "Offline/MC2Reco/inc/MC2TemplateFit.hh"
#include "TFile.h"
#include "TGraph.h"
#include "TH1.h"

#include <cmath>
#include <fstream>
#include <iostream>
#include <memory>
#include <sstream>
#include <unordered_map>

namespace mu2e {

  class MC2CaloRecoDigiMaker : public art::EDProducer {
    public:
      struct Config {
        using Name = fhicl::Name;
        using Comment = fhicl::Comment;

        // Input
        fhicl::Atom<art::InputTag> caloDigiCollection{Name("caloDigiCollection"), Comment("CaloDigi module label")};

        // Template and fit (see MC2TemplateFit.hh)
        fhicl::Atom<std::string> templateFile{Name("templateFile"),
            Comment("production-format template file (TGraph gtempl)")};
        fhicl::Atom<std::string> graphName{Name("graphName"), Comment("name of the template graph"), "gtempl"};
        fhicl::Atom<double> fitLo{Name("fitLo"), Comment("start of the fit window from the fitted peak (ns)"), -80.};
        fhicl::Atom<double> fitHi{Name("fitHi"), Comment("end of the fit window from the fitted peak (ns)"), -31.};
        fhicl::Atom<double> noise{Name("noise"), Comment("noise per sample (ADC), for chi2/ndf"), 3.6};
        fhicl::Atom<int> nPedSamples{Name("nPedSamples"), Comment("samples used for the pedestal"), 4};
        fhicl::Atom<double> maxBaseline{Name("maxBaseline"), Comment("limit of the residual baseline (ADC)"), 10.};
        fhicl::Atom<int> maxPre{Name("maxPre"),
            Comment("samples before the largest one allowed in the fit (<0: all)"), 16};
        fhicl::Atom<int> shapeMode{Name("shapeMode"),
            Comment("0 fixed template; amplitude model in templateFile: 1 shape tied to the fitted amplitude, "
                    "2 free shape, 3 constrained shape"), 0};
        fhicl::Atom<double> shapeSigma{Name("shapeSigma"),
            Comment("shapeMode 3: spread of the shape between hits (units of ln A)"), 0.5};
        fhicl::Atom<int> satADC{Name("satADC"), Comment("raw samples >= satADC are left out of the fit"), 4090};
        fhicl::Atom<bool> perChannel{Name("perChannel"),
            Comment("one template per SiPM from the histograms perChannelFormat of templateFile"), false};
        fhicl::Atom<std::string> perChannelFormat{Name("perChannelFormat"),
            Comment("printf format of the per-channel histogram names (board, channel)"), "tmpl_b%03d_c%02d"};

        // Selection of the reco digis
        fhicl::Atom<double> minAmplitude{Name("minAmplitude"), Comment("minimum fitted amplitude (ADC)"), 0.};
        fhicl::Atom<double> maxChi2Cut{Name("maxChi2Cut"), Comment("maximum chi2/ndf to keep the reco digi"), 1e9};

        // Time error and time offsets
        fhicl::Atom<double> timeErrA{Name("timeErrA"), Comment("time error term a/A (ns ADC)"), 110.};
        fhicl::Atom<double> timeErrC{Name("timeErrC"), Comment("time error constant term (ns)"), 0.18};
        fhicl::Atom<std::string> timeOffsetFile{Name("timeOffsetFile"),
            Comment("per-SiPM time offsets: SiPMID tcorr[ns]; empty = none"), ""};
        fhicl::Atom<bool> useCalCalibTime{Name("useCalCalibTimeOffset"),
            Comment("subtract the CalCalib time offset (CalTimeCalib table), as CaloRecoDigiMaker"), false};
        fhicl::Atom<art::InputTag> pbtTag{Name("ProtonBunchTimeTag"), Comment("ProtonBunchTime producer"),
            art::InputTag("EWMProducer")};
        fhicl::Atom<bool> usePBT{Name("UseProtonBunchTime"),
            Comment("subtract the proton bunch time as CaloRecoDigiMaker"), false};

        fhicl::Atom<int> diagLevel{Name("diagLevel"), Comment("diagnosis level"), 0};
      };

      explicit MC2CaloRecoDigiMaker(const art::EDProducer::Table<Config>& config);
      void produce(art::Event& event) override;

    private:
      static MC2TemplateFit::Config fitConfig(const Config& config);
      void readTimeOffsetFile(const std::string& fileName);
      void loadPerChannelTemplates(const art::Event& event);
      const MC2TemplateFit& fitterFor(int sipmId) const;
      double timeOffset(int sipmId, const CalCalib& calCalib) const;
      double protonBunchTimeOffset(const art::Event& event);

      // Input
      art::ProductToken<CaloDigiCollection> caloDigisToken_;
      art::ProductToken<ProtonBunchTime> pbtToken_;

      // Fit
      MC2TemplateFit::Config fitConfig_;
      MC2TemplateFit fitter_;                                       // common template
      bool perChannel_;
      std::string perChannelFormat_;
      std::vector<std::unique_ptr<MC2TemplateFit>> perChannelFitters_;   // by offline SiPM id; empty: common template
      bool perChannelLoaded_ = false;

      // Selection, errors, offsets
      double minAmplitude_;
      double maxChi2Cut_;
      double timeErrA_;
      double timeErrC_;
      bool useCalCalibTime_;
      bool usePBT_;
      std::unordered_map<int, double> timeOffsetFromFile_;   // SiPM id -> tcorr (ns)

      // Conditions
      ProditionsHandle<CalCalib> calCalibHandle_;
      ProditionsHandle<EventTiming> eventTimingHandle_;
      ProditionsHandle<CalDAQMap> daqMapHandle_;

      // Counters
      int diagLevel_;
      long nDigis_ = 0;
      long nFailedFits_ = 0;
      long nChi2Cut_ = 0;
  };

  MC2TemplateFit::Config MC2CaloRecoDigiMaker::fitConfig(const Config& config) {
    ConfigFileLookupPolicy resolve;
    MC2TemplateFit::Config fit;
    fit.templateFile = resolve(config.templateFile());
    fit.graphName = config.graphName();
    fit.fitLo = config.fitLo();
    fit.fitHi = config.fitHi();
    fit.noise = config.noise();
    fit.nPedSamples = config.nPedSamples();
    fit.maxBaseline = config.maxBaseline();
    fit.satADC = config.satADC();
    fit.maxPre = config.maxPre();
    fit.shapeMode = config.shapeMode();
    fit.shapeSigma = config.shapeSigma();
    return fit;
  }

  MC2CaloRecoDigiMaker::MC2CaloRecoDigiMaker(const art::EDProducer::Table<Config>& config) :
    EDProducer{config},
    caloDigisToken_(consumes<CaloDigiCollection>(config().caloDigiCollection())),
    fitConfig_(fitConfig(config())),
    fitter_(fitConfig_),
    perChannel_(config().perChannel()),
    perChannelFormat_(config().perChannelFormat()),
    minAmplitude_(config().minAmplitude()),
    maxChi2Cut_(config().maxChi2Cut()),
    timeErrA_(config().timeErrA()),
    timeErrC_(config().timeErrC()),
    useCalCalibTime_(config().useCalCalibTime()),
    usePBT_(config().usePBT()),
    diagLevel_(config().diagLevel())
  {
    produces<CaloRecoDigiCollection>();
    if (usePBT_) {
      pbtToken_ = consumes<ProtonBunchTime>(config().pbtTag());
    }

    const std::string offsetFile = config().timeOffsetFile();
    if (useCalCalibTime_ && !offsetFile.empty()) {
      throw cet::exception("MC2CALORECO")
          << "timeOffsetFile and useCalCalibTimeOffset together would subtract two offsets: use one\n";
    }
    if (!offsetFile.empty()) {
      readTimeOffsetFile(offsetFile);
    }
  }

  // Per-SiPM time offsets from a text file: "SiPMID tcorr" per line, # for comments.
  void MC2CaloRecoDigiMaker::readTimeOffsetFile(const std::string& fileName) {
    ConfigFileLookupPolicy resolve;
    const std::string path = resolve(fileName);
    std::ifstream in(path);
    if (!in) {
      throw cet::exception("MC2CALORECO") << "cannot open timeOffsetFile " << path << "\n";
    }
    std::string line;
    while (std::getline(in, line)) {
      if (line.empty() || line[0] == '#') {
        continue;
      }
      std::istringstream fields(line);
      int sipmId;
      double offset;
      if (fields >> sipmId >> offset) {
        timeOffsetFromFile_[sipmId] = offset;
      }
    }
    std::cout << "[MC2CaloRecoDigiMaker] " << timeOffsetFromFile_.size() << " SiPM time offsets from " << path
              << std::endl;
  }

  // Per-channel templates: for each offline SiPM, the histogram of its (board, channel) -> graph of the bin centres ->
  // its own fitter. SiPMs without a histogram (or with fewer than 10 bins) keep the common template.
  void MC2CaloRecoDigiMaker::loadPerChannelTemplates(const art::Event& event) {
    perChannelLoaded_ = true;
    if (!perChannel_) {
      return;
    }
    const CalDAQMap& daqMap = daqMapHandle_.get(event.id());
    std::unique_ptr<TFile> file(TFile::Open(fitConfig_.templateFile.c_str(), "READ"));
    if (!file || file->IsZombie()) {
      throw cet::exception("MC2CALORECO") << "cannot open " << fitConfig_.templateFile << "\n";
    }

    perChannelFitters_.resize(CaloConst::_nChannel);
    int nLoaded = 0;
    for (int sipmId = 0; sipmId < (int)CaloConst::_nChannel; ++sipmId) {
      const auto rawId = daqMap.rawId(CaloSiPMId(sipmId));
      char histogramName[64];
      snprintf(histogramName, sizeof(histogramName), perChannelFormat_.c_str(),
               int(rawId.dirac()), int(rawId.ROCchannel()));
      auto* histogram = dynamic_cast<TH1*>(file->Get(histogramName));
      if (!histogram || histogram->GetEntries() <= 0) {
        continue;
      }
      TGraph graph;
      for (int bin = 1; bin <= histogram->GetNbinsX(); ++bin) {
        graph.SetPoint(graph.GetN(), histogram->GetBinCenter(bin), histogram->GetBinContent(bin));
      }
      if (graph.GetN() < 10) {
        continue;
      }
      perChannelFitters_[sipmId] = std::make_unique<MC2TemplateFit>(fitConfig_, graph);
      ++nLoaded;
    }
    std::cout << "[MC2CaloRecoDigiMaker] per-channel templates for " << nLoaded << " SiPMs from "
              << fitConfig_.templateFile << " (common template for the others)" << std::endl;
  }

  const MC2TemplateFit& MC2CaloRecoDigiMaker::fitterFor(int sipmId) const {
    const bool hasOwnTemplate = sipmId >= 0 && sipmId < (int)perChannelFitters_.size() && perChannelFitters_[sipmId];
    return hasOwnTemplate ? *perChannelFitters_[sipmId] : fitter_;
  }

  double MC2CaloRecoDigiMaker::timeOffset(int sipmId, const CalCalib& calCalib) const {
    if (useCalCalibTime_) {
      return calCalib.timeoffset(sipmId);
    }
    const auto found = timeOffsetFromFile_.find(sipmId);
    return found == timeOffsetFromFile_.end() ? 0. : found->second;
  }

  // Proton bunch time minus the calo-tracker offset, as CaloRecoDigiMaker (0 without UseProtonBunchTime).
  double MC2CaloRecoDigiMaker::protonBunchTimeOffset(const art::Event& event) {
    if (!usePBT_) {
      return 0;
    }
    double offset = event.getValidHandle(pbtToken_)->pbtime_;
    offset -= eventTimingHandle_.get(event.id()).caloTrackerTimeOffset();
    return offset;
  }

  void MC2CaloRecoDigiMaker::produce(art::Event& event) {
    if (!perChannelLoaded_) {
      loadPerChannelTemplates(event);
    }
    const auto& digisHandle = event.getValidHandle(caloDigisToken_);
    const auto& digis = *digisHandle;
    auto recoDigis = std::make_unique<CaloRecoDigiCollection>();

    const double pbtOffset = protonBunchTimeOffset(event);
    const CalCalib& calCalib = calCalibHandle_.get(event.id());

    for (size_t i = 0; i < digis.size(); ++i) {
      const CaloDigi& digi = digis[i];
      ++nDigis_;
      const int sipmId = digi.SiPMID();

      // 1. fit
      const auto fit = fitterFor(sipmId).fit(digi.waveform());

      // 2. selection
      if (!fit.ok || fit.amplitude < minAmplitude_) {
        ++nFailedFits_;
        continue;
      }
      if (fit.chi2ndf > maxChi2Cut_) {
        ++nChi2Cut_;
        continue;
      }

      // 3. calibrated time, energy and errors
      const double time = digi.t0() * CaloConst::_digitizationPeriod + fit.time - pbtOffset
                        - timeOffset(sipmId, calCalib);
      const double adc2MeV = calCalib.ADC2MeV(sipmId);
      const double eDep = fit.amplitude * adc2MeV;
      const double eDepErr = std::sqrt(std::max(fit.amplitude, 1.)) * adc2MeV;
      const double timeErr = std::hypot(timeErrA_ / std::max(fit.amplitude, 1.), timeErrC_);

      art::Ptr<CaloDigi> digiPtr(digisHandle, i);
      recoDigis->emplace_back(CaloRecoDigi(digiPtr, eDep, eDepErr, time, timeErr, fit.chi2ndf * fit.ndf, fit.ndf,
                                           false));
      if (diagLevel_ > 2) {
        std::cout << "[MC2CaloRecoDigiMaker] SiPM " << sipmId << " A=" << fit.amplitude << " t=" << time
                  << " chi2/ndf=" << fit.chi2ndf << std::endl;
      }
    }

    if (diagLevel_ > 0 && event.id().event() % 10000 == 0) {
      std::cout << "[MC2CaloRecoDigiMaker] digis " << nDigis_ << ", failed fits " << nFailedFits_
                << ", chi2 cut " << nChi2Cut_ << std::endl;
    }
    event.put(std::move(recoDigis));
  }

}

DEFINE_ART_MODULE(mu2e::MC2CaloRecoDigiMaker)
