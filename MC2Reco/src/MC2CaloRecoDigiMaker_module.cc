//
// MC2CaloRecoDigiMaker: CaloDigi -> CaloRecoDigi with the MC2 template fit (MC2TemplateFit) and per-SiPM time offsets.
//
// It replaces mu2e::CaloRecoDigiMaker (Offline/CaloReco) for the MC2 cosmic-ray runs and writes the same data product, so
// that CaloHitMaker and everything downstream run unchanged: set CaloHitMaker.caloRecoDigiCollection to this module label.
// Both reconstructions can run in the same job (two collections) and be compared on the same events.
//
//   time    = 5 ns * t0(CaloDigi) + fitted peak time - [proton bunch time - calo-tracker offset] - tcorr(SiPMID)
//             (the bracket only if UseProtonBunchTime, exactly as in CaloRecoDigiMaker)
//   eDep    = A * CalCalib.ADC2MeV(SiPMID), as in CaloRecoDigiMaker
//   timeErr = sqrt((timeErrA / A)^2 + timeErrC^2), from the MC2 single-SiPM resolution on cosmics
//   tcorr   = per-SiPM offset in ns: with useCalCalibTimeOffset, the CalCalib time offset (Offline CalTimeCalib table, e.g. from
//             DbService.textFile), exactly as CaloRecoDigiMaker; otherwise from timeOffsetFile (text: "SiPMID tcorr" per line,
//             # for comments), 0 if the file is empty or the SiPM is not listed. Not both (exception).
//   perChannel: one template per SiPM from the histograms perChannelFormat (board, channel of the raw id, e.g. tmpl_b012_c07) of
//             templateFile, as written by mc2-ana laser_template_stream.C; SiPMs without one use the common graph graphName.
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
        using Name = fhicl::Name; using Comment = fhicl::Comment;
        fhicl::Atom<art::InputTag> caloDigiCollection { Name("caloDigiCollection"), Comment("CaloDigi module label") };
        fhicl::Atom<std::string>   templateFile       { Name("templateFile"),       Comment("production-format template file (TGraph gtempl)") };
        fhicl::Atom<std::string>   graphName          { Name("graphName"),          Comment("name of the template graph"), "gtempl" };
        fhicl::Atom<double>        fitLo              { Name("fitLo"),              Comment("start of the fit window from the fitted peak (ns)"), -80. };
        fhicl::Atom<double>        fitHi              { Name("fitHi"),              Comment("end of the fit window from the fitted peak (ns)"), -31. };
        fhicl::Atom<double>        noise              { Name("noise"),              Comment("noise per sample (ADC), for chi2/ndf"), 3.6 };
        fhicl::Atom<int>           nPedSamples        { Name("nPedSamples"),        Comment("samples used for the pedestal"), 4 };
        fhicl::Atom<double>        maxBaseline        { Name("maxBaseline"),        Comment("limit of the residual baseline (ADC)"), 10. };
        fhicl::Atom<int>           maxPre             { Name("maxPre"),             Comment("samples before the largest one allowed in the fit (<0: all)"), 16 };
        fhicl::Atom<int>           shapeMode          { Name("shapeMode"),          Comment("0 fixed template; amplitude model in templateFile: 1 shape tied to the fitted amplitude, 2 free shape, 3 constrained shape"), 0 };
        fhicl::Atom<double>        shapeSigma         { Name("shapeSigma"),         Comment("shapeMode 3: spread of the shape between hits (units of ln A)"), 0.5 };
        fhicl::Atom<int>           satADC             { Name("satADC"),             Comment("raw samples >= satADC are left out of the fit"), 4090 };
        fhicl::Atom<double>        minAmplitude       { Name("minAmplitude"),       Comment("minimum fitted amplitude (ADC)"), 0. };
        fhicl::Atom<double>        maxChi2Cut         { Name("maxChi2Cut"),         Comment("maximum chi2/ndf to keep the reco digi"), 1e9 };
        fhicl::Atom<double>        timeErrA           { Name("timeErrA"),           Comment("time error term a/A (ns ADC)"), 110. };
        fhicl::Atom<double>        timeErrC           { Name("timeErrC"),           Comment("time error constant term (ns)"), 0.18 };
        fhicl::Atom<std::string>   timeOffsetFile     { Name("timeOffsetFile"),     Comment("per-SiPM time offsets: SiPMID tcorr[ns]; empty = none"), "" };
        fhicl::Atom<bool>          useCalCalibTime    { Name("useCalCalibTimeOffset"), Comment("subtract the CalCalib time offset (CalTimeCalib table), as CaloRecoDigiMaker"), false };
        fhicl::Atom<art::InputTag> pbtTag             { Name("ProtonBunchTimeTag"), Comment("ProtonBunchTime producer"), art::InputTag("EWMProducer") };
        fhicl::Atom<bool>          usePBT             { Name("UseProtonBunchTime"), Comment("subtract the proton bunch time as CaloRecoDigiMaker"), false };
        fhicl::Atom<int>           diagLevel          { Name("diagLevel"),          Comment("diagnosis level"), 0 };
        fhicl::Atom<bool>          perChannel         { Name("perChannel"),         Comment("one template per SiPM from the histograms perChannelFormat of templateFile"), false };
        fhicl::Atom<std::string>   perChannelFormat   { Name("perChannelFormat"),   Comment("printf format of the per-channel histogram names (board, channel)"), "tmpl_b%03d_c%02d" };
      };

      explicit MC2CaloRecoDigiMaker(const art::EDProducer::Table<Config>& config);
      void produce(art::Event& event) override;

    private:
      art::ProductToken<CaloDigiCollection> caloDigisToken_;
      art::InputTag pbtTag_;
      bool usePBT_;
      art::ProductToken<ProtonBunchTime> pbtToken_;
      MC2TemplateFit fitter_;
      double minAmplitude_, maxChi2Cut_, timeErrA_, timeErrC_;
      int diagLevel_;
      bool useCalCalibTime_;
      std::unordered_map<int, double> tcorr_;
      ProditionsHandle<CalCalib> calCalibHandle_;
      ProditionsHandle<EventTiming> eventTimingHandle_;
      long nDigi_ = 0, nFail_ = 0, nCut_ = 0;
      static MC2TemplateFit::Config fitConfig(const Config& c);
      bool perChannel_; std::string perChFormat_;
      ProditionsHandle<CalDAQMap> daqMapHandle_;
      std::vector<std::unique_ptr<MC2TemplateFit>> perCh_;   // by offline SiPM id; empty = common template
      bool perChLoaded_ = false;
      void loadPerChannel(const art::Event& event, const MC2TemplateFit::Config& fc);
      MC2TemplateFit::Config fc_;
  };

  MC2TemplateFit::Config MC2CaloRecoDigiMaker::fitConfig(const Config& c) {
    MC2TemplateFit::Config f;
    ConfigFileLookupPolicy resolve;
    f.templateFile = resolve(c.templateFile()); f.graphName = c.graphName(); f.fitLo = c.fitLo(); f.fitHi = c.fitHi();
    f.noise = c.noise(); f.nPedSamples = c.nPedSamples(); f.maxBaseline = c.maxBaseline(); f.satADC = c.satADC(); f.maxPre = c.maxPre(); f.shapeMode = c.shapeMode(); f.shapeSigma = c.shapeSigma();
    return f;
  }

  MC2CaloRecoDigiMaker::MC2CaloRecoDigiMaker(const art::EDProducer::Table<Config>& config) :
    EDProducer{config},
    caloDigisToken_(consumes<CaloDigiCollection>(config().caloDigiCollection())),
    pbtTag_(config().pbtTag()), usePBT_(config().usePBT()),
    fitter_(fitConfig(config())),
    minAmplitude_(config().minAmplitude()), maxChi2Cut_(config().maxChi2Cut()),
    timeErrA_(config().timeErrA()), timeErrC_(config().timeErrC()), diagLevel_(config().diagLevel()), useCalCalibTime_(config().useCalCalibTime()),
    perChannel_(config().perChannel()), perChFormat_(config().perChannelFormat()), fc_(fitConfig(config()))
  {
    produces<CaloRecoDigiCollection>();
    if (usePBT_) pbtToken_ = consumes<ProtonBunchTime>(pbtTag_);
    if (useCalCalibTime_ && !config().timeOffsetFile().empty())
      throw cet::exception("MC2CALORECO") << "timeOffsetFile and useCalCalibTimeOffset together would subtract two offsets: use one\n";
    if (!config().timeOffsetFile().empty()) {
      ConfigFileLookupPolicy resolve;
      const std::string fn = resolve(config().timeOffsetFile());
      std::ifstream in(fn);
      if (!in) throw cet::exception("MC2CALORECO") << "cannot open timeOffsetFile " << fn << "\n";
      std::string line;
      while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::istringstream ss(line); int id; double t;
        if (ss >> id >> t) tcorr_[id] = t;
      }
      std::cout << "[MC2CaloRecoDigiMaker] " << tcorr_.size() << " SiPM time offsets from " << fn << std::endl;
    }
  }

  // per-channel templates: histogram of (board, channel) of each offline SiPM -> graph of the bin centres -> its own fitter
  void MC2CaloRecoDigiMaker::loadPerChannel(const art::Event& event, const MC2TemplateFit::Config& fc) {
    perChLoaded_ = true; if (!perChannel_) return;
    const CalDAQMap& dm = daqMapHandle_.get(event.id());
    std::unique_ptr<TFile> f(TFile::Open(fc.templateFile.c_str(), "READ"));
    if (!f || f->IsZombie()) throw cet::exception("MC2CALORECO") << "cannot open " << fc.templateFile << "\n";
    perCh_.resize(CaloConst::_nChannel); int n = 0;
    for (int i = 0; i < (int)CaloConst::_nChannel; ++i) {
      const auto raw = dm.rawId(CaloSiPMId(i)); char name[64];
      snprintf(name, sizeof(name), perChFormat_.c_str(), int(raw.dirac()), int(raw.ROCchannel()));
      auto* h = dynamic_cast<TH1*>(f->Get(name)); if (!h || h->GetEntries() <= 0) continue;
      TGraph g; for (int b = 1; b <= h->GetNbinsX(); ++b) g.SetPoint(g.GetN(), h->GetBinCenter(b), h->GetBinContent(b));
      if (g.GetN() < 10) continue;
      perCh_[i] = std::make_unique<MC2TemplateFit>(fc, g); ++n;
    }
    std::cout << "[MC2CaloRecoDigiMaker] per-channel templates for " << n << " SiPMs from " << fc.templateFile << " (common template for the others)" << std::endl;
  }

  void MC2CaloRecoDigiMaker::produce(art::Event& event) {
    if (!perChLoaded_) loadPerChannel(event, fc_);
    const auto& digisH = event.getValidHandle(caloDigisToken_);
    const auto& digis = *digisH;
    auto out = std::make_unique<CaloRecoDigiCollection>();

    double pbtOffset = 0;
    if (usePBT_) {
      pbtOffset = event.getValidHandle(pbtToken_)->pbtime_;
      pbtOffset -= eventTimingHandle_.get(event.id()).caloTrackerTimeOffset();
    }
    const CalCalib& calCalib = calCalibHandle_.get(event.id());

    for (size_t i = 0; i < digis.size(); ++i) {
      const CaloDigi& d = digis[i]; ++nDigi_;
      const int sid = d.SiPMID();
      const MC2TemplateFit& ft = (sid >= 0 && sid < (int)perCh_.size() && perCh_[sid]) ? *perCh_[sid] : fitter_;
      const auto r = ft.fit(d.waveform());
      if (!r.ok || r.amplitude < minAmplitude_) { ++nFail_; continue; }
      if (r.chi2ndf > maxChi2Cut_) { ++nCut_; continue; }
      const int id = d.SiPMID();
      const auto it = tcorr_.find(id);
      const double tcorr = useCalCalibTime_ ? calCalib.timeoffset(id) : (it == tcorr_.end() ? 0. : it->second);
      const double time = d.t0() * CaloConst::_digitizationPeriod + r.time - pbtOffset - tcorr;
      const double adc2MeV = calCalib.ADC2MeV(id);
      const double eDep = r.amplitude * adc2MeV, eDepErr = std::sqrt(std::max(r.amplitude, 1.)) * adc2MeV;
      const double timeErr = std::hypot(timeErrA_ / std::max(r.amplitude, 1.), timeErrC_);
      art::Ptr<CaloDigi> ptr(digisH, i);
      out->emplace_back(CaloRecoDigi(ptr, eDep, eDepErr, time, timeErr, r.chi2ndf * r.ndf, r.ndf, false));
      if (diagLevel_ > 2) std::cout << "[MC2CaloRecoDigiMaker] SiPM " << id << " A=" << r.amplitude << " t=" << time << " chi2/ndf=" << r.chi2ndf << std::endl;
    }
    if (diagLevel_ > 0 && event.id().event() % 10000 == 0)
      std::cout << "[MC2CaloRecoDigiMaker] digis " << nDigi_ << ", failed fits " << nFail_ << ", chi2 cut " << nCut_ << std::endl;
    event.put(std::move(out));
  }

}

DEFINE_ART_MODULE(mu2e::MC2CaloRecoDigiMaker)
