// MC2RawCheck: integrity of the raw calorimeter data, read directly from the DTC fragments of the art file (before any
// reconstruction), in the same way as Offline DAQ/src/CaloDigisFromDTCEvents_module.cc, but keeping also the hits that the
// decoder drops.
//
// Per event : event-window tag of every calorimeter sub-event (DTC) and of every ROC data block (link); monotonicity; ROC blocks
//             present; number of hits; laser events (hits >= laserMinHits).
// Per hit   : CaloDAQUtilities::isHitGood error code (the reason why Offline drops a hit) and the 4-bit ErrorFlags of the hit packet,
//             counted bit by bit, per channel and in time, for kept and for dropped hits;
//             in-payload event-window tag vs the sub-event tag; detector type (CALO / CAPHRI / LASER ...); waveform length, index of
//             the largest sample, saturation; pedestal = mean of the first nPed samples, its rms, and the Baseline field of the packet.
// Time      : slices of sliceMinutes from the event-window tag (windowMicroseconds per window, ASSUMPTION written in the output);
//             per channel and slice: pedestal sums (stability) and hit counts; in laser events, per board and slice: median sample of
//             the peak relative to the median of the event (whole-sample jumps of a board).
// Output    : histograms in the TFileService file, directory MC2RawCheck; read by Offline/MC2Reco/scripts/raw_check_report.py.
// Channel key = 20 * BoardID + ChannelID (raw); hKeyOffline maps a key to its offline SiPM id (+1; 0 = not in the map).
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art_root_io/TFileService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "fhiclcpp/types/Atom.h"
#include "artdaq-core-mu2e/Overlays/DTCEventFragment.hh"
#include "artdaq-core-mu2e/Overlays/Decoders/CalorimeterDataDecoder.hh"
#include "Offline/CaloConditions/inc/CalDAQMap.hh"
#include "Offline/ProditionsService/inc/ProditionsHandle.hh"
#include "Offline/DAQ/inc/CaloDAQUtilities.hh"
#include "Offline/DataProducts/inc/CaloRawSiPMId.hh"
#include "Offline/DataProducts/inc/CaloConst.hh"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include <algorithm>
#include <cmath>
#include <map>
#include <vector>

namespace mu2e {
class MC2RawCheck : public art::EDAnalyzer {
 public:
  struct Config {
    fhicl::Atom<double> windowMicroseconds{fhicl::Name("windowMicroseconds"), fhicl::Comment("length of one event window [us] (time axis)"), 100.};
    fhicl::Atom<double> sliceMinutes{fhicl::Name("sliceMinutes"), fhicl::Comment("time slice for the stability plots [min]"), 10.};
    fhicl::Atom<int> nSlices{fhicl::Name("nSlices"), fhicl::Comment("maximum number of time slices"), 144};
    fhicl::Atom<int> laserMinHits{fhicl::Name("laserMinHits"), fhicl::Comment("an event with at least this many hits is a laser pulse"), 1000};
    fhicl::Atom<int> nPed{fhicl::Name("nPed"), fhicl::Comment("samples used for the pedestal"), 4};
    fhicl::Atom<int> satADC{fhicl::Name("satADC"), fhicl::Comment("saturation [ADC]"), 4090};
  };
  using Parameters = art::EDAnalyzer::Table<Config>;
  explicit MC2RawCheck(Parameters const& c);
  void analyze(art::Event const& e) override;

 private:
  static constexpr int NB = 160, NC = 20, NK = NB * NC, NDTC = 32, NLINK = 8;
  double winUs_, sliceMin_; int nSl_, laserMin_, nPed_, sat_;
  ProditionsHandle<CalDAQMap> daqMap_h_;
  CaloDAQUtilities util_;
  long long tag0_ = -1, lastTag_ = -1; bool mapFilled_ = false;
  TH1D *hEvSlice_, *hLaserSlice_, *hTagStep_, *hNHits_, *hNBlocks_, *hKeyOffline_, *hOcc_, *hOccLaser_, *hPartner_, *hSat_, *hPeakIdx_, *hErrTot_, *hEvFailSlice_;
  TH2D *hBlocks_, *hBlockTagShift_, *hErr_, *hFlags_, *hHitTagShift_, *hDet_, *hNSamp_, *hPedSum_, *hPed2Sum_, *hNoiseSum_, *hPedN_, *hBaseDiff_, *hErrSlice_;
  TH2D *hFlagBitKept_, *hFlagBitDropped_, *hFlagBitSlice_; TH1D *hFlagsTot_, *hBaseline_, *hTagRange_; TH2D *hPartnerDt_, *hOccSlice_, *hRateFine_;
  TH3D* hBoardJump_;
  int slice(long long tag) const { if (tag0_ < 0) return 0; const double t = (tag - tag0_) * winUs_ * 1e-6 / 60. / sliceMin_; return std::clamp(int(std::floor(t)), 0, nSl_ - 1); }
};

MC2RawCheck::MC2RawCheck(Parameters const& c)
    : art::EDAnalyzer{c}, winUs_(c().windowMicroseconds()), sliceMin_(c().sliceMinutes()), nSl_(c().nSlices()), laserMin_(c().laserMinHits()),
      nPed_(c().nPed()), sat_(c().satADC()), util_("MC2RawCheck") {
  art::ServiceHandle<art::TFileService> tfs; auto d = tfs->mkdir("MC2RawCheck");
  const double sl = nSl_ * sliceMin_;
  hEvSlice_ = d.make<TH1D>("hEvSlice", ";time from the first event [min];events", nSl_, 0, sl);
  hLaserSlice_ = d.make<TH1D>("hLaserSlice", ";time [min];laser events", nSl_, 0, sl);
  hEvFailSlice_ = d.make<TH1D>("hEvFailSlice", ";time [min];events with at least one bad hit", nSl_, 0, sl);
  hTagStep_ = d.make<TH1D>("hTagStep", ";event-window tag - tag of the previous event;events", 41, -20.5, 20.5);
  hNHits_ = d.make<TH1D>("hNHits", ";hits per event;events", 300, 0, 3000);
  hNBlocks_ = d.make<TH1D>("hNBlocks", ";calorimeter ROC blocks per event;events", 200, -0.5, 199.5);
  hBlocks_ = d.make<TH2D>("hBlocks", ";DTC id;link;events with this ROC block", NDTC, -0.5, NDTC - 0.5, NLINK, -0.5, NLINK - 0.5);
  hBlockTagShift_ = d.make<TH2D>("hBlockTagShift", ";DTC id x 8 + link;ROC block tag - sub-event tag", NDTC * NLINK, -0.5, NDTC * NLINK - 0.5, 9, -4.5, 4.5);
  hKeyOffline_ = d.make<TH1D>("hKeyOffline", ";channel key;offline SiPM id + 1", NK, -0.5, NK - 0.5);
  hErr_ = d.make<TH2D>("hErr", ";channel key;isHitGood error code (0 = good)", NK, -0.5, NK - 0.5, 10, -0.5, 9.5);
  hErrTot_ = d.make<TH1D>("hErrTot", ";isHitGood error code;hits", 10, -0.5, 9.5);
  hErrSlice_ = d.make<TH2D>("hErrSlice", ";time [min];error code", nSl_, 0, sl, 10, -0.5, 9.5);
  hFlags_ = d.make<TH2D>("hFlags", ";channel key;ErrorFlags (4 bits)", NK, -0.5, NK - 0.5, 16, -0.5, 15.5);
  hHitTagShift_ = d.make<TH2D>("hHitTagShift", ";channel key;hit tag - sub-event tag (mod 2^16)", NK, -0.5, NK - 0.5, 9, -4.5, 4.5);
  hDet_ = d.make<TH2D>("hDet", ";channel key;DetectorID", NK, -0.5, NK - 0.5, 8, -0.5, 7.5);
  hOcc_ = d.make<TH1D>("hOcc", ";channel key;good hits, non-laser events", NK, -0.5, NK - 0.5);
  hOccLaser_ = d.make<TH1D>("hOccLaser", ";channel key;good hits, laser events", NK, -0.5, NK - 0.5);
  hPartner_ = d.make<TH1D>("hPartner", ";channel key;non-laser hits with a good hit of the other SiPM of the crystal in the same event", NK, -0.5, NK - 0.5);
  hPartnerDt_ = d.make<TH2D>("hPartnerDt", ";channel key;peak sample of this SiPM - of the other SiPM (raw, no channel offsets)", NK, -0.5, NK - 0.5, 121, -300.5, 300.5);
  hSat_ = d.make<TH1D>("hSat", ";channel key;saturated hits", NK, -0.5, NK - 0.5);
  hPeakIdx_ = d.make<TH1D>("hPeakIdx", ";index of the largest sample;hits", 128, -0.5, 127.5);
  hNSamp_ = d.make<TH2D>("hNSamp", ";channel key;samples in the waveform", NK, -0.5, NK - 0.5, 128, 0, 1024);
  hPedSum_ = d.make<TH2D>("hPedSum", ";channel key;time [min]", NK, -0.5, NK - 0.5, nSl_, 0, sl);
  hPed2Sum_ = d.make<TH2D>("hPed2Sum", ";channel key;time [min]", NK, -0.5, NK - 0.5, nSl_, 0, sl);
  hNoiseSum_ = d.make<TH2D>("hNoiseSum", ";channel key;time [min]", NK, -0.5, NK - 0.5, nSl_, 0, sl);   // sum over hits of the variance of the first nPed samples
  hFlagBitKept_ = d.make<TH2D>("hFlagBitKept", ";channel key;ErrorFlags bit set (hit kept by isHitGood)", NK, -0.5, NK - 0.5, 4, -0.5, 3.5);
  hFlagBitDropped_ = d.make<TH2D>("hFlagBitDropped", ";channel key;ErrorFlags bit set (hit dropped by isHitGood)", NK, -0.5, NK - 0.5, 4, -0.5, 3.5);
  hFlagBitSlice_ = d.make<TH2D>("hFlagBitSlice", ";time [min];ErrorFlags bit set", nSl_, 0, sl, 4, -0.5, 3.5);
  hFlagsTot_ = d.make<TH1D>("hFlagsTot", ";ErrorFlags value;hits", 16, -0.5, 15.5);
  hOccSlice_ = d.make<TH2D>("hOccSlice", ";channel key;time [min]", NK, -0.5, NK - 0.5, nSl_, 0, sl);   // good hits of non-laser events per channel and slice
  // counts in fine 10 s bins (up to 24 h), for the rate plots: 0 events, 1 laser events, 2 good hits disk 0, 3 good hits disk 1, 4 good LYSO hits, 5 dropped hits (non-laser events for 2-4)
  hRateFine_ = d.make<TH2D>("hRateFine", ";time [s] (10 s bins);category", 8640, 0, 86400, 6, -0.5, 5.5);
  hTagRange_ = d.make<TH1D>("hTagRange", "first and last event-window tag of the job (bin 1, 2)", 2, 0, 2);
  hBaseline_ = d.make<TH1D>("hBaseline", ";Baseline field of the hit packet [ADC];hits", 4096, -0.5, 4095.5);
  hPedN_ = d.make<TH2D>("hPedN", ";channel key;time [min]", NK, -0.5, NK - 0.5, nSl_, 0, sl);
  hBaseDiff_ = d.make<TH2D>("hBaseDiff", ";channel key;packet Baseline - mean of the first samples [ADC]", NK, -0.5, NK - 0.5, 81, -40.5, 40.5);
  hBoardJump_ = d.make<TH3D>("hBoardJump", ";board;time [min];peak sample - median of the laser event [samples]", NB, -0.5, NB - 0.5, nSl_, 0, sl, 41, -20.5, 20.5);
  TH1D* hp = d.make<TH1D>("hParams", "windowMicroseconds, sliceMinutes, nSlices, laserMinHits, nPed, satADC", 6, 0, 6);
  hp->SetBinContent(1, winUs_); hp->SetBinContent(2, sliceMin_); hp->SetBinContent(3, nSl_); hp->SetBinContent(4, laserMin_); hp->SetBinContent(5, nPed_); hp->SetBinContent(6, sat_);
}

void MC2RawCheck::analyze(art::Event const& event) {
  CalDAQMap const& daqMap = daqMap_h_.get(event.id());
  if (!mapFilled_) {
    for (int k = 0; k < NK; ++k) { CaloRawSiPMId raw(k / NC, k % NC); const auto off = daqMap.offlineId(raw).id(); if (off < CaloConst::_nChannel) hKeyOffline_->SetBinContent(k + 1, off + 1); }
    mapFilled_ = true;
  }
  auto frags = util_.getFragments(event);
  // first pass: the event-window tag of the event (first calorimeter sub-event) -> time slice
  long long evTag = -1;
  for (const auto& frag : frags) { DTCEventFragment ef(frag); for (auto& sub : ef.getSubsystemData(DTCLib::DTC_Subsystem::DTC_Subsystem_Calorimeter)) { evTag = (long long)sub.GetEventWindowTag().GetEventWindowTag(true); break; } if (evTag >= 0) break; }
  if (evTag < 0) return;
  if (tag0_ < 0) tag0_ = evTag;
  if (lastTag_ >= 0) hTagStep_->Fill(std::clamp<long long>(evTag - lastTag_, -20, 20));
  lastTag_ = evTag;
  if (hTagRange_->GetBinContent(1) == 0 || evTag < hTagRange_->GetBinContent(1)) hTagRange_->SetBinContent(1, double(evTag));
  if (evTag > hTagRange_->GetBinContent(2)) hTagRange_->SetBinContent(2, double(evTag));
  const double tmin = (slice(evTag) + 0.5) * sliceMin_;
  struct HitInfo { int key, peak; };
  std::vector<HitInfo> good; int nAll = 0, nBlocks = 0; bool anyBad = false;
  for (const auto& frag : frags) {
    DTCEventFragment ef(frag);
    for (auto& sub : ef.getSubsystemData(DTCLib::DTC_Subsystem::DTC_Subsystem_Calorimeter)) {
      const long long st = (long long)sub.GetEventWindowTag().GetEventWindowTag(true);
      CalorimeterDataDecoder dec(sub); const int dtc = sub.GetDTCID();
      for (size_t ib = 0; ib < dec.block_count(); ++ib) {
        auto blk = dec.dataAtBlockIndex(ib); if (!blk) continue; auto h = blk->GetHeader(); const int link = int(h->GetLinkID());
        ++nBlocks; hBlocks_->Fill(dtc, link);
        hBlockTagShift_->Fill(dtc * NLINK + link, std::clamp<long long>((long long)h->GetEventWindowTag().GetEventWindowTag(true) - st, -4, 4));
        auto data = dec.GetCalorimeterHitData(ib); if (!data) continue;
        for (auto& hp : *data) {
          ++nAll; const auto& p = hp.first; const auto& wf = hp.second; const int key = int(p.BoardID) * NC + int(p.ChannelID);
          const int code = int(util_.isHitGood(hp)); const int flags = int(p.ErrorFlags);
          hErrTot_->Fill(code); hErrSlice_->Fill(tmin, code); hFlagsTot_->Fill(flags); if (code) anyBad = true;
          for (int bit = 0; bit < 4; ++bit) if (flags & (1 << bit)) hFlagBitSlice_->Fill(tmin, bit);
          if (key < 0 || key >= NK) continue;
          hErr_->Fill(key, code); hFlags_->Fill(key, flags); hDet_->Fill(key, int(p.DetectorID));
          for (int bit = 0; bit < 4; ++bit) if (flags & (1 << bit)) (code ? hFlagBitDropped_ : hFlagBitKept_)->Fill(key, bit);
          int sh = int(p.InPayloadEventWindowTag) - int(st & 0xFFFF); if (sh > 32768) sh -= 65536; if (sh < -32768) sh += 65536;
          hHitTagShift_->Fill(key, std::clamp(sh, -4, 4));
          if (code) continue;
          hNSamp_->Fill(key, wf.size()); hPeakIdx_->Fill(p.IndexOfMaxDigitizerSample);
          int mx = 0; for (auto s : wf) mx = std::max<int>(mx, s); if (mx >= sat_) hSat_->Fill(key);
          if ((int)wf.size() > nPed_) {
            double s1 = 0, s2 = 0; for (int q = 0; q < nPed_; ++q) { s1 += wf[q]; s2 += double(wf[q]) * wf[q]; }
            const double m = s1 / nPed_, var = std::max(0., s2 / nPed_ - m * m);
            hBaseDiff_->Fill(key, double(p.Baseline) - m); hBaseline_->Fill(int(p.Baseline));
            hPedSum_->Fill(key, tmin, m); hPed2Sum_->Fill(key, tmin, m * m); hNoiseSum_->Fill(key, tmin, var); hPedN_->Fill(key, tmin);
          }
          good.push_back({key, int(p.Time) + int(p.IndexOfMaxDigitizerSample)});
        }
      }
    }
  }
  const double tsec = (evTag - tag0_) * winUs_ * 1e-6;
  hRateFine_->Fill(tsec, 0); hRateFine_->Fill(tsec, 5, double(nAll - (int)good.size()));
  hEvSlice_->Fill(tmin); hNHits_->Fill(nAll); hNBlocks_->Fill(nBlocks); if (anyBad) hEvFailSlice_->Fill(tmin);
  const bool laser = nAll >= laserMin_;
  if (laser) { hLaserSlice_->Fill(tmin); hRateFine_->Fill(tsec, 1); }
  std::map<int, int> peakOf; for (auto& h : good) peakOf[h.key] = h.peak;
  if (laser && !peakOf.empty()) {
    std::vector<int> pk; for (auto& kv : peakOf) pk.push_back(kv.second);
    std::nth_element(pk.begin(), pk.begin() + pk.size() / 2, pk.end()); const int med = pk[pk.size() / 2];
    std::map<int, std::vector<int>> byBoard; for (auto& kv : peakOf) byBoard[kv.first / NC].push_back(kv.second - med);
    for (auto& [b, v] : byBoard) { std::nth_element(v.begin(), v.begin() + v.size() / 2, v.end()); hBoardJump_->Fill(b, tmin, std::clamp(v[v.size() / 2], -20, 20)); }
  }
  std::map<int, int> peakOfOffline; for (auto& kv : peakOf) { const int off = int(hKeyOffline_->GetBinContent(kv.first + 1)) - 1; if (off >= 0) peakOfOffline[off] = kv.second; }
  for (auto& kv : peakOf) {
    (laser ? hOccLaser_ : hOcc_)->Fill(kv.first); if (!laser) hOccSlice_->Fill(kv.first, tmin);
    if (!laser) { const int o_ = int(hKeyOffline_->GetBinContent(kv.first + 1)) - 1; const bool ly = hDet_->GetBinContent(kv.first + 1, 2) > 0;
      if (ly) hRateFine_->Fill(tsec, 4); else if (o_ >= 0) hRateFine_->Fill(tsec, (o_ / 2) >= 674 ? 3 : 2); }
    if (laser) continue;
    const int off = int(hKeyOffline_->GetBinContent(kv.first + 1)) - 1; if (off < 0) continue;
    auto it = peakOfOffline.find(off ^ 1);   // the other SiPM of the crystal: offline ids 2c and 2c + 1
    if (it != peakOfOffline.end()) { hPartner_->Fill(kv.first); hPartnerDt_->Fill(kv.first, std::clamp(kv.second - it->second, -300, 300)); }
  }
}
}  // namespace mu2e
DEFINE_ART_MODULE(mu2e::MC2RawCheck)
