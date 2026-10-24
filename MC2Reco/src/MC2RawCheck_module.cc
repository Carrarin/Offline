//
// MC2RawCheck: integrity of the raw calorimeter data, read directly from the DTC fragments of the art file (before any
// reconstruction), in the same way as Offline DAQ/src/CaloDigisFromDTCEvents_module.cc, but keeping also the hits that
// the decoder drops.
//
// Input: the DTC fragments of the event (CaloDAQUtilities::getFragments) and the CalDAQMap condition (raw -> offline).
//
// Output: histograms in the TFileService file, directory MC2RawCheck, read by
// Offline/MC2Reco/scripts/raw_check_report.py (the names must not change).
// Channel key = 20 * BoardID + ChannelID (raw).
//
//   Per event : event-window tag of every calorimeter sub-event (DTC) and of every ROC data block (link);
//               monotonicity (hTagStep); ROC blocks present (hBlocks, hNBlocks, hBlockTagShift);
//               number of hits (hNHits); laser events, i.e. events with at least laserMinHits hits (hLaserSlice).
//   Per hit   : CaloDAQUtilities::isHitGood error code, i.e. the reason why Offline drops a hit
//               (hErr, hErrTot, hErrSlice), and the 4-bit ErrorFlags of the hit packet, counted bit by bit, per
//               channel and in time, for kept and for dropped hits
//               (hFlags, hFlagsTot, hFlagBitKept, hFlagBitDropped, hFlagBitSlice);
//               in-payload event-window tag vs the sub-event tag (hHitTagShift);
//               detector type CALO / CAPHRI / LASER ... (hDet); waveform length (hNSamp), index of the largest
//               sample (hPeakIdx), saturation (hSat); pedestal = mean of the first nPed samples, its rms, and the
//               Baseline field of the packet (hPedSum, hPed2Sum, hNoiseSum, hPedN, hBaseDiff, hBaseline).
//   Time      : slices of sliceMinutes from the event-window tag (windowMicroseconds per window: an ASSUMPTION,
//               written in the output in hParams). Per channel and slice: pedestal sums (stability) and hit
//               counts (hOccSlice); in laser events, per board and slice: median sample of the peak relative to
//               the median of the event (hBoardJump: whole-sample jumps of a board).
//               Fine 10 s bins for the rate plots (hRateFine).
//   Pairs     : non-laser hits with a good hit of the other SiPM of the same crystal (hPartner) and the
//               difference of their peak samples (hPartnerDt).
//   Maps      : hKeyOffline maps a channel key to its offline SiPM id + 1 (0 = not in the map);
//               hTagRange holds the first and the last event-window tag of the job.
//
// Steps per event:
//   1. at the first event, the map channel key -> offline SiPM id (hKeyOffline)
//   2. event-window tag of the event (first calorimeter sub-event); events without one are skipped
//   3. tag histograms (step from the previous event, first and last tag) and the time slice of the event
//   4. loop over the sub-events, ROC blocks and hits: block and hit checks, list of the good hits
//   5. per-event counts (events, hits, blocks, laser events, rates)
//   6. laser events: per-board jumps of the peak sample
//   7. occupancy, rates per category and pairs of SiPMs of the same crystal
//
// Parameters: windowMicroseconds, sliceMinutes, nSlices (time axis); laserMinHits; nPed (pedestal samples); satADC.
//
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
#include <cstdint>
#include <map>
#include <utility>
#include <vector>

namespace mu2e {

  class MC2RawCheck : public art::EDAnalyzer {
    public:
      struct Config {
        using Name = fhicl::Name;
        using Comment = fhicl::Comment;

        // Time axis
        fhicl::Atom<double> windowMicroseconds{Name("windowMicroseconds"),
            Comment("length of one event window [us] (time axis)"), 100.};
        fhicl::Atom<double> sliceMinutes{Name("sliceMinutes"),
            Comment("time slice for the stability plots [min]"), 10.};
        fhicl::Atom<int> nSlices{Name("nSlices"), Comment("maximum number of time slices"), 144};

        // Event and hit classification
        fhicl::Atom<int> laserMinHits{Name("laserMinHits"),
            Comment("an event with at least this many hits is a laser pulse"), 1000};
        fhicl::Atom<int> nPed{Name("nPed"), Comment("samples used for the pedestal"), 4};
        fhicl::Atom<int> satADC{Name("satADC"), Comment("saturation [ADC]"), 4090};
      };
      using Parameters = art::EDAnalyzer::Table<Config>;

      explicit MC2RawCheck(const Parameters& config);
      void analyze(const art::Event& event) override;

    private:
      // One raw hit as decoded from a ROC block: packet and waveform.
      using CaloHitData = std::pair<CalorimeterDataDecoder::CalorimeterHitDataPacket, std::vector<uint16_t>>;

      // A hit that passed isHitGood: channel key and peak sample (packet Time + index of the largest sample).
      struct GoodHit {
        int key;
        int peakSample;
      };

      // What the hit loop collects for the per-event histograms.
      struct EventHits {
        std::vector<GoodHit> goodHits;
        int nHits = 0;              // all hits, good and bad
        int nBlocks = 0;            // ROC blocks present
        bool anyBadHit = false;
      };

      // Raw channel numbering
      static constexpr int kNBoards = 160;
      static constexpr int kNChannelsPerBoard = 20;
      static constexpr int kNKeys = kNBoards * kNChannelsPerBoard;
      static constexpr int kNDTCs = 32;
      static constexpr int kNLinks = 8;

      void bookHistograms();
      void fillKeyToOfflineMap(const CalDAQMap& daqMap);
      long long eventWindowTag(const artdaq::Fragments& fragments) const;
      void fillTagHistograms(long long eventTag);
      int timeSlice(long long tag) const;
      EventHits checkRocBlocksAndHits(const artdaq::Fragments& fragments, double sliceTime);
      void checkHit(const CaloHitData& hit, long long subEventTag, double sliceTime, EventHits& eventHits);
      void fillPedestal(int key, const CalorimeterDataDecoder::CalorimeterHitDataPacket& packet,
                        const std::vector<uint16_t>& waveform, double sliceTime);
      bool fillEventCounts(const EventHits& eventHits, double sliceTime, double secondsFromStart);
      void fillBoardJumps(const std::map<int, int>& peakByKey, double sliceTime);
      void fillOccupancyAndPartners(const std::map<int, int>& peakByKey, bool isLaser, double sliceTime,
                                    double secondsFromStart);
      int offlineIdOfKey(int key) const;

      // Parameters
      double windowMicroseconds_;
      double sliceMinutes_;
      int nSlices_;
      int laserMinHits_;
      int nPedSamples_;
      int satADC_;

      // Conditions and decoding
      ProditionsHandle<CalDAQMap> daqMapHandle_;
      CaloDAQUtilities daqUtilities_;

      // State
      long long firstTag_ = -1;       // event-window tag of the first event: origin of the time axis
      long long lastTag_ = -1;        // tag of the previous event
      bool keyMapFilled_ = false;

      // Output histograms (member name = histogram name + _), see the header comment
      TH1D* hEvSlice_;
      TH1D* hLaserSlice_;
      TH1D* hEvFailSlice_;
      TH1D* hTagStep_;
      TH1D* hNHits_;
      TH1D* hNBlocks_;
      TH2D* hBlocks_;
      TH2D* hBlockTagShift_;
      TH1D* hKeyOffline_;
      TH2D* hErr_;
      TH1D* hErrTot_;
      TH2D* hErrSlice_;
      TH2D* hFlags_;
      TH2D* hHitTagShift_;
      TH2D* hDet_;
      TH1D* hOcc_;
      TH1D* hOccLaser_;
      TH1D* hPartner_;
      TH2D* hPartnerDt_;
      TH1D* hSat_;
      TH1D* hPeakIdx_;
      TH2D* hNSamp_;
      TH2D* hPedSum_;
      TH2D* hPed2Sum_;
      TH2D* hNoiseSum_;
      TH2D* hFlagBitKept_;
      TH2D* hFlagBitDropped_;
      TH2D* hFlagBitSlice_;
      TH1D* hFlagsTot_;
      TH2D* hOccSlice_;
      TH2D* hRateFine_;
      TH1D* hTagRange_;
      TH1D* hBaseline_;
      TH2D* hPedN_;
      TH2D* hBaseDiff_;
      TH3D* hBoardJump_;
  };

  MC2RawCheck::MC2RawCheck(const Parameters& config) :
    art::EDAnalyzer{config},
    windowMicroseconds_(config().windowMicroseconds()),
    sliceMinutes_(config().sliceMinutes()),
    nSlices_(config().nSlices()),
    laserMinHits_(config().laserMinHits()),
    nPedSamples_(config().nPed()),
    satADC_(config().satADC()),
    daqUtilities_("MC2RawCheck")
  {
    bookHistograms();
  }

  // All histograms, in the booking order of the output file; hParams stores the parameters for the report.
  void MC2RawCheck::bookHistograms() {
    art::ServiceHandle<art::TFileService> tfs;
    auto dir = tfs->mkdir("MC2RawCheck");
    const double maxMinutes = nSlices_ * sliceMinutes_;
    const int nK = kNKeys;

    // Events, tags and ROC blocks
    hEvSlice_ = dir.make<TH1D>("hEvSlice", ";time from the first event [min];events", nSlices_, 0, maxMinutes);
    hLaserSlice_ = dir.make<TH1D>("hLaserSlice", ";time [min];laser events", nSlices_, 0, maxMinutes);
    hEvFailSlice_ = dir.make<TH1D>("hEvFailSlice", ";time [min];events with at least one bad hit",
                                   nSlices_, 0, maxMinutes);
    hTagStep_ = dir.make<TH1D>("hTagStep", ";event-window tag - tag of the previous event;events", 41, -20.5, 20.5);
    hNHits_ = dir.make<TH1D>("hNHits", ";hits per event;events", 300, 0, 3000);
    hNBlocks_ = dir.make<TH1D>("hNBlocks", ";calorimeter ROC blocks per event;events", 200, -0.5, 199.5);
    hBlocks_ = dir.make<TH2D>("hBlocks", ";DTC id;link;events with this ROC block",
                              kNDTCs, -0.5, kNDTCs - 0.5, kNLinks, -0.5, kNLinks - 0.5);
    hBlockTagShift_ = dir.make<TH2D>("hBlockTagShift", ";DTC id x 8 + link;ROC block tag - sub-event tag",
                                     kNDTCs * kNLinks, -0.5, kNDTCs * kNLinks - 0.5, 9, -4.5, 4.5);
    hKeyOffline_ = dir.make<TH1D>("hKeyOffline", ";channel key;offline SiPM id + 1", nK, -0.5, nK - 0.5);

    // Hit errors, flags, tags and detector type
    hErr_ = dir.make<TH2D>("hErr", ";channel key;isHitGood error code (0 = good)",
                           nK, -0.5, nK - 0.5, 10, -0.5, 9.5);
    hErrTot_ = dir.make<TH1D>("hErrTot", ";isHitGood error code;hits", 10, -0.5, 9.5);
    hErrSlice_ = dir.make<TH2D>("hErrSlice", ";time [min];error code", nSlices_, 0, maxMinutes, 10, -0.5, 9.5);
    hFlags_ = dir.make<TH2D>("hFlags", ";channel key;ErrorFlags (4 bits)", nK, -0.5, nK - 0.5, 16, -0.5, 15.5);
    hHitTagShift_ = dir.make<TH2D>("hHitTagShift", ";channel key;hit tag - sub-event tag (mod 2^16)",
                                   nK, -0.5, nK - 0.5, 9, -4.5, 4.5);
    hDet_ = dir.make<TH2D>("hDet", ";channel key;DetectorID", nK, -0.5, nK - 0.5, 8, -0.5, 7.5);

    // Occupancy and pairs of SiPMs of the same crystal
    hOcc_ = dir.make<TH1D>("hOcc", ";channel key;good hits, non-laser events", nK, -0.5, nK - 0.5);
    hOccLaser_ = dir.make<TH1D>("hOccLaser", ";channel key;good hits, laser events", nK, -0.5, nK - 0.5);
    hPartner_ = dir.make<TH1D>("hPartner",
        ";channel key;non-laser hits with a good hit of the other SiPM of the crystal in the same event",
        nK, -0.5, nK - 0.5);
    hPartnerDt_ = dir.make<TH2D>("hPartnerDt",
        ";channel key;peak sample of this SiPM - of the other SiPM (raw, no channel offsets)",
        nK, -0.5, nK - 0.5, 121, -300.5, 300.5);

    // Waveforms and pedestals
    hSat_ = dir.make<TH1D>("hSat", ";channel key;saturated hits", nK, -0.5, nK - 0.5);
    hPeakIdx_ = dir.make<TH1D>("hPeakIdx", ";index of the largest sample;hits", 128, -0.5, 127.5);
    hNSamp_ = dir.make<TH2D>("hNSamp", ";channel key;samples in the waveform", nK, -0.5, nK - 0.5, 128, 0, 1024);
    hPedSum_ = dir.make<TH2D>("hPedSum", ";channel key;time [min]", nK, -0.5, nK - 0.5, nSlices_, 0, maxMinutes);
    hPed2Sum_ = dir.make<TH2D>("hPed2Sum", ";channel key;time [min]", nK, -0.5, nK - 0.5, nSlices_, 0, maxMinutes);
    // sum over hits of the variance of the first nPed samples
    hNoiseSum_ = dir.make<TH2D>("hNoiseSum", ";channel key;time [min]", nK, -0.5, nK - 0.5, nSlices_, 0, maxMinutes);

    // ErrorFlags bit by bit
    hFlagBitKept_ = dir.make<TH2D>("hFlagBitKept", ";channel key;ErrorFlags bit set (hit kept by isHitGood)",
                                   nK, -0.5, nK - 0.5, 4, -0.5, 3.5);
    hFlagBitDropped_ = dir.make<TH2D>("hFlagBitDropped", ";channel key;ErrorFlags bit set (hit dropped by isHitGood)",
                                      nK, -0.5, nK - 0.5, 4, -0.5, 3.5);
    hFlagBitSlice_ = dir.make<TH2D>("hFlagBitSlice", ";time [min];ErrorFlags bit set",
                                    nSlices_, 0, maxMinutes, 4, -0.5, 3.5);
    hFlagsTot_ = dir.make<TH1D>("hFlagsTot", ";ErrorFlags value;hits", 16, -0.5, 15.5);

    // good hits of non-laser events per channel and slice
    hOccSlice_ = dir.make<TH2D>("hOccSlice", ";channel key;time [min]", nK, -0.5, nK - 0.5, nSlices_, 0, maxMinutes);
    // counts in fine 10 s bins (up to 24 h), for the rate plots, per category: 0 events, 1 laser events,
    // 2 good hits disk 0, 3 good hits disk 1, 4 good LYSO hits, 5 dropped hits (categories 2-4 in non-laser events)
    hRateFine_ = dir.make<TH2D>("hRateFine", ";time [s] (10 s bins);category", 8640, 0, 86400, 6, -0.5, 5.5);
    hTagRange_ = dir.make<TH1D>("hTagRange", "first and last event-window tag of the job (bin 1, 2)", 2, 0, 2);
    hBaseline_ = dir.make<TH1D>("hBaseline", ";Baseline field of the hit packet [ADC];hits", 4096, -0.5, 4095.5);
    hPedN_ = dir.make<TH2D>("hPedN", ";channel key;time [min]", nK, -0.5, nK - 0.5, nSlices_, 0, maxMinutes);
    hBaseDiff_ = dir.make<TH2D>("hBaseDiff", ";channel key;packet Baseline - mean of the first samples [ADC]",
                                nK, -0.5, nK - 0.5, 81, -40.5, 40.5);
    hBoardJump_ = dir.make<TH3D>("hBoardJump", ";board;time [min];peak sample - median of the laser event [samples]",
                                 kNBoards, -0.5, kNBoards - 0.5, nSlices_, 0, maxMinutes, 41, -20.5, 20.5);

    TH1D* hParams = dir.make<TH1D>("hParams", "windowMicroseconds, sliceMinutes, nSlices, laserMinHits, nPed, satADC",
                                   6, 0, 6);
    hParams->SetBinContent(1, windowMicroseconds_);
    hParams->SetBinContent(2, sliceMinutes_);
    hParams->SetBinContent(3, nSlices_);
    hParams->SetBinContent(4, laserMinHits_);
    hParams->SetBinContent(5, nPedSamples_);
    hParams->SetBinContent(6, satADC_);
  }

  void MC2RawCheck::analyze(const art::Event& event) {
    const CalDAQMap& daqMap = daqMapHandle_.get(event.id());
    if (!keyMapFilled_) {
      fillKeyToOfflineMap(daqMap);
    }

    const auto fragments = daqUtilities_.getFragments(event);
    const long long eventTag = eventWindowTag(fragments);
    if (eventTag < 0) {
      return;
    }
    fillTagHistograms(eventTag);
    // centre of the time slice of the event [min]
    const double sliceTime = (timeSlice(eventTag) + 0.5) * sliceMinutes_;

    const EventHits eventHits = checkRocBlocksAndHits(fragments, sliceTime);

    const double secondsFromStart = (eventTag - firstTag_) * windowMicroseconds_ * 1e-6;
    const bool isLaser = fillEventCounts(eventHits, sliceTime, secondsFromStart);

    // peak sample per channel key (the last good hit of a channel in the event)
    std::map<int, int> peakByKey;
    for (const auto& hit : eventHits.goodHits) {
      peakByKey[hit.key] = hit.peakSample;
    }
    if (isLaser && !peakByKey.empty()) {
      fillBoardJumps(peakByKey, sliceTime);
    }
    fillOccupancyAndPartners(peakByKey, isLaser, sliceTime, secondsFromStart);
  }

  // hKeyOffline: channel key -> offline SiPM id + 1, for the keys in the DAQ map (bin content 0 for the others).
  void MC2RawCheck::fillKeyToOfflineMap(const CalDAQMap& daqMap) {
    for (int key = 0; key < kNKeys; ++key) {
      CaloRawSiPMId rawId(key / kNChannelsPerBoard, key % kNChannelsPerBoard);
      const auto offlineId = daqMap.offlineId(rawId).id();
      if (offlineId < CaloConst::_nChannel) {
        hKeyOffline_->SetBinContent(key + 1, offlineId + 1);
      }
    }
    keyMapFilled_ = true;
  }

  // Offline SiPM id of a channel key from hKeyOffline (-1 if not in the map).
  int MC2RawCheck::offlineIdOfKey(int key) const {
    return int(hKeyOffline_->GetBinContent(key + 1)) - 1;
  }

  // Event-window tag of the event: the tag of the first calorimeter sub-event (-1 if there is none).
  long long MC2RawCheck::eventWindowTag(const artdaq::Fragments& fragments) const {
    long long eventTag = -1;
    for (const auto& fragment : fragments) {
      DTCEventFragment dtcEvent(fragment);
      for (auto& subEvent : dtcEvent.getSubsystemData(DTCLib::DTC_Subsystem::DTC_Subsystem_Calorimeter)) {
        eventTag = (long long)subEvent.GetEventWindowTag().GetEventWindowTag(true);
        break;
      }
      if (eventTag >= 0) {
        break;
      }
    }
    return eventTag;
  }

  // Tag step from the previous event (monotonicity) and the first / last tag of the job.
  void MC2RawCheck::fillTagHistograms(long long eventTag) {
    if (firstTag_ < 0) {
      firstTag_ = eventTag;
    }
    if (lastTag_ >= 0) {
      hTagStep_->Fill(std::clamp<long long>(eventTag - lastTag_, -20, 20));
    }
    lastTag_ = eventTag;

    if (hTagRange_->GetBinContent(1) == 0 || eventTag < hTagRange_->GetBinContent(1)) {
      hTagRange_->SetBinContent(1, double(eventTag));
    }
    if (eventTag > hTagRange_->GetBinContent(2)) {
      hTagRange_->SetBinContent(2, double(eventTag));
    }
  }

  // Time slice of a tag, from the tag of the first event and windowMicroseconds per window, in [0, nSlices - 1].
  int MC2RawCheck::timeSlice(long long tag) const {
    if (firstTag_ < 0) {
      return 0;
    }
    const double slice = (tag - firstTag_) * windowMicroseconds_ * 1e-6 / 60. / sliceMinutes_;
    return std::clamp(int(std::floor(slice)), 0, nSlices_ - 1);
  }

  // Loop over the calorimeter sub-events (one per DTC) and their ROC blocks (one per link): block presence and tag,
  // then every hit of the block.
  MC2RawCheck::EventHits MC2RawCheck::checkRocBlocksAndHits(const artdaq::Fragments& fragments, double sliceTime) {
    EventHits eventHits;
    for (const auto& fragment : fragments) {
      DTCEventFragment dtcEvent(fragment);
      for (auto& subEvent : dtcEvent.getSubsystemData(DTCLib::DTC_Subsystem::DTC_Subsystem_Calorimeter)) {
        const long long subEventTag = (long long)subEvent.GetEventWindowTag().GetEventWindowTag(true);
        CalorimeterDataDecoder decoder(subEvent);
        const int dtcId = subEvent.GetDTCID();

        for (size_t blockIndex = 0; blockIndex < decoder.block_count(); ++blockIndex) {
          auto block = decoder.dataAtBlockIndex(blockIndex);
          if (!block) {
            continue;
          }
          auto header = block->GetHeader();
          const int link = int(header->GetLinkID());
          ++eventHits.nBlocks;
          hBlocks_->Fill(dtcId, link);
          const long long blockTag = (long long)header->GetEventWindowTag().GetEventWindowTag(true);
          hBlockTagShift_->Fill(dtcId * kNLinks + link, std::clamp<long long>(blockTag - subEventTag, -4, 4));

          auto hits = decoder.GetCalorimeterHitData(blockIndex);
          if (!hits) {
            continue;
          }
          for (const auto& hit : *hits) {
            checkHit(hit, subEventTag, sliceTime, eventHits);
          }
        }
      }
    }
    return eventHits;
  }

  // Checks of one hit. Error codes and flags are counted for every hit; the per-channel histograms only for valid
  // channel keys; waveform, pedestal and the list of good hits only for the hits kept by isHitGood.
  void MC2RawCheck::checkHit(const CaloHitData& hit, long long subEventTag, double sliceTime, EventHits& eventHits) {
    ++eventHits.nHits;
    const auto& packet = hit.first;
    const auto& waveform = hit.second;
    const int key = int(packet.BoardID) * kNChannelsPerBoard + int(packet.ChannelID);
    const int errorCode = int(daqUtilities_.isHitGood(hit));
    const int flags = int(packet.ErrorFlags);

    hErrTot_->Fill(errorCode);
    hErrSlice_->Fill(sliceTime, errorCode);
    hFlagsTot_->Fill(flags);
    if (errorCode) {
      eventHits.anyBadHit = true;
    }
    for (int bit = 0; bit < 4; ++bit) {
      if (flags & (1 << bit)) {
        hFlagBitSlice_->Fill(sliceTime, bit);
      }
    }
    if (key < 0 || key >= kNKeys) {
      return;
    }

    hErr_->Fill(key, errorCode);
    hFlags_->Fill(key, flags);
    hDet_->Fill(key, int(packet.DetectorID));
    for (int bit = 0; bit < 4; ++bit) {
      if (flags & (1 << bit)) {
        (errorCode ? hFlagBitDropped_ : hFlagBitKept_)->Fill(key, bit);
      }
    }

    // in-payload tag (16 bits) minus the sub-event tag, wrapped to [-32768, 32768]
    int tagShift = int(packet.InPayloadEventWindowTag) - int(subEventTag & 0xFFFF);
    if (tagShift > 32768) {
      tagShift -= 65536;
    }
    if (tagShift < -32768) {
      tagShift += 65536;
    }
    hHitTagShift_->Fill(key, std::clamp(tagShift, -4, 4));

    if (errorCode) {
      return;
    }

    hNSamp_->Fill(key, waveform.size());
    hPeakIdx_->Fill(packet.IndexOfMaxDigitizerSample);
    int maxSample = 0;
    for (auto sample : waveform) {
      maxSample = std::max<int>(maxSample, sample);
    }
    if (maxSample >= satADC_) {
      hSat_->Fill(key);
    }
    if ((int)waveform.size() > nPedSamples_) {
      fillPedestal(key, packet, waveform, sliceTime);
    }
    eventHits.goodHits.push_back({key, int(packet.Time) + int(packet.IndexOfMaxDigitizerSample)});
  }

  // Pedestal = mean of the first nPed samples, its variance, and the Baseline field of the packet.
  void MC2RawCheck::fillPedestal(int key, const CalorimeterDataDecoder::CalorimeterHitDataPacket& packet,
                                 const std::vector<uint16_t>& waveform, double sliceTime) {
    double sum = 0;
    double sumOfSquares = 0;
    for (int i = 0; i < nPedSamples_; ++i) {
      sum += waveform[i];
      sumOfSquares += double(waveform[i]) * waveform[i];
    }
    const double mean = sum / nPedSamples_;
    const double variance = std::max(0., sumOfSquares / nPedSamples_ - mean * mean);

    hBaseDiff_->Fill(key, double(packet.Baseline) - mean);
    hBaseline_->Fill(int(packet.Baseline));
    hPedSum_->Fill(key, sliceTime, mean);
    hPed2Sum_->Fill(key, sliceTime, mean * mean);
    hNoiseSum_->Fill(key, sliceTime, variance);
    hPedN_->Fill(key, sliceTime);
  }

  // Per-event counts: events, dropped hits, hits, blocks, events with bad hits, laser events. Returns true for a laser
  // event (at least laserMinHits hits, good or bad).
  bool MC2RawCheck::fillEventCounts(const EventHits& eventHits, double sliceTime, double secondsFromStart) {
    hRateFine_->Fill(secondsFromStart, 0);
    hRateFine_->Fill(secondsFromStart, 5, double(eventHits.nHits - (int)eventHits.goodHits.size()));
    hEvSlice_->Fill(sliceTime);
    hNHits_->Fill(eventHits.nHits);
    hNBlocks_->Fill(eventHits.nBlocks);
    if (eventHits.anyBadHit) {
      hEvFailSlice_->Fill(sliceTime);
    }

    const bool isLaser = eventHits.nHits >= laserMinHits_;
    if (isLaser) {
      hLaserSlice_->Fill(sliceTime);
      hRateFine_->Fill(secondsFromStart, 1);
    }
    return isLaser;
  }

  // Laser events: for every board, the median over its channels of (peak sample - median peak sample of the event).
  // A board that jumps by whole samples shows up as an entry away from 0.
  void MC2RawCheck::fillBoardJumps(const std::map<int, int>& peakByKey, double sliceTime) {
    std::vector<int> peaks;
    for (const auto& [key, peak] : peakByKey) {
      peaks.push_back(peak);
    }
    std::nth_element(peaks.begin(), peaks.begin() + peaks.size() / 2, peaks.end());
    const int eventMedian = peaks[peaks.size() / 2];

    std::map<int, std::vector<int>> shiftsByBoard;
    for (const auto& [key, peak] : peakByKey) {
      shiftsByBoard[key / kNChannelsPerBoard].push_back(peak - eventMedian);
    }
    for (auto& [board, shifts] : shiftsByBoard) {
      std::nth_element(shifts.begin(), shifts.begin() + shifts.size() / 2, shifts.end());
      hBoardJump_->Fill(board, sliceTime, std::clamp(shifts[shifts.size() / 2], -20, 20));
    }
  }

  // Occupancy per channel (laser and non-laser events), and for non-laser events: occupancy per slice, good-hit rates
  // per category (LYSO, disk 0, disk 1) and the pairs of SiPMs of the same crystal.
  void MC2RawCheck::fillOccupancyAndPartners(const std::map<int, int>& peakByKey, bool isLaser, double sliceTime,
                                             double secondsFromStart) {
    // peak sample per offline SiPM id
    std::map<int, int> peakByOfflineId;
    for (const auto& [key, peak] : peakByKey) {
      const int offlineId = offlineIdOfKey(key);
      if (offlineId >= 0) {
        peakByOfflineId[offlineId] = peak;
      }
    }

    for (const auto& [key, peak] : peakByKey) {
      (isLaser ? hOccLaser_ : hOcc_)->Fill(key);
      if (isLaser) {
        continue;
      }
      hOccSlice_->Fill(key, sliceTime);

      // rate category: LYSO (a channel that had hits with DetectorID 1, so far in the job), otherwise the disk
      // (offline crystal id >= 674 on disk 1)
      const int offlineId = offlineIdOfKey(key);
      const bool isLyso = hDet_->GetBinContent(key + 1, 2) > 0;
      if (isLyso) {
        hRateFine_->Fill(secondsFromStart, 4);
      } else if (offlineId >= 0) {
        hRateFine_->Fill(secondsFromStart, (offlineId / 2) >= 674 ? 3 : 2);
      }

      // the other SiPM of the crystal: offline ids 2c and 2c + 1
      if (offlineId < 0) {
        continue;
      }
      const auto partner = peakByOfflineId.find(offlineId ^ 1);
      if (partner != peakByOfflineId.end()) {
        hPartner_->Fill(key);
        hPartnerDt_->Fill(key, std::clamp(peak - partner->second, -300, 300));
      }
    }
  }

}  // namespace mu2e

DEFINE_ART_MODULE(mu2e::MC2RawCheck)
