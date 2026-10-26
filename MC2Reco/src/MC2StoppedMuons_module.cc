//
// MC2StoppedMuons: cosmic muons that stop in the calorimeter and their delayed deposits (decay electrons of the mu+,
// Michel spectrum), on calibrated CaloHits.
//
// It is a port of the mc2tools macro cosmic_delayed_coincidence.C (same selection, same tree), so that
// michel_scale.py reads its output as it is; in addition it fills the Michel spectra per disk and fits the lifetime
// at the end of the job.
//
// Input: CaloHitCollection (caloHitCollection, calibrated hits) and the CalCalib conditions (ADC2MeV).
//
// Definitions:
//   Crystal  : CaloHit with both SiPMs (nSiPMs == 2) and the amplitude of each SiPM >= minAmp ADC counts (energy /
//              CalCalib ADC2MeV); laser PIN diodes (id >= number of crystals) and the crystals of excludedCrystals
//              are left out;
//                time      = hit time (our template, after the time offsets of the conditions);
//                amplitude = mean of the two SiPM amplitudes;
//                energy    = CaloHit energy (MeV, mean of the two SiPMs, from the energy table of the conditions);
//                position  = crystal x, y (tracker frame of the Offline geometry).
//              More than one hit per crystal per event is allowed.
//   Cluster  : crystals in time order; a new cluster starts after a gap > clusGap ns.
//   Pair     : prompt cluster P and later cluster D of the same event window, dt = t_D - t_P in [dtMin, dtMax] us
//              (times of the first crystal of each cluster). P must start more than tEdge us after the window start and
//              more than dtMax before its end (window length windowNs), so that the acceptance in dt is full.
//                near (cls 1) = some crystal of D within dNear mm (same disk) of some crystal of P;
//                far  (cls 2) = all crystals of D farther than dFar mm (accidental control);
//              pairs in between are dropped.
//   Michel   : per disk of the delayed cluster, energy eD of the near pairs with dt in [michelDtMin, dtMax] us
//              (mu+ only: the mu- are captured in < 0.1 us in CsI), and of the far pairs; the accidentals are
//              subtracted with k = near/far in dt > dtNorm us (where the delayed rate is flat).
//   Lifetime : binned Poisson likelihood fit N0 exp(-t/tau) + B of the near pairs, dt in [fitMin, fitMax] us, all disks
//              and per disk.
//
// Output:
//   TTree "pairs", one entry per near or far pair (same branches as mc2tools pairs.root):
//     entry          art event number (not the tree entry)
//     dt             delay t_D - t_P [us]
//     dmin           minimum distance between a crystal of D and a crystal of P on the same disk [mm] (1e9 if none)
//     tP             time of the prompt cluster [us]
//     aP, aD, aDmax  sum of the crystal amplitudes of P, of D, largest crystal amplitude of D [ADC]
//     eP, eD, eDmax  same for the energies [MeV]
//     nP, nD         number of crystals of P and of D
//     disk           disk of the first crystal of D
//     diskP          disk of the crystals of P (2 if on both disks)
//     sameCry        1 if a crystal of D is also a crystal of P
//     cls            1 near, 2 far
//     nCl            number of clusters in the event
//     xD, yD         amplitude-weighted centroid of D [mm]
//     xP, yP         amplitude-weighted centroid of P [mm]
//     xLo, yLo       crystal of P with the lowest y; xHi, yHi the one with the highest y [mm]
//   Histograms:
//     hDtNear_<all|disk0|disk1>, hDtFar_<...>  delay of the near / far pairs (all, per disk of D)
//     hMichelNear_disk<d>, hMichelFar_disk<d>  eD of the near / far pairs with dt >= michelDtMin
//     hMichel_disk<d>                          near - k far (filled at the end of the job)
//     hNorm_disk<d>                            near (bin 0) and far (bin 1) pairs with dt >= dtNorm, for k
//   Text file summaryFile (also printed): counters, cuts, lifetime fits, Michel spectra after the subtraction.
//
// Steps per event: select the crystals, sort them in time, build the clusters, then for every prompt cluster in
// acceptance look at the later clusters, classify each pair (near / far), fill the tree and the histograms.
//
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art_root_io/TFileService.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"

#include "Offline/CaloConditions/inc/CalCalib.hh"
#include "Offline/CalorimeterGeom/inc/Calorimeter.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/ProditionsService/inc/ProditionsHandle.hh"
#include "Offline/RecoDataProducts/inc/CaloHit.hh"
#include "Offline/RecoDataProducts/inc/CaloRecoDigi.hh"
#include "TF1.h"
#include "TH1F.h"
#include "TTree.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <map>
#include <set>
#include <vector>

namespace mu2e {

  class MC2StoppedMuons : public art::EDAnalyzer {
    public:
      struct Config {
        using Name = fhicl::Name;
        using Comment = fhicl::Comment;

        // Input
        fhicl::Atom<art::InputTag> hitTag{Name("caloHitCollection"), Comment("calibrated calorimeter hits"),
            "MC2CaloHitMaker"};

        // Crystals and clusters
        fhicl::Atom<double> minAmp{Name("minAmp"), Comment("min amplitude of each SiPM [ADC]"), 100.};
        fhicl::Atom<double> clusGap{Name("clusGap"), Comment("time gap that starts a new cluster [ns]"), 30.};

        // Pairs
        fhicl::Atom<double> dtMin{Name("dtMin"), Comment("min delay [us]"), 0.1};
        fhicl::Atom<double> dtMax{Name("dtMax"), Comment("max delay [us]"), 20.};
        fhicl::Atom<double> tEdge{Name("tEdge"), Comment("min prompt time from the window start [us]"), 1.};
        fhicl::Atom<double> windowNs{Name("windowNs"), Comment("event window length [ns]"), 97495.};
        fhicl::Atom<double> dNear{Name("dNear"), Comment("near: max distance prompt-delayed crystal [mm]"), 80.};
        fhicl::Atom<double> dFar{Name("dFar"), Comment("far: min distance [mm]"), 300.};

        // Michel spectra and accidental subtraction
        fhicl::Atom<double> michelDtMin{Name("michelDtMin"), Comment("min delay for the Michel spectra [us]"), 1.};
        fhicl::Atom<double> dtNorm{Name("dtNorm"), Comment("near/far normalisation for dt > dtNorm [us]"), 15.};

        // Lifetime fit and summary
        fhicl::Atom<double> fitMin{Name("fitMin"), Comment("lifetime fit from [us]"), 1.};
        fhicl::Atom<double> fitMax{Name("fitMax"), Comment("lifetime fit to [us]"), 20.};
        fhicl::Atom<std::string> summaryFile{Name("summaryFile"), Comment("text summary"),
            "stopped_muons_summary.txt"};

        // Crystals left out of the analysis (e.g. the LYSO at operating voltage, whose intrinsic activity gives a
        // continuous rate of delayed deposits); empty: all crystals
        fhicl::Sequence<int> excludedCrystals{Name("excludedCrystals"), Comment("crystal ids left out"),
            std::vector<int>{}};
      };
      using Parameters = art::EDAnalyzer::Table<Config>;

      explicit MC2StoppedMuons(const Parameters& config);
      void beginJob() override;
      void analyze(const art::Event& event) override;
      void endJob() override;

    private:
      // A selected crystal (one CaloHit)
      struct CrystalHit {
        int disk;
        float x, y;          // tracker frame [mm]
        double time;         // [ns]
        double amplitude;    // mean of the two SiPMs [ADC]
        double energy;       // [MeV]
        int crystalId;
      };

      // A cluster: crystals [begin, end) of the time-ordered crystal list
      struct Cluster {
        int begin;
        int end;
      };

      // Sums over the prompt cluster
      struct PromptSums {
        double amplitude = 0;
        double energy = 0;
      };

      // Sums over the delayed cluster and its distance to the prompt cluster
      struct DelayedSums {
        double minDistance = 1e9;
        double amplitude = 0;
        double amplitudeMax = 0;
        double energy = 0;
        double energyMax = 0;
        int sameCrystal = 0;
      };

      // Branches of the "pairs" tree (names as the branch names)
      struct PairBranches {
        long long entry;
        float dt, dmin, tP, aP, aD, aDmax, xD, yD, xLo, yLo, xHi, yHi, xP, yP, eP, eD, eDmax;
        int nP, nD, disk, sameCry, cls, nCl, diskP;
      };

      static constexpr int kNear = 0;   // second index of the delay and Michel histograms
      static constexpr int kFar = 1;

      std::vector<CrystalHit> selectCrystals(const art::Event& event);
      bool meanSiPMAmplitude(const CaloHit& hit, const CalCalib& calib, double& meanAmplitude) const;
      std::vector<Cluster> buildClusters(const std::vector<CrystalHit>& crystals) const;
      bool promptInAcceptance(double promptTime) const;
      PromptSums sumPrompt(const std::vector<CrystalHit>& crystals, const Cluster& prompt) const;
      DelayedSums sumDelayed(const std::vector<CrystalHit>& crystals, const Cluster& delayed,
                             const Cluster& prompt) const;
      int pairClass(double minDistance) const;
      void fillPairTree(const art::Event& event, const std::vector<CrystalHit>& crystals,
                        const std::vector<Cluster>& clusters, size_t promptIndex, size_t delayedIndex,
                        double dt, double promptTime, const PromptSums& prompt, const DelayedSums& delayed,
                        int disk, int cls);
      void fillHistograms(int cls, int disk, double dt, double delayedEnergy);
      void fitLifetimes();
      void subtractAccidentals();
      template <typename... Args> void report(const char* format, Args... args);

      // Input
      art::InputTag hitTag_;
      ProditionsHandle<CalCalib> calibHandle_;

      // Parameters
      double minAmp_;
      double clusterGap_;
      double dtMin_;
      double dtMax_;
      double tEdge_;
      double windowNs_;
      double dNear_;
      double dFar_;
      double michelDtMin_;
      double dtNorm_;
      double fitMin_;
      double fitMax_;
      std::string summaryFileName_;
      std::set<int> excludedCrystals_;

      // Output
      TTree* pairTree_ = nullptr;
      PairBranches pair_;
      TH1F* hDelay_[3][2] = {{nullptr}};               // [all, disk 0, disk 1][near, far]
      TH1F* hMichel_[2][2] = {{nullptr}};              // [disk][near, far]
      TH1F* hMichelSubtracted_[2] = {nullptr, nullptr};   // [disk]: near - k far
      TH1F* hNormalisation_[2] = {nullptr, nullptr};      // [disk]: near, far pairs with dt >= dtNorm
      FILE* summaryFile_ = nullptr;

      // Counters
      long nEvents_ = 0;
      long nCrystals_ = 0;
      long nClusters_ = 0;
      long nPromptInAcceptance_ = 0;
  };

  MC2StoppedMuons::MC2StoppedMuons(const Parameters& config) :
    art::EDAnalyzer(config),
    hitTag_(config().hitTag()),
    minAmp_(config().minAmp()),
    clusterGap_(config().clusGap()),
    dtMin_(config().dtMin()),
    dtMax_(config().dtMax()),
    tEdge_(config().tEdge()),
    windowNs_(config().windowNs()),
    dNear_(config().dNear()),
    dFar_(config().dFar()),
    michelDtMin_(config().michelDtMin()),
    dtNorm_(config().dtNorm()),
    fitMin_(config().fitMin()),
    fitMax_(config().fitMax()),
    summaryFileName_(config().summaryFile())
  {
    const auto excluded = config().excludedCrystals();   // returned by value: one copy for begin and end
    excludedCrystals_.insert(excluded.begin(), excluded.end());
    consumes<CaloHitCollection>(hitTag_);
  }

  void MC2StoppedMuons::beginJob() {
    art::ServiceHandle<art::TFileService> tfs;

    pairTree_ = tfs->make<TTree>("pairs", "delayed coincidences (branches of mc2tools cosmic_delayed_coincidence.C)");
    pairTree_->Branch("entry", &pair_.entry);
    pairTree_->Branch("dt", &pair_.dt);
    pairTree_->Branch("dmin", &pair_.dmin);
    pairTree_->Branch("tP", &pair_.tP);
    pairTree_->Branch("aP", &pair_.aP);
    pairTree_->Branch("aD", &pair_.aD);
    pairTree_->Branch("aDmax", &pair_.aDmax);
    pairTree_->Branch("nP", &pair_.nP);
    pairTree_->Branch("nD", &pair_.nD);
    pairTree_->Branch("disk", &pair_.disk);
    pairTree_->Branch("sameCry", &pair_.sameCry);
    pairTree_->Branch("cls", &pair_.cls);
    pairTree_->Branch("nCl", &pair_.nCl);
    pairTree_->Branch("xD", &pair_.xD);
    pairTree_->Branch("yD", &pair_.yD);
    pairTree_->Branch("xLo", &pair_.xLo);
    pairTree_->Branch("yLo", &pair_.yLo);
    pairTree_->Branch("xHi", &pair_.xHi);
    pairTree_->Branch("yHi", &pair_.yHi);
    pairTree_->Branch("xP", &pair_.xP);
    pairTree_->Branch("yP", &pair_.yP);
    pairTree_->Branch("diskP", &pair_.diskP);
    pairTree_->Branch("eP", &pair_.eP);
    pairTree_->Branch("eD", &pair_.eD);
    pairTree_->Branch("eDmax", &pair_.eDmax);

    const char* selectionName[3] = {"all", "disk0", "disk1"};
    for (int k = 0; k < 3; ++k) {
      hDelay_[k][kNear] = tfs->make<TH1F>(Form("hDtNear_%s", selectionName[k]),
          Form("near pairs, %s;delay [#mus];pairs / 0.2 #mus", selectionName[k]), 100, 0, 20);
      hDelay_[k][kFar] = tfs->make<TH1F>(Form("hDtFar_%s", selectionName[k]),
          Form("far pairs, %s;delay [#mus];pairs / 0.2 #mus", selectionName[k]), 100, 0, 20);
    }
    for (int disk = 0; disk < 2; ++disk) {
      hMichel_[disk][kNear] = tfs->make<TH1F>(Form("hMichelNear_disk%d", disk),
          Form("disk %d, near pairs, delay > %.1f #mus;delayed energy [MeV];pairs / 2 MeV", disk, michelDtMin_),
          50, 0, 100);
      hMichel_[disk][kFar] = tfs->make<TH1F>(Form("hMichelFar_disk%d", disk),
          Form("disk %d, far pairs, delay > %.1f #mus;delayed energy [MeV];pairs / 2 MeV", disk, michelDtMin_),
          50, 0, 100);
      hMichelSubtracted_[disk] = tfs->make<TH1F>(Form("hMichel_disk%d", disk),
          Form("disk %d, near - k far (accidentals subtracted);delayed energy [MeV];pairs / 2 MeV", disk),
          50, 0, 100);
      hNormalisation_[disk] = tfs->make<TH1F>(Form("hNorm_disk%d", disk),
          Form("disk %d: near and far pairs with delay > %.0f #mus;0 near, 1 far;pairs", disk, dtNorm_),
          2, -0.5, 1.5);
    }
  }

  void MC2StoppedMuons::analyze(const art::Event& event) {
    ++nEvents_;

    // 1. crystals, in time order
    std::vector<CrystalHit> crystals = selectCrystals(event);
    nCrystals_ += crystals.size();
    if (crystals.size() < 2) {
      return;
    }
    std::sort(crystals.begin(), crystals.end(),
              [](const CrystalHit& first, const CrystalHit& second) { return first.time < second.time; });

    // 2. clusters
    const std::vector<Cluster> clusters = buildClusters(crystals);
    nClusters_ += clusters.size();
    if (clusters.size() < 2) {
      return;
    }

    // 3. pairs: every prompt cluster in acceptance with every later cluster in the delay range
    for (size_t p = 0; p < clusters.size(); ++p) {
      const double promptTime = crystals[clusters[p].begin].time;
      if (!promptInAcceptance(promptTime)) {
        continue;
      }
      ++nPromptInAcceptance_;
      const PromptSums prompt = sumPrompt(crystals, clusters[p]);

      for (size_t d = p + 1; d < clusters.size(); ++d) {
        const double dt = (crystals[clusters[d].begin].time - promptTime) * 1e-3;
        if (dt < dtMin_) {
          continue;
        }
        if (dt > dtMax_) {
          break;
        }
        const int disk = crystals[clusters[d].begin].disk;
        const DelayedSums delayed = sumDelayed(crystals, clusters[d], clusters[p]);
        const int cls = pairClass(delayed.minDistance);
        if (!cls) {
          continue;
        }
        fillPairTree(event, crystals, clusters, p, d, dt, promptTime, prompt, delayed, disk, cls);
        fillHistograms(cls, disk, dt, delayed.energy);
      }
    }
  }

  // Crystals of the event: hits with both SiPMs above minAmp.
  std::vector<MC2StoppedMuons::CrystalHit> MC2StoppedMuons::selectCrystals(const art::Event& event) {
    const Calorimeter& calorimeter = *GeomHandle<Calorimeter>();
    const CalCalib& calib = calibHandle_.get(event.id());
    const auto& hits = *event.getValidHandle<CaloHitCollection>(hitTag_);

    std::vector<CrystalHit> crystals;
    for (const auto& hit : hits) {
      if (hit.nSiPMs() < 2) {
        continue;
      }
      if (hit.crystalID() >= (int)calorimeter.nCrystals()) {   // laser PIN diodes: not crystals
        continue;
      }
      if (excludedCrystals_.count(hit.crystalID())) {
        continue;
      }
      double meanAmplitude;
      if (!meanSiPMAmplitude(hit, calib, meanAmplitude)) {
        continue;
      }
      const auto& crystal = calorimeter.crystal(hit.crystalID());
      const auto position = calorimeter.mu2eToTracker(crystal.position());
      crystals.push_back({int(crystal.diskID()), float(position.x()), float(position.y()), hit.time(), meanAmplitude,
                          hit.energyDep(), hit.crystalID()});
    }
    return crystals;
  }

  // Mean amplitude [ADC] of the SiPMs of the hit (energy / ADC2MeV, 0 if ADC2MeV <= 0). False if one SiPM is below
  // minAmp or the hit has no SiPM.
  bool MC2StoppedMuons::meanSiPMAmplitude(const CaloHit& hit, const CalCalib& calib, double& meanAmplitude) const {
    double amplitudeSum = 0;
    int nSiPMs = 0;
    for (const auto& recoDigi : hit.recoCaloDigis()) {
      const double adc2MeV = calib.ADC2MeV(recoDigi->SiPMID());
      const double amplitude = adc2MeV > 0 ? recoDigi->energyDep() / adc2MeV : 0;
      if (amplitude < minAmp_) {
        return false;
      }
      amplitudeSum += amplitude;
      ++nSiPMs;
    }
    if (nSiPMs == 0) {
      return false;
    }
    meanAmplitude = amplitudeSum / nSiPMs;
    return true;
  }

  // Clusters of the time-ordered crystals: a new cluster after a gap > clusGap.
  std::vector<MC2StoppedMuons::Cluster> MC2StoppedMuons::buildClusters(const std::vector<CrystalHit>& crystals) const {
    std::vector<Cluster> clusters;
    int clusterBegin = 0;
    for (size_t i = 1; i <= crystals.size(); ++i) {
      if (i == crystals.size() || crystals[i].time - crystals[i - 1].time > clusterGap_) {
        clusters.push_back({clusterBegin, (int)i});
        clusterBegin = i;
      }
    }
    return clusters;
  }

  // The prompt cluster must start after tEdge and at least dtMax before the end of the window [ns].
  bool MC2StoppedMuons::promptInAcceptance(double promptTime) const {
    return !(promptTime < tEdge_ * 1e3 || promptTime > windowNs_ - dtMax_ * 1e3);
  }

  MC2StoppedMuons::PromptSums MC2StoppedMuons::sumPrompt(const std::vector<CrystalHit>& crystals,
                                                         const Cluster& prompt) const {
    PromptSums sums;
    for (int i = prompt.begin; i < prompt.end; ++i) {
      sums.amplitude += crystals[i].amplitude;
      sums.energy += crystals[i].energy;
    }
    return sums;
  }

  // Sums over the delayed cluster, and its minimum distance to the prompt cluster (crystals on the same disk only).
  MC2StoppedMuons::DelayedSums MC2StoppedMuons::sumDelayed(const std::vector<CrystalHit>& crystals,
                                                           const Cluster& delayed, const Cluster& prompt) const {
    DelayedSums sums;
    for (int i = delayed.begin; i < delayed.end; ++i) {
      const CrystalHit& delayedCrystal = crystals[i];
      sums.amplitude += delayedCrystal.amplitude;
      sums.amplitudeMax = std::max(sums.amplitudeMax, delayedCrystal.amplitude);
      sums.energy += delayedCrystal.energy;
      sums.energyMax = std::max(sums.energyMax, delayedCrystal.energy);
      for (int j = prompt.begin; j < prompt.end; ++j) {
        const CrystalHit& promptCrystal = crystals[j];
        if (delayedCrystal.disk != promptCrystal.disk) {
          continue;
        }
        // hypot of two floats (float), as in the mc2tools macro
        sums.minDistance = std::min(sums.minDistance, double(std::hypot(delayedCrystal.x - promptCrystal.x,
                                                                        delayedCrystal.y - promptCrystal.y)));
        if (delayedCrystal.crystalId == promptCrystal.crystalId) {
          sums.sameCrystal = 1;
        }
      }
    }
    return sums;
  }

  // 1 near, 2 far, 0 neither (pair dropped).
  int MC2StoppedMuons::pairClass(double minDistance) const {
    return minDistance <= dNear_ ? 1 : (minDistance > dFar_ ? 2 : 0);
  }

  void MC2StoppedMuons::fillPairTree(const art::Event& event, const std::vector<CrystalHit>& crystals,
                                     const std::vector<Cluster>& clusters, size_t promptIndex, size_t delayedIndex,
                                     double dt, double promptTime, const PromptSums& prompt,
                                     const DelayedSums& delayed, int disk, int cls) {
    const Cluster& promptCluster = clusters[promptIndex];
    const Cluster& delayedCluster = clusters[delayedIndex];

    // amplitude-weighted centroid of the delayed cluster
    double sumAX = 0, sumAY = 0, sumA = 0;
    for (int i = delayedCluster.begin; i < delayedCluster.end; ++i) {
      sumAX += crystals[i].amplitude * crystals[i].x;
      sumAY += crystals[i].amplitude * crystals[i].y;
      sumA += crystals[i].amplitude;
    }

    // prompt cluster: lowest and highest crystal in y, amplitude-weighted centroid, disks
    int lowest = promptCluster.begin;
    int highest = promptCluster.begin;
    double promptSumAX = 0, promptSumAY = 0, promptSumA = 0;
    std::map<int, int> nCrystalsPerDisk;
    for (int i = promptCluster.begin; i < promptCluster.end; ++i) {
      if (crystals[i].y < crystals[lowest].y) {
        lowest = i;
      }
      if (crystals[i].y > crystals[highest].y) {
        highest = i;
      }
      promptSumAX += crystals[i].amplitude * crystals[i].x;
      promptSumAY += crystals[i].amplitude * crystals[i].y;
      promptSumA += crystals[i].amplitude;
      nCrystalsPerDisk[crystals[i].disk]++;
    }

    pair_.entry = event.event();
    pair_.dt = dt;
    pair_.dmin = delayed.minDistance;
    pair_.tP = promptTime * 1e-3;
    pair_.aP = prompt.amplitude;
    pair_.aD = delayed.amplitude;
    pair_.aDmax = delayed.amplitudeMax;
    pair_.nP = promptCluster.end - promptCluster.begin;
    pair_.nD = delayedCluster.end - delayedCluster.begin;
    pair_.disk = disk;
    pair_.sameCry = delayed.sameCrystal;
    pair_.cls = cls;
    pair_.nCl = clusters.size();
    pair_.xD = sumAX / sumA;
    pair_.yD = sumAY / sumA;
    pair_.xLo = crystals[lowest].x;
    pair_.yLo = crystals[lowest].y;
    pair_.xHi = crystals[highest].x;
    pair_.yHi = crystals[highest].y;
    pair_.xP = promptSumAX / promptSumA;
    pair_.yP = promptSumAY / promptSumA;
    pair_.diskP = nCrystalsPerDisk.size() > 1 ? 2 : nCrystalsPerDisk.begin()->first;
    pair_.eP = prompt.energy;
    pair_.eD = delayed.energy;
    pair_.eDmax = delayed.energyMax;
    pairTree_->Fill();
  }

  // Delay of the pair (all and per disk of the delayed cluster), Michel energy, near/far counts for the normalisation.
  void MC2StoppedMuons::fillHistograms(int cls, int disk, double dt, double delayedEnergy) {
    const int nearOrFar = cls == 1 ? kNear : kFar;
    const bool knownDisk = disk == 0 || disk == 1;
    hDelay_[0][nearOrFar]->Fill(dt);
    if (knownDisk) {
      hDelay_[1 + disk][nearOrFar]->Fill(dt);
    }
    if (knownDisk && dt >= michelDtMin_) {
      hMichel_[disk][nearOrFar]->Fill(delayedEnergy);
    }
    if (knownDisk && dt >= dtNorm_) {
      hNormalisation_[disk]->Fill(nearOrFar);
    }
  }

  // printf to the terminal and to the summary file.
  template <typename... Args> void MC2StoppedMuons::report(const char* format, Args... args) {
    printf(format, args...);
    if (summaryFile_) {
      fprintf(summaryFile_, format, args...);
    }
  }

  void MC2StoppedMuons::endJob() {
    summaryFile_ = fopen(summaryFileName_.c_str(), "w");
    report("MC2StoppedMuons: events %ld, crystals %ld, clusters %ld, prompt clusters in acceptance %ld\n",
           nEvents_, nCrystals_, nClusters_, nPromptInAcceptance_);
    report("cuts: minAmp %.0f ADC, clusGap %.0f ns, dt %.2f-%.1f us, tEdge %.1f us, near <= %.0f mm, far > %.0f mm\n",
           minAmp_, clusterGap_, dtMin_, dtMax_, tEdge_, dNear_, dFar_);
    fitLifetimes();
    subtractAccidentals();
    if (summaryFile_) {
      fclose(summaryFile_);
      summaryFile_ = nullptr;
    }
  }

  // Lifetime: binned likelihood fit N0 exp(-t/tau) + B of the near-pair delay, all disks and per disk.
  void MC2StoppedMuons::fitLifetimes() {
    const char* selectionName[3] = {"all", "disk 0", "disk 1"};
    for (int k = 0; k < 3; ++k) {
      TH1F* hNear = hDelay_[k][kNear];
      if (hNear->Integral() < 20) {
        report("lifetime %s: too few near pairs (%.0f)\n", selectionName[k], hNear->Integral());
        continue;
      }
      TF1 decay(Form("fexp%d", k), "[0]*exp(-x/[1])+[2]", fitMin_, fitMax_);
      decay.SetParameters(std::max(1., hNear->GetBinContent(hNear->FindBin(fitMin_)) * 1.5), 2.2,
                          std::max(0.1, hNear->GetBinContent(hNear->GetNbinsX())));
      decay.SetParLimits(1, 0.1, 20.);
      hNear->Fit(&decay, "LQR0");
      report("lifetime %s: tau = %.3f +- %.3f us, N0 = %.1f +- %.1f per bin, B = %.2f +- %.2f per bin, "
             "chi2/ndf = %.1f/%d, near pairs %.0f, far pairs %.0f\n",
             selectionName[k], decay.GetParameter(1), decay.GetParError(1), decay.GetParameter(0),
             decay.GetParError(0), decay.GetParameter(2), decay.GetParError(2), decay.GetChisquare(), decay.GetNDF(),
             hNear->Integral(), hDelay_[k][kFar]->Integral());
    }
  }

  // Michel spectra per disk: near - k far, with k = near/far for dt >= dtNorm (only accidentals there).
  void MC2StoppedMuons::subtractAccidentals() {
    for (int disk = 0; disk < 2; ++disk) {
      const double nNear = hNormalisation_[disk]->GetBinContent(1);
      const double nFar = hNormalisation_[disk]->GetBinContent(2);
      const double k = nFar > 0 ? nNear / nFar : 0;
      hMichelSubtracted_[disk]->Add(hMichel_[disk][kNear], hMichel_[disk][kFar], 1., -k);
      report("Michel disk %d: near %.0f, far %.0f (delay > %.1f us), accidental scale k = %.3f "
             "(near/far for delay > %.0f us: %.0f/%.0f); after subtraction %.0f pairs, mean %.1f MeV\n",
             disk, hMichel_[disk][kNear]->Integral(), hMichel_[disk][kFar]->Integral(), michelDtMin_, k, dtNorm_,
             nNear, nFar, hMichelSubtracted_[disk]->Integral(), hMichelSubtracted_[disk]->GetMean());
    }
  }

}

DEFINE_ART_MODULE(mu2e::MC2StoppedMuons)
