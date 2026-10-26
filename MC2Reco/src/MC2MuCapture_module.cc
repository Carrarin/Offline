//
// MC2MuCapture: search for the nuclear capture of stopped mu- in CsI, i.e. a second pulse in the crystal where a cosmic
// muon stops. It is a detector, not a time fit: the first pulse is removed and what is left is searched for "something
// with the shape of a pulse".
//
// Inputs (products):
//   caloTrackCollection     MC2CaloCosmicTracks (MC2CaloCosmicTrackFinder): the candidates
//   caloDigiCollection      CaloDigis: waveforms, and which crystals have a digi at all
//   caloHitCollection       CaloHits of the event: mu+ veto
//   caloRecoDigiCollection  CaloRecoDigis of the event: delayed deposits in the other crystals, delayed clusters, LYSO
//
// Samples (as mc2tools capture_double_pulse.C), from the calorimeter tracks of the job:
//   STOP    (kind 0): one-disk tracks with >= minCrystals crystals; the last crystal along the downward direction of
//                     the track in the disk plane, if the crystal one pitch (34.3 mm) further along that direction
//                     exists and has no CaloDigi on either SiPM;
//   THROUGH (kind 1): tracks with >= minThrough crystals, the crystal in the middle (the muon went through: no capture
//                     there).
// Waveforms: the CaloDigis of the two SiPMs of the crystal (through the CaloHit of the track); pedestal = mean of the
//   first 4 samples; rejected if shorter than minSamples, saturated (raw >= satADC) or with amplitude < minAmp ADC.
//
// Steps:
//   analyze (per event): find the STOP and THROUGH candidates; for each one store its two waveforms and everything that
//     needs the event (mu+ veto, delayed deposits, cluster deposits, tags), and write the tree "delayed" right away.
//   endJob (needs the whole run, for the template):
//     1. template from the THROUGH waveforms of the run (amplitude 300-1800 ADC), aligned at the 50 % crossing,
//        iterated nIterTemplate times (as mc2tools); 0.5 ns bins over [-60, xMax] ns. Optionally a stop template;
//     2. for every candidate: second-pulse search (below), tree "capture", delayed clusters (tree "cluster"), example
//        waveforms (tree "examples"), and the same search with an injected pulse (tree "inject").
//
// Second-pulse search, per crystal:
//   1. first pulse of each SiPM: amplitude and time from a template fit on [-20, +40] ns around the 50 % crossing
//      (leading edge and peak);
//   2. residual r = y - A1 T(t - t1) of each SiPM, normalised to its A1, the two SiPMs averaged on the common sample
//      grid (t0 of the CaloDigis): the light of a capture reaches both SiPMs, the electronic noise does not;
//   3. matched filter S(tau) = sum r T(t - t1 - tau) / (sigma_r sqrt(sum T^2)), tau scanned from dtMin to the end of
//      the waveform - 30 ns in 2.5 ns steps; maxS and its tau (delay of the second pulse), A2 = amplitude of the second
//      pulse at tau [ADC, mean of the SiPMs]; sigma_r from the noise per sample (noiseADC, white noise assumed in this
//      version: maxS is a ranking, the threshold comes from the THROUGH controls); tauF, A2F = the earliest delay with
//      S > sFirst and its amplitude;
//   4. tail charge excess Q = sum of r over t - t1 >= tailStart ns, in units of the first pulse (independent of the
//      filter and of tau).
//
// Efficiency (injection): every THROUGH crystal is processed a second time with an injected pulse a T(t - t1 - d) added
//   to both SiPMs, d uniform in [dtMin, injMaxDelay] ns, A2 uniform in [injMinA2, injMaxA2] ADC (fixed seed); tree
//   "inject" with the true and found values -> efficiency(delay, A2) for a given cut on maxS.
//
// mu+ veto: for every candidate, the CaloHits of the event in the same disk within vetoDist mm of the crystal and later
//   than its hit by [vetoMin, vetoMax] ns (decay electron of a mu+, tau = 2.2 us; below vetoMin the capture itself
//   lives): nMichel, the delay of the first one and their summed energy are stored; the candidates are not removed
//   (compare with and without the veto).
//
// Delayed deposits in the other crystals (tree "delayed"): for every candidate, the CaloRecoDigis of the event (single
//   SiPMs, time with the conditions' offsets, energy in MeV) in the same disk, not on the track and not the candidate
//   crystal, with dt = t_digi - t_hit of the candidate in [nbDtMin, nbDtMax] ns; ring 1 = within nbNear mm of the
//   candidate crystal, ring 2 = between nbFarMin and nbFarMax mm (accidental control). No first pulse to remove, and
//   no limit from the waveform length: capture (~85 ns), mu+ decays (2.2 us) and accidentals (flat) in one
//   distribution.
//
// Delayed cluster (tree "cluster"; decay electrons of the mu- in the muonic atom (DIO) and of early mu+, in the tail of
//   the stop): members = the second pulse of the candidate crystal (found when maxS > secS; energy = A2 x mean ADC2MeV
//   of its two SiPMs, time = t_hit + tau) and the crystals within clusR mm in the same disk, not on the track, with a
//   pulse on both SiPMs (CaloRecoDigis of the two SiPMs within pairDt ns; energy = mean of the two, MeV from the
//   conditions; time = mean), dt from the hit of the candidate in [dtMin, clusDtMax] ns; members in time order, a new
//   cluster after a gap > clusGap ns; one entry per cluster with E >= clusEmin MeV: dt (energy weighted), E, the part
//   in the candidate crystal, number of crystals, inWf (the cluster time is inside the candidate waveform, i.e. its
//   second pulse could have been seen). The threshold on E for the physics is eCut (counts printed at the end);
//   THROUGH candidates (kind 1) give the same clusters without a stop: tail and accidentals.
//
// Tags for the muonic X-rays (tree "capture", per candidate): the track (crystals, direction in the disk plane, crystal
//   one pitch after the STOP), the nearest LYSO of lysoCrystals in the disk (distance, on the track or not, per SiPM
//   the largest reco digi within |dt| < lysoWin of the candidate hit: energy, amplitude = energy / ADC2MeV, dt) and
//   the delayed activity around the candidate (cluster deposits with dt in [delMin, delMax] ns and E >= delEmin MeV:
//   number, summed energy, first dt): capture of a mu-.
//
// Stop template (stopTemplate: true): a second template built in the same way from the STOP candidates with a Michel
//   electron (mu+ that really stopped, no capture: most STOP candidates do not stop, so their mean is a crossing-muon
//   pulse; stopping muons have more slow light). Every candidate gets psd = chi2(stop template) / chi2(crossing
//   template) of its first pulse on [-20, 200] ns; a STOP candidate with psd < 1 (useTs) is searched with the stop
//   template, and injected (tree inject, kind 3).
//
// Output (TFileService):
//   tree "capture"  per candidate: kind (0 STOP, 1 THROUGH), cry, disk, run/subrun/event, amp (first pulse, ADC),
//                   t1, maxS, tau, A2, Q, ns (samples of the residual), nMichel, michelDt, michelE, tauF, A2F, and the
//                   tags (nTrk, ux, uy, cnext, stopE, LYSO, psd, useTs, delN, delE, delDt)
//   tree "inject"   the same search with an injected pulse (kind 2 THROUGH, 3 STOP), true injDelay and injA2
//   tree "delayed"  single-SiPM deposits in the other crystals of the disk (ring 1 near, 2 far)
//   tree "examples" example waveforms: cls 0 STOP found, 1 THROUGH found, 2 STOP not found
//   tree "cluster"  delayed clusters around the candidates
//   histograms htempl (THROUGH template) and htempl_stop (stop template, with stopTemplate)
//
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art_root_io/TFileService.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"

#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/CalorimeterGeom/inc/Calorimeter.hh"
#include "Offline/CaloConditions/inc/CalCalib.hh"
#include "Offline/ProditionsService/inc/ProditionsHandle.hh"
#include "Offline/DataProducts/inc/CaloConst.hh"
#include "Offline/RecoDataProducts/inc/CaloDigi.hh"
#include "Offline/RecoDataProducts/inc/CaloRecoDigi.hh"
#include "Offline/RecoDataProducts/inc/CaloHit.hh"
#include "Offline/MC2Reco/inc/MC2CaloCosmicTrack.hh"
#include "TH1D.h"
#include "TRandom3.h"
#include "TTree.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <map>
#include <set>
#include <vector>

namespace mu2e {

  namespace {

    constexpr double kSamplePeriod = 5.0;        // ns, digitizer period
    constexpr double kTemplateXMin = -60.;       // ns from the 50 % crossing, start of the template
    constexpr double kTemplateBinWidth = 0.5;    // ns
    constexpr double kCrystalPitch = 34.3;       // mm, to find the crystal after the STOP
    constexpr double kNextCrystalMaxDist = 20.;  // mm, max distance of that crystal from the extrapolated point
    constexpr double kEndMargin = 30.;           // ns, last delay searched / cluster "in the waveform": end - 30 ns
    constexpr double kScanStep = 2.5;            // ns, step of the matched-filter delay scan
    constexpr double kFirstPulseLo = -20.;       // ns, first-pulse fit window around the 50 % crossing
    constexpr double kFirstPulseHi = 40.;
    constexpr double kShapeHi = 200.;            // ns, end of the window of the pulse-shape chi2 (psd)
    constexpr double kTemplateMinAmp = 300.;     // ADC, amplitude range of the waveforms used for the templates
    constexpr double kTemplateMaxAmp = 1800.;
    constexpr size_t kMinTemplateWaveforms = 200;

    // Pulse template: bin contents on [kTemplateXMin, ...] in kTemplateBinWidth bins, linear interpolation between
    // the bin centres; the first bin value before the first centre, 0 after the last centre.
    struct PulseTemplate {
      std::vector<double> values;
      int nBins = 0;

      double at(double x) const {
        const double u = (x - kTemplateXMin) / kTemplateBinWidth - 0.5;
        if (u <= 0) {
          return values.empty() ? 0 : values.front();
        }
        const int i = int(u);
        if (i >= nBins - 1) {
          return 0.;
        }
        const double f = u - i;
        return values[i] * (1 - f) + values[i + 1] * f;
      }
    };

    // Waveform of one SiPM, pedestal subtracted and normalised to its amplitude.
    struct SiPMWaveform {
      std::vector<float> y;   // (raw - pedestal) / amp, per sample
      double t0;              // time of sample 0 [ns]
      float amp;              // largest sample - pedestal [ADC]
      double t50;             // 50 % crossing [ns from sample 0]
    };

    // A crystal deposit around a candidate.
    struct Deposit {
      float dt;   // from the hit of the candidate [ns]
      float e;    // [MeV]
      int cry;
    };

    // A STOP or THROUGH crystal, with everything taken from its event (the search itself is done in endJob).
    struct Candidate {
      int kind, cry, disk, run, sub, evt;
      SiPMWaveform w[2];
      // mu+ veto
      int nMichel = 0;
      float michelDt = -1;
      float michelE = 0;
      // calibration and timing of the crystal
      double a2m[2] = {0, 0};      // ADC2MeV of the two SiPMs
      double tEnd = 0;             // end of the shorter waveform from the hit time [ns]
      double stopE = 0;            // energy of the hit [MeV]
      std::vector<Deposit> nb;     // deposits of the other crystals for the delayed cluster
      // tags: track
      int nTrk = 0;
      int cnext = -1;
      float ux = 0;
      float uy = 0;
      // tags: nearest LYSO
      int lyCry = -1;
      int lyOnTrk = 0;
      float lyDist = -1;
      float lyE[2] = {0, 0};
      float lyA[2] = {0, 0};
      float lyDt[2] = {-999, -999};
      // tags: delayed activity
      int delN = 0;
      float delE = 0;
      float delDt = -1;
    };

    // Result of the second-pulse search on one candidate.
    struct PulseSearch {
      double amp = 0;          // first pulse, mean of the SiPMs [ADC]
      double t1 = 0;           // first pulse time, mean of the SiPMs [ns, absolute]
      double maxS = 0;         // largest matched-filter significance
      double tau = 0;          // delay of the second pulse at maxS [ns]
      double A2 = 0;           // amplitude of the second pulse at maxS [ADC]
      double Q = 0;            // tail charge excess (units of the first pulse)
      int ns = 0;              // samples of the residual
      double tauFirst = -1;    // earliest delay with S > sFirst [ns]
      double A2First = 0;      // its amplitude [ADC]
      double t1Sipm[2] = {0, 0};
      double A1Sipm[2] = {0, 0};   // first pulse amplitude of each SiPM [ADC]
    };

    // 50 % crossing of a normalised waveform: linear interpolation on the leading edge before the largest sample [ns].
    double halfMaxCrossingTime(const std::vector<float>& y) {
      const int im = int(std::max_element(y.begin(), y.end()) - y.begin());
      for (int q = im; q > 0; --q) {
        if (y[q - 1] < 0.5 && y[q] >= 0.5) {
          return (q - 1 + (0.5 - y[q - 1]) / (y[q] - y[q - 1])) * kSamplePeriod;
        }
      }
      return im * kSamplePeriod;
    }

    double lastSampleTime(const SiPMWaveform& w) {
      return w.t0 + kSamplePeriod * (w.y.size() - 1);
    }

    template <typename Point>
    double distanceInPlane(const Point& a, const Point& b) {
      return std::hypot(a.x() - b.x(), a.y() - b.y());
    }

    // Template from waveforms, each shifted by its alignment time: mean per bin, 5-bin moving average, normalised to
    // the peak, then shifted so that its 50 % crossing is at 0.
    PulseTemplate buildTemplate(const std::vector<const SiPMWaveform*>& waveforms,
                                const std::vector<double>& alignTimes, int nBins) {
      // mean per bin
      std::vector<double> sum(nBins, 0), count(nBins, 0);
      for (size_t i = 0; i < waveforms.size(); ++i) {
        for (size_t q = 0; q < waveforms[i]->y.size(); ++q) {
          const double x = q * kSamplePeriod - alignTimes[i];
          const int b = int((x - kTemplateXMin) / kTemplateBinWidth);
          if (b < 0 || b >= nBins) {
            continue;
          }
          sum[b] += waveforms[i]->y[q];
          count[b] += 1;
        }
      }
      std::vector<double> mean(nBins, 0);
      for (int b = 0; b < nBins; ++b) {
        mean[b] = count[b] > 0 ? sum[b] / count[b] : 0;
      }

      // smoothing over the filled bins of [b - 2, b + 2]
      PulseTemplate smoothed;
      smoothed.nBins = nBins;
      smoothed.values.assign(nBins, 0);
      for (int b = 0; b < nBins; ++b) {
        double total = 0;
        int nUsed = 0;
        for (int k = -2; k <= 2; ++k) {
          if (b + k >= 0 && b + k < nBins && count[b + k] > 0) {
            total += mean[b + k];
            ++nUsed;
          }
        }
        smoothed.values[b] = nUsed ? total / nUsed : 0;
      }

      // normalisation to the peak
      const double peak = *std::max_element(smoothed.values.begin(), smoothed.values.end());
      if (peak > 0) {
        for (auto& x : smoothed.values) {
          x /= peak;
        }
      }

      // 50 % crossing on the leading edge
      const int im = int(std::max_element(smoothed.values.begin(), smoothed.values.end()) - smoothed.values.begin());
      double x50 = 0;
      for (int b = im; b > 0; --b) {
        const auto& v = smoothed.values;
        if (v[b - 1] < 0.5 && v[b] >= 0.5) {
          x50 = kTemplateXMin + (b - 0.5 + (0.5 - v[b - 1]) / (v[b] - v[b - 1])) * kTemplateBinWidth;
          break;
        }
      }

      // shift: 50 % crossing at 0
      PulseTemplate shifted;
      shifted.nBins = nBins;
      shifted.values.resize(nBins);
      for (int b = 0; b < nBins; ++b) {
        const double x = kTemplateXMin + (b + 0.5) * kTemplateBinWidth + x50;
        shifted.values[b] = (x - kTemplateXMin) / kTemplateBinWidth - 0.5 >= nBins - 1 ? 0 : smoothed.at(x);
      }
      return shifted;
    }

    // chi2 (normalised units) of the template at time t on the samples with t_sample - t in [lo, hi]; the amplitude
    // is the linear least-squares one.
    double templateChi2(const SiPMWaveform& w, const PulseTemplate& T, double lo, double hi, double t,
                        double& amplitude) {
      double sty = 0, stt = 0;
      for (size_t q = 0; q < w.y.size(); ++q) {
        const double x = q * kSamplePeriod - t;
        if (x < lo || x > hi) {
          continue;
        }
        const double m = T.at(x);
        sty += m * w.y[q];
        stt += m * m;
      }
      amplitude = stt > 0 ? sty / stt : 0;
      double chi2 = 0;
      for (size_t q = 0; q < w.y.size(); ++q) {
        const double x = q * kSamplePeriod - t;
        if (x < lo || x > hi) {
          continue;
        }
        const double r = w.y[q] - amplitude * T.at(x);
        chi2 += r * r;
      }
      return chi2;
    }

    // First pulse: template fit on [lo, hi] around the time, scanning +-5 ns around the 50 % crossing in 0.1 ns steps.
    // Returns the chi2 (normalised units); bestTime and bestAmplitude (normalised to w.amp) of the best scan point.
    double fitFirstPulse(const SiPMWaveform& w, const PulseTemplate& T, double lo, double hi,
                         double& bestTime, double& bestAmplitude) {
      double best = 1e30, amplitude;
      bestTime = w.t50;
      for (double t = w.t50 - 5; t <= w.t50 + 5 + 1e-9; t += 0.1) {
        const double chi2 = templateChi2(w, T, lo, hi, t, amplitude);
        if (chi2 < best) {
          best = chi2;
          bestTime = t;
        }
      }
      templateChi2(w, T, lo, hi, bestTime, bestAmplitude);
      return best;
    }

  }

  class MC2MuCapture : public art::EDAnalyzer {
    public:
      struct Config {
        using Name = fhicl::Name;
        using Comment = fhicl::Comment;

        // Inputs
        fhicl::Atom<art::InputTag> trackTag{Name("caloTrackCollection"), Comment("MC2CaloCosmicTracks"),
            "MC2CaloCosmicTrackFinder"};
        fhicl::Atom<art::InputTag> digiTag{Name("caloDigiCollection"),
            Comment("CaloDigis (all waveforms of the event)"), "CaloDigis"};
        fhicl::Atom<art::InputTag> hitTag{Name("caloHitCollection"), Comment("CaloHits (mu+ veto)"),
            "MC2CaloHitMaker"};
        fhicl::Atom<art::InputTag> recoDigiTag{Name("caloRecoDigiCollection"),
            Comment("reco digis (delayed deposits in the other crystals)"), "MC2CaloRecoDigiMaker"};

        // Delayed deposits in the other crystals (tree delayed)
        fhicl::Atom<double> nbNear{Name("nbNear"), Comment("ring 1: max distance [mm]"), 80.};
        fhicl::Atom<double> nbFarMin{Name("nbFarMin"), Comment("ring 2: min distance [mm]"), 250.};
        fhicl::Atom<double> nbFarMax{Name("nbFarMax"), Comment("ring 2: max distance [mm]"), 400.};
        fhicl::Atom<double> nbDtMin{Name("nbDtMin"), Comment("min dt [ns]"), -1000.};
        fhicl::Atom<double> nbDtMax{Name("nbDtMax"), Comment("max dt [ns]"), 10000.};

        // mu+ veto
        fhicl::Atom<double> vetoDist{Name("vetoDist"), Comment("mu+ veto: max distance from the crystal [mm]"), 80.};
        fhicl::Atom<double> vetoMin{Name("vetoMin"), Comment("mu+ veto: min delay [ns]"), 1000.};
        fhicl::Atom<double> vetoMax{Name("vetoMax"), Comment("mu+ veto: max delay [ns]"), 20000.};

        // Candidates and waveforms
        fhicl::Atom<int> minCrystals{Name("minCrystals"), Comment("crystals of a one-disk track for a STOP"), 4};
        fhicl::Atom<int> minThrough{Name("minThrough"), Comment("crystals of a track for a THROUGH"), 6};
        fhicl::Atom<int> maxThrough{Name("maxThrough"), Comment("max THROUGH crystals kept"), 400000};
        fhicl::Atom<int> minSamples{Name("minSamples"), Comment("min waveform length [samples]"), 16};
        fhicl::Atom<double> minAmp{Name("minAmp"), Comment("min first-pulse amplitude [ADC]"), 100.};
        fhicl::Atom<int> satADC{Name("satADC"), Comment("saturation [raw ADC]"), 4090};

        // Templates
        fhicl::Atom<double> xMax{Name("xMax"), Comment("template range [ns]"), 1500.};
        fhicl::Atom<int> nIterTemplate{Name("nIterTemplate"), Comment("template iterations"), 2};
        fhicl::Atom<bool> stopTemplate{Name("stopTemplate"),
            Comment("search (and inject) the STOP crystals with a template built from the STOP waveforms"), false};

        // Second-pulse search
        fhicl::Atom<double> dtMin{Name("dtMin"), Comment("min delay of the second pulse [ns]"), 10.};
        fhicl::Atom<double> noiseADC{Name("noiseADC"), Comment("noise per sample [ADC]"), 3.75};
        fhicl::Atom<double> tailStart{Name("tailStart"), Comment("start of the tail charge [ns after t1]"), 60.};

        // Injection
        fhicl::Atom<double> injMaxDelay{Name("injMaxDelay"), Comment("injection: max delay [ns]"), 600.};
        fhicl::Atom<double> injMinA2{Name("injMinA2"), Comment("injection: min A2 [ADC]"), 20.};
        fhicl::Atom<double> injMaxA2{Name("injMaxA2"), Comment("injection: max A2 [ADC]"), 400.};
        fhicl::Atom<int> seed{Name("seed"), Comment("injection seed"), 12345};

        // Example waveforms (tree examples)
        fhicl::Atom<int> nExamples{Name("nExamples"),
            Comment("waveforms saved per class (tree examples): STOP found x5, THROUGH found, STOP not found"), 24};
        fhicl::Atom<double> exampleS{Name("exampleS"), Comment("found = maxS > exampleS and A2 > exampleA2"), 5.};
        fhicl::Atom<double> exampleA2{Name("exampleA2"), Comment("[ADC]"), 30.};

        // First second pulse
        fhicl::Atom<double> sFirst{Name("sFirst"),
            Comment("first second pulse: earliest delay with matched-filter significance above sFirst "
                    "(branches tauF, A2F)"), 8.};

        // Delayed clusters (tree cluster)
        fhicl::Atom<double> secS{Name("secS"), Comment("cluster: second pulse of the candidate crystal if maxS > secS"),
            5.};
        fhicl::Atom<double> clusR{Name("clusR"),
            Comment("cluster: max distance of the other crystals from the candidate [mm]"), 80.};
        fhicl::Atom<double> pairDt{Name("pairDt"), Comment("cluster: max |tL - tR| of the two SiPMs of a crystal [ns]"),
            10.};
        fhicl::Atom<double> clusGap{Name("clusGap"), Comment("cluster: new cluster after a gap > clusGap [ns]"), 10.};
        fhicl::Atom<double> clusDtMax{Name("clusDtMax"), Comment("cluster: max dt from the stop [ns]"), 20000.};
        fhicl::Atom<double> clusEmin{Name("clusEmin"), Comment("cluster: min energy written to the tree [MeV]"), 1.};
        fhicl::Atom<double> eCut{Name("eCut"),
            Comment("cluster: energy threshold of the physics sample [MeV] (summary counts)"), 15.};

        // Tags: LYSO and delayed activity
        fhicl::Sequence<int> lysoCrystals{Name("lysoCrystals"), Comment("crystal ids of the LYSO (CAPHRI)"),
            std::vector<int>{582, 609, 610, 637}};
        fhicl::Atom<double> lysoWin{Name("lysoWin"), Comment("LYSO deposits: max |dt| from the candidate hit [ns]"),
            50.};
        fhicl::Atom<double> delMin{Name("delMin"), Comment("delayed activity: min dt [ns]"), 20.};
        fhicl::Atom<double> delMax{Name("delMax"), Comment("delayed activity: max dt [ns]"), 400.};
        fhicl::Atom<double> delEmin{Name("delEmin"), Comment("delayed activity: min deposit energy [MeV]"), 1.};

        // Crystals left out of the analysis (e.g. the LYSO at operating voltage, whose intrinsic activity gives a
        // continuous rate of deposits): never a candidate, never in the mu+ veto, the delayed deposits or the
        // delayed clusters (the LYSO tag above still reads them); empty: all crystals
        fhicl::Sequence<int> excludedCrystals{Name("excludedCrystals"), Comment("crystal ids left out"),
            std::vector<int>{}};
      };

      using Parameters = art::EDAnalyzer::Table<Config>;
      explicit MC2MuCapture(const Parameters& config);
      void beginJob() override;
      void analyze(const art::Event& event) override;
      void endJob() override;

    private:
      // The products of one event used to build the candidates.
      struct EventData {
        const art::Event& event;
        const Calorimeter& cal;
        const CalCalib& calib;
        const CaloHitCollection& hits;
        const CaloRecoDigiCollection& recoDigis;
      };

      // Booking
      void bookSearchBranches(TTree* tree);
      void bookTagBranches();
      void bookDelayedTree();
      void bookExampleTree();
      void bookClusterTree();

      // analyze: candidates
      bool readWaveform(const CaloDigi& digi, SiPMWaveform& w) const;
      void addCandidate(const EventData& data, const MC2CaloCosmicTrack& track, double ux, double uy, int kind,
                        size_t hitIndex, int cnext);
      bool readCandidateWaveforms(const CaloHit& hit, Candidate& c) const;
      void applyMuPlusVeto(const EventData& data, const CaloHit& hit, Candidate& c) const;
      void fillDelayedDeposits(const EventData& data, const CaloHit& hit, const std::set<int>& onTrack,
                               const Candidate& c);
      void collectClusterDeposits(const EventData& data, const CaloHit& hit, const std::set<int>& onTrack,
                                  Candidate& c) const;
      void tagNearestLyso(const EventData& data, const CaloHit& hit, const std::set<int>& onTrack,
                          Candidate& c) const;
      void tagDelayedActivity(Candidate& c) const;
      int crystalAfterStop(const Calorimeter& cal, int lastCrystal, double ux, double uy) const;

      // endJob: templates and search
      PulseTemplate makeTemplate(const std::vector<const SiPMWaveform*>& waveforms, int nBins) const;
      PulseSearch searchSecondPulse(const Candidate& c, const PulseTemplate& T) const;
      void fillSearchBranches(const PulseSearch& search);
      void fillCaptureTree(const Candidate& c, const PulseSearch& search);
      void fillClusters(const Candidate& c, double maxS, double tau, double A2);
      void saveExample(const Candidate& c, const PulseSearch& search);
      void fillInjection(const Candidate& c, const PulseTemplate& T, const PulseTemplate& searchTemplate,
                         TRandom3& rnd);

      // Inputs
      art::InputTag trackTag_, digiTag_, hitTag_, recoDigiTag_;
      ProditionsHandle<CalCalib> calibHandle_;

      // Parameters: mu+ veto and delayed deposits
      double vetoDist_, vetoMin_, vetoMax_;
      double nbNear_, nbFarMin_, nbFarMax_, nbDtMin_, nbDtMax_;
      // Parameters: candidates and waveforms
      int minCrystals_, minThrough_, maxThrough_, minSamples_;
      double minAmp_;
      int satADC_;
      // Parameters: templates, search, injection, examples
      double xMax_;
      int nIterTemplate_;
      bool stopTemplate_ = false;
      double dtMin_, noiseADC_, tailStart_;
      double injMaxDelay_, injMinA2_, injMaxA2_;
      int seed_;
      int nExamples_;
      double exampleS_, exampleA2_;
      double sFirst_;
      // Parameters: clusters and tags
      double secS_, clusR_, pairDt_, clusGap_, clusDtMax_, clusEmin_, eCut_;
      std::vector<int> lysoCrystals_;
      double lysoWin_, delMin_, delMax_, delEmin_;
      std::set<int> excludedCrystals_;

      // State
      std::vector<Candidate> candidates_;
      long nThrough_ = 0;              // THROUGH candidates tried (limited by maxThrough)
      long nStopTried_ = 0;
      long nStopNextHasDigi_ = 0;      // STOP rejected: the crystal after it has a digi
      long nClustersAboveCut_[2] = {0, 0};   // clusters with E >= eCut, STOP and THROUGH
      int nExamplesSaved_[3] = {0, 0, 0};

      // Output: trees and templates
      TTree* captureTree_ = nullptr;
      TTree* injectTree_ = nullptr;
      TTree* examplesTree_ = nullptr;
      TTree* delayedTree_ = nullptr;
      TTree* clusterTree_ = nullptr;
      TH1D* hTemplate_ = nullptr;
      TH1D* hStopTemplate_ = nullptr;

      // Branches shared by capture, inject (and partly examples)
      int kind_, cry_, disk_, run_, sub_, evt_, ns_, nMi_;
      float amp_, t1_, maxS_, tau_, A2_, Q_, miDt_, miE_;
      float tTauF_ = -1, tA2F_ = 0;
      // inject only
      float injD_, injA2_;
      // capture only: tags
      int tNTrk_, tCnext_, tLyCry_, tLyOnTrk_, tDelN_, tUseTs_ = 0;
      float tPsd_ = -1;
      float tUx_, tUy_, tStopE_, tLyDist_, tLyE_[2], tLyA_[2], tLyDt_[2], tDelE_, tDelDt_;
      // delayed
      int dKind_, dCry_, dDisk_, dNMi_, dRing_, dNbCry_, dSiPM_;
      float dDt_, dE_, dDist_, dStopE_;
      // examples
      int exCls_;
      std::vector<float> exY0_, exY1_;
      float exT0a_, exT0b_, exA1a_, exA1b_, exT1a_, exT1b_, exTref_;
      // cluster
      int cKind_, cCry_, cDisk_, cRun_, cSub_, cEvt_, cN_, cInWf_, cNMi_;
      float cDt_, cE_, cEstop_, cStopE_, cMaxS_;
  };

  MC2MuCapture::MC2MuCapture(const Parameters& config) :
    art::EDAnalyzer(config),
    trackTag_(config().trackTag()),
    digiTag_(config().digiTag()),
    hitTag_(config().hitTag()),
    recoDigiTag_(config().recoDigiTag()),
    vetoDist_(config().vetoDist()),
    vetoMin_(config().vetoMin()),
    vetoMax_(config().vetoMax()),
    nbNear_(config().nbNear()),
    nbFarMin_(config().nbFarMin()),
    nbFarMax_(config().nbFarMax()),
    nbDtMin_(config().nbDtMin()),
    nbDtMax_(config().nbDtMax()),
    minCrystals_(config().minCrystals()),
    minThrough_(config().minThrough()),
    maxThrough_(config().maxThrough()),
    minSamples_(config().minSamples()),
    minAmp_(config().minAmp()),
    satADC_(config().satADC()),
    xMax_(config().xMax()),
    nIterTemplate_(config().nIterTemplate()),
    stopTemplate_(config().stopTemplate()),
    dtMin_(config().dtMin()),
    noiseADC_(config().noiseADC()),
    tailStart_(config().tailStart()),
    injMaxDelay_(config().injMaxDelay()),
    injMinA2_(config().injMinA2()),
    injMaxA2_(config().injMaxA2()),
    seed_(config().seed()),
    nExamples_(config().nExamples()),
    exampleS_(config().exampleS()),
    exampleA2_(config().exampleA2()),
    sFirst_(config().sFirst()),
    secS_(config().secS()),
    clusR_(config().clusR()),
    pairDt_(config().pairDt()),
    clusGap_(config().clusGap()),
    clusDtMax_(config().clusDtMax()),
    clusEmin_(config().clusEmin()),
    eCut_(config().eCut()),
    lysoCrystals_(config().lysoCrystals()),
    lysoWin_(config().lysoWin()),
    delMin_(config().delMin()),
    delMax_(config().delMax()),
    delEmin_(config().delEmin())
  {
    const auto excluded = config().excludedCrystals();   // returned by value: one copy for begin and end
    excludedCrystals_.insert(excluded.begin(), excluded.end());
    consumes<MC2CaloCosmicTrackCollection>(trackTag_);
    consumes<CaloDigiCollection>(digiTag_);
    consumes<CaloHitCollection>(hitTag_);
    consumes<CaloRecoDigiCollection>(recoDigiTag_);
  }

  //--------------------------------------------------------------------------------------------------------------------
  // Booking (the order of trees, branches and histograms is that of the output files)
  //--------------------------------------------------------------------------------------------------------------------

  void MC2MuCapture::beginJob() {
    art::ServiceHandle<art::TFileService> tfs;
    captureTree_ = tfs->make<TTree>("capture", "second-pulse search in STOP (kind 0) and THROUGH (kind 1) crystals");
    injectTree_ = tfs->make<TTree>("inject", "THROUGH crystals with an injected second pulse");
    for (TTree* tree : {captureTree_, injectTree_}) {
      bookSearchBranches(tree);
    }
    injectTree_->Branch("injDelay", &injD_);
    injectTree_->Branch("injA2", &injA2_);
    bookTagBranches();
    bookDelayedTree();
    bookExampleTree();
    bookClusterTree();

    const int nBins = int((xMax_ - kTemplateXMin) / kTemplateBinWidth);
    if (stopTemplate_) {
      hStopTemplate_ = tfs->make<TH1D>("htempl_stop",
          "template from the STOP crystals of the run;time - 50% crossing [ns];normalised", nBins, kTemplateXMin,
          xMax_);
    }
    hTemplate_ = tfs->make<TH1D>("htempl",
        "template from the THROUGH crystals of the run;time - 50% crossing [ns];normalised", nBins, kTemplateXMin,
        xMax_);
  }

  void MC2MuCapture::bookSearchBranches(TTree* tree) {
    tree->Branch("kind", &kind_);
    tree->Branch("cry", &cry_);
    tree->Branch("disk", &disk_);
    tree->Branch("run", &run_);
    tree->Branch("subrun", &sub_);
    tree->Branch("event", &evt_);
    tree->Branch("amp", &amp_);
    tree->Branch("t1", &t1_);
    tree->Branch("maxS", &maxS_);
    tree->Branch("tau", &tau_);
    tree->Branch("A2", &A2_);
    tree->Branch("Q", &Q_);
    tree->Branch("ns", &ns_);
    tree->Branch("nMichel", &nMi_);
    tree->Branch("michelDt", &miDt_);
    tree->Branch("michelE", &miE_);
    tree->Branch("tauF", &tTauF_);
    tree->Branch("A2F", &tA2F_);
  }

  void MC2MuCapture::bookTagBranches() {
    captureTree_->Branch("nTrk", &tNTrk_);
    captureTree_->Branch("ux", &tUx_);
    captureTree_->Branch("uy", &tUy_);
    captureTree_->Branch("cnext", &tCnext_);
    captureTree_->Branch("stopE", &tStopE_);
    captureTree_->Branch("lyCry", &tLyCry_);
    captureTree_->Branch("lyDist", &tLyDist_);
    captureTree_->Branch("lyOnTrk", &tLyOnTrk_);
    captureTree_->Branch("lyE", tLyE_, "lyE[2]/F");
    captureTree_->Branch("lyA", tLyA_, "lyA[2]/F");
    captureTree_->Branch("lyDt", tLyDt_, "lyDt[2]/F");
    captureTree_->Branch("psd", &tPsd_);
    captureTree_->Branch("useTs", &tUseTs_);
    captureTree_->Branch("delN", &tDelN_);
    captureTree_->Branch("delE", &tDelE_);
    captureTree_->Branch("delDt", &tDelDt_);
  }

  void MC2MuCapture::bookDelayedTree() {
    art::ServiceHandle<art::TFileService> tfs;
    delayedTree_ = tfs->make<TTree>("delayed",
        "single-SiPM deposits in the other crystals of the disk, dt from the candidate (ring 1 near, 2 far)");
    delayedTree_->Branch("kind", &dKind_);
    delayedTree_->Branch("cry", &dCry_);
    delayedTree_->Branch("disk", &dDisk_);
    delayedTree_->Branch("nMichel", &dNMi_);
    delayedTree_->Branch("ring", &dRing_);
    delayedTree_->Branch("nbCry", &dNbCry_);
    delayedTree_->Branch("sipm", &dSiPM_);
    delayedTree_->Branch("dt", &dDt_);
    delayedTree_->Branch("e", &dE_);
    delayedTree_->Branch("dist", &dDist_);
    delayedTree_->Branch("stopE", &dStopE_);
  }

  void MC2MuCapture::bookExampleTree() {
    art::ServiceHandle<art::TFileService> tfs;
    examplesTree_ = tfs->make<TTree>("examples",
        "example waveforms: cls 0 STOP found, 1 THROUGH found, 2 STOP not found");
    examplesTree_->Branch("cls", &exCls_);
    examplesTree_->Branch("kind", &kind_);
    examplesTree_->Branch("cry", &cry_);
    examplesTree_->Branch("disk", &disk_);
    examplesTree_->Branch("run", &run_);
    examplesTree_->Branch("subrun", &sub_);
    examplesTree_->Branch("event", &evt_);
    examplesTree_->Branch("amp", &amp_);
    examplesTree_->Branch("maxS", &maxS_);
    examplesTree_->Branch("tau", &tau_);
    examplesTree_->Branch("A2", &A2_);
    examplesTree_->Branch("Q", &Q_);
    examplesTree_->Branch("nMichel", &nMi_);
    examplesTree_->Branch("y0", &exY0_);   // raw - pedestal [ADC]
    examplesTree_->Branch("y1", &exY1_);
    examplesTree_->Branch("t0a", &exT0a_);   // t0 of the digis [ns]
    examplesTree_->Branch("t0b", &exT0b_);
    examplesTree_->Branch("A1a", &exA1a_);
    examplesTree_->Branch("A1b", &exA1b_);
    examplesTree_->Branch("t1a", &exT1a_);
    examplesTree_->Branch("t1b", &exT1b_);
    examplesTree_->Branch("tref", &exTref_);
  }

  void MC2MuCapture::bookClusterTree() {
    art::ServiceHandle<art::TFileService> tfs;
    clusterTree_ = tfs->make<TTree>("cluster",
        "delayed clusters around STOP (kind 0) and THROUGH (kind 1) crystals: "
        "second pulse of the crystal + neighbours");
    clusterTree_->Branch("kind", &cKind_);
    clusterTree_->Branch("cry", &cCry_);
    clusterTree_->Branch("disk", &cDisk_);
    clusterTree_->Branch("run", &cRun_);
    clusterTree_->Branch("subrun", &cSub_);
    clusterTree_->Branch("event", &cEvt_);
    clusterTree_->Branch("dt", &cDt_);
    clusterTree_->Branch("E", &cE_);
    clusterTree_->Branch("Estop", &cEstop_);
    clusterTree_->Branch("n", &cN_);
    clusterTree_->Branch("inWf", &cInWf_);
    clusterTree_->Branch("stopE", &cStopE_);
    clusterTree_->Branch("maxS", &cMaxS_);
    clusterTree_->Branch("nMichel", &cNMi_);
  }

  //--------------------------------------------------------------------------------------------------------------------
  // analyze: STOP and THROUGH candidates of the event
  //--------------------------------------------------------------------------------------------------------------------

  void MC2MuCapture::analyze(const art::Event& event) {
    const Calorimeter& cal = *GeomHandle<Calorimeter>();
    const auto& tracks = *event.getValidHandle<MC2CaloCosmicTrackCollection>(trackTag_);
    if (tracks.empty()) {
      return;
    }
    const auto& digis = *event.getValidHandle<CaloDigiCollection>(digiTag_);
    const auto& hits = *event.getValidHandle<CaloHitCollection>(hitTag_);
    const auto& recoDigis = *event.getValidHandle<CaloRecoDigiCollection>(recoDigiTag_);
    const CalCalib& calib = calibHandle_.get(event.id());
    const EventData data{event, cal, calib, hits, recoDigis};

    std::set<int> crystalsWithDigi;
    for (const auto& digi : digis) {
      crystalsWithDigi.insert(digi.SiPMID() / CaloConst::_nSiPMPerCrystal);
    }

    for (const auto& track : tracks) {
      if ((int)track.hits.size() < std::min(minCrystals_, minThrough_)) {
        continue;
      }

      // direction of the track in the disk plane, pointing down
      double ux = track.ux, uy = track.uy;
      const double norm = std::hypot(ux, uy);
      if (norm <= 0) {
        continue;
      }
      ux /= norm;
      uy /= norm;
      if (uy > 0) {
        ux = -ux;
        uy = -uy;
      }

      // crystals of the track ordered along the downward direction: (projection, hit index)
      std::vector<std::pair<double, size_t>> alongTrack;
      for (size_t i = 0; i < track.hits.size(); ++i) {
        const auto p = cal.mu2eToTracker(cal.crystal(track.hits[i]->crystalID()).position());
        alongTrack.push_back({p.x() * ux + p.y() * uy, i});
      }
      std::sort(alongTrack.begin(), alongTrack.end());

      // STOP: the last crystal of a one-disk track, if the next crystal along the track exists, is not on the track
      // and has no digi
      if ((track.type == 1 || track.type == 2) && (int)alongTrack.size() >= minCrystals_) {
        ++nStopTried_;
        const size_t lastHit = alongTrack.back().second;
        const int cnext = crystalAfterStop(cal, track.hits[lastHit]->crystalID(), ux, uy);
        bool nextOnTrack = false;
        for (const auto& hit : track.hits) {
          if (hit->crystalID() == cnext) {
            nextOnTrack = true;
          }
        }
        if (cnext >= 0 && !nextOnTrack) {
          if (crystalsWithDigi.count(cnext)) {
            ++nStopNextHasDigi_;
          } else {
            addCandidate(data, track, ux, uy, 0, lastHit, cnext);
          }
        }
      }

      // THROUGH: the crystal in the middle of a long track
      if ((int)alongTrack.size() >= minThrough_ && nThrough_ < maxThrough_) {
        addCandidate(data, track, ux, uy, 1, alongTrack[alongTrack.size() / 2].second, -1);
        ++nThrough_;
      }
    }
  }

  // The crystal one pitch after the last one along the downward direction, in the same disk (nearest centre within
  // 20 mm of the extrapolated point); -1 if none.
  int MC2MuCapture::crystalAfterStop(const Calorimeter& cal, int lastCrystal, double ux, double uy) const {
    const int disk = cal.crystal(lastCrystal).diskID();
    const auto pl = cal.mu2eToTracker(cal.crystal(lastCrystal).position());
    const double px = pl.x() + kCrystalPitch * ux, py = pl.y() + kCrystalPitch * uy;
    int cnext = -1;
    double bestDist = kNextCrystalMaxDist;
    for (unsigned k = 0; k < cal.nCrystals(); ++k) {
      if ((int)cal.crystal(k).diskID() != disk) {
        continue;
      }
      const auto q = cal.mu2eToTracker(cal.crystal(k).position());
      const double dist = std::hypot(q.x() - px, q.y() - py);
      if (dist < bestDist) {
        bestDist = dist;
        cnext = k;
      }
    }
    return cnext;
  }

  // Waveform of one digi: pedestal = mean of the first 4 samples; false if too short, saturated or too small.
  bool MC2MuCapture::readWaveform(const CaloDigi& digi, SiPMWaveform& w) const {
    const auto& adc = digi.waveform();
    if ((int)adc.size() < minSamples_) {
      return false;
    }
    double pedestal = 0;
    for (int q = 0; q < 4; ++q) {
      pedestal += adc[q];
    }
    pedestal /= 4.;
    int im = 0;
    for (size_t q = 0; q < adc.size(); ++q) {
      if (adc[q] >= satADC_) {
        return false;
      }
      if (adc[q] > adc[im]) {
        im = q;
      }
    }
    const double amp = adc[im] - pedestal;
    if (amp < minAmp_) {
      return false;
    }
    w.y.resize(adc.size());
    for (size_t q = 0; q < adc.size(); ++q) {
      w.y[q] = float((adc[q] - pedestal) / amp);
    }
    w.amp = amp;
    w.t0 = digi.t0() * kSamplePeriod;
    w.t50 = halfMaxCrossingTime(w.y);
    return true;
  }

  // Both SiPMs of the hit must give a good waveform.
  bool MC2MuCapture::readCandidateWaveforms(const CaloHit& hit, Candidate& c) const {
    int nGood = 0;
    for (const auto& recoDigi : hit.recoCaloDigis()) {
      const int s = recoDigi->SiPMID() % CaloConst::_nSiPMPerCrystal;
      if (readWaveform(*recoDigi->caloDigiPtr(), c.w[s])) {
        ++nGood;
      }
    }
    return nGood == 2;
  }

  // A candidate with everything taken from the event; the delayed deposits are written to their tree right away.
  void MC2MuCapture::addCandidate(const EventData& data, const MC2CaloCosmicTrack& track, double ux, double uy,
                                  int kind, size_t hitIndex, int cnext) {
    const Calorimeter& cal = data.cal;
    const auto& hit = track.hits[hitIndex];
    if (excludedCrystals_.count(hit->crystalID())) {
      return;
    }
    Candidate c;
    c.kind = kind;
    c.cry = hit->crystalID();
    c.disk = cal.crystal(c.cry).diskID();
    c.run = data.event.run();
    c.sub = data.event.subRun();
    c.evt = data.event.event();

    if (!readCandidateWaveforms(*hit, c)) {
      return;
    }
    for (const auto& recoDigi : hit->recoCaloDigis()) {
      c.a2m[recoDigi->SiPMID() % CaloConst::_nSiPMPerCrystal] = data.calib.ADC2MeV(recoDigi->SiPMID());
    }
    // approximate: digi t0 vs calibrated hit time
    c.tEnd = std::min(lastSampleTime(c.w[0]), lastSampleTime(c.w[1])) - hit->time();

    applyMuPlusVeto(data, *hit, c);

    std::set<int> onTrack;
    for (const auto& trackHit : track.hits) {
      onTrack.insert(trackHit->crystalID());
    }
    fillDelayedDeposits(data, *hit, onTrack, c);
    collectClusterDeposits(data, *hit, onTrack, c);
    c.stopE = hit->energyDep();

    // tags: track, nearest LYSO and its prompt deposits, delayed activity around the candidate
    c.nTrk = track.hits.size();
    c.ux = ux;
    c.uy = uy;
    c.cnext = cnext;
    tagNearestLyso(data, *hit, onTrack, c);
    tagDelayedActivity(c);

    candidates_.push_back(std::move(c));
  }

  // mu+ veto: later deposits near the crystal.
  void MC2MuCapture::applyMuPlusVeto(const EventData& data, const CaloHit& hit, Candidate& c) const {
    const Calorimeter& cal = data.cal;
    const auto pc = cal.mu2eToTracker(cal.crystal(c.cry).position());
    for (const auto& other : data.hits) {
      const double dt = other.time() - hit.time();
      if (dt < vetoMin_ || dt > vetoMax_ || other.crystalID() >= (int)cal.nCrystals()) {
        continue;
      }
      if (excludedCrystals_.count(other.crystalID())) {
        continue;
      }
      const auto& otherCrystal = cal.crystal(other.crystalID());
      if ((int)otherCrystal.diskID() != c.disk) {
        continue;
      }
      const auto pa = cal.mu2eToTracker(otherCrystal.position());
      if (distanceInPlane(pa, pc) > vetoDist_) {
        continue;
      }
      ++c.nMichel;
      c.michelE += other.energyDep();
      if (c.michelDt < 0 || dt < c.michelDt) {
        c.michelDt = dt;
      }
    }
  }

  // Tree "delayed": single-SiPM reco digis of the other crystals of the disk (not on the track; PIN diodes skipped),
  // in the near ring or in the far ring. Written here: no waveform needed.
  void MC2MuCapture::fillDelayedDeposits(const EventData& data, const CaloHit& hit, const std::set<int>& onTrack,
                                         const Candidate& c) {
    const Calorimeter& cal = data.cal;
    const auto pc = cal.mu2eToTracker(cal.crystal(c.cry).position());
    for (const auto& recoDigi : data.recoDigis) {
      const int crystal = recoDigi.SiPMID() / CaloConst::_nSiPMPerCrystal;
      if (crystal >= (int)cal.nCrystals() || crystal == c.cry || onTrack.count(crystal)) {
        continue;
      }
      if (excludedCrystals_.count(crystal)) {
        continue;
      }
      const auto& otherCrystal = cal.crystal(crystal);
      if ((int)otherCrystal.diskID() != c.disk) {
        continue;
      }
      const double dt = recoDigi.time() - hit.time();
      if (dt < nbDtMin_ || dt > nbDtMax_) {
        continue;
      }
      const auto pa = cal.mu2eToTracker(otherCrystal.position());
      const double dist = distanceInPlane(pa, pc);
      const int ring = dist <= nbNear_ ? 1 : ((dist > nbFarMin_ && dist <= nbFarMax_) ? 2 : 0);
      if (!ring) {
        continue;
      }
      dKind_ = c.kind;
      dCry_ = c.cry;
      dDisk_ = c.disk;
      dNMi_ = c.nMichel;
      dRing_ = ring;
      dNbCry_ = crystal;
      dSiPM_ = recoDigi.SiPMID();
      dDt_ = dt;
      dE_ = recoDigi.energyDep();
      dDist_ = dist;
      dStopE_ = hit.energyDep();
      delayedTree_->Fill();
    }
  }

  // Crystal deposits for the delayed cluster: crystals within clusR with a reco digi on both SiPMs within pairDt (for
  // each digi of SiPM 0, the closest in time of SiPM 1); time and energy = mean of the two.
  void MC2MuCapture::collectClusterDeposits(const EventData& data, const CaloHit& hit, const std::set<int>& onTrack,
                                            Candidate& c) const {
    const Calorimeter& cal = data.cal;
    const auto pc = cal.mu2eToTracker(cal.crystal(c.cry).position());
    std::map<int, std::vector<const CaloRecoDigi*>> bySiPM[2];   // [SiPM of the crystal][crystal] -> reco digis
    for (const auto& recoDigi : data.recoDigis) {
      const int crystal = recoDigi.SiPMID() / CaloConst::_nSiPMPerCrystal;
      if (crystal >= (int)cal.nCrystals() || crystal == c.cry || onTrack.count(crystal)) {
        continue;
      }
      if (excludedCrystals_.count(crystal)) {
        continue;
      }
      const auto& otherCrystal = cal.crystal(crystal);
      if ((int)otherCrystal.diskID() != c.disk) {
        continue;
      }
      const auto pa = cal.mu2eToTracker(otherCrystal.position());
      if (distanceInPlane(pa, pc) > clusR_) {
        continue;
      }
      bySiPM[recoDigi.SiPMID() % CaloConst::_nSiPMPerCrystal][crystal].push_back(&recoDigi);
    }

    for (const auto& [crystal, digisSiPM0] : bySiPM[0]) {
      const auto found = bySiPM[1].find(crystal);
      if (found == bySiPM[1].end()) {
        continue;
      }
      for (const auto* a : digisSiPM0) {
        const CaloRecoDigi* b = nullptr;
        for (const auto* x : found->second) {
          if (std::abs(x->time() - a->time()) < pairDt_
              && (!b || std::abs(x->time() - a->time()) < std::abs(b->time() - a->time()))) {
            b = x;
          }
        }
        if (!b) {
          continue;
        }
        const double dt = 0.5 * (a->time() + b->time()) - hit.time();
        if (dt < dtMin_ || dt > clusDtMax_) {
          continue;
        }
        c.nb.push_back({float(dt), float(0.5 * (a->energyDep() + b->energyDep())), crystal});
      }
    }
  }

  // Nearest LYSO crystal in the disk of the candidate and, per SiPM, its largest reco digi within lysoWin of the hit.
  void MC2MuCapture::tagNearestLyso(const EventData& data, const CaloHit& hit, const std::set<int>& onTrack,
                                    Candidate& c) const {
    const Calorimeter& cal = data.cal;
    const auto pc = cal.mu2eToTracker(cal.crystal(c.cry).position());
    for (int lyso : lysoCrystals_) {
      if (lyso < 0 || lyso >= (int)cal.nCrystals() || (int)cal.crystal(lyso).diskID() != c.disk) {
        continue;
      }
      const auto pl = cal.mu2eToTracker(cal.crystal(lyso).position());
      const float dist = distanceInPlane(pl, pc);
      if (c.lyCry < 0 || dist < c.lyDist) {
        c.lyCry = lyso;
        c.lyDist = dist;
      }
    }
    if (c.lyCry < 0) {
      return;
    }

    c.lyOnTrk = onTrack.count(c.lyCry) ? 1 : 0;
    for (const auto& recoDigi : data.recoDigis) {
      if (int(recoDigi.SiPMID() / CaloConst::_nSiPMPerCrystal) != c.lyCry) {
        continue;
      }
      const double dt = recoDigi.time() - hit.time();
      if (std::abs(dt) > lysoWin_) {
        continue;
      }
      const int s = recoDigi.SiPMID() % CaloConst::_nSiPMPerCrystal;
      if (recoDigi.energyDep() <= c.lyE[s]) {
        continue;
      }
      const double a2m = data.calib.ADC2MeV(recoDigi.SiPMID());
      c.lyE[s] = recoDigi.energyDep();
      c.lyA[s] = a2m > 0 ? recoDigi.energyDep() / a2m : -1;
      c.lyDt[s] = dt;
    }
  }

  // Delayed activity: cluster deposits with dt in [delMin, delMax] and E >= delEmin.
  void MC2MuCapture::tagDelayedActivity(Candidate& c) const {
    for (const auto& deposit : c.nb) {
      if (deposit.dt >= delMin_ && deposit.dt <= delMax_ && deposit.e >= delEmin_) {
        ++c.delN;
        c.delE += deposit.e;
        if (c.delDt < 0 || deposit.dt < c.delDt) {
          c.delDt = deposit.dt;
        }
      }
    }
  }

  //--------------------------------------------------------------------------------------------------------------------
  // Second-pulse search on one candidate
  //--------------------------------------------------------------------------------------------------------------------

  PulseSearch MC2MuCapture::searchSecondPulse(const Candidate& c, const PulseTemplate& T) const {
    PulseSearch result;

    // 1. first pulse of each SiPM
    double t1[2], A1[2];
    for (int s = 0; s < 2; ++s) {
      fitFirstPulse(c.w[s], T, kFirstPulseLo, kFirstPulseHi, t1[s], A1[s]);
    }
    for (int s = 0; s < 2; ++s) {
      result.t1Sipm[s] = t1[s];
      result.A1Sipm[s] = A1[s] * c.w[s].amp;
    }

    // 2. common grid: absolute times (t0 of the digi + sample), residual normalised to the first pulse, averaged over
    //    the SiPMs; times from tref, the mean first-pulse time
    const double ta = c.w[0].t0 + t1[0], tb = c.w[1].t0 + t1[1], tref = 0.5 * (ta + tb);
    const double lo = std::max(c.w[0].t0, c.w[1].t0), hi = std::min(lastSampleTime(c.w[0]), lastSampleTime(c.w[1]));
    std::vector<double> times, residual;
    for (double t = lo; t <= hi + 1e-6; t += kSamplePeriod) {
      double v = 0;
      for (int s = 0; s < 2; ++s) {
        const double qf = (t - c.w[s].t0) / kSamplePeriod;
        const int q = int(std::lround(qf));
        if (q < 0 || q >= (int)c.w[s].y.size()) {
          v = NAN;
          break;
        }
        const double tl = s == 0 ? ta : tb;
        v += (c.w[s].y[q] - A1[s] * T.at(t - tl)) / (A1[s] > 0 ? A1[s] : 1.);
      }
      if (std::isnan(v)) {
        continue;
      }
      times.push_back(t - tref);
      residual.push_back(0.5 * v);
    }
    result.amp = 0.5 * (A1[0] * c.w[0].amp + A1[1] * c.w[1].amp);
    result.t1 = tref;
    result.ns = residual.size();

    // noise of the averaged residual per sample, in units of the first pulse
    const double sigma = 0.5 * std::hypot(noiseADC_ / std::max(1e-3, A1[0] * c.w[0].amp),
                                          noiseADC_ / std::max(1e-3, A1[1] * c.w[1].amp));
    result.maxS = -1e9;
    result.tau = -1;
    result.A2 = 0;
    result.Q = 0;
    result.tauFirst = -1;
    result.A2First = 0;

    // 4. tail charge
    for (size_t i = 0; i < times.size(); ++i) {
      if (times[i] >= tailStart_) {
        result.Q += residual[i];
      }
    }
    if (times.empty()) {
      return result;
    }

    // 3. matched filter scanned over the delay
    const double tEnd = times.back() - kEndMargin;
    for (double d = dtMin_; d <= tEnd; d += kScanStep) {
      double sr = 0, s2 = 0;
      for (size_t i = 0; i < times.size(); ++i) {
        const double m = T.at(times[i] - d);
        if (m == 0) {
          continue;
        }
        sr += residual[i] * m;
        s2 += m * m;
      }
      if (s2 <= 0) {
        continue;
      }
      const double S = sr / (sigma * std::sqrt(s2));
      if (S > result.maxS) {
        result.maxS = S;
        result.tau = d;
        result.A2 = sr / s2 * result.amp;
      }
      if (result.tauFirst < 0 && S > sFirst_) {
        result.tauFirst = d;
        result.A2First = sr / s2 * result.amp;
      }
    }
    return result;
  }

  //--------------------------------------------------------------------------------------------------------------------
  // Delayed clusters of one candidate (tree "cluster")
  //--------------------------------------------------------------------------------------------------------------------

  void MC2MuCapture::fillClusters(const Candidate& c, double maxS, double tau, double A2) {
    // members: the deposits of the other crystals, plus the second pulse of the candidate crystal if found
    std::vector<Deposit> members = c.nb;
    double eStop2 = 0;
    if (maxS > secS_ && A2 > 0) {
      eStop2 = A2 * 0.5 * (c.a2m[0] + c.a2m[1]);
      members.push_back({float(tau), float(eStop2), c.cry});
    }
    std::sort(members.begin(), members.end(), [](const Deposit& a, const Deposit& b) { return a.dt < b.dt; });

    // clusters: runs of members in time order separated by gaps > clusGap
    for (size_t i = 0; i < members.size();) {
      size_t j = i + 1;
      while (j < members.size() && members[j].dt - members[j - 1].dt <= clusGap_) {
        ++j;
      }
      double e = 0, et = 0, es = 0;
      std::set<int> crystals;
      for (size_t k = i; k < j; ++k) {
        e += members[k].e;
        et += members[k].e * members[k].dt;
        crystals.insert(members[k].cry);
        if (members[k].cry == c.cry) {
          es += members[k].e;
        }
      }
      if (e >= clusEmin_ && e > 0) {
        cKind_ = c.kind;
        cCry_ = c.cry;
        cDisk_ = c.disk;
        cRun_ = c.run;
        cSub_ = c.sub;
        cEvt_ = c.evt;
        cDt_ = et / e;
        cE_ = e;
        cEstop_ = es;
        cN_ = crystals.size();
        cInWf_ = cDt_ < c.tEnd - kEndMargin;
        cStopE_ = c.stopE;
        cMaxS_ = maxS;
        cNMi_ = c.nMichel;
        clusterTree_->Fill();
        if (e >= eCut_ && c.kind <= 1) {
          ++nClustersAboveCut_[c.kind];
        }
      }
      i = j;
    }
  }

  //--------------------------------------------------------------------------------------------------------------------
  // endJob: templates, then the search on every candidate
  //--------------------------------------------------------------------------------------------------------------------

  // Template from a set of waveforms: first alignment at the 50 % crossing (spread uniformly over one sample by the
  // rank of its phase within the sample), then nIterTemplate iterations aligned with the first-pulse fit.
  PulseTemplate MC2MuCapture::makeTemplate(const std::vector<const SiPMWaveform*>& waveforms, int nBins) const {
    std::vector<double> phase;
    for (auto* w : waveforms) {
      phase.push_back(std::fmod(w->t50, kSamplePeriod));
    }
    std::vector<double> sortedPhase = phase;
    std::sort(sortedPhase.begin(), sortedPhase.end());
    std::vector<double> alignTimes(waveforms.size());
    for (size_t i = 0; i < waveforms.size(); ++i) {
      const auto rank = std::lower_bound(sortedPhase.begin(), sortedPhase.end(), phase[i]) - sortedPhase.begin();
      alignTimes[i] = waveforms[i]->t50 - phase[i] + kSamplePeriod * double(rank) / sortedPhase.size();
    }
    PulseTemplate T = buildTemplate(waveforms, alignTimes, nBins);
    for (int it = 0; it < nIterTemplate_; ++it) {
      for (size_t i = 0; i < waveforms.size(); ++i) {
        double bestTime, amplitude;
        fitFirstPulse(*waveforms[i], T, kFirstPulseLo, kFirstPulseHi, bestTime, amplitude);
        alignTimes[i] = bestTime;
      }
      T = buildTemplate(waveforms, alignTimes, nBins);
    }
    return T;
  }

  void MC2MuCapture::endJob() {
    const int nBins = int((xMax_ - kTemplateXMin) / kTemplateBinWidth);

    // 1. crossing-muon template from the THROUGH waveforms
    std::vector<const SiPMWaveform*> throughWaveforms;
    for (auto& c : candidates_) {
      if (c.kind == 1) {
        for (int s = 0; s < 2; ++s) {
          if (c.w[s].amp >= kTemplateMinAmp && c.w[s].amp <= kTemplateMaxAmp) {
            throughWaveforms.push_back(&c.w[s]);
          }
        }
      }
    }
    printf("MC2MuCapture: candidates %zu (STOP tried %ld, rejected because the next crystal has a digi %ld; "
           "THROUGH %ld); template from %zu waveforms\n",
           candidates_.size(), nStopTried_, nStopNextHasDigi_, nThrough_, throughWaveforms.size());
    if (throughWaveforms.size() < kMinTemplateWaveforms) {
      printf("MC2MuCapture: too few THROUGH waveforms for the template, nothing done\n");
      return;
    }
    PulseTemplate T = makeTemplate(throughWaveforms, nBins);
    for (int b = 0; b < nBins; ++b) {
      hTemplate_->SetBinContent(b + 1, T.values[b]);
    }

    // 2. stop template from the STOP candidates with a Michel electron (mu+ stops: really stopped, no capture)
    PulseTemplate Ts = T;
    bool haveTs = false;
    if (stopTemplate_) {
      std::vector<const SiPMWaveform*> stopWaveforms;
      for (auto& c : candidates_) {
        if (c.kind == 0 && c.nMichel > 0) {
          for (int s = 0; s < 2; ++s) {
            if (c.w[s].amp >= kTemplateMinAmp && c.w[s].amp <= kTemplateMaxAmp) {
              stopWaveforms.push_back(&c.w[s]);
            }
          }
        }
      }
      if (stopWaveforms.size() >= kMinTemplateWaveforms) {
        Ts = makeTemplate(stopWaveforms, nBins);
        haveTs = true;
        for (int b = 0; b < nBins; ++b) {
          hStopTemplate_->SetBinContent(b + 1, Ts.values[b]);
        }
      }
      printf("MC2MuCapture: stop template from %zu waveforms of STOP with Michel%s\n", stopWaveforms.size(),
             haveTs ? "" : " (too few, THROUGH template used)");
    }

    // 3. search on every candidate
    TRandom3 rnd(seed_);
    long nStop = 0, nThrough = 0;
    for (const auto& c : candidates_) {
      // pulse shape: chi2 of the two SiPMs on [-20, 200] ns with the stop template over the one with the
      // crossing-muon template; a STOP candidate is searched with the template that describes its first pulse better
      double chi2T = 0, chi2Ts = 0;
      if (haveTs) {
        for (int s = 0; s < 2; ++s) {
          double bestTime, amplitude;
          chi2T += fitFirstPulse(c.w[s], T, kFirstPulseLo, kShapeHi, bestTime, amplitude);
          chi2Ts += fitFirstPulse(c.w[s], Ts, kFirstPulseLo, kShapeHi, bestTime, amplitude);
        }
      }
      tPsd_ = (haveTs && chi2T > 0) ? float(chi2Ts / chi2T) : -1.f;
      tUseTs_ = (haveTs && c.kind == 0 && chi2Ts < chi2T) ? 1 : 0;
      const PulseTemplate& searchTemplate = tUseTs_ ? Ts : T;

      const PulseSearch search = searchSecondPulse(c, searchTemplate);
      fillCaptureTree(c, search);
      fillClusters(c, search.maxS, search.tau, search.A2);
      saveExample(c, search);
      (c.kind == 0 ? nStop : nThrough)++;

      // injection (THROUGH: kind 2; STOP: kind 3, with the search template chosen above)
      if (c.kind == 1 || (haveTs && c.kind == 0)) {
        fillInjection(c, T, searchTemplate, rnd);
      }
    }
    printf("MC2MuCapture: STOP crystals analysed %ld, THROUGH %ld (each also with an injected pulse)\n",
           nStop, nThrough);
    printf("MC2MuCapture: delayed clusters with E >= %.1f MeV: STOP %ld, THROUGH %ld\n", eCut_,
           nClustersAboveCut_[0], nClustersAboveCut_[1]);
  }

  // Branches of the search shared by the trees capture and inject.
  void MC2MuCapture::fillSearchBranches(const PulseSearch& search) {
    tTauF_ = float(search.tauFirst);
    tA2F_ = float(search.A2First);
    amp_ = search.amp;
    t1_ = search.t1;
    maxS_ = search.maxS;
    tau_ = search.tau;
    A2_ = search.A2;
    Q_ = search.Q;
    ns_ = search.ns;
  }

  void MC2MuCapture::fillCaptureTree(const Candidate& c, const PulseSearch& search) {
    fillSearchBranches(search);
    kind_ = c.kind;
    cry_ = c.cry;
    disk_ = c.disk;
    run_ = c.run;
    sub_ = c.sub;
    evt_ = c.evt;
    nMi_ = c.nMichel;
    miDt_ = c.michelDt;
    miE_ = c.michelE;
    tNTrk_ = c.nTrk;
    tUx_ = c.ux;
    tUy_ = c.uy;
    tCnext_ = c.cnext;
    tStopE_ = c.stopE;
    tLyCry_ = c.lyCry;
    tLyDist_ = c.lyDist;
    tLyOnTrk_ = c.lyOnTrk;
    for (int s = 0; s < 2; ++s) {
      tLyE_[s] = c.lyE[s];
      tLyA_[s] = c.lyA[s];
      tLyDt_[s] = c.lyDt[s];
    }
    tDelN_ = c.delN;
    tDelE_ = c.delE;
    tDelDt_ = c.delDt;
    captureTree_->Fill();
  }

  // Example waveforms (raw - pedestal): class 0 STOP found (5 x nExamples), 1 THROUGH found, 2 STOP not found.
  // Uses the branches of the capture tree just filled.
  void MC2MuCapture::saveExample(const Candidate& c, const PulseSearch& search) {
    const bool found = search.maxS > exampleS_ && search.A2 > exampleA2_;
    const int cls = c.kind == 0 ? (found ? 0 : 2) : (found ? 1 : -1);
    const int limit = cls == 0 ? 5 * nExamples_ : nExamples_;
    if (cls < 0 || nExamplesSaved_[cls] >= limit) {
      return;
    }
    ++nExamplesSaved_[cls];
    exCls_ = cls;
    exY0_.assign(c.w[0].y.begin(), c.w[0].y.end());
    for (auto& v : exY0_) {
      v *= c.w[0].amp;
    }
    exY1_.assign(c.w[1].y.begin(), c.w[1].y.end());
    for (auto& v : exY1_) {
      v *= c.w[1].amp;
    }
    exT0a_ = c.w[0].t0;
    exT0b_ = c.w[1].t0;
    exA1a_ = search.A1Sipm[0];
    exA1b_ = search.A1Sipm[1];
    exT1a_ = search.t1Sipm[0];
    exT1b_ = search.t1Sipm[1];
    exTref_ = search.t1;
    examplesTree_->Fill();
  }

  // Injection: a pulse of the crossing-muon template T, random delay and amplitude, added to both SiPMs after their
  // first pulse (fitted with the search template), then the same search. Uses the identification branches of the
  // capture tree just filled.
  void MC2MuCapture::fillInjection(const Candidate& c, const PulseTemplate& T, const PulseTemplate& searchTemplate,
                                   TRandom3& rnd) {
    Candidate injected = c;
    const double delay = rnd.Uniform(dtMin_, injMaxDelay_);
    const double a2 = rnd.Uniform(injMinA2_, injMaxA2_);
    for (int s = 0; s < 2; ++s) {
      double bestTime, amplitude;
      fitFirstPulse(c.w[s], searchTemplate, kFirstPulseLo, kFirstPulseHi, bestTime, amplitude);
      for (size_t q = 0; q < injected.w[s].y.size(); ++q) {
        injected.w[s].y[q] += float(a2 / c.w[s].amp * T.at(q * kSamplePeriod - bestTime - delay));
      }
    }
    const PulseSearch search = searchSecondPulse(injected, searchTemplate);
    fillSearchBranches(search);
    kind_ = c.kind == 1 ? 2 : 3;
    injD_ = delay;
    injA2_ = a2;
    injectTree_->Fill();
  }

}

DEFINE_ART_MODULE(mu2e::MC2MuCapture)
