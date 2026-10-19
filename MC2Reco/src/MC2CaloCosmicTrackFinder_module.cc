// MC2CaloCosmicTrackFinder: cosmic-ray tracks in the calorimeter alone, the algorithm of the mc2tools tagger
// (macros/cosmic_track_tagger.C) on CaloHits, so that it runs on any reconstruction (MC2, Offline, data or Monte Carlo).
//
//   1. crystal hits: CaloHits with both SiPMs (requireBothSiPMs) and the amplitude of each SiPM >= minSiPMAmp ADC counts
//      (amplitude = CaloRecoDigi energy / CalCalib ADC2MeV of that SiPM, i.e. the amplitude of the pulse fit; with rawAmplitude:
//      largest sample of the CaloDigi waveform - mean of its first 4 samples, independent of the fit);
//      position = crystal centre, tracker frame: x, y of the crystal, z = crystal front face + depthInCrystal + diskZOffset[disk];
//   2. time clusters: hits sorted in time, a new cluster after a gap > maxGap ns; clusters with >= minCrystals hits;
//   3. geometry: two-disk track (type 3) if >= min3DPerDisk crystals in each disk and the rms distance from a 3D line < rms3D;
//      one-disk track (type 1, 2) if >= minCrystals in one disk and < min3DPerDisk in the other, rms distance from a line in
//      the x-y plane < rms2D, and extent along the line >= minLength.
// The line is the principal axis of the crystal centres (power iteration from a fixed start vector, as in the tagger).
// Defaults = the tagger: 150 ADC, 4 samples = 20 ns, 20 mm, 40 mm, 100 mm, crystal centre (200/2 mm).
// Product: MC2CaloCosmicTrackCollection.
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/CalorimeterGeom/inc/Calorimeter.hh"
#include "Offline/CaloConditions/inc/CalCalib.hh"
#include "Offline/ProditionsService/inc/ProditionsHandle.hh"
#include "Offline/RecoDataProducts/inc/CaloHit.hh"
#include "Offline/RecoDataProducts/inc/CaloRecoDigi.hh"
#include "Offline/MC2Reco/inc/MC2CaloCosmicTrack.hh"
#include <algorithm>
#include <array>
#include <cmath>
#include <memory>
#include <vector>

namespace mu2e {
  namespace {
    // principal axis of a set of points; returns the rms distance from the line, fills centroid m and direction v
    double lineFit(const std::vector<std::array<double, 3>>& pt, double m[3], double v[3]) {
      for (int q = 0; q < 3; ++q) m[q] = 0;
      for (auto& p : pt) for (int q = 0; q < 3; ++q) m[q] += p[q];
      for (int q = 0; q < 3; ++q) m[q] /= pt.size();
      double cov[3][3] = {{0}};
      for (auto& p : pt) for (int a = 0; a < 3; ++a) for (int b = 0; b < 3; ++b) cov[a][b] += (p[a] - m[a]) * (p[b] - m[b]);
      v[0] = 0.31; v[1] = 0.52; v[2] = 0.79;
      for (int it = 0; it < 300; ++it) {
        double w[3] = {0, 0, 0}; for (int a = 0; a < 3; ++a) for (int b = 0; b < 3; ++b) w[a] += cov[a][b] * v[b];
        const double n = std::sqrt(w[0] * w[0] + w[1] * w[1] + w[2] * w[2]); if (n <= 0) break; for (int a = 0; a < 3; ++a) v[a] = w[a] / n;
      }
      double s = 0;
      for (auto& p : pt) { const double t = (p[0] - m[0]) * v[0] + (p[1] - m[1]) * v[1] + (p[2] - m[2]) * v[2];
        for (int a = 0; a < 3; ++a) { const double r = p[a] - m[a] - t * v[a]; s += r * r; } }
      return std::sqrt(s / pt.size());
    }
    struct CHit { size_t idx; int disk; double x, y, z, t, e; };
  }

  class MC2CaloCosmicTrackFinder : public art::EDProducer {
    public:
      struct Config {
        using Name = fhicl::Name; using Comment = fhicl::Comment;
        fhicl::Atom<art::InputTag> caloHitTag       { Name("caloHitCollection"), Comment("calorimeter hits"), "CaloHitMaker" };
        fhicl::Atom<bool>          requireBothSiPMs { Name("requireBothSiPMs"), Comment("crystal hit only if both SiPMs fired"), true };
        fhicl::Atom<double>        minSiPMAmp       { Name("minSiPMAmp"), Comment("minimum amplitude of each SiPM [ADC counts]"), 150. };
        fhicl::Atom<bool>          rawAmplitude     { Name("rawAmplitude"), Comment("cut on the raw waveform amplitude instead of the fit amplitude"), false };
        fhicl::Atom<double>        maxGap           { Name("maxGap"), Comment("time gap that starts a new cluster [ns]"), 20. };
        fhicl::Atom<int>           minCrystals      { Name("minCrystals"), Comment("crystals of a time cluster / of a one-disk track"), 4 };
        fhicl::Atom<int>           min3DPerDisk     { Name("min3DPerDisk"), Comment("crystals per disk for a two-disk track"), 3 };
        fhicl::Atom<double>        rms2D            { Name("rms2D"), Comment("max rms distance from the x-y line [mm]"), 20. };
        fhicl::Atom<double>        rms3D            { Name("rms3D"), Comment("max rms distance from the 3D line [mm]"), 40. };
        fhicl::Atom<double>        minLength        { Name("minLength"), Comment("min extent along the line, one-disk tracks [mm]"), 100. };
        fhicl::Atom<double>        depthInCrystal   { Name("depthInCrystal"), Comment("z of a hit behind the crystal front face [mm]"), 100. };
        fhicl::Sequence<double>    diskZOffset      { Name("diskZOffset"), Comment("added to the Offline z of disk 0 and 1 [mm]"), std::vector<double>{0., 0.} };
      };
      using Parameters = art::EDProducer::Table<Config>;
      explicit MC2CaloCosmicTrackFinder(const Parameters& conf);
      void produce(art::Event& event) override;
    private:
      art::InputTag hitTag_; bool both_; double minAmp_; bool raw_; double maxGap_; int minCry_, min3D_; double rms2D_, rms3D_, minLen_, depth_;
      std::vector<double> dzOff_; ProditionsHandle<CalCalib> calibH_;
  };

  MC2CaloCosmicTrackFinder::MC2CaloCosmicTrackFinder(const Parameters& conf) : art::EDProducer(conf),
    hitTag_(conf().caloHitTag()), both_(conf().requireBothSiPMs()), minAmp_(conf().minSiPMAmp()), raw_(conf().rawAmplitude()), maxGap_(conf().maxGap()),
    minCry_(conf().minCrystals()), min3D_(conf().min3DPerDisk()), rms2D_(conf().rms2D()), rms3D_(conf().rms3D()), minLen_(conf().minLength()),
    depth_(conf().depthInCrystal()), dzOff_(conf().diskZOffset()) {
    if (dzOff_.size() != 2) throw cet::exception("MC2CALOCOSMIC") << "diskZOffset needs two values\n";
    consumes<CaloHitCollection>(hitTag_); produces<MC2CaloCosmicTrackCollection>();
  }

  void MC2CaloCosmicTrackFinder::produce(art::Event& event) {
    auto out = std::make_unique<MC2CaloCosmicTrackCollection>();
    const Calorimeter& cal = *GeomHandle<Calorimeter>();
    const CalCalib& calib = calibH_.get(event.id());
    const auto hitH = event.getValidHandle<CaloHitCollection>(hitTag_); const auto& hits = *hitH;
    // 1. crystal hits
    std::vector<CHit> ch;
    for (size_t i = 0; i < hits.size(); ++i) {
      const auto& h = hits[i];
      if (both_ && h.nSiPMs() < 2) continue;
      bool ok = true;
      for (const auto& rd : h.recoCaloDigis()) {
        double a = -1;
        if (raw_) { const auto& w = rd->caloDigiPtr()->waveform();
          if (w.size() > 4) { const double ped = 0.25 * (w[0] + w[1] + w[2] + w[3]); a = *std::max_element(w.begin(), w.end()) - ped; } }
        else { const double a2m = calib.ADC2MeV(rd->SiPMID()); if (a2m > 0) a = rd->energyDep() / a2m; }
        if (a < minAmp_) { ok = false; break; } }
      if (!ok) continue;
      const auto& c = cal.crystal(h.crystalID()); const int d = c.diskID(); if (d < 0 || d > 1) continue;
      const auto p = cal.mu2eToTracker(c.position());
      ch.push_back({i, d, p.x(), p.y(), p.z() + depth_ + dzOff_[d], h.time(), h.energyDep()});
    }
    if ((int)ch.size() >= minCry_) {
      std::sort(ch.begin(), ch.end(), [](const CHit& a, const CHit& b) { return a.t < b.t; });
      // 2. time clusters
      std::vector<std::pair<int, int>> cl; int i0 = 0;
      for (size_t h = 1; h <= ch.size(); ++h)
        if (h == ch.size() || ch[h].t - ch[h - 1].t > maxGap_) { if ((int)h - i0 >= minCry_) cl.push_back({i0, (int)h}); i0 = h; }
      // 3. geometry
      for (auto& cc : cl) {
        std::vector<std::array<double, 3>> p3, p2[2]; int n[2] = {0, 0}; double td[2] = {0, 0}, zd[2] = {0, 0};
        for (int h = cc.first; h < cc.second; ++h) { const int d = ch[h].disk; ++n[d]; td[d] += ch[h].t; zd[d] = ch[h].z;
          p3.push_back({ch[h].x, ch[h].y, ch[h].z}); p2[d].push_back({ch[h].x, ch[h].y, 0.}); }
        int type = 0; double m[3], v[3], rms = 0;
        if (n[0] >= min3D_ && n[1] >= min3D_) { rms = lineFit(p3, m, v); if (rms < rms3D_) type = 3; }
        else { const int d = n[0] >= n[1] ? 0 : 1;
          if (n[d] >= minCry_ && n[1 - d] < min3D_) { rms = lineFit(p2[d], m, v); m[2] = zd[d]; if (rms < rms2D_) type = d + 1; } }
        if (!type) continue;
        MC2CaloCosmicTrack trk; double smin = 1e9, smax = -1e9, tm = 0, es = 0; int non = 0;
        for (int h = cc.first; h < cc.second; ++h) {
          if (type < 3 && ch[h].disk != type - 1) continue;
          const double q[3] = {ch[h].x, ch[h].y, type == 3 ? ch[h].z : m[2]};
          const double s = (q[0] - m[0]) * v[0] + (q[1] - m[1]) * v[1] + (q[2] - m[2]) * v[2]; smin = std::min(smin, s); smax = std::max(smax, s);
          tm += ch[h].t; es += ch[h].e; ++non; trk.hits.push_back(art::Ptr<CaloHit>(hitH, ch[h].idx));
        }
        if (type < 3 && smax - smin < minLen_) continue;
        tm /= non; double tr = 0;
        for (int h = cc.first; h < cc.second; ++h) if (type == 3 || ch[h].disk == type - 1) tr += (ch[h].t - tm) * (ch[h].t - tm);
        if (v[1] < 0) for (int a = 0; a < 3; ++a) v[a] = -v[a];   // pointing up
        trk.type = type; trk.n0 = n[0]; trk.n1 = n[1]; trk.nTimeClusters = cl.size(); trk.time = tm; trk.timeRms = std::sqrt(tr / non);
        trk.dtDisk = (type == 3) ? td[1] / n[1] - td[0] / n[0] : 0.;
        trk.cx = m[0]; trk.cy = m[1]; trk.cz = m[2]; trk.ux = v[0]; trk.uy = v[1]; trk.uz = v[2];
        trk.lineRms = rms; trk.length = smax - smin; trk.eSum = es;
        out->push_back(std::move(trk));
      }
    }
    event.put(std::move(out));
  }
}
DEFINE_ART_MODULE(mu2e::MC2CaloCosmicTrackFinder)
