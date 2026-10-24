//
// MC2CaloCosmicTrackFinder: cosmic-ray tracks in the calorimeter alone.
//
// It is the algorithm of the mc2tools tagger (macros/cosmic_track_tagger.C) run on CaloHits, so that it works on any
// reconstruction (MC2 or Offline, data or Monte Carlo).
//
// Input:  CaloHitCollection (caloHitCollection); CalCalib conditions (ADC2MeV, to turn the energies back into ADC
//         amplitudes); Calorimeter geometry.
// Output: MC2CaloCosmicTrackCollection (see MC2Reco/inc/MC2CaloCosmicTrack.hh), one entry per track found.
//
// Steps, for every event:
//
//   1. crystal hits: CaloHits with both SiPMs (requireBothSiPMs) and the amplitude of each SiPM >= minSiPMAmp ADC
//      counts. Amplitude = CaloRecoDigi energy / CalCalib ADC2MeV of that SiPM, i.e. the amplitude of the pulse fit;
//      with rawAmplitude: largest sample of the CaloDigi waveform - mean of its first 4 samples, independent of the
//      fit.
//      Position = crystal centre in the tracker frame: x, y of the crystal,
//      z = crystal front face + depthInCrystal + diskZOffset[disk].
//      Hits of the laser PIN diodes (crystal id beyond the calorimeter, read out in the KPP runs) are skipped.
//   2. time clusters: hits sorted in time, a new cluster starts after a gap > maxGap ns; only clusters with
//      >= minCrystals hits are kept.
//   3. geometry, for each time cluster:
//        two-disk track (type 3): >= min3DPerDisk crystals in each disk and rms distance from a 3D line < rms3D;
//        one-disk track (type 1 = disk 0, 2 = disk 1): >= minCrystals in one disk and < min3DPerDisk in the other,
//          rms distance from a line in the x-y plane < rms2D, and extent along the line >= minLength.
//      The line is the principal axis of the crystal centres (power iteration from a fixed start vector, as in the
//      tagger); its direction is then turned to point up (uy >= 0).
//
// Defaults = the tagger: 150 ADC, 4 samples = 20 ns, 20 mm, 40 mm, 100 mm, crystal centre (200/2 mm).
//
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "cetlib_except/exception.h"
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
#include <utility>
#include <vector>

namespace mu2e {

  namespace {

    // Principal axis of a set of points: fills the centroid and the unit direction, returns the rms distance of the
    // points from the line. The direction comes from a power iteration on the covariance matrix, started from a
    // fixed vector as in the tagger (300 iterations, or until the vector vanishes).
    double fitLine(const std::vector<std::array<double, 3>>& points, double centroid[3], double direction[3]) {
      for (int q = 0; q < 3; ++q) {
        centroid[q] = 0;
      }
      for (const auto& p : points) {
        for (int q = 0; q < 3; ++q) {
          centroid[q] += p[q];
        }
      }
      for (int q = 0; q < 3; ++q) {
        centroid[q] /= points.size();
      }

      double covariance[3][3] = {{0}};
      for (const auto& p : points) {
        for (int a = 0; a < 3; ++a) {
          for (int b = 0; b < 3; ++b) {
            covariance[a][b] += (p[a] - centroid[a]) * (p[b] - centroid[b]);
          }
        }
      }

      direction[0] = 0.31;
      direction[1] = 0.52;
      direction[2] = 0.79;
      for (int iteration = 0; iteration < 300; ++iteration) {
        double product[3] = {0, 0, 0};
        for (int a = 0; a < 3; ++a) {
          for (int b = 0; b < 3; ++b) {
            product[a] += covariance[a][b] * direction[b];
          }
        }
        const double norm = std::sqrt(product[0] * product[0] + product[1] * product[1] + product[2] * product[2]);
        if (norm <= 0) {
          break;
        }
        for (int a = 0; a < 3; ++a) {
          direction[a] = product[a] / norm;
        }
      }

      double sumOfSquares = 0;
      for (const auto& p : points) {
        const double alongLine = (p[0] - centroid[0]) * direction[0] + (p[1] - centroid[1]) * direction[1]
                               + (p[2] - centroid[2]) * direction[2];
        for (int a = 0; a < 3; ++a) {
          const double residual = p[a] - centroid[a] - alongLine * direction[a];
          sumOfSquares += residual * residual;
        }
      }
      return std::sqrt(sumOfSquares / points.size());
    }

    // A crystal hit selected for the track search.
    struct CrystalHit {
      size_t index;   // index in the CaloHitCollection
      int disk;
      double x, y, z;   // crystal centre, tracker frame [mm]
      double t;         // CaloHit time [ns]
      double e;         // CaloHit energy [MeV]
    };

    // A time cluster: the crystal hits [first, second) of the time-sorted vector.
    using TimeCluster = std::pair<int, int>;

    // The line through the crystals of a time cluster and what is needed to make the track from it.
    struct ClusterLine {
      int type = 0;                 // 0 no track, 1 disk 0, 2 disk 1, 3 both disks
      double centroid[3] = {0, 0, 0};
      double direction[3] = {0, 0, 0};
      double rms = 0;               // rms distance of the crystal centres from the line [mm]
      int nPerDisk[2] = {0, 0};     // crystals of the cluster in disk 0 and disk 1
      double timeSumPerDisk[2] = {0, 0};
    };

  }

  class MC2CaloCosmicTrackFinder : public art::EDProducer {
    public:
      struct Config {
        using Name = fhicl::Name;
        using Comment = fhicl::Comment;

        // Input
        fhicl::Atom<art::InputTag> caloHitTag{Name("caloHitCollection"), Comment("calorimeter hits"), "CaloHitMaker"};

        // 1. Crystal hits
        fhicl::Atom<bool> requireBothSiPMs{Name("requireBothSiPMs"), Comment("crystal hit only if both SiPMs fired"),
            true};
        fhicl::Atom<double> minSiPMAmp{Name("minSiPMAmp"), Comment("minimum amplitude of each SiPM [ADC counts]"),
            150.};
        fhicl::Atom<bool> rawAmplitude{Name("rawAmplitude"),
            Comment("cut on the raw waveform amplitude instead of the fit amplitude"), false};
        fhicl::Atom<double> depthInCrystal{Name("depthInCrystal"),
            Comment("z of a hit behind the crystal front face [mm]"), 100.};
        fhicl::Sequence<double> diskZOffset{Name("diskZOffset"),
            Comment("added to the Offline z of disk 0 and 1 [mm]"), std::vector<double>{0., 0.}};

        // 2. Time clusters
        fhicl::Atom<double> maxGap{Name("maxGap"), Comment("time gap that starts a new cluster [ns]"), 20.};
        fhicl::Atom<int> minCrystals{Name("minCrystals"), Comment("crystals of a time cluster / of a one-disk track"),
            4};

        // 3. Geometry
        fhicl::Atom<int> min3DPerDisk{Name("min3DPerDisk"), Comment("crystals per disk for a two-disk track"), 3};
        fhicl::Atom<double> rms2D{Name("rms2D"), Comment("max rms distance from the x-y line [mm]"), 20.};
        fhicl::Atom<double> rms3D{Name("rms3D"), Comment("max rms distance from the 3D line [mm]"), 40.};
        fhicl::Atom<double> minLength{Name("minLength"), Comment("min extent along the line, one-disk tracks [mm]"),
            100.};
      };
      using Parameters = art::EDProducer::Table<Config>;

      explicit MC2CaloCosmicTrackFinder(const Parameters& conf);
      void produce(art::Event& event) override;

    private:
      double sipmAmplitude(const CaloRecoDigi& recoDigi, const CalCalib& calCalib) const;
      std::vector<CrystalHit> selectCrystalHits(const CaloHitCollection& caloHits, const Calorimeter& calorimeter,
                                                const CalCalib& calCalib) const;
      std::vector<TimeCluster> findTimeClusters(const std::vector<CrystalHit>& crystalHits) const;
      ClusterLine fitClusterLine(const std::vector<CrystalHit>& crystalHits, const TimeCluster& cluster) const;
      bool makeTrack(const std::vector<CrystalHit>& crystalHits, const TimeCluster& cluster, ClusterLine& line,
                     int nTimeClusters, const art::ValidHandle<CaloHitCollection>& caloHitHandle,
                     MC2CaloCosmicTrack& track) const;

      // Input
      art::InputTag caloHitTag_;

      // Crystal hits
      bool requireBothSiPMs_;
      double minSiPMAmp_;
      bool rawAmplitude_;
      double depthInCrystal_;
      std::vector<double> diskZOffset_;   // [disk]

      // Time clusters
      double maxGap_;
      int minCrystals_;

      // Geometry
      int min3DPerDisk_;
      double rms2D_;
      double rms3D_;
      double minLength_;

      // Conditions
      ProditionsHandle<CalCalib> calCalibHandle_;
  };

  MC2CaloCosmicTrackFinder::MC2CaloCosmicTrackFinder(const Parameters& conf) :
    art::EDProducer(conf),
    caloHitTag_(conf().caloHitTag()),
    requireBothSiPMs_(conf().requireBothSiPMs()),
    minSiPMAmp_(conf().minSiPMAmp()),
    rawAmplitude_(conf().rawAmplitude()),
    depthInCrystal_(conf().depthInCrystal()),
    diskZOffset_(conf().diskZOffset()),
    maxGap_(conf().maxGap()),
    minCrystals_(conf().minCrystals()),
    min3DPerDisk_(conf().min3DPerDisk()),
    rms2D_(conf().rms2D()),
    rms3D_(conf().rms3D()),
    minLength_(conf().minLength())
  {
    if (diskZOffset_.size() != 2) {
      throw cet::exception("MC2CALOCOSMIC") << "diskZOffset needs two values\n";
    }
    consumes<CaloHitCollection>(caloHitTag_);
    produces<MC2CaloCosmicTrackCollection>();
  }

  void MC2CaloCosmicTrackFinder::produce(art::Event& event) {
    auto tracks = std::make_unique<MC2CaloCosmicTrackCollection>();
    const Calorimeter& calorimeter = *GeomHandle<Calorimeter>();
    const CalCalib& calCalib = calCalibHandle_.get(event.id());
    const auto caloHitHandle = event.getValidHandle<CaloHitCollection>(caloHitTag_);

    std::vector<CrystalHit> crystalHits = selectCrystalHits(*caloHitHandle, calorimeter, calCalib);

    if ((int)crystalHits.size() >= minCrystals_) {
      std::sort(crystalHits.begin(), crystalHits.end(),
                [](const CrystalHit& a, const CrystalHit& b) { return a.t < b.t; });

      const std::vector<TimeCluster> timeClusters = findTimeClusters(crystalHits);

      for (const auto& cluster : timeClusters) {
        ClusterLine line = fitClusterLine(crystalHits, cluster);
        if (!line.type) {
          continue;
        }
        MC2CaloCosmicTrack track;
        if (makeTrack(crystalHits, cluster, line, timeClusters.size(), caloHitHandle, track)) {
          tracks->push_back(std::move(track));
        }
      }
    }

    event.put(std::move(tracks));
  }

  // Amplitude of one SiPM [ADC counts]: fit amplitude (energy / ADC2MeV) or, with rawAmplitude, largest waveform sample
  // minus the mean of the first 4 samples. -1 if it cannot be computed (no calibration, too short a waveform).
  double MC2CaloCosmicTrackFinder::sipmAmplitude(const CaloRecoDigi& recoDigi, const CalCalib& calCalib) const {
    double amplitude = -1;
    if (rawAmplitude_) {
      const auto& waveform = recoDigi.caloDigiPtr()->waveform();
      if (waveform.size() > 4) {
        const double pedestal = 0.25 * (waveform[0] + waveform[1] + waveform[2] + waveform[3]);
        amplitude = *std::max_element(waveform.begin(), waveform.end()) - pedestal;
      }
    } else {
      const double adcToMeV = calCalib.ADC2MeV(recoDigi.SiPMID());
      if (adcToMeV > 0) {
        amplitude = recoDigi.energyDep() / adcToMeV;
      }
    }
    return amplitude;
  }

  // Step 1: the crystal hits used for the track search, in the order of the CaloHitCollection.
  std::vector<CrystalHit> MC2CaloCosmicTrackFinder::selectCrystalHits(const CaloHitCollection& caloHits,
                                                                      const Calorimeter& calorimeter,
                                                                      const CalCalib& calCalib) const {
    std::vector<CrystalHit> crystalHits;
    for (size_t i = 0; i < caloHits.size(); ++i) {
      const auto& hit = caloHits[i];
      if (requireBothSiPMs_ && hit.nSiPMs() < 2) {
        continue;
      }

      // every SiPM of the crystal above threshold
      bool allAboveThreshold = true;
      for (const auto& recoDigi : hit.recoCaloDigis()) {
        if (sipmAmplitude(*recoDigi, calCalib) < minSiPMAmp_) {
          allAboveThreshold = false;
          break;
        }
      }
      if (!allAboveThreshold) {
        continue;
      }

      // PIN diodes of the laser (read out in the KPP runs) have ids beyond the crystals
      if (hit.crystalID() >= int(calorimeter.nCrystals())) {
        continue;
      }
      const auto& crystal = calorimeter.crystal(hit.crystalID());
      const int disk = crystal.diskID();
      if (disk < 0 || disk > 1) {
        continue;
      }

      // crystal position() is the front face: move to depthInCrystal behind it, plus the disk offset
      const auto position = calorimeter.mu2eToTracker(crystal.position());
      crystalHits.push_back({i, disk, position.x(), position.y(), position.z() + depthInCrystal_ + diskZOffset_[disk],
                             hit.time(), hit.energyDep()});
    }
    return crystalHits;
  }

  // Step 2: time clusters of the time-sorted crystal hits. A cluster ends where the gap to the next hit is > maxGap_;
  // only clusters with >= minCrystals_ hits are returned.
  std::vector<TimeCluster> MC2CaloCosmicTrackFinder::findTimeClusters(const std::vector<CrystalHit>& crystalHits)
      const {
    std::vector<TimeCluster> clusters;
    int clusterStart = 0;
    for (size_t h = 1; h <= crystalHits.size(); ++h) {
      if (h == crystalHits.size() || crystalHits[h].t - crystalHits[h - 1].t > maxGap_) {
        if ((int)h - clusterStart >= minCrystals_) {
          clusters.push_back({clusterStart, (int)h});
        }
        clusterStart = h;
      }
    }
    return clusters;
  }

  // Step 3a: track type and line of a time cluster. Two-disk clusters get a 3D line; one-disk clusters a line in the
  // x-y plane of the disk, placed at the z of that disk's crystal centres. type stays 0 if the cluster is not a track.
  ClusterLine MC2CaloCosmicTrackFinder::fitClusterLine(const std::vector<CrystalHit>& crystalHits,
                                                       const TimeCluster& cluster) const {
    ClusterLine line;
    std::vector<std::array<double, 3>> points3D;
    std::vector<std::array<double, 3>> pointsXYPerDisk[2];   // z = 0
    double diskZ[2] = {0, 0};                                 // z of the (last) crystal of each disk
    for (int h = cluster.first; h < cluster.second; ++h) {
      const int disk = crystalHits[h].disk;
      ++line.nPerDisk[disk];
      line.timeSumPerDisk[disk] += crystalHits[h].t;
      diskZ[disk] = crystalHits[h].z;
      points3D.push_back({crystalHits[h].x, crystalHits[h].y, crystalHits[h].z});
      pointsXYPerDisk[disk].push_back({crystalHits[h].x, crystalHits[h].y, 0.});
    }

    const int* n = line.nPerDisk;
    if (n[0] >= min3DPerDisk_ && n[1] >= min3DPerDisk_) {
      line.rms = fitLine(points3D, line.centroid, line.direction);
      if (line.rms < rms3D_) {
        line.type = 3;
      }
    } else {
      const int disk = n[0] >= n[1] ? 0 : 1;
      if (n[disk] >= minCrystals_ && n[1 - disk] < min3DPerDisk_) {
        line.rms = fitLine(pointsXYPerDisk[disk], line.centroid, line.direction);
        line.centroid[2] = diskZ[disk];
        if (line.rms < rms2D_) {
          line.type = disk + 1;
        }
      }
    }
    return line;
  }

  // Step 3b: the track from the line of a cluster (type != 0). Uses the crystals of the cluster (only those of the disk
  // of the line for one-disk tracks): extent along the line, mean time and rms, energy sum. Returns false if a
  // one-disk track is shorter than minLength_.
  bool MC2CaloCosmicTrackFinder::makeTrack(const std::vector<CrystalHit>& crystalHits, const TimeCluster& cluster,
                                           ClusterLine& line, int nTimeClusters,
                                           const art::ValidHandle<CaloHitCollection>& caloHitHandle,
                                           MC2CaloCosmicTrack& track) const {
    const int type = line.type;
    const double* m = line.centroid;
    double* v = line.direction;

    double sMin = 1e9;
    double sMax = -1e9;
    double meanTime = 0;
    double energySum = 0;
    int nOnTrack = 0;
    for (int h = cluster.first; h < cluster.second; ++h) {
      if (type < 3 && crystalHits[h].disk != type - 1) {
        continue;
      }
      const double q[3] = {crystalHits[h].x, crystalHits[h].y, type == 3 ? crystalHits[h].z : m[2]};
      const double s = (q[0] - m[0]) * v[0] + (q[1] - m[1]) * v[1] + (q[2] - m[2]) * v[2];
      sMin = std::min(sMin, s);
      sMax = std::max(sMax, s);
      meanTime += crystalHits[h].t;
      energySum += crystalHits[h].e;
      ++nOnTrack;
      track.hits.push_back(art::Ptr<CaloHit>(caloHitHandle, crystalHits[h].index));
    }
    if (type < 3 && sMax - sMin < minLength_) {
      return false;
    }

    meanTime /= nOnTrack;
    double timeSumOfSquares = 0;
    for (int h = cluster.first; h < cluster.second; ++h) {
      if (type == 3 || crystalHits[h].disk == type - 1) {
        timeSumOfSquares += (crystalHits[h].t - meanTime) * (crystalHits[h].t - meanTime);
      }
    }

    // direction pointing up
    if (v[1] < 0) {
      for (int a = 0; a < 3; ++a) {
        v[a] = -v[a];
      }
    }

    const int* n = line.nPerDisk;
    const double* td = line.timeSumPerDisk;
    track.type = type;
    track.n0 = n[0];
    track.n1 = n[1];
    track.nTimeClusters = nTimeClusters;
    track.time = meanTime;
    track.timeRms = std::sqrt(timeSumOfSquares / nOnTrack);
    track.dtDisk = (type == 3) ? td[1] / n[1] - td[0] / n[0] : 0.;
    track.cx = m[0];
    track.cy = m[1];
    track.cz = m[2];
    track.ux = v[0];
    track.uy = v[1];
    track.uz = v[2];
    track.lineRms = line.rms;
    track.length = sMax - sMin;
    track.eSum = energySum;
    return true;
  }

}

DEFINE_ART_MODULE(mu2e::MC2CaloCosmicTrackFinder)
