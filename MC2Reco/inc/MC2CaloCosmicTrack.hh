#ifndef MC2Reco_MC2CaloCosmicTrack_hh
#define MC2Reco_MC2CaloCosmicTrack_hh
//
// MC2CaloCosmicTrack: a cosmic-ray track seen by the calorimeter alone (made by MC2CaloCosmicTrackFinder).
//
// A straight line through the centres of the crystals of one time cluster of CaloHits (both SiPMs above threshold):
//   type 1 = disk 0 only, type 2 = disk 1 only: line in the x-y plane of the disk, at z = the disk's crystal centre;
//   type 3 = both disks: line in 3D.
// Frame: tracker (detector) frame of the Offline geometry, mm, plus the disk z offsets given to the finder (fcl).
// Direction (ux, uy, uz): unit vector, pointing up (uy >= 0).
// Times in ns (CaloHit times, so with the per-SiPM offsets that the reconstruction applied, if any).
//
// This is a persisted data product (classes_def.xml): do not change the data members, their types or their order.
//
#include "canvas/Persistency/Common/Ptr.h"
#include "Offline/RecoDataProducts/inc/CaloHit.hh"
#include <vector>

namespace mu2e {

  struct MC2CaloCosmicTrack {
    int   type = 0;                        // 1 disk 0, 2 disk 1, 3 both
    int   n0 = 0, n1 = 0;                  // crystals of the time cluster in disk 0 and disk 1
    int   nTimeClusters = 0;               // time clusters with >= minCrystals (default 4) crystals in the event
    float time = 0, timeRms = 0;           // mean crystal time of the track and rms of the crystal times around it [ns]
    float dtDisk = 0;                      // mean time of disk 1 - disk 0 [ns] (type 3 only, 0 otherwise)
    float cx = 0, cy = 0, cz = 0;          // a point of the line (centroid of the crystal centres) [mm]
    float ux = 0, uy = 0, uz = 0;          // direction, unit vector, uy >= 0
    float lineRms = 0;                     // rms distance of the crystal centres from the line [mm]
    float length = 0;                      // extent of the crystals along the line [mm]
    float eSum = 0;                        // sum of the crystal energies [MeV]
    std::vector<art::Ptr<CaloHit>> hits;   // the crystals on the track (in the disk of the line for types 1 and 2)
  };

  typedef std::vector<MC2CaloCosmicTrack> MC2CaloCosmicTrackCollection;

}

#endif
