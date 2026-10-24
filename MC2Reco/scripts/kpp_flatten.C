// kpp_flatten.C: EventNtuple of the KPP combined run (mc2_kpp_trk_calo_ntuple.fcl) -> flat csv tables.
//
// Usage (muse environment, after "muse setup": the EventNtuple dictionaries and include path are needed):
//
//   root -l -b -q 'Offline/MC2Reco/scripts/kpp_flatten.C("nts.root","out/prefix")'
//
// Arguments: in = EventNtuple file (tree nt/ntuple with evtinfo, trk, trksegs, calohits, crvcoincs);
//            out = prefix of the csv files.
//
// Output (detector frame: mm, MeV/c, ns):
//   <out>_trk.csv : ev,itrk,status,goodfit,nhits,nactive,nplanes,ndof,chisq,fitcon,nseg     one row per track
//   <out>_seg.csv : ev,itrk,sid,inbounds,x,y,z,px,py,pz,time    one row per track segment (sid: Offline SurfaceId)
//   <out>_calo.csv: ev,crystal,time,eDep,x,y,z    one row per calorimeter hit (x y z: crystal front face)
//   <out>_crv.csv : ev,sector,x,y,z,time,PEs,nLayers,posCalc    one row per CRV coincidence
//   Events without tracks are skipped in all four tables. Prints the number of events, events with tracks, tracks.

#include <TFile.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>

#include <cstdio>
#include <vector>

#include "EventNtuple/inc/TrkInfo.hh"
#include "EventNtuple/inc/TrkSegInfo.hh"
#include "EventNtuple/inc/CaloHitInfo.hh"
#include "EventNtuple/inc/CrvHitInfoReco.hh"
#include "EventNtuple/inc/EventInfo.hh"

void kpp_flatten(const char* in, const char* out) {
  TFile inputFile(in);
  TTreeReader reader("nt/ntuple", &inputFile);
  TTreeReaderValue<mu2e::EventInfo> eventInfo(reader, "evtinfo");
  TTreeReaderValue<std::vector<mu2e::TrkInfo>> tracks(reader, "trk");
  TTreeReaderValue<std::vector<std::vector<mu2e::TrkSegInfo>>> trackSegments(reader, "trksegs");
  TTreeReaderValue<std::vector<mu2e::CaloHitInfo>> caloHits(reader, "calohits");
  TTreeReaderValue<std::vector<mu2e::CrvHitInfoReco>> crvCoincidences(reader, "crvcoincs");

  FILE* trackCsv = fopen(Form("%s_trk.csv", out), "w");
  FILE* segmentCsv = fopen(Form("%s_seg.csv", out), "w");
  FILE* caloCsv = fopen(Form("%s_calo.csv", out), "w");
  FILE* crvCsv = fopen(Form("%s_crv.csv", out), "w");
  fprintf(trackCsv, "ev,itrk,status,goodfit,nhits,nactive,nplanes,ndof,chisq,fitcon,nseg\n");
  fprintf(segmentCsv, "ev,itrk,sid,inbounds,x,y,z,px,py,pz,time\n");
  fprintf(caloCsv, "ev,crystal,time,eDep,x,y,z\n");
  fprintf(crvCsv, "ev,sector,x,y,z,time,PEs,nLayers,posCalc\n");

  long nEvents = 0;
  long nEventsWithTracks = 0;
  long nTracks = 0;
  while (reader.Next()) {
    ++nEvents;
    const int ev = eventInfo->event;
    if (tracks->empty()) {
      continue;
    }
    ++nEventsWithTracks;

    // tracks, and the segments of each track (the trksegs entry with the same index, if any)
    for (size_t i = 0; i < tracks->size(); ++i) {
      const auto& track = (*tracks)[i];
      ++nTracks;
      fprintf(trackCsv, "%d,%zu,%d,%d,%d,%d,%d,%d,%.4f,%.6g,%d\n", ev, i, track.status, track.goodfit, track.nhits,
              track.nactive, track.nplanes, track.ndof, track.chisq, track.fitcon, track.nseg);
      if (i < trackSegments->size()) {
        for (const auto& segment : (*trackSegments)[i]) {
          fprintf(segmentCsv, "%d,%zu,%d,%d,%.2f,%.2f,%.2f,%.3f,%.3f,%.3f,%.3f\n", ev, i, segment.sid,
                  (int)segment.inbounds, segment.pos.x(), segment.pos.y(), segment.pos.z(), segment.mom.x(),
                  segment.mom.y(), segment.mom.z(), segment.time);
        }
      }
    }

    // calorimeter hits and CRV coincidences of the event
    for (const auto& hit : *caloHits) {
      fprintf(caloCsv, "%d,%d,%.3f,%.3f,%.1f,%.1f,%.1f\n", ev, hit.crystalId_, hit.time_, hit.eDep_,
              hit.crystalPos_.x(), hit.crystalPos_.y(), hit.crystalPos_.z());
    }
    for (const auto& coincidence : *crvCoincidences) {
      fprintf(crvCsv, "%d,%d,%.1f,%.1f,%.1f,%.3f,%.1f,%d,%d\n", ev, coincidence.sectorType, coincidence.pos.x(),
              coincidence.pos.y(), coincidence.pos.z(), coincidence.time, coincidence.PEs, coincidence.nLayers,
              (int)coincidence.hitPosAndTimeCalculated);
    }
  }

  fclose(trackCsv);
  fclose(segmentCsv);
  fclose(caloCsv);
  fclose(crvCsv);
  printf("events %ld, with tracks %ld, tracks %ld\n", nEvents, nEventsWithTracks, nTracks);
}
