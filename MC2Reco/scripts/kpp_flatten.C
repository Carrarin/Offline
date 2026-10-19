// kpp_flatten.C -- EventNtuple of the KPP combined run (mc2_kpp_trk_calo_ntuple.fcl) -> flat text tables.
//   <out>_trk.csv : ev, itrk, status, goodfit, nhits, nactive, nplanes, ndof, chisq, fitcon, nseg
//   <out>_seg.csv : ev, itrk, sid, inbounds, x, y, z, px, py, pz, time   (detector frame, mm, MeV/c, ns; sid: Offline SurfaceId)
//   <out>_calo.csv: ev, crystal, time, eDep, x, y, z (crystal front face)  [events with >= 1 track]
//   <out>_crv.csv : ev, sector, x, y, z, time, PEs, nLayers, posCalc        [events with >= 1 track]
// Usage: root -l -b -q 'Offline/MC2Reco/scripts/kpp_flatten.C("nts.root","out/prefix")'   (muse environment, EventNtuple dictionaries)
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
  TFile f(in); TTreeReader r("nt/ntuple", &f);
  TTreeReaderValue<mu2e::EventInfo> ei(r, "evtinfo");
  TTreeReaderValue<std::vector<mu2e::TrkInfo>> trk(r, "trk");
  TTreeReaderValue<std::vector<std::vector<mu2e::TrkSegInfo>>> seg(r, "trksegs");
  TTreeReaderValue<std::vector<mu2e::CaloHitInfo>> ch(r, "calohits");
  TTreeReaderValue<std::vector<mu2e::CrvHitInfoReco>> cc(r, "crvcoincs");
  FILE* ft = fopen(Form("%s_trk.csv", out), "w"); FILE* fs = fopen(Form("%s_seg.csv", out), "w");
  FILE* fh = fopen(Form("%s_calo.csv", out), "w"); FILE* fc = fopen(Form("%s_crv.csv", out), "w");
  fprintf(ft, "ev,itrk,status,goodfit,nhits,nactive,nplanes,ndof,chisq,fitcon,nseg\n");
  fprintf(fs, "ev,itrk,sid,inbounds,x,y,z,px,py,pz,time\n");
  fprintf(fh, "ev,crystal,time,eDep,x,y,z\n"); fprintf(fc, "ev,sector,x,y,z,time,PEs,nLayers,posCalc\n");
  long n = 0, nt = 0, ntrk = 0;
  while (r.Next()) { ++n; const int ev = ei->event; if (trk->empty()) continue; ++nt;
    for (size_t i = 0; i < trk->size(); ++i) { const auto& t = (*trk)[i]; ++ntrk;
      fprintf(ft, "%d,%zu,%d,%d,%d,%d,%d,%d,%.4f,%.6g,%d\n", ev, i, t.status, t.goodfit, t.nhits, t.nactive, t.nplanes, t.ndof, t.chisq, t.fitcon, t.nseg);
      if (i < seg->size()) for (const auto& s : (*seg)[i])
        fprintf(fs, "%d,%zu,%d,%d,%.2f,%.2f,%.2f,%.3f,%.3f,%.3f,%.3f\n", ev, i, s.sid, (int)s.inbounds, s.pos.x(), s.pos.y(), s.pos.z(), s.mom.x(), s.mom.y(), s.mom.z(), s.time); }
    for (const auto& h : *ch) fprintf(fh, "%d,%d,%.3f,%.3f,%.1f,%.1f,%.1f\n", ev, h.crystalId_, h.time_, h.eDep_, h.crystalPos_.x(), h.crystalPos_.y(), h.crystalPos_.z());
    for (const auto& c : *cc) fprintf(fc, "%d,%d,%.1f,%.1f,%.1f,%.3f,%.1f,%d,%d\n", ev, c.sectorType, c.pos.x(), c.pos.y(), c.pos.z(), c.time, c.PEs, c.nLayers, (int)c.hitPosAndTimeCalculated); }
  fclose(ft); fclose(fs); fclose(fh); fclose(fc); printf("events %ld, with tracks %ld, tracks %ld\n", n, nt, ntrk);
}
