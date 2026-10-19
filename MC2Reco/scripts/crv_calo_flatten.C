// crv_calo_flatten.C -- EventNtuple (calo + CRV) -> two flat text tables for the coincidence analysis.
//   <out>_crv.csv : ev, sector, pos x y z (detector coordinates, mm), time, timeStart, timeEnd, PEs, nHits, nLayers, angle, posCalc
//   <out>_calo.csv: ev, crystalId, time, eDep, crystal x y z (detector coordinates, mm)   [only events with >= 1 CRV coincidence]
// Usage: root -l -b -q 'Offline/MC2Reco/scripts/crv_calo_flatten.C("nts.root","out/prefix")'
#include <TFile.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TTreeReaderValue.h>
#include <cstdio>
void crv_calo_flatten(const char* in, const char* out) {
  TFile f(in); TTreeReader r("nt/ntuple", &f);
  TTreeReaderArray<int> cs(r, "crvcoincs.sectorType"); TTreeReaderArray<float> cx(r, "crvcoincs.pos.fCoordinates.fX"), cy(r, "crvcoincs.pos.fCoordinates.fY"), cz(r, "crvcoincs.pos.fCoordinates.fZ");
  TTreeReaderArray<float> ct(r, "crvcoincs.time"), ct0(r, "crvcoincs.timeStart"), ct1(r, "crvcoincs.timeEnd"), cpe(r, "crvcoincs.PEs"), cang(r, "crvcoincs.angle");
  TTreeReaderArray<int> cnh(r, "crvcoincs.nHits"), cnl(r, "crvcoincs.nLayers"); TTreeReaderArray<bool> cpc(r, "crvcoincs.hitPosAndTimeCalculated");
  TTreeReaderArray<int> hid(r, "calohits.crystalId_"); TTreeReaderArray<float> ht(r, "calohits.time_"), he(r, "calohits.eDep_");
  TTreeReaderArray<float> hx(r, "calohits.crystalPos_.fCoordinates.fX"), hy(r, "calohits.crystalPos_.fCoordinates.fY"), hz(r, "calohits.crystalPos_.fCoordinates.fZ");
  TTreeReaderValue<int> evn(r, "evtinfo.event");
  FILE* fc = fopen(Form("%s_crv.csv", out), "w"); FILE* fh = fopen(Form("%s_calo.csv", out), "w");
  fprintf(fc, "ev,sector,x,y,z,time,timeStart,timeEnd,PEs,nHits,nLayers,angle,posCalc\n"); fprintf(fh, "ev,crystal,time,eDep,x,y,z\n");
  long n = 0, nc = 0;
  while (r.Next()) { ++n; if (cs.GetSize() == 0) continue; ++nc;
    for (size_t i = 0; i < cs.GetSize(); ++i) fprintf(fc, "%d,%d,%.1f,%.1f,%.1f,%.3f,%.3f,%.3f,%.1f,%d,%d,%.3f,%d\n", *evn, cs[i], cx[i], cy[i], cz[i], ct[i], ct0[i], ct1[i], cpe[i], cnh[i], cnl[i], cang[i], (int)cpc[i]);
    for (size_t i = 0; i < hid.GetSize(); ++i) fprintf(fh, "%d,%d,%.3f,%.3f,%.1f,%.1f,%.1f\n", *evn, hid[i], ht[i], he[i], hx[i], hy[i], hz[i]); }
  fclose(fc); fclose(fh); printf("events %ld, with CRV coincidences %ld\n", n, nc);
}
