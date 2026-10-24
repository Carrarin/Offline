// crv_calo_flatten.C: EventNtuple (calorimeter + CRV) -> two flat csv tables for the coincidence analysis
// (crv_calo_coinc.py, crv_calo_spatial.py).
//
// Usage (any ROOT: the branches are read by name, no EventNtuple dictionaries needed):
//
//   root -l -b -q 'Offline/MC2Reco/scripts/crv_calo_flatten.C("nts.root","out/prefix")'
//
// Arguments: in = EventNtuple file (tree nt/ntuple with crvcoincs, calohits, evtinfo); out = prefix of the csv files.
//
// Output (positions in detector coordinates, mm; times in ns):
//   <out>_crv.csv : ev,sector,x,y,z,time,timeStart,timeEnd,PEs,nHits,nLayers,angle,posCalc
//                   one row per CRV coincidence (sector = sectorType, posCalc = hitPosAndTimeCalculated)
//   <out>_calo.csv: ev,crystal,time,eDep,x,y,z
//                   one row per calorimeter hit (crystal = crystalId, x y z = crystal position), only in the events
//                   with at least one CRV coincidence
// Prints the number of events and of events with CRV coincidences.

#include <TFile.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TTreeReaderValue.h>

#include <cstdio>

void crv_calo_flatten(const char* in, const char* out) {
  TFile inputFile(in);
  TTreeReader reader("nt/ntuple", &inputFile);

  // CRV coincidences
  TTreeReaderArray<int> crvSector(reader, "crvcoincs.sectorType");
  TTreeReaderArray<float> crvX(reader, "crvcoincs.pos.fCoordinates.fX");
  TTreeReaderArray<float> crvY(reader, "crvcoincs.pos.fCoordinates.fY");
  TTreeReaderArray<float> crvZ(reader, "crvcoincs.pos.fCoordinates.fZ");
  TTreeReaderArray<float> crvTime(reader, "crvcoincs.time");
  TTreeReaderArray<float> crvTimeStart(reader, "crvcoincs.timeStart");
  TTreeReaderArray<float> crvTimeEnd(reader, "crvcoincs.timeEnd");
  TTreeReaderArray<float> crvPEs(reader, "crvcoincs.PEs");
  TTreeReaderArray<float> crvAngle(reader, "crvcoincs.angle");
  TTreeReaderArray<int> crvNHits(reader, "crvcoincs.nHits");
  TTreeReaderArray<int> crvNLayers(reader, "crvcoincs.nLayers");
  TTreeReaderArray<bool> crvPosCalculated(reader, "crvcoincs.hitPosAndTimeCalculated");

  // calorimeter hits
  TTreeReaderArray<int> caloCrystal(reader, "calohits.crystalId_");
  TTreeReaderArray<float> caloTime(reader, "calohits.time_");
  TTreeReaderArray<float> caloEDep(reader, "calohits.eDep_");
  TTreeReaderArray<float> caloX(reader, "calohits.crystalPos_.fCoordinates.fX");
  TTreeReaderArray<float> caloY(reader, "calohits.crystalPos_.fCoordinates.fY");
  TTreeReaderArray<float> caloZ(reader, "calohits.crystalPos_.fCoordinates.fZ");

  TTreeReaderValue<int> eventNumber(reader, "evtinfo.event");

  FILE* crvCsv = fopen(Form("%s_crv.csv", out), "w");
  FILE* caloCsv = fopen(Form("%s_calo.csv", out), "w");
  fprintf(crvCsv, "ev,sector,x,y,z,time,timeStart,timeEnd,PEs,nHits,nLayers,angle,posCalc\n");
  fprintf(caloCsv, "ev,crystal,time,eDep,x,y,z\n");

  long nEvents = 0;
  long nEventsWithCrv = 0;
  while (reader.Next()) {
    ++nEvents;
    if (crvSector.GetSize() == 0) {
      continue;
    }
    ++nEventsWithCrv;
    for (size_t i = 0; i < crvSector.GetSize(); ++i) {
      fprintf(crvCsv, "%d,%d,%.1f,%.1f,%.1f,%.3f,%.3f,%.3f,%.1f,%d,%d,%.3f,%d\n", *eventNumber, crvSector[i], crvX[i],
              crvY[i], crvZ[i], crvTime[i], crvTimeStart[i], crvTimeEnd[i], crvPEs[i], crvNHits[i], crvNLayers[i],
              crvAngle[i], (int)crvPosCalculated[i]);
    }
    for (size_t i = 0; i < caloCrystal.GetSize(); ++i) {
      fprintf(caloCsv, "%d,%d,%.3f,%.3f,%.1f,%.1f,%.1f\n", *eventNumber, caloCrystal[i], caloTime[i], caloEDep[i],
              caloX[i], caloY[i], caloZ[i]);
    }
  }

  fclose(crvCsv);
  fclose(caloCsv);
  printf("events %ld, with CRV coincidences %ld\n", nEvents, nEventsWithCrv);
}
