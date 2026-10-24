// tagger_offsets_to_tcorr.C: time offsets of the mc2tools tagger -> "SiPMID tcorr[ns]" file of MC2CaloRecoDigiMaker.
//
// Input: channel_time_offsets.csv of mc2tools (one header line, then board,chan,offset in samples, ...). The raw key of
// a channel is 20 * board + chan, the Offline raw SiPM id. The raw -> offline map is the hKeyOffline histogram written
// by MC2RawCheck (Offline CalDAQMap: bin k+1 = offline id + 1 of raw key k; 0 = no offline id).
// Same sign convention: the tagger uses t - offset, MC2CaloRecoDigiMaker t - tcorr.
//
// Usage (any ROOT):
//
//   root -l -b -q 'Offline/MC2Reco/scripts/tagger_offsets_to_tcorr.C("channel_time_offsets.csv","raw_check_<run>.root",
//                  "<run>_sipm_tcorr.txt")'
//
// Arguments: csv = tagger offsets; rawCheck = output of MC2RawCheck (histogram rawCheck/MC2RawCheck/hKeyOffline or
//            MC2RawCheck/hKeyOffline); out = file to write; nsPerSample = sample period (default 5 ns).
//
// Output: <out>, a comment line with the provenance, then one line "offlineSiPMId tcorr_ns" (4 decimals) per channel
// with an offline id, in the order of the csv. Prints the number of SiPMs written and of keys without an offline id.

#include <TFile.h>
#include <TH1.h>

#include <cstdio>
#include <fstream>
#include <sstream>
#include <string>

void tagger_offsets_to_tcorr(const char* csv, const char* rawCheck, const char* out, double nsPerSample = 5.) {
  // raw -> offline map, in either of the two directories where MC2RawCheck may have written it
  TFile rawCheckFile(rawCheck);
  TH1* keyToOffline = (TH1*)rawCheckFile.Get("rawCheck/MC2RawCheck/hKeyOffline");
  if (!keyToOffline) {
    keyToOffline = (TH1*)rawCheckFile.Get("MC2RawCheck/hKeyOffline");
  }
  if (!keyToOffline) {
    printf("ERROR: no hKeyOffline in %s\n", rawCheck);
    return;
  }

  // skip the header line of the csv
  std::ifstream offsetsCsv(csv);
  std::string line;
  std::getline(offsetsCsv, line);

  FILE* tcorrFile = fopen(out, "w");
  fprintf(tcorrFile, "# SiPMID tcorr[ns]  from %s (tagger offsets x %.1f ns), raw->offline map of %s\n", csv,
          nsPerSample, rawCheck);

  int nWritten = 0;
  int nMissing = 0;
  while (std::getline(offsetsCsv, line)) {
    // board, channel, offset in samples (first three fields)
    std::stringstream fields(line);
    std::string field;
    std::getline(fields, field, ',');
    const int board = std::stoi(field);
    std::getline(fields, field, ',');
    const int channel = std::stoi(field);
    std::getline(fields, field, ',');
    const double offsetSamples = std::stod(field);

    const int rawKey = 20 * board + channel;
    const int offlineId = int(keyToOffline->GetBinContent(rawKey + 1)) - 1;
    if (offlineId < 0) {
      ++nMissing;
      continue;
    }
    fprintf(tcorrFile, "%d %.4f\n", offlineId, offsetSamples * nsPerSample);
    ++nWritten;
  }

  fclose(tcorrFile);
  printf("written %s: %d SiPMs, %d keys without an offline id\n", out, nWritten, nMissing);
}
