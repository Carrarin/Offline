// tagger_offsets_to_tcorr.C -- mc2tools channel_time_offsets.csv (board, chan, offset in samples; key = 20 * board + chan, the
// Offline raw SiPM id) -> the "SiPMID tcorr[ns]" file of MC2CaloRecoDigiMaker. The raw -> offline map is the hKeyOffline histogram
// written by MC2RawCheck (Offline CalDAQMap, bin k+1 = offline id + 1 of raw key k).
// Same sign: the tagger uses t - offset, MC2CaloRecoDigiMaker t - tcorr.
// Usage: root -l -b -q 'Offline/MC2Reco/scripts/tagger_offsets_to_tcorr.C("channel_time_offsets.csv","raw_check_<run>.root","<run>_sipm_tcorr.txt")'
#include <TFile.h>
#include <TH1.h>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <string>
void tagger_offsets_to_tcorr(const char* csv, const char* rawCheck, const char* out, double nsPerSample = 5.) {
  TFile f(rawCheck); TH1* m = (TH1*)f.Get("rawCheck/MC2RawCheck/hKeyOffline"); if (!m) m = (TH1*)f.Get("MC2RawCheck/hKeyOffline");
  if (!m) { printf("ERROR: no hKeyOffline in %s\n", rawCheck); return; }
  std::ifstream in(csv); std::string line; std::getline(in, line); FILE* fo = fopen(out, "w");
  fprintf(fo, "# SiPMID tcorr[ns]  from %s (tagger offsets x %.1f ns), raw->offline map of %s\n", csv, nsPerSample, rawCheck);
  int n = 0, nmiss = 0;
  while (std::getline(in, line)) { std::stringstream ss(line); std::string t; int b, c; double o;
    std::getline(ss, t, ','); b = std::stoi(t); std::getline(ss, t, ','); c = std::stoi(t); std::getline(ss, t, ','); o = std::stod(t);
    const int key = 20 * b + c; const int off = int(m->GetBinContent(key + 1)) - 1;
    if (off < 0) { ++nmiss; continue; } fprintf(fo, "%d %.4f\n", off, o * nsPerSample); ++n; }
  fclose(fo); printf("written %s: %d SiPMs, %d keys without an offline id\n", out, n, nmiss);
}
