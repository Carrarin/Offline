// calo_tracks_vs_tagger.C -- event-by-event comparison of the art calorimeter tracks (MC2CaloCosmicTrackDump, tree dump/tracks)
// with the mc2tools tagger (tracks.root, tree tracks), on the same run.
// Events are matched by ew (art event number = roottuple currentDTCEventWindow). Crystals: art gives hitCrystal and hitKey (raw id of
// SiPM 0); the tagger gives hitKey of its L sensor. A tagger key is mapped to a crystal through the art (key, crystal) pairs of the run;
// a key never seen in art (the L sensor is the other SiPM of the crystal) is mapped by majority vote in the events where it is the only
// unknown key of a tagger track. A tagger track and an art track of the same event match if they share >= 50 % of the crystals of the
// smaller one (each art track used once).
// Output: printed summary (fractions matched per type, same type, angle between directions) and <out>.png.
// Usage: root -l -b -q 'Offline/MC2Reco/scripts/calo_tracks_vs_tagger.C+("calo_tracks_<run>.root","tracks.root","out/cmp_<run>")'
#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <cmath>
#include <cstdio>
#include <map>
#include <set>
#include <vector>
struct T { int type; double dx, dy, dz; std::set<int> cry; bool used = false; };
void calo_tracks_vs_tagger(const char* artFile, const char* tagFile, const char* out) {
  TFile fa(artFile), ft(tagFile); TTree* ta = (TTree*)fa.Get("dump/tracks"); TTree* tt = (TTree*)ft.Get("tracks");
  if (!ta || !tt) { printf("ERROR: missing trees\n"); return; }
  Long64_t ew; int type; double dx, dy, dz; std::vector<int>* key = nullptr; std::vector<int>* cry = nullptr;
  std::map<int, int> k2c;
  ta->SetBranchAddress("ew", &ew); ta->SetBranchAddress("type", &type); ta->SetBranchAddress("dx", &dx); ta->SetBranchAddress("dy", &dy); ta->SetBranchAddress("dz", &dz);
  ta->SetBranchAddress("hitKey", &key); ta->SetBranchAddress("hitCrystal", &cry);
  std::multimap<Long64_t, T> A;
  for (Long64_t i = 0; i < ta->GetEntries(); ++i) { ta->GetEntry(i); T x{type, dx, dy, dz, {}}; for (size_t j = 0; j < cry->size(); ++j) { x.cry.insert((*cry)[j]); k2c[(*key)[j]] = (*cry)[j]; } A.insert({ew, x}); }
  std::vector<int>* tkey = nullptr; tt->SetBranchAddress("ew", &ew); tt->SetBranchAddress("type", &type); tt->SetBranchAddress("dx", &dx); tt->SetBranchAddress("dy", &dy); tt->SetBranchAddress("dz", &dz); tt->SetBranchAddress("hitKey", &tkey);
  std::map<int, std::map<int, int>> vote;
  for (Long64_t i = 0; i < tt->GetEntries(); ++i) { tt->GetEntry(i); auto r = A.equal_range(ew); if (r.first == r.second) continue;
    std::set<int> known; std::vector<int> unk; for (int k : *tkey) { auto it = k2c.find(k); if (it != k2c.end()) known.insert(it->second); else unk.push_back(k); }
    if (unk.size() != 1) continue;
    for (auto it = r.first; it != r.second; ++it) for (int c : it->second.cry) if (!known.count(c)) vote[unk[0]][c]++; }
  int learned = 0; for (auto& kv : vote) { int best = -1, nb = 0; for (auto& cv : kv.second) if (cv.second > nb) { nb = cv.second; best = cv.first; } if (nb >= 3) { k2c[kv.first] = best; ++learned; } }
  TH1D hAng("hAng", "matched tracks;angle between the art and tagger directions [deg];tracks", 100, 0, 5);
  TH2D hType("hType", ";tagger type;art type", 3, 0.5, 3.5, 3, 0.5, 3.5);
  long nTag = 0, nMatch = 0, nNoEv = 0, nUnmapped = 0; long nTagT[4] = {0}, nMatT[4] = {0};
  for (Long64_t i = 0; i < tt->GetEntries(); ++i) { tt->GetEntry(i); ++nTag; ++nTagT[type];
    std::set<int> c; for (int k : *tkey) { auto it = k2c.find(k); if (it != k2c.end()) c.insert(it->second); else ++nUnmapped; }
    auto r = A.equal_range(ew); if (r.first == r.second) { ++nNoEv; continue; }
    T* best = nullptr; double bf = 0;
    for (auto it = r.first; it != r.second; ++it) { if (it->second.used) continue; int sh = 0; for (int x : c) sh += it->second.cry.count(x);
      const double fr = double(sh) / std::max<size_t>(1, std::min(c.size(), it->second.cry.size())); if (fr > bf) { bf = fr; best = &it->second; } }
    if (!best || bf < 0.5) continue;
    best->used = true; ++nMatch; ++nMatT[type]; hType.Fill(type, best->type);
    const double cs = std::fabs(dx * best->dx + dy * best->dy + dz * best->dz); hAng.Fill(std::acos(std::min(1., cs)) * 180. / M_PI); }
  long nArt = A.size(), nArtUsed = 0; for (auto& kv : A) nArtUsed += kv.second.used;
  printf("tagger tracks %ld (type 1/2/3: %ld %ld %ld); matched to an art track %ld (%.2f %%) [type 1/2/3: %.2f %.2f %.2f %%]\n", nTag, nTagT[1], nTagT[2], nTagT[3], nMatch, 100. * nMatch / nTag,
         100. * nMatT[1] / nTagT[1], 100. * nMatT[2] / nTagT[2], 100. * nMatT[3] / nTagT[3]);
  printf("art tracks %ld, matched %ld (%.2f %%); tagger tracks in events without art tracks %ld; tagger keys learned %d, unmapped hits %ld\n", nArt, nArtUsed, 100. * nArtUsed / nArt, nNoEv, learned, nUnmapped);
  printf("same type for matched tracks: %.2f %%; angle between directions: mean %.3f deg, rms %.3f deg, > 1 deg %.2f %%\n",
         100. * (hType.GetBinContent(1, 1) + hType.GetBinContent(2, 2) + hType.GetBinContent(3, 3)) / std::max(1., hType.GetEntries()), hAng.GetMean(), hAng.GetRMS(),
         100. * hAng.Integral(hAng.FindBin(1.0001), hAng.GetNbinsX() + 1) / std::max(1., hAng.GetEntries()));
  TCanvas cv("cv", "", 1200, 500); cv.Divide(2, 1); cv.cd(1)->SetLogy(); hAng.Draw(); cv.cd(2); hType.Draw("text colz"); cv.SaveAs(Form("%s.png", out));
}
