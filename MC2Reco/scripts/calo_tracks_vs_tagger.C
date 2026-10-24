// calo_tracks_vs_tagger.C: event-by-event comparison of the art calorimeter tracks with the mc2tools tagger, on the
// same run.
//
// Inputs: the art tracks (MC2CaloCosmicTrackDump, tree dump/tracks) and the tagger tracks (tracks.root, tree tracks).
// Branches used in both: ew, type (1, 2, 3), dx, dy, dz (direction), hitKey; art also has hitCrystal.
//
// Method:
//   - Events are matched by ew (art event number = roottuple currentDTCEventWindow).
//   - Crystals: art gives hitCrystal and hitKey (raw id of SiPM 0); the tagger gives the hitKey of its L sensor. A
//     tagger key is mapped to a crystal through the art (key, crystal) pairs of the run. A key never seen in art (the
//     L sensor is the other SiPM of the crystal) is mapped by majority vote (at least 3 votes) in the events where it
//     is the only unknown key of a tagger track: each vote is a crystal of the art tracks of the event not already
//     explained by the known keys.
//   - A tagger track and an art track of the same event match if they share >= 50 % of the crystals of the smaller
//     one (the art track with the largest shared fraction; each art track used once).
//
// Usage (any ROOT, compiled):
//
//   root -l -b -q 'Offline/MC2Reco/scripts/calo_tracks_vs_tagger.C+("calo_tracks_<run>.root","tracks.root",
//                  "out/cmp_<run>")'
//
// Arguments: artFile = file with dump/tracks; tagFile = tagger file with tracks; out = prefix of the figure.
//
// Output: printed summary (fractions matched per type, same type, angle between the directions) and <out>.png
// (angle between the directions of the matched tracks; tagger type vs art type).

#include <TCanvas.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TTree.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <map>
#include <set>
#include <vector>

// One art track: type, direction, crystals; used = already matched to a tagger track
struct ArtTrack {
  int type;
  double dx, dy, dz;
  std::set<int> crystals;
  bool used = false;
};
using ArtTracksByEvent = std::multimap<Long64_t, ArtTrack>;

// Branches of a track tree (the same names in art and in the tagger)
struct TrackBranches {
  Long64_t ew;
  int type;
  double dx, dy, dz;
  std::vector<int>* hitKey = nullptr;
};

void connectTrackBranches(TTree* tree, TrackBranches& branches) {
  tree->SetBranchAddress("ew", &branches.ew);
  tree->SetBranchAddress("type", &branches.type);
  tree->SetBranchAddress("dx", &branches.dx);
  tree->SetBranchAddress("dy", &branches.dy);
  tree->SetBranchAddress("dz", &branches.dz);
  tree->SetBranchAddress("hitKey", &branches.hitKey);
}

// All the art tracks by event, and the (key -> crystal) map seen in art
void readArtTracks(TTree* artTree, ArtTracksByEvent& artTracks, std::map<int, int>& keyToCrystal) {
  TrackBranches art;
  std::vector<int>* hitCrystal = nullptr;
  connectTrackBranches(artTree, art);
  artTree->SetBranchAddress("hitCrystal", &hitCrystal);
  for (Long64_t i = 0; i < artTree->GetEntries(); ++i) {
    artTree->GetEntry(i);
    ArtTrack track{art.type, art.dx, art.dy, art.dz, {}};
    for (size_t j = 0; j < hitCrystal->size(); ++j) {
      track.crystals.insert((*hitCrystal)[j]);
      keyToCrystal[(*art.hitKey)[j]] = (*hitCrystal)[j];
    }
    artTracks.insert({art.ew, track});
  }
  artTree->ResetBranchAddresses();   // the branch variables are local to this function
}

// Map the tagger keys never seen in art by majority vote; returns the number of keys learned
int learnUnknownKeys(TTree* taggerTree, TrackBranches& tagger, const ArtTracksByEvent& artTracks,
                     std::map<int, int>& keyToCrystal) {
  std::map<int, std::map<int, int>> votes;   // unknown key -> crystal -> votes
  for (Long64_t i = 0; i < taggerTree->GetEntries(); ++i) {
    taggerTree->GetEntry(i);
    auto eventTracks = artTracks.equal_range(tagger.ew);
    if (eventTracks.first == eventTracks.second) {
      continue;
    }
    std::set<int> knownCrystals;
    std::vector<int> unknownKeys;
    for (int key : *tagger.hitKey) {
      auto it = keyToCrystal.find(key);
      if (it != keyToCrystal.end()) {
        knownCrystals.insert(it->second);
      } else {
        unknownKeys.push_back(key);
      }
    }
    if (unknownKeys.size() != 1) {
      continue;
    }
    for (auto it = eventTracks.first; it != eventTracks.second; ++it) {
      for (int crystal : it->second.crystals) {
        if (!knownCrystals.count(crystal)) {
          votes[unknownKeys[0]][crystal]++;
        }
      }
    }
  }

  int learned = 0;
  for (auto& keyVotes : votes) {
    int bestCrystal = -1;
    int bestVotes = 0;
    for (auto& crystalVotes : keyVotes.second) {
      if (crystalVotes.second > bestVotes) {
        bestVotes = crystalVotes.second;
        bestCrystal = crystalVotes.first;
      }
    }
    if (bestVotes >= 3) {
      keyToCrystal[keyVotes.first] = bestCrystal;
      ++learned;
    }
  }
  return learned;
}

void calo_tracks_vs_tagger(const char* artFile, const char* tagFile, const char* out) {
  TFile artInput(artFile), taggerInput(tagFile);
  TTree* artTree = (TTree*)artInput.Get("dump/tracks");
  TTree* taggerTree = (TTree*)taggerInput.Get("tracks");
  if (!artTree || !taggerTree) {
    printf("ERROR: missing trees\n");
    return;
  }

  // 1. art tracks and key -> crystal map; 2. keys of the tagger unknown to art
  ArtTracksByEvent artTracks;
  std::map<int, int> keyToCrystal;
  readArtTracks(artTree, artTracks, keyToCrystal);
  TrackBranches tagger;
  connectTrackBranches(taggerTree, tagger);
  const int learned = learnUnknownKeys(taggerTree, tagger, artTracks, keyToCrystal);

  // 3. match every tagger track to the art track of the same event with the largest shared fraction of crystals
  TH1D hAng("hAng", "matched tracks;angle between the art and tagger directions [deg];tracks", 100, 0, 5);
  TH2D hType("hType", ";tagger type;art type", 3, 0.5, 3.5, 3, 0.5, 3.5);
  long nTag = 0, nMatch = 0, nNoEv = 0, nUnmapped = 0;
  long nTagT[4] = {0}, nMatT[4] = {0};   // per type
  for (Long64_t i = 0; i < taggerTree->GetEntries(); ++i) {
    taggerTree->GetEntry(i);
    ++nTag;
    ++nTagT[tagger.type];

    std::set<int> taggerCrystals;
    for (int key : *tagger.hitKey) {
      auto it = keyToCrystal.find(key);
      if (it != keyToCrystal.end()) {
        taggerCrystals.insert(it->second);
      } else {
        ++nUnmapped;
      }
    }

    auto eventTracks = artTracks.equal_range(tagger.ew);
    if (eventTracks.first == eventTracks.second) {
      ++nNoEv;
      continue;
    }

    ArtTrack* best = nullptr;
    double bestFraction = 0;
    for (auto it = eventTracks.first; it != eventTracks.second; ++it) {
      if (it->second.used) {
        continue;
      }
      int shared = 0;
      for (int crystal : taggerCrystals) {
        shared += it->second.crystals.count(crystal);
      }
      const double fraction =
          double(shared) / std::max<size_t>(1, std::min(taggerCrystals.size(), it->second.crystals.size()));
      if (fraction > bestFraction) {
        bestFraction = fraction;
        best = &it->second;
      }
    }
    if (!best || bestFraction < 0.5) {
      continue;
    }

    best->used = true;
    ++nMatch;
    ++nMatT[tagger.type];
    hType.Fill(tagger.type, best->type);
    const double cosAngle = std::fabs(tagger.dx * best->dx + tagger.dy * best->dy + tagger.dz * best->dz);
    hAng.Fill(std::acos(std::min(1., cosAngle)) * 180. / M_PI);
  }

  // 4. summary
  long nArt = artTracks.size(), nArtUsed = 0;
  for (auto& eventTrack : artTracks) {
    nArtUsed += eventTrack.second.used;
  }
  printf("tagger tracks %ld (type 1/2/3: %ld %ld %ld); matched to an art track %ld (%.2f %%) "
         "[type 1/2/3: %.2f %.2f %.2f %%]\n",
         nTag, nTagT[1], nTagT[2], nTagT[3], nMatch, 100. * nMatch / nTag, 100. * nMatT[1] / nTagT[1],
         100. * nMatT[2] / nTagT[2], 100. * nMatT[3] / nTagT[3]);
  printf("art tracks %ld, matched %ld (%.2f %%); tagger tracks in events without art tracks %ld; "
         "tagger keys learned %d, unmapped hits %ld\n",
         nArt, nArtUsed, 100. * nArtUsed / nArt, nNoEv, learned, nUnmapped);
  printf("same type for matched tracks: %.2f %%; "
         "angle between directions: mean %.3f deg, rms %.3f deg, > 1 deg %.2f %%\n",
         100. * (hType.GetBinContent(1, 1) + hType.GetBinContent(2, 2) + hType.GetBinContent(3, 3)) /
             std::max(1., hType.GetEntries()),
         hAng.GetMean(), hAng.GetRMS(),
         100. * hAng.Integral(hAng.FindBin(1.0001), hAng.GetNbinsX() + 1) / std::max(1., hAng.GetEntries()));

  // 5. figure: angle (log scale) and type vs type
  TCanvas cv("cv", "", 1200, 500);
  cv.Divide(2, 1);
  cv.cd(1)->SetLogy();
  hAng.Draw();
  cv.cd(2);
  hType.Draw("text colz");
  cv.SaveAs(Form("%s.png", out));
}
