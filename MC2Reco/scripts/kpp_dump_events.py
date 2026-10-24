"""Events of the KPP check of the geometry (3D display, residuals, times), from the EventNtuple of kpp_ntuple.fcl.

Coordinates and times: detector frame, mm, ns.
Per tracker track (KKLine) with at least one segment:
  - the line: point, direction and time of the middle segment; the direction is flipped to point down (y < 0),
    cosmics go down;
  - the straw hits: POCA and track time there (hits with state >= 0);
  - calorimeter hits (MC2 chain) with E > EMIN, not PIN diodes (crystalId >= 1348), not the hot crystal, whose time is
    within +-DTCAL of the track time at the crystal z plus CALOFF (calorimeter-tracker time offset);
  - CRV coincidences whose time minus the track time at the CRV height is within (CRVLO, CRVHI).

Inputs   one or more EventNtuple files (tree EventNtuple/ntuple).
Output   one JSON file: the cuts, the number of events and tracks, one entry per track
         (run, subrun, event, trk, p, d, t, nhits, hits, calo, crv).
Printed  "events N tracks N with calo N with crv N all N" (all = with calo and CRV).

Usage
  python3 kpp_dump_events.py out.json nts1.root [nts2.root ...]
  e.g. python3 Offline/MC2Reco/scripts/kpp_dump_events.py ev_125400_001.json nts_125400_001.root
"""
import argparse
import json

import ROOT

EMIN = 5.                    # MeV, calorimeter hit threshold
CALOFF = 908.                # ns, calorimeter time - track time at the crystal (calorimeter-tracker offset)
DTCAL = 150.                 # ns, half window around CALOFF
CRVLO, CRVHI = 0., 1000.     # ns, window of CRV time - track time at the CRV height
HOT = {1127}                 # hot crystal, excluded
N_CRYSTALS = 1348            # crystal ids below this; above: PIN diodes and other channels
FIRST_CRYSTAL_DISK1 = 674    # first crystal id of disk 1
SPEED_OF_LIGHT = 299.792458  # mm/ns


def parse_args():
    """Command line: output JSON and the EventNtuple files."""
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("output", help="output JSON file")
    parser.add_argument("ntuples", nargs="*", help="EventNtuple files of kpp_ntuple.fcl")
    return parser.parse_args()


class TrackLine:
    """Straight track through point p with unit direction d and time tref at p (the particle moves at c)."""

    def __init__(self, point, direction, tref):
        self.p = point
        self.d = direction
        self.tref = tref

    def at_y(self, y):
        """Point and time of the line at height y (tracker frame, y up)."""
        p, d = self.p, self.d
        k = (y - p[1]) / d[1]
        return [p[0] + k * d[0], y, p[2] + k * d[2]], self.tref + k / SPEED_OF_LIGHT

    def at_z(self, z):
        """Point and time of the line at the plane z; (None, None) if the line is parallel to it."""
        p, d = self.p, self.d
        if abs(d[2]) < 1e-6:
            return None, None
        k = (z - p[2]) / d[2]
        return [p[0] + k * d[0], p[1] + k * d[1], z], self.tref + k / SPEED_OF_LIGHT


def track_line(segment):
    """Line of the track from its middle segment, direction pointing down."""
    point = [segment.pos.x(), segment.pos.y(), segment.pos.z()]
    momentum = segment.mom.R()
    direction = [segment.mom.x() / momentum, segment.mom.y() / momentum, segment.mom.z() / momentum]
    if direction[1] > 0:          # cosmics go down: direction with y < 0
        direction = [-x for x in direction]
    return TrackLine(point, direction, segment.time)


def straw_hits(track_hits):
    """POCA and track time there of the active straw hits (state >= 0)."""
    return [[h.poca.x(), h.poca.y(), h.poca.z(), h.ptoca] for h in track_hits if getattr(h, "state", 0) >= 0]


def calo_hits_in_time(event, line):
    """Calorimeter crystals above EMIN, not PIN / hot, with dt = t_hit - t_track(z) within DTCAL of CALOFF."""
    calo = []
    for h in event.calohits:
        if h.eDep_ < EMIN or h.crystalId_ >= N_CRYSTALS or h.crystalId_ in HOT:
            continue
        position = [h.crystalPos_.x(), h.crystalPos_.y(), h.crystalPos_.z()]
        predicted, track_time = line.at_z(position[2])
        if predicted is None:
            continue
        dt = h.time_ - track_time
        if abs(dt - CALOFF) < DTCAL:
            calo.append(dict(id=h.crystalId_, disk=int(h.crystalId_ >= FIRST_CRYSTAL_DISK1), pos=position, E=h.eDep_,
                             t=h.time_, dt=dt, pred=predicted))
    return calo


def crv_in_time(event, line):
    """CRV coincidences with dt = t_CRV - t_track(y of the CRV) in (CRVLO, CRVHI)."""
    crv = []
    for c in event.crvcoincs:
        position = [c.pos.x(), c.pos.y(), c.pos.z()]
        predicted, track_time = line.at_y(position[1])
        dt = c.time - track_time
        if CRVLO < dt < CRVHI:
            crv.append(dict(type=int(c.sectorType), pos=position, PEs=c.PEs, t=c.time, dt=dt, pred=predicted))
    return crv


def dump_file(filename, out):
    """Add the events and tracks of one EventNtuple file to out."""
    f = ROOT.TFile.Open(filename)
    tree = f.Get("EventNtuple/ntuple")
    for event in tree:
        out["nevt"] += 1
        for itrk in range(len(event.trk)):
            segments = event.trksegs[itrk]
            if not len(segments):
                continue
            out["ntrk"] += 1
            line = track_line(segments[len(segments) // 2])
            hits = straw_hits(event.trkhits[itrk])
            calo = calo_hits_in_time(event, line)
            crv = crv_in_time(event, line)
            out["events"].append(dict(run=event.evtinfo.run, subrun=event.evtinfo.subrun, event=event.evtinfo.event,
                                      trk=itrk, p=line.p, d=line.d, t=line.tref, nhits=len(hits), hits=hits,
                                      calo=calo, crv=crv))


def main():
    """Dump all files, write the JSON, print the counts."""
    args = parse_args()
    out = dict(cuts=dict(EMIN=EMIN, CALOFF=CALOFF, DTCAL=DTCAL, CRVLO=CRVLO, CRVHI=CRVHI, hot=sorted(HOT)),
               events=[], ntrk=0, nevt=0)
    for filename in args.ntuples:
        dump_file(filename, out)
    json.dump(out, open(args.output, "w"))
    tracks = out["events"]
    print("events", out["nevt"], "tracks", out["ntrk"], "with calo", sum(1 for x in tracks if x["calo"]),
          "with crv", sum(1 for x in tracks if x["crv"]), "all", sum(1 for x in tracks if x["calo"] and x["crv"]))


if __name__ == "__main__":
    main()
