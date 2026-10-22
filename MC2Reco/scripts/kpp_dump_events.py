"""Events of the KPP check of the geometry (3D display, residuals, times), from the EventNtuple of kpp_ntuple.fcl (detector frame, mm, ns).
Per tracker track (KKLine): line from the segments, straw hit POCAs; calo hits (MC2 chain, E > EMIN, not PIN diodes, not the hot
crystal) within +-DTCAL of t0 + CALOFF; CRV coincidences within [CRVLO, CRVHI] of the track time at the CRV.
usage: python3 kpp_dump_events.py out.json nts1.root [nts2.root ...]"""
import json, sys, ROOT
EMIN, CALOFF, DTCAL, CRVLO, CRVHI, HOT = 5., 908., 150., 0., 1000., {1127}
out = dict(cuts=dict(EMIN=EMIN, CALOFF=CALOFF, DTCAL=DTCAL, CRVLO=CRVLO, CRVHI=CRVHI, hot=sorted(HOT)), events=[], ntrk=0, nevt=0)
for fn in sys.argv[2:]:
    f = ROOT.TFile.Open(fn); t = f.Get("EventNtuple/ntuple")
    for e in t:
        out["nevt"] += 1
        for it in range(len(e.trk)):
            segs = e.trksegs[it]
            if not len(segs): continue
            out["ntrk"] += 1
            s = segs[len(segs) // 2]; p = [s.pos.x(), s.pos.y(), s.pos.z()]; R = s.mom.R(); d = [s.mom.x() / R, s.mom.y() / R, s.mom.z() / R]
            if d[1] > 0: d = [-x for x in d]          # cosmics go down: direction with y < 0
            tref = s.time
            def at_y(y):   # point and time of the line at height y (tracker frame, y up)
                k = (y - p[1]) / d[1]; return [p[0] + k * d[0], y, p[2] + k * d[2]], tref + k / 299.792458
            def at_z(z):
                if abs(d[2]) < 1e-6: return None, None
                k = (z - p[2]) / d[2]; return [p[0] + k * d[0], p[1] + k * d[1], z], tref + k / 299.792458
            hits = [[h.poca.x(), h.poca.y(), h.poca.z(), h.ptoca] for h in e.trkhits[it] if getattr(h, "state", 0) >= 0]   # POCA and track time there
            calo = []
            for h in e.calohits:
                if h.eDep_ < EMIN or h.crystalId_ >= 1348 or h.crystalId_ in HOT: continue
                q = [h.crystalPos_.x(), h.crystalPos_.y(), h.crystalPos_.z()]; pr, tp = at_z(q[2])
                if pr is None: continue
                dt = h.time_ - tp
                if abs(dt - CALOFF) < DTCAL: calo.append(dict(id=h.crystalId_, disk=int(h.crystalId_ >= 674), pos=q, E=h.eDep_, t=h.time_, dt=dt, pred=pr))
            crv = []
            for c in e.crvcoincs:
                q = [c.pos.x(), c.pos.y(), c.pos.z()]; pr, tp = at_y(q[1]); dt = c.time - tp
                if CRVLO < dt < CRVHI: crv.append(dict(type=int(c.sectorType), pos=q, PEs=c.PEs, t=c.time, dt=dt, pred=pr))
            out["events"].append(dict(run=e.evtinfo.run, subrun=e.evtinfo.subrun, event=e.evtinfo.event, trk=it, p=p, d=d, t=tref,
                                      nhits=len(hits), hits=hits, calo=calo, crv=crv))
json.dump(out, open(sys.argv[1], "w")); ev = out["events"]
print("events", out["nevt"], "tracks", out["ntrk"], "with calo", sum(1 for x in ev if x["calo"]), "with crv", sum(1 for x in ev if x["crv"]), "all", sum(1 for x in ev if x["calo"] and x["crv"]))
