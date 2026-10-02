"""
Project 3 - track coverage diagnostic.

Run this ONCE before launching the 20 Stage 3 runs. It reports, for every
event and every lead time, whether the IBTrACS record reaches far enough back
to support that forecast origin.

Finding out now that three panels are unrunnable is a minute's work. Finding
out one at a time, three hours into a sequential loop, is not.

    python check_track_coverage.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
from climada.hazard import TCTracks

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

import config as cfg
from track_availability import (
    MAX_CLAMP_H, AvailabilityRecord, check_origin,
    InsufficientTrackHistory, record_unavailable,
)

IBTRACS_PROVIDER = "usa"
IBTRACS_BASIN = "NI"


def first_track_time(sid: str) -> np.datetime64:
    tracks = TCTracks.from_ibtracs_netcdf(
        provider=IBTRACS_PROVIDER, basin=IBTRACS_BASIN, storm_id=sid
    )
    if len(tracks.data) != 1:
        raise ValueError(f"Expected one IBTrACS track for {sid}, "
                         f"found {len(tracks.data)}")
    t = TCTracks(data=[tracks.data[0]])
    t.equal_timestep(time_step_h=1)
    return np.asarray(t.data[0].time.values)[0]


def main() -> int:
    events = cfg.load_events()
    out_path = cfg.OUTPUT_DIR / "track_availability.csv"

    print("=" * 76)
    print("PROJECT 3 - TRACK COVERAGE")
    print("=" * 76)
    print("\nIBTrACS begins when a system is first classified. A forecast")
    print("origin earlier than that is asking what the forecast was for a")
    print("storm that did not yet exist.\n")

    print(f"{'event':<9} {'first track':<17} {'landfall':<17} {'history':>9}")
    firsts = {}
    for name in cfg.EVENT_NAMES:
        ev = events[name]
        ft = first_track_time(ev.sid)
        firsts[name] = ft
        span = float((ev.landfall_time - ft) / np.timedelta64(1, "h"))
        print(f"{name:<9} {str(ft)[:16]:<17} {str(ev.landfall_time)[:16]:<17} "
              f"{span:8.2f}h")

    print(f"\n{'event':<9} {'lead':>6} {'origin':<17} {'status':<40}")
    n_unavailable = 0
    for name in cfg.EVENT_NAMES:
        ev = events[name]
        for lead in cfg.LEADS_H:
            try:
                origin, actual, rec = check_origin(
                    name, lead, ev.landfall_time, firsts[name], policy="skip"
                )
                status = "ok"
            except InsufficientTrackHistory as exc:
                n_unavailable += 1
                origin = exc.requested_origin
                clampable = exc.shortfall_h <= MAX_CLAMP_H
                status = (f"UNAVAILABLE, short by {exc.shortfall_h:.2f}h"
                          f"{' (clampable)' if clampable else ''}")
                rec = AvailabilityRecord(
                    event=name, lead_time_h=lead, available=False,
                    reason="forecast origin precedes the IBTrACS record; "
                           "system not yet classified",
                    shortfall_h=exc.shortfall_h,
                    first_track_time=str(firsts[name])[:16],
                    requested_origin=str(origin)[:16],
                    actual_lead_h=float(
                        (ev.landfall_time - firsts[name]) / np.timedelta64(1, "h")
                    ),
                )
                record_unavailable(out_path, rec)

            print(f"{name:<9} {'T-' + str(lead):>6} {str(origin)[:16]:<17} {status:<40}")

    print("\n" + "=" * 76)
    if n_unavailable:
        print(f"{n_unavailable} of {len(cfg.EVENT_NAMES)*len(cfg.LEADS_H)} "
              f"panels are unavailable.")
        print(f"Recorded in {out_path}")
        print("\nThese are a result, not a gap. A panel that is simply missing")
        print("looks like a forgotten run; a panel that is missing AND recorded")
        print("with its reason is a finding about when the decision clock")
        print("actually starts for that storm.")
    else:
        print("All panels have sufficient track history.")
    print("=" * 76)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
