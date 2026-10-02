"""
Project 3 - track availability.

Some (event, lead) combinations cannot be run, and not because of a bug.
IBTrACS begins when a system is first classified, so a forecast origin that
falls before the record starts is asking what the forecast was for a storm
that did not yet exist. Titli is the case in this set: its record begins
8 Oct 2018 00:00 UTC and its T-72 origin is 7 Oct 22:34 UTC, 1h26m earlier.

This is the same thing section 4.3 says about the cone phase, at a different
lead. An event response team at T-72h for Titli had nothing to respond to.
Reporting the panel as unavailable is the finding; back-extrapolating 86
minutes of track for a system that had not formed would be inventing the one
thing the project is about not inventing.

Two behaviours are offered and the choice is recorded per event:

    skip  (default)  Report the combination as unavailable, with the reason
                     and the shortfall in hours. Every reported T-72h then
                     means the same thing across all four events.

    clamp            Move the origin forward to the first available track
                     time and run at the SHORTER actual lead. Invents nothing,
                     but Titli's "T-72h" would really be T-70.6h, so the
                     panels stop being comparable and the landfall calibration
                     assertion no longer has a ladder entry to check against.

Default is skip. Clamp exists so the cost of the choice can be measured
rather than argued about.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path

import numpy as np


# Clamping is only ever considered inside this window. Beyond it the forecast
# origin is not "slightly early", it is a different scenario.
MAX_CLAMP_H = 3.0


class InsufficientTrackHistory(Exception):
    """The forecast origin precedes the start of the best-track record."""

    def __init__(self, event: str, lead_time_h: int, shortfall_h: float,
                 first_available: np.datetime64, requested_origin: np.datetime64):
        self.event = event
        self.lead_time_h = lead_time_h
        self.shortfall_h = float(shortfall_h)
        self.first_available = first_available
        self.requested_origin = requested_origin
        super().__init__(
            f"{event} T-{lead_time_h}h: forecast origin "
            f"{str(requested_origin)[:16]} precedes the start of the IBTrACS "
            f"record ({str(first_available)[:16]}) by {shortfall_h:.2f} h. "
            f"The system was not yet classified, so there was no storm to "
            f"forecast. This combination is unavailable, not broken."
        )


@dataclass(frozen=True)
class AvailabilityRecord:
    event: str
    lead_time_h: int
    available: bool
    reason: str
    shortfall_h: float
    first_track_time: str
    requested_origin: str
    actual_lead_h: float


def check_origin(event: str, lead_time_h: int,
                 landfall_time: np.datetime64,
                 first_track_time: np.datetime64,
                 policy: str = "skip") -> tuple[np.datetime64, float, AvailabilityRecord]:
    """Validate (and optionally clamp) a forecast origin against the record.

    Returns (origin, actual_lead_h, record). Raises InsufficientTrackHistory
    under the default 'skip' policy when the origin is unreachable.
    """
    if policy not in ("skip", "clamp"):
        raise ValueError(f"Unknown policy '{policy}'. Use 'skip' or 'clamp'.")

    requested = landfall_time - np.timedelta64(int(lead_time_h), "h")
    shortfall_h = float(
        (first_track_time - requested) / np.timedelta64(1, "h")
    )

    if shortfall_h <= 0:
        return requested, float(lead_time_h), AvailabilityRecord(
            event=event, lead_time_h=lead_time_h, available=True,
            reason="ok", shortfall_h=0.0,
            first_track_time=str(first_track_time)[:16],
            requested_origin=str(requested)[:16],
            actual_lead_h=float(lead_time_h),
        )

    if policy == "clamp" and shortfall_h <= MAX_CLAMP_H:
        actual = float((landfall_time - first_track_time) / np.timedelta64(1, "h"))
        return first_track_time, actual, AvailabilityRecord(
            event=event, lead_time_h=lead_time_h, available=True,
            reason=f"clamped to first track time, actual lead {actual:.2f} h",
            shortfall_h=shortfall_h,
            first_track_time=str(first_track_time)[:16],
            requested_origin=str(requested)[:16],
            actual_lead_h=actual,
        )

    raise InsufficientTrackHistory(
        event, lead_time_h, shortfall_h, first_track_time, requested
    )


def record_unavailable(path: Path, rec: AvailabilityRecord) -> None:
    """Append one availability record, so every skipped panel is auditable.

    A panel that is simply absent from the outputs looks like a run that was
    forgotten. A panel that is absent AND recorded with its reason is a result.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = list(AvailabilityRecord.__dataclass_fields__)
    exists = path.exists()

    rows = []
    if exists:
        with path.open("r", newline="", encoding="utf-8") as f:
            rows = [r for r in csv.DictReader(f)
                    if not (r["event"] == rec.event
                            and int(r["lead_time_h"]) == rec.lead_time_h)]

    rows.append({k: getattr(rec, k) for k in fields})
    rows.sort(key=lambda r: (str(r["event"]), -int(r["lead_time_h"])))

    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)
