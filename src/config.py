"""
Project 3 - Cyclone Event Response
Central configuration: paths, events, ladders and frozen decision rules.

Everything that more than one stage needs to agree on lives here. If two
stages disagree about a lead-time ladder or a file-naming convention, the
pipeline fails silently rather than loudly, which is the failure mode this
project exists to demonstrate.

Nothing in this module reads CLIMADA, so every downstream stage except
Stage 3 and Stage 4 can be unit-tested without it.
"""

from __future__ import annotations

import csv
import math
from dataclasses import dataclass
from pathlib import Path

import numpy as np


# =====================================================================
# Paths
# =====================================================================

SRC_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SRC_DIR.parent

DATA_DIR = PROJECT_ROOT / "data"
OUTPUT_DIR = PROJECT_ROOT / "outputs"

EVENTS_FILE = DATA_DIR / "events.csv"
ERROR_CSV = DATA_DIR / "imd_forecast_errors.csv"

HAZARD_DIR = OUTPUT_DIR / "hazard_benchmark"
IMPACT_DIR = OUTPUT_DIR / "impact_benchmark"
STAGE5_DIR = OUTPUT_DIR / "loss_distribution"
STAGE6_DIR = OUTPUT_DIR / "decision_layer"
STAGE7_DIR = OUTPUT_DIR / "variance_decomposition"
STAGE8_DIR = OUTPUT_DIR / "trigger_backtest"
STAGE9_DIR = OUTPUT_DIR / "event_briefs"

# External: Project 1 hazard grid. Override with P3_GRID_PATH if the repo
# lives somewhere else on another machine.
import os

P1_GRID = Path(
    os.environ.get(
        "P3_GRID_PATH",
        "/home/navneet/Projects/odisha-cyclone-risk/outputs/odisha_exposure_grid.gpkg",
    )
)

P2_LINKAGE = Path(
    os.environ.get(
        "P3_P2_LINKAGE",
        "/home/navneet/Projects/odisha-exposure-quality/"
        "outputs/aal_scenario_cell_linkage.csv",
    )
)

ALL_OUTPUT_DIRS = (
    HAZARD_DIR, IMPACT_DIR, STAGE5_DIR,
    STAGE6_DIR, STAGE7_DIR, STAGE8_DIR, STAGE9_DIR,
)


def ensure_output_dirs() -> None:
    for d in ALL_OUTPUT_DIRS:
        d.mkdir(parents=True, exist_ok=True)


# =====================================================================
# Events and ladders
# =====================================================================

EVENT_NAMES = ("Fani", "Phailin", "Titli", "Amphan")

# Loss-phase ladder. 96 h and 120 h are approach-path only: they carry no
# verified landfall-point error (FORECAST_ERROR_MODEL.md section 4), so they
# are not run as loss scenarios.
LEADS_H = (72, 48, 24, 12, 0)

N_MEMBERS = 500
BASE_SEED = 20260930

# Per-event narrative context. Amphan is the near-miss: its forecast cone
# included Odisha at long lead and its actual landfall did not. It is the
# case that makes the Stage 8 trigger backtest worth running.
EVENT_NOTES = {
    "Fani": "Direct hit near Puri, the highest-TIV part of the portfolio.",
    "Phailin": "Landfall at Gopalpur, south of the portfolio centre.",
    "Titli": "Landfall south of Gopalpur with a sharp post-landfall recurve.",
    "Amphan": "Near miss. Landfall in West Bengal, outside the Odisha portfolio.",
}

# Whether each event actually damaged THIS portfolio, used by Stage 8 as
# ground truth for scoring the notification trigger.
#
# These are a fallback only. Judging a "hit" by whether the storm made
# landfall on the Odisha coast is the wrong test: the portfolio has 39% of
# its TIV in a single cell near Puri, so a cyclone that lands 90 km away is
# a miss in loss terms however severe it was meteorologically. Scoring the
# trigger against landfall geography would mark the trigger as having missed
# events that caused this portfolio no loss at all.
#
# event_hit_portfolio() below derives the answer from the T-0 loss instead,
# which is the realised outcome with the forecast uncertainty gone.
EVENT_HIT_PORTFOLIO_FALLBACK = {
    "Fani": True,
    "Phailin": False,
    "Titli": False,
    "Amphan": False,
}

# A storm "hit" the portfolio if the realised loss would have reached the
# CAT XL attachment point. That is the same bar the trigger is being asked
# to predict, so the trigger is scored against the question it answers.
HIT_CRITERION_LEAD_H = 0


def event_hit_portfolio(event: str) -> bool:
    """Ground truth for Stage 8, derived from the T-0 loss where available.

    At T-0 the meteorological uncertainty is analysis uncertainty only, so the
    median T-0 loss is the realised outcome. Falls back to the table above if
    T-0 has not been run.
    """
    import numpy as np

    path = loss_file(event, HIT_CRITERION_LEAD_H)
    if not path.exists():
        return EVENT_HIT_PORTFOLIO_FALLBACK[event]
    losses = load_losses(event, HIT_CRITERION_LEAD_H)
    return bool(float(np.median(losses)) > ATTACHMENT)


class _HitPortfolio(dict):
    """Mapping interface so existing cfg.EVENT_HIT_PORTFOLIO[event] still works."""

    def __missing__(self, key):
        return event_hit_portfolio(key)

    def __getitem__(self, key):
        return event_hit_portfolio(key)

    def items(self):
        return ((e, event_hit_portfolio(e)) for e in EVENT_NAMES)


EVENT_HIT_PORTFOLIO = _HitPortfolio()


@dataclass(frozen=True)
class EventRecord:
    event: str
    sid: str
    landfall_time: np.datetime64
    landfall_lat: float
    landfall_lon: float
    landfall_method: str
    coastline_dataset: str
    imd_reference: np.datetime64
    imd_difference_min: float


def load_events(path: Path = EVENTS_FILE) -> dict[str, EventRecord]:
    """Load the canonical landfall table and validate it against IMD.

    The landfall definition is the geometric sea-to-land crossing of the
    interpolated best track (Option A). IMD's published reference is the
    independent check, not the input, so the agreement tolerance is asserted
    here rather than assumed.
    """
    if not path.exists():
        raise FileNotFoundError(f"Event table not found: {path}")

    events: dict[str, EventRecord] = {}
    with path.open("r", newline="", encoding="utf-8") as f:
        for row in csv.DictReader(f):
            rec = EventRecord(
                event=row["event"],
                sid=row["sid"],
                landfall_time=np.datetime64(row["landfall_time"]),
                landfall_lat=float(row["landfall_lat"]),
                landfall_lon=float(row["landfall_lon"]),
                landfall_method=row["landfall_method"],
                coastline_dataset=row["coastline_dataset"],
                imd_reference=np.datetime64(row["imd_reference"]),
                imd_difference_min=float(row["imd_difference_min"]),
            )
            events[rec.event] = rec

    missing = set(EVENT_NAMES) - set(events)
    if missing:
        raise AssertionError(f"Event table is missing: {sorted(missing)}")

    for rec in events.values():
        if rec.landfall_method != "geometric_sea_to_land":
            raise AssertionError(
                f"{rec.event}: landfall_method is '{rec.landfall_method}'. "
                "All four events must use the same rule (Option A); an "
                "event-dependent definition is exactly what was rejected."
            )
        # IMD cross-check. 1 h is roughly 20 km at typical translation speed,
        # inside both coastline representation error and one 27.75 km grid cell.
        if abs(rec.imd_difference_min) > MAX_IMD_DIFFERENCE_MIN:
            raise AssertionError(
                f"{rec.event}: derived landfall differs from the IMD reference "
                f"by {rec.imd_difference_min:.0f} min, outside the stated "
                f"+/-{MAX_IMD_DIFFERENCE_MIN:.0f} min tolerance."
            )

    return events


MAX_IMD_DIFFERENCE_MIN = 60.0


# =====================================================================
# Frozen decision rules (Stage 6)
# =====================================================================

RESERVE_PERCENTILE = 75

# ---------------------------------------------------------------------
# CAT XL layer, and why it is denominated the way it is.
#
# Project 1 prices its layer in USD on a USD exposure base:
#     $10B xs $12B on $144.86B of LitPop exposure.
#
# Project 3 computes losses in INR on Project 2's PRIMARY portfolio:
#     Rs 145.684B TIV.
#
# The layer is therefore NOT copied across as a number. It is rescaled to
# the same share of total exposure that Project 1's layer represented:
#
#     attachment share = 12.00 / 144.86 = 8.284%
#     limit share      = 10.00 / 144.86 = 6.903%
#
# Applied to Rs 145.684B this gives Rs 12.07B xs and Rs 10.06B of limit,
# which round to the figures below. The near-coincidence of the two
# portfolio magnitudes (144.86 and 145.684 in different currencies) means a
# naive currency copy would also have produced roughly these numbers. That
# is precisely why the derivation is written out instead of asserted.
# ---------------------------------------------------------------------

P1_EXPOSURE_USD_B = 144.86
P1_ATTACHMENT_USD_B = 12.00
P1_LIMIT_USD_B = 10.00

P2_PRIMARY_TIV_INR_B = 145.684

ATTACHMENT_SHARE = P1_ATTACHMENT_USD_B / P1_EXPOSURE_USD_B
LIMIT_SHARE = P1_LIMIT_USD_B / P1_EXPOSURE_USD_B

# Rounded to the nearest Rs 1B for a quotable layer, with the unrounded
# derivation retained so the rounding is visible.
ATTACHMENT_EXACT = P2_PRIMARY_TIV_INR_B * ATTACHMENT_SHARE * 1e9
LIMIT_EXACT = P2_PRIMARY_TIV_INR_B * LIMIT_SHARE * 1e9

ATTACHMENT = 12e9
LIMIT = 10e9

NOTIFICATION_THRESHOLD = 0.50

CURRENCY_SYMBOL = "Rs "
CURRENCY_CODE = "INR"


def describe_layer() -> str:
    return (
        f"CAT XL {CURRENCY_SYMBOL}{LIMIT/1e9:.0f}B xs "
        f"{CURRENCY_SYMBOL}{ATTACHMENT/1e9:.0f}B "
        f"(rescaled from Project 1's ${P1_LIMIT_USD_B:.0f}B xs "
        f"${P1_ATTACHMENT_USD_B:.0f}B at the same share of exposure: "
        f"{ATTACHMENT_SHARE:.3%} attachment, {LIMIT_SHARE:.3%} limit; "
        f"unrounded {CURRENCY_SYMBOL}{LIMIT_EXACT/1e9:.2f}B xs "
        f"{CURRENCY_SYMBOL}{ATTACHMENT_EXACT/1e9:.2f}B)"
    )


# =====================================================================
# File naming - one convention, used by every stage
# =====================================================================

def event_key(event: str) -> str:
    return event.lower()


def hazard_file(event: str, lead_h: int) -> Path:
    return HAZARD_DIR / f"{event_key(event)}_T{lead_h}_hazard.pkl"


def loss_file(event: str, lead_h: int) -> Path:
    return IMPACT_DIR / f"{event_key(event)}_T{lead_h}_event_losses.csv"


def loss_file_variant(event: str, lead_h: int, variant: str) -> Path:
    """Loss file for a perturbation variant (Stage 7 exposure scenarios)."""
    return IMPACT_DIR / f"{event_key(event)}_T{lead_h}_{variant}_event_losses.csv"


# =====================================================================
# Loss loading, with the validation every stage needs
# =====================================================================

def load_losses(
    event: str,
    lead_h: int,
    n_members: int = N_MEMBERS,
    path: Path | None = None,
) -> np.ndarray:
    """Load the Stage 4 simulated gross losses for one event and lead time.

    Member order is the Stage 2 ensemble order, which Stage 7 relies on to
    join losses back to their latent draws. That ordering is asserted here.
    """
    import pandas as pd

    p = path if path is not None else loss_file(event, lead_h)
    if not p.exists():
        raise FileNotFoundError(
            f"{event} T-{lead_h}: missing loss file {p}. "
            "Run Stage 4 for this event and lead before this stage."
        )

    df = pd.read_csv(p)

    required = {"event_id", "gross_loss"}
    missing = required - set(df.columns)
    if missing:
        raise AssertionError(f"{p.name}: missing columns {sorted(missing)}")

    if len(df) != n_members:
        raise AssertionError(
            f"{p.name}: expected {n_members} losses, found {len(df)}"
        )
    if not df["event_id"].is_unique:
        raise AssertionError(f"{p.name}: event_id values are not unique")

    # Stage 4 writes member_index: the Stage 2 ensemble member that produced
    # each loss. Stage 7 joins losses back to the latent draws, so the array
    # returned here must be ordered by member index, not by row order. Sorting
    # on it explicitly is what stops a reordered file from producing a
    # perfectly plausible but completely wrong decomposition.
    if "member_index" in df.columns:
        mi = df["member_index"].to_numpy(dtype=int)
        if not np.array_equal(np.sort(mi), np.arange(n_members)):
            raise AssertionError(
                f"{p.name}: member_index is not a permutation of 0..{n_members-1}. "
                "Stage 7 cannot join losses to their latent draws."
            )
        df = df.sort_values("member_index")

    losses = df["gross_loss"].to_numpy(dtype=float)

    if not np.all(np.isfinite(losses)):
        raise AssertionError(f"{p.name}: non-finite losses")
    if np.any(losses < 0):
        raise AssertionError(f"{p.name}: negative losses")

    return losses


def available_leads(event: str) -> tuple[int, ...]:
    return tuple(l for l in LEADS_H if loss_file(event, l).exists())


def available_events() -> tuple[str, ...]:
    return tuple(e for e in EVENT_NAMES if available_leads(e))


# =====================================================================
# Formatting
# =====================================================================

def fmt_money(value: float, decimals: int = 2) -> str:
    if not np.isfinite(value):
        return "n/a"
    return f"{CURRENCY_SYMBOL}{value/1e9:,.{decimals}f}B"


def fmt_lead(lead_h: int) -> str:
    return f"T-{lead_h}h"


# =====================================================================
# Self-check
# =====================================================================

if __name__ == "__main__":
    print("=" * 72)
    print("PROJECT 3 CONFIGURATION")
    print("=" * 72)
    ensure_output_dirs()

    events = load_events()
    print(f"\nEvents             : {', '.join(EVENT_NAMES)}")
    print(f"Loss-phase ladder  : {LEADS_H}")
    print(f"Members            : {N_MEMBERS}")
    print(f"Base seed          : {BASE_SEED}")
    print(f"\n{describe_layer()}")
    print(f"Reserve            : P{RESERVE_PERCENTILE}")
    print(f"Notification       : P(attach) >= {NOTIFICATION_THRESHOLD:.0%}")

    print("\nCanonical landfall table (geometric sea-to-land crossing):")
    print(f"{'event':<9} {'landfall (UTC)':<20} {'lat':>9} {'lon':>9} "
          f"{'IMD delta':>10}")
    for name in EVENT_NAMES:
        e = events[name]
        # Seconds are spurious precision: the best track is interpolated and
        # the coastline carries km-scale representation error.
        shown = str(e.landfall_time)[:16]
        print(f"{e.event:<9} {shown:<20} {e.landfall_lat:9.4f} "
              f"{e.landfall_lon:9.4f} {e.imd_difference_min:8.0f}m")

    print(f"\nAll four within +/-{MAX_IMD_DIFFERENCE_MIN:.0f} min of IMD: PASSED")
    print("\nConfiguration self-check PASSED")
