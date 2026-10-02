"""
Project 3 - Stage 4: impact and loss.

Takes each Stage 3 hazard set (500 perturbed wind fields on the exact
Project 1 525-cell grid), applies the Project 1 Odisha vulnerability function
and a Project 2 exposure scenario, and produces 500 event losses.

Generalised from the Fani-only benchmark
----------------------------------------
1. Event-aware. The benchmark hardcoded 'fani_T{lead}' in both the hazard
   filename and the output filename. That is the same defect that made
   Stage 5 and Stage 6 fail on Phailin.

2. Exposure-scenario aware. The Project 2 linkage file carries a 'scenario'
   column, so Stage 4 can run the same hazard against PRIMARY, P001, P003 and
   DIRTY. That is what lets Stage 7 compare meteorological variance against
   exposure-data variance, which is the link between Project 2 and Project 3.

3. Writes member_index explicitly. Stage 7 joins losses back to the latent
   draws that produced them. Relying on row order for that join would work
   silently until the day it did not, so the member index is parsed from the
   CLIMADA event name and written into the CSV, and its ordering is asserted.

Off-grid exposure
-----------------
PRIMARY must sit entirely on the Project 1 grid and is asserted to. Degraded
scenarios need not: the DIRTY portfolio contains mislocated records, which is
the Project 2 finding, not a bug here. Those cells cannot be modelled on a
hazard grid that does not cover them, so they are excluded, counted, and their
TIV reported. Excluding them silently would understate the dirty scenario's
loss without saying so.
"""

from __future__ import annotations

import argparse
import pickle
import re
import sys
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd

from climada.entity import Exposures
from climada.entity.impact_funcs import ImpactFunc, ImpactFuncSet
from climada.engine import ImpactCalc

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

import config as cfg


# Project 2 scenario linkage. Override with P3_P2_LINKAGE if the repo moves.
P2_LINKAGE = cfg.P2_LINKAGE

DEFAULT_SCENARIO = "PRIMARY"

# The PRIMARY portfolio is the project's modelling basis, so it is pinned to
# the exact Project 2 figure. A drift here means Project 2's output changed
# underneath Project 3 and every loss number is on a different footing.
PRIMARY_EXPECTED_TIV = 145_684_374_688.5043450280
PRIMARY_EXPECTED_CELLS = 17

IMPF_ID = 2
IMPF_NAME = "Odisha coastal buildings - OSDMA"

EVENT_NAME_RE = re.compile(r"_m(\d{4})$")

# Panels that Stage 3 recorded as unavailable, e.g. a forecast origin earlier
# than the IBTrACS record. Those are expected absences and Stage 4 skips them.
# A hazard file missing WITHOUT such a record is a Stage 3 run that did not
# happen, and that still fails loudly. Collapsing the two into one --skip-missing
# flag would let a forgotten run disappear as quietly as a genuine result.
AVAILABILITY_FILE = cfg.OUTPUT_DIR / "track_availability.csv"


def known_unavailable() -> dict[tuple[str, int], str]:
    """Load (event, lead) panels Stage 3 recorded as legitimately absent."""
    if not AVAILABILITY_FILE.exists():
        return {}
    df = pd.read_csv(AVAILABILITY_FILE)
    if "available" not in df.columns:
        return {}
    out: dict[tuple[str, int], str] = {}
    for _, r in df.iterrows():
        avail = r["available"]
        if isinstance(avail, str):
            avail = avail.strip().lower() in ("true", "1", "yes")
        if not bool(avail):
            out[(str(r["event"]), int(r["lead_time_h"]))] = str(
                r.get("reason", "recorded as unavailable")
            )
    return out


# ---------------------------------------------------------------------
# Vulnerability
# ---------------------------------------------------------------------

def build_odisha_impact_function() -> ImpactFunc:
    """Reproduce the validated Project 1 Odisha OSDMA vulnerability function.

    The regression checks below are the reason this is safe to copy between
    projects: the curve is pinned point by point, so a CLIMADA version change
    that altered interpolation would fail here rather than quietly shift every
    loss number in Project 3.
    """
    intensity = np.array(
        [0.0, 20.0, 24.0, 25.0, 30.0, 35.0,
         40.0, 45.0, 50.0, 55.0, 60.0, 70.0],
        dtype=float,
    )
    mdd = np.array(
        [0.0, 0.0, 0.0, 0.100, 0.155, 0.211,
         0.358, 0.569, 0.756, 0.922, 0.984, 1.000],
        dtype=float,
    )

    impf = ImpactFunc()
    impf.haz_type = "TC"
    impf.id = IMPF_ID
    impf.name = IMPF_NAME
    impf.intensity_unit = "m/s"
    impf.intensity = intensity
    impf.mdd = mdd
    impf.paa = np.ones_like(mdd)
    impf.check()

    for wind, expected in zip(intensity, mdd):
        actual = float(impf.calc_mdr(float(wind)))
        if not np.isclose(actual, float(expected), atol=1e-10):
            raise AssertionError(
                f"Impact-function regression failed at {wind} m/s: "
                f"{actual} != {expected}"
            )

    return impf


# ---------------------------------------------------------------------
# Exposure
# ---------------------------------------------------------------------

def load_exposure_on_p1_grid(scenario: str = DEFAULT_SCENARIO) -> tuple:
    """Load one Project 2 exposure scenario onto the exact P1 525-cell grid.

    Returns (Exposures, diagnostics dict).
    """
    scenario = scenario.upper()

    grid = gpd.read_file(cfg.P1_GRID)
    if len(grid) != 525:
        raise AssertionError(f"Expected 525 P1 cells, got {len(grid)}")
    if grid.crs is None or str(grid.crs).upper() != "EPSG:4326":
        raise AssertionError(f"P1 grid CRS is {grid.crs}, expected EPSG:4326")

    linkage = pd.read_csv(P2_LINKAGE)
    required = {"scenario", "grid_lon", "grid_lat", "total_tiv"}
    missing = required - set(linkage.columns)
    if missing:
        raise AssertionError(f"Missing linkage columns: {sorted(missing)}")

    available = sorted(linkage["scenario"].str.upper().unique())
    if scenario not in available:
        raise ValueError(
            f"Scenario '{scenario}' not in the Project 2 linkage file. "
            f"Available: {available}"
        )

    sel = linkage.loc[linkage["scenario"].str.upper() == scenario].copy()
    if sel.empty:
        raise AssertionError(f"Scenario {scenario} has no rows.")

    grid_lat = np.round(grid.geometry.y.to_numpy(), 8)
    grid_lon = np.round(grid.geometry.x.to_numpy(), 8)
    sel["grid_lat"] = sel["grid_lat"].round(8)
    sel["grid_lon"] = sel["grid_lon"].round(8)

    p1_keys = {(float(lo), float(la)) for lo, la in zip(grid_lon, grid_lat)}
    sel_keys = [(float(lo), float(la))
                for lo, la in zip(sel["grid_lon"], sel["grid_lat"])]

    on_grid = np.array([k in p1_keys for k in sel_keys])
    off_grid_tiv = float(sel.loc[~on_grid, "total_tiv"].sum())
    n_off_grid = int((~on_grid).sum())

    if scenario == DEFAULT_SCENARIO and n_off_grid:
        raise AssertionError(
            f"PRIMARY contains {n_off_grid} cells outside the P1 525-cell "
            f"grid, holding {off_grid_tiv:,.2f}. The accumulation-eligible "
            "portfolio must lie entirely on the hazard grid."
        )

    sel = sel.loc[on_grid].copy()

    exposure_df = pd.DataFrame({"latitude": grid_lat, "longitude": grid_lon})
    lookup = {
        (float(r.grid_lon), float(r.grid_lat)): float(r.total_tiv)
        for r in sel.itertuples()
    }
    exposure_df["value"] = [
        lookup.get((float(lo), float(la)), 0.0)
        for lo, la in zip(exposure_df["longitude"], exposure_df["latitude"])
    ]
    exposure_df["impf_"] = IMPF_ID

    modelled_tiv = float(exposure_df["value"].sum())
    occupied = int((exposure_df["value"] > 0).sum())

    if scenario == DEFAULT_SCENARIO:
        if not np.isclose(modelled_tiv, PRIMARY_EXPECTED_TIV, rtol=0, atol=1e-2):
            raise AssertionError(
                f"PRIMARY TIV mismatch: {modelled_tiv} != {PRIMARY_EXPECTED_TIV}. "
                "Project 2's output has changed under Project 3."
            )
        if occupied != PRIMARY_EXPECTED_CELLS:
            raise AssertionError(
                f"PRIMARY occupied cells: {occupied} != {PRIMARY_EXPECTED_CELLS}"
            )
    else:
        if modelled_tiv <= 0:
            raise AssertionError(f"{scenario}: modelled TIV is not positive.")

    exposure = Exposures(exposure_df)
    exposure.check()

    return exposure, {
        "scenario": scenario,
        "modelled_tiv": modelled_tiv,
        "occupied_cells": occupied,
        "off_grid_cells": n_off_grid,
        "off_grid_tiv": off_grid_tiv,
        "available_scenarios": ",".join(available),
    }


# ---------------------------------------------------------------------
# Hazard
# ---------------------------------------------------------------------

def load_hazard(event: str, lead_h: int, n_members: int = cfg.N_MEMBERS):
    path = cfg.hazard_file(event, lead_h)
    if not path.exists():
        raise FileNotFoundError(
            f"Hazard file not found: {path}. Run Stage 3 for {event} T-{lead_h}h."
        )

    with open(path, "rb") as f:
        hazard = pickle.load(f)

    if hazard.size != n_members:
        raise AssertionError(
            f"Expected {n_members} hazard events, got {hazard.size}"
        )
    if hazard.intensity.shape != (n_members, 525):
        raise AssertionError(f"Unexpected hazard shape: {hazard.intensity.shape}")
    if len(hazard.centroids.lat) != 525:
        raise AssertionError("Hazard centroids are not the 525-cell P1 grid.")

    return hazard


def member_indices(hazard, n_members: int = cfg.N_MEMBERS) -> np.ndarray:
    """Recover the Stage 2 member index for every hazard event, in order.

    Stage 7 joins losses to the latent draws that produced them. If the hazard
    events were ever reordered, or a member dropped by the plausibility screen,
    a positional join would still line up and every downstream number would
    still look plausible. Parsing the index out of the CLIMADA event name and
    asserting it is 0..N-1 in order is what prevents that.
    """
    names = [str(n) for n in hazard.event_name]
    idx = []
    for n in names:
        m = EVENT_NAME_RE.search(n)
        if not m:
            raise AssertionError(
                f"Hazard event name '{n}' does not carry a member index. "
                "Stage 2's to_climada_tctracks must name members "
                "'{event}_T{lead}_m{index:04d}'."
            )
        idx.append(int(m.group(1)))

    idx = np.asarray(idx, dtype=int)

    if len(idx) != n_members:
        raise AssertionError(f"Expected {n_members} members, got {len(idx)}")
    if not np.array_equal(idx, np.arange(n_members)):
        raise AssertionError(
            "Hazard member indices are not 0..N-1 in order. Either the "
            "plausibility screen rejected members or the events were "
            "reordered. Stage 7's join to the latent draws would be wrong."
        )

    return idx


# ---------------------------------------------------------------------
# Runner
# ---------------------------------------------------------------------

def run_one(event: str, lead_h: int, scenario: str, impf_set,
            exposure, exp_diag: dict) -> dict:
    print(f"\n{'-' * 72}")
    print(f"{event} T-{lead_h}h, exposure scenario {scenario}")
    print(f"{'-' * 72}")

    hazard = load_hazard(event, lead_h)
    idx = member_indices(hazard)

    hz_lat = np.asarray(hazard.centroids.lat)
    hz_lon = np.asarray(hazard.centroids.lon)
    ex_lat = np.asarray(exposure.latitude)
    ex_lon = np.asarray(exposure.longitude)

    if not (np.allclose(hz_lat, ex_lat, atol=1e-8)
            and np.allclose(hz_lon, ex_lon, atol=1e-8)):
        raise AssertionError(
            "Hazard and exposure are not on the same grid, cell for cell."
        )
    print("  grid alignment            : PASSED")
    print(f"  member ordering (0..{len(idx)-1}) : PASSED")

    impact = ImpactCalc(exposure, impf_set, hazard).impact(save_mat=True)

    losses = np.asarray(impact.at_event, dtype=float)
    if len(losses) != cfg.N_MEMBERS:
        raise AssertionError(f"Expected {cfg.N_MEMBERS} losses, got {len(losses)}")
    if not np.all(np.isfinite(losses)):
        raise AssertionError("Non-finite losses.")
    if np.any(losses < 0):
        raise AssertionError("Negative losses.")

    total_tiv = float(exposure.value.sum())
    if np.max(losses) > total_tiv + 1e-6:
        raise AssertionError(
            f"Maximum loss {np.max(losses):,.2f} exceeds total insured value "
            f"{total_tiv:,.2f}."
        )

    is_primary = scenario.upper() == DEFAULT_SCENARIO
    losses_path = (cfg.loss_file(event, lead_h) if is_primary
                   else cfg.loss_file_variant(event, lead_h, scenario.lower()))
    impact_path = losses_path.with_name(
        losses_path.name.replace("_event_losses.csv", "_impact.pkl")
    )

    cfg.IMPACT_DIR.mkdir(parents=True, exist_ok=True)

    pd.DataFrame({
        "event_id": np.arange(len(losses)),
        "member_index": idx,
        "gross_loss": losses,
    }).to_csv(losses_path, index=False)

    with open(impact_path, "wb") as f:
        pickle.dump(impact, f)

    zero_frac = float(np.mean(losses <= 1e3))

    print(f"  mean gross loss           : {cfg.fmt_money(float(np.mean(losses)))}")
    print(f"  median gross loss         : {cfg.fmt_money(float(np.median(losses)))}")
    print(f"  maximum gross loss        : {cfg.fmt_money(float(np.max(losses)))}")
    print(f"  members with no loss      : {zero_frac:.1%}")
    print(f"  -> {losses_path.name}")

    return {
        "event": event,
        "lead_time_h": lead_h,
        "scenario": scenario.upper(),
        "members": len(losses),
        "modelled_tiv": exp_diag["modelled_tiv"],
        "occupied_cells": exp_diag["occupied_cells"],
        "off_grid_cells": exp_diag["off_grid_cells"],
        "off_grid_tiv": exp_diag["off_grid_tiv"],
        "mean_gross_loss": float(np.mean(losses)),
        "median_gross_loss": float(np.median(losses)),
        "max_gross_loss": float(np.max(losses)),
        "p_zero_loss": zero_frac,
        "losses_file": str(losses_path),
    }


def main() -> int:
    ap = argparse.ArgumentParser(description="Project 3 Stage 4: impact and loss.")
    ap.add_argument("--event", choices=cfg.EVENT_NAMES, action="append",
                    help="Event(s). Default: all four.")
    ap.add_argument("--lead", type=int, choices=cfg.LEADS_H, action="append",
                    help="Lead time(s). Default: the whole loss-phase ladder.")
    ap.add_argument("--scenario", default=DEFAULT_SCENARIO,
                    help="Project 2 exposure scenario. Default PRIMARY.")
    ap.add_argument("--all-scenarios", action="store_true",
                    help="Run every scenario in the Project 2 linkage file. "
                         "Enables the Stage 7 exposure comparison.")
    ap.add_argument("--skip-missing", action="store_true",
                    help="Also skip hazards that are missing with no recorded "
                         "reason. Use only when you know why they are absent.")
    args = ap.parse_args()

    events = tuple(args.event) if args.event else cfg.EVENT_NAMES
    leads = tuple(args.lead) if args.lead else cfg.LEADS_H

    print("=" * 72)
    print("STAGE 4 - IMPACT AND LOSS")
    print("=" * 72)

    print("\n[1] Building the Project 1 Odisha impact function...")
    impf = build_odisha_impact_function()
    impf_set = ImpactFuncSet()
    impf_set.append(impf)
    print(f"  id {IMPF_ID}, {IMPF_NAME}")
    print("  point-by-point regression against the P1 curve : PASSED")

    if args.all_scenarios:
        probe = pd.read_csv(P2_LINKAGE)
        scenarios = tuple(sorted(probe["scenario"].str.upper().unique()))
        print(f"\n  Running all Project 2 scenarios: {', '.join(scenarios)}")
    else:
        scenarios = (args.scenario.upper(),)

    unavailable = known_unavailable()
    if unavailable:
        print(f"\n  Known-unavailable panels (recorded by Stage 3):")
        for (ev, ld), why in sorted(unavailable.items()):
            print(f"    {ev} T-{ld}h : {why}")

    rows = []
    skipped = []
    for scenario in scenarios:
        print(f"\n[2] Loading exposure scenario {scenario} on the P1 grid...")
        exposure, diag = load_exposure_on_p1_grid(scenario)
        print(f"  modelled TIV   : {cfg.fmt_money(diag['modelled_tiv'])}")
        print(f"  occupied cells : {diag['occupied_cells']} of 525")
        if diag["off_grid_cells"]:
            print(f"  off-grid cells : {diag['off_grid_cells']} holding "
                  f"{cfg.fmt_money(diag['off_grid_tiv'])}, EXCLUDED")
            print("                   (mislocated exposure cannot be modelled "
                  "on a grid that does not cover it)")

        for event in events:
            for lead in leads:
                key = (event, lead)
                if key in unavailable:
                    print(f"\n  {event} T-{lead}h, scenario {scenario}: "
                          f"skipped, panel is unavailable")
                    print(f"    reason: {unavailable[key]}")
                    skipped.append({"event": event, "lead_time_h": lead,
                                    "scenario": scenario,
                                    "reason": unavailable[key]})
                    continue
                try:
                    rows.append(run_one(event, lead, scenario, impf_set,
                                        exposure, diag))
                except FileNotFoundError as exc:
                    if args.skip_missing:
                        print(f"\n  skipped {event} T-{lead}h "
                              f"(--skip-missing): {exc}")
                        skipped.append({"event": event, "lead_time_h": lead,
                                        "scenario": scenario,
                                        "reason": "hazard file absent, "
                                                  "no recorded reason"})
                        continue
                    print(f"\nFAILED: {exc}")
                    print("\nThis hazard is missing with no recorded reason, so "
                          "it is a Stage 3 run that did not happen rather than a "
                          "panel that cannot exist. Run Stage 3 for it, or pass "
                          "--skip-missing if you know why it is absent.")
                    return 1

    if not rows:
        # Nothing ran for one of two very different reasons. Either every
        # requested panel is recorded as unavailable, which is a result and
        # exits 0, or a hazard file is missing with no recorded reason, which
        # is a Stage 3 run that did not happen and exits 1. Conflating them is
        # how one correctly unavailable panel stops the whole pipeline.
        all_unavailable = bool(skipped) and all(
            (s["event"], s["lead_time_h"]) in unavailable for s in skipped
            )
        if all_unavailable:
            skip_path = cfg.IMPACT_DIR / "impact_skipped.csv"
            pd.DataFrame(skipped).to_csv(skip_path, index=False)
            print("\nNo panel ran because every requested panel is recorded as "
                    "unavailable. That is a result, not a failure. The reasons "
                    f"are in {cfg.OUTPUT_DIR / 'track_availability.csv'}.")
            print(f"Recorded in {skip_path}")
            return 0
        print("\nNothing ran, and no requested panel is recorded as "
                "unavailable, so a Stage 3 hazard file is missing with no "
                "reason. Run Stage 3 for it, or pass --skip-missing if you "
                "know why it is absent.")
        return 1

    summary = pd.DataFrame(rows)
    path = cfg.IMPACT_DIR / "impact_summary.csv"
    summary.to_csv(path, index=False)

    if skipped:
        skip_path = cfg.IMPACT_DIR / "impact_skipped.csv"
        pd.DataFrame(skipped).to_csv(skip_path, index=False)

    print("\n" + "=" * 72)
    print(f"STAGE 4 COMPLETE: {len(rows)} run(s), {len(skipped)} skipped")
    print("=" * 72)
    print(f"Summary: {path}")
    if skipped:
        print(f"Skipped: {skip_path}")
        print("\nSkipped panels are recorded rather than merely absent. A panel")
        print("that is simply missing from the outputs looks like a forgotten")
        print("run; a panel that is missing with its reason attached is a result.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
