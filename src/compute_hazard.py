
"""
Project 3 — Stage 3 hazard benchmark.

Generic four-event Stage 3 runner.

Event:
    Fani / Phailin / Titli / Amphan

Pipeline:
    canonical geometric landfall
        -> historical IBTrACS physical track
        -> 500-member forecast ensemble
        -> CLIMADA TCTracks
        -> exact P1 525-cell hazard grid
        -> TropCyclone wind fields

This script intentionally benchmarks ONE event/lead at a time.
No exposure or impact calculation is performed here.
"""

from pathlib import Path
import sys
import time
import pickle
import csv
import argparse

import config as cfg
import geopandas as gpd
import numpy as np
import cartopy.io.shapereader as shpreader
from pyproj import Geod
from shapely.ops import nearest_points
from climada.hazard import Centroids, TropCyclone, TCTracks


# ---------------------------------------------------------------------
# Local imports
# ---------------------------------------------------------------------

HERE = Path(__file__).resolve().parent
REPO = HERE.parent

if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from error_model import load_error_models, landfall_sampler_parameters
from generate_ensemble import (
    generate_ensemble,
    to_climada_tctracks,
    EXTENSION_SIGMAS,
)

# ---------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------

N_MEMBERS = 500
BASE_SEED = 20260930

VALID_LEADS_H = (72, 48, 24, 12, 0)

# One source of truth, with the P3_GRID_PATH override, so the repo runs on a
# machine where Project 1 lives somewhere else.
P1_GRID = cfg.P1_GRID

OUTPUT_DIR = REPO / "outputs" / "hazard_benchmark"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

EVENTS_FILE = REPO / "data" / "events.csv"

IBTRACS_PROVIDER = "usa"
IBTRACS_BASIN = "NI"


# ---------------------------------------------------------------------
# Arguments
# ---------------------------------------------------------------------

def parse_args():
    parser = argparse.ArgumentParser(
        description="Project 3 Stage 3 hazard benchmark."
    )

    parser.add_argument(
        "--event",
        required=True,
        choices=("Fani", "Phailin", "Titli", "Amphan"),
        help="Cyclone event.",
    )

    parser.add_argument(
        "--lead",
        type=int,
        default=48,
        choices=VALID_LEADS_H,
        help="Forecast lead time in hours.",
    )

    return parser.parse_args()


ARGS = parse_args()
EVENT_NAME = ARGS.event
LEAD_TIME_H = ARGS.lead


# ---------------------------------------------------------------------
# Event metadata
# ---------------------------------------------------------------------

def load_event_metadata(event_name):
    rows = []

    with EVENTS_FILE.open(
        "r",
        newline="",
        encoding="utf-8",
    ) as f:
        reader = csv.DictReader(f)

        for row in reader:
            if row["event"] == event_name:
                rows.append(row)

    if len(rows) != 1:
        raise ValueError(
            f"Expected exactly one row for {event_name}, "
            f"found {len(rows)}."
        )

    row = rows[0]

    return {
        "event": row["event"],
        "sid": row["sid"],
        "landfall_time": np.datetime64(row["landfall_time"]),
        "landfall_lat": float(row["landfall_lat"]),
        "landfall_lon": float(row["landfall_lon"]),
        "landfall_method": row["landfall_method"],
        "coastline_dataset": row["coastline_dataset"],
    }


EVENT = load_event_metadata(EVENT_NAME)

LANDFALL_TIME = EVENT["landfall_time"]

OUTPUT_FILE = (
    OUTPUT_DIR
    / f"{EVENT_NAME.lower()}_T{LEAD_TIME_H}_hazard.pkl"
)


# ---------------------------------------------------------------------
# Historical IBTrACS track
# ---------------------------------------------------------------------

def load_historical_track(sid):
    print(
        f"\nLoading IBTrACS track: {sid}"
    )

    tracks = TCTracks.from_ibtracs_netcdf(
        provider=IBTRACS_PROVIDER,
        basin=IBTRACS_BASIN,
        storm_id=sid,
    )

    if len(tracks.data) != 1:
        raise ValueError(
            f"Expected exactly one IBTrACS track for {sid}, "
            f"found {len(tracks.data)}."
        )

    track = tracks.data[0]

    tracks_1h = TCTracks(data=[track])
    tracks_1h.equal_timestep(time_step_h=1)

    track = tracks_1h.data[0]

    return track


# ---------------------------------------------------------------------
# Interpolate historical physical track to canonical landfall clock
# ---------------------------------------------------------------------

def build_best_track(track, landfall_time, lead_time_h):
    """
    Construct the Stage 2 BestTrack on an hourly forecast horizon.

    The canonical geometric landfall timestamp is preserved exactly.
    Historical IBTrACS fields are linearly interpolated onto the
    corresponding physical timestamps.

    For positive leads the grid continues PAST landfall by as many hours as are
    needed to cover EXTENSION_SIGMAS times sigma_a at landfall, measured as
    path length along the real track. Without that extension the array ends at
    landfall, so a member displaced backward along track has every one of its
    positions offshore, the wind field never comes ashore, and the loss is
    exactly zero for a reason that has nothing to do with where the storm might
    have gone. See section 6.5.9.

    For T-0 one physical landfall point is returned and there is no extension:
    a single point admits no along-track displacement, and the T-0 panels
    contained no zero-loss members.

    Returns (best_track, target_times, extension_info).
    """

    from generate_ensemble import BestTrack
    from track_availability import check_origin

    geod_local = Geod(ellps="WGS84")

    time_raw = np.asarray(track.time.values)

    lat_raw = np.asarray(track.lat.values, dtype=float)
    lon_raw = np.asarray(track.lon.values, dtype=float)
    wind_raw = np.asarray(
        track.max_sustained_wind.values,
        dtype=float,
    )

    forecast_origin = (
        landfall_time
        - np.timedelta64(lead_time_h, "h")
    )

    source_seconds = (
        (time_raw - time_raw[0])
        / np.timedelta64(1, "s")
    ).astype(float)

    # If the requested forecast origin predates the first IBTrACS observation,
    # this panel is genuinely unavailable. Do not extrapolate or clamp.
    check_origin(
        EVENT_NAME,
        lead_time_h,
        landfall_time,
        time_raw[0],
        policy="skip",
    )

    landfall_seconds = float(
        (landfall_time - time_raw[0]) / np.timedelta64(1, "s")
    )

    if landfall_seconds > source_seconds.max():
        raise ValueError(
            f"{EVENT_NAME}: canonical landfall exceeds "
            "available IBTrACS data."
        )

    def interp_at(times):
        ts = (
            (times - time_raw[0])
            / np.timedelta64(1, "s")
        ).astype(float)
        return (
            np.interp(ts, source_seconds, lat_raw),
            np.interp(ts, source_seconds, lon_raw),
            np.interp(ts, source_seconds, wind_raw),
        )

    # -----------------------------------------------------------------
    # Size the post-landfall extension.
    #
    # The requirement is a DISTANCE, so the search walks path length and
    # converts to hours, rather than assuming a translation speed. A slow
    # storm needs more hours to cover the same ground.
    # -----------------------------------------------------------------

    sigma_a_landfall_km, _ = landfall_sampler_parameters(lead_time_h)
    required_km = float(EXTENSION_SIGMAS) * float(sigma_a_landfall_km)

    hours_available = int(
        np.floor(
            (time_raw[-1] - landfall_time) / np.timedelta64(1, "h")
        )
    )
    hours_available = max(hours_available, 0)

    extension_h = 0
    available_km = 0.0

    if lead_time_h > 0 and hours_available > 0:
        probe_times = (
            landfall_time
            + np.arange(hours_available + 1)
            * np.timedelta64(1, "h")
        )
        p_lat, p_lon, p_wind = interp_at(probe_times)

        _, _, seg_m = geod_local.inv(
            p_lon[:-1], p_lat[:-1], p_lon[1:], p_lat[1:]
        )
        cum_km = np.cumsum(np.asarray(seg_m, dtype=float) / 1000.0)

        # Stop before intensity reaches zero. BestTrack rejects a
        # non-positive intensity outright, and CLIMADA cannot use one.
        bad = np.flatnonzero(p_wind[1:] <= 0)
        max_h = int(bad[0]) if len(bad) else hours_available

        reached = np.flatnonzero(cum_km >= required_km)
        needed_h = int(reached[0]) + 1 if len(reached) else hours_available

        extension_h = int(min(needed_h, max_h))
        available_km = float(cum_km[extension_h - 1]) if extension_h else 0.0

    extension_info = {
        "event": EVENT_NAME,
        "lead_time_h": int(lead_time_h),
        "sigma_a_landfall_km": float(sigma_a_landfall_km),
        "required_extension_km": required_km,
        "available_extension_km": available_km,
        "extension_hours": int(extension_h),
        "hours_available_in_record": int(hours_available),
        "fully_covered": bool(
            lead_time_h == 0 or available_km >= required_km
        ),
    }

    # -----------------------------------------------------------------
    # Build the grid: origin -> landfall -> extension.
    # -----------------------------------------------------------------

    n_points = lead_time_h + 1 + extension_h

    target_times = (
        forecast_origin
        + np.arange(n_points)
        * np.timedelta64(1, "h")
    )

    lat, lon, wind = interp_at(target_times)

    # The grid is hourly from the origin, so landfall is at this index.
    i_landfall = int(lead_time_h)

    assert target_times[0] == forecast_origin
    assert target_times[i_landfall] == landfall_time

    horizon_h = (
        (target_times - target_times[0])
        / np.timedelta64(1, "h")
    ).astype(float)

    assert horizon_h[0] == 0
    assert horizon_h[i_landfall] == lead_time_h
    assert len(horizon_h) == n_points

    assert np.isfinite(lat).all()
    assert np.isfinite(lon).all()
    assert np.isfinite(wind).all()
    assert np.all(wind > 0)

    # -----------------------------------------------------------------
    # Endpoint should reproduce the canonical geometric landfall
    # position to within the interpolation/source tolerance.
    #
    # Indexed at LANDFALL, not at the end of the array, which is now
    # well inland.
    #
    # This is a diagnostic only. The canonical event table remains
    # authoritative for the landfall definition.
    # -----------------------------------------------------------------

    _, _, endpoint_error_m = geod_local.inv(
        lon[i_landfall],
        lat[i_landfall],
        EVENT["landfall_lon"],
        EVENT["landfall_lat"],
    )

    endpoint_error_km = endpoint_error_m / 1000.0

    print(f"\n{EVENT_NAME} canonical landfall:")
    print(f"  Time       : {str(landfall_time)}")
    print(
        f"  Event table: "
        f"{EVENT['landfall_lat']:.6f} N, "
        f"{EVENT['landfall_lon']:.6f} E"
    )
    print(
        f"  IBTrACS int: "
        f"{lat[i_landfall]:.6f} N, "
        f"{lon[i_landfall]:.6f} E"
    )
    print(f"  Endpoint difference: {endpoint_error_km:.3f} km")

    print(f"\n{EVENT_NAME} post-landfall extension (section 6.5.9):")
    print(f"  sigma_a at landfall : {sigma_a_landfall_km:.2f} km")
    print(
        f"  Required ({EXTENSION_SIGMAS:.0f} sigma) : "
        f"{required_km:.1f} km"
    )
    print(
        f"  Provided            : "
        f"{available_km:.1f} km over {extension_h} h"
    )
    if not extension_info["fully_covered"]:
        print(
            "  SHORT of the 4-sigma target. The best track record ends "
            "before it. validate_ensemble will check whether the shortfall "
            "actually bit the members drawn."
        )

    best_track = BestTrack(
        event_name=EVENT_NAME,
        horizon_h=horizon_h,
        lat=lat,
        lon=lon,
        max_wind_kt=wind,
        landfall_index=i_landfall,
    )

    return best_track, target_times, extension_info


# ---------------------------------------------------------------------
# Physical reference track for CLIMADA
# ---------------------------------------------------------------------

def build_reference_physical_track(
    track,
    landfall_time,
    lead_time_h,
    target_times=None,
):
    """Physical fields CLIMADA needs, on the ensemble's own time axis.

    target_times comes from build_best_track so that the two tracks share one
    grid by construction. to_climada_tctracks compares the timestamp arrays for
    exact equality, and deriving the grid independently in two places is how
    they come to disagree.
    """

    import xarray as xr

    time_raw = np.asarray(
        track.time.values
    )

    if lead_time_h == 0:

        previous_time = (
            landfall_time
            - np.timedelta64(
                1,
                "h",
            )
        )

        target_times = np.array(
            [
                previous_time,
                landfall_time,
            ],
            dtype="datetime64[ns]",
        )

    else:

        if target_times is None:
            raise ValueError(
                "target_times is required for a positive lead time. Pass the "
                "grid returned by build_best_track so the reference track and "
                "the ensemble share one time axis, including the "
                "post-landfall extension."
            )

        target_times = np.asarray(target_times)

    source_seconds = (
        (time_raw - time_raw[0])
        / np.timedelta64(
            1,
            "s",
        )
    ).astype(float)

    target_seconds = (
        (target_times - time_raw[0])
        / np.timedelta64(
            1,
            "s",
        )
    ).astype(float)


    if target_seconds.max() > source_seconds.max():
        raise ValueError(
            f"{EVENT_NAME}: reference segment ends "
            "after available IBTrACS data."
        )

    interpolated = {}

    for var in [
        "lat",
        "lon",
        "max_sustained_wind",
        "central_pressure",
        "environmental_pressure",
        "radius_max_wind",
    ]:

        values = np.asarray(
            track[var].values,
            dtype=float,
        )

        if not np.isfinite(values).all():
            raise ValueError(
                f"{EVENT_NAME} reference {var} "
                "contains non-finite values."
            )

        interpolated[var] = np.interp(
            target_seconds,
            source_seconds,
            values,
        )

    reference = xr.Dataset(
        {
            "lat": (
                ("time",),
                interpolated["lat"],
            ),
            "lon": (
                ("time",),
                interpolated["lon"],
            ),
            "max_sustained_wind": (
                ("time",),
                interpolated[
                    "max_sustained_wind"
                ],
            ),
            "central_pressure": (
                ("time",),
                interpolated[
                    "central_pressure"
                ],
            ),
            "environmental_pressure": (
                ("time",),
                interpolated[
                    "environmental_pressure"
                ],
            ),
            "radius_max_wind": (
                ("time",),
                interpolated[
                    "radius_max_wind"
                ],
            ),
        },
        coords={
            "time": target_times,
        },
    )

    assert len(reference.time) == len(target_times)

    assert len(reference.lat) == len(
        reference.time
    )

    assert len(reference.lon) == len(
        reference.time
    )

    assert len(
        reference.max_sustained_wind
    ) == len(reference.time)

    assert len(
        reference.central_pressure
    ) == len(reference.time)

    assert len(
        reference.environmental_pressure
    ) == len(reference.time)

    assert len(
        reference.radius_max_wind
    ) == len(reference.time)

    return reference


# ---------------------------------------------------------------------
# Stage 3 — exact P1 grid
# ---------------------------------------------------------------------

print("=" * 72)
print("PROJECT 3 — STAGE 3 HAZARD BENCHMARK")
print(f"Event: {EVENT_NAME}")
print(f"Forecast lead: T-{LEAD_TIME_H} h")
print("=" * 72)

print("\n[1] Loading exact P1 525-cell grid...")

grid = gpd.read_file(P1_GRID)

assert len(grid) == 525
assert str(grid.crs) == "EPSG:4326"
assert grid.geometry.notna().all()

grid_lon = grid.geometry.x.to_numpy()
grid_lat = grid.geometry.y.to_numpy()

assert np.isfinite(grid_lon).all()
assert np.isfinite(grid_lat).all()

print(f"Grid cells : {len(grid)}")
print(f"CRS        : {grid.crs}")
print(
    f"Bounds     : "
    f"[{grid_lon.min():.2f}, {grid_lat.min():.2f}] "
    f"to "
    f"[{grid_lon.max():.2f}, {grid_lat.max():.2f}]"
)


# ---------------------------------------------------------------------
# CLIMADA centroids
# ---------------------------------------------------------------------

centroids_odisha = Centroids(
    lat=grid_lat,
    lon=grid_lon,
)

assert len(centroids_odisha.lat) == 525
assert len(centroids_odisha.lon) == 525

print("CLIMADA centroids: 525")


# ---------------------------------------------------------------------
# P1-compatible distance-to-coast
# ---------------------------------------------------------------------

print(
    "\nCalculating WGS84 distance-to-coast "
    "from Natural Earth 10m coastline..."
)

shpfilename = shpreader.natural_earth(
    resolution="10m",
    category="physical",
    name="coastline",
)

coastline = gpd.read_file(shpfilename)

assert coastline.crs.to_epsg() == 4326
assert len(coastline) > 0

geod = Geod(ellps="WGS84")


def distance_to_coast(point, coastline):
    distances = coastline.geometry.distance(point)
    nearest_idx = distances.idxmin()

    coast_point = nearest_points(
        point,
        coastline.geometry.loc[nearest_idx],
    )[1]

    _, _, distance_m = geod.inv(
        point.x,
        point.y,
        coast_point.x,
        coast_point.y,
    )

    return distance_m


centroids_odisha.gdf["dist_coast"] = (
    centroids_odisha.gdf.geometry.apply(
        lambda point: distance_to_coast(
            point,
            coastline,
        )
    )
)

assert len(centroids_odisha.gdf["dist_coast"]) == 525
assert np.isfinite(
    centroids_odisha.gdf["dist_coast"].to_numpy()
).all()
assert (
    centroids_odisha.gdf["dist_coast"] >= 0
).all()

print(
    "Distance-to-coast calculated for",
    len(centroids_odisha.gdf),
    "centroids",
)

print(
    f"Minimum distance : "
    f"{centroids_odisha.gdf['dist_coast'].min():.2f} m"
)

print(
    f"Maximum distance : "
    f"{centroids_odisha.gdf['dist_coast'].max():.2f} m"
)


# ---------------------------------------------------------------------
# Stage 1 / Stage 2 historical track
# ---------------------------------------------------------------------

print(
    f"\n[2] Loading {EVENT_NAME} historical IBTrACS track..."
)

historical_track = load_historical_track(
    EVENT["sid"]
)

print(
    f"[2] Building {EVENT_NAME} T-{LEAD_TIME_H} BestTrack..."
)

from track_availability import InsufficientTrackHistory, record_unavailable, AvailabilityRecord

try:
    best_track, target_times, extension_info = build_best_track(
        historical_track, LANDFALL_TIME, LEAD_TIME_H
    )
except InsufficientTrackHistory as exc:
    print(f"\n{exc}\n")
    record_unavailable(REPO / "outputs" / "track_availability.csv",
                       AvailabilityRecord(
                           event=EVENT_NAME, lead_time_h=LEAD_TIME_H,
                           available=False,
                           reason="origin precedes IBTrACS record",
                           shortfall_h=exc.shortfall_h,
                           first_track_time=str(exc.first_available)[:16],
                           requested_origin=str(exc.requested_origin)[:16],
                           actual_lead_h=float("nan")))
    raise SystemExit(0)      # clean exit, the bash loop continues

print(
    f"BestTrack points : "
    f"{len(best_track.horizon_h)} "
    f"({best_track.i_landfall + 1} forecast + "
    f"{best_track.n_extension_points} extension)"
)

print(
    f"Lead time        : "
    f"{best_track.landfall_horizon_h:.0f} h"
)


# ---------------------------------------------------------------------
# Reference physical track for CLIMADA
# ---------------------------------------------------------------------

print(
    "\nPreparing physical reference track for Stage 3..."
)

reference_physical_track = build_reference_physical_track(
    historical_track,
    LANDFALL_TIME,
    LEAD_TIME_H,
    target_times=target_times,
)

print(
    "Reference physical fields :",
    "central_pressure,",
    "environmental_pressure,",
    "radius_max_wind",
)

print(
    f"Reference points           : "
    f"{len(reference_physical_track.time)}"
)

assert best_track.event_name == EVENT_NAME
assert best_track.horizon_h[0] == 0
assert int(best_track.landfall_horizon_h) == LEAD_TIME_H
assert best_track.i_landfall == LEAD_TIME_H
assert len(best_track.horizon_h) == (
    LEAD_TIME_H + 1 + best_track.n_extension_points
)
if LEAD_TIME_H > 0:
    assert best_track.n_extension_points > 0, (
        "Positive-lead panel has no post-landfall extension. Backward "
        "along-track members would end offshore and report zero loss. "
        "See section 6.5.9."
    )


# ---------------------------------------------------------------------
# Stage 2 ensemble
# ---------------------------------------------------------------------

models = load_error_models()

print(
    "\n[3] Generating "
    f"{N_MEMBERS}-member Stage 2 ensemble..."
)

ensemble = generate_ensemble(
    best_track,
    lead_time_h=LEAD_TIME_H,
    models=models,
    n_members=N_MEMBERS,
    base_seed=BASE_SEED,
)

assert ensemble.n_members == N_MEMBERS

n_accepted = int(
    np.count_nonzero(ensemble.accepted)
)

n_rejected = (
    ensemble.n_members
    - n_accepted
)

print(f"Members generated : {ensemble.n_members}")
print(f"Members accepted  : {n_accepted}")
print(f"Members rejected  : {n_rejected}")
print(
    f"Rejection fraction: "
    f"{ensemble.rejection_fraction:.4%}"
)

print(
    "Rejection reasons :",
    ensemble.rejection_reasons,
)

if n_rejected > 0:
    raise RuntimeError(
        f"{EVENT_NAME} T-{LEAD_TIME_H}: "
        f"Stage 2 rejected {n_rejected}/{N_MEMBERS} members. "
        f"Reasons: {ensemble.rejection_reasons}"
    )

assert n_accepted == N_MEMBERS

# ---------------------------------------------------------------------
# CLIMADA tracks
# ---------------------------------------------------------------------

print(
    "\n[4] Converting ensemble to CLIMADA TCTracks..."
)

tracks = to_climada_tctracks(
    ensemble,
    landfall_time=LANDFALL_TIME,
    accepted_only=True,
    reference_physical_track=reference_physical_track,
)

assert len(tracks.data) == N_MEMBERS

print(
    f"CLIMADA tracks : {len(tracks.data)}"
)


# ---------------------------------------------------------------------
# Stage 3 hazard calculation
# ---------------------------------------------------------------------

print(
    "\n[5] Running TropCyclone.from_tracks()..."
)

print(
    "    This is the benchmarked operation."
)

t0 = time.perf_counter()

hazard = TropCyclone.from_tracks(
    tracks,
    centroids_odisha,
)

elapsed = time.perf_counter() - t0


# ---------------------------------------------------------------------
# Hazard validation
# ---------------------------------------------------------------------

print("\n[6] Validating hazard output...")

assert len(hazard.event_name) == N_MEMBERS

assert hazard.intensity.shape[0] == N_MEMBERS
assert hazard.intensity.shape[1] == 525

intensity = hazard.intensity

assert np.isfinite(
    intensity.data
).all()

max_wind = float(
    intensity.max()
)

nonzero = int(
    np.count_nonzero(
        intensity.data
    )
)

assert max_wind > 0.0

print(
    f"Hazard events       : "
    f"{len(hazard.event_name)}"
)

print(
    f"Hazard centroids    : "
    f"{intensity.shape[1]}"
)

print(
    f"Intensity shape     : "
    f"{intensity.shape}"
)

print(
    f"Maximum wind        : "
    f"{max_wind:.4f} m/s"
)

print(
    f"Non-zero entries    : "
    f"{nonzero:,}"
)

print(
    f"Runtime             : "
    f"{elapsed:.3f} s"
)


# ---------------------------------------------------------------------
# Save
# ---------------------------------------------------------------------

print("\n[7] Saving benchmark hazard...")

with open(
    OUTPUT_FILE,
    "wb",
) as f:
    pickle.dump(
        hazard,
        f,
    )

print(
    f"Saved to: {OUTPUT_FILE}"
)

# ---------------------------------------------------------------------
# Record the post-landfall extension (section 6.5.9)
# ---------------------------------------------------------------------

EXTENSION_CSV = REPO / "outputs" / "track_extension.csv"

_ext_fields = list(extension_info.keys())
_ext_rows = []

if EXTENSION_CSV.exists():
    with EXTENSION_CSV.open("r", newline="", encoding="utf-8") as f:
        _ext_rows = [
            r for r in csv.DictReader(f)
            if not (
                r["event"] == EVENT_NAME
                and int(r["lead_time_h"]) == LEAD_TIME_H
            )
        ]

_ext_rows.append({k: extension_info[k] for k in _ext_fields})
_ext_rows.sort(key=lambda r: (str(r["event"]), -int(r["lead_time_h"])))

with EXTENSION_CSV.open("w", newline="", encoding="utf-8") as f:
    w = csv.DictWriter(f, fieldnames=_ext_fields)
    w.writeheader()
    w.writerows(_ext_rows)

print(f"Extension recorded in: {EXTENSION_CSV}")

# ---------------------------------------------------------------------
# Final summary
# ---------------------------------------------------------------------

print(
    "\n" + "=" * 72
)

print(
    f"STAGE 3 {EVENT_NAME.upper()} "
    f"T-{LEAD_TIME_H} BENCHMARK PASSED"
)

print(
    "=" * 72
)

print(
    f"Events              : "
    f"{len(hazard.event_name)}"
)

print(
    f"Centroids           : "
    f"{len(hazard.centroids.gdf)}"
)

print(
    f"Intensity shape     : "
    f"{hazard.intensity.shape}"
)

print(
    f"Runtime              : "
    f"{elapsed:.3f} s"
)

print(
    f"Output               : "
    f"{OUTPUT_FILE}"
)

print("=" * 72)