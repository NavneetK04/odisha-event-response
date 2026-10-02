"""
Project 3 - Cyclone Event Response
Stage 2: forecast-track ensemble generation.

What this implements
--------------------
The persistent-deviate scheme locked in FORECAST_ERROR_MODEL.md section 6.5.3:

    For member m of a forecast issued at lead time T:
        draw z_m = (z_a, z_c, z_v), each ~ N(0, 1), ONCE per member

    For each track point at forecast horizon h in [0, T]:
        along-track offset   a(h) = sigma_a(h) * z_a
        cross-track offset   c(h) = sigma_c(h) * z_c
        intensity offset     v(h) = sigma_v(h) * z_v
        rotate (a, c) into geographic coordinates by the local track heading

A member draws its wrongness once and keeps it. Correlation between the offsets
at any two points on a member's track is exactly 1; only the magnitude grows.

Why one ladder serves both the approach path and the landfall point
-------------------------------------------------------------------
Section 7 requires that general track error and landfall-point error never be
conflated. The section 6.5 calibration satisfies both simultaneously at each
horizon h:

    E[ sqrt(a^2 + c^2) ] = DPE(h)      <- 2-D position error, approach path
    E[ |c| ]             = LPE(h)      <- 1-D coastal displacement, landfall

so the total displacement along the approach path reproduces the published
track error, while the cross-track component alone reproduces the published
landfall error. At h = T the storm is at landfall by construction, so the
landfall spread of the ensemble is the published landfall statistic without
any further work. The separation is structural, not bookkeeping.

Horizons beyond 72 h (the cone phase, section 6.5.7) have no calibrated
landfall error, so they are drawn isotropically from DPE.

Coordinate handling
-------------------
Offsets are applied with a geodesic forward solution (pyproj.Geod, WGS84)
rather than a flat-earth degrees approximation. At 245 km offsets and 17-22N
the flat-earth error is small but it is free to avoid, and the project has
already been bitten once by a degree-space shortcut (Project 2, section 10.7).

Conventions
-----------
- horizon h : hours ahead of the forecast origin. h = 0 is the origin,
              h = T is landfall. Increasing h means further into the future.
- along-track : positive = ahead of the unperturbed position (storm running early)
- cross-track : positive = to the RIGHT of the direction of motion
- heading     : bearing in degrees clockwise from true north
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Mapping
import math
import hashlib
import numpy as np
from pyproj import Geod

import error_model as em


GEOD = Geod(ellps="WGS84")

# A member's position offset at the forecast origin is not zero: the origin
# position is known only as well as the analysis allows. That is exactly the
# T-0 analysis uncertainty in section 5.1, so the ladder is used as-is at h = 0.
ORIGIN_HORIZON_H = 0

# Plausibility bounds. Under a persistent-deviate scheme tracks are smooth
# deformations of the best track, so the rejection rate should be essentially
# zero. Any rejection is a signal that the perturbation scheme is wrong, not
# noise to tolerate. See the tightened threshold in MAX_REJECTION_FRACTION.
MIN_INTENSITY_KT = 5.0
MAX_TRANSLATION_SPEED_KMH = 120.0
MAX_REJECTION_FRACTION = 0.002

DEFAULT_SEED = 42

# ---------------------------------------------------------------------
# Track-relative frame
# ---------------------------------------------------------------------
# Section 6.5.3 rotates each member's offset into geographic coordinates by
# "the local track heading". Using the LOCAL heading at every point makes the
# frame rotate with the best track, and a member sitting far off that track
# gets swung by the rotation. At T-72h a 3-sigma member is roughly 530 km out,
# so a 13-degree heading change between two consecutive hourly points moves it
# 120 km in one hour. That is not the storm travelling; it is the frame
# turning, and the plausibility screen correctly rejects it as an implausible
# translation speed.
#
# The fix is to use ONE frame for the whole forecast, anchored to the heading
# on approach to landfall. Nothing is lost:
#
#   - The position-error MAGNITUDE (DPE) is frame-invariant, so the approach
#     path stays calibrated at every horizon.
#   - The along/cross SPLIT is only ever calibrated at landfall, where
#     sigma_c is the displacement along the coast. Anchoring the frame there
#     makes that the exact quantity IMD publishes.
#   - Frame rotation, and therefore the artefact, is removed entirely rather
#     than merely reduced.
#
# "local" is retained so the two can be compared and the difference reported.
FRAME_MODE = "landfall"
VALID_FRAME_MODES = ("landfall", "local", "smoothed")
SMOOTH_WINDOW_H = 6.0

# Post-landfall extension (section 6.5.9).
#
# The ensemble track must continue past landfall far enough that a member
# displaced BACKWARD along track still reaches the coast. 4 sigma leaves about
# 3 draws in 100,000 uncovered on the backward side, so at N = 500 the expected
# number of still-truncated members is under 0.02.
EXTENSION_SIGMAS = 4.0

# Over the extension the best track is decaying inland, so best-track wind plus
# a negative intensity deviate can reach zero or below. CLIMADA needs a positive
# intensity, so the extension is clipped at this floor. Inside the forecast
# window a below-floor intensity remains a rejection, not something to clip.
EXTENSION_WIND_FLOOR_KT = 10.0


# ---------------------------------------------------------------------
# Inputs and outputs
# ---------------------------------------------------------------------

@dataclass(frozen=True)
class BestTrack:
    """A best track segment from the forecast origin to landfall, plus an
    optional post-landfall extension.

    horizon_h must be strictly increasing and start at 0 (the forecast origin).
    landfall_index marks which point is landfall. Points after it are the
    extension, and they exist for one reason: so that a member displaced
    backward along track can still reach the coast.

    Without the extension the array ended AT landfall, so a backward along-track
    displacement left the member's last position offshore. CLIMADA then saw a
    track that never came ashore and the loss was exactly zero. That is not a
    cross-track miss, it is the array running out, and it is not symmetric: only
    backward displacements are affected. The signature is a zero-loss set that
    is one-sided in z_a, which is what outputs/zero_loss_mechanism.csv reports
    and what section 6.5.9 records.
    """

    event_name: str
    horizon_h: np.ndarray        # (P,) hours ahead of origin, ascending, [0] == 0
    lat: np.ndarray              # (P,) degrees north
    lon: np.ndarray              # (P,) degrees east
    max_wind_kt: np.ndarray      # (P,) knots
    landfall_index: int = -1     # -1 means the last point, i.e. no extension

    def __post_init__(self) -> None:
        p = len(self.horizon_h)
        if p < 1:
            raise ValueError("A best track segment needs at least one point.")
        for nm, arr in (("lat", self.lat), ("lon", self.lon), ("max_wind_kt", self.max_wind_kt)):
            if len(arr) != p:
                raise ValueError(f"{nm} has length {len(arr)}, expected {p}.")
        if not np.all(np.diff(self.horizon_h) > 0):
            raise ValueError("horizon_h must be strictly increasing.")
        if self.horizon_h[0] != ORIGIN_HORIZON_H:
            raise ValueError("horizon_h must start at 0 (the forecast origin).")
        if not np.all(np.isfinite(self.lat)) or not np.all(np.isfinite(self.lon)):
            raise ValueError("Best track contains non-finite coordinates.")
        if np.any(self.max_wind_kt <= 0):
            raise ValueError("Best track contains non-positive intensity.")
        if not (-1 <= int(self.landfall_index) < p):
            raise ValueError(
                f"landfall_index {self.landfall_index} is out of range for "
                f"{p} track points."
            )

    @property
    def i_landfall(self) -> int:
        """Index of the landfall point. Points beyond it are the extension."""
        p = len(self.horizon_h)
        return p - 1 if int(self.landfall_index) < 0 else int(self.landfall_index)

    @property
    def landfall_horizon_h(self) -> float:
        return float(self.horizon_h[self.i_landfall])

    @property
    def n_extension_points(self) -> int:
        return len(self.horizon_h) - 1 - self.i_landfall


@dataclass(frozen=True)
class Ensemble:
    event_name: str
    lead_time_h: int
    n_members: int
    seed: int

    horizon_h: np.ndarray        # (P,)
    lat: np.ndarray              # (M, P)
    lon: np.ndarray              # (M, P)
    max_wind_kt: np.ndarray      # (M, P)

    # Diagnostics, kept so the assertions can be checked without re-deriving.
    along_km: np.ndarray         # (M, P)
    cross_km: np.ndarray         # (M, P)
    sigma_a_km: np.ndarray       # (P,)
    sigma_c_km: np.ndarray       # (P,)
    sigma_v_kt: np.ndarray       # (P,)
    frame_heading_deg: np.ndarray  # (P,) the frame the offsets were applied in
    frame_mode: str
    z_a: np.ndarray              # (M,)
    z_c: np.ndarray              # (M,)
    z_v: np.ndarray              # (M,)

    accepted: np.ndarray         # (M,) bool
    i_landfall: int = -1         # -1 means the last column
    rejection_reasons: dict[str, int] = field(default_factory=dict)

    @property
    def rejection_fraction(self) -> float:
        return float((~self.accepted).mean())

    @property
    def i_lf(self) -> int:
        """Resolved landfall column. Every landfall quantity reads this, not -1.

        Before the post-landfall extension existed, -1 and landfall were the
        same column and the distinction did not matter. They are no longer the
        same column, and a stray -1 would silently report an inland quantity as
        a landfall one.
        """
        return self.lat.shape[1] - 1 if int(self.i_landfall) < 0 else int(self.i_landfall)

    @property
    def landfall_cross_km(self) -> np.ndarray:
        """Cross-track displacement at landfall: the landfall-point error."""
        return self.cross_km[:, self.i_lf]

    @property
    def landfall_along_km(self) -> np.ndarray:
        """Along-track displacement at landfall: the timing error, in distance."""
        return self.along_km[:, self.i_lf]


# ---------------------------------------------------------------------
# Sigma ladder
# ---------------------------------------------------------------------

def build_sigma_ladder(
    models: Mapping[int, "em.ForecastError"],
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Assemble (leads, sigma_a, sigma_c) over every horizon the model covers.

    Anisotropic where section 6.5 calibrated it (0 to 72 h). Isotropic beyond,
    drawn from DPE, per section 6.5.7: the fitted anisotropy ratio reaches ~1 by
    96 h and a ratio below 1 is unphysical, so isotropy is the honest treatment
    there rather than an extrapolation.
    """
    leads: list[int] = []
    sig_a: list[float] = []
    sig_c: list[float] = []

    for lead in sorted(set(em.LANDFALL_SAMPLER_LEADS_H) | set(em.TRACK_LADDER_H)):
        if lead in em.LANDFALL_SAMPLER_LEADS_H:
            a, c = em.landfall_sampler_parameters(lead)
        else:
            dpe = models[lead].track_error_km
            if not np.isfinite(dpe):
                raise ValueError(f"No track error available at {lead} h.")
            a = c = em.mean_error_to_rayleigh_sigma(dpe)
        leads.append(lead)
        sig_a.append(a)
        sig_c.append(c)

    leads_arr = np.asarray(leads, dtype=float)
    if not np.all(np.diff(leads_arr) > 0):
        raise AssertionError("Sigma ladder leads are not strictly increasing.")
    return leads_arr, np.asarray(sig_a), np.asarray(sig_c)


def interpolate_sigma(
    horizon_h: np.ndarray,
    ladder_leads: np.ndarray,
    ladder_values: np.ndarray,
) -> np.ndarray:
    """Linear interpolation in horizon between calibrated points.

    Linear is a stated assumption (section 6.5.8). It matters most across the
    0 to 12 h gap, where sigma_a rises roughly thirteenfold with no published
    data in between. Extrapolation beyond the ladder is refused rather than
    silently flat-extended.
    """
    h = np.asarray(horizon_h, dtype=float)
    if h.min() < ladder_leads[0] - 1e-9 or h.max() > ladder_leads[-1] + 1e-9:
        raise ValueError(
            f"Horizon range [{h.min()}, {h.max()}] falls outside the calibrated "
            f"ladder [{ladder_leads[0]}, {ladder_leads[-1]}]. Refusing to extrapolate."
        )
    return np.interp(h, ladder_leads, ladder_values)


def build_intensity_sigma(
    horizon_h: np.ndarray,
    models: Mapping[int, "em.ForecastError"],
) -> np.ndarray:
    """Intensity sampler width from published RMSE, not MAE.

    RMSE is the second-moment quantity the sampler needs; MAE understates the
    tail. The gap is a factor of ~2 at 72 h, which is where the loss fan is
    widest (section 3.3).
    """
    leads, vals = [], []
    for lead in sorted(models):
        rmse = models[lead].intensity_rmse_kt
        if rmse is not None and np.isfinite(rmse):
            leads.append(float(lead))
            vals.append(float(rmse))
    if not leads:
        raise ValueError("No intensity RMSE available in the error model.")

    # Intensity error at the origin is best-track intensity uncertainty, not
    # zero (section 5.2). The ladder starts at 12 h, so anchor h = 0 explicitly.
    if leads[0] > 0:
        leads.insert(0, 0.0)
        vals.insert(0, em.T0_INTENSITY_ERROR_KT)

    h = np.asarray(horizon_h, dtype=float)
    if h.max() > leads[-1] + 1e-9:
        raise ValueError(
            f"Horizon {h.max()} h exceeds the calibrated intensity ladder "
            f"(max {leads[-1]} h). Refusing to extrapolate intensity error."
        )
    return np.interp(h, np.asarray(leads), np.asarray(vals))


def interpolate_sigma_extended(
    horizon_h: np.ndarray,
    ladder_leads: np.ndarray,
    ladder_values: np.ndarray,
    i_landfall: int,
) -> np.ndarray:
    """Ladder interpolation up to landfall, held constant beyond it.

    The ladder is calibrated against forecast lead time and stops at 72 h, so
    there is nothing to interpolate past landfall and interpolate_sigma would
    refuse, correctly. Freezing sigma over the extension makes the post-landfall
    segment a RIGID TRANSLATION of the real track by each member's landfall
    offset, which is precisely what a persistent along-track error means: the
    storm arrives late, then follows its own path inland.

    The two alternatives are both wrong. Growing sigma past landfall fabricates
    forecast error for a horizon the forecast was never verified at. Shrinking
    it pulls late members back towards the best track and quietly undoes the
    displacement whose effect the panel exists to measure.
    """
    h = np.asarray(horizon_h, dtype=float)
    out = np.empty_like(h)
    cut = int(i_landfall) + 1
    out[:cut] = interpolate_sigma(h[:cut], ladder_leads, ladder_values)
    if cut < len(h):
        out[cut:] = out[cut - 1]
    return out


def build_intensity_sigma_extended(
    horizon_h: np.ndarray,
    models: Mapping[int, "em.ForecastError"],
    i_landfall: int,
) -> np.ndarray:
    """Intensity sampler width, held constant past landfall. See above."""
    h = np.asarray(horizon_h, dtype=float)
    out = np.empty_like(h)
    cut = int(i_landfall) + 1
    out[:cut] = build_intensity_sigma(h[:cut], models)
    if cut < len(h):
        out[cut:] = out[cut - 1]
    return out

# ---------------------------------------------------------------------
# Geometry
# ---------------------------------------------------------------------

def track_heading_deg(lat: np.ndarray, lon: np.ndarray) -> np.ndarray:
    """Bearing of motion at each point, degrees clockwise from north.

    Forward difference at the first point, backward at the last, centred
    elsewhere via the geodesic azimuth between neighbouring points.
    """
    p = len(lat)
    if p == 1:
    # T-0 has no track direction. Position uncertainty is isotropic,
    # so heading is mathematically irrelevant.
        return np.array([0.0])
    if p < 1:
        raise ValueError("Need at least one point.")
    az_seg, _, _ = GEOD.inv(lon[:-1], lat[:-1], lon[1:], lat[1:])  # (P-1,)
    az_seg = np.asarray(az_seg, dtype=float)

    heading = np.empty(p, dtype=float)
    heading[0] = az_seg[0]
    heading[-1] = az_seg[-1]
    if p > 2:
        # Circular mean of the two adjacent segment bearings.
        a = np.deg2rad(az_seg[:-1])
        b = np.deg2rad(az_seg[1:])
        heading[1:-1] = np.rad2deg(
            np.arctan2(0.5 * (np.sin(a) + np.sin(b)), 0.5 * (np.cos(a) + np.cos(b)))
        )
    return np.mod(heading, 360.0)


def track_frame_heading(
    lat: np.ndarray,
    lon: np.ndarray,
    mode: str = FRAME_MODE,
    horizon_h: np.ndarray | None = None,
    i_landfall: int = -1,
) -> np.ndarray:
    """Heading array defining the track-relative frame, one value per point.

    landfall : a single bearing for the whole forecast, taken on approach to
               landfall. No rotation, so no lever-arm artefact. Default.
    local    : the instantaneous heading at each point. The original scheme,
               retained for comparison; it produces the artefact above.
    smoothed : a centred circular mean over SMOOTH_WINDOW_H hours. Reduces the
               artefact without removing it.
    """
    if mode not in VALID_FRAME_MODES:
        raise ValueError(f"Unknown frame mode '{mode}'. Use one of {VALID_FRAME_MODES}.")

    local = track_heading_deg(lat, lon)
    p = len(local)
    i_lf = p - 1 if int(i_landfall) < 0 else int(i_landfall)

    if mode == "local" or p == 1:
        return local

    if mode == "landfall":
        # The bearing on approach to LANDFALL, not at the end of the array.
        # With a post-landfall extension the last point is well inland on a
        # decaying, often turning track, and using its heading would rotate
        # the frame for every member and move every panel.
        return np.full(p, float(local[i_lf]))

    h = np.arange(p, dtype=float) if horizon_h is None else np.asarray(horizon_h, float)
    rad = np.deg2rad(local)
    sin_s, cos_s = np.sin(rad), np.cos(rad)
    out = np.empty(p, dtype=float)
    for i in range(p):
        w = np.abs(h - h[i]) <= SMOOTH_WINDOW_H / 2.0
        out[i] = np.rad2deg(np.arctan2(sin_s[w].mean(), cos_s[w].mean()))
    return np.mod(out, 360.0)


def apply_offsets(
    lat: np.ndarray,
    lon: np.ndarray,
    heading_deg: np.ndarray,
    along_km: np.ndarray,
    cross_km: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Displace each point by (along, cross) in the track-relative frame.

    Cross-track positive is 90 degrees to the RIGHT of motion, so the offset
    bearing is heading + atan2(cross, along).
    """
    r_km = np.hypot(along_km, cross_km)
    bearing = np.mod(heading_deg + np.rad2deg(np.arctan2(cross_km, along_km)), 360.0)

    flat_lon = np.broadcast_to(lon, along_km.shape).ravel()
    flat_lat = np.broadcast_to(lat, along_km.shape).ravel()

    new_lon, new_lat, _ = GEOD.fwd(
        flat_lon, flat_lat, bearing.ravel(), r_km.ravel() * 1000.0
    )
    return (
        np.asarray(new_lat).reshape(along_km.shape),
        np.asarray(new_lon).reshape(along_km.shape),
    )


# ---------------------------------------------------------------------
# Plausibility filter
# ---------------------------------------------------------------------

def screen_members(
    lat: np.ndarray,
    lon: np.ndarray,
    wind_kt: np.ndarray,
    horizon_h: np.ndarray,
    i_landfall: int = -1,
) -> tuple[np.ndarray, dict[str, int]]:
    """Physical plausibility screen. Returns (accepted mask, reason counts).

    Screening covers the forecast window [origin, landfall] only. The
    post-landfall extension is real best-track data rigidly translated, so it
    has nothing to screen. Including it would reject nearly everything: the
    best track decays inland, so wind + dv drops under MIN_INTENSITY_KT for a
    large share of members and the intensity_below_floor rule would fire on
    them all. That would turn a zero-loss artifact into a rejection artifact.
    """
    m = lat.shape[0]
    i_lf = lat.shape[1] - 1 if int(i_landfall) < 0 else int(i_landfall)
    lat = lat[:, : i_lf + 1]
    lon = lon[:, : i_lf + 1]
    wind_kt = wind_kt[:, : i_lf + 1]
    horizon_h = np.asarray(horizon_h, dtype=float)[: i_lf + 1]
    accepted = np.ones(m, dtype=bool)
    reasons: dict[str, int] = {}

    def fail(mask: np.ndarray, label: str) -> None:
        nonlocal accepted
        newly = mask & accepted
        n = int(newly.sum())
        if n:
            reasons[label] = reasons.get(label, 0) + n
        accepted = accepted & ~mask

    fail(~np.isfinite(lat).all(axis=1) | ~np.isfinite(lon).all(axis=1), "non_finite_position")
    fail((np.abs(lat) > 89.0).any(axis=1), "pole_crossing")
    fail((wind_kt < MIN_INTENSITY_KT).any(axis=1), "intensity_below_floor")

    dt = np.diff(horizon_h)
    if np.any(dt <= 0):
        raise ValueError("horizon_h must be strictly increasing.")
    _, _, seg_m = GEOD.inv(lon[:, :-1], lat[:, :-1], lon[:, 1:], lat[:, 1:])
    speed_kmh = (np.asarray(seg_m) / 1000.0) / dt[None, :]
    fail((speed_kmh > MAX_TRANSLATION_SPEED_KMH).any(axis=1), "implausible_translation_speed")

    return accepted, reasons

def best_track_with_extension(
    event_name: str,
    horizon_h: np.ndarray,
    lat: np.ndarray,
    lon: np.ndarray,
    max_wind_kt: np.ndarray,
    lead_time_h: int,
    sigma_a_landfall_km: float,
    extension_sigmas: float = EXTENSION_SIGMAS,
) -> tuple[BestTrack, dict]:
    """Slice [origin, landfall + extension] out of a FULL best track.

    Pass the whole IBTrACS track with horizon_h measured in hours relative to
    the forecast origin, so horizon 0 is the origin and horizon lead_time_h is
    landfall. Points before the origin and after the extension are dropped.

    The extension is sized in PATH LENGTH, not hours, because what has to be
    covered is a displacement in kilometres: extension_sigmas times the
    along-track sigma at landfall. A storm moving slowly needs more hours of
    track to cover the same distance, and hours would get that backwards.

    The walk stops early if the best track ends, or if intensity reaches zero.
    An event whose record ends soon after landfall cannot be fully covered, and
    the returned dict says so. That shortfall is a documented limitation, not a
    silent one: validate_ensemble asserts on it.
    """
    h = np.asarray(horizon_h, dtype=float)
    lat = np.asarray(lat, dtype=float)
    lon = np.asarray(lon, dtype=float)
    wind = np.asarray(max_wind_kt, dtype=float)

    i_origin = int(np.argmin(np.abs(h - 0.0)))
    if abs(h[i_origin]) > 1e-6:
        raise ValueError(
            "Full track contains no point at horizon 0 (the forecast origin)."
        )
    i_lf = int(np.argmin(np.abs(h - float(lead_time_h))))
    if abs(h[i_lf] - float(lead_time_h)) > 1e-6:
        raise ValueError(
            f"Full track contains no point at horizon {lead_time_h} h (landfall)."
        )

    required_km = float(extension_sigmas) * float(sigma_a_landfall_km)
    acc_km, i_end = 0.0, i_lf
    while (
        i_end + 1 < len(h)
        and acc_km < required_km
        and np.isfinite(wind[i_end + 1])
        and wind[i_end + 1] > 0
    ):
        _, _, d_m = GEOD.inv(lon[i_end], lat[i_end], lon[i_end + 1], lat[i_end + 1])
        acc_km += float(d_m) / 1000.0
        i_end += 1

    info = dict(
        event=event_name,
        lead_time_h=int(lead_time_h),
        sigma_a_landfall_km=float(sigma_a_landfall_km),
        required_extension_km=required_km,
        available_extension_km=acc_km,
        extension_points=int(i_end - i_lf),
        hours_past_landfall=float(h[i_end] - h[i_lf]),
        fully_covered=bool(acc_km >= required_km),
    )

    seg = BestTrack(
        event_name=event_name,
        horizon_h=h[i_origin : i_end + 1] - h[i_origin],
        lat=lat[i_origin : i_end + 1],
        lon=lon[i_origin : i_end + 1],
        max_wind_kt=wind[i_origin : i_end + 1],
        landfall_index=int(i_lf - i_origin),
    )
    return seg, info

# ---------------------------------------------------------------------
# The generator
# ---------------------------------------------------------------------

def member_seed(base_seed: int, event_name: str, lead_time_h: int) -> int:
    """Stable deterministic seed for one (event, lead-time) ensemble.

    Uses SHA-256 rather than Python's built-in hash(), because Python
    randomizes hash() between interpreter processes.
    """

    key = (
        f"{base_seed}|{event_name}|{int(lead_time_h)}"
    ).encode("utf-8")

    digest = hashlib.sha256(key).digest()

    # Use the first 4 bytes to create a stable 32-bit seed.
    seed = int.from_bytes(
        digest[:4],
        byteorder="big",
        signed=False,
    )

    # Keep the seed inside NumPy's positive 32-bit range.
    return seed % (2**32 - 1)


def generate_ensemble(
    best_track: BestTrack,
    lead_time_h: int,
    models: Mapping[int, "em.ForecastError"],
    n_members: int = 500,
    base_seed: int = DEFAULT_SEED,
    latent_draws: tuple[np.ndarray, np.ndarray, np.ndarray] | None = None,
    frame_mode: str = FRAME_MODE,
) -> Ensemble:
    """Generate an N-member perturbed-track ensemble for one forecast issue."""
    if n_members <= 0:
        raise ValueError("n_members must be positive.")
    if abs(best_track.landfall_horizon_h - lead_time_h) > 1e-6:
        raise ValueError(
            f"Best track landfall point sits at horizon "
            f"{best_track.landfall_horizon_h} h but the forecast lead time is "
            f"{lead_time_h} h. landfall_index must point at the landfall point, "
            "not at the end of the array."
        )
    if lead_time_h == 0 and len(best_track.horizon_h) != 1:
            raise ValueError("T-0 BestTrack must contain exactly one point.")

    if lead_time_h > 0 and len(best_track.horizon_h) < 2:
        raise ValueError("Positive-lead BestTrack must contain at least two points.")

    horizon = np.asarray(best_track.horizon_h, dtype=float)
    p = len(horizon)

    i_lf = best_track.i_landfall

    leads, lad_a, lad_c = build_sigma_ladder(models)
    sigma_a = interpolate_sigma_extended(horizon, leads, lad_a, i_lf)
    sigma_c = interpolate_sigma_extended(horizon, leads, lad_c, i_lf)
    sigma_v = build_intensity_sigma_extended(horizon, models, i_lf)

    if latent_draws is None:
        rng = np.random.default_rng(member_seed(base_seed, best_track.event_name, lead_time_h))
        z_a = rng.normal(size=n_members)
        z_c = rng.normal(size=n_members)
        z_v = rng.normal(size=n_members)
    else:
        z_a, z_c, z_v = (np.asarray(z, dtype=float) for z in latent_draws)
        if not (len(z_a) == len(z_c) == len(z_v) == n_members):
            raise ValueError("latent_draws must each have length n_members.")

    # One draw per member, applied at every horizon: correlation 1 by construction.
    along = z_a[:, None] * sigma_a[None, :]
    cross = z_c[:, None] * sigma_c[None, :]
    dv = z_v[:, None] * sigma_v[None, :]

    heading = track_frame_heading(
        best_track.lat, best_track.lon, frame_mode,
        horizon_h=horizon, i_landfall=i_lf,
    )
    lat, lon = apply_offsets(best_track.lat, best_track.lon, heading, along, cross)
    wind = best_track.max_wind_kt[None, :] + dv

    # Clip the extension only. See EXTENSION_WIND_FLOOR_KT.
    if i_lf + 1 < wind.shape[1]:
        np.clip(
            wind[:, i_lf + 1 :], EXTENSION_WIND_FLOOR_KT, None,
            out=wind[:, i_lf + 1 :],
        )

    accepted, reasons = screen_members(lat, lon, wind, horizon, i_landfall=i_lf)

    return Ensemble(
        event_name=best_track.event_name,
        lead_time_h=int(lead_time_h),
        n_members=n_members,
        seed=member_seed(base_seed, best_track.event_name, lead_time_h),
        horizon_h=horizon,
        lat=lat,
        lon=lon,
        max_wind_kt=wind,
        along_km=along,
        cross_km=cross,
        sigma_a_km=sigma_a,
        sigma_c_km=sigma_c,
        sigma_v_kt=sigma_v,
        frame_heading_deg=heading,
        frame_mode=frame_mode,
        z_a=z_a,
        z_c=z_c,
        z_v=z_v,
        accepted=accepted,
        i_landfall=int(i_lf),
        rejection_reasons=reasons
    )


# ---------------------------------------------------------------------
# Executable assertions
# ---------------------------------------------------------------------

def validate_ensemble(
    ens: Ensemble,
    best_track: BestTrack,
    models: Mapping[int, "em.ForecastError"],
    verbose: bool = True,
) -> dict[str, float]:
    """Assert every property section 6.5 and blueprint section 8.2 require.

    Every Monte Carlo tolerance is derived from N, never chosen. A tolerance
    picked as a round number is how a correct sampler gets debugged for a day.
    """

    m, p = ens.lat.shape
    i_lf = ens.i_lf
    n_ext = p - 1 - i_lf
    diag: dict[str, float] = {"extension_points": float(n_ext)}
    if i_lf > 0:
        nz = np.flatnonzero(ens.sigma_a_km[: i_lf + 1] > 0)
        i, j = int(nz[0]), int(nz[-1])
        r_a = float(np.corrcoef(ens.along_km[:, i], ens.along_km[:, j])[0, 1])
        r_c = float(np.corrcoef(ens.cross_km[:, i], ens.cross_km[:, j])[0, 1])
        assert r_a > 1 - 1e-10, f"Along-track offsets not persistent along the track (r={r_a})."
        assert r_c > 1 - 1e-10, f"Cross-track offsets not persistent along the track (r={r_c})."

    # --- growth: sigma must be non-decreasing AND actually grow ------------
    # Non-decreasing alone is not enough: a rigid translation (sigma frozen at
    # the landfall value) satisfies it, and that is the second failure mode in
    # section 6.5.1. Pin both ends to the ladder instead.
    
    assert np.all(np.diff(ens.sigma_a_km[: i_lf + 1]) >= -1e-9), (
        "sigma_a decreases along the forecast."
    )
    assert np.all(np.diff(ens.sigma_c_km[: i_lf + 1]) >= -1e-9), (
        "sigma_c decreases along the forecast."
    )

    # Over the extension sigma must be EXACTLY flat, not merely non-decreasing.
    # Flat makes the post-landfall segment a rigid translation by the landfall
    # offset. Any drift there would either fabricate forecast error past the
    # calibrated horizon or pull late members back onto the best track.
    if n_ext:
        assert np.allclose(ens.sigma_a_km[i_lf:], ens.sigma_a_km[i_lf], atol=1e-9), (
            "sigma_a is not constant across the post-landfall extension."
        )
        assert np.allclose(ens.sigma_c_km[i_lf:], ens.sigma_c_km[i_lf], atol=1e-9), (
            "sigma_c is not constant across the post-landfall extension."
        )

    a0, c0 = em.landfall_sampler_parameters(ORIGIN_HORIZON_H)
    assert abs(ens.sigma_a_km[0] - a0) < 1e-6, (
        f"sigma_a at the forecast origin is {ens.sigma_a_km[0]:.4f} km, expected the "
        f"T-0 analysis value {a0:.4f} km. Error must grow from the origin, not be "
        "applied rigidly."
    )
    assert abs(ens.sigma_c_km[0] - c0) < 1e-6, (
        f"sigma_c at the forecast origin is {ens.sigma_c_km[0]:.4f} km, expected {c0:.4f} km."
    )
    if T_ := ens.lead_time_h:
        if T_ in em.LANDFALL_SAMPLER_LEADS_H:
            aT, cT = em.landfall_sampler_parameters(T_)
            assert abs(ens.sigma_a_km[i_lf] - aT) < 1e-6, (
                f"sigma_a at landfall is {ens.sigma_a_km[i_lf]:.4f} km, expected the "
                f"lead-{T_}h value {aT:.4f} km."
            )
            assert abs(ens.sigma_c_km[-1] - cT) < 1e-6, (
                f"sigma_c at landfall is {ens.sigma_c_km[-1]:.4f} km, expected {cT:.4f} km."
            )
    if i_lf > 0:
        assert ens.sigma_a_km[i_lf] > ens.sigma_a_km[0], (
            "sigma_a does not grow along the forecast."
        )
        assert ens.sigma_c_km[i_lf] > ens.sigma_c_km[0], (
            "sigma_c does not grow along the forecast."
        )
    # --- unbiasedness: ensemble mean position is the best track ------------
    # Tolerance from the Monte Carlo SE of the mean offset, not a round number.
    mean_along = np.abs(ens.along_km.mean(axis=0))
    mean_cross = np.abs(ens.cross_km.mean(axis=0))
    se_along = ens.sigma_a_km / math.sqrt(m)
    se_cross = ens.sigma_c_km / math.sqrt(m)
    tol_a = np.maximum(4.0 * se_along, 1e-9)
    tol_c = np.maximum(4.0 * se_cross, 1e-9)
    assert np.all(mean_along <= tol_a), (
        f"Ensemble is biased along-track: max |mean| {mean_along.max():.3f} km "
        f"against 4 SE {tol_a.max():.3f} km."
    )
    assert np.all(mean_cross <= tol_c), (
        f"Ensemble is biased cross-track: max |mean| {mean_cross.max():.3f} km "
        f"against 4 SE {tol_c.max():.3f} km."
    )

    # --- landfall spread reproduces the published LPE ----------------------
    T = ens.lead_time_h
    if T in em.EXPECTED_LANDFALL_POINT_KM:
        lpe_pub = em.EXPECTED_LANDFALL_POINT_KM[T]

        # The published landfall figure is a 1-D coastal displacement at
        # 12-72 h, so it is compared against the cross-track component alone.
        # At T-0 it is NOT: the 5 km is analysis uncertainty about where the
        # centre is, which is a 2-D isotropic position error with no preferred
        # direction. Comparing it against the cross component would expect
        # 5 km from a quantity whose mean is 5 x sqrt(2/pi) = 3.19 km, and the
        # check would fail on a sampler that is behaving exactly as specified.
        if T in em.ISOTROPIC_LANDFALL_LEADS_H:
            mag = np.hypot(ens.landfall_along_km, ens.landfall_cross_km)
            sampled = float(mag.mean())
            se = float(mag.std(ddof=1)) / math.sqrt(m)
            basis = "2-D magnitude"
        else:
            sampled = float(np.mean(np.abs(ens.landfall_cross_km)))
            sc = ens.sigma_c_km[i_lf]
            se = sc * math.sqrt(1.0 - 2.0 / math.pi) / math.sqrt(m)
            basis = "1-D cross-track"

        diag["landfall_point_error_km"] = sampled
        assert abs(sampled - lpe_pub) <= 3 * se, (
            f"Landfall spread {sampled:.3f} km ({basis}) does not match "
            f"published {lpe_pub:.3f} km at {T} h (3 SE = {3*se:.3f})."
        )

    # --- total position spread reproduces the published DPE ----------------
    if T in em.EXPECTED_TRACK_KM:
        dpe_pub = em.EXPECTED_TRACK_KM[T]
        mag = np.hypot(ens.along_km[:, i_lf], ens.cross_km[:, i_lf])
        sampled = float(mag.mean())
        se = float(mag.std(ddof=1)) / math.sqrt(m)
        diag["position_error_km"] = sampled
        assert abs(sampled - dpe_pub) <= 3 * se, (
            f"Position spread {sampled:.3f} km does not match published DPE "
            f"{dpe_pub:.3f} km at {T} h (3 SE = {3*se:.3f})."
        )

    # --- the extension covers the displacement actually drawn ---------------
    # This is the regression test for the truncation defect. If the extension
    # is shorter in path length than the most backward member drawn, that
    # member's track still ends offshore, CLIMADA never brings it to the coast,
    # and its zero loss is an artifact of array length rather than a statement
    # about where the storm might have gone.
    #
    # The defect was invisible for a long time because a zero loss is a
    # perfectly ordinary value. It was found by reconstructing the latent draws
    # and noticing the zeros were one-sided in z_a, up to 14 sigma. An assertion
    # is cheaper than a diagnostic, so it lives here now.
    if n_ext:
        _, _, seg_m = GEOD.inv(
            best_track.lon[i_lf:-1], best_track.lat[i_lf:-1],
            best_track.lon[i_lf + 1 :], best_track.lat[i_lf + 1 :],
        )
        extension_km = float(np.sum(np.asarray(seg_m)) / 1000.0)
        worst_backward_km = float(max(0.0, -ens.landfall_along_km.min()))
        diag["extension_km"] = extension_km
        diag["worst_backward_km"] = worst_backward_km
        assert extension_km >= worst_backward_km, (
            f"Post-landfall extension is {extension_km:.1f} km but the most "
            f"backward member is displaced {worst_backward_km:.1f} km along "
            "track. Members beyond the extension never reach the coast, so "
            "their zero loss is an artifact, not a miss. Lengthen the "
            "extension or raise EXTENSION_SIGMAS (section 6.5.9)."
        )
    elif ens.lead_time_h > 0:
        raise AssertionError(
            f"{ens.event_name} T-{ens.lead_time_h}h has no post-landfall "
            "extension. Backward along-track members will end offshore and "
            "report zero loss. See section 6.5.9."
        )

    # --- geometry round-trip: recovered offsets match the ones applied -----
    # Guards the heading convention and the geodesic call. If cross-track were
    # applied to the left instead of the right, or along/cross were swapped,
    # this is what would catch it.
    az, _, dist_m = GEOD.inv(
        np.broadcast_to(best_track.lon, ens.lon.shape).ravel(),
        np.broadcast_to(best_track.lat, ens.lat.shape).ravel(),
        ens.lon.ravel(),
        ens.lat.ravel(),
    )
    # Use the frame the ensemble actually used. Recomputing it here would let
    # the generator and the check drift apart, and the round-trip would then
    # fail for a reason that has nothing to do with the geometry.
    heading = ens.frame_heading_deg
    rel = np.deg2rad(np.mod(np.asarray(az).reshape(ens.lon.shape)
                            - np.broadcast_to(heading, ens.lon.shape), 360.0))
    r_km = np.asarray(dist_m).reshape(ens.lon.shape) / 1000.0
    rec_along = r_km * np.cos(rel)
    rec_cross = r_km * np.sin(rel)
    max_err = float(np.max(np.abs(rec_along - ens.along_km) + np.abs(rec_cross - ens.cross_km)))
    diag["geometry_roundtrip_max_km"] = max_err
    assert max_err < 0.5, (
        f"Offset round-trip mismatch {max_err:.4f} km: the heading convention or "
        "the geodesic displacement is wrong."
    )

    # --- intensity sampler uses RMSE and stays physical --------------------
    assert np.all(ens.sigma_v_kt > 0), "Intensity sigma must be positive everywhere."
    if 72 in models and models[72].intensity_rmse_kt is not None:
        assert ens.sigma_v_kt.max() <= models[72].intensity_rmse_kt + 1e-9, (
            "Intensity sigma exceeds the published 72 h RMSE."
        )

    # --- plausibility filter: rejections are a signal, not noise -----------
    diag["rejection_fraction"] = ens.rejection_fraction
    assert ens.rejection_fraction <= MAX_REJECTION_FRACTION, (
        f"Rejection fraction {ens.rejection_fraction:.4f} exceeds "
        f"{MAX_REJECTION_FRACTION}. Under a persistent-deviate scheme tracks are "
        f"smooth by construction, so this means the perturbation scheme is wrong, "
        f"not that the filter is strict. Reasons: {ens.rejection_reasons}"
    )

    # --- reported, not asserted: implied landfall timing error -------------
    # Independent of everything used to calibrate sigma_a and sigma_c. A storm
    # whose translation speed differs from the ~23 km/h population average will
    # not match the published LTE, so this is a diagnostic rather than a test.
    if p > 1:
        _, _, seg_m = GEOD.inv(
        best_track.lon[:-1],
        best_track.lat[:-1],
        best_track.lon[1:],
        best_track.lat[1:],
    )
        total_km = float(np.sum(seg_m) / 1000.0)
        span_h = float(best_track.horizon_h[-1] - best_track.horizon_h[0])

        if span_h > 0:
            v_kmh = total_km / span_h
            diag["storm_translation_speed_kmh"] = v_kmh

            if v_kmh > 0:
                diag["implied_landfall_time_error_h"] = (
                    float(np.mean(np.abs(ens.landfall_along_km))) / v_kmh
                )
    if verbose:
        print(f"  Ensemble {ens.event_name} T-{ens.lead_time_h}h: "
              f"M={m}, P={p}, rejected={ens.rejection_fraction:.4%}")
        for k, v in diag.items():
            print(f"    {k:38s} {v:10.4f}")

    return diag


# ---------------------------------------------------------------------
# CLIMADA handoff (lazy import: the geometry above is testable without it)
# ---------------------------------------------------------------------

def to_climada_tctracks(
    ens: Ensemble,
    landfall_time: "np.datetime64",
    accepted_only: bool = True,
    reference_physical_track=None,
):
    """Convert to a CLIMADA TCTracks object for Stage 3.

    Imported lazily so that everything above can be unit-tested in an
    environment without CLIMADA installed.

    For T-0, Stage 2 contains a genuine one-point ensemble at landfall.
    CLIMADA's H08 windfield implementation requires at least two track
    positions because translational speed and pressure change are computed
    using backward differences. Therefore, T-0 receives one immediately
    preceding physical Fani point solely as computational context for the
    CLIMADA handoff. The uncertainty itself remains applied only to the
    T-0 endpoint.
    """
    import xarray as xr
    from climada.hazard import TCTracks

    idx = (
        np.flatnonzero(ens.accepted)
        if accepted_only
        else np.arange(ens.n_members)
    )

    # ------------------------------------------------------------------
    # Detect genuine T-0 ensemble
    # ------------------------------------------------------------------
    is_t0 = (
        len(ens.horizon_h) == 1
        and ens.lead_time_h == 0
    )

    if reference_physical_track is None:
        raise ValueError(
            "reference_physical_track is required for the "
            "Stage 3 CLIMADA handoff."
        )

    # ------------------------------------------------------------------
    # T-0: provide one preceding physical point as CLIMADA context.
    #
    # Stage 2 remains genuinely one-point. Only this computational
    # handoff gets the preceding physical state required by H08.
    # ------------------------------------------------------------------
    if is_t0:
        if len(reference_physical_track.time) < 2:
            raise ValueError(
                "T-0 CLIMADA handoff requires at least two "
                "reference physical points."
            )

        reference_for_handoff = reference_physical_track.isel(
            time=slice(
                len(reference_physical_track.time) - 2,
                None,
            )
        )

        times = reference_for_handoff.time.values

        step_h = float(
            (
                times[1] - times[0]
            ) / np.timedelta64(1, "h")
        )

        if step_h <= 0:
            raise ValueError(
                "T-0 physical context interval must be positive."
            )

    # ------------------------------------------------------------------
    # Normal forecast leads: existing behaviour
    # ------------------------------------------------------------------
    else:
        # Anchored on the LANDFALL horizon, not the end of the array. With a
        # post-landfall extension the last horizon is hours past landfall, and
        # subtracting it would stamp the entire track that many hours early.
        # The hazard would still compute and nothing would complain.
        hours_before_landfall = (
            ens.horizon_h - ens.horizon_h[ens.i_lf]
        )

        times = landfall_time + (
            hours_before_landfall * 3600e9
        ).astype("timedelta64[ns]")

        step_h = (
            float(np.diff(ens.horizon_h)[0])
            if len(ens.horizon_h) > 1
            else 1.0
        )

        reference_for_handoff = reference_physical_track

    # ------------------------------------------------------------------
    # Validate physical reference fields
    # ------------------------------------------------------------------
    required_vars = [
        "central_pressure",
        "environmental_pressure",
        "radius_max_wind",
    ]

    for var in required_vars:
        if var not in reference_for_handoff:
            raise ValueError(
                f"Reference physical track is missing '{var}'."
            )

    if len(reference_for_handoff.time) != len(times):
        raise ValueError(
            "Reference physical track length does not match "
            "the CLIMADA handoff timestamps."
        )

    if not np.array_equal(
        reference_for_handoff.time.values,
        times,
    ):
        raise ValueError(
            "Reference physical track timestamps do not match "
            "the CLIMADA handoff timestamps."
        )

    # ------------------------------------------------------------------
    # Build CLIMADA tracks
    # ------------------------------------------------------------------
    data = []

    for m in idx:

        # --------------------------------------------------------------
        # T-0:
        #
        # Point 0 = immediately preceding physical Fani point
        # Point 1 = perturbed T-0 landfall point
        #
        # Only the endpoint receives Stage-2 uncertainty.
        # --------------------------------------------------------------
        if is_t0:
            lat_values = np.array(
                [
                    reference_for_handoff.lat.values[0],
                    ens.lat[m, 0],
                ],
                dtype=float,
            )

            lon_values = np.array(
                [
                    reference_for_handoff.lon.values[0],
                    ens.lon[m, 0],
                ],
                dtype=float,
            )

            wind_values = np.array(
                [
                    reference_for_handoff.max_sustained_wind.values[0],
                    ens.max_wind_kt[m, 0],
                ],
                dtype=float,
            )

        # --------------------------------------------------------------
        # T-72 / T-48 / T-24 / T-12:
        # Existing Stage-2 ensemble is passed through unchanged.
        # --------------------------------------------------------------
        else:
            lat_values = np.asarray(
                ens.lat[m],
                dtype=float,
            )

            lon_values = np.asarray(
                ens.lon[m],
                dtype=float,
            )

            wind_values = np.asarray(
                ens.max_wind_kt[m],
                dtype=float,
            )

        ds = xr.Dataset(
            {
                "lat": (
                    "time",
                    lat_values,
                ),

                "lon": (
                    "time",
                    lon_values,
                ),

                "max_sustained_wind": (
                    "time",
                    wind_values,
                ),

                "time_step": (
                    "time",
                    np.full(
                        len(times),
                        step_h,
                        dtype=float,
                    ),
                ),

                "central_pressure": (
                    "time",
                    np.asarray(
                        reference_for_handoff.central_pressure.values,
                        dtype=float,
                    ),
                ),

                "environmental_pressure": (
                    "time",
                    np.asarray(
                        reference_for_handoff.environmental_pressure.values,
                        dtype=float,
                    ),
                ),

                "radius_max_wind": (
                    "time",
                    np.asarray(
                        reference_for_handoff.radius_max_wind.values,
                        dtype=float,
                    ),
                ),

                "basin": (
                    "time",
                    np.full(
                        len(times),
                        "NI",
                        dtype="U2",
                    ),
                ),
            },

            coords={
                "time": times,
            },

            attrs={
                "max_sustained_wind_unit": "kn",
                "central_pressure_unit": "mb",
                "sid": (
                    f"{ens.event_name}_T"
                    f"{ens.lead_time_h}_m{int(m):04d}"
                ),
                "name": (
                    f"{ens.event_name}_T"
                    f"{ens.lead_time_h}_m{int(m):04d}"
                ),
                "orig_event_flag": False,
                "category": None,
            },
        )

        data.append(ds)

    tracks = TCTracks()
    tracks.data = data

    return tracks