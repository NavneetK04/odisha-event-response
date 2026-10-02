"""
Project 3 — Cyclone Event Response
Forecast-error model and executable assertions.

Purpose
-------
Load the fixed IMD 2020–2024 forecast-error model, validate every
required quantity, and expose runtime parameters for Stage 2.

Important modelling distinction
--------------------------------
- General track/DPE error -> approach-path uncertainty
- Landfall-point error   -> landfall-point uncertainty

These are deliberately kept separate.

The CSV contains audit columns for Rayleigh sigma. The sampler does
NOT trust those columns: sigma is recalculated at runtime and the
stored values are checked against the calculation.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import math
import numpy as np
import pandas as pd


# ---------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------

PROJECT_ROOT = Path(__file__).resolve().parents[1]
ERROR_CSV = PROJECT_ROOT / "data" / "imd_forecast_errors.csv"


# ---------------------------------------------------------------------
# Locked model ladders
# ---------------------------------------------------------------------

TRACK_LADDER_H = (12, 24, 48, 72, 96, 120)
LANDFALL_LADDER_H = (0, 12, 24, 48, 72)

# 36 h and 60 h are retained in the CSV as evidence of the growth curve,
# but are NOT scenario lead times.
LANDFALL_SUPPORTING_H = (36, 60)


# ---------------------------------------------------------------------
# Tolerances
# ---------------------------------------------------------------------

VALUE_TOL_KM = 0.02
VALUE_TOL_H = 0.02
SIGMA_TOL_KM = 0.01

# Worst observed residual between hoyt_mean_magnitude(sigma_a, sigma_c) and the
# published DPE is 0.004 km. A 1.0 km tolerance would admit a sigma_a wrong by
# roughly 1 per cent.
DPE_CONSISTENCY_TOL_KM = 0.01

# Single seed for every seeded check in this module.
SAMPLER_SEED = 42

# Published/reconstructed values can be represented at different
# rounding levels in the CSV, so the 12 h RMSE check allows the
# published 5.0 kt representation of the pooled 4.960727 kt result.
INTENSITY_TOL_KT = 0.05


# ---------------------------------------------------------------------
# Expected locked values
# ---------------------------------------------------------------------

EXPECTED_LANDFALL_POINT_KM = {
    0: 5.0,
    12: 10.53,
    24: 16.19,
    48: 39.32,
    72: 69.52,
}

EXPECTED_LANDFALL_TIME_H = {
    0: 0.0,
    12: 1.85,
    24: 2.89,
    48: 4.19,
    72: 7.45,
}

EXPECTED_TRACK_KM = {
    12: 44.866834,
    24: 71.969578,
    48: 110.943644,
    72: 154.028889,
    96: 182.694545,
    120: 245.225000,
}

EXPECTED_TRACK_N = {
    12: 398,
    24: 332,
    48: 236,
    72: 135,
    96: 55,
    120: 28,
}

EXPECTED_INTENSITY_MAE_KT = {
    12: 3.384171,
    24: 5.9,
    48: 8.3,
    72: 9.8,
}

EXPECTED_INTENSITY_RMSE_KT = {
    12: 4.960727,
    24: 7.9,
    48: 11.0,
    72: 19.2,
}


# ---------------------------------------------------------------------
# Runtime representation
# ---------------------------------------------------------------------

@dataclass(frozen=True)
class ForecastError:
    lead_time_h: int

    # General track / approach-path
    track_error_km: float
    track_sigma_km: float
    track_n: int | None
    track_provenance: str

    # Landfall point
    landfall_point_error_km: float | None
    landfall_sigma_km: float | None
    landfall_sigma_a_km: float | None
    landfall_sigma_c_km: float | None
    landfall_n: int | None
    landfall_provenance: str

    # Landfall timing
    landfall_time_error_h: float | None

    # Intensity
    intensity_mae_kt: float | None
    intensity_rmse_kt: float | None
    intensity_n: int | None
    intensity_provenance: str


# ---------------------------------------------------------------------
# Mathematical helpers
# ---------------------------------------------------------------------

def mean_error_to_rayleigh_sigma(mean_error: float) -> float:
    """
    Convert mean Rayleigh error magnitude to Rayleigh sigma.

    For R ~ Rayleigh(sigma):

        E[R] = sigma * sqrt(pi / 2)

    Therefore:

        sigma = E[R] * sqrt(2 / pi)
    """
    if mean_error < 0:
        raise ValueError("Mean error cannot be negative.")

    return mean_error * math.sqrt(2.0 / math.pi)


# ---------------------------------------------------------------------
# §6.5 Anisotropic landfall-error model
# ---------------------------------------------------------------------

ANISOTROPIC_LANDFALL_LEADS_H = (12, 24, 48, 72)

# T-0 is analysis uncertainty, not forecast error. The anisotropy in §6.5
# exists because along-track error IS timing error; at analysis time there is
# no trajectory to be early or late along, so centre-fixing ignorance has no
# preferred direction. The ratio therefore drops to 1 discontinuously at T-0
# rather than decaying to it. This is physical, not an extrapolation.
ISOTROPIC_LANDFALL_LEADS_H = (0,)
T0_POSITION_ERROR_KM = 5.0

# Section 5.2: best-track intensity uncertainty at analysis time. A judgement,
# below the published 24 h MAE of 5.9 kt, not a published statistic.
T0_INTENSITY_ERROR_KT = 5.0
_T0_SIGMA_KM = T0_POSITION_ERROR_KM * math.sqrt(2.0 / math.pi)   # 3.9894 km

# Every lead in LANDFALL_LADDER_H must appear here or the sampler has a hole.
LANDFALL_SAMPLER_LEADS_H = ISOTROPIC_LANDFALL_LEADS_H + ANISOTROPIC_LANDFALL_LEADS_H

# Locked §6.5 source-calibrated component scales.
EXPECTED_LANDFALL_SIGMA_A_KM = {0: _T0_SIGMA_KM, 12: 52.40, 24: 84.47, 48: 120.25, 72: 154.00}
EXPECTED_LANDFALL_SIGMA_C_KM = {0: _T0_SIGMA_KM, 12: 13.20, 24: 20.29, 48: 49.28, 72: 87.13}

# Angular quadrature avoids the origin cusp of Cartesian Gauss-Hermite.
HOYT_ANGULAR_QUADRATURE_ORDER = 256

def landfall_sigma_c_from_lpe(landfall_point_error_km: float) -> float:
    if landfall_point_error_km < 0:
        raise ValueError("Landfall point error cannot be negative.")
    return landfall_point_error_km * math.sqrt(math.pi / 2.0)

def hoyt_mean_magnitude(sigma_a_km: float, sigma_c_km: float, quadrature_order: int = HOYT_ANGULAR_QUADRATURE_ORDER) -> float:
    """Deterministic mean magnitude for independent normal components."""
    if sigma_a_km <= 0 or sigma_c_km <= 0:
        raise ValueError("Hoyt component sigmas must be positive.")
    nodes, weights = np.polynomial.legendre.leggauss(quadrature_order)
    theta = (nodes + 1.0) * (math.pi / 4.0)
    w = weights * (math.pi / 4.0)
    angular = np.sqrt(sigma_a_km**2 * np.cos(theta)**2 + sigma_c_km**2 * np.sin(theta)**2)
    angular_mean = (2.0 / math.pi) * float(np.sum(w * angular))
    return math.sqrt(math.pi / 2.0) * angular_mean

def sigma_a_for_dpe(dpe_km: float, sigma_c_km: float, quadrature_order: int = HOYT_ANGULAR_QUADRATURE_ORDER) -> float:
    """Diagnostic numerical solution for sigma_a at a requested DPE."""
    if dpe_km <= 0 or sigma_c_km <= 0:
        raise ValueError("DPE and sigma_c must be positive.")
    lo, hi = 0.0, max(3.0 * dpe_km, 2.0 * sigma_c_km)
    while hoyt_mean_magnitude(hi, sigma_c_km, quadrature_order) < dpe_km:
        hi *= 2.0
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        if hoyt_mean_magnitude(mid, sigma_c_km, quadrature_order) < dpe_km:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)

def landfall_sampler_parameters(lead_time_h: int) -> tuple[float, float]:
    """Return the locked (sigma_a, sigma_c) pair for every landfall scenario lead.

    Anisotropic at 12-72 h (§6.5.4). Isotropic at T-0 (§6.5, analysis uncertainty).
    """
    if lead_time_h not in LANDFALL_SAMPLER_LEADS_H:
        raise ValueError(f"No landfall sampler parameters are defined at {lead_time_h} h.")
    return EXPECTED_LANDFALL_SIGMA_A_KM[lead_time_h], EXPECTED_LANDFALL_SIGMA_C_KM[lead_time_h]


# Backwards-compatible alias.
anisotropic_landfall_parameters = landfall_sampler_parameters

def sample_anisotropic_landfall_errors(sigma_a_km: float, sigma_c_km: float, n_members: int, random_seed: int = 42, latent_draws: tuple[np.ndarray, np.ndarray] | None = None) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Sample persistent standardized along/cross-track errors."""
    if n_members <= 0 or sigma_a_km <= 0 or sigma_c_km <= 0:
        raise ValueError("Invalid sampler parameters.")
    if latent_draws is None:
        rng = np.random.default_rng(random_seed)
        za, zc = rng.normal(size=n_members), rng.normal(size=n_members)
    else:
        za, zc = latent_draws
        if len(za) != n_members or len(zc) != n_members:
            raise ValueError("latent_draws length must equal n_members.")
    along, cross = sigma_a_km * za, sigma_c_km * zc
    return za, zc, np.hypot(along, cross)


# ---------------------------------------------------------------------
# Assertion helpers
# ---------------------------------------------------------------------

def assert_close(
    actual: float,
    expected: float,
    tolerance: float,
    label: str,
) -> None:
    if not math.isclose(
        actual,
        expected,
        rel_tol=0.0,
        abs_tol=tolerance,
    ):
        raise AssertionError(
            f"{label}: expected {expected}, got {actual}"
        )


def assert_equal(
    actual,
    expected,
    label: str,
) -> None:
    if actual != expected:
        raise AssertionError(
            f"{label}: expected {expected}, got {actual}"
        )


# ---------------------------------------------------------------------
# Load CSV
# ---------------------------------------------------------------------

def load_raw_error_table(
    csv_path: Path = ERROR_CSV,
) -> pd.DataFrame:
    if not csv_path.exists():
        raise FileNotFoundError(
            f"Forecast-error CSV not found: {csv_path}"
        )

    df = pd.read_csv(csv_path)

    required_columns = {
        "lead_time_h",
        "track_error_km",
        "track_n",
        "landfall_point_error_km",
        "landfall_n",
        "landfall_time_error_h",
        "intensity_mae_kt",
        "intensity_rmse_kt",
        "intensity_n",
        "track_sigma_km",
        "landfall_sigma_km",
        "source_period",
        "track_provenance",
        "landfall_provenance",
        "intensity_provenance",
    }

    missing = required_columns - set(df.columns)

    if missing:
        raise AssertionError(
            f"CSV is missing required columns: {sorted(missing)}"
        )

    return df


# ---------------------------------------------------------------------
# Basic table assertions
# ---------------------------------------------------------------------

def validate_schema(df: pd.DataFrame) -> None:
    assert df["lead_time_h"].is_unique, (
        "Each lead time must appear exactly once."
    )

    assert df["lead_time_h"].notna().all(), (
        "Lead time cannot contain missing values."
    )

    assert (df["lead_time_h"] >= 0).all(), (
        "Lead time cannot be negative."
    )

    # All primary published rows must use the locked period.
    published_rows = df[
        (
            df["track_provenance"] == "imd_published"
        )
        | (
            df["landfall_provenance"] == "imd_published"
        )
        | (
            df["intensity_provenance"] == "imd_published"
        )
    ]

    assert (
        published_rows["source_period"] == "2020-2024"
    ).all(), (
        "Published primary rows must use source period 2020-2024."
    )


def validate_provenance(df: pd.DataFrame) -> None:
    allowed = {
        "imd_published",
        "assumed",
        "rejected_small_sample",
        "no_data",
        "REQUIRED_EXTRACTION",
    }

    for column in (
        "track_provenance",
        "landfall_provenance",
        "intensity_provenance",
    ):
        unknown = set(df[column].dropna()) - allowed

        assert not unknown, (
            f"Unknown provenance in {column}: {sorted(unknown)}"
        )


# ---------------------------------------------------------------------
# Required ladder assertions
# ---------------------------------------------------------------------

def validate_track_ladder(df: pd.DataFrame) -> None:
    indexed = df.set_index("lead_time_h")

    for lead in TRACK_LADDER_H:
        assert lead in indexed.index, (
            f"Required track lead time {lead} h is missing."
        )

        row = indexed.loc[lead]

        assert pd.notna(row["track_error_km"]), (
            f"Track error missing at {lead} h."
        )

        assert row["track_provenance"] == "imd_published", (
            f"Track lead {lead} h must be imd_published."
        )

    # Locked sample counts.
    for lead, expected_n in EXPECTED_TRACK_N.items():
        actual_n = int(indexed.loc[lead, "track_n"])
        assert_equal(
            actual_n,
            expected_n,
            f"Track n at {lead} h",
        )

    # Locked error values.
    for lead, expected in EXPECTED_TRACK_KM.items():
        actual = float(indexed.loc[lead, "track_error_km"])

        assert_close(
            actual,
            expected,
            VALUE_TOL_KM,
            f"Track error at {lead} h",
        )

    # Strictly increasing approach-path uncertainty.
    values = [
        float(indexed.loc[lead, "track_error_km"])
        for lead in TRACK_LADDER_H
    ]

    assert all(
        later > earlier
        for earlier, later in zip(values, values[1:])
    ), (
        "General track error must strictly increase "
        "across 12 -> 24 -> 48 -> 72 -> 96 -> 120 h."
    )

    # Explicitly protect the long-lead decision.
    assert int(indexed.loc[120, "track_n"]) >= 10, (
        "120 h track error cannot be used with n < 10."
    )


def validate_landfall_ladder(df: pd.DataFrame) -> None:
    indexed = df.set_index("lead_time_h")

    for lead in LANDFALL_LADDER_H:
        assert lead in indexed.index, (
            f"Required landfall lead time {lead} h is missing."
        )

        row = indexed.loc[lead]

        assert pd.notna(row["landfall_point_error_km"]), (
            f"Landfall point error missing at {lead} h."
        )

        assert row["landfall_provenance"] in {
            "imd_published",
            "assumed",
        }, (
            f"Landfall lead {lead} h has unusable provenance: "
            f"{row['landfall_provenance']}"
        )

    # Locked point errors.
    for lead, expected in EXPECTED_LANDFALL_POINT_KM.items():
        actual = float(
            indexed.loc[lead, "landfall_point_error_km"]
        )

        assert_close(
            actual,
            expected,
            VALUE_TOL_KM,
            f"Landfall point error at {lead} h",
        )

    # Locked time errors.
    for lead, expected in EXPECTED_LANDFALL_TIME_H.items():
        actual = float(
            indexed.loc[lead, "landfall_time_error_h"]
        )

        assert_close(
            actual,
            expected,
            VALUE_TOL_H,
            f"Landfall time error at {lead} h",
        )

    # Strictly increasing point uncertainty.
    values = [
        float(indexed.loc[lead, "landfall_point_error_km"])
        for lead in LANDFALL_LADDER_H
    ]

    assert all(
        later > earlier
        for earlier, later in zip(values, values[1:])
    ), (
        "Landfall point error must strictly increase "
        "across 0 -> 12 -> 24 -> 48 -> 72 h."
    )

    # Blueprint distinction:
    # 96/120 h landfall values must NOT be used.
    for lead in (96, 120):
        if lead in indexed.index:
            assert pd.isna(
                indexed.loc[lead, "landfall_point_error_km"]
            ), (
                f"Landfall error at {lead} h must remain unavailable; "
                "96/120 h are approach-path only."
            )


# ---------------------------------------------------------------------
# Intensity assertions
# ---------------------------------------------------------------------

def validate_intensity(df: pd.DataFrame) -> None:
    indexed = df.set_index("lead_time_h")

    for lead in EXPECTED_INTENSITY_MAE_KT:
        assert lead in indexed.index, (
            f"Required intensity lead time {lead} h is missing."
        )

        row = indexed.loc[lead]

        assert row["intensity_provenance"] == "imd_published", (
            f"Intensity lead {lead} h must be imd_published."
        )

        mae = float(row["intensity_mae_kt"])
        rmse = float(row["intensity_rmse_kt"])

        assert_close(
            mae,
            EXPECTED_INTENSITY_MAE_KT[lead],
            INTENSITY_TOL_KT,
            f"Intensity MAE at {lead} h",
        )

        assert_close(
            rmse,
            EXPECTED_INTENSITY_RMSE_KT[lead],
            INTENSITY_TOL_KT,
            f"Intensity RMSE at {lead} h",
        )

        assert rmse >= mae, (
            f"Intensity RMSE must be >= MAE at {lead} h."
        )

    # Explicit pooled 12 h reconstruction.
    assert_close(
        float(indexed.loc[12, "intensity_mae_kt"]),
        3.384171,
        INTENSITY_TOL_KT,
        "12 h intensity MAE",
    )

    assert_close(
        float(indexed.loc[12, "intensity_rmse_kt"]),
        4.960727,
        INTENSITY_TOL_KT,
        "12 h pooled intensity RMSE",
    )


# ---------------------------------------------------------------------
# Rayleigh sigma assertions
# ---------------------------------------------------------------------

def validate_sigma_columns(df: pd.DataFrame) -> None:
    for _, row in df.iterrows():

        # Track sigma
        if pd.notna(row["track_error_km"]):
            expected = mean_error_to_rayleigh_sigma(
                float(row["track_error_km"])
            )

            stored = float(row["track_sigma_km"])

            assert_close(
                stored,
                expected,
                SIGMA_TOL_KM,
                (
                    f"Track sigma at "
                    f"{int(row['lead_time_h'])} h"
                ),
            )

        # Landfall sigma is dimension-aware: LPE is a 1-D coastal displacement,
        # so the sampler scale is sigma_c = 1.25331 * LPE, NOT 0.79788 * LPE.
        # The CSV column landfall_sigma_km holds the legacy 2-D audit value;
        # both are checked so neither can drift unnoticed.
        lead = int(row["lead_time_h"])
        if pd.notna(row["landfall_point_error_km"]):
            lpe = float(row["landfall_point_error_km"])

            assert_close(
                float(row["landfall_sigma_km"]),
                mean_error_to_rayleigh_sigma(lpe),
                SIGMA_TOL_KM,
                f"Legacy landfall audit sigma at {lead} h (2-D basis)",
            )

            if lead in ANISOTROPIC_LANDFALL_LEADS_H:
                assert_close(
                    EXPECTED_LANDFALL_SIGMA_C_KM[lead],
                    landfall_sigma_c_from_lpe(lpe),
                    0.05,
                    f"Landfall sigma_c at {lead} h (1-D basis)",
                )


def validate_anisotropic_landfall_model(df: pd.DataFrame) -> None:
    """Validate the §6.5 landfall model: calibration, sampling, and persistence."""
    indexed = df.set_index("lead_time_h")

    # Assertion 12: exact Hoyt/Rayleigh limiting identity.
    # Cheapest assertion in the suite; it is what catches a wrong polar Jacobian.
    sigma = 50.0
    expected = sigma * math.sqrt(math.pi / 2.0)
    assert_close(hoyt_mean_magnitude(sigma, sigma), expected, 1e-10, "Hoyt equal-sigma Rayleigh limit")

    # Assertion: every landfall scenario lead has usable sampler parameters.
    # T-0 previously fell through this gap with sigma_a = sigma_c = None.
    for lead in LANDFALL_LADDER_H:
        sa, sc = landfall_sampler_parameters(lead)
        assert sa is not None and sc is not None and sa > 0 and sc > 0, (
            f"Landfall scenario lead {lead} h has no usable sampler parameters."
        )

    # T-0 must be isotropic and must match the assumed 5 km position error.
    sa0, sc0 = landfall_sampler_parameters(0)
    assert_close(sa0, sc0, 1e-12, "T-0 landfall sampler must be isotropic")
    assert_close(sa0, mean_error_to_rayleigh_sigma(T0_POSITION_ERROR_KM), 1e-12,
                 "T-0 sigma must equal the 2-D conversion of the assumed 5 km position error")
    assert_close(float(indexed.loc[0, "landfall_point_error_km"]), T0_POSITION_ERROR_KM,
                 VALUE_TOL_KM, "T-0 landfall point error")

    for lead in ANISOTROPIC_LANDFALL_LEADS_H:
        dpe = float(indexed.loc[lead, "track_error_km"])
        lpe = float(indexed.loc[lead, "landfall_point_error_km"])
        sigma_a, sigma_c = landfall_sampler_parameters(lead)

        assert_close(sigma_a, EXPECTED_LANDFALL_SIGMA_A_KM[lead], 0.05, f"§6.5 sigma_a at {lead} h")
        assert_close(sigma_c, EXPECTED_LANDFALL_SIGMA_C_KM[lead], 0.05, f"§6.5 sigma_c at {lead} h")

        # Dimension-aware: LPE is a 1-D coastal displacement (half-normal).
        implied_lpe = sigma_c * math.sqrt(2.0 / math.pi)
        assert_close(implied_lpe, lpe, 0.05, f"§6.5 LPE calibration at {lead} h")

        # DPE is a 2-D magnitude (Hoyt). Tolerance sized to the real residual
        # (worst observed 0.004 km), not to a round number.
        implied_dpe = hoyt_mean_magnitude(sigma_a, sigma_c)
        assert abs(implied_dpe - dpe) <= DPE_CONSISTENCY_TOL_KM, (
            f"§6.5 DPE consistency failed at {lead} h: "
            f"implied={implied_dpe:.4f}, target={dpe:.4f}, tol={DPE_CONSISTENCY_TOL_KM}."
        )

        # Seeded Monte Carlo, per-lead generator so adding a lead cannot
        # silently shift every subsequent result.
        n = 500
        _, zc, mag = sample_anisotropic_landfall_errors(
            sigma_a, sigma_c, n, random_seed=SAMPLER_SEED + lead
        )
        sampled_lpe = float(np.mean(np.abs(sigma_c * zc)))
        sampled_dpe = float(np.mean(mag))
        lpe_se = sigma_c * math.sqrt(1.0 - 2.0 / math.pi) / math.sqrt(n)
        dpe_se = float(np.std(mag, ddof=1)) / math.sqrt(n)
        assert abs(sampled_lpe - lpe) <= 3 * lpe_se, f"§6.5 LPE sampling failed at {lead} h"
        assert abs(sampled_dpe - implied_dpe) <= 3 * dpe_se, f"§6.5 DPE sampling failed at {lead} h"

    # -----------------------------------------------------------------
    # Persistence (§6.5.3): a member draws its wrongness once and keeps it.
    #
    # This MUST route through the sampler. Correlating two scalar multiples
    # of the same array is true by algebra and would pass even if the
    # sampler redrew independently at every call.
    # -----------------------------------------------------------------
    n = 500

    # (a) same seed across leads -> identical latent draws
    za_ref, zc_ref, _ = sample_anisotropic_landfall_errors(
        *landfall_sampler_parameters(12), n, random_seed=SAMPLER_SEED
    )
    for lead in (24, 48, 72):
        za, zc, _ = sample_anisotropic_landfall_errors(
            *landfall_sampler_parameters(lead), n, random_seed=SAMPLER_SEED
        )
        assert float(np.corrcoef(za_ref, za)[0, 1]) > 1 - 1e-12, (
            f"Along-track latent draws are not persistent between 12 h and {lead} h."
        )
        assert float(np.corrcoef(zc_ref, zc)[0, 1]) > 1 - 1e-12, (
            f"Cross-track latent draws are not persistent between 12 h and {lead} h."
        )

    # (b) explicit latent_draws path, which is what Stage 2 uses along a track
    rng = np.random.default_rng(SAMPLER_SEED)
    shared = (rng.normal(size=n), rng.normal(size=n))
    offsets = {}
    for lead in ANISOTROPIC_LANDFALL_LEADS_H:
        sa, sc = landfall_sampler_parameters(lead)
        za, zc, _ = sample_anisotropic_landfall_errors(sa, sc, n, latent_draws=shared)
        offsets[lead] = (sa * za, sc * zc)
    for lead in (24, 48, 72):
        assert float(np.corrcoef(offsets[12][0], offsets[lead][0])[0, 1]) > 1 - 1e-12, (
            f"Along-track offsets not perfectly correlated between 12 h and {lead} h."
        )
        assert float(np.corrcoef(offsets[12][1], offsets[lead][1])[0, 1]) > 1 - 1e-12, (
            f"Cross-track offsets not perfectly correlated between 12 h and {lead} h."
        )

    # (c) along and cross components must be independent within a member
    za, zc, _ = sample_anisotropic_landfall_errors(
        *landfall_sampler_parameters(48), 20000, random_seed=SAMPLER_SEED
    )
    assert abs(float(np.corrcoef(za, zc)[0, 1])) < 0.05, (
        "Along-track and cross-track latent draws must be independent."
    )

    # (d) determinism
    za1, zc1, _ = sample_anisotropic_landfall_errors(1, 1, n, random_seed=SAMPLER_SEED)
    za2, zc2, _ = sample_anisotropic_landfall_errors(1, 1, n, random_seed=SAMPLER_SEED)
    assert np.array_equal(za1, za2) and np.array_equal(zc1, zc2), "Sampler is not deterministic."

    # (e) 96/120 h must never acquire landfall sampler parameters
    assert 96 not in LANDFALL_SAMPLER_LEADS_H and 120 not in LANDFALL_SAMPLER_LEADS_H


# ---------------------------------------------------------------------
# Sampler safety assertion
# ---------------------------------------------------------------------

def validate_sampler_safety(df: pd.DataFrame) -> None:
    """
    Ensure unusable provenance tags cannot silently enter the model.

    Each uncertainty component has its own scenario ladder:

    - Track / approach path: 12, 24, 48, 72, 96, 120 h
    - Landfall point:        0, 12, 24, 48, 72 h
    - Intensity:             12, 24, 48, 72 h

    A `no_data` landfall field at 96/120 h is intentional and must
    remain in the CSV for auditability.
    """

    forbidden = {
        "rejected_small_sample",
        "no_data",
        "REQUIRED_EXTRACTION",
    }

    # -------------------------------------------------------------
    # Track / approach-path component
    # -------------------------------------------------------------

    track_rows = df[
        df["lead_time_h"].isin(TRACK_LADDER_H)
    ]

    bad_track = track_rows["track_provenance"].isin(forbidden)

    assert not bad_track.any(), (
        "Forbidden provenance reached the track sampler: "
        f"{track_rows.loc[bad_track, 'track_provenance'].tolist()}"
    )

    # -------------------------------------------------------------
    # Landfall-point component
    # -------------------------------------------------------------

    landfall_rows = df[
        df["lead_time_h"].isin(LANDFALL_LADDER_H)
    ]

    bad_landfall = (
        landfall_rows["landfall_provenance"].isin(forbidden)
    )

    assert not bad_landfall.any(), (
        "Forbidden provenance reached the landfall sampler: "
        f"{landfall_rows.loc[bad_landfall, 'landfall_provenance'].tolist()}"
    )

    # -------------------------------------------------------------
    # Intensity component
    # -------------------------------------------------------------

    intensity_rows = df[
        df["lead_time_h"].isin(
            EXPECTED_INTENSITY_MAE_KT.keys()
        )
    ]

    bad_intensity = (
        intensity_rows["intensity_provenance"].isin(forbidden)
    )

    assert not bad_intensity.any(), (
        "Forbidden provenance reached the intensity sampler: "
        f"{intensity_rows.loc[bad_intensity, 'intensity_provenance'].tolist()}"
    )

    # -------------------------------------------------------------
    # Explicitly verify that excluded 96/120 h landfall fields
    # remain excluded rather than being accidentally consumed.
    # -------------------------------------------------------------

    excluded_landfall_rows = df[
        df["lead_time_h"].isin((96, 120))
    ]

    for _, row in excluded_landfall_rows.iterrows():

        assert row["landfall_provenance"] == "no_data", (
            f"Expected no_data landfall provenance at "
            f"{int(row['lead_time_h'])} h."
        )

        assert pd.isna(
            row["landfall_point_error_km"]
        ), (
            f"Landfall point error must remain unavailable at "
            f"{int(row['lead_time_h'])} h."
        )

    # If we reach here, no excluded landfall field has been
    # accidentally admitted into the sampler.


# ---------------------------------------------------------------------
# Monte Carlo sampler calibration assertion
# ---------------------------------------------------------------------

def validate_rayleigh_sampler_calibration(
    df: pd.DataFrame,
    n_members: int = 500,
    random_seed: int = 42,
    verbose: bool = True,
) -> None:
    """
    Validate the Rayleigh position-error sampler.

    The published IMD statistic is the mean error magnitude.
    For a Rayleigh distribution:

        E[R] = sigma * sqrt(pi / 2)

    therefore:

        sigma = published_mean * sqrt(2 / pi)

    Two checks are performed:

    1. Analytic calibration:
       sigma must equal the value implied by the published mean.

    2. Seeded Monte Carlo check:
       the sampled mean must fall within 3 standard errors of the
       published mean.

    The Monte Carlo tolerance is derived from the Rayleigh variance,
    rather than being an arbitrary percentage.

    This is intentionally a validation of the sampler construction,
    not a claim that a finite 500-member ensemble must reproduce the
    published mean exactly.
    """

    assert n_members > 0, "n_members must be positive."

    rayleigh_factor = np.sqrt(2.0 / np.pi)

    # Rayleigh standard deviation:
    #
    # SD(R) = sigma * sqrt((4 - pi) / 2)
    rayleigh_sd_factor = np.sqrt((4.0 - np.pi) / 2.0)

    # -------------------------------------------------------------
    # Test the approach-path track ladder
    # -------------------------------------------------------------

    for lead in TRACK_LADDER_H:

        row = df.loc[
            df["lead_time_h"] == lead
        ]

        assert len(row) == 1, (
            f"Expected exactly one CSV row at {lead} h."
        )

        row = row.iloc[0]

        published_mean = float(
            row["track_error_km"]
        )

        sigma_expected = (
            published_mean * rayleigh_factor
        )

        # ---------------------------------------------------------
        # Analytic calibration
        # ---------------------------------------------------------

        sigma_from_csv = float(
            row["track_sigma_km"]
        )

        assert np.isclose(
            sigma_from_csv,
            sigma_expected,
            rtol=1e-5,
            atol=1e-6,
        ), (
            f"Track sigma calibration failed at {lead} h: "
            f"CSV={sigma_from_csv:.6f}, "
            f"expected={sigma_expected:.6f}."
        )

        # ---------------------------------------------------------
        # Construct the sampler directly from the calibrated sigma
        # ---------------------------------------------------------

        # Per-lead generator: adding or removing a lead time must not
        # silently change the draws at every other lead.
        rng = np.random.default_rng(random_seed + lead)

        dx = rng.normal(loc=0.0, scale=sigma_expected, size=n_members)
        dy = rng.normal(loc=0.0, scale=sigma_expected, size=n_members)

        sampled_magnitude = np.sqrt(
            dx**2 + dy**2
        )

        sampled_mean = float(
            np.mean(sampled_magnitude)
        )

        # ---------------------------------------------------------
        # Monte Carlo standard error
        # ---------------------------------------------------------

        population_sd = (
            sigma_expected * rayleigh_sd_factor
        )

        standard_error = (
            population_sd / np.sqrt(n_members)
        )

        tolerance = 3.0 * standard_error

        assert abs(
            sampled_mean - published_mean
        ) < tolerance, (
            f"Rayleigh Monte Carlo calibration failed at "
            f"{lead} h: "
            f"sample_mean={sampled_mean:.3f} km, "
            f"published_mean={published_mean:.3f} km, "
            f"3SE={tolerance:.3f} km, "
            f"N={n_members}, seed={random_seed}."
        )

    # The landfall ladder is NOT validated here. §6.5 supersedes the isotropic
    # Rayleigh treatment of landfall error: LPE is a 1-D coastal displacement
    # taking sigma_c = 1.25331 * LPE, not the 2-D factor 0.79788 * LPE used
    # above for DPE. Those differ by pi/2 (36%). Landfall calibration and
    # sampling are validated in validate_anisotropic_landfall_model().

    if verbose:
        print(
            f"Track (Rayleigh) sampler calibration PASSED "
            f"(N={n_members}, seed={random_seed}, "
            f"tolerance=3 SE)."
        )

# ---------------------------------------------------------------------
# Main validation entry point
# ---------------------------------------------------------------------

def validate_error_model(
    csv_path: Path = ERROR_CSV,
) -> pd.DataFrame:

    df = load_raw_error_table(csv_path)

    validate_schema(df)
    validate_provenance(df)
    validate_track_ladder(df)
    validate_landfall_ladder(df)
    validate_intensity(df)
    validate_sigma_columns(df)
    validate_sampler_safety(df)
    validate_rayleigh_sampler_calibration(df)
    validate_anisotropic_landfall_model(df)

    return df


# ---------------------------------------------------------------------
# Runtime model loader
# ---------------------------------------------------------------------

def load_error_models(
    csv_path: Path = ERROR_CSV,
) -> dict[int, ForecastError]:

    df = load_raw_error_table(csv_path)

    validate_schema(df)
    validate_provenance(df)
    validate_track_ladder(df)
    validate_landfall_ladder(df)
    validate_intensity(df)
    validate_sigma_columns(df)
    validate_sampler_safety(df)
    validate_rayleigh_sampler_calibration(
        df,
        verbose=False,
    )
    validate_anisotropic_landfall_model(df)
    indexed = df.set_index("lead_time_h")

    models: dict[int, ForecastError] = {}

    for lead in sorted(
        set(TRACK_LADDER_H) | set(LANDFALL_LADDER_H)
    ):

        row = indexed.loc[lead]

        track_error = float(row["track_error_km"])
        track_sigma = mean_error_to_rayleigh_sigma(
            track_error
        )

        landfall_error = None
        landfall_sigma = None
        landfall_sigma_a = None
        landfall_sigma_c = None

        if pd.notna(row["landfall_point_error_km"]):
            landfall_error = float(
                row["landfall_point_error_km"]
            )

            # Legacy Rayleigh value retained for audit compatibility.
            landfall_sigma = mean_error_to_rayleigh_sigma(
                landfall_error
            )

            # §6.5: anisotropic at 12-72 h, isotropic at T-0.
            if lead in LANDFALL_SAMPLER_LEADS_H:
                landfall_sigma_a, landfall_sigma_c = landfall_sampler_parameters(lead)

        models[lead] = ForecastError(
            lead_time_h=lead,

            track_error_km=track_error,
            track_sigma_km=track_sigma,
            track_n=(
                int(row["track_n"])
                if pd.notna(row["track_n"])
                else None
            ),
            track_provenance=str(
                row["track_provenance"]
            ),

            landfall_point_error_km=landfall_error,
            landfall_sigma_km=landfall_sigma,
            landfall_sigma_a_km=landfall_sigma_a,
            landfall_sigma_c_km=landfall_sigma_c,
            landfall_n=(
                int(row["landfall_n"])
                if pd.notna(row["landfall_n"])
                else None
            ),
            landfall_provenance=str(
                row["landfall_provenance"]
            ),

            landfall_time_error_h=(
                float(row["landfall_time_error_h"])
                if pd.notna(row["landfall_time_error_h"])
                else None
            ),

            intensity_mae_kt=(
                float(row["intensity_mae_kt"])
                if pd.notna(row["intensity_mae_kt"])
                else None
            ),
            intensity_rmse_kt=(
                float(row["intensity_rmse_kt"])
                if pd.notna(row["intensity_rmse_kt"])
                else None
            ),
            intensity_n=(
                int(row["intensity_n"])
                if pd.notna(row["intensity_n"])
                else None
            ),
            intensity_provenance=str(
                row["intensity_provenance"]
            ),
        )

    return models


# ---------------------------------------------------------------------
# Standalone execution
# ---------------------------------------------------------------------

if __name__ == "__main__":

    df = validate_error_model()

    models = load_error_models()

    print()
    print("P3 FORECAST ERROR MODEL ASSERTIONS PASSED")
    print("§6.5 anisotropic landfall model PASSED")
    print()
    print("Approach-path track ladder:")
    print("Lead (h) | DPE (km) | n | sigma (km)")
    print("--------------------------------------")

    for lead in TRACK_LADDER_H:
        m = models[lead]

        print(
            f"{lead:8d} | "
            f"{m.track_error_km:8.3f} | "
            f"{m.track_n:3d} | "
            f"{m.track_sigma_km:10.3f}"
        )

    print()
    print("Landfall-point ladder (§6.5):")
    print("Lead (h) | LPE (km) | n | sigma_a | sigma_c")
    print("-----------------------------------------------")

    for lead in LANDFALL_LADDER_H:
        m = models[lead]

        sigma_a = (
            f"{m.landfall_sigma_a_km:7.3f}"
            if m.landfall_sigma_a_km is not None
            else "   None"
        )
        sigma_c = (
            f"{m.landfall_sigma_c_km:7.3f}"
            if m.landfall_sigma_c_km is not None
            else "   None"
        )

        print(
            f"{lead:8d} | "
            f"{m.landfall_point_error_km:8.3f} | "
            f"{str(m.landfall_n):>3} | "
            f"{sigma_a} | {sigma_c}"
        )

    print()
    print("12 h intensity:")
    m = models[12]

    print(
        f"MAE  = {m.intensity_mae_kt:.6f} kt"
    )
    print(
        f"RMSE = {m.intensity_rmse_kt:.6f} kt"
    )

    print()