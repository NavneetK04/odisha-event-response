# Forecast Error Model

**Project 3 — Cyclone Event Response: Forecast Uncertainty and the Decision Clock**

This document derives every forecast-error number the project uses, states where
each one came from, and names the two that are assumptions rather than published
statistics.

It exists because the entire project rests on this model. If the error growth
curve is wrong, every loss fan chart, every attachment probability and the whole
variance decomposition inherits the error, and none of it announces itself. The
arithmetic is therefore written out rather than left in code.

Companion file: `data/imd_forecast_errors.csv`.

---

## 1. Source

India Meteorological Department, RSMC New Delhi, official verification of
tropical cyclone forecasts over the North Indian Ocean.

| What | Where |
|---|---|
| Landfall point error, Table 1, 2003–2025 | [landfall-forecast.php](https://rsmcnewdelhi.imd.gov.in/landfall-forecast.php) |
| Landfall time error, Table 2, 2003–2025 | [landfall-forecast.php](https://rsmcnewdelhi.imd.gov.in/landfall-forecast.php) |
| Track forecast error (DPE), Tables 1–2 | [track-forecast.php](https://rsmcnewdelhi.imd.gov.in/track-forecast.php) |
| Intensity forecast error | [intensity-forecast.php](https://rsmcnewdelhi.imd.gov.in/intensity-forecast.php) |
| Annual Verification Report 2024 | [Annual_Veri_2024.pdf](https://rsmcnewdelhi.imd.gov.in/uploads/Annual_Veri_2024.pdf) |
| Method review, lead-time structure | [TC Forecast Verification by IMD for NIO: a Review](https://www.sciencedirect.com/science/article/pii/S2225603218301140) |

IMD publishes these as **annual means with the number of verified cases in
parentheses**, by lead time, in 12-hour steps from 12 h to 120 h. Track and
landfall forecasts have been issued to 72 h since 2009 and to 120 h since 2013.

The case counts matter and are carried through this document. A cell reading
`69.5 (10)` is a mean over ten verified landfalls, not a population statistic.

---

## 2. Method, and the check that validates it

For mean-error quantities (landfall point/time, track DPE and absolute intensity
error), multi-year figures are **sample-weighted means** across the annual values:

```
                Σ_years ( annual_mean × n_cases )
period_mean  =  ─────────────────────────────────
                       Σ_years n_cases
```

**This method is not assumed — it is verified.** Applying it to IMD's annual
table for 2020–2024 reproduces the multi-year figures IMD itself publishes:

| Quantity | Lead | Reconstructed | IMD published | Match |
|---|---:|---:|---:|:--:|
| Landfall point error | 24 h | 16.19 km | 16.2 km | ✓ |
| Landfall point error | 48 h | 39.32 km | 39.3 km | ✓ |
| Landfall point error | 72 h | 69.52 km | 69.5 km | ✓ |
| Landfall time error | 24 h | 2.887 h | 2.9 h | ✓ |
| Landfall time error | 48 h | 4.189 h | 4.2 h | ✓ |
| Landfall time error | 72 h | 7.450 h | 7.5 h | ✓ |

All six reconciled to the published precision. Every other number in this
document is produced by the same arithmetic on the same table, so it carries the
same status as the published values — it is a recomputation, not an estimate.

A simple unweighted mean across years does **not** reproduce the landfall/track
figures. RMSE is handled separately because it is a square-root second-moment
quantity; §3.3 gives the correct pooled-RMSE reconstruction. If a future edit
changes either aggregation rule, the reconciliation checks will fail, which is
the point of recording them.

---

## 3. Primary error model — 2020 to 2024

### 3.1 Landfall point and time error

| Lead | Point error | n | Time error | n |
|---:|---:|---:|---:|---:|
| 0 h | **5.0 km** (assumed, §5) | — | 0.0 h (by construction) | — |
| 12 h | **10.53 km** | 24 | **1.85 h** | 24 |
| 24 h | **16.19 km** | 23 | **2.89 h** | 23 |
| 36 h | 29.36 km | 21 | 3.76 h | 21 |
| 48 h | **39.32 km** | 19 | **4.19 h** | 19 |
| 60 h | 49.24 km | 16 | 6.42 h | 16 |
| 72 h | **69.52 km** | 10 | **7.45 h** | 10 |

Bold rows are the project's lead-time ladder. 36 h and 60 h are recorded because
they are published and because they are the best available evidence that the
growth curve is smooth; they are not run as scenarios.

The landfall **time** error is not used by the sampler. It is used in §6.5.6 as a
held-out check on the along-track widths derived from position statistics alone.

### 3.2 General track error (direct position error)

Used for the approach path, **not** for the landfall point. See §7. The 2020–2024
12 h, 96 h and 120 h values below are sample-weighted recomputations from IMD's
annual track table.

| Lead | Mean DPE | n | Status |
|---:|---:|---:|---|
| 12 h | **44.87 km** | 398 | `imd_published` |
| 24 h | **71.97 km** | 332 | `imd_published` |
| 48 h | **110.94 km** | 236 | `imd_published` |
| 72 h | **154.03 km** | 135 | `imd_published` |
| 96 h | **182.69 km** | 55 | `imd_published` |
| 120 h | **245.23 km** | 28 | `imd_published` |

For comparison with the rounded multi-year figures reported by IMD, these are
72.0, 111.0 and 154.0 km at 24/48/72 h. The unrounded values are retained here
so the arithmetic is auditable.

The 96 h and 120 h track values are **not** used as landfall errors. They remain
valid approach-path inputs because both have n ≥ 10 and the 12–120 h weighted
sequence is monotonic: 44.87 → 71.97 → 110.94 → 154.03 → 182.69 → 245.23 km.

### 3.3 Intensity error

| Lead | MAE | RMSE | n |
|---:|---:|---:|---:|
| 12 h | **3.384 kt** | **4.961 kt ≈ 5.0 kt** | 398 |
| 24 h | 5.9 kt | 7.9 kt | not recorded |
| 48 h | 8.3 kt | 11.0 kt | not recorded |
| 72 h | 9.8 kt | 19.2 kt | not recorded |

The 12 h values are reconstructed from IMD's **Table 1a** (absolute intensity
error) and **Table 1b** (root-mean-square intensity error), using the 2020–2024
annual case counts shown in those tables. MAE is aggregated by the case-weighted
mean. RMSE is **not** averaged linearly: because RMSE is a square-root second-moment
quantity, the pooled RMSE is reconstructed as the square root of the
case-weighted mean of the annual squared RMSE values. The result, 4.96073 kt,
rounds to IMD's published 2020–2024 value of 5.0 kt. The published annual inputs are:

| Year | 12 h MAE (kt) | n | 12 h RMSE (kt) | n |
|---:|---:|---:|---:|---:|
| 2020 | 5.0 | 62 | 6.8 | 62 |
| 2021 | 3.5 | 68 | 4.9 | 68 |
| 2022 | 2.4 | 75 | 3.3 | 75 |
| 2023 | 3.7 | 125 | 5.4 | 125 |
| 2024 | 2.3 | 68 | 3.5 | 68 |
| **Total** | | **398** | | **398** |

Thus:

```
12 h MAE  = [5.0×62 + 3.5×68 + 2.4×75 + 3.7×125 + 2.3×68] / 398
          = 3.38417 kt

12 h RMSE = sqrt([6.8²×62 + 4.9²×68 + 3.3²×75 + 5.4²×125 + 3.5²×68] / 398)
          = 4.96073 kt ≈ 5.0 kt (IMD published 2020–2024 value)
```

**The 24, 48 and 72 h intensity figures are taken as published, not
recomputed.** The annual case counts behind them were not transcribed, so the
§2 weighting method cannot be verified for those three rows the way it was for
12 h. They are tagged `imd_published` because they are IMD's own multi-year
values, but they do not carry the reconciliation guarantee the 12 h row does.

Sampler width uses **RMSE, not MAE**, because the sampler needs a second-moment
quantity and MAE understates the tail. The distinction is worth roughly a factor
of two at 72 h, where RMSE is 19.2 kt against an MAE of 9.8 kt — a gap that says
the 72 h intensity error distribution is strongly skewed by a few large misses.

---

## 4. Why the landfall ladder stops at 72 hours

The blueprint specifies lead times `{120, 96, 72, 48, 24, 12, 0}` for the forecast
origin and approach-path reconstruction. **Only the landfall-point error ladder
stops at 72 h.** This distinction is deliberate: the IMD landfall table becomes
sparsely sampled beyond 72 h, while the general track table still has adequate
2020–2024 track cases at 96 h and 120 h.

### 4.1 There is no data at 120 h

Verified landfall forecasts at 120 h, 2020–2024: **zero**. Across the entire
2003–2025 record there is exactly one, in 2015.

No distribution can be fitted to that. Interpolating a value would place a
fabricated number at the most visually prominent end of every fan chart in the
project.

### 4.2 The long-lead landfall cells that do exist are not measurements

Extending the 2020–2024 **landfall** column past 72 h:

| Lead | Point error | n |
|---:|---:|---:|
| 72 h | 69.5 km | 10 |
| 84 h | 49.9 km | 7 |
| 96 h | 55.4 km | 6 |
| 108 h | 28.1 km | 2 |
| 120 h | — | 0 |

The landfall error **falls** as lead time grows. That is not a credible forecast-skill
curve, and the samples are too small to support those cells.

Two things are happening, and both disqualify the landfall cells:

**Sample size.** n drops from 10 to 2 across four rows. At n=2 a single easy
storm sets the mean.

**Selection.** IMD issues a landfall forecast only once a landfall point can be
identified at all. The handful of storms that received a 96 h landfall forecast
were those already on an obvious track four days out. The long-lead sample is
therefore biased toward the easiest cases, which is exactly why the mean falls.

The 2015–2019 window makes the instability unmistakable:

| Lead | 2020–2024 | n | 2015–2019 | n | Ratio |
|---:|---:|---:|---:|---:|---:|
| 72 h | 69.5 km | 10 | 109.3 km | 7 | 1.6× |
| 84 h | 49.9 km | 7 | 268.9 km | 2 | 5.4× |
| 96 h | 55.4 km | 6 | 279.4 km | 2 | 5.0× |
| 108 h | 28.1 km | 2 | 403.6 km | 2 | 14.4× |
| 120 h | — | 0 | 435.9 km | 1 | — |

At 12–72 h the two windows differ by a factor of 1.6–2.8, which is a credible
skill improvement over five years. At 84–108 h they differ by 5–14×, on samples
of one or two. That is noise, not skill.

**These cells are recorded here and deliberately kept out of the data file.**
`validate_landfall_ladder` asserts that `landfall_point_error_km` is null at
96 h and 120 h, and the run fails otherwise:

### 4.3 What this becomes in the write-up

The project therefore keeps **two lead-time structures** rather than forcing one
error statistic to do two different jobs:

- **Approach-path track error:** `{120, 96, 72, 48, 24, 12, 0}`.
- **Landfall-point error:** `{72, 48, 24, 12, 0}`.

At 120/96 h the hindcast can therefore represent **where the cyclone is expected
to travel**, but it does not claim a separately verified 120/96 h landfall-point
uncertainty, and no loss is computed at those leads.

The loss phase therefore runs five lead times rather than seven:
`500 × 5 × 4 = 10,000` wind fields in principle, and **9,500 in practice**,
because Titli's T-72h forecast origin falls 1.43 h before the start of its
best-track record and that panel cannot exist. The unavailable panel is recorded
in `outputs/track_availability.csv` rather than silently omitted.

---

## 5. The two assumptions

Everything above is published. These two are not, and they are labelled
`assumed` in the CSV.

### 5.1 T-0h position error: 5 km

**There is no published value and there cannot be one.** T-0 is analysis time,
not a forecast. What exists is best-track position uncertainty, which IMD does
not publish as a verification statistic.

**Chosen value: 5.0 km.** Reasoning:

- It must be below the 12 h error of 10.53 km, which is the hard published ceiling.
- The best individual years run 5.4 km (2024) and 6.8 km (2021) at 12 h, so
  sub-10 km centre fixing is demonstrably achievable in current practice.
- At landfall the system is inside coastal Doppler radar coverage, which is the
  best-constrained the centre position ever gets.
- 5 km is roughly half the 12 h value, consistent with the observed growth ratio
  between adjacent lead times.

It is applied **isotropically**, as a two-dimensional position uncertainty, which
is why T-0 is the one lead where σ_a = σ_c. See §6.5.7.

**Sensitivity: not run.** The sensitivity this calls for is 0, 5 and 10.53 km,
reporting whether the conclusions hold across that range. It was not performed
and is recorded as future work in the README.

**Do not set this to zero.** The project's hypothesis is that at short lead times
the dominant uncertainty stops being meteorological and becomes exposure data
quality. Setting meteorological error to zero at T-0 guarantees that result by
construction at the endpoint, making the final panel of the variance
decomposition circular. It is the same self-selection failure the blueprint
warns about in §8.2, wearing a different costume.

For the same reason the T-0 row of `imd_forecast_errors.csv` carries 5.0 km in
**both** the track and landfall columns. An earlier version of that file carried
zero in the track column, contradicting this paragraph.

### 5.2 T-0h intensity error: 5 kt

Same logic. Best-track intensity over the NIO rests largely on Dvorak-technique
estimates, which carry meaningful uncertainty even at analysis time. 5 kt is a
judgement, below the 24 h MAE of 5.9 kt. It is recorded in the RMSE column,
because RMSE is the quantity the sampler consumes (§3.3).

**Sensitivity: not run.** The sensitivity this calls for is 0 and 5 kt. Recorded
as future work.

### 5.3 T-0h landfall time error: 0 h

This one *is* zero by construction. At T-0 the storm is making landfall; there is
no timing left to forecast. No sensitivity needed.

---

## 6. From a published mean to a sampler

IMD publishes **mean error magnitudes**. The ensemble needs distributions. The
conversions are written out here because a reviewer will ask for them, and
because there are two of them rather than one.

### 6.1 Two conversions, not one

This is the easiest arithmetic in the project to get wrong, and the result of
getting it wrong is a table that looks entirely reasonable.

**Direct position error (DPE) is a two-dimensional magnitude.** It is the length
of the error vector between forecast and actual position. For independent
zero-mean Gaussian components of equal width σ, that magnitude is **Rayleigh**
distributed:

```
mean   = σ · √(π/2) = 1.25331 · σ      ⟹   σ = 0.79788 · DPE
RMS    = σ · √2     = 1.41421 · σ
median = σ · √(2 ln 2) = 1.17741 · σ
```

**Landfall point error (LPE) is a one-dimensional displacement.** It measures how
far along the coast the forecast landfall point fell from the actual one. A
one-dimensional zero-mean Gaussian has a **half-normal** magnitude:

```
mean = σ · √(2/π) = 0.79788 · σ        ⟹   σ_c = 1.25331 · LPE
```

**The two factors are reciprocal.** Applying the Rayleigh conversion to LPE
understates the cross-track width by 36%: at 12 h it gives 8.40 km where the
correct value is 13.20 km. Applying the half-normal conversion to DPE overstates
by 57% in the other direction. Both errors produce a monotone, plausible-looking
ladder, and neither announces itself.

The one exception is T-0. The assumed 5 km there is a two-dimensional position
uncertainty, not a coastal displacement, so the Rayleigh conversion applies to it
in both the track and the landfall columns. See §6.5.7.

### 6.2 Worked values

Landfall point, 2020–2024. `imd_forecast_errors.csv` carries **three** sigma
columns for the landfall quantity, and only two of them reach the sampler:

| Column | Basis | Formula | Used by |
|---|---|---|---|
| `landfall_sigma_km` | 2-D Rayleigh, as if LPE were a position magnitude | 0.79788 × LPE | nothing; audit only |
| `landfall_sigma_c_km` | 1-D half-normal, the correct basis for a coastal displacement | 1.25331 × LPE | the sampler |
| `landfall_sigma_a_km` | solved jointly against DPE | §6.5.4 | the sampler |

| Lead | LPE | `landfall_sigma_km` *(audit)* | σ_c *(sampler)* | σ_a *(sampler)* |
|---:|---:|---:|---:|---:|
| 0 h | 5.00 km *(assumed)* | 3.99 km | 3.99 km | 3.99 km |
| 12 h | 10.53 km | 8.40 km | **13.20 km** | 52.40 km |
| 24 h | 16.19 km | 12.92 km | **20.29 km** | 84.47 km |
| 48 h | 39.32 km | 31.37 km | **49.28 km** | 120.25 km |
| 72 h | 69.52 km | 55.47 km | **87.13 km** | 154.00 km |

**The audit column is retained deliberately.** It is the value the first version
of this model used, before the dimension of LPE was examined, and it is 36% below
the correct σ_c because it applies a 2-D conversion to a 1-D quantity. Keeping it
and asserting it on every run means the wrong value cannot quietly reappear: if
someone recomputes the landfall sigma with the Rayleigh factor again, the
assertion on the audit column still passes, the assertion on σ_c fails, and the
run stops. Deleting the column would remove that tripwire.

At T-0 all three agree, because the assumed 5 km there is a two-dimensional
position uncertainty rather than a coastal displacement, so the Rayleigh
conversion is the correct one and the error is isotropic. See §6.5.7.

The `landfall_sigma_a_km` and `landfall_sigma_c_km` columns are documentation.
The loader computes both from `landfall_sampler_parameters()` rather than reading
them, so the CSV cannot silently change what the sampler does; the columns exist
so a reader can see the ladder without running the code.

General track error, approach path:

| Lead | Published mean | σ per component |
|---:|---:|---:|
| 12 h | 44.87 km | 35.80 km |
| 24 h | 71.97 km | 57.42 km |
| 48 h | 110.94 km | 88.52 km |
| 72 h | 154.03 km | 122.90 km |
| 96 h | 182.69 km | 145.77 km |
| 120 h | 245.23 km | 195.66 km |

### 6.3 Sampling procedure

For the approach path, where only DPE is available and no decomposition is
possible:

```
σ      = 0.79788 × DPE(T)
dx     ~ N(0, σ)        # along-track component
dy     ~ N(0, σ)        # cross-track component
offset = (dx, dy) rotated into the track-relative frame
```

Equivalently, draw magnitude `r ~ Rayleigh(σ)` and bearing `θ ~ Uniform(0, 2π)`.
The two formulations are identical under equal component widths.

For the landfall-phase ensemble, which is what every loss figure in the project
is built from, this is replaced by the anisotropic scheme in §6.5.

### 6.4 One caveat that survives

**The model assumes the forecast is unbiased.** A Rayleigh or half-normal
magnitude arises from zero-mean components. IMD's published errors include any
systematic bias, so if NIO landfall forecasts carry a consistent directional
bias, this model spreads that bias symmetrically instead of reproducing it.
Stated, not corrected — the published tables give no bias decomposition.

**Isotropy, and why it was abandoned.** The first version of this model drew
position error isotropically, with a single σ per lead, and recorded the
resulting understatement of loss-relevant spread as a known weakness pending an
external decomposition. That weakness is now closed. IMD's own DPE and LPE at the
same lead turn out to constrain the along-track and cross-track components
jointly, so the split can be derived from published data rather than assumed or
imported. §6.5 replaces the isotropic model for the landfall phase. The approach
path in §6.3 remains isotropic, because only DPE is published beyond 72 h and
there is nothing to decompose it against.

### 6.5 The anisotropic landfall model

#### 6.5.1 The scheme, and its two failure modes

Each ensemble member draws three standard normal deviates **once** and carries
them across the whole forecast:

```
z_a, z_c, z_v  ~  N(0,1), drawn once per member

for each forecast horizon h in [0, T]:
    along(h) = sigma_a(h) * z_a
    cross(h) = sigma_c(h) * z_c
    dv(h)    = sigma_v(h) * z_v
```

Correlation between a member's offsets at any two horizons is exactly 1 by
construction, and the displacement grows as σ grows. A member is consistently
left or right of forecast and consistently early or late, which is what forecast
error looks like.

Two things can go wrong, and both produce output that looks correct.

**Independent draws at each horizon.** Reproduces the marginal spread at every
lead and nothing about the structure. A track 80 km left at 48 h and 80 km right
at 24 h is not a forecast error, it is noise, and it would wander across the
forecast corridor crossing the best track repeatedly.

**σ frozen at its landfall value.** Makes every member a rigid translation of the
best track, so the ensemble is a bundle of parallel tracks with no fan at all.
This one is subtle because it satisfies a naive monotonicity check: σ
non-decreasing along the forecast is true of a constant. The assertion therefore
pins both ends, σ at the origin to the T-0 value and σ at landfall to the ladder,
rather than only checking the direction of travel.

#### 6.5.2 Why one σ is not enough

A single σ per lead treats along-track and cross-track error as equal. They are
not equivalent for loss:

- **Cross-track error moves the landfall point along the coast.** It decides
  which exposure the footprint intersects.
- **Along-track error is timing error.** It mostly decides when the storm
  arrives, and moves the landfall point only along its own approach bearing.

Treating them as equal understates the quantity that drives loss and overstates
the one that does not.

#### 6.5.3 The track-relative frame

The offsets are applied in a frame defined by the direction of motion, with
cross-track positive 90 degrees to the right. Three modes exist in the code:

| Mode | Frame heading | Phailin T-72h rejections | Max implied speed |
|---|---|---|---|
| `local` | instantaneous heading at each point | **144 of 500** | 407 km/h |
| `smoothed` | centred circular mean over 6 h | intermediate | intermediate |
| `landfall` | the heading on approach to landfall, applied to every point | **0 of 500** | 24 km/h |

**`landfall` is the default.** The failure in `local` is geometric rather than
statistical. Where the track curves, the frame rotates between successive points,
so a member carrying a large constant along-track deviate is swung about a lever
arm: two adjacent displaced positions end up far apart even though the underlying
best-track positions are one hour of travel apart. Phailin's T-72h approach has
the sharpest curvature in the event set, and the physical plausibility screen
correctly rejected 28.8% of its members at implied translation speeds up to
407 km/h. Anchoring the frame on the landfall heading took that to zero.

This was diagnosed by running it, not by arguing it.

Where a post-landfall extension is present (§6.5.9) the frame is anchored on the
**landfall index**, not the last point of the array. The last point is inland on
a decaying and often turning track, and using its heading would rotate the frame
for every member in every panel.

#### 6.5.4 Deriving the split from published data

IMD publishes two quantities at the same lead time that measure different
projections of the same error:

- **LPE** is a displacement along the coast. It is one-dimensional and
  half-normal, so it constrains σ_c alone (§6.1).
- **DPE** is the magnitude of the full two-dimensional error vector. It
  constrains σ_a and σ_c jointly.

Two equations, two unknowns. With independent zero-mean Gaussian components of
unequal width, the magnitude follows the Hoyt (Nakagami-q) distribution, whose
mean has no elementary closed form:

```
sigma_c                =  LPE * sqrt(pi/2)                                  (1)

E[ sqrt(A^2 + C^2) ]   =  DPE,   A ~ N(0, sigma_a^2),  C ~ N(0, sigma_c^2)  (2)
```

(1) is direct. (2) is solved numerically for σ_a by Gauss-Legendre quadrature
over the angular integral, with σ_c fixed by (1).

**The equal-σ limit is asserted, not assumed.** Forcing σ_a = σ_c must reduce (2)
to the isotropic Rayleigh result, mean = σ·√(π/2). This single check is what
caught a dropped Jacobian factor in the angular integral: the integrand carries
one power of r from the magnitude and one from the polar area element, and losing
either leaves a solver that converges, satisfies both constraints to tolerance,
and produces a table inflated by roughly 25%.

#### 6.5.5 The derived ladder

| Lead | DPE (km) | n | LPE (km) | n | σ_a (km) | σ_c (km) | σ_a/σ_c |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 0 h | 5.00 *(assumed)* | — | 5.00 *(assumed)* | — | 3.99 | 3.99 | 1.00 |
| 12 h | 44.87 | 398 | 10.53 | 24 | 52.40 | 13.20 | 3.97 |
| 24 h | 71.97 | 332 | 16.19 | 23 | 84.47 | 20.29 | 4.16 |
| 48 h | 110.94 | 236 | 39.32 | 19 | 120.25 | 49.28 | 2.44 |
| 72 h | 154.03 | 135 | 69.52 | 10 | 154.00 | 87.13 | 1.77 |

The anisotropy is strongly lead-dependent and peaks at 24 h. **A 24-hour forecast
knows where a storm will come ashore about four times better than it knows when.**
The ratio falls towards 72 h, where the forecast is losing its grip on placement
as well as timing.

#### 6.5.6 An independent check the derivation never used

σ_a is a distance, and dividing it by the forecast's timing uncertainty should
give a translation speed. IMD publishes landfall **time** error separately, and
it was not used anywhere in §6.5.4:

| Lead | σ_a (km) | LTE (h) | σ_t = LTE·√(π/2) | implied speed |
|---:|---:|---:|---:|---:|
| 12 h | 52.40 | 1.85 | 2.32 h | 22.6 km/h |
| 24 h | 84.47 | 2.89 | 3.62 h | 23.3 km/h |
| 48 h | 120.25 | 4.19 | 5.25 h | 22.9 km/h |
| 72 h | 154.00 | 7.45 | 9.34 h | 16.5 km/h |

16 to 23 km/h, which is the right range for North Indian Ocean tropical cyclone
translation, and roughly constant across three of the four leads. Nothing forced
that. The along-track widths were solved from position statistics alone, and the
timing statistics they reproduce were held out.

This is also the diagnostic that first exposed the Jacobian error in §6.5.4. The
inflated σ_a implied about 30 km/h, which is fast for these storms. That was
initially explained away rather than treated as a signal, which is the mistake
worth recording: a physically implausible by-product of a solver that reports
convergence is still evidence against the solver.

#### 6.5.7 T-0 is isotropic, and that is not an oversight

Every positive lead has σ_a > σ_c. T-0 has σ_a = σ_c = 3.99 km.

**Anisotropy is a property of forecast error.** The along and cross split exists
because a forecast constrains a storm's position across its track better than
along it. T-0 is analysis time: the uncertainty is about where the centre is
right now, and it has no preferred direction. Carrying the 12 h anisotropy ratio
down to T-0 would apply a forecast property to a quantity that is not a forecast.

This has a consequence for the assertions. At 12 to 72 h the published landfall
figure is a one-dimensional coastal displacement and is checked against the
cross-track component alone. At T-0 it is a two-dimensional position uncertainty
and must be checked against the magnitude. Comparing the T-0 figure against the
cross component would expect 5.0 km from a quantity whose mean is
5.0 × √(2/π) = 3.19 km, and would fail a sampler behaving exactly as specified.

#### 6.5.8 Interpolation between calibrated leads

The ladder is calibrated at 0, 12, 24, 48 and 72 h. Every forecast horizon in
between takes a **linear** interpolation of σ_a and σ_c, and extrapolation beyond
the ladder is refused rather than silently flat-extended.

Linear is a stated assumption, not a derived one. It matters most across the 0 to
12 h gap, where σ_a rises from 3.99 km to 52.40 km, a factor of thirteen, with no
published data anywhere in between. A power law or an exponential would give
materially different widths at 4 and 8 hours. There is nothing to choose between
them on the evidence, so the simplest is used and recorded here.

#### 6.5.9 Post-landfall track extension

**Status: `derived`. Corrects a defect present in the first full pipeline run.**

§6.5.1 places each member at horizons 0 to T, where T is landfall. Implemented
literally, the ensemble track array ends at landfall. A member displaced backward
along track therefore has every one of its positions offshore, the wind field
never comes ashore, and the loss is exactly zero.

In a loss file that zero is indistinguishable from the zero produced by a member
that made landfall too far from the portfolio to cause damage. The two have
opposite meaning. One is a statement about forecast uncertainty; the other is a
statement about the length of an array.

They are distinguishable in the latent draws. Cross-track displacement does not
depend on z_a at all, so a genuine miss is symmetric in z_a. Truncation affects
only backward displacements, so it is one-sided. Reconstructing the draws from
the seed and testing the sign balance of z_a among zero-loss members gave:

| Event | Leads tested | Fraction of zeros with z_a < 0 | Deviation from symmetric |
|---|---|---|---|
| Fani | 72, 48, 24, 12 | 0.948, **1.000**, **1.000**, **1.000** | 11.1, 11.7, 7.1, 5.2 σ |
| Phailin | 72, 48, 24, 12 | 0.882, 0.916, 0.990, **1.000** | 11.8, 12.5, 13.7, 11.4 σ |
| Titli | 48, 24, 12, 0 | 0.589, 0.679, 0.582, 0.631 | 3.3, 6.5, 2.9, 3.2 σ |
| Amphan | 72, 48 | 0.479, 0.516 | 0.9, 0.7 σ |

Four panels at exactly 1.000. Amphan is the control: symmetric in z_a, with mean
|z_c| of 0.64 among its zero-loss members against 1.44 among its loss-producing
ones, which is the correct signature for a storm that can only reach the
portfolio through a large cross-track displacement.

**The fix.** Each ensemble track continues past landfall along the real best
track, by enough *path length* to cover 4 σ_a of backward displacement. Sizing in
path length rather than hours matters, because a slow storm needs more hours to
cover the same ground.

σ is held **exactly constant** across the extension, making the post-landfall
segment a rigid translation of the real track by the member's landfall offset.
That is what a persistent along-track error means: the storm arrives late and
then follows its own path inland. Growing σ past landfall would fabricate
forecast error for a horizon that was never verified; shrinking it would pull
late members back onto the best track and quietly undo the displacement being
tested. Note that this is the frozen-σ failure mode of §6.5.1 applied
deliberately and only beyond landfall, so the assertion is split: σ must grow
strictly over the forecast window and be exactly flat over the extension.

Four σ leaves about three draws in 100,000 uncovered on the backward side, so the
expected number of still-truncated members at N = 500 is under 0.02 per panel.
Realised coverage is reported in `outputs/track_extension.csv` and
`outputs/extension_coverage.csv`: fourteen of fifteen positive-lead panels are
fully covered, and Phailin T-72h retains one member of 500 displaced 463 km
against the 422 km its best-track record can supply.

**The sign test is a screen, not a verdict.** One-sidedness is necessary for
truncation but not sufficient, because along-track displacement is a real
displacement. Where a storm's direction of motion points towards the exposure,
forward and backward shifts have genuinely different loss consequences and the
zeros are one-sided for physical reasons. Phailin, which made landfall moving
northwest up the coast towards the portfolio, is exactly that case. The
unconfounded test compares each panel's extension length against the largest
backward displacement actually drawn.

**What this cost.** Zero fractions were overstated, attachment probabilities
understated at every lead, and the along-track share of the variance
decomposition inflated, because the truncation switch is driven entirely by z_a.
Fani's T-24h along-track share fell from 0.926 to 0.714 after the fix. All Stage 4
to Stage 9 outputs were regenerated. The superseded diagnostic is retained at
`outputs/zero_loss_mechanism_PRE_FIX.csv`.

**What it says about the assertions.** All fourteen checks in
`generate_ensemble.py` passed on the broken version, because every one of them
tests a property of the ensemble at or before landfall, and the ensemble was
correct there. They verified that the sampler produced the specified
distribution. They could not verify that the specification was adequate for what
came next. `validate_ensemble` now asserts the extension against the realised
draw and refuses any positive-lead panel with no extension at all.

---

## 7. Two error models, not one

Restating blueprint §4.1 because it is the easiest thing in this project to get
wrong, and the consequence is invisible:

- **Landfall point** perturbation uses the **landfall point error**
  (10.5 / 16.2 / 39.3 / 69.5 km), decomposed into σ_a and σ_c by §6.5.
- **Approach path** positions use the **general track error**
  (— / 72 / 111 / 154 km), isotropically.

Forecasters constrain landfall far better than open-ocean position, and the gap
is widest exactly where the project's argument lives:

| Lead | Landfall error | Track error | Ratio |
|---:|---:|---:|---:|
| 24 h | 16.19 km | 72 km | 0.22 |
| 48 h | 39.32 km | 111 km | 0.35 |
| 72 h | 69.52 km | 154 km | 0.45 |

Using the general track error to perturb the landfall point would inflate
uncertainty by 4.4× at 24 h and 2.2× at 72 h. Because the distortion is largest
at short lead, it would not merely add noise — it would **manufacture the
project's hypothesis**, by making meteorological variance appear to dominate at
precisely the lead time where the project claims it does not.

**Assertion to code:** ensemble landfall spread at each lead time must match the
landfall error in this document to within a stated tolerance, dimension-aware per
§6.5.7. Fail the run otherwise. (Blueprint §8.2.) Implemented in
`validate_ensemble`.

---

## 8. Period choice and the contemporaneous problem

**Primary model: 2020–2024**, as specified in the blueprint. Recent enough to
reflect current skill, long enough to be stable at 12–72 h.

But three of the four events predate that window, and IMD's skill improved
sharply over the period:

| Event | Year | In primary window? |
|---|---:|---|
| Phailin | 2013 | No — precedes both stated windows |
| Titli | 2018 | No — 2015–2019 |
| Fani | 2019 | No — 2015–2019 |
| Amphan | 2020 | Yes |

The three windows, landfall point error:

| Lead | 2011–2015 | n | 2015–2019 | n | 2020–2024 | n |
|---:|---:|---:|---:|---:|---:|---:|
| 12 h | 36.53 km | 12 | 25.39 km | 17 | 10.53 km | 24 |
| 24 h | 56.34 km | 12 | 44.73 km | 17 | 16.19 km | 23 |
| 48 h | 93.49 km | 11 | 69.36 km | 13 | 39.32 km | 19 |
| 72 h | 105.70 km | 6 | 109.34 km | 7 | 69.52 km | 10 |

**This is not a cosmetic sensitivity.** At 24 h the contemporaneous error for
Phailin is 3.5× the primary model's. Running Phailin under 2020–2024 skill models
a forecast nobody could have issued in 2013.

**Decision:** all reported results use 2020–2024 for all four events, so the
events are compared on a common basis and the fan charts mean the same thing
across panels.

**The contemporaneous run was not performed.** It would use 2011–2015 for
Phailin, 2015–2019 for Titli and Fani, and 2020–2024 for Amphan. It is recorded
as future work in the README, and it is the most interesting of the four
outstanding sensitivities: Phailin's attachment probability peaks at 0.308 under
2020–2024 skill, and under 2013 skill its fan would be roughly three times wider
at 24 h, which could well carry it across the 0.50 notification threshold. That
would reframe the project's single false alarm as a result about forecast skill
improvement rather than a footnote.

State plainly which is which. The common-basis run answers *"how does forecast
uncertainty propagate to loss?"*; a contemporaneous run would answer *"what would
an event response team actually have faced at the time?"* Both are legitimate;
conflating them is not, and claiming to have run one when only the other exists
is worse than either.

---

## 9. Provenance

Every row of `data/imd_forecast_errors.csv` carries one of:

| Tag | Meaning | Where it appears |
|---|---|---|
| `imd_published` | Sample-weighted recomputation of IMD's published annual table (method validated, §2) | CSV |
| `assumed` | No published value exists; chosen and justified in §5 | CSV, T-0 row |
| `no_data` | No value is available to the sampler at this lead, either because zero cases were verified (120 h landfall) or because the published cells were rejected as unusable (96 h landfall, §4.2) | CSV |
| `rejected_small_sample` | Published but unusable. **Documented in §4.2, not carried in the CSV.** | this document |
| `REQUIRED_EXTRACTION` | Known gap, not yet filled | no longer any |

The two long-lead landfall cells are a deliberate asymmetry worth stating. The
120 h cell is `no_data` because IMD verified zero landfall forecasts at that lead
in 2020–2024. The 96 h cell had six, at 55.4 km, and that figure exists; it is
rejected on the grounds in §4.2 and recorded there rather than in the data file,
because `validate_landfall_ladder` requires it to be absent. Both therefore read
`no_data` in the CSV, and the distinction between "never measured" and "measured
and rejected" lives in §4.2 where it can be read and reasoned about.

This is the one place where the project records an exclusion in prose rather than
in data. The trade is deliberate: a value that is present and tagged can be read
by a future consumer that forgets to check the tag, and a value that is absent
and asserted cannot.

---

## 10. Known weaknesses

**10.1 The 72 h figure rests on ten cases, and one year dominates it.**
2021 contributes 158.5 km on n=2. Removing 2021 drops the 2020–2024 72 h mean
from 69.5 km to roughly 47 km, which would propagate to σ_c at 72 h and widen or
narrow every T-72h fan chart accordingly. **This sensitivity was not run.**
Recorded as future work. The same thin-sample fragility was recorded in Project 2
§10.9; it has not gone away, it has moved.

**10.2 Case counts are not equal across lead times.** n falls from 24 at 12 h to
10 at 72 h. The 72 h error statistic therefore carries more estimation
uncertainty than the 12 h one, and the ensemble spread arguably ought to be
widened at long lead to reflect uncertainty in the error parameter itself. Not
implemented. Stated.

**10.3 2025 is excluded.** The 2025 row is a single case with a 12 h error of
71 km against a 10.5 km five-year baseline. Including it would distort every
cell. The blueprint's 2020–2024 choice stands.

**10.4 Landfall errors are conditioned on a landfall forecast existing.** Every
figure here comes from storms IMD judged confident enough to forecast a landfall
for. Storms that never got one are absent from the denominator. This is the same
selection effect that disqualifies the long leads (§4.2); it is weaker at 12–72 h
but it is not zero.

**10.5 No bias decomposition.** Covered in §6.4. The anisotropy weakness
previously recorded here is closed by §6.5.

**10.6 The approach path is still isotropic.** §6.5 decomposes the landfall
phase only, because LPE exists only to 72 h. The 96 h and 120 h approach-path
positions are drawn with equal component widths. No loss is computed at those
leads, so this does not reach any reported result, but the approach-path
visualisation is isotropic where the landfall phase is not.

**10.7 Four sensitivities are specified and none has been run.** §5.1 (T-0
position at 0, 5, 10.53 km), §5.2 (T-0 intensity at 0 and 5 kt), §8
(contemporaneous skill windows) and §10.1 (72 h figure without 2021). All four
are recorded as future work in the README rather than reported as results.

---

## 11. Transcription status: closed

All previously open transcription items are complete. No published value from the
IMD error model is intentionally left unresolved for the primary approach-path and
landfall ladders.

- [x] **12 h general track error (DPE), 2020–2024:** 44.866834 km, n=398.
- [x] **96 h and 120 h general track error:** 182.694545 km (n=55) and
      245.225000 km (n=28); both satisfy the n≥10 screening rule and the
      approach-path error grows monotonically.
- [x] **12 h intensity error (MAE and RMSE), 2020–2024:** 3.384171 kt and
      4.960727 kt reconstructed from the annual table, consistent with IMD's
      published 3.4 kt and 5.0 kt values, n=398.
- [x] **Along/cross decomposition:** derived in §6.5.4 from published DPE and LPE
      jointly, closing the open item recorded in the original §6.4.

The track and MAE values were case-weighted as described in §2; the 12 h RMSE was
pooled in squared-error space as described in §3.3. The 24/48/72 h track and
landfall figures still reproduce IMD's published multi-year values, so the
aggregation logic remains reconciled.

The implementation and assertions this document calls for are in
`src/error_model.py` and `src/generate_ensemble.py`, and run as Stage 1 of the
pipeline.

---

## 12. Assertions this document implies

Implemented as executable checks, not left as prose.

1. Reconstructed 2020–2024 landfall point error at 24/48/72 h equals
   16.2 / 39.3 / 69.5 km to within 0.1 km. *(Guards the weighting method.)*
2. Reconstructed 2020–2024 landfall time error at 24/48/72 h equals
   2.9 / 4.2 / 7.5 h to within 0.1 h.
3. Every required approach-path and landfall lead time has a row in
   `imd_forecast_errors.csv` with a non-null value and a provenance tag.
4. The landfall ladder is empty beyond 72 h. `landfall_point_error_km` must be
   null at 96 h and 120 h, and a run that finds a value there fails rather than
   proceeding. *(This is the structural form of "no rejected cell ever reaches
   the sampler". It is stronger than tagging, because an absent value cannot be
   read by a consumer that forgets to check a tag, and it fails on the edit that
   would introduce the problem rather than on the run that suffers from it.)*
5. Cross-track width σ_c is strictly increasing across the landfall ladder
   0 → 12 → 24 → 48 → 72 h, and so is σ_a.
6. General track error is strictly increasing across the approach-path ladder
   12 → 24 → 48 → 72 → 96 → 120 h.
7. Sampled ensemble landfall spread at each landfall lead time matches the value
   in this document to within 3 Monte Carlo standard errors, compared against the
   cross-track component at 12–72 h and against the 2-D magnitude at T-0
   (§6.5.7). *(Blueprint §8.2.)*
8. Sampled total position spread at landfall matches the published DPE to within
   3 Monte Carlo standard errors.
9. Ensemble mean position at lead time T equals the best-track position at T,
   within 4 Monte Carlo standard errors. *(Guards against a biased sampler.)*
10. Positional σ values satisfy the conversion appropriate to their dimension,
    and both landfall bases are checked so neither can drift unnoticed:
    `track_sigma_km = 0.79788 × track_error_km`; the audit column
    `landfall_sigma_km = 0.79788 × LPE`; and the sampler's
    σ_c = 1.25331 × LPE at 12–72 h. T-0 is the exception, where the assumed
    5 km is a 2-D position uncertainty and all three coincide at 3.99 km (§6.1,
    §6.5.7).
11. Forcing σ_a = σ_c in the anisotropic solver reproduces the isotropic Rayleigh
    mean, σ·√(π/2), to within tolerance. *(Guards the angular quadrature; this is
    the check that caught the dropped Jacobian, §6.5.4.)*
12. σ grows strictly across the forecast window and is exactly constant across
    the post-landfall extension (§6.5.9).
13. The post-landfall extension is at least as long, in path length, as the
    largest backward along-track displacement actually drawn in that panel. Any
    positive-lead panel with no extension is refused outright.
14. Offsets are persistent along the track: correlation between a member's
    displacement at the first and last forecast horizon is 1.
15. Each member's applied offset is recovered from its resulting coordinates by
    inverse geodesic. *(Guards the heading convention; catches a left/right or
    along/cross swap.)*
16. The 12 h intensity reconstruction gives MAE = 3.384171 kt and pooled
    RMSE = 4.960727 kt, consistent with IMD's published 3.4/5.0 kt values.

Assertion 4 is the one that matters most. The failure this project is about is
wrong data producing confident, well-formatted, entirely plausible output, and a
rejected 96 h landfall cell reaching the sampler is exactly that failure: it
would widen every long-lead fan chart using a mean over six storms that were easy
enough to forecast four days out. The assertion does not test that the sampler
skips such a value. It tests that the value is not there to be skipped.

Assertion 13 exists because fourteen passing assertions did not catch the
truncation defect in §6.5.9. Every one of them tested the ensemble at or before
landfall, where it was correct.

---

## 13. Change log

| Date | Change |
|---|---|
| 2026-09-28 | Initial derivation. Added 12 h landfall statistics from IMD tables. Recorded T-0 assumptions (§5). Added 2011–2015 window for Phailin (§8). |
| 2026-09-28 | Completed 12 h general-track and intensity MAE/RMSE extraction; retained 96/120 h for the approach-path model after n≥10 and monotonic-growth checks; restricted the landfall ladder to 72 h and below. |
| 2026-09-30 | Derived the along/cross anisotropy from published DPE and LPE jointly (§6.5), closing the isotropy weakness recorded in the original §6.4. Corrected the landfall conversion from the 2-D Rayleigh factor to the 1-D half-normal factor (§6.1); the previous `landfall_sigma_km` column understated σ_c by 36%. Added the held-out translation-speed check (§6.5.6). |
| 2026-09-30 | Fixed a dropped Jacobian factor in the Hoyt angular integral, found by the equal-σ Rayleigh-limit assertion after an implied translation speed of 30 km/h was noticed and initially explained away. σ_a had been inflated by about 25%. |
| 2026-10-01 | Adopted the geometric sea-to-land landfall definition for all four events; all within 60 min of the IMD reference. Superseded the earlier Fani 06:00 UTC definition and the 20.2 N inland landfall point. Anchored the perturbation frame on the landfall heading after the local-heading frame produced 144 of 500 implausible members on Phailin T-72h (§6.5.3). |
| 2026-10-02 | Added §6.5.9, post-landfall track extension, after a one-sided z_a test showed zero-loss members were truncated rather than missing. Regenerated all Stage 4 to Stage 9 outputs. |
| 2026-10-02 | Corrected the T-0 row of `imd_forecast_errors.csv`, which carried zero track error in contradiction of §5.1. Retained the 96 h landfall cell as `rejected_small_sample` so assertion 4 has something to guard. Marked the four unrun sensitivities (§5.1, §5.2, §8, §10.1) as future work rather than reported results. Corrected the §4.3 wind-field count from 14,000 to 9,500. |
