"""
Project 3 - zero-loss mechanism diagnostic.

At long lead a large share of ensemble members produce exactly zero loss:
31% of Fani members at T-72, 47% of Phailin members at T-72. Two mechanisms
could produce that, and they have opposite implications.

  CROSS-TRACK MISS (a result)
      The perturbed track makes landfall far enough from the portfolio that
      no cell sees wind above the vulnerability threshold. Zeros should then
      be roughly SYMMETRIC in z_a (about half negative) and concentrated at
      LARGE |z_c|.

  ALONG-TRACK TRUNCATION (a bug)
      A backward along-track displacement leaves the track ending before it
      reaches the coast, so no landfall is simulated at all. Zeros would then
      be ONE-SIDED in z_a, piling up at one sign, and would show little or no
      dependence on |z_c|.

The first is the finding the project is built on. The second would invalidate
every long-lead loss distribution in it, because the zeros would be an
artifact of how the ensemble was constructed rather than a statement about
where the storm might have gone.

The test is cheap because Stage 7 already reconstructs the latent draws
exactly from the seed, so no re-run of Stage 3 or Stage 4 is needed.

    python src/check_zero_mechanism.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

import config as cfg
from stage7_variance_decomposition import recover_latent_draws, validate_draw_order

# A zero set that is this lopsided in z_a is not a cross-track miss pattern.
ONE_SIDED_THRESHOLD = 0.70

# Panels with fewer zeros than this cannot support the test.
MIN_ZEROS = 25


def analyse(event: str, lead_h: int) -> dict | None:
    path = cfg.loss_file(event, lead_h)
    if not path.exists():
        return None

    df = pd.read_csv(path)
    if "member_index" in df.columns:
        df = df.sort_values("member_index")
    losses = df["gross_loss"].to_numpy(dtype=float)

    draws = recover_latent_draws(event, lead_h, n_members=len(losses))
    z_a = np.asarray(draws["along"], dtype=float)
    z_c = np.asarray(draws["cross"], dtype=float)

    is_zero = losses == 0.0
    n_zero = int(is_zero.sum())
    if n_zero < MIN_ZEROS or n_zero == len(losses):
        return dict(event=event, lead_h=lead_h, n_zero=n_zero, verdict="too few zeros")

    frac_neg = float((z_a[is_zero] < 0).mean())
    # Binomial standard error under the symmetric null.
    se = 0.5 / np.sqrt(n_zero)
    lopsided = max(frac_neg, 1.0 - frac_neg)

    abs_c_zero = float(np.abs(z_c[is_zero]).mean())
    abs_c_nonzero = float(np.abs(z_c[~is_zero]).mean())
    abs_a_zero = float(np.abs(z_a[is_zero]).mean())
    abs_a_nonzero = float(np.abs(z_a[~is_zero]).mean())

    if lopsided > ONE_SIDED_THRESHOLD and abs_c_zero < 1.15 * abs_c_nonzero:
        verdict = "TRUNCATION SUSPECTED"
    elif abs_c_zero > 1.25 * abs_c_nonzero:
        verdict = "cross-track miss"
    else:
        verdict = "unclear, inspect"

    return dict(
        event=event, lead_h=lead_h, n_zero=n_zero,
        frac_z_a_negative=frac_neg,
        sigmas_from_symmetric=abs(frac_neg - 0.5) / se,
        mean_abs_z_c_zero=abs_c_zero,
        mean_abs_z_c_nonzero=abs_c_nonzero,
        mean_abs_z_a_zero=abs_a_zero,
        mean_abs_z_a_nonzero=abs_a_nonzero,
        verdict=verdict,
    )


def main() -> int:
    validate_draw_order()

    print("=" * 78)
    print("PROJECT 3 - ZERO-LOSS MECHANISM")
    print("=" * 78)
    print("\nUnder a cross-track miss, zeros sit at large |z_c| and are split")
    print("evenly between negative and positive z_a. Under along-track")
    print("truncation, zeros pile up at one sign of z_a and |z_c| looks")
    print("the same as it does for the members that did produce loss.\n")

    rows = []
    for event in cfg.EVENT_NAMES:
        for lead in cfg.LEADS_H:
            rec = analyse(event, lead)
            if rec is not None:
                rows.append(rec)

    if not rows:
        print("No panels with enough zeros to test.")
        return 0

    df = pd.DataFrame(rows)
    testable = df[df["verdict"] != "too few zeros"].copy()

    if not testable.empty:
        show = testable[[
            "event", "lead_h", "n_zero", "frac_z_a_negative",
            "sigmas_from_symmetric", "mean_abs_z_c_zero",
            "mean_abs_z_c_nonzero", "verdict",
        ]]
        print(show.to_string(index=False, float_format=lambda x: f"{x:.3f}"))

    skipped = df[df["verdict"] == "too few zeros"]
    if not skipped.empty:
        print("\nNot testable (fewer than "
              f"{MIN_ZEROS} zeros, or all members zero):")
        for _, r in skipped.iterrows():
            print(f"  {r['event']} T-{r['lead_h']}: {r['n_zero']} zeros")

    out = cfg.OUTPUT_DIR / "zero_loss_mechanism.csv"
    df.to_csv(out, index=False)

    print("\n" + "=" * 78)
    suspect = testable[testable["verdict"] == "TRUNCATION SUSPECTED"]
    if len(suspect):
        print(f"{len(suspect)} panel(s) look one-sided in z_a. The zeros may be")
        print("an artifact of track length rather than a cross-track miss.")
        print("Before writing anything up, plot one zero-loss member's track")
        print("against the best track and check whether it reaches the coast.")
    else:
        print("No panel shows the truncation signature. The zeros behave like")
        print("cross-track misses, which is what the project claims they are.")
    print(f"\nWritten to {out}")
    print("=" * 78)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
