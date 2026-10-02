"""
Project 3 - Stage 5: loss distribution.

Turns each Stage 4 set of 500 simulated losses into a distribution, for every
event and every lead time in the loss-phase ladder.

Two things this does that the Fani-only version did not
-------------------------------------------------------
1. It is event-aware. The previous version hardcoded 'fani_T{lead}' in the
   filename and 'Fani' in the summary row, which is the only reason Phailin
   failed. Nothing about Phailin was wrong.

2. It reports RELATIVE interval width alongside absolute width.

   The absolute P95-P5 width barely narrows toward landfall, because the
   distribution shifts right as members converge on the true track at the
   same time as it narrows. Measured in absolute rupees the collapse looks
   like 18%; normalised by the median it is a factor of five. Both are true
   and they answer different questions: the relative width is "how much does
   the uncertainty collapse", the absolute width is "how much must an insurer
   reserve against". Reporting only the first would overstate the result and
   reporting only the second would hide it.

P95/P5 is undefined whenever P5 = 0, which happens at long lead because more
than 5% of plausible tracks miss the portfolio entirely. Rather than invent an
infinite ratio, that case is reported as P(loss = 0), which is the more
informative quantity and is the cone question expressed in loss terms.
"""

from __future__ import annotations

import argparse

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import config as cfg


QUANTILES = (5, 25, 50, 75, 95)

# A member is treated as producing no loss if its loss is below this. Exactly
# zero is the normal case (storm misses the portfolio); the floor guards
# against a vanishing but non-zero wind field at the edge of the grid.
ZERO_LOSS_FLOOR = 1e3


def summarize(event: str, lead_h: int, losses: np.ndarray) -> dict:
    p5, p25, p50, p75, p95 = np.percentile(losses, QUANTILES)

    if not (p5 <= p25 <= p50 <= p75 <= p95):
        raise AssertionError(
            f"{event} T-{lead_h}: percentiles are not monotonic."
        )
    if p5 < 0:
        raise AssertionError(f"{event} T-{lead_h}: negative P5.")

    width = p95 - p5
    zero_fraction = float(np.mean(losses <= ZERO_LOSS_FLOOR))

    # Relative width. Undefined if the median is zero, which means more than
    # half the ensemble misses the portfolio - itself a reportable result.
    rel_width = width / p50 if p50 > 0 else np.nan

    # Retained for continuity with the earlier Fani run, but it is undefined
    # whenever P5 = 0 and zero_fraction is the better statistic.
    p95_p5 = p95 / p5 if p5 > 0 else np.nan

    return {
        "event": event,
        "lead_time_h": lead_h,
        "n_members": int(len(losses)),
        "p5": float(p5),
        "p25": float(p25),
        "p50": float(p50),
        "p75": float(p75),
        "p95": float(p95),
        "mean": float(np.mean(losses)),
        "minimum": float(np.min(losses)),
        "maximum": float(np.max(losses)),
        "width_p95_p5": float(width),
        "relative_width": float(rel_width) if np.isfinite(rel_width) else np.nan,
        "p_zero_loss": zero_fraction,
        "p95_p5_ratio": float(p95_p5) if np.isfinite(p95_p5) else np.nan,
    }


def report_convergence(df: pd.DataFrame, event: str) -> dict:
    """Report how the interval narrows, absolutely and relatively.

    The blueprint asks us to TEST monotonic narrowing, not to enforce it. A
    non-monotonic absolute width is a diagnostic, not an error, because the
    distribution is shifting as well as narrowing.
    """
    sub = df[df["event"] == event].sort_values("lead_time_h", ascending=False)
    widths = sub["width_p95_p5"].to_numpy()
    rel = sub["relative_width"].to_numpy()
    leads = sub["lead_time_h"].to_numpy()

    abs_monotonic = bool(np.all(np.diff(widths) <= 1e-9))
    finite = np.isfinite(rel)
    rel_monotonic = bool(np.all(np.diff(rel[finite]) <= 1e-9)) if finite.sum() > 1 else False

    abs_collapse = widths[0] / widths[-1] if widths[-1] > 0 else np.nan
    rel_collapse = (
        rel[finite][0] / rel[finite][-1]
        if finite.sum() > 1 and rel[finite][-1] > 0
        else np.nan
    )

    print(f"\n  [{event}] interval width, T-{int(leads[0])}h -> T-{int(leads[-1])}h")
    print(f"  {'lead':>6} {'abs width':>16} {'median':>16} {'rel width':>10} {'P(loss=0)':>10}")
    for _, r in sub.iterrows():
        rw = f"{r['relative_width']:10.2f}" if np.isfinite(r["relative_width"]) else f"{'n/a':>10}"
        print(f"  {'T-' + str(int(r['lead_time_h'])):>6} "
              f"{cfg.fmt_money(r['width_p95_p5']):>16} "
              f"{cfg.fmt_money(r['p50']):>16} {rw} {r['p_zero_loss']:10.3f}")

    print(f"  absolute narrowing  : {abs_collapse:.2f}x "
          f"({'monotonic' if abs_monotonic else 'NOT monotonic'})")
    if np.isfinite(rel_collapse):
        print(f"  relative narrowing  : {rel_collapse:.2f}x "
              f"({'monotonic' if rel_monotonic else 'NOT monotonic'})")
    else:
        print("  relative narrowing  : undefined (median reaches zero)")

    return {
        "event": event,
        "absolute_collapse": float(abs_collapse) if np.isfinite(abs_collapse) else np.nan,
        "relative_collapse": float(rel_collapse) if np.isfinite(rel_collapse) else np.nan,
        "absolute_monotonic": abs_monotonic,
        "relative_monotonic": rel_monotonic,
    }


def fan_chart(df: pd.DataFrame, event: str):
    sub = df[df["event"] == event].sort_values("lead_time_h", ascending=False)
    x = np.arange(len(sub))

    fig, ax = plt.subplots(figsize=(9, 5.5))
    ax.fill_between(x, sub["p5"] / 1e9, sub["p95"] / 1e9, alpha=0.20, label="P5 to P95")
    ax.fill_between(x, sub["p25"] / 1e9, sub["p75"] / 1e9, alpha=0.35, label="P25 to P75")
    ax.plot(x, sub["p50"] / 1e9, marker="o", linewidth=2, label="Median")
    ax.plot(x, sub["mean"] / 1e9, marker="s", linestyle="--", linewidth=1.4, label="Mean")

    ax.set_xticks(x)
    ax.set_xticklabels([f"T-{int(v)}h" for v in sub["lead_time_h"]])
    ax.set_xlabel("Forecast lead time")
    ax.set_ylabel(f"Gross loss ({cfg.CURRENCY_CODE} billion)")
    ax.set_title(f"{event}: forecast loss distribution by lead time")
    ax.grid(True, alpha=0.25)
    ax.legend()
    fig.tight_layout()

    path = cfg.STAGE5_DIR / f"{cfg.event_key(event)}_loss_fan_chart.png"
    fig.savefig(path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    return path


def combined_chart(df: pd.DataFrame):
    """All four events on one relative-width axis: the convergence result."""
    fig, ax = plt.subplots(figsize=(9, 5.5))

    for event in sorted(df["event"].unique()):
        sub = df[df["event"] == event].sort_values("lead_time_h", ascending=False)
        rel = sub["relative_width"].to_numpy()
        finite = np.isfinite(rel)
        if finite.sum() < 2:
            continue
        ax.plot(np.arange(len(sub))[finite], rel[finite],
                marker="o", linewidth=1.8, label=event)

    leads = sorted(df["lead_time_h"].unique(), reverse=True)
    ax.set_xticks(np.arange(len(leads)))
    ax.set_xticklabels([f"T-{int(v)}h" for v in leads])
    ax.set_xlabel("Forecast lead time")
    ax.set_ylabel("(P95 - P5) / median")
    ax.set_title("Relative loss uncertainty collapses as the decision clock runs down")
    ax.grid(True, alpha=0.25)
    ax.legend()
    fig.tight_layout()

    path = cfg.STAGE5_DIR / "all_events_relative_width.png"
    fig.savefig(path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    return path


def run(events: tuple[str, ...] | None = None) -> pd.DataFrame:
    cfg.ensure_output_dirs()

    events = events or cfg.available_events()
    if not events:
        raise RuntimeError(
            "No Stage 4 loss files found. Run Stage 4 before Stage 5."
        )

    print("=" * 72)
    print("STAGE 5 - LOSS DISTRIBUTION")
    print("=" * 72)

    rows = []
    for event in events:
        leads = cfg.available_leads(event)
        if not leads:
            print(f"\n[{event}] no loss files, skipped")
            continue
        print(f"\n[{event}] leads: {', '.join('T-' + str(l) for l in leads)}")
        for lead in leads:
            losses = cfg.load_losses(event, lead)
            rows.append(summarize(event, lead, losses))

    df = pd.DataFrame(rows)

    convergence = [report_convergence(df, e) for e in df["event"].unique()]

    summary_path = cfg.STAGE5_DIR / "loss_distribution_summary.csv"
    df.to_csv(summary_path, index=False)

    conv_path = cfg.STAGE5_DIR / "convergence_summary.csv"
    pd.DataFrame(convergence).to_csv(conv_path, index=False)

    charts = [fan_chart(df, e) for e in df["event"].unique()]
    charts.append(combined_chart(df))

    print("\n[OUTPUTS]")
    print(f"  {summary_path}")
    print(f"  {conv_path}")
    for c in charts:
        print(f"  {c}")

    print("\n" + "=" * 72)
    print("STAGE 5 COMPLETE")
    print("=" * 72)
    return df


def main():
    ap = argparse.ArgumentParser(description="Project 3 Stage 5.")
    ap.add_argument("--event", choices=cfg.EVENT_NAMES, action="append",
                    help="Restrict to one or more events. Default: all available.")
    args = ap.parse_args()
    run(tuple(args.event) if args.event else None)


if __name__ == "__main__":
    main()
