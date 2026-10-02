"""
Project 3 - Stage 6: decision layer.

Converts each loss distribution into the three numbers an event response team
actually produces:

    1. A reserve recommendation (P75 of the loss distribution).
    2. The probability the CAT XL layer attaches.
    3. A notification decision: tell the reinsurer, or wait.

Attachment probability is computed from the 500 underlying losses directly,
never inferred from percentiles. P(loss > attachment) read off a fitted
distribution would be a second model layered on top of the first.

On the layer denomination
-------------------------
Project 1 prices its CAT XL in USD on a USD exposure base. Project 3 computes
losses in INR on Project 2's portfolio. The layer is rescaled by share of
exposure rather than copied as a number; config.describe_layer() prints the
derivation. The two portfolio magnitudes happen to be numerically close in
their respective currencies, so a naive currency copy would have produced
roughly the right answer for entirely the wrong reason. That is why the
derivation is in the code and printed on every run.
"""

from __future__ import annotations

import argparse

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import config as cfg


def reserve(losses: np.ndarray) -> float:
    return float(np.percentile(losses, cfg.RESERVE_PERCENTILE))


def attachment_probability(losses: np.ndarray) -> float:
    """Empirical P(loss > attachment). Strict '>' is intentional: a loss
    exactly at the attachment point recovers nothing."""
    return float(np.mean(losses > cfg.ATTACHMENT))


def recoveries(losses: np.ndarray) -> np.ndarray:
    """CAT XL recovery per member: min(max(loss - attachment, 0), limit)."""
    return np.minimum(np.maximum(losses - cfg.ATTACHMENT, 0.0), cfg.LIMIT)


def exhaustion_probability(losses: np.ndarray) -> float:
    return float(np.mean(losses >= cfg.ATTACHMENT + cfg.LIMIT))


def summarize(event: str, lead_h: int, losses: np.ndarray) -> dict:
    rec = recoveries(losses)
    p_attach = attachment_probability(losses)

    return {
        "event": event,
        "lead_time_h": lead_h,
        "members": int(len(losses)),
        "reserve_p75": reserve(losses),
        "attachment_point": cfg.ATTACHMENT,
        "limit": cfg.LIMIT,
        "attachment_probability": p_attach,
        "attachment_count": int(np.sum(losses > cfg.ATTACHMENT)),
        "exhaustion_probability": exhaustion_probability(losses),
        "notification_threshold": cfg.NOTIFICATION_THRESHOLD,
        "notification_trigger": bool(p_attach >= cfg.NOTIFICATION_THRESHOLD),
        "mean_gross_loss": float(np.mean(losses)),
        "median_gross_loss": float(np.median(losses)),
        "mean_cat_xl_recovery": float(np.mean(rec)),
        "median_cat_xl_recovery": float(np.median(rec)),
        "max_cat_xl_recovery": float(np.max(rec)),
    }


def validate(df: pd.DataFrame) -> None:
    if np.any(df["members"].to_numpy() != cfg.N_MEMBERS):
        raise AssertionError("Not every row has the expected member count.")

    p = df["attachment_probability"].to_numpy()
    if np.any((p < 0) | (p > 1)):
        raise AssertionError("Attachment probability outside [0, 1].")

    # The probability must be exactly count / members. If these ever disagree
    # it means the probability came from somewhere other than the raw losses.
    implied = df["attachment_count"].to_numpy(dtype=float) / cfg.N_MEMBERS
    if not np.allclose(p, implied):
        raise AssertionError(
            "Attachment probability does not equal count / members, so it was "
            "not computed from the underlying losses."
        )

    if np.any(df["reserve_p75"].to_numpy() < 0):
        raise AssertionError("Negative reserve.")

    rec = df["mean_cat_xl_recovery"].to_numpy()
    if np.any(rec < 0) or np.any(rec > cfg.LIMIT):
        raise AssertionError("Mean recovery outside [0, limit].")

    if np.any(df["exhaustion_probability"].to_numpy()
              > df["attachment_probability"].to_numpy() + 1e-12):
        raise AssertionError(
            "Exhaustion probability exceeds attachment probability, which is "
            "impossible: exhaustion implies attachment."
        )

    print("  members per lead              : PASSED")
    print("  attachment probability in 0-1 : PASSED")
    print("  probability = count / members : PASSED")
    print("  recovery within layer         : PASSED")
    print("  exhaustion <= attachment      : PASSED")


def attachment_chart(df: pd.DataFrame):
    fig, ax = plt.subplots(figsize=(9, 5.5))

    leads = sorted(df["lead_time_h"].unique(), reverse=True)
    x = np.arange(len(leads))

    for event in sorted(df["event"].unique()):
        sub = df[df["event"] == event].set_index("lead_time_h").reindex(leads)
        ax.plot(x, sub["attachment_probability"].to_numpy(),
                marker="o", linewidth=1.8, label=event)

    ax.axhline(cfg.NOTIFICATION_THRESHOLD, linestyle="--", linewidth=1.4,
               color="0.3",
               label=f"Notification threshold ({cfg.NOTIFICATION_THRESHOLD:.0%})")

    ax.set_xticks(x)
    ax.set_xticklabels([f"T-{int(v)}h" for v in leads])
    ax.set_ylim(0, 1)
    ax.set_xlabel("Forecast lead time")
    ax.set_ylabel(f"P(loss > {cfg.fmt_money(cfg.ATTACHMENT, 0)})")
    ax.set_title("CAT XL attachment probability by lead time")
    ax.grid(True, alpha=0.25)
    ax.legend()
    fig.tight_layout()

    path = cfg.STAGE6_DIR / "attachment_probability.png"
    fig.savefig(path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    return path


def reserve_chart(df: pd.DataFrame):
    fig, ax = plt.subplots(figsize=(9, 5.5))
    leads = sorted(df["lead_time_h"].unique(), reverse=True)
    x = np.arange(len(leads))

    for event in sorted(df["event"].unique()):
        sub = df[df["event"] == event].set_index("lead_time_h").reindex(leads)
        ax.plot(x, sub["reserve_p75"].to_numpy() / 1e9,
                marker="o", linewidth=1.8, label=event)

    ax.axhline(cfg.ATTACHMENT / 1e9, linestyle=":", linewidth=1.4, color="0.3",
               label="CAT XL attachment")
    ax.set_xticks(x)
    ax.set_xticklabels([f"T-{int(v)}h" for v in leads])
    ax.set_xlabel("Forecast lead time")
    ax.set_ylabel(f"P{cfg.RESERVE_PERCENTILE} reserve ({cfg.CURRENCY_CODE} billion)")
    ax.set_title(f"P{cfg.RESERVE_PERCENTILE} reserve recommendation by lead time")
    ax.grid(True, alpha=0.25)
    ax.legend()
    fig.tight_layout()

    path = cfg.STAGE6_DIR / "reserve_by_lead.png"
    fig.savefig(path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    return path


def run(events: tuple[str, ...] | None = None) -> pd.DataFrame:
    cfg.ensure_output_dirs()
    events = events or cfg.available_events()
    if not events:
        raise RuntimeError("No Stage 4 loss files found. Run Stage 4 first.")

    print("=" * 72)
    print("STAGE 6 - DECISION LAYER")
    print("=" * 72)
    print("\n[FROZEN DECISION RULES]")
    print(f"  Reserve percentile     : P{cfg.RESERVE_PERCENTILE}")
    print(f"  {cfg.describe_layer()}")
    print("  Attachment condition   : loss > attachment (strict)")
    print(f"  Notification threshold : {cfg.NOTIFICATION_THRESHOLD:.0%}")

    rows = []
    for event in events:
        leads = cfg.available_leads(event)
        if not leads:
            continue
        print(f"\n[{event}]")
        print(f"  {'lead':>6} {'P75 reserve':>16} {'P(attach)':>10} "
              f"{'attached':>9} {'mean recov':>16} {'notify':>7}")
        for lead in leads:
            s = summarize(event, lead, cfg.load_losses(event, lead))
            rows.append(s)
            print(f"  {'T-' + str(lead):>6} {cfg.fmt_money(s['reserve_p75']):>16} "
                  f"{s['attachment_probability']:10.4f} "
                  f"{s['attachment_count']:>4}/{s['members']:<4} "
                  f"{cfg.fmt_money(s['mean_cat_xl_recovery']):>16} "
                  f"{'YES' if s['notification_trigger'] else 'no':>7}")

    df = pd.DataFrame(rows).sort_values(
        ["event", "lead_time_h"], ascending=[True, False]
    ).reset_index(drop=True)

    print("\n[VALIDATION]")
    validate(df)

    summary_path = cfg.STAGE6_DIR / "decision_layer_summary.csv"
    df.to_csv(summary_path, index=False)

    charts = [attachment_chart(df), reserve_chart(df)]

    print("\n[OUTPUTS]")
    print(f"  {summary_path}")
    for c in charts:
        print(f"  {c}")

    print("\n" + "=" * 72)
    print("STAGE 6 COMPLETE")
    print("=" * 72)
    return df


def main():
    ap = argparse.ArgumentParser(description="Project 3 Stage 6.")
    ap.add_argument("--event", choices=cfg.EVENT_NAMES, action="append")
    args = ap.parse_args()
    run(tuple(args.event) if args.event else None)


if __name__ == "__main__":
    main()
