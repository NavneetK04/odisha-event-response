"""
Project 3 - Stage 8: notification trigger backtest.

Stage 6 fixes a notification rule: tell the reinsurer when P(loss >
attachment) reaches 50%. This stage asks whether that rule would actually have
worked, across all four events.

Why this stage is the one that matters
--------------------------------------
Fani, Phailin and Titli all made landfall inside the modelled Odisha
portfolio. Amphan did not: its forecast cone included Odisha at long lead and
its track resolved north into West Bengal. So Amphan is the only case that can
produce a false alarm, and a trigger that fires on all four events is not a
decision rule, it is a constant.

Two behaviours are scored separately:

    Warning time   - how early does the trigger fire on events that did hit?
    False alarms   - does it fire on the event that did not?

These trade off against each other, so the threshold sweep reports both rather
than collapsing them into one score. Collapsing them would bury the judgement
in an objective function and present it as an optimisation.

Whipsaw
-------
A trigger that fires at long lead and then stops firing as the forecast
resolves is the most expensive kind of wrong: the notification has already
gone out. That pattern is detected and reported per event, because it is
exactly what an event response team fears and exactly what Amphan is in the
set to demonstrate.
"""

from __future__ import annotations

import argparse

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import config as cfg
import stage6_decision_layer as s6


THRESHOLD_GRID = np.round(np.arange(0.05, 0.96, 0.05), 2)


def attachment_profile(event: str) -> pd.DataFrame:
    """P(attach) at every available lead for one event, ordered by lead."""
    rows = []
    for lead in cfg.available_leads(event):
        losses = cfg.load_losses(event, lead)
        rows.append({
            "event": event,
            "lead_time_h": lead,
            "attachment_probability": s6.attachment_probability(losses),
            "mean_loss": float(np.mean(losses)),
            "median_loss": float(np.median(losses)),
            "p_zero_loss": float(np.mean(losses <= 1e3)),
            "hit_portfolio": cfg.EVENT_HIT_PORTFOLIO[event],
        })
    return pd.DataFrame(rows).sort_values("lead_time_h", ascending=False)


def evaluate_threshold(profiles: dict[str, pd.DataFrame],
                       threshold: float) -> dict:
    """Score one threshold across all events.

    A 'warning' is the earliest lead at which the trigger fires. For an event
    that hit, earlier is better. For an event that did not hit, firing at all
    is a false alarm.
    """
    warned_hits, missed_hits, false_alarms, correct_silence = [], [], [], []
    warning_leads = {}

    for event, prof in profiles.items():
        fired = prof["attachment_probability"].to_numpy() >= threshold
        leads = prof["lead_time_h"].to_numpy()
        hit = bool(prof["hit_portfolio"].iloc[0])

        if fired.any():
            first = int(leads[np.argmax(fired)])
            warning_leads[event] = first
            (warned_hits if hit else false_alarms).append(event)
        else:
            warning_leads[event] = None
            (missed_hits if hit else correct_silence).append(event)

    hit_events = [e for e, p in profiles.items() if bool(p["hit_portfolio"].iloc[0])]
    lead_values = [warning_leads[e] for e in warned_hits if warning_leads[e] is not None]

    return {
        "threshold": float(threshold),
        "warned_hits": len(warned_hits),
        "missed_hits": len(missed_hits),
        "false_alarms": len(false_alarms),
        "correct_silence": len(correct_silence),
        "hit_detection_rate": len(warned_hits) / len(hit_events) if hit_events else np.nan,
        "mean_warning_lead_h": float(np.mean(lead_values)) if lead_values else np.nan,
        "min_warning_lead_h": float(np.min(lead_values)) if lead_values else np.nan,
        "false_alarm_events": ",".join(sorted(false_alarms)),
        "missed_events": ",".join(sorted(missed_hits)),
    }


def detect_whipsaw(prof: pd.DataFrame, threshold: float) -> dict:
    """Did the trigger fire and then stop firing as the forecast resolved?

    Leads run from long to short, so a True followed later by a False is a
    notification that was sent and then contradicted.
    """
    fired = prof["attachment_probability"].to_numpy() >= threshold
    leads = prof["lead_time_h"].to_numpy()

    whipsaw, first_fire, stopped_at = False, None, None
    for i, f in enumerate(fired):
        if f and first_fire is None:
            first_fire = int(leads[i])
        if first_fire is not None and not f:
            whipsaw, stopped_at = True, int(leads[i])
            break

    return {
        "event": str(prof["event"].iloc[0]),
        "threshold": float(threshold),
        "first_fire_lead_h": first_fire,
        "whipsaw": whipsaw,
        "stopped_firing_at_lead_h": stopped_at,
        "hit_portfolio": bool(prof["hit_portfolio"].iloc[0]),
    }


def profile_chart(profiles: dict[str, pd.DataFrame]):
    fig, ax = plt.subplots(figsize=(9.5, 5.5))
    leads = sorted({l for p in profiles.values() for l in p["lead_time_h"]},
                   reverse=True)
    x = np.arange(len(leads))

    for event, prof in sorted(profiles.items()):
        sub = prof.set_index("lead_time_h").reindex(leads)
        hit = bool(prof["hit_portfolio"].iloc[0])
        ax.plot(x, sub["attachment_probability"].to_numpy(),
                marker="o" if hit else "^",
                linestyle="-" if hit else "--",
                linewidth=2.0 if hit else 1.6,
                label=(
    f"{event} (struck portfolio)" if event == "Fani"
    else f"{event} (near miss)" if event == "Phailin"
    else f"{event} (missed)"
))

    ax.axhline(cfg.NOTIFICATION_THRESHOLD, color="0.3", linestyle=":",
               linewidth=1.4,
               label=f"Threshold ({cfg.NOTIFICATION_THRESHOLD:.0%})")
    ax.set_xticks(x)
    ax.set_xticklabels([f"T-{int(v)}h" for v in leads])
    ax.set_ylim(0, 1)
    ax.set_xlabel("Forecast lead time")
    ax.set_ylabel("P(loss > attachment)")
    ax.set_title("Notification trigger backtest: Phailin is the test")
    ax.grid(True, alpha=0.25)
    ax.legend(fontsize=9)
    fig.tight_layout()

    path = cfg.STAGE8_DIR / "trigger_profiles.png"
    fig.savefig(path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    return path


def sweep_chart(sweep: pd.DataFrame):
    fig, ax1 = plt.subplots(figsize=(9, 5.5))

    ax1.plot(
        sweep["threshold"],
        sweep["mean_warning_lead_h"],
        marker="o",
        linewidth=1.8,
        label="Mean warning lead time",
    )
    ax1.set_xlabel("Notification threshold, P(attach)")
    ax1.set_ylabel("Mean warning lead time (hours)")
    ax1.set_ylim(
        max(0, sweep["mean_warning_lead_h"].min() - 6),
        sweep["mean_warning_lead_h"].max() + 6,
    )
    ax1.grid(True, alpha=0.25)

    ax2 = ax1.twinx()
    ax2.plot(
        sweep["threshold"],
        sweep["false_alarms"],
        marker="s",
        linestyle="--",
        linewidth=1.6,
        color="0.35",
        label="False alarms",
    )
    ax2.set_ylabel("False alarms (count)")
    ax2.set_ylim(
        -0.1,
        max(1.2, sweep["false_alarms"].max() + 0.3),
    )

    ax1.axvline(
        cfg.NOTIFICATION_THRESHOLD,
        color="0.3",
        linestyle=":",
        linewidth=1.4,
        label=f"Frozen rule ({cfg.NOTIFICATION_THRESHOLD:.0%})",
    )

    lines = ax1.get_lines() + ax2.get_lines()
    ax1.legend(
        lines,
        [line.get_label() for line in lines],
        loc="center right",
        fontsize=9,
    )

    ax1.set_title("Threshold sweep: warning lead time against false alarms")
    fig.tight_layout()

    path = cfg.STAGE8_DIR / "threshold_sweep.png"
    fig.savefig(path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    return path

def run(events: tuple[str, ...] | None = None) -> pd.DataFrame:
    cfg.ensure_output_dirs()
    events = events or cfg.available_events()
    if not events:
        raise RuntimeError("No Stage 4 loss files found. Run Stage 4 first.")

    print("=" * 72)
    print("STAGE 8 - NOTIFICATION TRIGGER BACKTEST")
    print("=" * 72)

    profiles = {e: attachment_profile(e) for e in events
                if len(cfg.available_leads(e)) > 0}

    hits = [e for e in profiles if cfg.EVENT_HIT_PORTFOLIO[e]]
    misses = [e for e in profiles if not cfg.EVENT_HIT_PORTFOLIO[e]]

    print(f"\nEvents that hit the portfolio : {', '.join(hits) or 'none'}")
    print(f"Near misses                   : {', '.join(misses) or 'none'}")
    if not misses:
        print("\n  WARNING: no near-miss event is present. Without one, a")
        print("  trigger that fires on everything scores perfectly and the")
        print("  backtest cannot distinguish a rule from a constant.")

    print(f"\n[ATTACHMENT PROBABILITY BY LEAD, threshold "
          f"{cfg.NOTIFICATION_THRESHOLD:.0%}]")
    for event, prof in sorted(profiles.items()):
        tag = "" if cfg.EVENT_HIT_PORTFOLIO[event] else "  (near miss)"
        print(f"\n  {event}{tag}")
        print(f"    {'lead':>6} {'P(attach)':>10} {'fires':>7} {'P(loss=0)':>10}")
        for _, r in prof.iterrows():
            fires = r["attachment_probability"] >= cfg.NOTIFICATION_THRESHOLD
            print(f"    {'T-' + str(int(r['lead_time_h'])):>6} "
                  f"{r['attachment_probability']:10.4f} "
                  f"{'YES' if fires else 'no':>7} {r['p_zero_loss']:10.3f}")

    # Whipsaw at the frozen threshold.
    whip = pd.DataFrame([detect_whipsaw(p, cfg.NOTIFICATION_THRESHOLD)
                         for p in profiles.values()])
    print("\n[WHIPSAW CHECK at the frozen threshold]")
    print("  A trigger that fires at long lead and then stops firing is a")
    print("  notification already sent and then contradicted.")
    for _, r in whip.iterrows():
        if r["whipsaw"]:
            print(f"  {r['event']:<9} fired at "
                  f"T-{int(r['first_fire_lead_h'])}h, stopped by "
                  f"T-{int(r['stopped_firing_at_lead_h'])}h  <-- WHIPSAW"
                  f"{'' if r['hit_portfolio'] else ' (on the near miss)'}")
        elif pd.notna(r["first_fire_lead_h"]):
            print(f"  {r['event']:<9} fired at "
                  f"T-{int(r['first_fire_lead_h'])}h and stayed fired")
        else:
            print(f"  {r['event']:<9} never fired")

    # Threshold sweep.
    sweep = pd.DataFrame([evaluate_threshold(profiles, t) for t in THRESHOLD_GRID])

    print("\n[THRESHOLD SWEEP]")
    print(f"  {'thresh':>7} {'detected':>9} {'missed':>7} {'false':>6} "
          f"{'mean warn':>10}  false alarm on")
    for _, r in sweep.iterrows():
        warn = (f"{r['mean_warning_lead_h']:8.1f}h"
                if np.isfinite(r["mean_warning_lead_h"]) else f"{'n/a':>9}")
        print(f"  {r['threshold']:7.2f} {r['warned_hits']:>4}/{len(hits):<4} "
              f"{r['missed_hits']:>7} {r['false_alarms']:>6} {warn}  "
              f"{r['false_alarm_events']}")

    clean = sweep[(sweep["false_alarms"] == 0) & (sweep["missed_hits"] == 0)]
    print("\n[READING]")
    if len(clean):
        best = clean.loc[clean["mean_warning_lead_h"].idxmax()]
        print(f"  Thresholds from {clean['threshold'].min():.2f} to "
              f"{clean['threshold'].max():.2f} detect every hit with no false "
              f"alarm.")
        print(f"  Earliest mean warning in that band: "
              f"{best['mean_warning_lead_h']:.1f} h at threshold "
              f"{best['threshold']:.2f}.")
    else:
        print("  No threshold both catches every hit and avoids every false")
        print("  alarm. That is a result, not a failure: it means the trigger")
        print("  cannot separate this near miss from the direct hits on loss")
        print("  probability alone, and the rule needs another input.")

    prof_path = cfg.STAGE8_DIR / "attachment_profiles.csv"
    pd.concat(profiles.values(), ignore_index=True).to_csv(prof_path, index=False)

    sweep_path = cfg.STAGE8_DIR / "threshold_sweep.csv"
    sweep.to_csv(sweep_path, index=False)

    whip_path = cfg.STAGE8_DIR / "whipsaw_check.csv"
    whip.to_csv(whip_path, index=False)

    charts = [profile_chart(profiles), sweep_chart(sweep)]

    print("\n[OUTPUTS]")
    for p in (prof_path, sweep_path, whip_path, *charts):
        print(f"  {p}")

    print("\n" + "=" * 72)
    print("STAGE 8 COMPLETE")
    print("=" * 72)
    return sweep


def main():
    ap = argparse.ArgumentParser(description="Project 3 Stage 8.")
    ap.add_argument("--event", choices=cfg.EVENT_NAMES, action="append")
    args = ap.parse_args()
    run(tuple(args.event) if args.event else None)


if __name__ == "__main__":
    main()
