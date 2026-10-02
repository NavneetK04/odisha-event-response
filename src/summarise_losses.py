"""
Project 3 - compact loss summary for review.

Reduces every PRIMARY loss file to one row so the whole Stage 4 result can be
checked without moving 9,500 rows around. Computed straight from the raw
losses, so comparing this against loss_distribution_summary.csv and
decision_layer_summary.csv is a real cross-check: those come from Stage 5 and
Stage 6, this does not.

Also lists which exposure scenario variants exist per panel, which is what
Stage 7 discovers by glob, so a missing variant shows up here rather than as a
quietly thinner decomposition.

    python src/summarise_losses.py
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

import config as cfg

PRIMARY_RE = re.compile(r"^([a-z]+)_T(\d+)_event_losses\.csv$")
VARIANT_RE = re.compile(r"^([a-z]+)_T(\d+)_(.+)_event_losses\.csv$")


def main() -> int:
    rows, variants = [], {}

    for p in sorted(cfg.IMPACT_DIR.glob("*_event_losses.csv")):
        mv = VARIANT_RE.match(p.name)
        if mv:
            key = (mv.group(1).capitalize(), int(mv.group(2)))
            variants.setdefault(key, []).append(mv.group(3))
            continue

        mp = PRIMARY_RE.match(p.name)
        if not mp:
            continue
        event, lead = mp.group(1).capitalize(), int(mp.group(2))

        df = pd.read_csv(p)
        if "member_index" in df.columns:
            df = df.sort_values("member_index")
        L = df["gross_loss"].to_numpy(dtype=float)

        rows.append({
            "event": event,
            "lead_time_h": lead,
            "n_members": len(L),
            "p_zero": float((L == 0).mean()),
            "mean_B": L.mean() / 1e9,
            "p5_B": np.percentile(L, 5) / 1e9,
            "p50_B": np.median(L) / 1e9,
            "p75_B": np.percentile(L, 75) / 1e9,
            "p95_B": np.percentile(L, 95) / 1e9,
            "max_B": L.max() / 1e9,
            "p_attach": float((L > cfg.ATTACHMENT).mean()),
            "n_attach": int((L > cfg.ATTACHMENT).sum()),
            "p_exhaust": float((L > cfg.ATTACHMENT + cfg.LIMIT).mean()),
        })

    if not rows:
        print(f"No PRIMARY loss files found in {cfg.IMPACT_DIR}")
        return 1

    df = pd.DataFrame(rows).sort_values(
        ["event", "lead_time_h"], ascending=[True, False]
    )
    df["variants"] = [
        ";".join(sorted(variants.get((r.event, r.lead_time_h), [])))
        for r in df.itertuples()
    ]

    out = cfg.OUTPUT_DIR / "loss_summary_check.csv"
    df.to_csv(out, index=False)

    pd.set_option("display.width", 220)
    print("=" * 110)
    print("PROJECT 3 - LOSS SUMMARY (computed from raw member losses)")
    print("=" * 110)
    print(cfg.describe_layer())
    print()
    print(
        df[["event", "lead_time_h", "p_zero", "mean_B", "p50_B", "p75_B",
            "max_B", "p_attach"]].to_string(
            index=False, float_format=lambda x: f"{x:.3f}"
        )
    )

    print("\nExposure scenario variants present per panel:")
    for r in df.itertuples():
        n = len(variants.get((r.event, r.lead_time_h), []))
        print(f"  {r.event:<8} T-{r.lead_time_h:<3} {n} variant(s): {r.variants}")

    print(f"\nWritten to {out}")
    print("=" * 110)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
