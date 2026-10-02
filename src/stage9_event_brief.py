"""
Project 3 - Stage 9: event response brief.

The showcase deliverable. Produces the one-page brief an event response
analyst would actually send: what we think the loss is, how uncertain that
is, whether the CAT XL layer is in play, and what we recommend doing right
now.

Design rule for this page
-------------------------
Every number on it must trace to a committed CSV. Nothing is retyped and
nothing is rounded into the HTML by hand, because a brief is exactly the kind
of artefact where a number gets prettified once and then quoted forever.

The brief also carries its own uncertainty statement. A loss estimate issued
at T-72h with a plausible range spanning an order of magnitude is not the same
product as the same number issued at T-12h, and a brief that does not say so
is misleading even when every figure in it is correct.
"""

from __future__ import annotations

import argparse
import html
from datetime import datetime, timezone

import numpy as np
import pandas as pd

import config as cfg
import stage5_loss_distribution as s5
import stage6_decision_layer as s6


def _recommendation(dist: dict, dec: dict, lead_h: int) -> tuple[str, str]:
    """Recommendation and its one-line rationale."""
    p = dec["attachment_probability"]
    rel = dist["relative_width"]

    if dec["notification_trigger"]:
        action = "NOTIFY REINSURER"
        why = (f"Attachment probability {p:.0%} is at or above the "
               f"{cfg.NOTIFICATION_THRESHOLD:.0%} notification threshold.")
    elif p >= cfg.NOTIFICATION_THRESHOLD * 0.5:
        action = "MONITOR CLOSELY"
        why = (f"Attachment probability {p:.0%} is material but below the "
               f"{cfg.NOTIFICATION_THRESHOLD:.0%} threshold.")
    else:
        action = "NO ACTION"
        why = (f"Attachment probability {p:.0%} is well below the "
               f"{cfg.NOTIFICATION_THRESHOLD:.0%} threshold.")

    if np.isfinite(rel) and rel > 3.0 and lead_h >= 48:
        why += (f" Note the plausible range is still {rel:.1f} times the "
                f"central estimate at this lead time.")
    return action, why


def build_brief(event: str, lead_h: int,
                dist_df: pd.DataFrame, dec_df: pd.DataFrame,
                var_df: pd.DataFrame | None = None) -> str:
    ev = cfg.load_events()[event]

    d = dist_df[(dist_df["event"] == event) &
                (dist_df["lead_time_h"] == lead_h)]
    k = dec_df[(dec_df["event"] == event) &
               (dec_df["lead_time_h"] == lead_h)]
    if d.empty or k.empty:
        raise ValueError(f"No Stage 5/6 results for {event} T-{lead_h}h.")

    dist, dec = d.iloc[0].to_dict(), k.iloc[0].to_dict()
    action, why = _recommendation(dist, dec, lead_h)

    landfall = str(ev.landfall_time)[:16].replace("T", " ")
    origin = str(ev.landfall_time - np.timedelta64(lead_h, "h"))[:16].replace("T", " ")

    rel = (f"{dist['relative_width']:.2f}x"
           if np.isfinite(dist["relative_width"]) else "undefined")

    var_rows = ""
    if var_df is not None:
        v = var_df[(var_df["event"] == event) &
                   (var_df["lead_time_h"] == lead_h)]
        for _, r in v.sort_values("omega_sq", ascending=False).iterrows():
            var_rows += (f"<tr><td>{html.escape(str(r['label']))}</td>"
                         f"<td class='n'>{r['omega_sq']:.1%}</td></tr>")

    generated = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")

    return f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>{html.escape(event)} T-{lead_h}h Event Response Brief</title>
<style>
 :root {{ --ink:#15181d; --muted:#5b6472; --line:#d9dee6; --bg:#ffffff;
          --accent:#1d4ed8; --warn:#b45309; }}
 * {{ box-sizing:border-box; }}
 body {{ margin:0; padding:28px 16px; background:var(--bg); color:var(--ink);
        font:15px/1.55 -apple-system,BlinkMacSystemFont,"Segoe UI",Roboto,sans-serif; }}
 .page {{ max-width:860px; margin:0 auto; }}
 header {{ border-bottom:2px solid var(--ink); padding-bottom:12px; margin-bottom:20px; }}
 h1 {{ font-size:22px; margin:0 0 4px; letter-spacing:-0.01em; }}
 .sub {{ color:var(--muted); font-size:13px; }}
 .action {{ border:2px solid var(--accent); border-radius:6px; padding:14px 16px;
            margin:18px 0; }}
 .action .label {{ font-size:12px; text-transform:uppercase; letter-spacing:.08em;
                   color:var(--muted); }}
 .action .verdict {{ font-size:20px; font-weight:700; color:var(--accent);
                     margin:4px 0 6px; }}
 h2 {{ font-size:13px; text-transform:uppercase; letter-spacing:.08em;
       color:var(--muted); margin:22px 0 8px; border-bottom:1px solid var(--line);
       padding-bottom:5px; }}
 table {{ width:100%; border-collapse:collapse; font-size:14px; }}
 td {{ padding:6px 0; border-bottom:1px solid var(--line); }}
 td.n {{ text-align:right; font-variant-numeric:tabular-nums; font-weight:600; }}
 .grid {{ display:grid; grid-template-columns:repeat(auto-fit,minmax(200px,1fr));
          gap:14px; margin:14px 0; }}
 .stat {{ border:1px solid var(--line); border-radius:5px; padding:11px 13px; }}
 .stat .k {{ font-size:11px; text-transform:uppercase; letter-spacing:.06em;
             color:var(--muted); }}
 .stat .v {{ font-size:19px; font-weight:700; font-variant-numeric:tabular-nums;
             margin-top:3px; }}
 .caveat {{ background:#fffbeb; border-left:3px solid var(--warn);
            padding:10px 14px; font-size:13px; margin:16px 0; }}
 footer {{ margin-top:26px; padding-top:12px; border-top:1px solid var(--line);
           font-size:11.5px; color:var(--muted); }}
 @media print {{ body {{ padding:0; }} }}
</style></head><body><div class="page">

<header>
 <h1>{html.escape(event)} &mdash; Event Response Brief, T-{lead_h}h</h1>
 <div class="sub">Forecast issued {origin} UTC &nbsp;&middot;&nbsp;
  Projected landfall {landfall} UTC &nbsp;&middot;&nbsp;
  {cfg.N_MEMBERS}-member forecast ensemble</div>
</header>

<div class="action">
 <div class="label">Recommendation</div>
 <div class="verdict">{html.escape(action)}</div>
 <div>{html.escape(why)}</div>
</div>

<h2>Loss estimate</h2>
<div class="grid">
 <div class="stat"><div class="k">Central estimate (median)</div>
   <div class="v">{cfg.fmt_money(dist['p50'])}</div></div>
 <div class="stat"><div class="k">Mean</div>
   <div class="v">{cfg.fmt_money(dist['mean'])}</div></div>
 <div class="stat"><div class="k">P{cfg.RESERVE_PERCENTILE} reserve</div>
   <div class="v">{cfg.fmt_money(dec['reserve_p75'])}</div></div>
</div>

<table>
 <tr><td>P5 &ndash; P95 range</td>
     <td class="n">{cfg.fmt_money(dist['p5'])} &ndash; {cfg.fmt_money(dist['p95'])}</td></tr>
 <tr><td>P25 &ndash; P75 range</td>
     <td class="n">{cfg.fmt_money(dist['p25'])} &ndash; {cfg.fmt_money(dist['p75'])}</td></tr>
 <tr><td>Range relative to central estimate</td><td class="n">{rel}</td></tr>
 <tr><td>Probability of no loss to this portfolio</td>
     <td class="n">{dist['p_zero_loss']:.1%}</td></tr>
</table>

<h2>Reinsurance position</h2>
<table>
 <tr><td>Layer</td>
     <td class="n">{cfg.fmt_money(cfg.LIMIT,0)} xs {cfg.fmt_money(cfg.ATTACHMENT,0)}</td></tr>
 <tr><td>Probability of attachment</td>
     <td class="n">{dec['attachment_probability']:.1%}</td></tr>
 <tr><td>Probability of exhaustion</td>
     <td class="n">{dec['exhaustion_probability']:.1%}</td></tr>
 <tr><td>Expected recovery</td>
     <td class="n">{cfg.fmt_money(dec['mean_cat_xl_recovery'])}</td></tr>
 <tr><td>Members attaching</td>
     <td class="n">{dec['attachment_count']} of {dec['members']}</td></tr>
</table>

{"<h2>Where this uncertainty comes from</h2><table>" + var_rows + "</table>" if var_rows else ""}

<div class="caveat">
 <strong>Basis.</strong> Forecast uncertainty is sampled from IMD's published
 2020&ndash;2024 verification statistics, not from issued advisories. This is a
 hindcast under a published error model. Losses are modelled economic losses on
 a synthetic exposure portfolio, wind peril only: storm surge and rainfall
 flooding are excluded and are material for this basin.
</div>

<footer>
 Project 3, Cyclone Event Response. Landfall reference
 {ev.landfall_lat:.4f}N {ev.landfall_lon:.4f}E, derived by
 {html.escape(ev.landfall_method.replace('_',' '))} on
 {html.escape(ev.coastline_dataset.replace('_',' '))}, within
 {abs(ev.imd_difference_min):.0f} min of the IMD published reference.
 Generated {generated}. Every figure traces to a committed CSV.
</footer>

</div></body></html>"""


def run(event: str | None = None, lead_h: int | None = None,
        all_briefs: bool = False) -> list:
    cfg.ensure_output_dirs()

    dist_path = cfg.STAGE5_DIR / "loss_distribution_summary.csv"
    dec_path = cfg.STAGE6_DIR / "decision_layer_summary.csv"
    var_path = cfg.STAGE7_DIR / "variance_decomposition.csv"

    for p in (dist_path, dec_path):
        if not p.exists():
            raise FileNotFoundError(
                f"Missing {p}. Run Stages 5 and 6 before Stage 9."
            )

    dist_df = pd.read_csv(dist_path)
    dec_df = pd.read_csv(dec_path)
    var_df = pd.read_csv(var_path) if var_path.exists() else None

    print("=" * 72)
    print("STAGE 9 - EVENT RESPONSE BRIEF")
    print("=" * 72)

    targets = []
    if all_briefs:
        for _, r in dec_df.iterrows():
            targets.append((r["event"], int(r["lead_time_h"])))
    else:
        targets.append((event or "Fani", lead_h if lead_h is not None else 48))

    written = []
    for ev, lh in targets:
        htm = build_brief(ev, lh, dist_df, dec_df, var_df)
        path = cfg.STAGE9_DIR / f"event_response_brief_{cfg.event_key(ev)}_T{lh}.html"
        path.write_text(htm, encoding="utf-8")
        written.append(path)
        print(f"  {ev} T-{lh}h -> {path.name}")

    print("\n" + "=" * 72)
    print(f"STAGE 9 COMPLETE ({len(written)} brief(s))")
    print("=" * 72)
    return written


def main():
    ap = argparse.ArgumentParser(description="Project 3 Stage 9.")
    ap.add_argument("--event", choices=cfg.EVENT_NAMES, default="Fani")
    ap.add_argument("--lead", type=int, choices=cfg.LEADS_H, default=48)
    ap.add_argument("--all", action="store_true",
                    help="Write a brief for every event and lead.")
    args = ap.parse_args()
    run(args.event, args.lead, args.all)


if __name__ == "__main__":
    main()
