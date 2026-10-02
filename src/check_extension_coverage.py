"""
Project 3 - did the post-landfall extension actually cover every member?

This is the decisive test for the truncation defect, and it is the one that
check_zero_mechanism.py cannot perform. That script tests whether the zero-loss
set is one-sided in z_a. One-sidedness is necessary for truncation but it is NOT
sufficient, because along-track displacement is a real displacement: if a
storm's direction of motion points towards the portfolio, then forward and
backward shifts have genuinely different loss consequences and the zeros will be
one-sided for reasons that have nothing to do with array length.

Phailin is the case in this set. It made landfall moving northwest, up the coast
towards the exposure block, so a forward shift carries it closer and a backward
shift carries it out to sea. Its zeros are one-sided by geometry.

The test that is not confounded compares, per panel, the largest BACKWARD
along-track displacement actually drawn against the path length of the
extension that was provided. If the extension is longer, no member ran off the
end of its track, and every zero is a statement about where the storm went.

Needs outputs/track_extension.csv, written by Stage 3.

    python src/check_extension_coverage.py
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
import error_model as em
from stage7_variance_decomposition import recover_latent_draws


def main() -> int:
    ext_path = cfg.OUTPUT_DIR / "track_extension.csv"
    if not ext_path.exists():
        print(f"Missing {ext_path}. Run Stage 3 first.")
        return 1

    ext = pd.read_csv(ext_path)
    ext_km = {
        (r.event, int(r.lead_time_h)): float(r.available_extension_km)
        for r in ext.itertuples()
    }

    rows = []
    for event in cfg.EVENT_NAMES:
        for lead in cfg.LEADS_H:
            if lead == 0:
                continue  # one point, no along-track displacement possible
            loss_path = cfg.loss_file(event, lead)
            if not loss_path.exists():
                continue
            if (event, lead) not in ext_km:
                print(f"  {event} T-{lead}h has no extension record, skipped.")
                continue

            df = pd.read_csv(loss_path)
            if "member_index" in df.columns:
                df = df.sort_values("member_index")
            losses = df["gross_loss"].to_numpy(dtype=float)

            draws = recover_latent_draws(event, lead, n_members=len(losses))
            z_a = np.asarray(draws["along"], dtype=float)

            sigma_a_lf, _ = em.landfall_sampler_parameters(lead)
            along_km = z_a * float(sigma_a_lf)

            available = ext_km[(event, lead)]
            worst_backward = float(max(0.0, -along_km.min()))
            beyond = along_km < -available
            n_beyond = int(beyond.sum())
            n_beyond_zero = int((beyond & (losses == 0.0)).sum())

            rows.append({
                "event": event,
                "lead_time_h": lead,
                "sigma_a_landfall_km": float(sigma_a_lf),
                "extension_km": available,
                "extension_sigmas": available / float(sigma_a_lf),
                "worst_backward_km": worst_backward,
                "worst_backward_sigmas": worst_backward / float(sigma_a_lf),
                "members_beyond_extension": n_beyond,
                "of_which_zero_loss": n_beyond_zero,
                "n_zero_total": int((losses == 0.0).sum()),
                "covered": bool(n_beyond == 0),
            })

    if not rows:
        print("No panels to check.")
        return 1

    out = pd.DataFrame(rows)
    pd.set_option("display.width", 220)

    print("=" * 108)
    print("PROJECT 3 - POST-LANDFALL EXTENSION COVERAGE")
    print("=" * 108)
    print("\nA panel is covered when no member's backward along-track")
    print("displacement exceeds the extension provided. Covered means every")
    print("zero loss in that panel is a miss, not a truncated track.\n")

    print(out[[
        "event", "lead_time_h", "extension_km", "extension_sigmas",
        "worst_backward_km", "worst_backward_sigmas",
        "members_beyond_extension", "of_which_zero_loss", "n_zero_total",
        "covered",
    ]].to_string(index=False, float_format=lambda x: f"{x:.2f}"))

    path = cfg.OUTPUT_DIR / "extension_coverage.csv"
    out.to_csv(path, index=False)

    bad = out[~out["covered"]]
    print("\n" + "=" * 108)
    if bad.empty:
        print("Every panel is covered. No member ran off the end of its track,")
        print("so every zero loss in the project is a miss. Any remaining")
        print("one-sidedness in z_a among the zeros is track geometry: the")
        print("direction of motion relative to the exposure block.")
    else:
        total = int(bad["members_beyond_extension"].sum())
        print(f"{len(bad)} panel(s) left {total} member(s) beyond the extension.")
        print("Those members' zero losses are artifacts. If the count is small")
        print("relative to N it is a documented residual rather than a defect,")
        print("but the count belongs in the limitations either way:")
        for _, r in bad.iterrows():
            print(
                f"  {r['event']} T-{int(r['lead_time_h'])}h: "
                f"{int(r['members_beyond_extension'])} of {cfg.N_MEMBERS} "
                f"({int(r['members_beyond_extension']) / cfg.N_MEMBERS:.2%}), "
                f"extension {r['extension_km']:.0f} km = "
                f"{r['extension_sigmas']:.2f} sigma, worst draw "
                f"{r['worst_backward_km']:.0f} km"
            )
    print(f"\nWritten to {path}")
    print("=" * 108)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
