"""
Project 3 - Stage 7: variance decomposition.

Answers the question the whole project was built to answer: at each point on
the decision clock, where does the uncertainty in the loss estimate actually
come from?

How this works without re-running CLIMADA
-----------------------------------------
Every ensemble member is fully described by three latent draws:

    z_a  along-track  (timing)
    z_c  cross-track  (where it lands along the coast)
    z_v  intensity

and the loss for that member is a deterministic function of them. The draws
depend only on (base_seed, event, lead_time, n_members) through
generate_ensemble.member_seed, NOT on the best track, so they can be
reproduced exactly without loading IBTrACS, CLIMADA or the hazard files. That
makes this stage cheap and independently testable.

The estimator
-------------
Loss is strongly non-linear in position: a member either hits the portfolio or
misses it. A linear regression would badly understate the position terms, so
first-order sensitivity indices are estimated variance-wise rather than by
correlation:

    S_i = Var( E[ loss | z_i ] ) / Var( loss )

estimated from the existing 500 members by binning each z into quantile
groups and measuring between-group variance. The naive correlation ratio
(eta-squared) is biased upward with finite samples, so omega-squared is
reported as the primary figure and eta-squared alongside it. Spearman rank
correlation is reported as a third, assumption-free cross-check.

S_a + S_c + S_v will not sum to 1. The remainder is interaction: members that
are both badly placed and badly sized. That residual is reported rather than
distributed, because attributing it would be a modelling choice presented as
a measurement.

The result to look for
----------------------
sigma_a exceeds sigma_c by a factor of 1.8 to 4.2 at every lead. If S_cross
nonetheless dominates S_along, that independently confirms the section 6.5
decomposition from the loss side, using nothing that went into calibrating
it: along-track error moves the landfall TIME, cross-track error moves the
landfall PLACE, and only the second one moves the loss.
"""

from __future__ import annotations

import argparse

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import config as cfg
from generate_ensemble import member_seed


N_BINS = 10

SOURCE_LABELS = {
    "along": "Track position, along-track (timing)",
    "cross": "Track position, cross-track (landfall place)",
    "intensity": "Intensity",
}


def recover_latent_draws(
    event: str,
    lead_h: int,
    n_members: int = cfg.N_MEMBERS,
    base_seed: int = cfg.BASE_SEED,
) -> dict[str, np.ndarray]:
    """Reproduce a Stage 2 ensemble's latent draws exactly.

    This mirrors generate_ensemble.generate_ensemble's draw order precisely.
    If that order ever changes, this silently desynchronises, so the order is
    asserted against the generator in validate_draw_order().
    """
    rng = np.random.default_rng(member_seed(base_seed, event, lead_h))
    z_a = rng.normal(size=n_members)
    z_c = rng.normal(size=n_members)
    z_v = rng.normal(size=n_members)
    return {"along": z_a, "cross": z_c, "intensity": z_v}


def validate_draw_order() -> None:
    """Guard against the generator's draw order drifting away from this stage.

    Stage 7 joins losses to latent draws by member index. If the generator
    ever drew z_c before z_a, or added a fourth draw, every index would still
    line up and every number would still look plausible. This check is the
    only thing standing between that and a silently wrong decomposition.
    """
    import inspect
    import generate_ensemble as ge

    src = inspect.getsource(ge.generate_ensemble)
    try:
        i_a = src.index("z_a = rng.normal")
        i_c = src.index("z_c = rng.normal")
        i_v = src.index("z_v = rng.normal")
    except ValueError as exc:
        raise AssertionError(
            "generate_ensemble no longer draws z_a, z_c, z_v in the expected "
            "form. Stage 7 reproduces those draws by seed and would "
            "desynchronise silently. Re-check the draw order."
        ) from exc

    if not (i_a < i_c < i_v):
        raise AssertionError(
            "generate_ensemble draws the latent variables in a different "
            "order than Stage 7 reproduces them."
        )


def _binned_indices(z: np.ndarray, n_bins: int) -> np.ndarray:
    """Quantile bins, so each group holds roughly the same number of members."""
    edges = np.quantile(z, np.linspace(0, 1, n_bins + 1))
    edges[0] -= 1e-12
    edges[-1] += 1e-12
    return np.clip(np.digitize(z, edges[1:-1]), 0, n_bins - 1)


def first_order_index(y: np.ndarray, z: np.ndarray, n_bins: int = N_BINS) -> dict:
    """Estimate Var(E[Y|Z]) / Var(Y) from existing samples.

    Returns eta-squared (naive, biased up), omega-squared (bias-corrected,
    primary) and the Spearman rank correlation squared (assumption-free
    cross-check that catches monotone relationships only).
    """
    n = len(y)
    groups = _binned_indices(z, n_bins)
    grand = float(np.mean(y))
    ss_total = float(np.sum((y - grand) ** 2))

    if ss_total <= 0:
        return {"eta_sq": 0.0, "omega_sq": 0.0, "spearman_sq": 0.0,
                "n_bins": n_bins, "degenerate": True}

    ss_between = 0.0
    k_used = 0
    for g in np.unique(groups):
        mask = groups == g
        if mask.sum() == 0:
            continue
        k_used += 1
        ss_between += mask.sum() * (float(np.mean(y[mask])) - grand) ** 2

    ss_within = ss_total - ss_between
    df_within = max(n - k_used, 1)
    ms_within = ss_within / df_within

    eta_sq = ss_between / ss_total
    omega_sq = (ss_between - (k_used - 1) * ms_within) / (ss_total + ms_within)
    omega_sq = float(max(omega_sq, 0.0))

    # Spearman: rank-based, so immune to the shape of the loss response.
    rank_y = pd.Series(y).rank().to_numpy()
    rank_z = pd.Series(z).rank().to_numpy()
    if np.std(rank_y) == 0 or np.std(rank_z) == 0:
        spearman = 0.0
    else:
        spearman = float(np.corrcoef(rank_y, rank_z)[0, 1])

    return {
        "eta_sq": float(eta_sq),
        "omega_sq": omega_sq,
        "spearman_sq": float(spearman ** 2),
        "n_bins": k_used,
        "degenerate": False,
    }


def decompose(event: str, lead_h: int) -> list[dict]:
    # load_losses returns losses ordered by Stage 2 member index, so position
    # i here is member i, which is the same member recover_latent_draws
    # reproduces at position i.
    losses = cfg.load_losses(event, lead_h)
    draws = recover_latent_draws(event, lead_h, n_members=len(losses))

    loss_var = float(np.var(losses, ddof=1))
    degenerate = loss_var <= 0.0

    rows = []
    for source, z in draws.items():
        idx = first_order_index(losses, z)
        idx["degenerate"] = idx["degenerate"] or degenerate
        rows.append({
            "event": event,
            "lead_time_h": lead_h,
            "source": source,
            "label": SOURCE_LABELS[source],
            "loss_variance": float(np.var(losses, ddof=1)),
            "loss_mean": float(np.mean(losses)),
            **idx,
        })
    return rows


def add_residuals(df: pd.DataFrame) -> pd.DataFrame:
    """Add the unexplained share for each (event, lead).

    This is interaction between sources plus higher-order structure. It is
    reported, never distributed across the named sources.
    """
    out = []
    for (event, lead), sub in df.groupby(["event", "lead_time_h"]):
        explained = float(sub["omega_sq"].sum())
        # A loss distribution with no variance has nothing to attribute.
        # Reporting residual = 1.0 there would read as "100% unexplained"
        # when the truth is "every member produced the same loss".
        degenerate = bool(sub["degenerate"].any())
        out.append({
            "event": event,
            "lead_time_h": lead,
            "source": "interaction_residual",
            "label": "Interaction and higher order (not attributed)",
            "loss_variance": float(sub["loss_variance"].iloc[0]),
            "loss_mean": float(sub["loss_mean"].iloc[0]),
            "eta_sq": np.nan,
            "omega_sq": np.nan if degenerate else float(max(1.0 - explained, 0.0)),
            "spearman_sq": np.nan,
            "n_bins": np.nan,
            "degenerate": degenerate,
        })
    return pd.concat([df, pd.DataFrame(out)], ignore_index=True)


def discover_exposure_variants(event: str, lead_h: int) -> dict[str, "Path"]:
    """Find every Project 2 exposure scenario Stage 4 produced for this panel.

    Scenario names come from the Project 2 linkage file and are not known in
    advance (COMBINED_SENSITIVITY, P003_SENSITIVITY and so on), so they are
    discovered by globbing rather than hardcoded. A hardcoded list silently
    finds nothing the moment Project 2 renames a scenario, and the comparison
    then reports "not available" instead of "the names do not match".
    """
    from pathlib import Path

    found: dict[str, Path] = {}

    primary = cfg.loss_file(event, lead_h)
    if primary.exists():
        found["PRIMARY"] = primary

    stem = f"{cfg.event_key(event)}_T{lead_h}_"
    for p in sorted(cfg.IMPACT_DIR.glob(f"{stem}*_event_losses.csv")):
        variant = p.name[len(stem):-len("_event_losses.csv")]
        if variant:
            found[variant.upper()] = p

    return found


def exposure_comparison(event: str, lead_h: int,
                        variants: tuple[str, ...] | None = None) -> dict | None:
    """Compare meteorological spread against exposure-scenario spread.

    Meteorological spread is the variance across the 500 ensemble members for
    one fixed exposure scenario. Exposure spread is how far the MEAN loss moves
    when the exposure scenario changes with the hazard held fixed.

    Caveat worth stating in the write-up: the scenarios are not equiprobable.
    PRIMARY is the treated portfolio and DIRTY is the untreated one; they are
    not two draws from a distribution. So this is a spread across defensible
    alternatives, not a probability-weighted variance, and the range is
    reported alongside the variance for that reason.
    """
    found = discover_exposure_variants(event, lead_h)
    if len(found) < 2:
        return None

    means = {
        name: float(np.mean(cfg.load_losses(event, lead_h, path=p)))
        for name, p in found.items()
    }

    base = cfg.load_losses(event, lead_h)
    met_var = float(np.var(base, ddof=1))
    values = np.array(list(means.values()), dtype=float)
    exposure_var = float(np.var(values, ddof=1))

    return {
        "event": event,
        "lead_time_h": lead_h,
        "meteorological_sd": float(np.sqrt(met_var)),
        "exposure_scenario_sd": float(np.sqrt(exposure_var)),
        "exposure_scenario_range": float(values.max() - values.min()),
        "exposure_share_of_total_variance": (
            exposure_var / (met_var + exposure_var)
            if (met_var + exposure_var) > 0 else np.nan
        ),
        "n_scenarios": len(found),
        "scenarios": ",".join(sorted(found)),
        **{f"mean_loss_{k}": v for k, v in sorted(means.items())},
    }


def stacked_chart(df: pd.DataFrame, event: str):
    sub = df[df["event"] == event]
    leads = sorted(sub["lead_time_h"].unique(), reverse=True)
    order = ["cross", "along", "intensity", "interaction_residual"]

    shares = {
        s: [float(sub[(sub["lead_time_h"] == l) & (sub["source"] == s)]["omega_sq"].sum())
            for l in leads]
        for s in order
    }

    fig, ax = plt.subplots(figsize=(9, 5.5))
    x = np.arange(len(leads))
    bottom = np.zeros(len(leads))
    labels = {
        "cross": "Cross-track (where it lands)",
        "along": "Along-track (timing)",
        "intensity": "Intensity",
        "interaction_residual": "Interaction / unattributed",
    }
    for s in order:
        vals = np.asarray(shares[s])
        ax.bar(x, vals, bottom=bottom, label=labels[s])
        bottom += vals

    ax.set_xticks(x)
    ax.set_xticklabels([f"T-{int(v)}h" for v in leads])
    ax.set_ylim(0, 1)
    ax.set_xlabel("Forecast lead time")
    ax.set_ylabel("Share of loss variance")
    ax.set_title(f"{event}: where the loss uncertainty comes from")
    ax.legend(loc="upper right", fontsize=9)
    ax.grid(True, axis="y", alpha=0.25)
    fig.tight_layout()

    path = cfg.STAGE7_DIR / f"{cfg.event_key(event)}_variance_decomposition.png"
    fig.savefig(path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    return path


def report_anisotropy_check(df: pd.DataFrame) -> pd.DataFrame:
    """The independent confirmation of section 6.5, from the loss side.

    sigma_a > sigma_c at every lead. If cross-track nonetheless explains more
    loss variance than along-track, the along/cross interpretation is doing
    real work rather than fitting the data it was calibrated on.
    """
    import error_model as em

    rows = []
    for (event, lead), sub in df.groupby(["event", "lead_time_h"]):
        if bool(sub["degenerate"].any()):
            continue
        s_cross = float(sub[sub["source"] == "cross"]["omega_sq"].sum())
        s_along = float(sub[sub["source"] == "along"]["omega_sq"].sum())
        try:
            sa, sc = em.landfall_sampler_parameters(int(lead))
            sigma_ratio = sa / sc
        except ValueError:
            sigma_ratio = np.nan

        rows.append({
            "event": event,
            "lead_time_h": lead,
            "sigma_a_over_sigma_c": sigma_ratio,
            "S_cross": s_cross,
            "S_along": s_along,
            "variance_ratio_cross_over_along": (
                s_cross / s_along if s_along > 0 else np.inf
            ),
            "confirms_decomposition": bool(s_cross > s_along),
        })
    return pd.DataFrame(rows).sort_values(
        ["event", "lead_time_h"], ascending=[True, False]
    ).reset_index(drop=True)


def run(events: tuple[str, ...] | None = None,
        variants: tuple[str, ...] | None = None) -> pd.DataFrame:
    cfg.ensure_output_dirs()
    validate_draw_order()

    events = events or cfg.available_events()
    if not events:
        raise RuntimeError("No Stage 4 loss files found. Run Stage 4 first.")

    print("=" * 72)
    print("STAGE 7 - VARIANCE DECOMPOSITION")
    print("=" * 72)
    print("\nLatent draw order verified against the Stage 2 generator.")
    print(f"Estimator: first-order variance index, {N_BINS} quantile bins,")
    print("omega-squared (bias corrected) reported as primary.")

    rows = []
    for event in events:
        leads = cfg.available_leads(event)
        if not leads:
            continue
        print(f"\n[{event}]")
        print(f"  {'lead':>6} {'cross':>8} {'along':>8} {'intensity':>10} "
              f"{'residual':>9}")
        for lead in leads:
            r = decompose(event, lead)
            rows.extend(r)
            if any(x["degenerate"] for x in r):
                print(f"  {'T-' + str(lead):>6} {'-':>8} {'-':>8} {'-':>10} "
                      f"{'-':>9}   no loss variance to decompose")
                continue
            d = {x["source"]: x["omega_sq"] for x in r}
            resid = max(1.0 - sum(d.values()), 0.0)
            print(f"  {'T-' + str(lead):>6} {d['cross']:8.3f} {d['along']:8.3f} "
                  f"{d['intensity']:10.3f} {resid:9.3f}")

    df = add_residuals(pd.DataFrame(rows))

    aniso = report_anisotropy_check(df)
    print("\n[ANISOTROPY CROSS-CHECK]")
    print("  sigma_a exceeds sigma_c everywhere. Does cross-track still")
    print("  dominate the loss variance?")
    print(f"  {'event':<9} {'lead':>6} {'sig_a/sig_c':>12} {'S_cross':>9} "
          f"{'S_along':>9} {'confirms':>9}")
    for _, r in aniso.iterrows():
        ratio = (f"{r['sigma_a_over_sigma_c']:12.2f}"
                 if np.isfinite(r["sigma_a_over_sigma_c"]) else f"{'n/a':>12}")
        print(f"  {r['event']:<9} {'T-' + str(int(r['lead_time_h'])):>6} {ratio} "
              f"{r['S_cross']:9.3f} {r['S_along']:9.3f} "
              f"{'YES' if r['confirms_decomposition'] else 'no':>9}")

    n_conf = int(aniso["confirms_decomposition"].sum())
    print(f"\n  Cross-track dominates along-track in {n_conf}/{len(aniso)} cases.")

    # Optional exposure comparison.
    exposure_rows = []
    for event in events:
        for lead in cfg.available_leads(event):
            res = exposure_comparison(event, lead, variants)
            if res:
                exposure_rows.append(res)

    if exposure_rows:
        ex_df = pd.DataFrame(exposure_rows)
        ex_path = cfg.STAGE7_DIR / "exposure_vs_meteorological.csv"
        ex_df.to_csv(ex_path, index=False)
        print("\n[EXPOSURE SCENARIO COMPARISON]")
        print(ex_df.to_string(index=False))
    else:
        print("\n[EXPOSURE SCENARIO COMPARISON]")
        print("  Skipped: fewer than two exposure scenarios found in")
        print(f"  {cfg.IMPACT_DIR}")
        print("  Run: python compute_impact.py --all-scenarios")

    summary_path = cfg.STAGE7_DIR / "variance_decomposition.csv"
    df.to_csv(summary_path, index=False)

    aniso_path = cfg.STAGE7_DIR / "anisotropy_cross_check.csv"
    aniso.to_csv(aniso_path, index=False)

    charts = [stacked_chart(df, e) for e in df["event"].unique()]

    print("\n[OUTPUTS]")
    print(f"  {summary_path}")
    print(f"  {aniso_path}")
    for c in charts:
        print(f"  {c}")

    print("\n" + "=" * 72)
    print("STAGE 7 COMPLETE")
    print("=" * 72)
    return df


def main():
    ap = argparse.ArgumentParser(description="Project 3 Stage 7.")
    ap.add_argument("--event", choices=cfg.EVENT_NAMES, action="append")
    args = ap.parse_args()
    run(tuple(args.event) if args.event else None)


if __name__ == "__main__":
    main()
