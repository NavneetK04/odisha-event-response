"""
Project 3 - pipeline runner.

Runs the whole chain in order and stops at the first stage that fails, rather
than carrying a broken intermediate forward into a stage that will happily
produce confident output from it.

    Stage 1  error model validation          (error_model.py)
    Stage 3  hazard, per event and lead      (compute_hazard.py, needs CLIMADA)
    Stage 4  impact, per event and lead      (compute_impact.py, needs CLIMADA)
    Stage 5  loss distribution
    Stage 6  decision layer
    Stage 7  variance decomposition
    Stage 8  trigger backtest
    Stage 9  event response briefs

Stages 3 and 4 are run as subprocesses, one event and lead at a time, because
they are the expensive ones and a crash halfway through should not lose the
work already done. Everything from Stage 5 onward runs in-process.

Usage
-----
    python run_pipeline.py --from 5              # post-processing only
    python run_pipeline.py --all                 # everything, all events
    python run_pipeline.py --all --event Fani    # everything, one event
    python run_pipeline.py --from 5 --skip-missing
"""

from __future__ import annotations

import argparse
import subprocess
import sys
import time
from pathlib import Path

import config as cfg


HERE = Path(__file__).resolve().parent


def banner(text: str) -> None:
    print("\n" + "#" * 72)
    print(f"# {text}")
    print("#" * 72)


def run_script(script: str, args: list[str]) -> bool:
    cmd = [sys.executable, str(HERE / script), *args]
    print(f"\n$ {' '.join(cmd[1:])}")
    t0 = time.perf_counter()
    result = subprocess.run(cmd, cwd=HERE)
    dt = time.perf_counter() - t0
    ok = result.returncode == 0
    print(f"  -> {'ok' if ok else 'FAILED'} in {dt:.1f}s")
    return ok


def stage1() -> bool:
    banner("STAGE 1 - forecast error model validation")
    return run_script("error_model.py", [])


def stage3_4(events: tuple[str, ...], leads: tuple[int, ...],
             skip_missing: bool) -> bool:
    ok = True
    for stage, script in ((3, "compute_hazard.py"), (4, "compute_impact.py")):
        banner(f"STAGE {stage} - {script}")
        if not (HERE / script).exists():
            msg = f"  {script} not found in {HERE}"
            if skip_missing:
                print(msg + " (skipped)")
                continue
            print(msg)
            return False
        for event in events:
            for lead in leads:
                # Stage 4 must produce every exposure scenario in the same
                # invocation. Running PRIMARY now and the variants later leaves
                # Stage 7 comparing loss files generated from different
                # ensembles, and it has no way to notice.
                extra = ["--all-scenarios"] if stage == 4 else []
                args_ = ["--event", event, "--lead", str(lead), *extra]
                if not run_script(script, args_):
                    ok = False
                    if not skip_missing:
                        return False
    return ok


def main() -> int:
    ap = argparse.ArgumentParser(description="Project 3 pipeline runner.")
    ap.add_argument("--from", dest="from_stage", type=int, default=5,
                    choices=(1, 3, 5, 6, 7, 8, 9),
                    help="First stage to run. Default 5 (post-processing).")
    ap.add_argument("--all", action="store_true",
                    help="Run from Stage 1, including hazard and impact.")
    ap.add_argument("--event", choices=cfg.EVENT_NAMES, action="append")
    ap.add_argument("--lead", type=int, choices=cfg.LEADS_H, action="append")
    ap.add_argument("--skip-missing", action="store_true",
                    help="Continue past missing scripts or failed event/lead runs.")
    ap.add_argument("--brief-event", choices=cfg.EVENT_NAMES, default="Fani")
    ap.add_argument("--brief-lead", type=int, choices=cfg.LEADS_H, default=48)
    args = ap.parse_args()

    cfg.ensure_output_dirs()

    events = tuple(args.event) if args.event else cfg.EVENT_NAMES
    leads = tuple(args.lead) if args.lead else cfg.LEADS_H
    start = 1 if args.all else args.from_stage

    print("=" * 72)
    print("PROJECT 3 - CYCLONE EVENT RESPONSE PIPELINE")
    print("=" * 72)
    print(f"Events      : {', '.join(events)}")
    print(f"Leads       : {', '.join('T-' + str(l) for l in leads)}")
    print(f"From stage  : {start}")
    print(f"Members     : {cfg.N_MEMBERS}   Base seed: {cfg.BASE_SEED}")

    t0 = time.perf_counter()

    if start <= 1 and not stage1():
        print("\nStage 1 failed. The error model is the spine of the project; "
              "nothing downstream is meaningful until it passes.")
        return 1

    if start <= 3 and not stage3_4(events, leads, args.skip_missing):
        print("\nStage 3/4 failed. Stopping before post-processing.")
        return 1

    ev = tuple(e for e in events if cfg.available_leads(e))
    if not ev:
        print("\nNo Stage 4 loss files found for the requested events.")
        print(f"Expected files like: {cfg.loss_file(events[0], leads[0])}")
        return 1

    if start <= 5:
        banner("STAGE 5 - loss distribution")
        import stage5_loss_distribution as s5
        s5.run(ev)

    if start <= 6:
        banner("STAGE 6 - decision layer")
        import stage6_decision_layer as s6
        s6.run(ev)

    if start <= 7:
        banner("STAGE 7 - variance decomposition")
        import stage7_variance_decomposition as s7
        s7.run(ev)

    if start <= 8:
        banner("STAGE 8 - trigger backtest")
        import stage8_trigger_backtest as s8
        s8.run(ev)

    if start <= 9:
        banner("STAGE 9 - event response briefs")
        import stage9_event_brief as s9
        s9.run(all_briefs=True)

    dt = time.perf_counter() - t0
    print("\n" + "=" * 72)
    print(f"PIPELINE COMPLETE in {dt:.1f}s")
    print("=" * 72)
    print(f"\nOutputs under {cfg.OUTPUT_DIR}")
    for d in cfg.ALL_OUTPUT_DIRS:
        files = sorted(d.glob("*"))
        if files:
            print(f"\n  {d.name}/")
            for f in files:
                print(f"    {f.name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
