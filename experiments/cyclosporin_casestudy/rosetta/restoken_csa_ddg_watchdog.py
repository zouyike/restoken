#!/usr/bin/env python3
"""Auto-generated watchdog: restoken_csa_ddg"""
import json
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, "/public/home/genesis/agent_exchange/shared_scripts")
from cron_watchdog import run_watchdog, WatchdogResult

OUTPUT_DIR = Path("/scratch/genesis/NCAA_tokenization/experiments/cyclosporin_casestudy/rosetta")
DONE_FILE = OUTPUT_DIR / "DONE"
LOG_FILE = OUTPUT_DIR / "batch_r3.log"


def check_fn():
    # Check if process is alive
    try:
        pid = 1023563
        import os
        os.kill(pid, 0)
        running = True
    except (OSError, ProcessLookupError):
        running = False

    # Check done file
    if DONE_FILE.exists():
        return WatchdogResult(
            report="restoken_csa_ddg COMPLETE. Check output in /scratch/genesis/NCAA_tokenization/experiments/cyclosporin_casestudy/rosetta",
            is_complete=True, is_error=False,
        )

    # Check log tail
    log_tail = ""
    if LOG_FILE.exists():
        lines = LOG_FILE.read_text().strip().split("\n")
        log_tail = "\n".join(lines[-3:])

    if not running and not DONE_FILE.exists():
        return WatchdogResult(
            report=f"restoken_csa_ddg STOPPED (PID 987136 gone, no DONE file).\nLast log: {log_tail}",
            is_complete=True, is_error=True,
        )

    return WatchdogResult(
        report=f"restoken_csa_ddg running (PID 987136).\nLast log: {log_tail}",
        is_complete=False, is_error=False,
    )


if __name__ == "__main__":
    run_watchdog("restoken_csa_ddg", check_fn, max_runs=48)
