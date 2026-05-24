#!/usr/bin/env python3
"""Re-run Experiment 2 (EDIT/FROZEN controllability) for API models.

Reason: previous results invalid due to 3-digit ID truncation parser bug,
now fixed in run_benchmark.py (commit f586e12).

Models: gemini-2.5-pro, gemini-2.5-flash, gpt-4o
"""

import os
import sys
import time
from pathlib import Path

# Setup paths
PROJECT_ROOT = Path("/scratch/genesis/NCAA_tokenization")
sys.path.insert(0, str(PROJECT_ROOT))
sys.path.insert(0, str(PROJECT_ROOT / "experiments" / "script"))

# Proxy for API calls
os.environ.setdefault("HTTPS_PROXY", "http://127.0.0.1:48890")
os.environ.setdefault("HTTP_PROXY", "http://127.0.0.1:48890")

from restoken.src.library import BlockLibrary
from run_benchmark import GeminiBackend, OpenAIBackend, run_exp2

N_PARENTS = 10
N_VARIANTS = 20
LENGTH = 6

MODELS = [
    ("gemini", "gemini-2.5-pro"),
    ("gemini", "gemini-2.5-flash"),
    ("openai", "gpt-4o"),
]


def main():
    print(f"=== Exp2 Re-run (parser fix) — {time.strftime('%Y-%m-%d %H:%M:%S')} ===")
    print(f"  n_parents={N_PARENTS}, n_variants={N_VARIANTS}, length={LENGTH}")
    print()

    lib = BlockLibrary()
    print(f"Library loaded: {lib.size} blocks\n")

    for backend_type, model_name in MODELS:
        print(f"\n{'='*60}")
        print(f"  Model: {model_name} (backend: {backend_type})")
        print(f"{'='*60}")

        if backend_type == "gemini":
            backend = GeminiBackend(model_name)
        elif backend_type == "openai":
            backend = OpenAIBackend(model_name)

        try:
            summary = run_exp2(backend, lib, N_PARENTS, N_VARIANTS, LENGTH)
            print(f"\n  [OK] {model_name} complete.")
        except Exception as e:
            print(f"\n  [FAIL] {model_name}: {e}")
            import traceback
            traceback.print_exc()

        # Brief pause between models
        time.sleep(5)

    print(f"\n=== All exp2 re-runs complete — {time.strftime('%Y-%m-%d %H:%M:%S')} ===")


if __name__ == "__main__":
    main()
