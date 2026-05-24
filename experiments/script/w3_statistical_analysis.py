#!/usr/bin/env python3
"""
W3 Statistical Analysis for ResToken JCIM Paper.

Computes:
  1. Bootstrap 95% CIs for exp1 validity rates (per model x representation)
  2. Chi-square tests: ResToken vs SMILES, ResToken vs HELM validity (per model)
  3. Aggregate cross-model statistics for exp1, exp2, exp3
  4. Effect-size enrichment of ResToken constraint satisfaction over random baseline
  5. JSON summary for manuscript use

Usage:
    python w3_statistical_analysis.py
"""

import json
import os
import sys
from pathlib import Path
from collections import defaultdict

import numpy as np
from scipy import stats

# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------
BASE = Path("/scratch/genesis/NCAA_tokenization/experiments/output")
OUT_JSON = BASE / "w3_statistical_summary.json"

MODELS = [
    "gemini-2.5-pro",
    "gemini-2.5-flash",
    "gpt-4o",
    "Qwen3.5-9B",
    "gemma-3-12b-it-bnb-4bit",
]
REPRESENTATIONS = ["restoken", "smiles", "helm"]
PROFILES = ["permeable", "charged_binder", "rigid_scaffold"]

# Random baselines for enrichment analysis (exp3)
RANDOM_BASELINES = {
    "permeable": 0.0207,
    "charged_binder": 0.0182,
    "rigid_scaffold": 0.0081,
}

N_BOOTSTRAP = 10_000
SEED = 42
ALPHA = 0.05  # for 95% CI

# Display names (shorter)
MODEL_SHORT = {
    "gemini-2.5-pro": "Gemini-Pro",
    "gemini-2.5-flash": "Gemini-Flash",
    "gpt-4o": "GPT-4o",
    "Qwen3.5-9B": "Qwen3.5-9B",
    "gemma-3-12b-it-bnb-4bit": "Gemma-12B",
}


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def load_json(path):
    """Load JSON, return None if missing or empty."""
    if not path.exists():
        return None
    with open(path) as f:
        data = json.load(f)
    return data


def get_exp1_results(model, rep):
    """Load per-sequence results for exp1 validity."""
    fname = f"{model}_{rep}_unconstrained_len6_results.json"
    return load_json(BASE / "exp1_validity" / fname)


def get_exp1_summary(model, rep):
    fname = f"{model}_{rep}_unconstrained_len6_summary.json"
    return load_json(BASE / "exp1_validity" / fname)


def get_exp2_summary(model):
    fname = f"{model}_exp2_controllability_summary.json"
    return load_json(BASE / "exp2_controllability" / fname)


def get_exp2_results(model):
    fname = f"{model}_exp2_controllability_results.json"
    return load_json(BASE / "exp2_controllability" / fname)


def get_exp3_summary(model, rep, profile):
    fname = f"{model}_{rep}_{profile}_len6_summary.json"
    return load_json(BASE / "exp3_property_constrained" / fname)


def get_exp3_results(model, rep, profile):
    fname = f"{model}_{rep}_{profile}_len6_results.json"
    return load_json(BASE / "exp3_property_constrained" / fname)


def bootstrap_ci(binary_array, n_boot=N_BOOTSTRAP, alpha=ALPHA, rng=None):
    """Bootstrap 95% CI for a proportion (binary array)."""
    if rng is None:
        rng = np.random.default_rng(SEED)
    arr = np.asarray(binary_array, dtype=float)
    n = len(arr)
    if n == 0:
        return {"point": 0.0, "ci_lower": 0.0, "ci_upper": 0.0, "n": 0}
    point = float(arr.mean())
    # Generate bootstrap resamples
    boot_indices = rng.integers(0, n, size=(n_boot, n))
    boot_means = arr[boot_indices].mean(axis=1)
    ci_lower = float(np.percentile(boot_means, 100 * alpha / 2))
    ci_upper = float(np.percentile(boot_means, 100 * (1 - alpha / 2)))
    return {"point": point, "ci_lower": ci_lower, "ci_upper": ci_upper, "n": n}


def extract_validity_binary(results, rep):
    """Extract binary validity vector from per-sequence results."""
    if results is None or len(results) == 0:
        return np.array([])
    return np.array([1.0 if r.get("overall_valid", False) else 0.0 for r in results])


def cramers_v(chi2, n, k=2):
    """Cramer's V for 2x2 table (k = min(rows,cols) = 2)."""
    return float(np.sqrt(chi2 / (n * (k - 1)))) if n > 0 and chi2 > 0 else 0.0


def chi_square_2x2(n1_success, n1_total, n2_success, n2_total):
    """Chi-square test for two proportions (2x2 contingency table).
    Returns chi2, p, cramers_v.
    """
    table = np.array([
        [n1_success, n1_total - n1_success],
        [n2_success, n2_total - n2_success],
    ])
    # Check minimum expected frequency
    row_sums = table.sum(axis=1)
    col_sums = table.sum(axis=0)
    total = table.sum()
    if total == 0 or any(row_sums == 0) or any(col_sums == 0):
        return 0.0, 1.0, 0.0
    # Use chi2_contingency with Yates correction
    chi2, p, dof, expected = stats.chi2_contingency(table, correction=True)
    v = cramers_v(chi2, total)
    return float(chi2), float(p), float(v)


# ---------------------------------------------------------------------------
# Section 1: Bootstrap CIs for exp1 validity
# ---------------------------------------------------------------------------
def compute_bootstrap_cis():
    print("=" * 80)
    print("SECTION 1: Bootstrap 95% Confidence Intervals for Exp1 Validity Rates")
    print("=" * 80)
    rng = np.random.default_rng(SEED)
    results = {}

    # Header
    header = f"{'Model':<20} {'Rep':<10} {'N':>5} {'Rate':>8} {'95% CI Lower':>14} {'95% CI Upper':>14}"
    print(header)
    print("-" * len(header))

    for model in MODELS:
        for rep in REPRESENTATIONS:
            res = get_exp1_results(model, rep)
            binary = extract_validity_binary(res, rep)
            ci = bootstrap_ci(binary, rng=rng)
            key = f"{model}|{rep}"
            results[key] = ci

            mname = MODEL_SHORT.get(model, model)
            print(f"{mname:<20} {rep:<10} {ci['n']:>5} {ci['point']:>8.3f} "
                  f"{ci['ci_lower']:>14.3f} {ci['ci_upper']:>14.3f}")
        print()

    return results


# ---------------------------------------------------------------------------
# Section 2: Chi-square tests (ResToken vs SMILES, ResToken vs HELM)
# ---------------------------------------------------------------------------
def compute_chi_square_tests(bootstrap_results):
    print("\n" + "=" * 80)
    print("SECTION 2: Chi-Square Tests — ResToken vs SMILES / HELM (Exp1 Validity)")
    print("=" * 80)

    n_tests = len(MODELS) * 2  # 2 comparisons per model
    bonferroni_alpha = ALPHA / n_tests
    print(f"Bonferroni-corrected alpha = {ALPHA} / {n_tests} = {bonferroni_alpha:.4f}\n")

    results = {}
    header = (f"{'Model':<20} {'Comparison':<22} {'chi2':>8} {'p-value':>12} "
              f"{'Cramers V':>10} {'Signif':>8}")
    print(header)
    print("-" * len(header))

    for model in MODELS:
        rt_data = get_exp1_results(model, "restoken")
        sm_data = get_exp1_results(model, "smiles")
        hm_data = get_exp1_results(model, "helm")

        rt_bin = extract_validity_binary(rt_data, "restoken")
        sm_bin = extract_validity_binary(sm_data, "smiles")
        hm_bin = extract_validity_binary(hm_data, "helm")

        mname = MODEL_SHORT.get(model, model)

        for comp_name, comp_bin in [("vs SMILES", sm_bin), ("vs HELM", hm_bin)]:
            n1_total = len(rt_bin)
            n1_success = int(rt_bin.sum()) if len(rt_bin) > 0 else 0
            n2_total = len(comp_bin)
            n2_success = int(comp_bin.sum()) if len(comp_bin) > 0 else 0

            chi2, p, v = chi_square_2x2(n1_success, n1_total, n2_success, n2_total)
            sig = "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < bonferroni_alpha else "ns"

            key = f"{model}|restoken {comp_name}"
            results[key] = {
                "chi2": chi2,
                "p_value": p,
                "p_value_bonferroni_threshold": bonferroni_alpha,
                "cramers_v": v,
                "significant": p < bonferroni_alpha,
                "n_restoken": n1_total,
                "n_valid_restoken": n1_success,
                "n_comparator": n2_total,
                "n_valid_comparator": n2_success,
            }

            print(f"{mname:<20} {comp_name:<22} {chi2:>8.2f} {p:>12.2e} "
                  f"{v:>10.3f} {sig:>8}")
        print()

    return results


# ---------------------------------------------------------------------------
# Section 3: Aggregate cross-model statistics
# ---------------------------------------------------------------------------
def compute_aggregate_stats():
    print("\n" + "=" * 80)
    print("SECTION 3: Aggregate Cross-Model Statistics")
    print("=" * 80)

    # --- 3a: Exp1 validity ---
    print("\n--- 3a. Exp1: Mean +/- Std Validity Across Models ---")
    exp1_rates = defaultdict(list)
    for model in MODELS:
        for rep in REPRESENTATIONS:
            summary = get_exp1_summary(model, rep)
            if summary and summary.get("n_parsed", 0) > 0:
                exp1_rates[rep].append(summary["overall_valid"])

    header = f"{'Representation':<15} {'N Models':>10} {'Mean':>8} {'Std':>8} {'Min':>8} {'Max':>8}"
    print(header)
    print("-" * len(header))
    exp1_agg = {}
    for rep in REPRESENTATIONS:
        vals = exp1_rates[rep]
        if vals:
            m, s = np.mean(vals), np.std(vals, ddof=1) if len(vals) > 1 else 0.0
            exp1_agg[rep] = {
                "mean": float(m), "std": float(s),
                "min": float(min(vals)), "max": float(max(vals)),
                "n_models": len(vals), "values": [float(v) for v in vals],
            }
            print(f"{rep:<15} {len(vals):>10} {m:>8.3f} {s:>8.3f} "
                  f"{min(vals):>8.3f} {max(vals):>8.3f}")
        else:
            exp1_agg[rep] = {"mean": 0, "std": 0, "n_models": 0}
            print(f"{rep:<15} {'N/A':>10}")

    # --- 3b: Exp2 controllability ---
    print("\n--- 3b. Exp2: Edit & Frozen Compliance (ResToken only) ---")
    exp2_metrics = defaultdict(list)
    for model in MODELS:
        summary = get_exp2_summary(model)
        if summary and summary.get("n_total_variants", 0) > 0:
            for k in ["edit_compliance", "frozen_compliance", "id_valid", "length_match"]:
                exp2_metrics[k].append(summary[k])

    header = f"{'Metric':<22} {'N Models':>10} {'Mean':>8} {'Std':>8}"
    print(header)
    print("-" * len(header))
    exp2_agg = {}
    for k in ["edit_compliance", "frozen_compliance", "id_valid", "length_match"]:
        vals = exp2_metrics[k]
        if vals:
            m, s = np.mean(vals), np.std(vals, ddof=1) if len(vals) > 1 else 0.0
            exp2_agg[k] = {"mean": float(m), "std": float(s), "n_models": len(vals),
                           "values": [float(v) for v in vals]}
            print(f"{k:<22} {len(vals):>10} {m:>8.3f} {s:>8.3f}")

    # --- 3c: Exp3 constraint satisfaction by profile (ResToken only) ---
    print("\n--- 3c. Exp3: Constraint Satisfaction by Profile (ResToken) ---")
    exp3_agg = {}
    for profile in PROFILES:
        vals = []
        for model in MODELS:
            summary = get_exp3_summary(model, "restoken", profile)
            if summary and "constraint_satisfaction" in summary:
                vals.append(summary["constraint_satisfaction"])
        if vals:
            m, s = np.mean(vals), np.std(vals, ddof=1) if len(vals) > 1 else 0.0
            exp3_agg[profile] = {
                "mean": float(m), "std": float(s), "n_models": len(vals),
                "values": [float(v) for v in vals],
            }
            print(f"  {profile:<22} n={len(vals):>2}  mean={m:.3f} +/- {s:.3f}  "
                  f"range=[{min(vals):.3f}, {max(vals):.3f}]")
        else:
            exp3_agg[profile] = {"mean": 0, "std": 0, "n_models": 0}

    # --- 3d: Exp3 overall validity by representation (cross-model) ---
    print("\n--- 3d. Exp3: Overall Validity by Representation x Profile ---")
    exp3_validity = {}
    for profile in PROFILES:
        exp3_validity[profile] = {}
        for rep in REPRESENTATIONS:
            vals = []
            for model in MODELS:
                summary = get_exp3_summary(model, rep, profile)
                if summary and summary.get("n_parsed", 0) > 0:
                    vals.append(summary.get("overall_valid", 0.0))
            if vals:
                m = np.mean(vals)
                s = np.std(vals, ddof=1) if len(vals) > 1 else 0.0
                exp3_validity[profile][rep] = {
                    "mean": float(m), "std": float(s), "n_models": len(vals),
                    "values": [float(v) for v in vals],
                }
                print(f"  {profile:<20} {rep:<10} n={len(vals):>2}  "
                      f"mean={m:.3f} +/- {s:.3f}")

    return {
        "exp1_validity_by_rep": exp1_agg,
        "exp2_controllability": exp2_agg,
        "exp3_constraint_satisfaction_restoken": exp3_agg,
        "exp3_validity_by_rep_profile": exp3_validity,
    }


# ---------------------------------------------------------------------------
# Section 4: Enrichment over random baseline
# ---------------------------------------------------------------------------
def compute_enrichment():
    print("\n" + "=" * 80)
    print("SECTION 4: ResToken Constraint Satisfaction Enrichment over Random Baseline")
    print("=" * 80)

    results = {}
    header = (f"{'Model':<20} {'Profile':<18} {'Observed':>10} {'Random':>10} "
              f"{'Enrichment':>12} {'N valid':>8}")
    print(header)
    print("-" * len(header))

    for model in MODELS:
        for profile in PROFILES:
            summary = get_exp3_summary(model, "restoken", profile)
            if summary and "constraint_satisfaction" in summary:
                obs = summary["constraint_satisfaction"]
                baseline = RANDOM_BASELINES[profile]
                enrichment = obs / baseline if baseline > 0 else float("inf") if obs > 0 else 0.0
                n_valid = summary.get("n_valid", 0)

                key = f"{model}|{profile}"
                results[key] = {
                    "observed": obs,
                    "random_baseline": baseline,
                    "enrichment_fold": float(enrichment),
                    "n_valid": n_valid,
                    "n_constraint_pass": summary.get("n_constraint_pass", 0),
                }

                mname = MODEL_SHORT.get(model, model)
                print(f"{mname:<20} {profile:<18} {obs:>10.3f} {baseline:>10.4f} "
                      f"{enrichment:>11.1f}x {n_valid:>8}")

    # Aggregate enrichment per profile
    print("\n--- Aggregate Enrichment per Profile ---")
    agg = {}
    for profile in PROFILES:
        enrichments = []
        for model in MODELS:
            key = f"{model}|{profile}"
            if key in results and results[key]["observed"] > 0:
                enrichments.append(results[key]["enrichment_fold"])
        if enrichments:
            m = np.mean(enrichments)
            agg[profile] = {
                "mean_enrichment": float(m),
                "n_models_with_nonzero": len(enrichments),
            }
            print(f"  {profile:<20} mean enrichment = {m:.1f}x  "
                  f"(n={len(enrichments)} models with >0 satisfaction)")
        else:
            agg[profile] = {"mean_enrichment": 0.0, "n_models_with_nonzero": 0}
            print(f"  {profile:<20} no models achieved >0 satisfaction")

    results["aggregate"] = agg
    return results


# ---------------------------------------------------------------------------
# Section 5: Detailed exp2 per-model table
# ---------------------------------------------------------------------------
def print_exp2_detail():
    print("\n" + "=" * 80)
    print("DETAIL: Exp2 Controllability Per-Model Results (ResToken)")
    print("=" * 80)
    header = (f"{'Model':<20} {'Edit Compl':>12} {'Frozen Compl':>14} "
              f"{'ID Valid':>10} {'Len Match':>10} {'N Variants':>12}")
    print(header)
    print("-" * len(header))
    for model in MODELS:
        summary = get_exp2_summary(model)
        if summary:
            mname = MODEL_SHORT.get(model, model)
            print(f"{mname:<20} {summary['edit_compliance']:>12.3f} "
                  f"{summary['frozen_compliance']:>14.3f} "
                  f"{summary['id_valid']:>10.3f} "
                  f"{summary['length_match']:>10.3f} "
                  f"{summary['n_total_variants']:>12}")


# ---------------------------------------------------------------------------
# Section 6: Detailed exp3 per-model table
# ---------------------------------------------------------------------------
def print_exp3_detail():
    print("\n" + "=" * 80)
    print("DETAIL: Exp3 Property-Constrained Per-Model Results")
    print("=" * 80)

    # ResToken constraint satisfaction
    print("\n--- ResToken Constraint Satisfaction ---")
    header = f"{'Model':<20} {'permeable':>12} {'charged':>12} {'rigid':>12}"
    print(header)
    print("-" * len(header))
    for model in MODELS:
        mname = MODEL_SHORT.get(model, model)
        vals = []
        for profile in PROFILES:
            summary = get_exp3_summary(model, "restoken", profile)
            if summary and "constraint_satisfaction" in summary:
                vals.append(f"{summary['constraint_satisfaction']:.3f}")
            else:
                vals.append("N/A")
        print(f"{mname:<20} {vals[0]:>12} {vals[1]:>12} {vals[2]:>12}")

    # Overall validity across representations for 'permeable' (most data)
    print("\n--- Overall Validity (permeable profile) ---")
    header = f"{'Model':<20} {'restoken':>12} {'smiles':>12} {'helm':>12}"
    print(header)
    print("-" * len(header))
    for model in MODELS:
        mname = MODEL_SHORT.get(model, model)
        vals = []
        for rep in REPRESENTATIONS:
            summary = get_exp3_summary(model, rep, "permeable")
            if summary and summary.get("n_parsed", 0) > 0:
                vals.append(f"{summary.get('overall_valid', 0.0):.3f}")
            else:
                vals.append("N/A")
        print(f"{mname:<20} {vals[0]:>12} {vals[1]:>12} {vals[2]:>12}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    print("ResToken W2 Benchmark — Statistical Analysis")
    print(f"Models: {', '.join(MODEL_SHORT.values())}")
    print(f"Bootstrap resamples: {N_BOOTSTRAP}")
    print(f"Random seed: {SEED}")
    print()

    # Section 1
    bootstrap_cis = compute_bootstrap_cis()

    # Section 2
    chi_sq_results = compute_chi_square_tests(bootstrap_cis)

    # Section 3
    aggregate = compute_aggregate_stats()

    # Section 4
    enrichment = compute_enrichment()

    # Detail tables
    print_exp2_detail()
    print_exp3_detail()

    # ---------------------------------------------------------------------------
    # Assemble JSON output
    # ---------------------------------------------------------------------------
    output = {
        "metadata": {
            "models": MODELS,
            "representations": REPRESENTATIONS,
            "profiles": PROFILES,
            "n_bootstrap": N_BOOTSTRAP,
            "seed": SEED,
            "alpha": ALPHA,
            "random_baselines": RANDOM_BASELINES,
        },
        "exp1_bootstrap_ci": {},
        "exp1_chi_square": {},
        "exp2_controllability_per_model": {},
        "exp3_constraint_satisfaction_per_model": {},
        "aggregate": aggregate,
        "enrichment": enrichment,
    }

    # Restructure bootstrap CIs
    for model in MODELS:
        output["exp1_bootstrap_ci"][model] = {}
        for rep in REPRESENTATIONS:
            key = f"{model}|{rep}"
            if key in bootstrap_cis:
                output["exp1_bootstrap_ci"][model][rep] = bootstrap_cis[key]

    # Restructure chi-square tests
    for model in MODELS:
        output["exp1_chi_square"][model] = {}
        for comp in ["vs SMILES", "vs HELM"]:
            key = f"{model}|restoken {comp}"
            if key in chi_sq_results:
                output["exp1_chi_square"][model][comp] = chi_sq_results[key]

    # Exp2 per model
    for model in MODELS:
        summary = get_exp2_summary(model)
        if summary and summary.get("n_total_variants", 0) > 0:
            output["exp2_controllability_per_model"][model] = {
                k: summary[k] for k in
                ["edit_compliance", "frozen_compliance", "id_valid",
                 "length_match", "n_total_variants"]
            }

    # Exp3 per model (restoken constraint satisfaction)
    for model in MODELS:
        output["exp3_constraint_satisfaction_per_model"][model] = {}
        for profile in PROFILES:
            summary = get_exp3_summary(model, "restoken", profile)
            if summary and "constraint_satisfaction" in summary:
                output["exp3_constraint_satisfaction_per_model"][model][profile] = {
                    "constraint_satisfaction": summary["constraint_satisfaction"],
                    "n_constraint_pass": summary.get("n_constraint_pass", 0),
                    "n_valid": summary.get("n_valid", 0),
                    "n_parsed": summary.get("n_parsed", 0),
                }

    # Write JSON
    with open(OUT_JSON, "w") as f:
        json.dump(output, f, indent=2)

    print(f"\n{'=' * 80}")
    print(f"JSON summary written to: {OUT_JSON}")
    print(f"{'=' * 80}")


if __name__ == "__main__":
    main()
