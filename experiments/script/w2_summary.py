#!/usr/bin/env python3
"""W2 Benchmark Summary — aggregate all results into tables and figures.

Usage:
    python w2_summary.py                  # Print tables to stdout
    python w2_summary.py --figures        # Also generate figures
    python w2_summary.py --latex          # Output LaTeX tables
"""

import argparse
import json
import sys
from pathlib import Path
from collections import defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

OUTPUT_ROOT = Path(__file__).resolve().parents[1] / "output"
FIG_DIR = Path(__file__).resolve().parents[2] / "manuscript" / "figures" / "generated"
TABLE_DIR = Path(__file__).resolve().parents[2] / "manuscript" / "tables"

MODELS = [
    "gemini-2.5-pro",
    "gemini-2.5-flash",
    "gpt-4o",
    "claude-sonnet-4-6",
    "Qwen3.5-9B",
    "gemma-3-12b-it-bnb-4bit",
    "txgemma-9b-chat",
]

MODEL_DISPLAY = {
    "gemini-2.5-pro": "Gemini 2.5 Pro",
    "gemini-2.5-flash": "Gemini 2.5 Flash",
    "gpt-4o": "GPT-4o",
    "claude-sonnet-4-6": "Claude Sonnet 4.6",
    "Qwen3.5-9B": "Qwen 3.5 9B",
    "gemma-3-12b-it-bnb-4bit": "Gemma 3 12B",
    "txgemma-9b-chat": "TxGemma 9B",
}

REPS = ["restoken", "smiles", "helm"]
PROFILES = ["permeable", "charged_binder", "rigid_scaffold"]


def load_summary(exp_dir, pattern):
    p = OUTPUT_ROOT / exp_dir / pattern
    matches = list((OUTPUT_ROOT / exp_dir).glob(pattern))
    if not matches:
        return None
    try:
        return json.loads(matches[0].read_text())
    except Exception:
        return None


def exp1_table():
    """Exp1: Zero-shot validity rates by model × representation."""
    print("\n" + "=" * 80)
    print("EXPERIMENT 1: Zero-Shot Validity (n=200 requested)")
    print("=" * 80)

    header = f"{'Model':<25} {'ResToken':>12} {'SMILES':>12} {'HELM':>12} {'ResToken':>10} {'SMILES':>10} {'HELM':>10}"
    print(f"\n{'':<25} {'--- n_parsed ---':>36} {'--- validity% ---':>30}")
    print(header)
    print("-" * 95)

    rows = []
    for model in MODELS:
        row = {"model": model}
        for rep in REPS:
            d = load_summary("exp1_validity",
                             f"{model}_{rep}_unconstrained_len6_summary.json")
            if d and d.get("n_parsed", 0) > 0:
                n = d.get("n_parsed", 0)
                n_valid = d.get("n_valid", 0)
                rate = n_valid / n * 100 if n > 0 else 0
                row[f"{rep}_n"] = n
                row[f"{rep}_rate"] = rate
            else:
                row[f"{rep}_n"] = 0
                row[f"{rep}_rate"] = 0
        rows.append(row)

        name = MODEL_DISPLAY.get(model, model)
        rn = row.get("restoken_n", 0)
        sn = row.get("smiles_n", 0)
        hn = row.get("helm_n", 0)
        rr = row.get("restoken_rate", 0)
        sr = row.get("smiles_rate", 0)
        hr = row.get("helm_rate", 0)
        print(f"{name:<25} {rn:>12} {sn:>12} {hn:>12} {rr:>9.1f}% {sr:>9.1f}% {hr:>9.1f}%")

    return rows


def exp1_diversity_table():
    """Exp1: Diversity metrics (uniqueness) by model × representation."""
    print("\n" + "=" * 80)
    print("EXPERIMENT 1: Sequence Diversity (uniqueness = unique/total)")
    print("=" * 80)

    print(f"\n{'Model':<25} {'ResToken':>12} {'SMILES':>12} {'HELM':>12}")
    print("-" * 65)

    for model in MODELS:
        name = MODEL_DISPLAY.get(model, model)
        vals = []
        for rep in REPS:
            d = load_summary("exp1_validity",
                             f"{model}_{rep}_unconstrained_len6_summary.json")
            if d:
                u = d.get("uniqueness", d.get("diversity", {}).get("uniqueness", None))
                if u is not None:
                    vals.append(f"{u:.3f}")
                else:
                    vals.append("—")
            else:
                vals.append("—")
        print(f"{name:<25} {vals[0]:>12} {vals[1]:>12} {vals[2]:>12}")


def exp2_table():
    """Exp2: Controllability metrics."""
    print("\n" + "=" * 80)
    print("EXPERIMENT 2: EDIT/FROZEN Controllability (ResToken, 10 parents × 20 variants)")
    print("=" * 80)

    print(f"\n{'Model':<25} {'n_variants':>12} {'length%':>10} {'frozen%':>10} {'edit%':>10} {'id_valid%':>10}")
    print("-" * 80)

    for model in MODELS:
        d = load_summary("exp2_controllability",
                         f"{model}_exp2_controllability_summary.json")
        name = MODEL_DISPLAY.get(model, model)
        if d and d.get("n_total_variants", 0) > 0:
            n = d["n_total_variants"]
            lm = d.get("length_match", 0) * 100
            fc = d.get("frozen_compliance", 0) * 100
            ec = d.get("edit_compliance", 0) * 100
            iv = d.get("id_valid", 0) * 100
            print(f"{name:<25} {n:>12} {lm:>9.1f}% {fc:>9.1f}% {ec:>9.1f}% {iv:>9.1f}%")
        else:
            print(f"{name:<25} {'0':>12} {'—':>10} {'—':>10} {'—':>10} {'—':>10}")


def exp3_table():
    """Exp3: Property-constrained generation."""
    print("\n" + "=" * 80)
    print("EXPERIMENT 3: Property-Constrained Generation (ResToken, n=100 requested)")
    print("=" * 80)

    print(f"\n{'Model':<25} {'Profile':<20} {'n_parsed':>10} {'n_valid':>10} {'constraint%':>12}")
    print("-" * 80)

    for model in MODELS:
        name = MODEL_DISPLAY.get(model, model)
        first = True
        for profile in PROFILES:
            d = load_summary("exp3_property_constrained",
                             f"{model}_restoken_{profile}_len6_summary.json")
            label = name if first else ""
            first = False
            if d:
                n_parsed = d.get("n_parsed", 0)
                n_valid = d.get("n_valid", 0)
                cs = d.get("constraint_satisfaction", 0) * 100
                print(f"{label:<25} {profile:<20} {n_parsed:>10} {n_valid:>10} {cs:>11.1f}%")
            else:
                print(f"{label:<25} {profile:<20} {'—':>10} {'—':>10} {'—':>12}")


def fig_exp1_validity(rows):
    """Bar chart: validity rate by model × representation."""
    FIG_DIR.mkdir(parents=True, exist_ok=True)

    models_with_data = [r for r in rows if any(r.get(f"{rep}_n", 0) > 0 for rep in REPS)]
    if not models_with_data:
        print("  [SKIP] No data for exp1 figure")
        return

    n_models = len(models_with_data)
    x = np.arange(n_models)
    width = 0.25

    fig, ax = plt.subplots(figsize=(7, 4))
    colors = {"restoken": "#4C72B0", "smiles": "#C44E52", "helm": "#8172B2"}

    for i, rep in enumerate(REPS):
        vals = [r.get(f"{rep}_rate", 0) for r in models_with_data]
        ax.bar(x + (i - 1) * width, vals, width, label=rep.upper(),
               color=colors[rep], edgecolor="white", linewidth=0.5)

    ax.set_ylabel("Validity Rate (%)")
    ax.set_title("Experiment 1: Zero-Shot Sequence Validity")
    ax.set_xticks(x)
    ax.set_xticklabels([MODEL_DISPLAY.get(r["model"], r["model"]) for r in models_with_data],
                       rotation=30, ha="right", fontsize=7)
    ax.legend(frameon=False)
    ax.set_ylim(0, 105)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    out = FIG_DIR / "w2_exp1_validity.png"
    fig.savefig(out, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out}")


def fig_exp2_compliance():
    """Grouped bar: edit and frozen compliance by model."""
    FIG_DIR.mkdir(parents=True, exist_ok=True)

    data = []
    for model in MODELS:
        d = load_summary("exp2_controllability",
                         f"{model}_exp2_controllability_summary.json")
        if d and d.get("n_total_variants", 0) > 0:
            data.append({
                "model": model,
                "frozen": d.get("frozen_compliance", 0) * 100,
                "edit": d.get("edit_compliance", 0) * 100,
                "id_valid": d.get("id_valid", 0) * 100,
            })

    if not data:
        print("  [SKIP] No data for exp2 figure")
        return

    x = np.arange(len(data))
    width = 0.25
    fig, ax = plt.subplots(figsize=(6, 4))

    ax.bar(x - width, [d["frozen"] for d in data], width, label="Frozen Compliance",
           color="#55A868", edgecolor="white", linewidth=0.5)
    ax.bar(x, [d["edit"] for d in data], width, label="Edit Compliance",
           color="#4C72B0", edgecolor="white", linewidth=0.5)
    ax.bar(x + width, [d["id_valid"] for d in data], width, label="ID Valid",
           color="#DD8452", edgecolor="white", linewidth=0.5)

    ax.set_ylabel("Rate (%)")
    ax.set_title("Experiment 2: EDIT/FROZEN Controllability")
    ax.set_xticks(x)
    ax.set_xticklabels([MODEL_DISPLAY.get(d["model"], d["model"]) for d in data],
                       rotation=30, ha="right", fontsize=7)
    ax.legend(frameon=False, fontsize=7)
    ax.set_ylim(0, 105)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    out = FIG_DIR / "w2_exp2_controllability.png"
    fig.savefig(out, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out}")


def fig_exp3_constraints():
    """Heatmap: constraint satisfaction by model × profile."""
    FIG_DIR.mkdir(parents=True, exist_ok=True)

    models_with_data = []
    matrix = []
    for model in MODELS:
        row_vals = []
        has_data = False
        for profile in PROFILES:
            d = load_summary("exp3_property_constrained",
                             f"{model}_restoken_{profile}_len6_summary.json")
            if d and d.get("n_valid", 0) > 0:
                row_vals.append(d.get("constraint_satisfaction", 0) * 100)
                has_data = True
            else:
                row_vals.append(np.nan)
        if has_data:
            models_with_data.append(model)
            matrix.append(row_vals)

    if not matrix:
        print("  [SKIP] No data for exp3 figure")
        return

    matrix = np.array(matrix)
    fig, ax = plt.subplots(figsize=(5, max(3, len(models_with_data) * 0.6)))

    cmap = plt.cm.RdYlGn.copy()
    cmap.set_bad(color="#f0f0f0")
    im = ax.imshow(matrix, cmap=cmap, aspect="auto", vmin=0, vmax=100)

    ax.set_xticks(range(len(PROFILES)))
    ax.set_xticklabels([p.replace("_", "\n") for p in PROFILES], fontsize=7)
    ax.set_yticks(range(len(models_with_data)))
    ax.set_yticklabels([MODEL_DISPLAY.get(m, m) for m in models_with_data], fontsize=7)

    for i in range(len(models_with_data)):
        for j in range(len(PROFILES)):
            val = matrix[i, j]
            if np.isnan(val):
                ax.text(j, i, "—", ha="center", va="center", fontsize=7, color="#999")
            else:
                color = "white" if val > 60 else "black"
                ax.text(j, i, f"{val:.0f}%", ha="center", va="center",
                        fontsize=7, fontweight="bold", color=color)

    ax.set_title("Experiment 3: Constraint Satisfaction (ResToken)", fontsize=9)
    fig.colorbar(im, ax=ax, shrink=0.7, label="Satisfaction %")

    out = FIG_DIR / "w2_exp3_constraints.png"
    fig.savefig(out, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out}")


def _esc(s):
    return s.replace("_", r"\_").replace("%", r"\%").replace("&", r"\&")


def latex_exp1_validity():
    """LaTeX Table 2: Zero-shot validity by model × representation."""
    lines = [
        r"\begin{table}[htbp]",
        r"\centering",
        r"\caption{Zero-shot sequence validity across models and representations "
        r"(Experiment~1, $n=200$ requested per condition). "
        r"Validity = fraction of parsed sequences passing all seven constraint checks. "
        r"TxGemma~9B failed to produce any parseable output and is omitted.}",
        r"\label{tab:w2-exp1-validity}",
        r"\small",
        r"\begin{tabular}{l *{3}{r} *{3}{r}}",
        r"\toprule",
        r"& \multicolumn{3}{c}{\textbf{Parsed ($n$)}} & \multicolumn{3}{c}{\textbf{Validity (\%)}} \\",
        r"\cmidrule(lr){2-4} \cmidrule(lr){5-7}",
        r"Model & ResToken & SMILES & HELM & ResToken & SMILES & HELM \\",
        r"\midrule",
    ]

    best = {"restoken": 0, "smiles": 0, "helm": 0}
    rows_data = []
    for model in MODELS:
        row = {"model": model}
        for rep in REPS:
            d = load_summary("exp1_validity",
                             f"{model}_{rep}_unconstrained_len6_summary.json")
            if d and d.get("n_parsed", 0) > 0:
                n = d["n_parsed"]
                nv = d.get("n_valid", 0)
                rate = nv / n * 100
                row[f"{rep}_n"] = n
                row[f"{rep}_rate"] = rate
                if rate > best[rep]:
                    best[rep] = rate
            else:
                row[f"{rep}_n"] = 0
                row[f"{rep}_rate"] = None
        rows_data.append(row)

    for row in rows_data:
        if row.get("restoken_n", 0) == 0 and row.get("smiles_n", 0) == 0 and row.get("helm_n", 0) == 0:
            continue
        name = MODEL_DISPLAY.get(row["model"], row["model"])
        cells = [_esc(name)]
        for rep in REPS:
            n = row.get(f"{rep}_n", 0)
            cells.append(str(n) if n > 0 else "---")
        for rep in REPS:
            r = row.get(f"{rep}_rate")
            if r is not None:
                fmt = f"{r:.1f}"
                if r == best[rep]:
                    fmt = r"\textbf{" + fmt + "}"
                cells.append(fmt)
            else:
                cells.append("---")
        lines.append(" & ".join(cells) + r" \\")

    lines += [
        r"\bottomrule",
        r"\end{tabular}",
        r"\end{table}",
    ]
    return "\n".join(lines)


def latex_exp1_diversity():
    """LaTeX table: Sequence diversity (uniqueness) by model × representation."""
    lines = [
        r"\begin{table}[htbp]",
        r"\centering",
        r"\caption{Sequence diversity (uniqueness = unique valid sequences / total valid) "
        r"by model and representation. Higher values indicate greater generation diversity. "
        r"Models with no valid output are omitted.}",
        r"\label{tab:w2-exp1-diversity}",
        r"\small",
        r"\begin{tabular}{l rrr}",
        r"\toprule",
        r"Model & ResToken & SMILES & HELM \\",
        r"\midrule",
    ]

    for model in MODELS:
        vals = []
        any_data = False
        for rep in REPS:
            d = load_summary("exp1_validity",
                             f"{model}_{rep}_unconstrained_len6_summary.json")
            if d:
                u = d.get("uniqueness", d.get("diversity", {}).get("uniqueness"))
                if u is not None:
                    vals.append(f"{u:.3f}")
                    any_data = True
                else:
                    vals.append("---")
            else:
                vals.append("---")
        if not any_data:
            continue
        name = _esc(MODEL_DISPLAY.get(model, model))
        lines.append(f"{name} & {vals[0]} & {vals[1]} & {vals[2]}" + r" \\")

    lines += [
        r"\bottomrule",
        r"\end{tabular}",
        r"\end{table}",
    ]
    return "\n".join(lines)


def latex_exp2():
    """LaTeX table: Exp2 controllability metrics."""
    lines = [
        r"\begin{table}[htbp]",
        r"\centering",
        r"\caption{EDIT/FROZEN controllability (Experiment~2, ResToken only, "
        r"10~parent sequences $\times$ 20~variants each). "
        r"Length = fraction matching parent length; "
        r"Frozen = fraction of FROZEN positions unchanged; "
        r"Edit = fraction of EDIT positions changed; "
        r"ID~Valid = fraction with all IDs in library.}",
        r"\label{tab:w2-exp2-controllability}",
        r"\small",
        r"\begin{tabular}{l r rrrr}",
        r"\toprule",
        r"Model & $n$ & Length (\%) & Frozen (\%) & Edit (\%) & ID Valid (\%) \\",
        r"\midrule",
    ]

    for model in MODELS:
        d = load_summary("exp2_controllability",
                         f"{model}_exp2_controllability_summary.json")
        if not d or d.get("n_total_variants", 0) == 0:
            continue
        name = _esc(MODEL_DISPLAY.get(model, model))
        n = d["n_total_variants"]
        lm = d.get("length_match", 0) * 100
        fc = d.get("frozen_compliance", 0) * 100
        ec = d.get("edit_compliance", 0) * 100
        iv = d.get("id_valid", 0) * 100
        lines.append(
            f"{name} & {n} & {lm:.1f} & {fc:.1f} & {ec:.1f} & {iv:.1f}" + r" \\"
        )

    lines += [
        r"\bottomrule",
        r"\end{tabular}",
        r"\end{table}",
    ]
    return "\n".join(lines)


def latex_exp3():
    """LaTeX table: Exp3 property-constrained generation."""
    profile_display = {
        "permeable": "Permeable",
        "charged_binder": "Charged binder",
        "rigid_scaffold": "Rigid scaffold",
    }
    lines = [
        r"\begin{table}[htbp]",
        r"\centering",
        r"\caption{Property-constrained generation (Experiment~3, ResToken, "
        r"$n=100$ requested per condition). "
        r"Constraint satisfaction = fraction of valid sequences meeting all profile constraints. "
        r"Random baseline acceptance rates: permeable 2.07\%, charged 1.82\%, rigid 0.81\%.}",
        r"\label{tab:w2-exp3-constraints}",
        r"\small",
        r"\begin{tabular}{l l rrr}",
        r"\toprule",
        r"Model & Profile & Parsed & Valid & Constr.\ Sat.\ (\%) \\",
        r"\midrule",
    ]

    for model in MODELS:
        name = _esc(MODEL_DISPLAY.get(model, model))
        any_data = False
        model_lines = []
        first = True
        for profile in PROFILES:
            d = load_summary("exp3_property_constrained",
                             f"{model}_restoken_{profile}_len6_summary.json")
            label = name if first else ""
            first = False
            if d and d.get("n_parsed", 0) > 0:
                any_data = True
                np_ = d["n_parsed"]
                nv = d.get("n_valid", 0)
                cs = d.get("constraint_satisfaction", 0) * 100
                cs_fmt = f"\\textbf{{{cs:.1f}}}" if cs >= 80 else f"{cs:.1f}"
                model_lines.append(
                    f"{label} & {profile_display[profile]} & {np_} & {nv} & {cs_fmt}" + r" \\"
                )
            else:
                model_lines.append(
                    f"{label} & {profile_display[profile]} & --- & --- & ---" + r" \\"
                )
        if any_data:
            lines.extend(model_lines)
            lines.append(r"\addlinespace")

    if lines[-1] == r"\addlinespace":
        lines.pop()

    lines += [
        r"\bottomrule",
        r"\end{tabular}",
        r"\end{table}",
    ]
    return "\n".join(lines)


def latex_combined():
    """Combined supplementary table (Table 2 per manuscript scaffold)."""
    lines = [
        r"\begin{table*}[htbp]",
        r"\centering",
        r"\caption{Complete benchmark results across all models, representations, and experiments. "
        r"Validity = fraction passing all constraint checks; "
        r"Uniq = unique/total ratio; "
        r"Edit/Frozen = compliance in EDIT/FROZEN controllability test (ResToken only); "
        r"Constr = mean constraint satisfaction across three profiles (ResToken only). "
        r"Models with zero parseable output in all conditions are omitted.}",
        r"\label{tab:w2-combined}",
        r"\footnotesize",
        r"\begin{tabular}{l *{3}{r} *{3}{r} rr r}",
        r"\toprule",
        r"& \multicolumn{3}{c}{Validity (\%)} & \multicolumn{3}{c}{Uniqueness} "
        r"& \multicolumn{2}{c}{Controllability (\%)} & Constr. \\",
        r"\cmidrule(lr){2-4} \cmidrule(lr){5-7} \cmidrule(lr){8-9} \cmidrule(lr){10-10}",
        r"Model & RT & SM & HE & RT & SM & HE & Edit & Frozen & Sat (\%) \\",
        r"\midrule",
    ]

    for model in MODELS:
        vals_v, vals_u = [], []
        any_data = False
        for rep in REPS:
            d = load_summary("exp1_validity",
                             f"{model}_{rep}_unconstrained_len6_summary.json")
            if d and d.get("n_parsed", 0) > 0:
                any_data = True
                n = d["n_parsed"]
                nv = d.get("n_valid", 0)
                vals_v.append(f"{nv/n*100:.1f}")
                u = d.get("uniqueness", d.get("diversity", {}).get("uniqueness"))
                vals_u.append(f"{u:.2f}" if u is not None else "---")
            else:
                vals_v.append("---")
                vals_u.append("---")

        if not any_data:
            continue

        d2 = load_summary("exp2_controllability",
                          f"{model}_exp2_controllability_summary.json")
        if d2 and d2.get("n_total_variants", 0) > 0:
            ec = f"{d2.get('edit_compliance', 0)*100:.1f}"
            fc = f"{d2.get('frozen_compliance', 0)*100:.1f}"
        else:
            ec, fc = "---", "---"

        cs_vals = []
        for profile in PROFILES:
            d3 = load_summary("exp3_property_constrained",
                              f"{model}_restoken_{profile}_len6_summary.json")
            if d3 and d3.get("n_valid", 0) > 0:
                cs_vals.append(d3.get("constraint_satisfaction", 0) * 100)
        mean_cs = f"{np.mean(cs_vals):.1f}" if cs_vals else "---"

        name = _esc(MODEL_DISPLAY.get(model, model))
        cells = [name] + vals_v + vals_u + [ec, fc, mean_cs]
        lines.append(" & ".join(cells) + r" \\")

    lines += [
        r"\bottomrule",
        r"\end{tabular}",
        r"\end{table*}",
    ]
    return "\n".join(lines)


def write_latex_tables():
    TABLE_DIR.mkdir(parents=True, exist_ok=True)
    tables = {
        "w2_table2_validity.tex": latex_exp1_validity(),
        "w2_table_diversity.tex": latex_exp1_diversity(),
        "w2_table_controllability.tex": latex_exp2(),
        "w2_table_constraints.tex": latex_exp3(),
        "w2_table_combined.tex": latex_combined(),
    }
    for fname, content in tables.items():
        out = TABLE_DIR / fname
        out.write_text(content)
        print(f"  Saved: {out}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--figures", action="store_true", help="Generate figures")
    parser.add_argument("--latex", action="store_true", help="Output LaTeX tables")
    args = parser.parse_args()

    rows = exp1_table()
    exp1_diversity_table()
    exp2_table()
    exp3_table()

    if args.latex:
        print("\n\nGenerating LaTeX tables...")
        write_latex_tables()

    if args.figures:
        print("\n\nGenerating figures...")
        fig_exp1_validity(rows)
        fig_exp2_compliance()
        fig_exp3_constraints()

    print("\nDone.")


if __name__ == "__main__":
    main()
