#!/usr/bin/env python3
"""W2 Benchmark — Publication-quality figures for JCIM paper.

Generates:
  Fig A: Exp1 validity grouped bars (model × rep)
  Fig B: Exp1 diversity (uniqueness) grouped bars
  Fig C: Exp2 controllability stacked metrics
  Fig D: Exp3 constraint satisfaction heatmap
  Fig E: Combined 2×2 panel figure
  Fig F: Radar chart per model (overall capability)
"""

import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

OUTPUT_ROOT = Path(__file__).resolve().parents[1] / "output"
FIG_DIR = Path(__file__).resolve().parents[2] / "manuscript" / "figures" / "generated"
FIG_DIR.mkdir(parents=True, exist_ok=True)

# ─── Style ────────────────────────────────────────────────────────────────────
plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 8,
    "axes.titlesize": 10,
    "axes.labelsize": 9,
    "xtick.labelsize": 7.5,
    "ytick.labelsize": 7.5,
    "legend.fontsize": 7.5,
    "figure.dpi": 300,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
    "axes.linewidth": 0.6,
    "xtick.major.width": 0.6,
    "ytick.major.width": 0.6,
    "axes.spines.top": False,
    "axes.spines.right": False,
})

# Exclude txgemma (complete failure) from main figures
MODELS = [
    "gemini-2.5-pro",
    "gemini-2.5-flash",
    "gpt-4o",
    "Qwen3.5-9B",
    "gemma-3-12b-it-bnb-4bit",
]

MODEL_SHORT = {
    "gemini-2.5-pro": "Gemini\n2.5 Pro",
    "gemini-2.5-flash": "Gemini\n2.5 Flash",
    "gpt-4o": "GPT-4o",
    "Qwen3.5-9B": "Qwen 3.5\n9B",
    "gemma-3-12b-it-bnb-4bit": "Gemma 3\n12B",
}

MODEL_COLORS = {
    "gemini-2.5-pro": "#4285F4",
    "gemini-2.5-flash": "#34A853",
    "gpt-4o": "#10A37F",
    "Qwen3.5-9B": "#7C3AED",
    "gemma-3-12b-it-bnb-4bit": "#EA4335",
}

REP_COLORS = {
    "restoken": "#2563EB",
    "smiles": "#DC2626",
    "helm": "#7C3AED",
}

REPS = ["restoken", "smiles", "helm"]
PROFILES = ["permeable", "charged_binder", "rigid_scaffold"]
PROFILE_DISPLAY = {"permeable": "Permeable", "charged_binder": "Charged\nBinder", "rigid_scaffold": "Rigid\nScaffold"}


def load(exp_dir, pattern):
    matches = list((OUTPUT_ROOT / exp_dir).glob(pattern))
    if not matches:
        return None
    try:
        return json.loads(matches[0].read_text())
    except Exception:
        return None


# ─── Fig A: Exp1 Validity ─────────────────────────────────────────────────────

def fig_exp1_validity():
    fig, ax = plt.subplots(figsize=(7, 3.5))
    x = np.arange(len(MODELS))
    width = 0.22
    offsets = [-width, 0, width]

    for i, rep in enumerate(REPS):
        vals = []
        for model in MODELS:
            d = load("exp1_validity", f"{model}_{rep}_unconstrained_len6_summary.json")
            if d and d.get("n_parsed", 0) > 0:
                n_valid = d.get("n_valid", 0)
                n_parsed = d.get("n_parsed", 1)
                vals.append(n_valid / n_parsed * 100)
            else:
                vals.append(0)
        bars = ax.bar(x + offsets[i], vals, width * 0.9, label=rep.upper(),
                      color=REP_COLORS[rep], edgecolor="white", linewidth=0.3, alpha=0.85)
        for bar, v in zip(bars, vals):
            if v > 0:
                ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 1,
                        f"{v:.0f}", ha="center", va="bottom", fontsize=5.5, color="#333")

    ax.set_ylabel("Validity Rate (%)")
    ax.set_title("(a) Zero-Shot Sequence Validity", fontweight="bold", loc="left")
    ax.set_xticks(x)
    ax.set_xticklabels([MODEL_SHORT[m] for m in MODELS], fontsize=7)
    ax.set_ylim(0, 115)
    ax.legend(frameon=False, ncol=3, loc="upper right")
    ax.axhline(y=100, color="#ccc", linewidth=0.5, linestyle="--", zorder=0)

    out = FIG_DIR / "w2_fig_a_validity.pdf"
    fig.savefig(out)
    fig.savefig(out.with_suffix(".png"))
    plt.close(fig)
    print(f"  {out.name}")
    return fig


# ─── Fig B: Exp1 Diversity ────────────────────────────────────────────────────

def fig_exp1_diversity():
    fig, ax = plt.subplots(figsize=(7, 3.5))
    x = np.arange(len(MODELS))
    width = 0.22
    offsets = [-width, 0, width]

    for i, rep in enumerate(REPS):
        vals = []
        for model in MODELS:
            d = load("exp1_validity", f"{model}_{rep}_unconstrained_len6_summary.json")
            if d:
                u = d.get("uniqueness", 0)
                vals.append(u * 100)
            else:
                vals.append(0)
        bars = ax.bar(x + offsets[i], vals, width * 0.9, label=rep.upper(),
                      color=REP_COLORS[rep], edgecolor="white", linewidth=0.3, alpha=0.85)
        for bar, v in zip(bars, vals):
            if v > 0:
                ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 1,
                        f"{v:.0f}", ha="center", va="bottom", fontsize=5.5, color="#333")

    ax.set_ylabel("Uniqueness (%)")
    ax.set_title("(b) Sequence Diversity (Unique / Total)", fontweight="bold", loc="left")
    ax.set_xticks(x)
    ax.set_xticklabels([MODEL_SHORT[m] for m in MODELS], fontsize=7)
    ax.set_ylim(0, 115)
    ax.legend(frameon=False, ncol=3, loc="upper right")
    ax.axhline(y=100, color="#ccc", linewidth=0.5, linestyle="--", zorder=0)

    out = FIG_DIR / "w2_fig_b_diversity.pdf"
    fig.savefig(out)
    fig.savefig(out.with_suffix(".png"))
    plt.close(fig)
    print(f"  {out.name}")


# ─── Fig C: Exp2 Controllability ──────────────────────────────────────────────

def fig_exp2_controllability():
    fig, ax = plt.subplots(figsize=(6, 3.5))

    data = []
    labels = []
    for model in MODELS:
        d = load("exp2_controllability", f"{model}_exp2_controllability_summary.json")
        if d and d.get("n_total_variants", 0) > 0:
            data.append(d)
            labels.append(MODEL_SHORT[model])

    x = np.arange(len(data))
    width = 0.2

    metrics = [
        ("edit_compliance", "Edit Compliance", "#2563EB"),
        ("frozen_compliance", "Frozen Compliance", "#16A34A"),
        ("id_valid", "ID Valid", "#F59E0B"),
    ]

    for i, (key, label, color) in enumerate(metrics):
        vals = [d.get(key, 0) * 100 for d in data]
        bars = ax.bar(x + (i - 1) * width, vals, width * 0.85, label=label,
                      color=color, edgecolor="white", linewidth=0.3, alpha=0.85)
        for bar, v in zip(bars, vals):
            if v > 0:
                ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 1,
                        f"{v:.0f}", ha="center", va="bottom", fontsize=5.5, color="#333")

    ax.set_ylabel("Rate (%)")
    ax.set_title("(c) EDIT/FROZEN Controllability (ResToken)", fontweight="bold", loc="left")
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=7)
    ax.set_ylim(0, 115)
    ax.legend(frameon=False, ncol=3, fontsize=7)
    ax.axhline(y=100, color="#ccc", linewidth=0.5, linestyle="--", zorder=0)

    out = FIG_DIR / "w2_fig_c_controllability.pdf"
    fig.savefig(out)
    fig.savefig(out.with_suffix(".png"))
    plt.close(fig)
    print(f"  {out.name}")


# ─── Fig D: Exp3 Constraint Heatmap ──────────────────────────────────────────

def fig_exp3_heatmap():
    models_with_data = []
    matrix = []
    n_valid_matrix = []

    for model in MODELS:
        row = []
        nv_row = []
        has_any = False
        for profile in PROFILES:
            d = load("exp3_property_constrained",
                     f"{model}_restoken_{profile}_len6_summary.json")
            if d and d.get("n_valid", 0) > 0:
                row.append(d.get("constraint_satisfaction", 0) * 100)
                nv_row.append(d.get("n_valid", 0))
                has_any = True
            elif d and d.get("n_parsed", 0) > 0:
                row.append(0)
                nv_row.append(d.get("n_parsed", 0))
                has_any = True
            else:
                row.append(np.nan)
                nv_row.append(0)
        if has_any:
            models_with_data.append(model)
            matrix.append(row)
            n_valid_matrix.append(nv_row)

    matrix = np.array(matrix)
    fig, ax = plt.subplots(figsize=(4.5, 3.5))

    cmap = plt.cm.RdYlGn.copy()
    cmap.set_bad(color="#f5f5f5")
    im = ax.imshow(matrix, cmap=cmap, aspect="auto", vmin=0, vmax=100)

    ax.set_xticks(range(len(PROFILES)))
    ax.set_xticklabels([PROFILE_DISPLAY[p] for p in PROFILES], fontsize=8)
    ax.set_yticks(range(len(models_with_data)))
    ax.set_yticklabels([MODEL_SHORT[m].replace("\n", " ") for m in models_with_data], fontsize=7.5)

    for i in range(len(models_with_data)):
        for j in range(len(PROFILES)):
            val = matrix[i, j]
            nv = n_valid_matrix[i][j]
            if np.isnan(val):
                ax.text(j, i, "N/A", ha="center", va="center", fontsize=7, color="#aaa")
            else:
                color = "white" if val > 55 else ("black" if val > 0 else "#666")
                txt = f"{val:.0f}%\n(n={nv})"
                ax.text(j, i, txt, ha="center", va="center", fontsize=6.5,
                        fontweight="bold", color=color)

    ax.set_title("(d) Constraint Satisfaction (ResToken)", fontweight="bold", loc="left", fontsize=10)
    cbar = fig.colorbar(im, ax=ax, shrink=0.8, pad=0.08)
    cbar.set_label("Satisfaction %", fontsize=8)
    cbar.ax.tick_params(labelsize=7)

    out = FIG_DIR / "w2_fig_d_constraints.pdf"
    fig.savefig(out)
    fig.savefig(out.with_suffix(".png"))
    plt.close(fig)
    print(f"  {out.name}")


# ─── Fig E: Combined 2×2 Panel ───────────────────────────────────────────────

def fig_combined_panel():
    fig = plt.figure(figsize=(14, 10))
    gs = fig.add_gridspec(2, 2, hspace=0.35, wspace=0.3)

    # (a) Validity
    ax1 = fig.add_subplot(gs[0, 0])
    x = np.arange(len(MODELS))
    width = 0.22
    offsets = [-width, 0, width]
    for i, rep in enumerate(REPS):
        vals = []
        for model in MODELS:
            d = load("exp1_validity", f"{model}_{rep}_unconstrained_len6_summary.json")
            if d and d.get("n_parsed", 0) > 0:
                vals.append(d.get("n_valid", 0) / d.get("n_parsed", 1) * 100)
            else:
                vals.append(0)
        ax1.bar(x + offsets[i], vals, width * 0.9, label=rep.upper(),
                color=REP_COLORS[rep], edgecolor="white", linewidth=0.3, alpha=0.85)
    ax1.set_ylabel("Validity Rate (%)")
    ax1.set_title("(a) Zero-Shot Validity", fontweight="bold", loc="left")
    ax1.set_xticks(x)
    ax1.set_xticklabels([MODEL_SHORT[m] for m in MODELS], fontsize=6.5)
    ax1.set_ylim(0, 112)
    ax1.legend(frameon=False, ncol=3, fontsize=6.5)
    ax1.axhline(y=100, color="#ccc", linewidth=0.4, linestyle="--", zorder=0)

    # (b) Diversity
    ax2 = fig.add_subplot(gs[0, 1])
    for i, rep in enumerate(REPS):
        vals = []
        for model in MODELS:
            d = load("exp1_validity", f"{model}_{rep}_unconstrained_len6_summary.json")
            vals.append(d.get("uniqueness", 0) * 100 if d else 0)
        ax2.bar(x + offsets[i], vals, width * 0.9, label=rep.upper(),
                color=REP_COLORS[rep], edgecolor="white", linewidth=0.3, alpha=0.85)
    ax2.set_ylabel("Uniqueness (%)")
    ax2.set_title("(b) Sequence Diversity", fontweight="bold", loc="left")
    ax2.set_xticks(x)
    ax2.set_xticklabels([MODEL_SHORT[m] for m in MODELS], fontsize=6.5)
    ax2.set_ylim(0, 112)
    ax2.legend(frameon=False, ncol=3, fontsize=6.5)
    ax2.axhline(y=100, color="#ccc", linewidth=0.4, linestyle="--", zorder=0)

    # (c) Controllability
    ax3 = fig.add_subplot(gs[1, 0])
    ctrl_data = []
    ctrl_labels = []
    for model in MODELS:
        d = load("exp2_controllability", f"{model}_exp2_controllability_summary.json")
        if d and d.get("n_total_variants", 0) > 0:
            ctrl_data.append(d)
            ctrl_labels.append(MODEL_SHORT[model])
    cx = np.arange(len(ctrl_data))
    cw = 0.2
    metrics = [("edit_compliance", "Edit", "#2563EB"),
               ("frozen_compliance", "Frozen", "#16A34A"),
               ("id_valid", "ID Valid", "#F59E0B")]
    for i, (key, label, color) in enumerate(metrics):
        vals = [d.get(key, 0) * 100 for d in ctrl_data]
        ax3.bar(cx + (i-1)*cw, vals, cw*0.85, label=label,
                color=color, edgecolor="white", linewidth=0.3, alpha=0.85)
    ax3.set_ylabel("Rate (%)")
    ax3.set_title("(c) EDIT/FROZEN Controllability", fontweight="bold", loc="left")
    ax3.set_xticks(cx)
    ax3.set_xticklabels(ctrl_labels, fontsize=6.5)
    ax3.set_ylim(0, 112)
    ax3.legend(frameon=False, ncol=3, fontsize=6.5)
    ax3.axhline(y=100, color="#ccc", linewidth=0.4, linestyle="--", zorder=0)

    # (d) Constraints heatmap
    ax4 = fig.add_subplot(gs[1, 1])
    models_wd = []
    mat = []
    nv_mat = []
    for model in MODELS:
        row, nv = [], []
        has = False
        for profile in PROFILES:
            d = load("exp3_property_constrained",
                     f"{model}_restoken_{profile}_len6_summary.json")
            if d and d.get("n_valid", 0) > 0:
                row.append(d.get("constraint_satisfaction", 0) * 100)
                nv.append(d.get("n_valid", 0))
                has = True
            elif d and d.get("n_parsed", 0) > 0:
                row.append(0)
                nv.append(d.get("n_parsed", 0))
                has = True
            else:
                row.append(np.nan)
                nv.append(0)
        if has:
            models_wd.append(model)
            mat.append(row)
            nv_mat.append(nv)

    mat = np.array(mat)
    cmap = plt.cm.RdYlGn.copy()
    cmap.set_bad(color="#f5f5f5")
    im = ax4.imshow(mat, cmap=cmap, aspect="auto", vmin=0, vmax=100)
    ax4.set_xticks(range(len(PROFILES)))
    ax4.set_xticklabels([PROFILE_DISPLAY[p] for p in PROFILES], fontsize=7)
    ax4.set_yticks(range(len(models_wd)))
    ax4.set_yticklabels([MODEL_SHORT[m].replace("\n", " ") for m in models_wd], fontsize=7)
    for i in range(len(models_wd)):
        for j in range(len(PROFILES)):
            val = mat[i, j]
            nv = nv_mat[i][j]
            if np.isnan(val):
                ax4.text(j, i, "N/A", ha="center", va="center", fontsize=6, color="#aaa")
            else:
                c = "white" if val > 55 else ("black" if val > 0 else "#666")
                ax4.text(j, i, f"{val:.0f}%\n(n={nv})", ha="center", va="center",
                         fontsize=6, fontweight="bold", color=c)
    ax4.set_title("(d) Constraint Satisfaction", fontweight="bold", loc="left")
    fig.colorbar(im, ax=ax4, shrink=0.7, pad=0.06).set_label("Sat. %", fontsize=7)

    out = FIG_DIR / "w2_fig_combined_panel.pdf"
    fig.savefig(out)
    fig.savefig(out.with_suffix(".png"))
    plt.close(fig)
    print(f"  {out.name}")


# ─── Fig F: Radar Chart ──────────────────────────────────────────────────────

def fig_radar():
    categories = [
        "ResToken\nValidity",
        "SMILES\nValidity",
        "HELM\nValidity",
        "Diversity\n(ResToken)",
        "Edit\nCompliance",
        "Constraint\nSat.",
    ]
    N = len(categories)
    angles = np.linspace(0, 2 * np.pi, N, endpoint=False).tolist()
    angles += angles[:1]

    fig, ax = plt.subplots(figsize=(6, 6), subplot_kw=dict(polar=True))

    for model in MODELS:
        vals = []
        # ResToken validity
        d = load("exp1_validity", f"{model}_restoken_unconstrained_len6_summary.json")
        vals.append(d.get("n_valid", 0) / max(d.get("n_parsed", 1), 1) * 100 if d else 0)
        # SMILES validity
        d = load("exp1_validity", f"{model}_smiles_unconstrained_len6_summary.json")
        vals.append(d.get("n_valid", 0) / max(d.get("n_parsed", 1), 1) * 100 if d else 0)
        # HELM validity
        d = load("exp1_validity", f"{model}_helm_unconstrained_len6_summary.json")
        vals.append(d.get("n_valid", 0) / max(d.get("n_parsed", 1), 1) * 100 if d else 0)
        # Diversity (restoken uniqueness)
        d = load("exp1_validity", f"{model}_restoken_unconstrained_len6_summary.json")
        vals.append(d.get("uniqueness", 0) * 100 if d else 0)
        # Edit compliance
        d = load("exp2_controllability", f"{model}_exp2_controllability_summary.json")
        vals.append(d.get("edit_compliance", 0) * 100 if d else 0)
        # Avg constraint satisfaction across 3 profiles
        cs_vals = []
        for profile in PROFILES:
            d = load("exp3_property_constrained",
                     f"{model}_restoken_{profile}_len6_summary.json")
            if d and d.get("n_valid", 0) > 0:
                cs_vals.append(d.get("constraint_satisfaction", 0) * 100)
        vals.append(np.mean(cs_vals) if cs_vals else 0)

        vals += vals[:1]
        ax.plot(angles, vals, linewidth=1.5, label=MODEL_SHORT[model].replace("\n", " "),
                color=MODEL_COLORS[model], alpha=0.8)
        ax.fill(angles, vals, color=MODEL_COLORS[model], alpha=0.08)

    ax.set_xticks(angles[:-1])
    ax.set_xticklabels(categories, fontsize=7)
    ax.set_ylim(0, 105)
    ax.set_yticks([25, 50, 75, 100])
    ax.set_yticklabels(["25", "50", "75", "100"], fontsize=6, color="#888")
    ax.set_title("Overall Model Capability Profile", fontweight="bold", fontsize=10, pad=20)
    ax.legend(loc="upper right", bbox_to_anchor=(1.35, 1.1), frameon=False, fontsize=7)

    out = FIG_DIR / "w2_fig_f_radar.pdf"
    fig.savefig(out)
    fig.savefig(out.with_suffix(".png"))
    plt.close(fig)
    print(f"  {out.name}")


# ─── Main ─────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    print("Generating W2 benchmark figures...")
    fig_exp1_validity()
    fig_exp1_diversity()
    fig_exp2_controllability()
    fig_exp3_heatmap()
    fig_combined_panel()
    fig_radar()
    print(f"\nAll figures saved to {FIG_DIR}/")
