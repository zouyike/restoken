#!/usr/bin/env python
"""
Figure 2: Benchmark Results for ResToken manuscript (JCIM).
Three-panel figure: (a) Exp1 Validity, (b) Exp2 Controllability Heatmap,
(c) Exp3 Constraint Satisfaction.
"""

import json
import pathlib

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np

# ── paths ──────────────────────────────────────────────────────────────
BASE = pathlib.Path("/scratch/genesis/NCAA_tokenization")
SUMMARY = BASE / "experiments/output/w3_statistical_summary.json"
OUT_DIR = BASE / "manuscript/figures/generated"
OUT_DIR.mkdir(parents=True, exist_ok=True)

# ── load data ──────────────────────────────────────────────────────────
with open(SUMMARY) as f:
    data = json.load(f)

# ── model name mapping ─────────────────────────────────────────────────
MODEL_ORDER = [
    "gemini-2.5-pro",
    "gemini-2.5-flash",
    "gpt-4o",
    "Qwen3.5-9B",
    "gemma-3-12b-it-bnb-4bit",
]
DISPLAY_NAMES = {
    "gemini-2.5-pro": "Gemini Pro",
    "gemini-2.5-flash": "Gemini Flash",
    "gpt-4o": "GPT-4o",
    "Qwen3.5-9B": "Qwen 3.5",
    "gemma-3-12b-it-bnb-4bit": "Gemma 3",
}

# ── colors ─────────────────────────────────────────────────────────────
REP_COLORS = {
    "restoken": "#2196F3",
    "smiles": "#FF9800",
    "helm": "#4CAF50",
}
REP_LABELS = {
    "restoken": "ResToken",
    "smiles": "SMILES",
    "helm": "HELM",
}
REPS = ["restoken", "smiles", "helm"]

# Qualitative palette for per-model bars in panel (c)
MODEL_COLORS = {
    "gemini-2.5-pro": "#2196F3",
    "gemini-2.5-flash": "#FF9800",
    "gpt-4o": "#4CAF50",
    "Qwen3.5-9B": "#E91E63",
    "gemma-3-12b-it-bnb-4bit": "#9C27B0",
}

PROFILE_DISPLAY = {
    "permeable": "Permeable",
    "charged_binder": "Charged\nBinder",
    "rigid_scaffold": "Rigid\nScaffold",
}
PROFILES = ["permeable", "charged_binder", "rigid_scaffold"]

# ── style ──────────────────────────────────────────────────────────────
plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 8,
    "axes.labelsize": 8,
    "axes.titlesize": 8,
    "xtick.labelsize": 7,
    "ytick.labelsize": 7,
    "legend.fontsize": 7,
    "figure.dpi": 300,
    "savefig.dpi": 300,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "figure.facecolor": "white",
    "axes.facecolor": "white",
})

# ── figure setup ───────────────────────────────────────────────────────
fig = plt.figure(figsize=(7, 2.8))
gs = fig.add_gridspec(1, 3, width_ratios=[3, 1.3, 2.5], wspace=0.45)
ax_a = fig.add_subplot(gs[0])
ax_b = fig.add_subplot(gs[1])
ax_c = fig.add_subplot(gs[2])

# ═══════════════════════════════════════════════════════════════════════
# Panel (a): Exp1 Validity — grouped bar chart
# ═══════════════════════════════════════════════════════════════════════
exp1 = data["exp1_bootstrap_ci"]
n_models = len(MODEL_ORDER)
n_reps = len(REPS)
bar_width = 0.22
x = np.arange(n_models)

for j, rep in enumerate(REPS):
    vals, lo, hi = [], [], []
    for model in MODEL_ORDER:
        entry = exp1.get(model, {}).get(rep)
        if entry:
            pt = entry["point"] * 100
            vals.append(pt)
            lo.append(pt - entry["ci_lower"] * 100)
            hi.append(entry["ci_upper"] * 100 - pt)
        else:
            vals.append(0)
            lo.append(0)
            hi.append(0)
    offset = (j - 1) * bar_width
    ax_a.bar(x + offset, vals, bar_width,
             yerr=[lo, hi], capsize=2, linewidth=0.5,
             color=REP_COLORS[rep], label=REP_LABELS[rep],
             edgecolor="white", error_kw={"linewidth": 0.7})

ax_a.set_xticks(x)
ax_a.set_xticklabels([DISPLAY_NAMES[m] for m in MODEL_ORDER], rotation=25,
                      ha="right")
ax_a.set_ylabel("Validity Rate (%)")
ax_a.set_ylim(0, 110)
ax_a.yaxis.set_major_locator(mticker.MultipleLocator(20))
ax_a.legend(loc="upper right", frameon=False, ncol=1)
ax_a.set_title("(a) Validity", fontweight="bold", loc="left")

# ═══════════════════════════════════════════════════════════════════════
# Panel (b): Exp2 Controllability — heatmap
# ═══════════════════════════════════════════════════════════════════════
exp2 = data["exp2_controllability_per_model"]
metrics = ["edit_compliance", "frozen_compliance"]
metric_labels = ["Edit\nCompl.", "Frozen\nCompl."]
heatmap_data = np.full((n_models, 2), np.nan)
for i, model in enumerate(MODEL_ORDER):
    entry = exp2.get(model)
    if entry:
        heatmap_data[i, 0] = entry["edit_compliance"] * 100
        heatmap_data[i, 1] = entry["frozen_compliance"] * 100

im = ax_b.imshow(heatmap_data, cmap="Blues", aspect="auto",
                 vmin=0, vmax=100)
# Annotate cells
for i in range(n_models):
    for j in range(2):
        val = heatmap_data[i, j]
        if not np.isnan(val):
            text_color = "white" if val > 55 else "black"
            ax_b.text(j, i, f"{val:.0f}%", ha="center", va="center",
                      fontsize=7, color=text_color, fontweight="bold")

ax_b.set_xticks(np.arange(2))
ax_b.set_xticklabels(metric_labels, fontsize=7)
ax_b.set_yticks(np.arange(n_models))
ax_b.set_yticklabels([DISPLAY_NAMES[m] for m in MODEL_ORDER], fontsize=7)
ax_b.set_title("(b) Controllability", fontweight="bold", loc="left")
# Remove spines for heatmap — they look odd; use thin border instead
for spine in ax_b.spines.values():
    spine.set_visible(False)
ax_b.tick_params(length=0)

# ═══════════════════════════════════════════════════════════════════════
# Panel (c): Exp3 Constraint Satisfaction — grouped bar chart
# ═══════════════════════════════════════════════════════════════════════
enrichment = data["enrichment"]

# Determine which models have any n_valid > 0 across profiles (ResToken only)
models_c = []
for model in MODEL_ORDER:
    has_valid = False
    for profile in PROFILES:
        key = f"{model}|{profile}"
        entry = enrichment.get(key, {})
        if entry.get("n_valid", 0) > 0:
            has_valid = True
            break
    if has_valid:
        models_c.append(model)

n_profiles = len(PROFILES)
n_models_c = len(models_c)
bar_width_c = 0.8 / max(n_models_c, 1)
x_c = np.arange(n_profiles)

for j, model in enumerate(models_c):
    vals = []
    for profile in PROFILES:
        key = f"{model}|{profile}"
        entry = enrichment.get(key, {})
        nv = entry.get("n_valid", 0)
        if nv > 0:
            vals.append(entry.get("observed", 0) * 100)
        else:
            vals.append(0)
    offset = (j - (n_models_c - 1) / 2) * bar_width_c
    ax_c.bar(x_c + offset, vals, bar_width_c,
             color=MODEL_COLORS[model],
             label=DISPLAY_NAMES[model],
             edgecolor="white", linewidth=0.5)

ax_c.set_xticks(x_c)
ax_c.set_xticklabels([PROFILE_DISPLAY[p] for p in PROFILES])
ax_c.set_ylabel("Constraint Satisfaction (%)")
ax_c.set_ylim(0, 115)
ax_c.yaxis.set_major_locator(mticker.MultipleLocator(20))
ax_c.legend(loc="upper right", frameon=False, fontsize=6, ncol=1)
ax_c.set_title("(c) Constraint Satisfaction", fontweight="bold", loc="left")

# ── save ───────────────────────────────────────────────────────────────
fig.subplots_adjust(left=0.07, right=0.98, bottom=0.22, top=0.90)
for fmt in ("png", "pdf"):
    out = OUT_DIR / f"fig2_benchmark_results.{fmt}"
    fig.savefig(out, dpi=300, bbox_inches="tight", facecolor="white")
    print(f"Saved: {out}")

plt.close(fig)
print("Done.")
