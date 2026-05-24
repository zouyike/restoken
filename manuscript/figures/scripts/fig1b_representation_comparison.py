#!/usr/bin/env python
"""
Figure 1b: Comparison of three peptide representations — SMILES, HELM, ResToken.

Demonstrates why ResToken is superior for LLM applications:
shorter, human-readable, and semantically encoded.

Output: manuscript/figures/generated/fig1b_representation_comparison.{png,pdf}
"""

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch

# ── ACS / JCIM style settings ──────────────────────────────────────────────
plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["DejaVu Sans"],
    "font.size": 7,
    "axes.linewidth": 0.6,
    "axes.labelsize": 7,
    "figure.dpi": 300,
    "savefig.dpi": 300,
    "pdf.fonttype": 42,       # TrueType fonts in PDF (ACS requirement)
    "ps.fonttype": 42,
})

# ── Peptide data ────────────────────────────────────────────────────────────
residues = [
    {"token": "A01", "backbone": "alpha", "chirality": "L", "bulk": "large",
     "smiles": "CCC(C)(C)[C@H](NC(C)=O)C(=O)NC"},
    {"token": "N15", "backbone": "beta",  "chirality": "D", "bulk": "small",
     "smiles": "CNC(=O)C[C@@H](C)NC(C)=O"},
    {"token": "s01", "backbone": "gamma", "chirality": "L", "bulk": "small",
     "smiles": "CNC(=O)C[C@@H]1CCN(C(C)=O)C1"},
    {"token": "K08", "backbone": "beta",  "chirality": "L", "bulk": "large",
     "smiles": "CNC(=O)C[C@H](CCCC[NH3+])NC(C)=O"},
    {"token": "E09", "backbone": "alpha", "chirality": "D", "bulk": "large",
     "smiles": "CNC(=O)[C@@H](CCC(C)C)NC(C)=O"},
    {"token": "U04", "backbone": "gamma", "chirality": "D", "bulk": "large",
     "smiles": "CNC(=O)[C@@H]1C[C@H]1[C@@H](NC(C)=O)c1ccccc1"},
]

# Full-peptide SMILES (comma-separated per-residue, typical LLM input)
smiles_str = ",".join(r["smiles"] for r in residues)

# HELM
helm_str = "PEPTIDE1{[A01].[N15].[s01].[K08].[E09].[U04]}$$$$V2.0"
helm_str_display = helm_str.replace("$", r"\$")

# ResToken
restoken_str = "A01-N15-s01-K08-E09-U04"

# Colors
COLORS = {
    "alpha": "#2B6CB0",   # blue
    "beta":  "#38A169",   # green
    "gamma": "#DD6B20",   # orange
}
BG_PANEL = "#F8F8F8"
BORDER   = "#D0D0D0"
TEXT_DIM = "#777777"
TEXT_DARK = "#1A1A1A"
SMILES_RED = "#A0522D"    # sienna - readable dark red-brown

# ── Figure layout ───────────────────────────────────────────────────────────
fig_w = 3.25   # single-column JCIM
fig_h = 4.2

fig = plt.figure(figsize=(fig_w, fig_h))

# Use gridspec for precise control
gs = fig.add_gridspec(3, 1, height_ratios=[1.0, 0.55, 1.7],
                      hspace=0.12, left=0.02, right=0.98,
                      top=0.97, bottom=0.02)

axes = [fig.add_subplot(gs[i]) for i in range(3)]

for ax in axes:
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")


def draw_panel_bg(ax, color=BG_PANEL, edge=BORDER):
    """Draw a rounded background spanning the full axes."""
    rect = FancyBboxPatch(
        (0.0, 0.0), 1.0, 1.0,
        boxstyle="round,pad=0.02",
        facecolor=color, edgecolor=edge, linewidth=0.5,
        transform=ax.transAxes, zorder=0,
        clip_on=False
    )
    ax.add_patch(rect)


# ═══════════════════════════════════════════════════════════════════════════
# Panel (i): SMILES
# ═══════════════════════════════════════════════════════════════════════════
ax = axes[0]
draw_panel_bg(ax)

ax.text(0.04, 0.90, "(i)", fontsize=8, fontweight="bold", color=TEXT_DARK,
        va="top", transform=ax.transAxes)
ax.text(0.12, 0.90, "SMILES", fontsize=7.5, fontweight="bold", color=TEXT_DARK,
        va="top", transform=ax.transAxes)

# Wrap SMILES into 3 lines to fit within panel
chars_per_line = 55
lines = []
s = smiles_str
while s:
    if len(s) <= chars_per_line:
        lines.append(s)
        break
    # Find a comma near the target length to split cleanly
    cut = s.rfind(",", 0, chars_per_line + 1)
    if cut == -1:
        cut = chars_per_line
    lines.append(s[:cut + 1])
    s = s[cut + 1:]

y_text = 0.72
for line in lines:
    ax.text(0.04, y_text, line, fontsize=4.2, fontfamily="monospace",
            color=SMILES_RED, va="top", transform=ax.transAxes)
    y_text -= 0.17

# Annotation
ax.text(0.04, 0.08,
        f"{len(smiles_str)} characters  |  No semantic encoding  |  Not human-readable",
        fontsize=5.0, color=TEXT_DIM, va="bottom", style="italic",
        transform=ax.transAxes)


# ═══════════════════════════════════════════════════════════════════════════
# Panel (ii): HELM
# ═══════════════════════════════════════════════════════════════════════════
ax = axes[1]
draw_panel_bg(ax)

ax.text(0.04, 0.88, "(ii)", fontsize=8, fontweight="bold", color=TEXT_DARK,
        va="top", transform=ax.transAxes)
ax.text(0.13, 0.88, "HELM", fontsize=7.5, fontweight="bold", color=TEXT_DARK,
        va="top", transform=ax.transAxes)

ax.text(0.04, 0.48, helm_str_display, fontsize=4.8, fontfamily="monospace",
        color="#2E7D32", va="center", transform=ax.transAxes)

ax.text(0.04, 0.10,
        f"{len(helm_str)} characters  |  Monomer IDs only  |  Requires external library",
        fontsize=5.0, color=TEXT_DIM, va="bottom", style="italic",
        transform=ax.transAxes)


# ═══════════════════════════════════════════════════════════════════════════
# Panel (iii): ResToken (the hero panel)
# ═══════════════════════════════════════════════════════════════════════════
ax = axes[2]
draw_panel_bg(ax)

ax.text(0.04, 0.97, "(iii)", fontsize=8, fontweight="bold", color=TEXT_DARK,
        va="top", transform=ax.transAxes)
ax.text(0.13, 0.97, "ResToken", fontsize=7.5, fontweight="bold", color=TEXT_DARK,
        va="top", transform=ax.transAxes)
ax.text(0.30, 0.97, "(ours)", fontsize=6.5, color=TEXT_DIM,
        va="top", style="italic", transform=ax.transAxes)

# ── Draw color-coded tokens ─────────────────────────────────────────────
token_y = 0.82
token_w = 0.105
token_h = 0.075
gap = 0.04        # space between tokens (includes hyphen)
x_start = 0.06

x = x_start
for i, res in enumerate(residues):
    color = COLORS[res["backbone"]]

    # Colored pill background
    pill = FancyBboxPatch(
        (x - 0.005, token_y - token_h / 2 - 0.005),
        token_w + 0.01, token_h + 0.01,
        boxstyle="round,pad=0.008",
        facecolor=color, edgecolor="none", alpha=0.13,
        transform=ax.transAxes, zorder=1
    )
    ax.add_patch(pill)

    # Token text
    ax.text(x + token_w / 2, token_y, res["token"],
            fontsize=7.5, fontfamily="monospace", fontweight="bold",
            color=color, ha="center", va="center",
            transform=ax.transAxes, zorder=2)

    # Hyphen separator
    if i < len(residues) - 1:
        ax.text(x + token_w + gap / 2, token_y, "-",
                fontsize=7, fontfamily="monospace", color="#999999",
                ha="center", va="center", transform=ax.transAxes)

    x += token_w + gap

# ── Property annotations below tokens ──────────────────────────────────
annot_top = 0.64
line_spacing = 0.055

x = x_start
for i, res in enumerate(residues):
    color = COLORS[res["backbone"]]
    cx = x + token_w / 2

    # Backbone type (bold, colored)
    ax.text(cx, annot_top, res["backbone"],
            fontsize=5, color=color, ha="center", va="top",
            fontweight="bold", transform=ax.transAxes)
    # Chirality
    ax.text(cx, annot_top - line_spacing, f"{res['chirality']}-config",
            fontsize=4.5, color=TEXT_DIM, ha="center", va="top",
            transform=ax.transAxes)
    # Bulk
    ax.text(cx, annot_top - 2 * line_spacing, res["bulk"],
            fontsize=4.5, color=TEXT_DIM, ha="center", va="top",
            transform=ax.transAxes)

    x += token_w + gap

# ── Encoding legend ─────────────────────────────────────────────────────
legend_top = 0.38
lx = 0.06

ax.text(lx, legend_top, "Token encoding scheme:",
        fontsize=5.5, fontweight="bold", color=TEXT_DARK,
        va="top", transform=ax.transAxes)

# Divider line
ax.plot([lx, 0.94], [legend_top - 0.03, legend_top - 0.03],
        color=BORDER, linewidth=0.4, transform=ax.transAxes, zorder=1)

# Left column: prefix letter rules
col1_x = lx + 0.01
col1_top = legend_top - 0.06

rules_l = [
    ("Prefix letter  →  backbone type", TEXT_DARK, "bold"),
    ("  A–E = alpha", COLORS["alpha"], "normal"),
    ("  K–O = beta", COLORS["beta"], "normal"),
    ("  U–Y = gamma", COLORS["gamma"], "normal"),
]

y = col1_top
for text, col, weight in rules_l:
    ax.text(col1_x, y, text, fontsize=4.5, fontfamily="monospace",
            color=col, va="top", fontweight=weight, transform=ax.transAxes)
    y -= 0.045

# Right column: case + suffix
col2_x = 0.50
y = col1_top

rules_r = [
    ("Case  →  chirality", TEXT_DARK, "bold"),
    ("  UPPER (A01) = L-config", TEXT_DARK, "normal"),
    ("  lower (s01) = D-config", TEXT_DARK, "normal"),
    ("Suffix digits  →  side-chain ID", TEXT_DARK, "bold"),
]

for text, col, weight in rules_r:
    ax.text(col2_x, y, text, fontsize=4.5, fontfamily="monospace",
            color=col, va="top", fontweight=weight, transform=ax.transAxes)
    y -= 0.045

# Backbone color dots legend (far right)
dot_x = 0.88
dot_y = col1_top + 0.01
for btype, col in COLORS.items():
    circle = mpatches.Circle((dot_x, dot_y), 0.009, color=col,
                             transform=ax.transAxes, zorder=3,
                             clip_on=False)
    ax.add_patch(circle)
    ax.text(dot_x + 0.025, dot_y, btype, fontsize=4.5, color=col,
            va="center", fontweight="bold", transform=ax.transAxes)
    dot_y -= 0.05

# Bottom annotation
ax.text(0.04, 0.02,
        f"{len(restoken_str)} characters  |  Semantically encoded  |  Human-readable",
        fontsize=5.5, color=TEXT_DARK, va="bottom", style="italic",
        fontweight="bold", transform=ax.transAxes)

# ── Save ────────────────────────────────────────────────────────────────
out_base = "manuscript/figures/generated/fig1b_representation_comparison"
fig.savefig(f"{out_base}.png", dpi=300, bbox_inches="tight",
            facecolor="white", pad_inches=0.06)
fig.savefig(f"{out_base}.pdf", bbox_inches="tight",
            facecolor="white", pad_inches=0.06)
plt.close(fig)

print(f"Saved: {out_base}.png")
print(f"Saved: {out_base}.pdf")
print(f"SMILES length: {len(smiles_str)}")
print(f"HELM length:   {len(helm_str)}")
print(f"ResToken length: {len(restoken_str)}")
