"""Generate publication-quality figures for the ResToken JCIM paper.

Usage:
    python plot_figures.py --fig all        # generate all available figures
    python plot_figures.py --fig tokenization
    python plot_figures.py --fig treemap
    python plot_figures.py --fig baseline
    python plot_figures.py --fig comparison
"""

import argparse
import csv
import json
import sys
from pathlib import Path
from collections import Counter

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.gridspec import GridSpec
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from restoken.src.library import BlockLibrary

FIG_DIR = Path(__file__).resolve().parents[2] / "manuscript" / "figures" / "generated"
BASELINE_DIR = Path(__file__).resolve().parents[2] / "experiments" / "output" / "random_baseline"

# ACS style
plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 8,
    "axes.titlesize": 9,
    "axes.labelsize": 8,
    "xtick.labelsize": 7,
    "ytick.labelsize": 7,
    "legend.fontsize": 7,
    "figure.dpi": 300,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
    "axes.linewidth": 0.6,
    "xtick.major.width": 0.6,
    "ytick.major.width": 0.6,
})

COLORS = {
    "alpha": "#4C72B0",
    "beta": "#55A868",
    "gamma": "#C44E52",
    "restoken": "#4C72B0",
    "smiles": "#C44E52",
    "helm": "#8172B2",
    "random": "#999999",
}


def fig_tokenization(lib: BlockLibrary):
    """Figure 1a: Three example NCAAs with 2D structures and property cards."""
    from rdkit import Chem
    from rdkit.Chem import Draw, rdDepictor
    from rdkit.Chem.Draw import rdMolDraw2D
    from io import BytesIO
    from PIL import Image

    rdDepictor.SetPreferCoordGen(True)

    examples = [
        ("A01", "Alpha-L, Isoleucine analog"),
        ("N03", "Beta-D, small bulk"),
        ("S05", "Gamma-L, Ile-like, NCY-free"),
    ]

    fig, axes = plt.subplots(1, 3, figsize=(7.5, 3.0))

    for ax, (bid, desc) in zip(axes, examples):
        block = lib[bid]
        mol = Chem.MolFromSmiles(block.aa_smiles)
        if mol is None:
            ax.text(0.5, 0.5, f"{bid}: invalid SMILES", ha="center", va="center")
            continue

        rdDepictor.Compute2DCoords(mol)

        drawer = rdMolDraw2D.MolDraw2DCairo(400, 280)
        opts = drawer.drawOptions()
        opts.bondLineWidth = 1.5
        opts.padding = 0.15
        drawer.DrawMolecule(mol)
        drawer.FinishDrawing()
        png = drawer.GetDrawingText()

        img = Image.open(BytesIO(png))
        ax.imshow(img, aspect="auto")
        ax.axis("off")

        props_text = (
            f"{bid}  [{block.aa_class}]  {block.chirality}\n"
            f"{block.mc_type} | {block.mc_nmod} | {block.charge_label}\n"
            f"bulk={block.sc_bulk} | HBD={block.sc_hbd} | HBA={block.sc_hba} | rot={block.rot_total}"
        )

        ax.set_title(f"{bid}: {desc}", fontsize=8, fontweight="bold", pad=3)
        ax.text(
            0.5, -0.02, props_text, transform=ax.transAxes,
            ha="center", va="top", fontsize=6.5,
            bbox=dict(boxstyle="round,pad=0.3", facecolor="#f0f0f0", edgecolor="#cccccc", alpha=0.9),
            family="monospace",
        )

    plt.subplots_adjust(wspace=0.05)
    out = FIG_DIR / "fig1a_tokenization_scheme.png"
    fig.savefig(out, dpi=300)
    plt.close()
    print(f"Saved: {out}")
    return out


def fig_treemap(lib: BlockLibrary):
    """Figure 5: Library coverage — backbone × class × chirality as a grouped bar chart."""
    classes_order = ["I", "F", "A", "T", "P", "V", "S", "H", "K", "L", "M", "W",
                     "D", "Y", "E", "G", "C", "N", "Q", "R", "ncaa"]

    data = {bt: [] for bt in ["alpha", "beta", "gamma"]}
    for cls in classes_order:
        for bt in ["alpha", "beta", "gamma"]:
            count = len(lib.blocks_by_property(mc_type=bt, aa_class=cls))
            data[bt].append(count)

    x = np.arange(len(classes_order))
    width = 0.25

    fig, ax = plt.subplots(figsize=(7.5, 3.0))

    for i, (bt, color) in enumerate([("alpha", COLORS["alpha"]),
                                       ("beta", COLORS["beta"]),
                                       ("gamma", COLORS["gamma"])]):
        offset = (i - 1) * width
        bars = ax.bar(x + offset, data[bt], width, label=f"{bt} ({sum(data[bt])})",
                      color=color, edgecolor="white", linewidth=0.3)
        for bar, val in zip(bars, data[bt]):
            if val > 0:
                ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.3,
                        str(val), ha="center", va="bottom", fontsize=5)

    ax.set_xticks(x)
    ax.set_xticklabels(classes_order, fontsize=7)
    ax.set_xlabel("Functional Class (canonical AA analog)")
    ax.set_ylabel("Number of Building Blocks")
    ax.set_title("ResToken v11 Library: 400 Blocks by Backbone Type × Functional Class")
    ax.legend(loc="upper right", framealpha=0.9)
    ax.set_ylim(0, max(max(v) for v in data.values()) + 8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    out = FIG_DIR / "fig5_library_coverage.png"
    fig.savefig(out, dpi=300)
    plt.close()
    print(f"Saved: {out}")
    return out


def fig_baseline_distributions():
    """Baseline property distributions across 4 constraint profiles."""
    profiles = ["unconstrained", "permeable", "charged_binder", "rigid_scaffold"]
    profile_labels = ["Unconstrained\n(100%)", "Permeable\n(2.07%)",
                      "Charged Binder\n(1.82%)", "Rigid Scaffold\n(0.81%)"]

    fig, axes = plt.subplots(2, 2, figsize=(7.5, 5.5))
    props = [
        ("net_charge", "Net Charge"),
        ("total_hbd", "Total HBD"),
        ("total_hba", "Total HBA"),
        ("total_rot", "Total Rotatable Bonds"),
    ]

    for ax_idx, (prop, ylabel) in enumerate(props):
        ax = axes[ax_idx // 2, ax_idx % 2]
        all_data = []
        all_labels = []

        for i, (profile, label) in enumerate(zip(profiles, profile_labels)):
            csv_files = sorted(BASELINE_DIR.glob(f"baseline_{profile}_len6_*.csv"))
            if not csv_files:
                continue
            csv_file = csv_files[0]
            vals = []
            with open(csv_file) as f:
                reader = csv.DictReader(f)
                for row in reader:
                    vals.append(float(row[prop]))
            all_data.append(vals)
            all_labels.append(label)

        if not all_data:
            ax.text(0.5, 0.5, "No data", ha="center", va="center")
            continue

        colors_list = [COLORS["alpha"], COLORS["beta"], COLORS["gamma"], "#DD8452"]
        bp = ax.boxplot(all_data, labels=all_labels, patch_artist=True,
                        widths=0.6, showfliers=False, medianprops=dict(color="black", linewidth=1))
        for patch, color in zip(bp["boxes"], colors_list):
            patch.set_facecolor(color)
            patch.set_alpha(0.7)

        ax.set_ylabel(ylabel)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

        for i, vals in enumerate(all_data):
            mean_val = np.mean(vals)
            ax.plot(i + 1, mean_val, "D", color="white", markersize=3, markeredgecolor="black", markeredgewidth=0.5)

    fig.suptitle("Random Baseline Property Distributions (6-mer, N=1000 each)", fontsize=9, fontweight="bold")
    plt.tight_layout(rect=[0, 0, 1, 0.96])

    out = FIG_DIR / "fig_baseline_distributions.png"
    fig.savefig(out, dpi=300)
    plt.close()
    print(f"Saved: {out}")
    return out


def fig_comparison(lib: BlockLibrary):
    """Figure 1b: Same peptide in SMILES vs HELM vs ResToken."""
    seq = ["A01", "K03", "N12", "S05", "E02", "a07"]
    blocks = [lib[tid] for tid in seq]

    smiles_parts = [b.aa_smiles for b in blocks]
    smiles_full = ".".join(smiles_parts)

    helm_raw = "PEPTIDE1{" + ".".join(f"[{tid}]" for tid in seq) + "}$PEPTIDE1,PEPTIDE1,1:R1-6:R2$$$V2.0"
    helm_str = helm_raw.replace("$", r"\$")

    restoken_str = "-".join(seq)
    restoken_props = "  ".join(
        f"{b.aa_class}/{b.chirality}/{b.charge_label}" for b in blocks
    )
    restoken_bb = "  ".join(f"{b.mc_type:>5s}" for b in blocks)

    fig, ax = plt.subplots(figsize=(7.5, 3.5))
    ax.axis("off")

    y = 0.92
    dy = 0.12

    ax.text(0.02, y, "SMILES", fontsize=9, fontweight="bold", color=COLORS["smiles"],
            transform=ax.transAxes, va="top")
    # Truncate SMILES to show it's long and unreadable
    smiles_display = smiles_full[:90] + "..." if len(smiles_full) > 90 else smiles_full
    ax.text(0.02, y - 0.04, smiles_display, fontsize=5.5, family="monospace",
            transform=ax.transAxes, va="top", color="#555555",
            bbox=dict(boxstyle="round,pad=0.3", facecolor="#fff0f0", edgecolor=COLORS["smiles"], alpha=0.3))
    ax.text(0.98, y - 0.02, f"{len(smiles_full)} chars", fontsize=6, ha="right",
            transform=ax.transAxes, color=COLORS["smiles"], fontstyle="italic")

    y -= dy + 0.06
    ax.text(0.02, y, "HELM", fontsize=9, fontweight="bold", color=COLORS["helm"],
            transform=ax.transAxes, va="top")
    ax.text(0.02, y - 0.04, helm_str, fontsize=6, family="monospace",
            transform=ax.transAxes, va="top", color="#555555",
            bbox=dict(boxstyle="round,pad=0.3", facecolor="#f0f0ff", edgecolor=COLORS["helm"], alpha=0.3))
    ax.text(0.98, y - 0.02, f"{len(helm_raw)} chars", fontsize=6, ha="right",
            transform=ax.transAxes, color=COLORS["helm"], fontstyle="italic")

    y -= dy + 0.06
    ax.text(0.02, y, "ResToken", fontsize=9, fontweight="bold", color=COLORS["restoken"],
            transform=ax.transAxes, va="top")
    ax.text(0.02, y - 0.04, restoken_str, fontsize=8, family="monospace", fontweight="bold",
            transform=ax.transAxes, va="top", color="#333333",
            bbox=dict(boxstyle="round,pad=0.3", facecolor="#f0f5ff", edgecolor=COLORS["restoken"], alpha=0.3))
    ax.text(0.98, y - 0.02, f"{len(restoken_str)} chars", fontsize=6, ha="right",
            transform=ax.transAxes, color=COLORS["restoken"], fontstyle="italic")

    y -= 0.05
    ax.text(0.02, y - 0.04, f"Class:     {restoken_props}", fontsize=6, family="monospace",
            transform=ax.transAxes, va="top", color="#666666")
    ax.text(0.02, y - 0.09, f"Backbone: {restoken_bb}", fontsize=6, family="monospace",
            transform=ax.transAxes, va="top", color="#666666")

    charge_str = "  ".join(f"{b.charge:+d}" if b.charge != 0 else " 0" for b in blocks)
    net_chg = sum(b.charge for b in blocks)
    ax.text(0.02, y - 0.14, f"Charge:   {charge_str}     net = {net_chg:+d}",
            fontsize=6, family="monospace", transform=ax.transAxes, va="top", color="#666666")

    ax.set_title("Same 6-mer NCAA cyclic peptide in three representations",
                 fontsize=9, fontweight="bold", pad=10)

    out = FIG_DIR / "fig1b_representation_comparison.png"
    fig.savefig(out, dpi=300)
    plt.close()
    print(f"Saved: {out}")
    return out


def fig_chirality_charge(lib: BlockLibrary):
    """Chirality and charge distribution stacked bars by backbone type."""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7.5, 2.8))

    bts = ["alpha", "beta", "gamma"]
    bt_labels = [f"Alpha\n(n=250)", f"Beta\n(n=101)", f"Gamma\n(n=49)"]

    # Chirality
    chiral_data = {}
    for ch in ["L", "D", "A"]:
        chiral_data[ch] = [len(lib.blocks_by_property(mc_type=bt, chirality=ch)) for bt in bts]

    bottom = np.zeros(3)
    ch_colors = {"L": "#4C72B0", "D": "#C44E52", "A": "#999999"}
    ch_labels = {"L": "L-config", "D": "D-config", "A": "Achiral"}
    for ch in ["L", "D", "A"]:
        ax1.bar(bt_labels, chiral_data[ch], bottom=bottom, label=ch_labels[ch],
                color=ch_colors[ch], edgecolor="white", linewidth=0.5)
        for i, v in enumerate(chiral_data[ch]):
            if v > 5:
                ax1.text(i, bottom[i] + v / 2, str(v), ha="center", va="center", fontsize=6, color="white")
        bottom += chiral_data[ch]

    ax1.set_ylabel("Number of Blocks")
    ax1.set_title("Chirality Distribution")
    ax1.legend(loc="upper right", fontsize=6)
    ax1.spines["top"].set_visible(False)
    ax1.spines["right"].set_visible(False)

    # Charge
    charge_data = {}
    for cl, label in [("neu", "Neutral"), ("pos", "Positive"), ("neg", "Negative")]:
        charge_data[label] = [len(lib.blocks_by_property(mc_type=bt, charge_label=cl)) for bt in bts]

    bottom = np.zeros(3)
    cl_colors = {"Neutral": "#999999", "Positive": "#4C72B0", "Negative": "#C44E52"}
    for label in ["Neutral", "Positive", "Negative"]:
        ax2.bar(bt_labels, charge_data[label], bottom=bottom, label=label,
                color=cl_colors[label], edgecolor="white", linewidth=0.5)
        for i, v in enumerate(charge_data[label]):
            if v > 5:
                ax2.text(i, bottom[i] + v / 2, str(v), ha="center", va="center", fontsize=6, color="white")
        bottom += charge_data[label]

    ax2.set_ylabel("Number of Blocks")
    ax2.set_title("Charge Distribution")
    ax2.legend(loc="upper right", fontsize=6)
    ax2.spines["top"].set_visible(False)
    ax2.spines["right"].set_visible(False)

    plt.tight_layout()
    out = FIG_DIR / "fig_chirality_charge_dist.png"
    fig.savefig(out, dpi=300)
    plt.close()
    print(f"Saved: {out}")
    return out


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--fig", default="all",
                        choices=["all", "tokenization", "treemap", "baseline",
                                 "comparison", "chirality"])
    args = parser.parse_args()

    FIG_DIR.mkdir(parents=True, exist_ok=True)
    lib = BlockLibrary()

    figs = args.fig
    generated = []

    if figs in ("all", "tokenization"):
        generated.append(fig_tokenization(lib))

    if figs in ("all", "treemap"):
        generated.append(fig_treemap(lib))

    if figs in ("all", "baseline"):
        generated.append(fig_baseline_distributions())

    if figs in ("all", "comparison"):
        generated.append(fig_comparison(lib))

    if figs in ("all", "chirality"):
        generated.append(fig_chirality_charge(lib))

    print(f"\nGenerated {len(generated)} figures in {FIG_DIR}")


if __name__ == "__main__":
    main()
