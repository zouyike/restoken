"""Additional diverse figures for the ResToken JCIM paper.

Usage:
    python plot_figures_extra.py --fig all
    python plot_figures_extra.py --fig heatmap
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
import matplotlib.colors as mcolors
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from restoken.src.library import BlockLibrary

FIG_DIR = Path(__file__).resolve().parents[2] / "manuscript" / "figures" / "generated"
BASELINE_DIR = Path(__file__).resolve().parents[2] / "experiments" / "output" / "random_baseline"

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
})


def fig_heatmap(lib: BlockLibrary):
    """Heatmap: Functional Class × Backbone type coverage."""
    classes = ["I", "F", "A", "T", "P", "V", "S", "H", "K", "L", "M",
               "W", "D", "Y", "E", "G", "C", "N", "Q", "R", "ncaa"]
    bts = ["alpha", "beta", "gamma"]

    matrix = np.zeros((len(classes), len(bts)))
    for i, cls in enumerate(classes):
        for j, bt in enumerate(bts):
            matrix[i, j] = len(lib.blocks_by_property(mc_type=bt, aa_class=cls))

    fig, ax = plt.subplots(figsize=(3.5, 6.0))

    cmap = plt.cm.Blues.copy()
    cmap.set_under("white")
    im = ax.imshow(matrix, cmap=cmap, aspect="auto", vmin=0.5, vmax=matrix.max(),
                   interpolation="nearest")

    ax.set_xticks(range(len(bts)))
    ax.set_xticklabels([f"{bt}\n({int(matrix[:,j].sum())})" for j, bt in enumerate(bts)])
    ax.set_yticks(range(len(classes)))
    ax.set_yticklabels(classes)
    ax.set_xlabel("Backbone Type")
    ax.set_ylabel("Functional Class")

    for i in range(len(classes)):
        for j in range(len(bts)):
            val = int(matrix[i, j])
            if val == 0:
                ax.text(j, i, "—", ha="center", va="center", fontsize=6, color="#cccccc")
            else:
                color = "white" if val > matrix.max() * 0.6 else "black"
                ax.text(j, i, str(val), ha="center", va="center", fontsize=6,
                        fontweight="bold", color=color)

    ax.set_title("Building Block Coverage\n(Class × Backbone)", fontsize=9, fontweight="bold")

    cbar = fig.colorbar(im, ax=ax, shrink=0.6, pad=0.08)
    cbar.set_label("Number of blocks", fontsize=7)

    plt.tight_layout()
    out = FIG_DIR / "fig_heatmap_class_backbone.png"
    fig.savefig(out, dpi=300)
    plt.close()
    print(f"Saved: {out}")
    return out


def fig_scatter_properties(lib: BlockLibrary):
    """Scatter: HBD vs Rotatable bonds for all 400 blocks, colored by backbone."""
    fig, ax = plt.subplots(figsize=(5, 4))

    colors = {"alpha": "#4C72B0", "beta": "#55A868", "gamma": "#C44E52"}
    markers = {"alpha": "o", "beta": "s", "gamma": "^"}
    sizes = {"small": 15, "med": 40, "large": 80}

    for bt in ["alpha", "beta", "gamma"]:
        blocks = lib.blocks_by_property(mc_type=bt)
        x = [b.rot_total + np.random.uniform(-0.2, 0.2) for b in blocks]
        y = [b.sc_hbd + np.random.uniform(-0.15, 0.15) for b in blocks]
        s = [sizes.get(b.sc_bulk, 30) for b in blocks]
        ax.scatter(x, y, c=colors[bt], marker=markers[bt], s=s,
                   alpha=0.6, edgecolors="white", linewidth=0.3,
                   label=f"{bt} (n={len(blocks)})")

    ax.set_xlabel("Total Rotatable Bonds")
    ax.set_ylabel("Side Chain HBD")
    ax.set_title("Property Space: Flexibility vs H-Bond Donors", fontweight="bold")
    ax.legend(loc="upper right")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    size_legend = [
        plt.scatter([], [], s=15, c="gray", alpha=0.5, marker="o", label="Small"),
        plt.scatter([], [], s=40, c="gray", alpha=0.5, marker="o", label="Medium"),
        plt.scatter([], [], s=80, c="gray", alpha=0.5, marker="o", label="Large"),
    ]
    leg2 = ax.legend(handles=size_legend, title="Bulk", loc="center right",
                     fontsize=6, title_fontsize=7)
    ax.add_artist(ax.legend(loc="upper right"))

    out = FIG_DIR / "fig_scatter_hbd_rot.png"
    fig.savefig(out, dpi=300)
    plt.close()
    print(f"Saved: {out}")
    return out


def fig_radar(lib: BlockLibrary):
    """Radar chart: Property fingerprints of 4 representative blocks."""
    examples = [
        ("A01", "Alpha-L, Ile analog"),
        ("K08", "Beta-D, Lys(+) analog"),
        ("a07", "Alpha-D, Pro(NCY)"),
        ("S05", "Gamma-L, Ile analog"),
    ]

    categories = ["Charge", "HBD", "HBA", "Rot Bonds", "Bulk\n(encoded)"]

    def encode_bulk(b):
        return {"small": 1, "med": 2, "large": 3}.get(b, 0)

    fig, ax = plt.subplots(figsize=(5, 5), subplot_kw=dict(polar=True))

    n = len(categories)
    angles = np.linspace(0, 2 * np.pi, n, endpoint=False).tolist()
    angles += angles[:1]

    colors = ["#4C72B0", "#C44E52", "#55A868", "#DD8452"]

    maxvals = [1, 10, 14, 12, 3]

    for i, (bid, desc) in enumerate(examples):
        b = lib[bid]
        values = [
            abs(b.charge) / maxvals[0],
            b.sc_hbd / maxvals[1],
            b.sc_hba / maxvals[2],
            b.rot_total / maxvals[3],
            encode_bulk(b.sc_bulk) / maxvals[4],
        ]
        values += values[:1]

        ax.plot(angles, values, "o-", linewidth=1.5, color=colors[i],
                label=f"{bid} ({desc})", markersize=4)
        ax.fill(angles, values, alpha=0.1, color=colors[i])

    ax.set_xticks(angles[:-1])
    ax.set_xticklabels(categories, fontsize=7)
    ax.set_ylim(0, 1.1)
    ax.set_yticks([0.25, 0.5, 0.75, 1.0])
    ax.set_yticklabels(["25%", "50%", "75%", "max"], fontsize=5, color="gray")

    ax.set_title("Property Fingerprints of Representative Blocks",
                 fontsize=9, fontweight="bold", pad=20)
    ax.legend(loc="lower right", bbox_to_anchor=(1.3, -0.05), fontsize=6.5)

    out = FIG_DIR / "fig_radar_property_fingerprints.png"
    fig.savefig(out, dpi=300)
    plt.close()
    print(f"Saved: {out}")
    return out


def fig_donut(lib: BlockLibrary):
    """Double donut: inner ring = backbone type, outer ring = N-modification."""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7, 3.5))

    # Donut 1: Backbone composition
    bt_counts = Counter(lib[bid].mc_type for bid in lib.all_ids)
    bt_order = ["alpha", "beta", "gamma"]
    bt_vals = [bt_counts[bt] for bt in bt_order]
    bt_colors = ["#4C72B0", "#55A868", "#C44E52"]
    bt_labels = [f"Alpha\n{bt_vals[0]}", f"Beta\n{bt_vals[1]}", f"Gamma\n{bt_vals[2]}"]

    wedges1, texts1, autotexts1 = ax1.pie(
        bt_vals, labels=bt_labels, colors=bt_colors, autopct="%1.0f%%",
        pctdistance=0.75, startangle=90, textprops={"fontsize": 7},
        wedgeprops=dict(width=0.4, edgecolor="white", linewidth=2))
    for t in autotexts1:
        t.set_fontsize(7)
        t.set_fontweight("bold")
    ax1.set_title("Backbone Type", fontsize=9, fontweight="bold")

    # Donut 2: N-modification
    nmod_counts = Counter(lib[bid].mc_nmod for bid in lib.all_ids)
    nmod_order = ["NO", "NCY", "NME"]
    nmod_vals = [nmod_counts.get(nm, 0) for nm in nmod_order]
    nmod_colors = ["#999999", "#8172B2", "#DD8452"]
    nmod_labels = [f"None\n{nmod_vals[0]}", f"N-cyclic\n{nmod_vals[1]}", f"N-methyl\n{nmod_vals[2]}"]

    wedges2, texts2, autotexts2 = ax2.pie(
        nmod_vals, labels=nmod_labels, colors=nmod_colors, autopct="%1.0f%%",
        pctdistance=0.75, startangle=90, textprops={"fontsize": 7},
        wedgeprops=dict(width=0.4, edgecolor="white", linewidth=2))
    for t in autotexts2:
        t.set_fontsize(7)
        t.set_fontweight("bold")
    ax2.set_title("N-Modification", fontsize=9, fontweight="bold")

    fig.suptitle("Library Composition (N=400)", fontsize=10, fontweight="bold", y=1.02)
    plt.tight_layout()

    out = FIG_DIR / "fig_donut_composition.png"
    fig.savefig(out, dpi=300)
    plt.close()
    print(f"Saved: {out}")
    return out


def fig_acceptance_rate():
    """Horizontal bar: Random baseline acceptance rates with difficulty annotations."""
    profiles = [
        ("Unconstrained", 100.0, "No constraints beyond valid IDs + length 6"),
        ("Permeable", 2.07, "charge=0, HBD≤1, ≥2 NMe/NCY, ≥3 large"),
        ("Charged Binder", 1.82, "charge=+2, HBD≥2, ≥1 aromatic"),
        ("Rigid Scaffold", 0.81, "rot≤18, ≥3 beta backbone"),
    ]

    fig, ax = plt.subplots(figsize=(7, 3))

    names = [p[0] for p in profiles]
    rates = [p[1] for p in profiles]
    descs = [p[2] for p in profiles]
    colors = ["#55A868", "#4C72B0", "#8172B2", "#C44E52"]

    y_pos = range(len(profiles))
    bars = ax.barh(y_pos, rates, color=colors, height=0.6, edgecolor="white", linewidth=0.5)

    ax.set_yticks(y_pos)
    ax.set_yticklabels(names, fontsize=8)
    ax.set_xlabel("Random Acceptance Rate (%)")
    ax.set_title("Constraint Difficulty: What Fraction of Random Sequences Pass?",
                 fontsize=9, fontweight="bold")
    ax.set_xscale("log")
    ax.set_xlim(0.5, 200)

    for i, (bar, rate, desc) in enumerate(zip(bars, rates, descs)):
        if rate > 10:
            ax.text(rate - 2, i, f"{rate:.0f}%", ha="right", va="center",
                    fontsize=8, fontweight="bold", color="white")
        else:
            ax.text(rate + rate * 0.3, i, f"{rate:.2f}%", ha="left", va="center",
                    fontsize=8, fontweight="bold", color=colors[i])

        ax.text(150, i, desc, ha="right", va="center", fontsize=5.5,
                color="#888888", fontstyle="italic")

    ax.axvline(x=5, color="gray", linestyle="--", linewidth=0.5, alpha=0.5)
    ax.text(5.5, len(profiles) - 0.3, "5% threshold", fontsize=5, color="gray", alpha=0.7)

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.invert_yaxis()

    out = FIG_DIR / "fig_acceptance_rates.png"
    fig.savefig(out, dpi=300)
    plt.close()
    print(f"Saved: {out}")
    return out


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--fig", default="all",
                        choices=["all", "heatmap", "scatter", "radar", "donut", "acceptance"])
    args = parser.parse_args()

    FIG_DIR.mkdir(parents=True, exist_ok=True)
    lib = BlockLibrary()

    generated = []
    if args.fig in ("all", "heatmap"):
        generated.append(fig_heatmap(lib))
    if args.fig in ("all", "scatter"):
        generated.append(fig_scatter_properties(lib))
    if args.fig in ("all", "radar"):
        generated.append(fig_radar(lib))
    if args.fig in ("all", "donut"):
        generated.append(fig_donut(lib))
    if args.fig in ("all", "acceptance"):
        generated.append(fig_acceptance_rate())

    print(f"\nGenerated {len(generated)} figures in {FIG_DIR}")


if __name__ == "__main__":
    main()
