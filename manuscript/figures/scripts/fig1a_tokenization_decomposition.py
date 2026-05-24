#!/usr/bin/env python
"""
Figure 1a: ResToken Tokenization Decomposition
Shows how NCAA chemical structures are decomposed into semantic token properties.

3 representative blocks:
  - A01: alpha-L, large bulk (tert-leucine/isoleucine analog)
  - N15: beta-D, small bulk (beta-alanine analog)
  - s01: gamma-L, N-cyclic (gamma-proline analog)

Output: fig1a_tokenization_decomposition.png/.pdf
"""

import json
import io
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.gridspec import GridSpec
import numpy as np
from PIL import Image

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import rdMolDraw2D


# ── Configuration ───────────────────────────────────────────────────────
FONT = 'DejaVu Sans'

DATA_DIR = Path("/scratch/genesis/NCAA_tokenization/restoken/data")
OUT_DIR = Path("/scratch/genesis/NCAA_tokenization/manuscript/figures/generated")
OUT_DIR.mkdir(parents=True, exist_ok=True)

BLOCK_IDS = ["A01", "N15", "s01"]

BB_COLORS = {
    "alpha": "#2C6FBB",
    "beta":  "#1E8449",
    "gamma": "#D4740E",
}

BB_COLORS_LIGHT = {
    "alpha": "#DAEAF8",
    "beta":  "#D5F5E3",
    "gamma": "#FDEBD0",
}

ROW_LABELS = {
    "A01": r"$\alpha$-amino acid",
    "N15": r"$\beta$-amino acid",
    "s01": r"$\gamma$-amino acid",
}

ROW_DESCRIPTIONS = {
    "A01": "isoleucine analog",
    "N15": r"$\beta$-alanine analog",
    "s01": r"$\gamma$-proline analog",
}

PROP_KEYS = [
    ("class",        "class"),
    ("chirality",    "chirality"),
    ("charge_label", "charge"),
    ("mc_type",      "backbone"),
    ("mc_nmod",      "N-mod"),
    ("sc_bulk",      "bulk"),
    ("polarity_bin", "polarity"),
    ("flex_bin",     "flexibility"),
]


# ── Data loading ────────────────────────────────────────────────────────
with open(DATA_DIR / "bb_dict_backend_v11.json") as f:
    backend = json.load(f)
with open(DATA_DIR / "bb_dict_llm_v11.json") as f:
    llm_data = json.load(f)

blocks_b = {b['id']: b for b in backend['blocks']}
blocks_l = {b['id']: b for b in llm_data['blocks']}


# ── Helpers ─────────────────────────────────────────────────────────────
def render_mol_2d(smiles, size=(600, 450)):
    """Render a 2D molecule image from SMILES, return PIL Image."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Cannot parse SMILES: {smiles}")
    AllChem.Compute2DCoords(mol)

    drawer = rdMolDraw2D.MolDraw2DCairo(size[0], size[1])
    opts = drawer.drawOptions()
    opts.bondLineWidth = 2.8
    opts.minFontSize = 18
    opts.maxFontSize = 24
    opts.padding = 0.15
    opts.backgroundColour = (1, 1, 1, 0)
    drawer.DrawMolecule(mol)
    drawer.FinishDrawing()
    png_data = drawer.GetDrawingText()
    return Image.open(io.BytesIO(png_data))


def draw_property_card(ax, props, bb_type, block_id):
    """Draw a property card on the given axes."""
    color = BB_COLORS[bb_type]
    color_light = BB_COLORS_LIGHT[bb_type]

    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_aspect('auto')
    ax.axis('off')

    n_props = len(PROP_KEYS)
    card_top = 0.97
    card_bottom = 0.03
    card_left = 0.02
    card_right = 0.98
    card_height = card_top - card_bottom
    title_h = 0.08
    row_area_top = card_top - title_h
    row_height = (row_area_top - card_bottom) / n_props

    # Card background
    card_bg = mpatches.FancyBboxPatch(
        (card_left, card_bottom), card_right - card_left, card_height,
        boxstyle="round,pad=0.012",
        facecolor='white', edgecolor=color, linewidth=2.0
    )
    ax.add_patch(card_bg)

    # Title bar
    title_rect = mpatches.FancyBboxPatch(
        (card_left, row_area_top), card_right - card_left, title_h,
        boxstyle="round,pad=0.01",
        facecolor=color, edgecolor=color, linewidth=1.5
    )
    ax.add_patch(title_rect)
    # Mask bottom rounding of title
    mask_rect = mpatches.Rectangle(
        (card_left + 0.008, row_area_top), card_right - card_left - 0.016, title_h * 0.45,
        facecolor=color, edgecolor='none'
    )
    ax.add_patch(mask_rect)

    ax.text(0.50, card_top - title_h * 0.45, f"[{block_id}]",
            ha='center', va='center', fontsize=10.5, fontweight='bold',
            color='white', fontfamily=FONT)

    # Property rows
    for i, (key, label) in enumerate(PROP_KEYS):
        y_center = row_area_top - (i + 0.5) * row_height
        y_bottom = row_area_top - (i + 1) * row_height
        val = props.get(key, "?")

        # Alternating row background
        if i % 2 == 0:
            row_bg = mpatches.Rectangle(
                (card_left + 0.012, y_bottom + 0.002),
                card_right - card_left - 0.024,
                row_height - 0.004,
                facecolor=color_light, edgecolor='none', alpha=0.55
            )
            ax.add_patch(row_bg)

        # Separator line
        if i > 0:
            ax.plot(
                [card_left + 0.04, card_right - 0.04],
                [row_area_top - i * row_height, row_area_top - i * row_height],
                color='#DDDDDD', linewidth=0.4, zorder=2
            )

        # Label
        ax.text(card_left + 0.06, y_center, label, ha='left', va='center',
                fontsize=8.5, fontfamily=FONT, color='#555555')

        # Value -- color-coded
        val_color = '#222222'
        fw = 'bold'
        if key == 'mc_type':
            val_color = color
        elif key == 'mc_nmod' and val != 'NO':
            val_color = '#C0392B'
            highlight = mpatches.FancyBboxPatch(
                (card_right - 0.20, y_center - 0.028), 0.16, 0.056,
                boxstyle="round,pad=0.006",
                facecolor='#FADBD8', edgecolor='#E74C3C', linewidth=0.7
            )
            ax.add_patch(highlight)
        elif key == 'chirality':
            val_color = '#8E44AD' if val == 'D' else '#2471A3'

        ax.text(card_right - 0.06, y_center, str(val), ha='right', va='center',
                fontsize=8.5, fontfamily=FONT, fontweight=fw,
                color=val_color, zorder=5)


# ── Main Figure ─────────────────────────────────────────────────────────
def create_figure():
    fig_width = 7.0
    fig_height = 6.2
    fig = plt.figure(figsize=(fig_width, fig_height), dpi=300, facecolor='white')

    gs = GridSpec(
        3, 3,
        figure=fig,
        width_ratios=[0.40, 0.07, 0.53],
        height_ratios=[1, 1, 1],
        hspace=0.22,
        wspace=0.02,
        left=0.03, right=0.97,
        top=0.88, bottom=0.03
    )

    for row_idx, block_id in enumerate(BLOCK_IDS):
        b_data = blocks_b[block_id]
        l_data = blocks_l[block_id]
        bb_type = l_data['mc_type']
        color = BB_COLORS[bb_type]
        smiles = b_data['structure']['aa_smiles']

        # ── Left panel: 2D structure ──
        ax_mol = fig.add_subplot(gs[row_idx, 0])
        mol_img = render_mol_2d(smiles, size=(600, 450))
        ax_mol.imshow(mol_img, aspect='auto', interpolation='lanczos')
        ax_mol.axis('off')

        # Colored left accent bar
        mol_pos = ax_mol.get_position()
        fig.patches.append(mpatches.FancyBboxPatch(
            (mol_pos.x0 - 0.006, mol_pos.y0 + 0.004),
            0.009, mol_pos.height - 0.008,
            boxstyle="round,pad=0.002",
            facecolor=color, edgecolor='none',
            transform=fig.transFigure, clip_on=False, zorder=5
        ))

        # Row label below molecule
        ax_mol.text(
            0.5, -0.06,
            f'{ROW_LABELS[block_id]}  ({ROW_DESCRIPTIONS[block_id]})',
            ha='center', va='top',
            fontsize=8, fontfamily=FONT,
            color=color, fontstyle='italic',
            transform=ax_mol.transAxes
        )

        # ── Middle: arrow ──
        ax_arrow = fig.add_subplot(gs[row_idx, 1])
        ax_arrow.set_xlim(0, 1)
        ax_arrow.set_ylim(0, 1)
        ax_arrow.axis('off')

        ax_arrow.annotate(
            '', xy=(0.92, 0.5), xytext=(0.08, 0.5),
            arrowprops=dict(
                arrowstyle='-|>,head_width=0.4,head_length=0.3',
                color=color, lw=2.5,
                connectionstyle='arc3,rad=0'
            ),
        )

        # ── Right panel: property card ──
        ax_card = fig.add_subplot(gs[row_idx, 2])
        draw_property_card(ax_card, l_data, bb_type, block_id)

    # ── Figure title ──
    fig.text(0.50, 0.965,
             'Semantic Decomposition of NCAA Building Blocks',
             ha='center', va='top', fontsize=12, fontweight='bold',
             fontfamily=FONT, color='#1A1A1A')

    # ── Column headers ──
    fig.text(0.22, 0.92,
             'Chemical Structure',
             ha='center', va='center', fontsize=9.5, fontweight='bold',
             fontfamily=FONT, color='#444444')
    fig.text(0.735, 0.92,
             'Semantic Token Properties',
             ha='center', va='center', fontsize=9.5, fontweight='bold',
             fontfamily=FONT, color='#444444')

    # Separator line
    fig.patches.append(mpatches.FancyBboxPatch(
        (0.04, 0.898), 0.92, 0.0012,
        boxstyle="round,pad=0.0005",
        facecolor='#CCCCCC', edgecolor='none',
        transform=fig.transFigure, clip_on=False
    ))

    # Save
    for ext in ['png', 'pdf']:
        out_path = OUT_DIR / f"fig1a_tokenization_decomposition.{ext}"
        fig.savefig(out_path, dpi=300, bbox_inches='tight',
                    facecolor='white', edgecolor='none')
        print(f"Saved: {out_path}")

    plt.close(fig)


if __name__ == "__main__":
    create_figure()
