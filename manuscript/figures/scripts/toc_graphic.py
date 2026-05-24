"""TOC/Abstract graphic for JCIM submission.
Shows: NCAA structure → ResToken (semantic properties) → LLM → valid peptide sequences
ACS TOC graphic: 3.25 x 1.75 inches, 300 DPI
"""
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
import numpy as np
import sys, os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', '..'))

fig, ax = plt.subplots(1, 1, figsize=(3.25, 1.75), dpi=300)
ax.set_xlim(0, 10)
ax.set_ylim(0, 5)
ax.axis('off')

font = 'DejaVu Sans'

# --- Left panel: NCAA building block as token ---
token_x, token_y = 1.2, 2.7
token_w, token_h = 2.0, 3.0

# Token box
token_box = FancyBboxPatch((token_x - token_w/2, token_y - token_h/2),
                            token_w, token_h,
                            boxstyle="round,pad=0.1",
                            facecolor='#E3F2FD', edgecolor='#1565C0',
                            linewidth=1.2)
ax.add_patch(token_box)

# Token ID
ax.text(token_x, token_y + 1.15, 'A01', fontsize=9, fontweight='bold',
        ha='center', va='center', fontfamily=font, color='#1565C0')

# Property labels (compact)
props = [
    ('alpha', '#4CAF50'),
    ('L-chiral', '#FF9800'),
    ('neutral', '#9E9E9E'),
    ('large', '#E91E63'),
    ('low pol.', '#00BCD4'),
]
for i, (label, color) in enumerate(props):
    y_pos = token_y + 0.55 - i * 0.45
    pill = FancyBboxPatch((token_x - 0.75, y_pos - 0.15), 1.5, 0.32,
                           boxstyle="round,pad=0.05",
                           facecolor=color, edgecolor='none', alpha=0.2)
    ax.add_patch(pill)
    ax.text(token_x, y_pos, label, fontsize=5, ha='center', va='center',
            fontfamily=font, color=color, fontweight='bold')

ax.text(token_x, token_y - token_h/2 - 0.05, 'ResToken', fontsize=6,
        ha='center', va='top', fontfamily=font, fontweight='bold', color='#1565C0')

# --- Arrow 1: Token → LLM ---
ax.annotate('', xy=(3.6, 2.5), xytext=(2.5, 2.5),
            arrowprops=dict(arrowstyle='->', color='#555', lw=1.2))

# --- Middle: LLM box ---
llm_x, llm_y = 4.6, 2.5
llm_box = FancyBboxPatch((llm_x - 0.9, llm_y - 0.6), 1.8, 1.2,
                           boxstyle="round,pad=0.12",
                           facecolor='#FFF3E0', edgecolor='#E65100',
                           linewidth=1.2)
ax.add_patch(llm_box)
ax.text(llm_x, llm_y + 0.15, 'LLM', fontsize=9, fontweight='bold',
        ha='center', va='center', fontfamily=font, color='#E65100')
ax.text(llm_x, llm_y - 0.25, 'zero-shot', fontsize=5,
        ha='center', va='center', fontfamily=font, color='#E65100', style='italic')

# --- Arrow 2: LLM → Output ---
ax.annotate('', xy=(6.5, 2.5), xytext=(5.7, 2.5),
            arrowprops=dict(arrowstyle='->', color='#555', lw=1.2))

# --- Right: Output sequences with validity ---
out_x = 8.0
seq_data = [
    ('A01-K03-N12-S05-E02-a07', True),
    ('A15-N33-K22-a01-E08-S11', True),
    ('T09-F12-K03-D01-A15-s03', True),
]
ax.text(out_x, 4.2, 'Generated Peptides', fontsize=5.5, ha='center', va='center',
        fontfamily=font, fontweight='bold', color='#333')

for i, (seq, valid) in enumerate(seq_data):
    y = 3.5 - i * 0.7
    color = '#2E7D32' if valid else '#C62828'
    bg = '#E8F5E9' if valid else '#FFEBEE'
    seq_box = FancyBboxPatch((out_x - 1.7, y - 0.22), 3.4, 0.44,
                              boxstyle="round,pad=0.05",
                              facecolor=bg, edgecolor=color,
                              linewidth=0.6, alpha=0.8)
    ax.add_patch(seq_box)
    ax.text(out_x, y, seq, fontsize=4.2, ha='center', va='center',
            fontfamily='DejaVu Sans Mono', color=color)

# --- Bottom: key result banner ---
banner = FancyBboxPatch((1.0, 0.15), 8.0, 0.55,
                         boxstyle="round,pad=0.08",
                         facecolor='#1565C0', edgecolor='none', alpha=0.9)
ax.add_patch(banner)
ax.text(5.0, 0.42, '90% validity  |  100% positional control  |  123× enrichment',
        fontsize=5, ha='center', va='center', fontfamily=font,
        color='white', fontweight='bold')

# --- "×400" annotation on token ---
ax.text(token_x, token_y + token_h/2 + 0.2, '400 blocks', fontsize=5,
        ha='center', va='bottom', fontfamily=font, color='#1565C0', style='italic')

plt.tight_layout(pad=0.1)
outdir = os.path.join(os.path.dirname(__file__), '..', 'generated')
for ext in ['png', 'pdf']:
    fig.savefig(os.path.join(outdir, f'toc_graphic.{ext}'),
                dpi=300, bbox_inches='tight', facecolor='white')
print("TOC graphic saved.")
plt.close()
