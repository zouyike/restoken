import io, json
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager as fm
from scipy.stats import spearmanr
from PIL import Image, ImageDraw, ImageFont

WD = "/scratch/genesis/NCAA_tokenization/experiments/permeability_validation"
OUT = f"{WD}/output"
NAVY="#002060"; BROWN="#B36B3C"; TEAL="#2E7F7A"; GRAYBLUE="#4A5D73"; DARKRED="#990033"
plt.rcParams.update({"font.size": 13, "axes.edgecolor": NAVY, "axes.labelcolor": NAVY,
                     "xtick.color": NAVY, "ytick.color": NAVY, "axes.linewidth": 1.2})


def fig_to_image(fig, target_w):
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=150, facecolor="white", bbox_inches="tight", pad_inches=0.12)
    plt.close(fig); buf.seek(0)
    img = Image.open(buf).convert("RGB")
    if img.width != target_w:
        img = img.resize((target_w, int(img.height*target_w/img.width)), Image.LANCZOS)
    return img


def panel_calibration(w):
    d = pd.read_csv(f"{OUT}/test_predictions.csv")
    m = json.load(open(f"{OUT}/test_metrics.json"))
    fig, ax = plt.subplots(figsize=(5.2, 5.0))
    ax.scatter(d.y_true, d.y_pred, s=14, c=NAVY, alpha=0.35, edgecolors="none")
    lo, hi = -10.5, -3.5
    ax.plot([lo, hi], [lo, hi], color=DARKRED, lw=1.6, ls="--", zorder=1)
    ax.set_xlim(lo, hi); ax.set_ylim(lo, hi)
    ax.set_xlabel("Experimental log P$_{app}$"); ax.set_ylabel("Predicted log P$_{app}$")
    ax.set_title("Held-out test calibration", color=NAVY, fontweight="bold", pad=10)
    txt = (f"Spearman $\\rho$ = {m['spearman']:.2f}\n"
           f"Pearson $r$ = {m['pearson']:.2f}\n"
           f"RMSE = {m['rmse']:.2f}\n$n$ = {len(d)}")
    ax.text(0.04, 0.96, txt, transform=ax.transAxes, va="top", ha="left", fontsize=12,
            color=NAVY, bbox=dict(boxstyle="round,pad=0.4", fc="white", ec=GRAYBLUE, lw=1.0))
    ax.grid(alpha=0.25)
    return fig_to_image(fig, w)


def panel_ecdf(w):
    df = pd.read_csv(f"{OUT}/generation_scores.csv")
    def g(src, con=None):
        m = df.source == src
        if con: m &= df.constraint == con
        return df.loc[m, "perm"].values
    groups = [("ResToken, permeable-constrained", g("restoken", "permeable"), NAVY),
              ("ResToken, unconstrained", g("restoken", "unconstrained"), TEAL),
              ("Random NCAA library", g("random"), BROWN)]
    fig, ax = plt.subplots(figsize=(5.4, 5.0))
    for label, v, c in groups:
        xs = np.sort(v); ys = np.arange(1, len(xs)+1)/len(xs)
        ax.plot(xs, ys, color=c, lw=2.4, label=f"{label} ($n$={len(v)})")
    ax.set_xlabel("Predicted log P$_{app}$  (higher = more permeable)")
    ax.set_ylabel("Cumulative fraction")
    ax.set_title("Permeability of generated peptides", color=NAVY, fontweight="bold", pad=10)
    ax.legend(fontsize=10.5, loc="upper left", frameon=False)
    ax.set_xlim(-11, -4); ax.grid(alpha=0.25)
    return fig_to_image(fig, w)


def panel_effects(w):
    res = json.load(open(f"{OUT}/comparison_stats.json"))
    labels = ["vs. ResToken\nunconstrained", "vs. random\nlibrary (all)", "vs. random,\npermeable-filtered"]
    order = ["restoken_permeable > restoken_unconstrained",
             "restoken_permeable > random_all",
             "restoken_permeable > random_permeable"]
    by = {r["contrast"]: r for r in res}
    fig, ax = plt.subplots(figsize=(5.4, 5.0))
    ys = np.arange(len(order))[::-1]
    for y, key, lab in zip(ys, order, labels):
        r = by[key]
        lo, hi = r["median_diff_95ci"]; mid = (lo+hi)/2
        ax.plot([lo, hi], [y, y], color=NAVY, lw=2.6, solid_capstyle="round")
        ax.plot(mid, y, "o", color=NAVY, ms=9)
        ax.text(hi+0.03, y, f"$\\delta$={r['cliffs_delta']:.2f}", va="center", ha="left",
                fontsize=11.5, color=DARKRED)
    ax.axvline(0, color=GRAYBLUE, lw=1.2, ls=":")
    ax.set_yticks(ys); ax.set_yticklabels(labels, fontsize=11.5)
    ax.set_xlabel("$\\Delta$ median log P$_{app}$  (ResToken permeable $-$ baseline)")
    ax.set_title("Effect size (95% CI)", color=NAVY, fontweight="bold", pad=10)
    ax.set_xlim(-0.15, 1.45); ax.grid(axis="x", alpha=0.25)
    return fig_to_image(fig, w)


def add_label(img, letter):
    im = img.copy(); d = ImageDraw.Draw(im)
    try:
        font = ImageFont.truetype("/usr/share/fonts/dejavu/DejaVuSans-Bold.ttf", 54)
    except Exception:
        font = ImageFont.load_default()
    d.text((14, 8), letter, fill=NAVY, font=font)
    return im


def hstack(imgs, gap=30):
    h = max(i.height for i in imgs); W = sum(i.width for i in imgs) + gap*(len(imgs)-1)
    canvas = Image.new("RGB", (W, h), "white"); x = 0
    for i in imgs:
        canvas.paste(i, (x, 0)); x += i.width + gap
    return canvas


def main():
    sub = 760
    a = add_label(panel_calibration(sub), "a")
    b = add_label(panel_ecdf(sub), "b")
    c = add_label(panel_effects(sub), "c")
    row = hstack([a, b, c], gap=34)
    M = 36
    canvas = Image.new("RGB", (row.width + 2*M, row.height + 2*M), "white")
    canvas.paste(row, (M, M))
    out = f"{OUT}/fig_permeability_validation.png"
    canvas.save(out, dpi=(300, 300))
    print("wrote", out, canvas.size)


if __name__ == "__main__":
    main()
