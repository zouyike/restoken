import os, sys, json, numpy as np, pandas as pd
import torch
from scipy.stats import mannwhitneyu

WD = "/scratch/genesis/NCAA_tokenization/experiments/permeability_validation"
sys.path.insert(0, f"{WD}/script")
from train_perm import Regressor, get_tok, DS, collate, MODEL_DIR
from torch.utils.data import DataLoader

PALETTE = {"navy": "#002060", "brown": "#B36B3C", "teal": "#2E7F7A",
           "grayblue": "#4A5D73", "darkred": "#990033"}


def score_smiles(smiles, model, tok, dev, mu, sd, maxlen, bs=64):
    ds = DS(smiles, None, tok, maxlen)
    dl = DataLoader(ds, batch_size=bs, shuffle=False,
                    collate_fn=lambda b: collate(b, tok.pad_token_id))
    model.eval(); out = []
    with torch.no_grad():
        for b in dl:
            p = model(b["input_ids"].to(dev), b["attention_mask"].to(dev)).cpu().numpy()
            out.append(np.atleast_1d(p))
    return np.concatenate(out) * sd + mu


def cliffs_delta(a, b):
    a = np.asarray(a); b = np.asarray(b)
    gt = sum((a[:, None] > b[None, :]).sum(axis=1))
    lt = sum((a[:, None] < b[None, :]).sum(axis=1))
    return (gt - lt) / (len(a) * len(b))


def boot_median_diff(a, b, n=5000, seed=0):
    rng = np.random.RandomState(seed)
    a = np.asarray(a); b = np.asarray(b)
    d = [np.median(rng.choice(a, len(a))) - np.median(rng.choice(b, len(b))) for _ in range(n)]
    return float(np.percentile(d, 2.5)), float(np.percentile(d, 97.5))


def compare(name, a, b):
    u, p = mannwhitneyu(a, b, alternative="greater")
    lo, hi = boot_median_diff(a, b)
    return {"contrast": name, "n_a": len(a), "n_b": len(b),
            "median_a": float(np.median(a)), "median_b": float(np.median(b)),
            "mannwhitney_U": float(u), "p_one_sided": float(p),
            "cliffs_delta": float(cliffs_delta(a, b)),
            "median_diff_95ci": [lo, hi]}


def main():
    dev = "cuda" if torch.cuda.is_available() else "cpu"
    ck = torch.load(f"{WD}/output/perm_model.pt", map_location=dev)
    mu, sd, maxlen = ck["mu"], ck["sd"], ck["maxlen"]
    tok = get_tok()
    model = Regressor(MODEL_DIR).to(dev); model.load_state_dict(ck["model"])

    df = pd.read_csv(f"{WD}/data/assembled_smiles.csv")
    df["perm"] = score_smiles(df.smiles.tolist(), model, tok, dev, mu, sd, maxlen)
    df.to_csv(f"{WD}/output/generation_scores.csv", index=False)

    def grp(source, constraint=None):
        m = df.source == source
        if constraint:
            m &= df.constraint == constraint
        return df.loc[m, "perm"].values

    rt_perm = grp("restoken", "permeable")
    rt_unc = grp("restoken", "unconstrained")
    rnd_perm = grp("random", "permeable")
    rnd_all = grp("random")

    results = [
        compare("restoken_permeable > restoken_unconstrained", rt_perm, rt_unc),
        compare("restoken_permeable > random_permeable", rt_perm, rnd_perm),
        compare("restoken_permeable > random_all", rt_perm, rnd_all),
    ]
    json.dump(results, open(f"{WD}/output/comparison_stats.json", "w"), indent=2)
    print(json.dumps(results, indent=2))

    # figure: ECDF of permeability by group
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    groups = [("ResToken permeable", rt_perm, PALETTE["navy"]),
              ("ResToken unconstrained", rt_unc, PALETTE["teal"]),
              ("Random permeable", rnd_perm, PALETTE["brown"])]
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.5))
    for label, vals, col in groups:
        xs = np.sort(vals); ys = np.arange(1, len(xs) + 1) / len(xs)
        ax1.plot(xs, ys, color=col, lw=2, label=f"{label} (n={len(vals)})")
    ax1.set_xlabel("Predicted log Papp (higher = more permeable)")
    ax1.set_ylabel("Cumulative fraction"); ax1.set_title("Permeability ECDF")
    ax1.legend(fontsize=8); ax1.grid(alpha=0.3)
    parts = ax2.violinplot([g[1] for g in groups], showmedians=True)
    for pc, (_, _, col) in zip(parts["bodies"], groups):
        pc.set_facecolor(col); pc.set_alpha(0.6)
    ax2.set_xticks([1, 2, 3]); ax2.set_xticklabels([g[0] for g in groups], rotation=15, fontsize=8)
    ax2.set_ylabel("Predicted log Papp"); ax2.set_title("Permeability distribution")
    ax2.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(f"{WD}/output/permeability_comparison.png", dpi=200)
    print("wrote output/permeability_comparison.png")


if __name__ == "__main__":
    main()
