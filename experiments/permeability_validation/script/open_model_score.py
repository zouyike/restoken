"""Cross-model ResToken agnosticism analysis.

For the open-source-LLM sweep (experiments/output/open_model_bench): assemble SMILES
from each model's restoken generations, score permeability with the held-out-validated
PeptideCLM judge (perm_model.pt, test rho=0.74), and produce:
  1. validity/constraint-pass table across models  -> open_model_validity.csv
  2. per-model permeable-vs-unconstrained permeability shift figure + stats
     -> open_model_permeability.png / open_model_perm_stats.json

Run on a GPU node (esm3 env). Imports the shared assembler + scoring model — no copies.
"""
import os, sys, json, glob, argparse
from collections import defaultdict
import numpy as np
import pandas as pd
import torch
from scipy.stats import mannwhitneyu

WD = "/scratch/genesis/NCAA_tokenization/experiments/permeability_validation"
sys.path.insert(0, f"{WD}/script")
sys.path.insert(0, "/scratch/genesis/NCAA_tokenization/restoken/src")
from cyclic_assembler import CyclicPeptideAssembler
from train_perm import Regressor, get_tok, DS, collate, MODEL_DIR
from torch.utils.data import DataLoader

PALETTE = {"navy": "#002060", "brown": "#B36B3C", "teal": "#2E7F7A",
           "grayblue": "#4A5D73", "darkred": "#990033"}
CONSTRAINTS = ["unconstrained", "permeable", "charged_binder", "rigid_scaffold"]
# display order / nice names (sorted by descending unconstrained validity for the figure)
MODELS = ["gemma-3-4b-it", "gemma-3-12b-it", "Qwen3-32B",
          "Llama-3.3-70B", "Mistral-Small-24B", "Mistral-Nemo-12B"]


def restoken_seqs(path):
    out = []
    for block in json.load(open(path)):
        for line in str(block).splitlines():
            line = line.strip()
            if line:
                out.append(line)
    return out


def assemble(in_dir):
    """-> DataFrame[model, constraint, token_seq, smiles] (dedup within model+constraint)."""
    asm = CyclicPeptideAssembler()
    rows, seen = [], set()
    for path in sorted(glob.glob(f"{in_dir}/*_restoken_*_len6_raw.json")):
        base = os.path.basename(path)
        model = base.split("_restoken_")[0]
        constraint = base.split("_restoken_")[1].split("_len6")[0]
        for seq in restoken_seqs(path):
            try:
                smi = asm.assemble(seq)
            except Exception:
                smi = None
            if not smi:
                continue
            k = (model, constraint, smi)
            if k in seen:
                continue
            seen.add(k)
            rows.append((model, constraint, seq, smi))
    return pd.DataFrame(rows, columns=["model", "constraint", "token_seq", "smiles"])


def validity_table(in_dir):
    rows = []
    for path in sorted(glob.glob(f"{in_dir}/*_restoken_*_len6_results.json")):
        d = json.load(open(path))
        rows.append({
            "model": d["model"], "constraint": d["profile"],
            "n_parsed": d.get("n_parsed"), "n_valid": d.get("n_valid"),
            "validity_rate": d.get("validity_rate"),
            "constraint_pass_rate": d.get("constraint_pass_rate"),
            "uniqueness": d.get("uniqueness"),
        })
    return pd.DataFrame(rows)


def score(df):
    dev = "cuda" if torch.cuda.is_available() else "cpu"
    ck = torch.load(f"{WD}/output/perm_model.pt", map_location=dev)
    mu, sd, maxlen = ck["mu"], ck["sd"], ck["maxlen"]
    tok = get_tok()
    model = Regressor(MODEL_DIR).to(dev)
    model.load_state_dict(ck["model"])
    model.eval()
    ds = DS(df.smiles.tolist(), None, tok, maxlen)
    dl = DataLoader(ds, batch_size=64, shuffle=False,
                    collate_fn=lambda b: collate(b, tok.pad_token_id))
    out = []
    with torch.no_grad():
        for b in dl:
            p = model(b["input_ids"].to(dev), b["attention_mask"].to(dev)).cpu().numpy()
            out.append(np.atleast_1d(p))
    return np.concatenate(out) * sd + mu


def cliffs_delta(a, b):
    a, b = np.asarray(a), np.asarray(b)
    gt = sum((a[:, None] > b[None, :]).sum(axis=1))
    lt = sum((a[:, None] < b[None, :]).sum(axis=1))
    return (gt - lt) / (len(a) * len(b))


def boot_ci(a, b, n=5000, seed=0):
    rng = np.random.RandomState(seed)
    a, b = np.asarray(a), np.asarray(b)
    d = [np.median(rng.choice(a, len(a))) - np.median(rng.choice(b, len(b))) for _ in range(n)]
    return float(np.percentile(d, 2.5)), float(np.percentile(d, 97.5))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in_dir",
                    default="/scratch/genesis/NCAA_tokenization/experiments/output/open_model_bench")
    ap.add_argument("--out_prefix", default=f"{WD}/output/open_model")
    args = ap.parse_args()

    # ---- 1. validity / constraint-pass table ----
    vt = validity_table(args.in_dir)
    vt.to_csv(f"{args.out_prefix}_validity.csv", index=False)
    print("=== validity table ===")
    print(vt.pivot(index="model", columns="constraint", values="validity_rate").to_string())

    # ---- 2. assemble + score ----
    df = assemble(args.in_dir)
    print(f"\nassembled {len(df)} unique SMILES across {df.model.nunique()} models")
    print(df.groupby(["model", "constraint"]).size().to_string())
    df["perm"] = score(df)
    df.to_csv(f"{args.out_prefix}_scores.csv", index=False)

    # ---- 3. per-model permeable vs unconstrained stats ----
    stats = []
    for m in df.model.unique():
        a = df[(df.model == m) & (df.constraint == "permeable")]["perm"].values
        b = df[(df.model == m) & (df.constraint == "unconstrained")]["perm"].values
        rec = {"model": m, "n_permeable": len(a), "n_unconstrained": len(b),
               "median_permeable": float(np.median(a)) if len(a) else None,
               "median_unconstrained": float(np.median(b)) if len(b) else None}
        if len(a) >= 10 and len(b) >= 10:
            u, p = mannwhitneyu(a, b, alternative="greater")
            lo, hi = boot_ci(a, b)
            rec.update({"delta_median": float(np.median(a) - np.median(b)),
                        "cliffs_delta": float(cliffs_delta(a, b)),
                        "p_one_sided": float(p), "median_diff_95ci": [lo, hi]})
        stats.append(rec)
    json.dump(stats, open(f"{args.out_prefix}_perm_stats.json", "w"), indent=2)
    print("\n=== permeable vs unconstrained (per model) ===")
    print(json.dumps(stats, indent=2))

    # ---- 4. figure: per-model permeable vs unconstrained medians ----
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    order = [m for m in MODELS if m in df.model.unique()]
    fig, ax = plt.subplots(figsize=(10, 5))
    x = np.arange(len(order)); w = 0.38
    med_unc = [np.median(df[(df.model == m) & (df.constraint == "unconstrained")]["perm"]) for m in order]
    med_perm = [np.median(df[(df.model == m) & (df.constraint == "permeable")]["perm"]) for m in order]
    ax.bar(x - w/2, med_unc, w, color=PALETTE["grayblue"], label="unconstrained")
    ax.bar(x + w/2, med_perm, w, color=PALETTE["navy"], label="permeable")
    for i, m in enumerate(order):
        n_p = (df[(df.model == m) & (df.constraint == "permeable")]).shape[0]
        ax.text(i + w/2, med_perm[i], f"n={n_p}", ha="center", va="bottom", fontsize=7)
    ax.set_xticks(x); ax.set_xticklabels(order, rotation=20, ha="right", fontsize=9)
    ax.set_ylabel("Median predicted log Papp (higher = more permeable)")
    ax.set_title("ResToken permeability constraint replicates across open LLMs\n"
                 "(permeable instruction shifts generations more permeable)")
    ax.legend(); ax.grid(axis="y", alpha=0.3)
    plt.tight_layout()
    plt.savefig(f"{args.out_prefix}_permeability.png", dpi=200)
    print(f"\nwrote {args.out_prefix}_permeability.png")


if __name__ == "__main__":
    main()
