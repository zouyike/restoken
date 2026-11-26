#!/usr/bin/env python3
"""Train PeptideCLM (or PeptideCLM-2) permeability regressor on a specified split
(random / cluster / source) and re-score all Gemini (8,129) and open-weight model macrocycles.
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import torch
import torch.nn as nn
from scipy.stats import mannwhitneyu, pearsonr, spearmanr
from torch.utils.data import DataLoader, Dataset
from transformers import AutoModel, AutoTokenizer, RoFormerModel

ROOT = Path("/scratch/genesis/NCAA_tokenization")
PERM_DIR = ROOT / "experiments" / "permeability_validation"
SPLIT_ROOT = ROOT / "experiments" / "dd_revision" / "perm_splits"


def cliffs_delta(a, b):
    a = np.asarray(a)
    b = np.asarray(b)
    gt = 0
    lt = 0
    for i in range(0, len(a), 500):
        chunk = a[i : i + 500, None]
        gt += int((chunk > b[None, :]).sum())
        lt += int((chunk < b[None, :]).sum())
    return (gt - lt) / (len(a) * len(b))


def boot_ci(a, b, n=3000, seed=0):
    rng = np.random.RandomState(seed)
    a = np.asarray(a)
    b = np.asarray(b)
    d = [np.median(rng.choice(a, len(a))) - np.median(rng.choice(b, len(b))) for _ in range(n)]
    return [float(np.percentile(d, 2.5)), float(np.percentile(d, 97.5))]


def load_tokenizer(model_type, model_dir):
    if model_type == "peptideclm":
        sys.path.insert(0, str(model_dir))
        from tokenizer.my_tokenizers import SMILES_SPE_Tokenizer
        return SMILES_SPE_Tokenizer(
            str(Path(model_dir) / "tokenizer" / "new_vocab.txt"),
            str(Path(model_dir) / "tokenizer" / "new_splits.txt"),
        )
    else:
        return AutoTokenizer.from_pretrained(model_dir, trust_remote_code=True)


class DS(Dataset):
    def __init__(self, smiles, y, tok, maxlen):
        self.s = list(smiles)
        self.y = None if y is None else np.asarray(y, dtype=np.float32)
        self.tok = tok
        self.maxlen = maxlen

    def __len__(self):
        return len(self.s)

    def __getitem__(self, i):
        e = self.tok(self.s[i], truncation=True, max_length=self.maxlen)
        item = {
            "input_ids": torch.tensor(e["input_ids"], dtype=torch.long),
            "attention_mask": torch.tensor(e["attention_mask"], dtype=torch.long),
        }
        if self.y is not None:
            item["y"] = torch.tensor(self.y[i], dtype=torch.float32)
        return item


def collate(batch, pad_id):
    m = max(len(b["input_ids"]) for b in batch)
    ids = torch.full((len(batch), m), pad_id, dtype=torch.long)
    am = torch.zeros((len(batch), m), dtype=torch.long)
    ys = []
    for j, b in enumerate(batch):
        L = len(b["input_ids"])
        ids[j, :L] = b["input_ids"]
        am[j, :L] = b["attention_mask"]
        if "y" in b:
            ys.append(b["y"])
    out = {"input_ids": ids, "attention_mask": am}
    if ys:
        out["y"] = torch.stack(ys)
    return out


class Regressor(nn.Module):
    def __init__(self, model_type, model_dir, drop=0.1):
        super().__init__()
        if model_type == "peptideclm":
            self.enc = RoFormerModel.from_pretrained(model_dir)
        else:
            self.enc = AutoModel.from_pretrained(model_dir, trust_remote_code=True)
        h = getattr(self.enc.config, "hidden_size", None) or getattr(self.enc.config, "embed_dim")
        self.head = nn.Sequential(
            nn.Dropout(drop),
            nn.Linear(h, 256),
            nn.GELU(),
            nn.Dropout(drop),
            nn.Linear(256, 1),
        )

    def forward(self, input_ids, attention_mask):
        o = self.enc(input_ids=input_ids, attention_mask=attention_mask).last_hidden_state
        mask = attention_mask.unsqueeze(-1).float()
        pooled = (o * mask).sum(1) / mask.sum(1).clamp(min=1e-6)
        return self.head(pooled).squeeze(-1)


def predict_loader(model, dl, dev, mu, sd):
    model.eval()
    P, Y = [], []
    with torch.no_grad():
        for b in dl:
            p = model(b["input_ids"].to(dev), b["attention_mask"].to(dev)).cpu().numpy()
            P.append(np.atleast_1d(p))
            if "y" in b:
                Y.append(b["y"].numpy())
    P = np.concatenate(P) * sd + mu
    if Y:
        Y = np.concatenate(Y) * sd + mu
        ss_res = float(((Y - P) ** 2).sum())
        ss_tot = float(((Y - Y.mean()) ** 2).sum())
        metrics = {
            "r2": float(1.0 - ss_res / ss_tot),
            "rmse": float(np.sqrt(np.mean((Y - P) ** 2))),
            "mae": float(np.mean(np.abs(Y - P))),
            "spearman": float(spearmanr(Y, P).correlation),
            "pearson": float(pearsonr(Y, P)[0]),
        }
        return metrics, P, Y
    return None, P, None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--split", required=True, choices=["random", "cluster", "source"])
    ap.add_argument("--model_type", default="peptideclm", choices=["peptideclm", "hf_auto"])
    ap.add_argument("--model_dir", default=str(PERM_DIR / "model"))
    ap.add_argument("--tag", default="peptideclm")
    ap.add_argument("--epochs", type=int, default=40)
    ap.add_argument("--bs", type=int, default=32)
    ap.add_argument("--lr", type=float, default=2e-5)
    ap.add_argument("--maxlen", type=int, default=512)
    ap.add_argument("--patience", type=int, default=6)
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()

    torch.manual_seed(args.seed)
    np.random.seed(args.seed)
    dev = "cuda" if torch.cuda.is_available() else "cpu"
    print(f"tag={args.tag} split={args.split} dev={dev}", flush=True)

    sdir = SPLIT_ROOT / args.split
    out_dir = SPLIT_ROOT / f"{args.tag}_{args.split}"
    out_dir.mkdir(parents=True, exist_ok=True)

    trn = pd.read_csv(sdir / "train.csv")
    va = pd.read_csv(sdir / "val.csv")
    te = pd.read_csv(sdir / "test.csv")
    mu = float(trn["PAMPA"].mean())
    sd = float(trn["PAMPA"].std())
    print(f"train={len(trn)} val={len(va)} test={len(te)} mu={mu:.3f} sd={sd:.3f}", flush=True)

    tok = load_tokenizer(args.model_type, args.model_dir)
    pad = tok.pad_token_id if tok.pad_token_id is not None else 0

    def mk(df_in, shuffle=False):
        ys = None if "PAMPA" not in df_in.columns else (df_in["PAMPA"].values - mu) / sd
        smis = df_in["SMILES"].tolist() if "SMILES" in df_in.columns else df_in["smiles"].tolist()
        return DataLoader(
            DS(smis, ys, tok, args.maxlen),
            batch_size=args.bs if shuffle else 64,
            shuffle=shuffle,
            collate_fn=lambda b: collate(b, pad),
        )

    dl_tr = mk(trn, shuffle=True)
    dl_va = mk(va, shuffle=False)
    dl_te = mk(te, shuffle=False)

    model = Regressor(args.model_type, args.model_dir).to(dev)
    opt = torch.optim.AdamW(model.parameters(), lr=args.lr, weight_decay=0.01)
    lossfn = nn.MSELoss()

    best_rho = -1e9
    bad = 0
    best_ep = 0
    ckpt_path = out_dir / "best_model.pt"

    for ep in range(1, args.epochs + 1):
        model.train()
        tl = 0.0
        for b in dl_tr:
            opt.zero_grad()
            p = model(b["input_ids"].to(dev), b["attention_mask"].to(dev))
            loss = lossfn(p, b["y"].to(dev))
            loss.backward()
            opt.step()
            tl += float(loss.item())
        vm, _, _ = predict_loader(model, dl_va, dev, mu, sd)
        print(
            f"ep{ep:02d} trloss={tl/len(dl_tr):.4f} val_r2={vm['r2']:.3f} val_rho={vm['spearman']:.3f} val_rmse={vm['rmse']:.3f}",
            flush=True,
        )
        if vm["spearman"] > best_rho:
            best_rho = vm["spearman"]
            best_ep = ep
            bad = 0
            torch.save({"model": model.state_dict(), "mu": mu, "sd": sd, "maxlen": args.maxlen}, ckpt_path)
        else:
            bad += 1
            if bad >= args.patience:
                print(f"Early stopping at epoch {ep} (best ep {best_ep}, val_rho={best_rho:.3f})", flush=True)
                break

    ck = torch.load(ckpt_path, map_location=dev)
    model.load_state_dict(ck["model"])
    tm, P, Y = predict_loader(model, dl_te, dev, mu, sd)
    tm["best_epoch"] = best_ep
    tm["best_val_spearman"] = float(best_rho)
    print("TEST METRICS:", json.dumps(tm, indent=2), flush=True)
    (out_dir / "test_metrics.json").write_text(json.dumps(tm, indent=2))
    pd.DataFrame({"smiles": te["SMILES"], "y_true": Y, "y_pred": P}).to_csv(out_dir / "test_predictions.csv", index=False)

    # Re-score Gemini 8,129 assembled macrocycles
    gem = pd.read_csv(PERM_DIR / "data" / "assembled_smiles.csv")
    _, gem_pred, _ = predict_loader(model, mk(gem, shuffle=False), dev, mu, sd)
    gem["perm"] = gem_pred
    gem.to_csv(out_dir / "gemini_generation_scores.csv", index=False)

    def grp(source, constraint=None):
        m = gem["source"] == source
        if constraint:
            m &= gem["constraint"] == constraint
        return gem.loc[m, "perm"].values

    rt_p = grp("restoken", "permeable")
    rt_u = grp("restoken", "unconstrained")
    rn_p = grp("random", "permeable")
    rn_a = grp("random")
    gem_contrasts = []
    for cname, a, b in [
        ("restoken_permeable > restoken_unconstrained", rt_p, rt_u),
        ("restoken_permeable > random_permeable", rt_p, rn_p),
        ("restoken_permeable > random_all", rt_p, rn_a),
    ]:
        u, pval = mannwhitneyu(a, b, alternative="greater")
        gem_contrasts.append({
            "contrast": cname,
            "n_a": int(len(a)),
            "n_b": int(len(b)),
            "median_a": float(np.median(a)),
            "median_b": float(np.median(b)),
            "delta_median": float(np.median(a) - np.median(b)),
            "cliffs_delta": float(cliffs_delta(a, b)),
            "p_one_sided": float(pval),
            "median_diff_95ci": boot_ci(a, b),
        })
    (out_dir / "gemini_contrasts.json").write_text(json.dumps(gem_contrasts, indent=2))
    print("GEMINI CONTRASTS:", json.dumps(gem_contrasts, indent=2), flush=True)

    # Re-score Open-Weight Models assembled macrocycles
    opn = pd.read_csv(PERM_DIR / "output" / "open_model_scores.csv")
    _, opn_pred, _ = predict_loader(model, mk(opn, shuffle=False), dev, mu, sd)
    opn["perm"] = opn_pred
    opn.to_csv(out_dir / "open_model_scores.csv", index=False)

    open_stats = []
    for mname in opn["model"].unique():
        a = opn[(opn["model"] == mname) & (opn["constraint"] == "permeable")]["perm"].values
        b = opn[(opn["model"] == mname) & (opn["constraint"] == "unconstrained")]["perm"].values
        rec = {
            "model": mname,
            "n_permeable": int(len(a)),
            "n_unconstrained": int(len(b)),
            "median_permeable": float(np.median(a)) if len(a) else None,
            "median_unconstrained": float(np.median(b)) if len(b) else None,
        }
        if len(a) >= 10 and len(b) >= 10:
            u, pval = mannwhitneyu(a, b, alternative="greater")
            rec.update({
                "delta_median": float(np.median(a) - np.median(b)),
                "cliffs_delta": float(cliffs_delta(a, b)),
                "p_one_sided": float(pval),
                "median_diff_95ci": boot_ci(a, b),
            })
        open_stats.append(rec)
    (out_dir / "open_model_contrasts.json").write_text(json.dumps(open_stats, indent=2))
    print("OPEN MODEL CONTRASTS:", json.dumps(open_stats, indent=2), flush=True)


if __name__ == "__main__":
    main()
