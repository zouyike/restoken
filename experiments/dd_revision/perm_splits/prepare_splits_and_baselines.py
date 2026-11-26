#!/usr/bin/env python3
"""CycPeptMPDB leakage-free splits, simple/trivial baselines, dataset/generation bias,
and ResToken 8-descriptor degeneracy analysis for Digital Discovery Major Revision.

Addresses:
- R2 (1b, 1c, 2a, 2b, 4b, 5a, 5b, 5c)
- R3 (2)
- R1 (4)
"""
import json
import math
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors
from rdkit.ML.Cluster import Butina
from scipy.stats import mannwhitneyu, pearsonr, spearmanr
from sklearn.ensemble import RandomForestRegressor
from sklearn.linear_model import Ridge

RDLogger.DisableLog("rdApp.*")

ROOT = Path("/scratch/genesis/NCAA_tokenization")
PERM_DIR = ROOT / "experiments" / "permeability_validation"
OUT_DIR = ROOT / "experiments" / "dd_revision" / "perm_splits"
OUT_DIR.mkdir(parents=True, exist_ok=True)


def cliffs_delta(a, b):
    a = np.asarray(a)
    b = np.asarray(b)
    # Vectorized chunked to avoid memory blowup
    gt = 0
    lt = 0
    for i in range(0, len(a), 500):
        chunk = a[i : i + 500, None]
        gt += int((chunk > b[None, :]).sum())
        lt += int((chunk < b[None, :]).sum())
    return (gt - lt) / (len(a) * len(b))


def boot_median_diff(a, b, n=3000, seed=0):
    rng = np.random.RandomState(seed)
    a = np.asarray(a)
    b = np.asarray(b)
    diffs = [np.median(rng.choice(a, len(a))) - np.median(rng.choice(b, len(b))) for _ in range(n)]
    return [float(np.percentile(diffs, 2.5)), float(np.percentile(diffs, 97.5))]


def reg_metrics(y_true, y_pred):
    y_true = np.asarray(y_true, dtype=float)
    y_pred = np.asarray(y_pred, dtype=float)
    ss_res = float(((y_true - y_pred) ** 2).sum())
    ss_tot = float(((y_true - y_true.mean()) ** 2).sum())
    return {
        "r2": float(1.0 - ss_res / ss_tot) if ss_tot > 0 else 0.0,
        "rmse": float(np.sqrt(np.mean((y_true - y_pred) ** 2))),
        "mae": float(np.mean(np.abs(y_true - y_pred))),
        "spearman": float(spearmanr(y_true, y_pred).correlation),
        "pearson": float(pearsonr(y_true, y_pred)[0]),
    }


def compute_fp_and_desc(smiles_list):
    fps_bv = []
    fps_arr = np.zeros((len(smiles_list), 2048), dtype=np.float32)
    descs = np.zeros((len(smiles_list), 8), dtype=np.float32)
    can_iso = []
    can_flat = []
    for i, s in enumerate(smiles_list):
        m = Chem.MolFromSmiles(s)
        bv = AllChem.GetMorganFingerprintAsBitVect(m, 2, nBits=2048)
        fps_bv.append(bv)
        arr = np.zeros((2048,), dtype=np.int8)
        DataStructs.ConvertToNumpyArray(bv, arr)
        fps_arr[i] = arr
        descs[i] = [
            Descriptors.MolWt(m),
            Descriptors.MolLogP(m),
            Descriptors.TPSA(m),
            rdMolDescriptors.CalcNumHBD(m),
            rdMolDescriptors.CalcNumHBA(m),
            rdMolDescriptors.CalcNumRotatableBonds(m),
            rdMolDescriptors.CalcNumRings(m),
            rdMolDescriptors.CalcFractionCSP3(m),
        ]
        can_iso.append(Chem.MolToSmiles(m, isomericSmiles=True))
        can_flat.append(Chem.MolToSmiles(m, isomericSmiles=False))
    return fps_bv, fps_arr, descs, can_iso, can_flat


def build_butina_clusters(fps_bv, cutoff=0.6):
    """Butina clustering using distance = 1 - Tanimoto <= 1 - cutoff."""
    n = len(fps_bv)
    dists = []
    for i in range(1, n):
        sims = DataStructs.BulkTanimotoSimilarity(fps_bv[i], fps_bv[:i])
        dists.extend([1.0 - x for x in sims])
    clusters = Butina.ClusterData(dists, n, 1.0 - cutoff, isDistData=True)
    return clusters


def assign_groups_to_splits(groups, n_total, target_val=0.10, target_test=0.10, seed=42):
    """Assign list of index-lists (groups) into train/val/test without splitting any group.
    Any mega-group larger than 40% of the val/test budget goes directly to train so val/test
    contain multiple independent clusters/sources and stay close to 80/10/10.
    """
    rng = np.random.RandomState(seed)
    n_test_target = int(round(n_total * target_test))
    n_val_target = int(round(n_total * target_val))
    max_holdout_group = max(250, int(0.45 * n_test_target))

    train_idx, val_idx, test_idx = [], [], []
    small_groups = []
    for g in groups:
        if len(g) > max_holdout_group:
            train_idx.extend(g)
        else:
            small_groups.append(list(g))

    order = list(range(len(small_groups)))
    rng.shuffle(order)
    for g_i in order:
        g = small_groups[g_i]
        if len(test_idx) < n_test_target and (len(test_idx) + len(g) <= n_test_target * 1.15 or len(test_idx) < n_test_target * 0.7):
            test_idx.extend(g)
        elif len(val_idx) < n_val_target and (len(val_idx) + len(g) <= n_val_target * 1.15 or len(val_idx) < n_val_target * 0.7):
            val_idx.extend(g)
        else:
            train_idx.extend(g)
    return sorted(train_idx), sorted(val_idx), sorted(test_idx)


def verify_split(df, tr_idx, va_idx, te_idx, fps_bv, can_iso, can_flat, name):
    tr_iso = {can_iso[i] for i in tr_idx}
    va_iso = {can_iso[i] for i in va_idx}
    te_iso = {can_iso[i] for i in te_idx}
    tr_flat = {can_flat[i] for i in tr_idx}
    va_flat = {can_flat[i] for i in va_idx}
    te_flat = {can_flat[i] for i in te_idx}
    tr_src = set(df.iloc[tr_idx]["Source"].dropna())
    te_src = set(df.iloc[te_idx]["Source"].dropna())
    tr_fps = [fps_bv[i] for i in tr_idx]
    te_nn = [max(DataStructs.BulkTanimotoSimilarity(fps_bv[i], tr_fps)) for i in te_idx]
    return {
        "split": name,
        "n_train": len(tr_idx),
        "n_val": len(va_idx),
        "n_test": len(te_idx),
        "iso_overlap_tr_te": len(tr_iso & te_iso),
        "iso_overlap_tr_va": len(tr_iso & va_iso),
        "flat_overlap_tr_te": len(tr_flat & te_flat),
        "flat_overlap_tr_va": len(tr_flat & va_flat),
        "source_overlap_tr_te": len(tr_src & te_src),
        "test_nn_tanimoto_median": float(np.median(te_nn)),
        "test_nn_tanimoto_mean": float(np.mean(te_nn)),
        "test_nn_tanimoto_max": float(np.max(te_nn)),
        "test_nn_tanimoto_p95": float(np.percentile(te_nn, 95)),
    }


def analyze_library_degeneracy():
    llm_d = json.loads((ROOT / "restoken" / "data" / "bb_dict_llm_v11.json").read_text())
    back_d = json.loads((ROOT / "restoken" / "data" / "bb_dict_backend_v11.json").read_text())
    back_by_id = {b["id"]: b for b in back_d["blocks"]}
    blocks = llm_d["blocks"]

    nmod_counts = Counter(b["mc_nmod"] for b in blocks)
    case_nmod = Counter((b["id"][0].islower(), b["mc_nmod"]) for b in blocks)

    clusters = defaultdict(list)
    for b in blocks:
        prof = (
            b["class"], b["chirality"], b["charge_label"],
            b["mc_type"], b["mc_nmod"], b["sc_bulk"],
            b["polarity_bin"], b["flex_bin"],
        )
        clusters[prof].append(b["id"])

    # Compute Morgan fingerprints of all 400 blocks (AA-SMILES)
    fp_map = {}
    for bid, b in back_by_id.items():
        smi = b.get("structure", {}).get("aa_smiles") or b.get("aa_smiles") or b.get("smiles")
        if not smi:
            continue
        m = Chem.MolFromSmiles(smi)
        if m is not None:
            fp_map[bid] = AllChem.GetMorganFingerprintAsBitVect(m, 2, nBits=2048)

    intra_sims = []
    cluster_rows = []
    for prof, ids in sorted(clusters.items(), key=lambda kv: -len(kv[1])):
        if len(ids) > 1:
            sims = []
            for i in range(len(ids)):
                for j in range(i + 1, len(ids)):
                    if ids[i] in fp_map and ids[j] in fp_map:
                        s = DataStructs.TanimotoSimilarity(fp_map[ids[i]], fp_map[ids[j]])
                        sims.append(s)
                        intra_sims.append(s)
            cluster_rows.append({
                "profile": "|".join(str(x) for x in prof),
                "size": len(ids),
                "ids_sample": ids[:6],
                "mean_tanimoto": float(np.mean(sims)) if sims else None,
                "min_tanimoto": float(np.min(sims)) if sims else None,
            })

    # Random background pairwise Tanimoto across 400 blocks
    all_ids = sorted(fp_map.keys())
    rng = np.random.RandomState(42)
    bg_sims = []
    for _ in range(5000):
        i, j = rng.choice(len(all_ids), 2, replace=False)
        bg_sims.append(DataStructs.TanimotoSimilarity(fp_map[all_ids[i]], fp_map[all_ids[j]]))

    out = {
        "n_blocks": len(blocks),
        "n_mod_counts": dict(nmod_counts),
        "lowercase_by_nmod": {f"lower={k[0]}|{k[1]}": v for k, v in case_nmod.items()},
        "n_unique_8prop_profiles": len(clusters),
        "n_singleton_profiles": sum(1 for v in clusters.values() if len(v) == 1),
        "n_colliding_profiles": sum(1 for v in clusters.values() if len(v) > 1),
        "n_blocks_in_colliding_profiles": sum(len(v) for v in clusters.values() if len(v) > 1),
        "frac_blocks_in_colliding_profiles": sum(len(v) for v in clusters.values() if len(v) > 1) / len(blocks),
        "max_cluster_size": max(len(v) for v in clusters.values()),
        "intra_cluster_tanimoto_mean": float(np.mean(intra_sims)),
        "intra_cluster_tanimoto_median": float(np.median(intra_sims)),
        "intra_cluster_tanimoto_p10": float(np.percentile(intra_sims, 10)),
        "intra_cluster_tanimoto_p90": float(np.percentile(intra_sims, 90)),
        "background_pairwise_tanimoto_mean": float(np.mean(bg_sims)),
        "background_pairwise_tanimoto_median": float(np.median(bg_sims)),
        "top_10_colliding_clusters": cluster_rows[:10],
    }
    (OUT_DIR / "library_degeneracy_stats.json").write_text(json.dumps(out, indent=2))
    return out


def analyze_generation_bias():
    """Analyze monomer usage bias, entropy, and library skew in generated sequences (R2-1c)."""
    llm_d = json.loads((ROOT / "restoken" / "data" / "bb_dict_llm_v11.json").read_text())
    lib_blocks = {b["id"]: b for b in llm_d["blocks"]}
    n_lib = len(lib_blocks)
    max_ent = math.log2(n_lib)

    rows = []
    # 1. Gemini / GPT-4o / Qwen / Gemma Task 1 exp1_validity
    exp1_dir = ROOT / "experiments" / "output" / "exp1_validity"
    for p in sorted(exp1_dir.glob("*_restoken_unconstrained_len6_results.json")):
        model = p.name.split("_restoken_")[0]
        data = json.loads(p.read_text())
        seqs = [
            (r.get("tokens") or r.get("raw", "").split("-"))
            for r in data
            if (r.get("overall_valid") or r.get("valid")) and (r.get("tokens") or r.get("raw"))
        ]
        if not seqs:
            continue
        toks = [t for s in seqs for t in s if t in lib_blocks]
        c = Counter(toks)
        total = len(toks)
        probs = np.array([v / total for v in c.values()])
        ent = float(-(probs * np.log2(probs)).sum())
        top10_share = sum(v for _, v in c.most_common(10)) / total
        mc_dist = Counter(lib_blocks[t]["mc_type"] for t in toks)
        rows.append({
            "suite": "Task1_unconstrained",
            "model": model,
            "n_valid_seqs": len(seqs),
            "n_unique_seqs": len({tuple(s) for s in seqs}),
            "n_unique_blocks_used": len(c),
            "library_coverage_pct": 100.0 * len(c) / n_lib,
            "shannon_entropy_bits": round(ent, 3),
            "normalized_entropy": round(ent / max_ent, 3),
            "top10_block_share_pct": round(100.0 * top10_share, 2),
            "alpha_pct": round(100.0 * mc_dist["alpha"] / total, 1),
            "beta_pct": round(100.0 * mc_dist["beta"] / total, 1),
            "gamma_pct": round(100.0 * mc_dist["gamma"] / total, 1),
            "top5_blocks": c.most_common(5),
        })

    # 2. Gemini assembled_smiles.csv (unconstrained vs permeable vs random)
    gem = pd.read_csv(PERM_DIR / "data" / "assembled_smiles.csv")
    for (src, cons), g in gem.groupby(["source", "constraint"]):
        seqs = [s.split("-") for s in g["token_seq"].dropna() if isinstance(s, str)]
        toks = [t for s in seqs for t in s if t in lib_blocks]
        if not toks:
            continue
        c = Counter(toks)
        total = len(toks)
        probs = np.array([v / total for v in c.values()])
        ent = float(-(probs * np.log2(probs)).sum())
        top10_share = sum(v for _, v in c.most_common(10)) / total
        mc_dist = Counter(lib_blocks[t]["mc_type"] for t in toks)
        rows.append({
            "suite": f"perm_eval_{src}",
            "model": cons,
            "n_valid_seqs": len(seqs),
            "n_unique_seqs": len({tuple(s) for s in seqs}),
            "n_unique_blocks_used": len(c),
            "library_coverage_pct": 100.0 * len(c) / n_lib,
            "shannon_entropy_bits": round(ent, 3),
            "normalized_entropy": round(ent / max_ent, 3),
            "top10_block_share_pct": round(100.0 * top10_share, 2),
            "alpha_pct": round(100.0 * mc_dist["alpha"] / total, 1),
            "beta_pct": round(100.0 * mc_dist["beta"] / total, 1),
            "gamma_pct": round(100.0 * mc_dist["gamma"] / total, 1),
            "top5_blocks": c.most_common(5),
        })

    (OUT_DIR / "generation_bias_stats.json").write_text(json.dumps(rows, indent=2))
    return rows


def main():
    print("=== 1. Library degeneracy & N-mod check ===", flush=True)
    deg = analyze_library_degeneracy()
    print(f"Unique 8-prop profiles: {deg['n_unique_8prop_profiles']}/400, "
          f"colliding blocks: {deg['n_blocks_in_colliding_profiles']} ({deg['frac_blocks_in_colliding_profiles']:.1%}), "
          f"intra-cluster Tanimoto mean={deg['intra_cluster_tanimoto_mean']:.3f} vs bg={deg['background_pairwise_tanimoto_mean']:.3f}", flush=True)

    print("\n=== 2. Generation bias check ===", flush=True)
    gen_bias = analyze_generation_bias()
    print(f"Computed bias metrics for {len(gen_bias)} model/condition groups.", flush=True)

    print("\n=== 3. CycPeptMPDB funnel & bias ===", flush=True)
    raw = pd.read_csv(PERM_DIR / "data" / "CycPeptMPDB_Peptide_All.csv", low_memory=False)
    n_total = len(raw)
    raw["PAMPA"] = pd.to_numeric(raw["PAMPA"], errors="coerce")
    d1 = raw[raw["PAMPA"].notna()].copy()
    n_pampa = len(d1)
    d2 = d1.dropna(subset=["SMILES"]).copy()
    n_nonnull_smi = len(d2)
    d3 = d2.drop_duplicates(subset=["SMILES"]).copy()
    n_dedup_smi = len(d3)
    ok_mask = d3["SMILES"].map(lambda s: Chem.MolFromSmiles(s) is not None)
    d = d3[ok_mask].reset_index(drop=True)
    n_rdkit = len(d)

    src_counts = d["Source"].value_counts()
    funnel = {
        "database": "CycPeptMPDB v1.2",
        "access_date": "2026-06-01",
        "year_range": [int(raw["Year"].min()), int(raw["Year"].max())],
        "n_total_raw": n_total,
        "n_removed_missing_pampa": n_total - n_pampa,
        "n_with_pampa": n_pampa,
        "n_removed_missing_smiles": n_pampa - n_nonnull_smi,
        "n_removed_duplicate_smiles": n_nonnull_smi - n_dedup_smi,
        "n_removed_rdkit_invalid": n_dedup_smi - n_rdkit,
        "n_final_clean": n_rdkit,
        "n_unique_sources": int(d["Source"].nunique()),
        "top_5_sources": {str(k): int(v) for k, v in src_counts.head(5).items()},
        "top_3_sources_share_pct": round(100.0 * src_counts.head(3).sum() / n_rdkit, 2),
        "pampa_mean": float(d["PAMPA"].mean()),
        "pampa_std": float(d["PAMPA"].std()),
        "pampa_median": float(d["PAMPA"].median()),
        "monomer_length_range": [int(d["Monomer_Length"].min()), int(d["Monomer_Length"].max())],
    }
    print(json.dumps(funnel, indent=2), flush=True)

    print("\n=== 4. Fingerprints, descriptors, and leakage-free splits ===", flush=True)
    fps_bv, fps_arr, descs, can_iso, can_flat = compute_fp_and_desc(d["SMILES"].tolist())
    d["can_iso"] = can_iso
    d["can_flat"] = can_flat
    funnel["n_unique_canonical_isomeric"] = len(set(can_iso))
    funnel["n_unique_canonical_flat"] = len(set(can_flat))

    # Split A: exact original Random split (seed=42 for 90/10, seed=0 for 90/10 val)
    rng42 = np.random.RandomState(42)
    idx42 = rng42.permutation(len(d))
    ntest = int(len(d) * 0.1)
    te_rand = sorted(idx42[:ntest].tolist())
    trva_rand = idx42[ntest:]
    rng0 = np.random.RandomState(0)
    idx0 = rng0.permutation(len(trva_rand))
    nval = int(len(trva_rand) * 0.1)
    va_rand = sorted(trva_rand[idx0[:nval]].tolist())
    tr_rand = sorted(trva_rand[idx0[nval:]].tolist())

    # Split B: Butina Cluster split (Tanimoto cutoff = 0.60, merged by flat graph to guarantee 0 stereoisomer leakage)
    print("Running Butina clustering (cutoff=0.60) on 6,889 peptides...", flush=True)
    clusters = build_butina_clusters(fps_bv, cutoff=0.60)
    funnel["n_butina_clusters_0p6"] = len(clusters)
    funnel["top5_butina_cluster_sizes"] = [len(c) for c in clusters[:5]]
    funnel["singleton_clusters_0p6"] = sum(1 for c in clusters if len(c) == 1)

    # Union-Find to merge any clusters sharing the same flat-graph SMILES
    parent = list(range(len(clusters)))
    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[rb] = ra

    flat_to_cid = {}
    for cid, members in enumerate(clusters):
        for m_idx in members:
            fg = can_flat[m_idx]
            if fg in flat_to_cid:
                union(cid, flat_to_cid[fg])
            else:
                flat_to_cid[fg] = cid
    merged = defaultdict(list)
    for cid, members in enumerate(clusters):
        merged[find(cid)].extend(members)
    cluster_groups = list(merged.values())
    tr_clu, va_clu, te_clu = assign_groups_to_splits(cluster_groups, len(d), 0.10, 0.10, seed=42)

    # Split C: Source / Literature split (disjoint studies; any flat-graph shared across sources assigned to train)
    flat_sources = defaultdict(set)
    for i, (fg, src) in enumerate(zip(can_flat, d["Source"].fillna("UNKNOWN"))):
        flat_sources[fg].add(src)
    # Build source groups; merge sources that share any flat graph? Better: hold out pure-source studies
    # with zero flat-graph overlap with remaining sources!
    src_to_idxs = defaultdict(list)
    for i, src in enumerate(d["Source"].fillna("UNKNOWN")):
        src_to_idxs[src].append(i)
    # Union sources that share a flat graph so zero flat-graph leakage is guaranteed
    src_list = sorted(src_to_idxs.keys())
    src_id = {s: i for i, s in enumerate(src_list)}
    sp = list(range(len(src_list)))
    def sfind(x):
        while sp[x] != x:
            sp[x] = sp[sp[x]]
            x = sp[x]
        return x
    def sunion(a, b):
        ra, rb = sfind(a), sfind(b)
        if ra != rb:
            sp[rb] = ra
    for fg, srcs in flat_sources.items():
        if len(srcs) > 1:
            s0 = next(iter(srcs))
            for s in srcs:
                sunion(src_id[s0], src_id[s])
    src_merged = defaultdict(list)
    for s, idxs in src_to_idxs.items():
        src_merged[sfind(src_id[s])].extend(idxs)
    source_groups = list(src_merged.values())
    tr_src, va_src, te_src = assign_groups_to_splits(source_groups, len(d), 0.10, 0.10, seed=7)

    splits = {
        "random": (tr_rand, va_rand, te_rand),
        "cluster": (tr_clu, va_clu, te_clu),
        "source": (tr_src, va_src, te_src),
    }

    split_reports = {}
    for sname, (tr_i, va_i, te_i) in splits.items():
        sdir = OUT_DIR / sname
        sdir.mkdir(exist_ok=True)
        cols = ["CycPeptMPDB_ID", "Source", "Year", "SMILES", "PAMPA", "Monomer_Length"]
        d.iloc[tr_i][cols].to_csv(sdir / "train.csv", index=False)
        d.iloc[va_i][cols].to_csv(sdir / "val.csv", index=False)
        d.iloc[te_i][cols].to_csv(sdir / "test.csv", index=False)
        rep = verify_split(d, tr_i, va_i, te_i, fps_bv, can_iso, can_flat, sname)
        split_reports[sname] = rep
        print(f"Split [{sname}]:", json.dumps(rep), flush=True)

    (OUT_DIR / "cycpeptmpdb_funnel_and_splits.json").write_text(
        json.dumps({"funnel": funnel, "splits": split_reports}, indent=2)
    )

    print("\n=== 5. Simple/Trivial Baselines across all 3 splits (R2-4b) ===", flush=True)
    y_all = d["PAMPA"].values.astype(np.float32)

    # Precompute assembled macrocycle fingerprints for re-scoring with RF baselines
    gem = pd.read_csv(PERM_DIR / "data" / "assembled_smiles.csv")
    gem_fps = np.zeros((len(gem), 2048), dtype=np.float32)
    for i, s in enumerate(gem["smiles"]):
        m = Chem.MolFromSmiles(s)
        if m is not None:
            bv = AllChem.GetMorganFingerprintAsBitVect(m, 2, nBits=2048)
            arr = np.zeros((2048,), dtype=np.int8)
            DataStructs.ConvertToNumpyArray(bv, arr)
            gem_fps[i] = arr

    baseline_results = {}
    rf_perm_contrasts = {}

    for sname, (tr_i, va_i, te_i) in splits.items():
        trva_i = tr_i + va_i
        X_tr, y_tr = fps_arr[trva_i], y_all[trva_i]
        X_te, y_te = fps_arr[te_i], y_all[te_i]
        D_tr, D_te = descs[trva_i], descs[te_i]

        res = {}
        # 1. Mean baseline
        res["mean_baseline"] = reg_metrics(y_te, np.full_like(y_te, y_tr.mean()))

        # 2. 1-NN Tanimoto
        tr_bvs = [fps_bv[i] for i in trva_i]
        nn_pred = []
        for i in te_i:
            sims = DataStructs.BulkTanimotoSimilarity(fps_bv[i], tr_bvs)
            nn_pred.append(y_tr[int(np.argmax(sims))])
        res["tanimoto_1nn"] = reg_metrics(y_te, np.array(nn_pred))

        # 3. ECFP4 + Ridge (tune alpha on val)
        best_a, best_rho = 10.0, -9.0
        for alpha in [1.0, 10.0, 50.0, 100.0, 500.0]:
            m_r = Ridge(alpha=alpha).fit(fps_arr[tr_i], y_all[tr_i])
            p_v = m_r.predict(fps_arr[va_i])
            rho = float(spearmanr(y_all[va_i], p_v).correlation)
            if rho > best_rho:
                best_rho, best_a = rho, alpha
        m_ridge = Ridge(alpha=best_a).fit(X_tr, y_tr)
        res["ecfp4_ridge"] = {**reg_metrics(y_te, m_ridge.predict(X_te)), "alpha": best_a}

        # 4. RDKit 2D descriptors + Random Forest
        m_desc_rf = RandomForestRegressor(n_estimators=200, random_state=42, n_jobs=-1)
        m_desc_rf.fit(D_tr, y_tr)
        res["rdkit2d_rf"] = reg_metrics(y_te, m_desc_rf.predict(D_te))

        # 5. ECFP4 + Random Forest
        m_rf = RandomForestRegressor(n_estimators=200, max_features="sqrt", random_state=42, n_jobs=-1)
        m_rf.fit(X_tr, y_tr)
        res["ecfp4_rf"] = reg_metrics(y_te, m_rf.predict(X_te))

        baseline_results[sname] = res
        print(f"\nBaselines [{sname}]:", flush=True)
        for bname, bmet in res.items():
            print(f"  {bname:16s} rho={bmet['spearman']:.3f} r={bmet['pearson']:.3f} R2={bmet['r2']:.3f} RMSE={bmet['rmse']:.3f}", flush=True)

        # Re-score Gemini 8,129 assembled macrocycles with this split's ECFP4-RF judge!
        preds = m_rf.predict(gem_fps)
        def grp(source, constraint=None):
            mask = (gem["source"] == source)
            if constraint:
                mask &= (gem["constraint"] == constraint)
            return preds[mask.values]
        rt_p = grp("restoken", "permeable")
        rt_u = grp("restoken", "unconstrained")
        rn_p = grp("random", "permeable")
        rn_a = grp("random")
        contrasts = []
        for cname, a, b in [
            ("restoken_permeable > restoken_unconstrained", rt_p, rt_u),
            ("restoken_permeable > random_permeable", rt_p, rn_p),
            ("restoken_permeable > random_all", rt_p, rn_a),
        ]:
            u, pval = mannwhitneyu(a, b, alternative="greater")
            contrasts.append({
                "contrast": cname,
                "n_a": int(len(a)),
                "n_b": int(len(b)),
                "median_a": float(np.median(a)),
                "median_b": float(np.median(b)),
                "delta_median": float(np.median(a) - np.median(b)),
                "cliffs_delta": float(cliffs_delta(a, b)),
                "p_one_sided": float(pval),
                "median_diff_95ci": boot_median_diff(a, b),
            })
        rf_perm_contrasts[sname] = contrasts

    (OUT_DIR / "simple_baselines_results.json").write_text(
        json.dumps({"test_metrics": baseline_results, "rf_generation_contrasts": rf_perm_contrasts}, indent=2)
    )
    print("\nALL DONE! Outputs in:", OUT_DIR, flush=True)


if __name__ == "__main__":
    main()
