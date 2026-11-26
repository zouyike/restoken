"""Applicability-domain check for the PeptideCLM permeability judge.

Reviewer question: the judge was trained on CycPeptMPDB peptides, but ResToken
molecules are heavily N-methylated NCAAs — is the judge just extrapolating
out-of-distribution? We answer it structurally: for every assembled ResToken
molecule, compute the max Morgan-Tanimoto similarity to the 6,201 training
peptides (its nearest training neighbour), and compare that distribution to the
held-out TEST set's nearest-neighbour distribution. The test set is exactly
where the judge attains rho=0.74, so it defines the "in-domain" regime. If
ResToken's NN-similarity overlaps the test set's, the judge is operating in the
same chemical neighbourhood, not extrapolating.

Env: ai (rdkit). CPU only, no GPU.
"""
import numpy as np
import pandas as pd
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem
from rdkit import RDLogger
RDLogger.DisableLog("rdApp.*")

BASE = "/scratch/genesis/NCAA_tokenization/experiments/permeability_validation"
OUT = f"{BASE}/output/applicability_domain.json"

def fps(smiles_iter):
    out = []
    for s in smiles_iter:
        if not isinstance(s, str) or not s:
            continue
        m = Chem.MolFromSmiles(s)
        if m is None:
            continue
        out.append(AllChem.GetMorganFingerprintAsBitVect(m, 2, nBits=2048))
    return out


def nn_to_train(query_fps, train_fps):
    """Max Tanimoto of each query to the training set."""
    return np.array([max(DataStructs.BulkTanimotoSimilarity(q, train_fps))
                     for q in query_fps])


def summ(name, arr):
    q = np.percentile(arr, [5, 25, 50, 75, 95])
    return {
        "set": name, "n": int(len(arr)),
        "mean": float(arr.mean()),
        "p5": float(q[0]), "p25": float(q[1]), "median": float(q[2]),
        "p75": float(q[3]), "p95": float(q[4]),
        "frac_lt_0.30": float((arr < 0.30).mean()),  # clearly OOD
        "frac_lt_0.20": float((arr < 0.20).mean()),
    }


def main():
    train = pd.read_csv(f"{BASE}/data/pampa_train.csv")
    test = pd.read_csv(f"{BASE}/data/pampa_test.csv")
    gem = pd.read_csv(f"{BASE}/data/assembled_smiles.csv")
    opn = pd.read_csv(f"{BASE}/output/open_model_scores.csv")

    print("building train fps...", flush=True)
    train_fps = fps(train["SMILES"])
    print("  train fps:", len(train_fps), flush=True)

    print("test NN...", flush=True)
    test_nn = nn_to_train(fps(test["SMILES"]), train_fps)

    # ResToken: dedup assembled SMILES across both benchmarks
    rt_smiles = pd.concat([gem["smiles"], opn["smiles"]]).dropna().drop_duplicates()
    print("ResToken unique assembled:", len(rt_smiles), flush=True)
    rt_nn = nn_to_train(fps(rt_smiles), train_fps)

    # Also split ResToken by permeable vs the rest (sanity: same domain?)
    perm_smiles = pd.concat([
        gem.loc[gem["constraint"] == "permeable", "smiles"],
        opn.loc[opn["constraint"] == "permeable", "smiles"],
    ]).dropna().drop_duplicates()
    perm_nn = nn_to_train(fps(perm_smiles), train_fps)

    res = [summ("test (in-domain reference)", test_nn),
           summ("ResToken all", rt_nn),
           summ("ResToken permeable-only", perm_nn)]
    for r in res:
        print(r, flush=True)

    import json
    json.dump(res, open(OUT, "w"), indent=2)
    print("wrote", OUT)


if __name__ == "__main__":
    main()
