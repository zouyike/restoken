import os, sys, json, glob, csv

WD = "/scratch/genesis/NCAA_tokenization/experiments/permeability_validation"
W3 = "/scratch/genesis/NCAA_tokenization/experiments/output/w3_scaleup"
RB = "/scratch/genesis/NCAA_tokenization/experiments/output/random_baseline"
sys.path.insert(0, "/scratch/genesis/NCAA_tokenization/restoken/src")
from cyclic_assembler import CyclicPeptideAssembler


def restoken_seqs(path):
    """raw.json is list[str]; each str is newline-separated dash-token sequences."""
    out = []
    for block in json.load(open(path)):
        for line in str(block).splitlines():
            line = line.strip()
            if line:
                out.append(line)
    return out


def random_seqs(path):
    out = []
    with open(path) as f:
        for row in csv.DictReader(f):
            s = (row.get("sequence") or "").strip()
            if s:
                out.append(s)
    return out


def main():
    asm = CyclicPeptideAssembler()
    rows = []          # (source, model, constraint, token_seq, smiles)
    seen = set()       # dedup on (source, smiles)

    # ResToken generations (both models, all constraints)
    for path in sorted(glob.glob(f"{W3}/*_restoken_*_len6_raw.json")):
        base = os.path.basename(path)
        model = "pro" if "pro" in base else "flash"
        constraint = base.split("_restoken_")[1].split("_len6")[0]
        for seq in restoken_seqs(path):
            try:
                smi = asm.assemble(seq)
            except Exception:
                smi = None
            if not smi:
                continue
            k = ("restoken", smi)
            if k in seen:
                continue
            seen.add(k)
            rows.append(("restoken", model, constraint, seq, smi))

    # Random baselines
    for path in sorted(glob.glob(f"{RB}/baseline_*_len6_n1000.csv")):
        base = os.path.basename(path)
        constraint = base.split("baseline_")[1].split("_len6")[0]
        for seq in random_seqs(path):
            try:
                smi = asm.assemble(seq)
            except Exception:
                smi = None
            if not smi:
                continue
            k = ("random", smi)
            if k in seen:
                continue
            seen.add(k)
            rows.append(("random", "na", constraint, seq, smi))

    out = f"{WD}/data/assembled_smiles.csv"
    with open(out, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["source", "model", "constraint", "token_seq", "smiles"])
        w.writerows(rows)

    # summary
    from collections import Counter
    c = Counter((r[0], r[2]) for r in rows)
    print(f"wrote {len(rows)} unique assembled SMILES -> {out}")
    for (src, con), n in sorted(c.items()):
        print(f"  {src:9s} {con:16s} {n}")


if __name__ == "__main__":
    main()
