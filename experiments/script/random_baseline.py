"""Random baseline generator for ResToken benchmark.

Rejection sampling: random sequence → validate → keep/reject.
Establishes the baseline for "how hard is it to randomly generate valid sequences?"

Usage:
    python random_baseline.py --n 1000 --length 6 --profile unconstrained
    python random_baseline.py --n 1000 --length 6 --profile permeable
    python random_baseline.py --n 1000 --length 8 --profile charged_binder
"""

import argparse
import csv
import json
import random
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from restoken.src.library import BlockLibrary
from restoken.src.validator import SequenceValidator
from restoken.src.reconstructor import SMILESReconstructor

CONSTRAINT_PROFILES = {
    "unconstrained": {},
    "permeable": {
        "target_charge": 0,
        "max_hbd": 1,
        "max_rot": None,
        "custom": lambda blocks: (
            sum(1 for b in blocks if b.mc_nmod in ("NME", "NCY")) >= 2
            and sum(1 for b in blocks if b.sc_bulk == "large") >= 3
        ),
    },
    "charged_binder": {
        "target_charge": 2,
        "custom": lambda blocks: (
            sum(1 for b in blocks if b.sc_hbd > 0) >= 2  # proxy for total HBD >= 2
            and any(b.aa_class in ("F", "Y", "W", "H") for b in blocks)
        ),
    },
    "rigid_scaffold": {
        "max_rot": 18,
        "custom": lambda blocks: (
            sum(1 for b in blocks if b.mc_type == "beta") >= 3
        ),
    },
}


def generate_random_sequence(lib: BlockLibrary, length: int) -> list[str]:
    all_ids = list(lib.all_ids)
    return [random.choice(all_ids) for _ in range(length)]


def check_custom(lib: BlockLibrary, seq: list[str], custom_fn) -> bool:
    blocks = [lib[tid] for tid in seq]
    return custom_fn(blocks)


def main():
    parser = argparse.ArgumentParser(description="Random baseline generator")
    parser.add_argument("--n", type=int, default=1000, help="Number of valid sequences to generate")
    parser.add_argument("--length", type=int, default=6)
    parser.add_argument("--profile", default="unconstrained",
                        choices=list(CONSTRAINT_PROFILES.keys()))
    parser.add_argument("--max-attempts", type=int, default=1_000_000)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--output-dir", default=None)
    args = parser.parse_args()

    random.seed(args.seed)
    lib = BlockLibrary()

    profile = CONSTRAINT_PROFILES[args.profile]
    custom_fn = profile.pop("custom", None)

    sv = SequenceValidator(
        library=lib,
        allowed_lengths=[args.length],
        target_charge=profile.get("target_charge"),
        max_hbd=profile.get("max_hbd"),
        max_hba=profile.get("max_hba"),
        max_rot=profile.get("max_rot"),
    )

    if args.output_dir:
        out_dir = Path(args.output_dir)
    else:
        out_dir = Path(__file__).resolve().parent.parent / "output" / "random_baseline"
    out_dir.mkdir(parents=True, exist_ok=True)

    valid_sequences = []
    attempts = 0
    t0 = time.time()

    while len(valid_sequences) < args.n and attempts < args.max_attempts:
        seq = generate_random_sequence(lib, args.length)
        attempts += 1

        result = sv.validate(seq)
        if not result.valid:
            continue

        if custom_fn and not check_custom(lib, seq, custom_fn):
            continue

        valid_sequences.append({
            "sequence": "-".join(seq),
            "tokens": seq,
            **result.properties,
        })

    elapsed = time.time() - t0

    # Compute diversity (if rdkit available)
    diversity_stats = {}
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
        from rdkit import DataStructs
        import numpy as np

        recon = SMILESReconstructor(library=lib)
        fps = []
        for entry in valid_sequences:
            smiles_list = recon.residue_smiles(entry["tokens"])
            combined = ".".join(smiles_list)
            mol = Chem.MolFromSmiles(combined)
            if mol:
                fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)
                fps.append(fp)

        if len(fps) >= 2:
            tanimotos = []
            for i in range(min(len(fps), 500)):
                for j in range(i + 1, min(len(fps), 500)):
                    tanimotos.append(DataStructs.TanimotoSimilarity(fps[i], fps[j]))
            diversity_stats = {
                "mean_tanimoto": float(np.mean(tanimotos)),
                "std_tanimoto": float(np.std(tanimotos)),
                "mean_distance": float(1 - np.mean(tanimotos)),
                "n_fps": len(fps),
            }
    except ImportError:
        pass

    # Uniqueness
    unique_seqs = set(e["sequence"] for e in valid_sequences)

    # Report
    print(f"Profile: {args.profile}")
    print(f"Target: {args.n} valid sequences of length {args.length}")
    print(f"Generated: {len(valid_sequences)} valid / {attempts} attempts")
    print(f"Acceptance rate: {len(valid_sequences)/attempts*100:.2f}%")
    print(f"Unique: {len(unique_seqs)}/{len(valid_sequences)}")
    print(f"Time: {elapsed:.1f}s")

    if valid_sequences:
        charges = [e["net_charge"] for e in valid_sequences]
        hbds = [e["total_hbd"] for e in valid_sequences]
        hbas = [e["total_hba"] for e in valid_sequences]
        rots = [e["total_rot"] for e in valid_sequences]
        print(f"\nProperty distributions (mean ± std):")
        print(f"  Charge: {sum(charges)/len(charges):.2f}")
        print(f"  HBD: {sum(hbds)/len(hbds):.2f}")
        print(f"  HBA: {sum(hbas)/len(hbas):.2f}")
        print(f"  Rot: {sum(rots)/len(rots):.2f}")

    if diversity_stats:
        print(f"\nDiversity (Morgan FP Tanimoto):")
        print(f"  Mean similarity: {diversity_stats['mean_tanimoto']:.4f}")
        print(f"  Mean distance: {diversity_stats['mean_distance']:.4f}")

    # Save sequences
    csv_path = out_dir / f"baseline_{args.profile}_len{args.length}_n{len(valid_sequences)}.csv"
    with open(csv_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["sequence", "net_charge", "total_hbd", "total_hba", "total_rot"])
        writer.writeheader()
        for e in valid_sequences:
            writer.writerow({
                "sequence": e["sequence"],
                "net_charge": e["net_charge"],
                "total_hbd": e["total_hbd"],
                "total_hba": e["total_hba"],
                "total_rot": e["total_rot"],
            })

    # Save summary stats
    stats = {
        "profile": args.profile,
        "length": args.length,
        "target_n": args.n,
        "generated": len(valid_sequences),
        "attempts": attempts,
        "acceptance_rate": len(valid_sequences) / attempts,
        "unique": len(unique_seqs),
        "elapsed_s": elapsed,
        "diversity": diversity_stats,
    }
    json_path = out_dir / f"baseline_{args.profile}_len{args.length}_stats.json"
    with open(json_path, "w") as f:
        json.dump(stats, f, indent=2)

    print(f"\nSaved: {csv_path}")
    print(f"Saved: {json_path}")


if __name__ == "__main__":
    main()
