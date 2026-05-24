"""
Demo 2: Guideline-based cell-penetrating cyclic peptide design.

Pipeline:
    1. Load guideline generator with CPP design rules
    2. Generate candidates via constrained sampling (CPP-like mode)
    3. Filter by hard constraints (charge, cationic count, ring size)
    4. Assemble top candidates into macrocyclic SMILES
    5. Compute physicochemical properties
    6. Render 2D structures

No LLM API key required — uses the built-in constrained sampler.

Usage:
    python cpp_design_demo.py
    python cpp_design_demo.py --ring-size 8 --n 100 --output cpp_hits.png
"""

import argparse
import csv
import sys
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import Crippen, Descriptors, rdDepictor, Draw

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent))
from restoken.src.guideline_generator import GuidelineGenerator
from restoken.src.cyclic_assembler import CyclicPeptideAssembler

rdDepictor.SetPreferCoordGen(True)

_GUIDELINE = Path(__file__).resolve().parent.parent / "guidelines" / "cell_penetrating_cyclic_peptide.md"


def main():
    parser = argparse.ArgumentParser(description="CPP cyclic peptide design demo")
    parser.add_argument("--ring-size", type=int, default=7,
                        help="Ring size (default: 7)")
    parser.add_argument("-n", type=int, default=100,
                        help="Number of candidates to generate (default: 100)")
    parser.add_argument("--top", type=int, default=6,
                        help="Number of top hits to render (default: 6)")
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--output", default="cpp_design_results.png")
    parser.add_argument("--csv-out", default=None, help="CSV output path")
    parser.add_argument("--guideline", default=None,
                        help="Custom guideline file (default: built-in CPP guideline)")
    parser.add_argument("--exclude-exotic", action="store_true",
                        help="Exclude blocks with exotic flags (high MW, poly-ring, etc.)")
    parser.add_argument("--max-exotic", type=int, default=None,
                        help="Max exotic blocks per sequence")
    args = parser.parse_args()

    guideline = args.guideline or str(_GUIDELINE)
    if Path(guideline).exists():
        gen = GuidelineGenerator(guideline)
    else:
        gen = GuidelineGenerator.from_params()

    print(f"Generating {args.n} CPP-like candidates (ring_size={args.ring_size})...")
    candidates = gen.generate(n=args.n, mode="cpp_like",
                              ring_size=args.ring_size, seed=args.seed,
                              exclude_exotic=args.exclude_exotic,
                              max_exotic_per_seq=args.max_exotic)
    filtered = gen.filter(candidates, mode="cpp_like")
    print(f"  {len(candidates)} generated, {len(filtered)} passed hard filters")

    asm = CyclicPeptideAssembler()
    results = []
    for c in filtered:
        try:
            smi = asm.assemble(c["sequence"])
            if smi is None:
                continue
            mol = Chem.MolFromSmiles(smi)
            if mol is None:
                continue
            alogp = Crippen.MolLogP(mol)
            mw = Descriptors.MolWt(mol)
            hbd = Descriptors.NumHDonors(mol)
            hba = Descriptors.NumHAcceptors(mol)
            tpsa = Descriptors.TPSA(mol)
            results.append({
                **c, "smiles": smi, "mol": mol,
                "alogp": alogp, "mw": mw, "hbd": hbd, "hba": hba, "tpsa": tpsa,
            })
        except Exception as e:
            print(f"  Assembly error for {c['sequence']}: {e}")

    if not results:
        print("No valid molecules produced.")
        sys.exit(1)

    print(f"\n{'Sequence':<45} {'Score':>5} {'Chg':>4} {'Cat':>4} {'Aro':>4} "
          f"{'MW':>6} {'AlogP':>6} {'TPSA':>5} {'HBD':>4}")
    print("-" * 110)
    for r in results[:20]:
        print(f"  {r['sequence']:<43} {r['score']:5.1f} {r['net_charge']:+4d} "
              f"{r['cationic_count']:4d} {r['aromatic_hydrophobe_count']:4d} "
              f"{r['mw']:6.0f} {r['alogp']:6.2f} {r['tpsa']:5.0f} {r['hbd']:4d}")

    best = results[:args.top]
    mols, legends = [], []
    for r in best:
        rdDepictor.Compute2DCoords(r["mol"])
        mols.append(r["mol"])
        legends.append(
            f"{r['sequence']}\n"
            f"score={r['score']:.1f}  charge={r['net_charge']:+d}  "
            f"MW={r['mw']:.0f}  AlogP={r['alogp']:.1f}"
        )

    img = Draw.MolsToGridImage(mols, molsPerRow=3, subImgSize=(1000, 800),
                                legends=legends, useSVG=False)
    img.save(args.output)
    print(f"\nRendered top {len(best)} structures: {args.output}")

    print(f"\nSummary ({len(results)} assembled):")
    print(f"  Net charge range: {min(r['net_charge'] for r in results):+d} to "
          f"{max(r['net_charge'] for r in results):+d}")
    print(f"  MW range: {min(r['mw'] for r in results):.0f} - "
          f"{max(r['mw'] for r in results):.0f}")
    print(f"  AlogP range: {min(r['alogp'] for r in results):.2f} - "
          f"{max(r['alogp'] for r in results):.2f}")
    print(f"  Avg score: {sum(r['score'] for r in results)/len(results):.1f}")

    if args.csv_out:
        fields = ["sequence", "smiles", "score", "net_charge", "cationic_count",
                  "arg_like_count", "aromatic_hydrophobe_count", "ring_size",
                  "mw", "alogp", "tpsa", "hbd", "hba", "priority"]
        with open(args.csv_out, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
            w.writeheader()
            for r in results:
                w.writerow(r)
        print(f"  Saved: {args.csv_out}")


if __name__ == "__main__":
    main()
