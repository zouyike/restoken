"""
LUNA18 semantic-analogue design demo with ResToken.

Question answered: "Can we design analogues of LUNA18 by editing its ResToken
sequence the way you edit words in a sentence?"

What this IS:
    A reproducible, local (no-API) demonstration that, because ResToken encodes
    LUNA18 as a discrete token sequence

        A18-a31-a02-A10-i01-A170-a01-I03-A13-a31-a31

    each residue is an editable *token*. Swapping one token for a
    *semantically equivalent* block (same backbone type, chirality, charge, and
    side-chain role, but a different side chain) yields a synthesizable,
    structurally close analogue. Each analogue here is a SINGLE-token edit on
    LUNA18, and every edit preserves LUNA18's monocyclic topology and its
    permeability-defining backbone N-methylation pattern.

What this is NOT (honest scope):
    - NOT validated KRAS / target binders. ResToken matches chemistry, not
      affinity. Binding claims would route through the AF3-surrogate pipeline.
    - The baseline is LUNA18's *ResToken reconstruction*, NOT the literal
      drug. The reconstruction carries one known side-chain coverage gap:
      residue 6 (A170) is LUNA18's penta-substituted FLUORO-aryl, but the
      400-block library has no fluorinated aromatic, so it matches a
      penta-HYDROXY-aryl block (+5 spurious OH -> reconstruction HBD = 8 vs the
      drug's true 3). Every analogue inherits the same gap identically, so the
      reported deltas isolate the effect of each semantic edit rather than the
      coverage gap. See git e45090a / test_nmod_decompose.py.

The five semantic edits (each a one-token swap, chosen for medchem meaning):
    1. pos4  A10 -> A84   3,4-xylyl-methyl  -> phenyl        (des-dimethyl; trim lipophilicity)
    2. pos4  A10 -> A86   3,4-xylyl-methyl  -> 2-naphthyl    (aromatic ring-grow; deeper hydrophobic reach)
    3. pos9  A13 -> A19   1-aminocyclopentyl -> adamantyl    (rigidified lipophilic anchor)
    4. pos1  A18 -> A48   sec-butyl (Ile-like) -> 1-Me-cyclopropyl (metabolically soft isostere)
    5. pos3  a02 -> a01   azetidine-2-COOH  -> proline        (backbone ring homologation; conformational probe)

Usage:
    python luna18_semantic_analogue_demo.py
    python luna18_semantic_analogue_demo.py --output hits.png --csv_out hits.csv
"""

import argparse
import copy
import csv
import sys
from pathlib import Path

from rdkit import RDLogger, Chem
from rdkit.Chem import Crippen, Descriptors, rdMolDescriptors, rdDepictor, Draw

RDLogger.DisableLog("rdApp.*")

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent))
from restoken.src.decomposer import PeptideDecomposer
from restoken.src.library import BlockLibrary
from restoken.src.cyclic_assembler import CyclicPeptideAssembler

rdDepictor.SetPreferCoordGen(True)

# Published LUNA18 SMILES (same string used in test_nmod_decompose.py).
LUNA18 = (
    "CC[C@H](C)[C@H]1C(=O)N([C@H](C(=O)N2CC[C@H]2C(=O)N([C@H](C(=O)N(CC(=O)N"
    "[C@H](C(=O)N3CCC[C@H]3C(=O)NC4(CCCC4)C(=O)N([C@H](C(=O)N([C@@H](CC(=O)N"
    "([C@H](C(=O)N1)CC(C)C)C)C(=O)N(C)C)C)C5CCCC5)C)CCC6=CC(=C(C(=C6)F)C(F)"
    "(F)F)F)C)CC7=CC=C(C=C7)C)CC)C)C"
)

# (residue index 0-based, new block id, short rationale tag)
SEMANTIC_EDITS = [
    (3, "A84", "A10->A84 phenyl (des-Me2)"),
    (3, "A86", "A10->A86 2-naphthyl (ring-grow)"),
    (8, "A19", "A13->A19 adamantyl (rigid)"),
    (0, "A48", "A18->A48 Me-cyclopropyl (soft)"),
    (2, "a01", "a02->a01 proline (ring homolog)"),
]


def measure(mol):
    return {
        "mw": Descriptors.MolWt(mol),
        "alogp": Crippen.MolLogP(mol),
        "tpsa": rdMolDescriptors.CalcTPSA(mol),
        "hbd": rdMolDescriptors.CalcNumHBD(mol),
        "hba": rdMolDescriptors.CalcNumHBA(mol),
        "rotb": rdMolDescriptors.CalcNumRotatableBonds(mol),
        "charge": Chem.GetFormalCharge(mol),
        "rings": rdMolDescriptors.CalcNumRings(mol),
        "heavy": mol.GetNumHeavyAtoms(),
    }


def edited_sequence(base_tokens, pos, new_block):
    toks = list(base_tokens)
    toks[pos] = f"[{new_block}]"  # mark the edited token
    return "-".join(toks)


def main():
    ap = argparse.ArgumentParser(description="LUNA18 semantic-analogue demo (ResToken)")
    ap.add_argument("--output", default=str(Path(__file__).parent / "luna18_semantic_analogue_results.png"))
    ap.add_argument("--csv_out", default=str(Path(__file__).parent / "luna18_semantic_analogue_results.csv"))
    ap.add_argument("--backend", default=None)
    args = ap.parse_args()

    dec = PeptideDecomposer(library=BlockLibrary())
    asm = CyclicPeptideAssembler(backend_path=args.backend) if args.backend else CyclicPeptideAssembler()

    # 1. Tokenize LUNA18 -> discrete sequence.
    res = dec.tokenize(LUNA18)
    base_tokens = res["sequence"].split("-")
    print(f"LUNA18 ResToken sequence: {res['sequence']}")
    print(f"  {res['n_residues']} residues, {res['n_crosslinks']} crosslinks, "
          f"{res['n_nmod']} N-modified")

    # Baseline = LUNA18's own ResToken reconstruction (shares the A170 gap).
    ref_smi = dec.reconstruct(res, assembler=asm, nmod_correct=True)
    ref_mol = Chem.MolFromSmiles(ref_smi)
    ref = measure(ref_mol)
    drug_mol = Chem.MolFromSmiles(LUNA18)
    drug = measure(drug_mol)
    print(f"\nLUNA18 drug      : MW={drug['mw']:.0f}  AlogP={drug['alogp']:.2f}  "
          f"HBD={drug['hbd']}  HBA={drug['hba']}  rings={drug['rings']}")
    print(f"LUNA18 reconstr. : MW={ref['mw']:.0f}  AlogP={ref['alogp']:.2f}  "
          f"HBD={ref['hbd']}  HBA={ref['hba']}  rings={ref['rings']}   <- analogue baseline")
    print("  (reconstr. HBD 8 = 3 true + 5 from the A170 fluoro->hydroxy coverage gap)")

    # 2-4. Apply each semantic edit, reconstruct (N-mod aware), measure.
    analogues = []
    for pos, new_block, tag in SEMANTIC_EDITS:
        r2 = copy.deepcopy(res)
        from_block = r2["residues"][pos]["block"]
        r2["residues"][pos]["block"] = new_block
        smi = dec.reconstruct(r2, assembler=asm, nmod_correct=True)
        mol = Chem.MolFromSmiles(smi) if smi else None
        if mol is None:
            print(f"  WARN: edit {tag} failed to assemble", file=sys.stderr)
            continue
        env = measure(mol)
        seq = edited_sequence(base_tokens, pos, new_block)
        analogues.append({
            "pos": pos + 1, "from": from_block, "to": new_block, "tag": tag,
            "seq": seq, "smi": smi, "mol": mol,
            "d_mw": env["mw"] - ref["mw"], "d_alogp": env["alogp"] - ref["alogp"],
            "d_hbd": env["hbd"] - ref["hbd"], "topology_ok": env["rings"] >= ref["rings"],
            **env,
        })

    print(f"\n{len(analogues)} semantic analogues (single-token edits), Δ vs LUNA18 reconstruction:")
    print(f"  {'edit':32s} {'MW':>6s} {'ΔMW':>6s} {'ΔAlogP':>7s} {'rings':>5s} {'topo':>5s}")
    for a in analogues:
        topo = "ok" if a["topology_ok"] else "BROKEN"
        print(f"  {a['tag']:32s} {a['mw']:6.0f} {a['d_mw']:+6.0f} {a['d_alogp']:+7.2f} "
              f"{a['rings']:5d} {topo:>5s}")

    # 5. Render: LUNA18 reconstruction first, then the 5 analogues.
    rdDepictor.Compute2DCoords(ref_mol)
    mols = [ref_mol]
    legends = [f"LUNA18 (ResToken reconstruction)  MW={ref['mw']:.0f} AlogP={ref['alogp']:.2f}"]
    for a in analogues:
        rdDepictor.Compute2DCoords(a["mol"])
        mols.append(a["mol"])
        legends.append(f"{a['tag']}  MW={a['mw']:.0f} ({a['d_mw']:+.0f})  AlogP={a['alogp']:.2f}")
    img = Draw.MolsToGridImage(mols, molsPerRow=3, subImgSize=(1000, 800),
                               legends=legends, useSVG=False)
    img.save(args.output)
    print(f"\nSaved figure: {args.output}")

    with open(args.csv_out, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["label", "token_sequence", "edit_pos", "from_block", "to_block",
                    "smiles", "mw", "alogp", "tpsa", "hbd", "hba", "rings", "heavy",
                    "d_mw", "d_alogp", "d_hbd"])
        w.writerow(["LUNA18_drug", res["sequence"], "", "", "", LUNA18,
                    f"{drug['mw']:.0f}", f"{drug['alogp']:.2f}", f"{drug['tpsa']:.0f}",
                    drug["hbd"], drug["hba"], drug["rings"], drug["heavy"], "", "", ""])
        w.writerow(["LUNA18_reconstruction", res["sequence"], "", "", "", ref_smi,
                    f"{ref['mw']:.0f}", f"{ref['alogp']:.2f}", f"{ref['tpsa']:.0f}",
                    ref["hbd"], ref["hba"], ref["rings"], ref["heavy"], "0", "0.00", "0"])
        for a in analogues:
            w.writerow([f"analogue_p{a['pos']}_{a['to']}", a["seq"], a["pos"],
                        a["from"], a["to"], a["smi"], f"{a['mw']:.0f}",
                        f"{a['alogp']:.2f}", f"{a['tpsa']:.0f}", a["hbd"], a["hba"],
                        a["rings"], a["heavy"], f"{a['d_mw']:+.0f}",
                        f"{a['d_alogp']:+.2f}", f"{a['d_hbd']:+d}"])
    print(f"Saved CSV: {args.csv_out}")


if __name__ == "__main__":
    main()
