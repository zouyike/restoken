"""
Semantic-analogue design for any cyclic peptide via ResToken token editing.

Generalizes the LUNA18 worked example: decode a macrocyclic-peptide SMILES into
its ResToken sequence, then design analogues by swapping ONE token for a
semantically equivalent block (same backbone type, chirality, charge, and
side-chain role -- different side chain), the way you edit a word in a sentence.
Each edit reconstructs N-mod-aware, so the parent's backbone N-methylation
(its permeability driver) and macrocyclic topology are preserved.

Honest baseline: analogues are scored against the parent's own ResToken
*reconstruction*, NOT the literal drug. Any side-chain the 400-block library
cannot represent (e.g. a fluoroaromatic, or cyclosporine's MeBmt side chain)
appears as a fixed coverage gap in the reconstruction; because every analogue
inherits it identically, the reported deltas isolate each edit's effect.

These are structure/property variants, NOT validated target binders -- ResToken
matches chemistry, not affinity. Binding claims route through the AF3-surrogate
pipeline.

Usage:
    python semantic_analogue_demo.py --name LUNA18 --smiles "<SMILES>" \
        --edits "4:A84,4:A86,9:A19,1:A48,3:a01" \
        --out_prefix luna18_semantic_analogue
    # --edits are 1-based: "<position>:<new_block_id>", comma-separated.
    # Find candidate blocks per position with --list_neighbors.
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


def semantic_neighbors(lib, block, strict=True):
    """Blocks sharing backbone class/chirality/charge/type (+ polarity/bulk if
    strict), with a different side chain -- i.e. semantically interchangeable."""
    out = []
    for oid in lib.all_ids:
        o = lib.get(oid)
        if o.id == block.id or o.sc_smiles == block.sc_smiles:
            continue
        if (o.aa_class != block.aa_class or o.chirality != block.chirality
                or o.charge != block.charge or o.mc_type != block.mc_type):
            continue
        if strict and (o.polarity_bin != block.polarity_bin
                       or o.sc_bulk != block.sc_bulk):
            continue
        out.append(oid)
    return out


def main():
    ap = argparse.ArgumentParser(description="Semantic-analogue design (ResToken)")
    ap.add_argument("--name", required=True, help="Parent compound name")
    ap.add_argument("--smiles", required=True, help="Parent macrocyclic-peptide SMILES")
    ap.add_argument("--edits", default="",
                    help="Comma-separated 1-based 'pos:new_block' single-token edits")
    ap.add_argument("--out_prefix", required=True,
                    help="Output prefix (writes <prefix>_results.{png,csv})")
    ap.add_argument("--list_neighbors", action="store_true",
                    help="Just print semantic neighbors per position and exit")
    ap.add_argument("--backend", default=None)
    args = ap.parse_args()

    lib = BlockLibrary()
    dec = PeptideDecomposer(library=lib)
    asm = CyclicPeptideAssembler(backend_path=args.backend) if args.backend else CyclicPeptideAssembler()

    res = dec.tokenize(args.smiles)
    base_tokens = res["sequence"].split("-")
    print(f"{args.name} ResToken sequence: {res['sequence']}")
    print(f"  {res['n_residues']} residues, {res['n_crosslinks']} crosslinks, "
          f"{res['n_nmod']} N-modified, mean Tanimoto {res['mean_tanimoto']:.3f}")

    if args.list_neighbors:
        for i, bid in enumerate(base_tokens):
            b = lib.get(bid)
            ns = semantic_neighbors(lib, b, strict=True)
            ex = [(x, lib.get(x).sc_smiles) for x in ns[:5]]
            print(f"  pos{i+1:2d} {bid:5s} aa={b.aa_class} chir={b.chirality} "
                  f"bulk={b.sc_bulk:5s} pol={b.polarity_bin:4s} sc={b.sc_smiles:20s} "
                  f"| {len(ns)} strict: {ex}")
        return

    ref_smi = dec.reconstruct(res, assembler=asm, nmod_correct=True)
    ref = measure(Chem.MolFromSmiles(ref_smi))
    drug = measure(Chem.MolFromSmiles(args.smiles))
    print(f"\n{args.name} drug       : MW={drug['mw']:.0f} AlogP={drug['alogp']:.2f} "
          f"HBD={drug['hbd']} HBA={drug['hba']} rings={drug['rings']}")
    print(f"{args.name} reconstr.  : MW={ref['mw']:.0f} AlogP={ref['alogp']:.2f} "
          f"HBD={ref['hbd']} HBA={ref['hba']} rings={ref['rings']}   <- analogue baseline")

    analogues = []
    for tok in [e for e in args.edits.split(",") if e.strip()]:
        pos_s, new_block = tok.split(":")
        pos = int(pos_s) - 1
        r2 = copy.deepcopy(res)
        from_block = r2["residues"][pos]["block"]
        r2["residues"][pos]["block"] = new_block
        smi = dec.reconstruct(r2, assembler=asm, nmod_correct=True)
        mol = Chem.MolFromSmiles(smi) if smi else None
        if mol is None:
            print(f"  WARN: edit pos{pos+1} {from_block}->{new_block} failed", file=sys.stderr)
            continue
        env = measure(mol)
        seq = "-".join(f"[{new_block}]" if j == pos else t
                       for j, t in enumerate(base_tokens))
        analogues.append({
            "pos": pos + 1, "from": from_block, "to": new_block,
            "tag": f"pos{pos+1} {from_block}->{new_block}", "seq": seq,
            "smi": smi, "mol": mol, "d_mw": env["mw"] - ref["mw"],
            "d_alogp": env["alogp"] - ref["alogp"], "d_hbd": env["hbd"] - ref["hbd"],
            "topology_ok": env["rings"] >= ref["rings"], **env,
        })

    print(f"\n{len(analogues)} semantic analogues (Δ vs {args.name} reconstruction):")
    for a in analogues:
        topo = "ok" if a["topology_ok"] else "BROKEN"
        print(f"  {a['tag']:20s} MW={a['mw']:6.0f} ΔMW={a['d_mw']:+5.0f} "
              f"ΔAlogP={a['d_alogp']:+6.2f} rings={a['rings']} {topo}")

    rdDepictor.SetPreferCoordGen(True)
    ref_mol = Chem.MolFromSmiles(ref_smi)
    rdDepictor.Compute2DCoords(ref_mol)
    mols = [ref_mol]
    legends = [f"{args.name} (reconstruction)  MW={ref['mw']:.0f} AlogP={ref['alogp']:.2f}"]
    for a in analogues:
        rdDepictor.Compute2DCoords(a["mol"])
        mols.append(a["mol"])
        legends.append(f"{a['tag']}  MW={a['mw']:.0f} ({a['d_mw']:+.0f}) AlogP={a['alogp']:.2f}")
    png = f"{args.out_prefix}_results.png"
    Draw.MolsToGridImage(mols, molsPerRow=3, subImgSize=(1000, 800),
                         legends=legends, useSVG=False).save(png)
    print(f"\nSaved figure: {png}")

    csv_path = f"{args.out_prefix}_results.csv"
    with open(csv_path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["label", "token_sequence", "edit_pos", "from_block", "to_block",
                    "smiles", "mw", "alogp", "tpsa", "hbd", "hba", "rings", "heavy",
                    "d_mw", "d_alogp", "d_hbd"])
        w.writerow([f"{args.name}_drug", res["sequence"], "", "", "", args.smiles,
                    f"{drug['mw']:.0f}", f"{drug['alogp']:.2f}", f"{drug['tpsa']:.0f}",
                    drug["hbd"], drug["hba"], drug["rings"], drug["heavy"], "", "", ""])
        w.writerow([f"{args.name}_reconstruction", res["sequence"], "", "", "", ref_smi,
                    f"{ref['mw']:.0f}", f"{ref['alogp']:.2f}", f"{ref['tpsa']:.0f}",
                    ref["hbd"], ref["hba"], ref["rings"], ref["heavy"], "0", "0.00", "0"])
        for a in analogues:
            w.writerow([f"analogue_p{a['pos']}_{a['to']}", a["seq"], a["pos"],
                        a["from"], a["to"], a["smi"], f"{a['mw']:.0f}",
                        f"{a['alogp']:.2f}", f"{a['tpsa']:.0f}", a["hbd"], a["hba"],
                        a["rings"], a["heavy"], f"{a['d_mw']:+.0f}",
                        f"{a['d_alogp']:+.2f}", f"{a['d_hbd']:+d}"])
    print(f"Saved CSV: {csv_path}")


if __name__ == "__main__":
    main()
