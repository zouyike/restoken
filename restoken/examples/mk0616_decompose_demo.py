"""
MK-0616 (enlicitide) decompose -> ResToken -> reconstruct round trip.

Question answered: "Can ResToken tokenize a real macrocyclic drug, and what
does the tokenization capture vs. lose?"

Pipeline (all local, no API):
    1. Parse published MK-0616 SMILES.
    2. Decompose it with PeptideDecomposer: locate backbone amides
       (alpha/beta/gamma aware), break side-chain crosslinks, re-cap each
       residue, and match it to the nearest of the 400 library blocks by
       chirality-aware Morgan-fingerprint Tanimoto.
    3. Reconstruct a head-to-tail macrocycle from the recovered block IDs.
    4. Compare descriptor envelopes (original vs reconstructed) to quantify
       exactly what head-to-tail ResToken cannot encode.

Honest scope:
    MK-0616 is *bicyclic*: a head-to-tail macrocycle PLUS a large aromatic
    crosslink bridge carrying its quaternary-ammonium charge. ResToken's
    vocabulary is head-to-tail only, so the bridge mass, the +1 charge, and
    the exotic linkers fall outside the block set. They surface here as
    low-Tanimoto "bridge-bearing" residues, detected crosslinks, filtered
    linker debris, and a quantified MW / charge / ring deficit in the
    reconstruction. This is a measurement of the framework's reach, not a
    claim that the reconstruction reproduces MK-0616.

Outputs:
    - console: per-residue match table + crosslink/debris/envelope accounting
    - CSV: per-residue records + envelope comparison
    - CDXML (optional, --cdxml): ACS-1996 ChemDraw of MK-0616 + reconstruction
      + the recovered residue blocks, via the shared cdxml-generation skill
    - PNG (optional, --png): RDKit quick-look grid

Usage:
    python mk0616_decompose_demo.py
    python mk0616_decompose_demo.py --cdxml out.cdxml --png out.png --csv out.csv
"""

import argparse
import csv
import subprocess
import sys
from pathlib import Path

from rdkit import RDLogger, Chem
from rdkit.Chem import Crippen, Descriptors, rdMolDescriptors, rdDepictor, Draw

RDLogger.DisableLog("rdApp.*")

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent))
from restoken.src.library import BlockLibrary
from restoken.src.cyclic_assembler import CyclicPeptideAssembler
from restoken.src.decomposer import PeptideDecomposer

# Published MK-0616 (enlicitide) SMILES — same string used in the manuscript.
ENLICITIDE_SMILES = (
    "C[C@H]1C(=O)N[C@@H](C(=O)N[C@H]2CC3=CC(=CC=C3)CNC(=O)CO[C@H]4"
    "CCN5[C@@H]4C(=O)N[C@H](C(=O)N[C@H](C(=O)N6CCC[C@]6(C(=O)NCCC7"
    "=CC=C(CN(CCCCCCN8C=C(C[C@@H](C5=O)NC2=O)C9=C8C=CC(=C9)F)C(=O)"
    "CCC(=O)N1)C=C7)C)CC1=CC=C(C=C1)OC)[C@@H](C)O)CNC(=O)CCCCC[N+](C)(C)C"
)

CDXML_HELPER = (
    Path.home() / "agent_exchange/shared_skills/helpers/smiles_to_cdxml_grid.py"
)

# Unified manuscript palette (navy / brown / teal / grayblue / darkred).
NAVY = "#002060"


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


def main():
    ap = argparse.ArgumentParser(description="MK-0616 ResToken decompose/reconstruct")
    here = Path(__file__).parent
    ap.add_argument("--bridge_tol", type=float, default=0.60,
                    help="Tanimoto below which a residue is flagged bridge-bearing")
    ap.add_argument("--csv", default=str(here / "mk0616_decompose_results.csv"))
    ap.add_argument("--cdxml", nargs="?", const=str(here / "mk0616_decompose.cdxml"),
                    default=None, help="Write ACS-1996 ChemDraw CDXML (borrowed skill)")
    ap.add_argument("--png", nargs="?", const=str(here / "mk0616_decompose.png"),
                    default=None, help="Write RDKit quick-look grid PNG")
    ap.add_argument("--backend", default=None)
    args = ap.parse_args()

    lib = BlockLibrary(backend_path=args.backend)
    dec = PeptideDecomposer(library=lib)
    asm = CyclicPeptideAssembler(backend_path=args.backend)

    ref_mol = Chem.MolFromSmiles(ENLICITIDE_SMILES)
    if ref_mol is None:
        print("ERROR: could not parse MK-0616 SMILES", file=sys.stderr)
        sys.exit(1)
    target = measure(ref_mol)

    # 1-2. Decompose.
    res = dec.tokenize(ENLICITIDE_SMILES)
    print("=== MK-0616 -> ResToken decomposition ===")
    print(f"  backbone residues recovered : {res['n_residues']}")
    print(f"  side-chain crosslinks broken : {res['n_crosslinks']}  "
          f"(head-to-tail ResToken cannot encode these)")
    print(f"  linker debris stubs filtered : {res['n_debris']}")
    print(f"  mean residue Tanimoto        : {res['mean_tanimoto']:.3f}")
    print(f"  recovered sequence           : {res['sequence']}")
    print("\n  per-residue match:")
    rows = []
    for i, r in enumerate(res["residues"], 1):
        b = lib.get(r["block"])
        cl = f"{b.aa_class}/{b.mc_type}/{b.chirality}/{b.mc_nmod}"
        flag = "bridge-bearing" if r["tanimoto"] < args.bridge_tol else "clean"
        print(f"   {i:2d}. {r['block']:5s}  T={r['tanimoto']:.3f}  {cl:22s} {flag}")
        rows.append({"pos": i, "block": r["block"], "tanimoto": r["tanimoto"],
                     "class": cl, "match": flag, "capped_smiles": r["capped"]})

    # 3. Reconstruct head-to-tail.
    seq = res["sequence"]
    recon = asm.assemble(seq)
    print("\n=== reconstruction (head-to-tail) ===")
    if not recon:
        print("  reconstruction failed", file=sys.stderr)
        recon_env = None
    else:
        recon_mol = Chem.MolFromSmiles(recon)
        recon_env = measure(recon_mol)
        print(f"  reconstructed SMILES: {recon}")
        print("\n  envelope comparison (what tokenization loses):")
        print(f"  {'prop':8s}{'MK-0616':>12s}{'ResToken':>12s}{'delta':>12s}")
        for k in ("mw", "charge", "rings", "heavy", "alogp", "tpsa", "hbd", "hba"):
            d = recon_env[k] - target[k]
            print(f"  {k:8s}{target[k]:>12.1f}{recon_env[k]:>12.1f}{d:>12.1f}")
        print(f"\n  => the {target['mw']-recon_env['mw']:.0f} Da / "
              f"{target['charge']-recon_env['charge']:+.0f} charge / "
              f"{target['rings']-recon_env['rings']:.0f} ring deficit is the "
              f"bicyclic crosslink bridge ResToken cannot represent.")

    # 4. Can the crosslink primitive (assemble_advanced) close that deficit?
    # It re-forms a bridge as a single bond between two *side-chain* atoms of
    # otherwise head-to-tail residues. That only works when each detected
    # crosslink is atom-disjoint from the backbone. Report the verdict honestly.
    print("\n=== crosslink-primitive re-encodability ===")
    rep = dec.crosslink_reencodability(ENLICITIDE_SMILES)
    print(f"  detected side-chain crosslinks : {rep['n_crosslinks']}")
    print(f"  disjoint (re-encodable)        : {rep['n_disjoint']}")
    print(f"  fused with backbone (cannot)   : {rep['n_fused']}")
    for d in rep["detail"]:
        kind = "FUSED (shares backbone atom %s)" % d["shared_atoms"] if d["fused"] \
            else "disjoint side-chain staple"
        print(f"   bond {d['bond'][0]}-{d['bond'][1]}: {kind}")
    if rep["reencodable"]:
        print("  => bridges are disjoint staples; assemble_advanced can re-close them.")
    else:
        print("  => MK-0616's bridges emanate from backbone-amide carbons, so the\n"
              "     macrocycle and bridge share atoms (a fused polycycle). This is\n"
              "     outside the head-to-tail + disjoint-staple model even with the\n"
              "     crosslink primitive: it is a measured representational limit,\n"
              "     not an assembler bug. (The primitive itself is proven on a clean\n"
              "     designed staple in tests/test_crosslink_assembler.py.)")

    # CSV.
    with open(args.csv, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["pos", "block", "tanimoto", "class", "match", "capped_smiles"])
        for r in rows:
            w.writerow([r["pos"], r["block"], r["tanimoto"], r["class"],
                        r["match"], r["capped_smiles"]])
        w.writerow([])
        w.writerow(["envelope", "property", "mk0616", "restoken", "delta"])
        if recon_env:
            for k in ("mw", "charge", "rings", "heavy", "alogp", "tpsa", "hbd", "hba"):
                w.writerow(["", k, f"{target[k]:.1f}", f"{recon_env[k]:.1f}",
                            f"{recon_env[k]-target[k]:.1f}"])
    print(f"\nSaved CSV: {args.csv}")

    # PNG quick-look (RDKit).
    if args.png and recon:
        mols, legs = [ref_mol, Chem.MolFromSmiles(recon)], \
            [f"MK-0616 (MW {target['mw']:.0f}, q {target['charge']:+d})",
             f"ResToken reconstruction (MW {recon_env['mw']:.0f}, q {recon_env['charge']:+d})"]
        for r in res["residues"]:
            m = Chem.MolFromSmiles(r["capped"])
            if m:
                rdDepictor.Compute2DCoords(m)
                mols.append(m)
                legs.append(f"{r['block']}  T={r['tanimoto']:.2f}")
        for m in mols[:2]:
            rdDepictor.Compute2DCoords(m)
        img = Draw.MolsToGridImage(mols, molsPerRow=3, subImgSize=(700, 560),
                                   legends=legs, useSVG=False)
        img.save(args.png)
        print(f"Saved PNG: {args.png}")

    # CDXML (borrowed cdxml-generation skill). One grid: the two whole
    # molecules followed by every recovered residue block, labelled with its
    # block id, Tanimoto, and clean/bridge-bearing call. Uses the helper's
    # JSON --input mode so labels never collide with comma splitting.
    if args.cdxml and recon:
        import json
        import tempfile
        entries = [
            {"smiles": ENLICITIDE_SMILES,
             "label": f"MK-0616  (MW {target['mw']:.0f}  q{target['charge']:+d})"},
            {"smiles": recon,
             "label": f"ResToken recon  (MW {recon_env['mw']:.0f}  q{recon_env['charge']:+d})"},
        ]
        for i, r in enumerate(res["residues"], 1):
            flag = "bridge" if r["tanimoto"] < args.bridge_tol else "clean"
            entries.append({"smiles": r["capped"],
                            "label": f"{i}. {r['block']}  T={r['tanimoto']:.2f}  ({flag})"})
        tmpdir = Path("/public/home/genesis/.local/tmp")
        tmpdir.mkdir(parents=True, exist_ok=True)
        with tempfile.NamedTemporaryFile("w", suffix=".json", dir=tmpdir,
                                         delete=False) as jf:
            json.dump(entries, jf)
            json_path = jf.name
        cmd = [sys.executable, str(CDXML_HELPER),
               "--input", json_path, "--cols", "3", "--output", args.cdxml]
        if args.backend:
            cmd += ["--backend", args.backend]
        try:
            subprocess.run(cmd, check=True)
            print(f"Saved CDXML: {args.cdxml}  ({len(entries)} structures)")
        except (subprocess.CalledProcessError, FileNotFoundError) as e:
            print(f"CDXML generation skipped ({e})", file=sys.stderr)
        finally:
            Path(json_path).unlink(missing_ok=True)


if __name__ == "__main__":
    main()
