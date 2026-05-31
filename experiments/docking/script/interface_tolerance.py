"""
Params-free interface-tolerance readout for a grafted LUNA18 analogue in the
crystal bound pose. Chemistry-agnostic (works for any NCAA side chain), needs no
Rosetta params. Quantifies, for the edited residue:

    buried_SASA  = SASA(side chain in isolated peptide) - SASA(in complex)
                   -> how much of the new side chain KRAS actually buries
    kras_contacts= heavy-atom contacts < 4.0 A to KRAS
    clashes      = heavy-atom pairs < 2.2 A to KRAS (steric tolerance flag)

A swap that buries MORE hydrophobic surface without clashing is "tolerated and
possibly favorable"; one that clashes or loses burial is flagged. This is a
structure/tolerance readout, NOT an affinity prediction.

Usage:
    python interface_tolerance.py --complex <complex.pdb> --resnum 9 --label adamantyl
"""
import argparse
import numpy as np
from Bio.PDB import PDBParser
from Bio.PDB.SASA import ShrakeRupley

BB = {"N", "CA", "C", "O"}


def load(path):
    return PDBParser(QUIET=True).get_structure("x", path)[0]


def sasa(model):
    ShrakeRupley(probe_radius=1.4, n_points=200).compute(model, level="A")


def atoms_of(model, chain, resnum):
    out = []
    for r in model[chain]:
        if r.id[1] == resnum:
            out.extend(r.get_atoms())
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--complex", required=True)
    ap.add_argument("--resnum", type=int, required=True)
    ap.add_argument("--label", required=True)
    ap.add_argument("--pep_chain", default="A")
    ap.add_argument("--tgt_chain", default="B")
    args = ap.parse_args()

    # complex (peptide A + KRAS B)
    cx = load(args.complex)
    sasa(cx)
    res_atoms_cx = atoms_of(cx, args.pep_chain, args.resnum)
    side_cx = {a.name: a.sasa for a in res_atoms_cx if a.name not in BB}
    pep_sasa_cx = sum(a.sasa for ch in cx if ch.id == args.pep_chain for a in ch.get_atoms())

    # isolated peptide: reload and drop target chain
    iso = load(args.complex)
    iso.detach_child(args.tgt_chain)
    sasa(iso)
    res_atoms_iso = atoms_of(iso, args.pep_chain, args.resnum)
    side_iso = {a.name: a.sasa for a in res_atoms_iso if a.name not in BB}
    pep_sasa_iso = sum(a.sasa for ch in iso for a in ch.get_atoms())

    buried_side = sum(side_iso[k] - side_cx.get(k, 0) for k in side_iso)
    buried_total = pep_sasa_iso - pep_sasa_cx
    side_exposed_iso = sum(side_iso.values())

    # contacts / clashes to target
    tgt = np.array([a.coord for ch in cx if ch.id == args.tgt_chain for a in ch.get_atoms()])
    side_coords = [(a.name, a.coord) for a in res_atoms_cx if a.name not in BB]
    contacts = clashes = 0
    closest = 99.0
    for nm, c in side_coords:
        d = np.linalg.norm(tgt - c, axis=1)
        closest = min(closest, d.min())
        if (d < 4.0).any():
            contacts += 1
        if (d < 2.2).any():
            clashes += 1

    print(f"=== {args.label}  (res {args.resnum}) ===")
    print(f"  side-chain heavy atoms        : {len(side_iso)}")
    print(f"  side-chain SASA (isolated)    : {side_exposed_iso:7.1f} A^2")
    print(f"  side-chain buried by KRAS     : {buried_side:7.1f} A^2  "
          f"({100*buried_side/max(side_exposed_iso,1e-6):4.0f}% of side chain)")
    print(f"  whole-peptide interface burial: {buried_total:7.1f} A^2")
    print(f"  side-chain atoms contacting KRAS(<4A): {contacts}/{len(side_coords)}")
    print(f"  side-chain clashes (<2.2A)    : {clashes}   closest approach {closest:.2f} A")


if __name__ == "__main__":
    main()
