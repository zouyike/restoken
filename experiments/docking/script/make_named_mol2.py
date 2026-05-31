#!/usr/bin/env python3
"""Produce a mol2 for an NCAA residue whose heavy-atom names MATCH a reference
PDB residue (e.g. the crystal). This lets molfile_to_params_polymer_v4
--use_original_atom_name emit a .params whose ATOM names line up with the input
PDB, so the residue can be read straight out of an NCAA-peptide complex.

Strategy: the capped SMILES (ACE-X-NME) supplies correct topology, caps and
bond orders; we substructure-match the crystal residue's heavy-atom graph onto
the capped molecule and copy the crystal atom names across. Cap atoms and all H
get fresh, collision-free names (write_M_session renames the caps to ACE/NME
later anyway).

Usage:
  make_named_mol2.py --smiles "<ACE-X-NME>" --res-pdb <residue.pdb> \
                     --out <named.mol2> [--obabel <obabel_bin>]
The residue.pdb must contain ONLY the target residue's ATOM/HETATM records.
"""
import argparse, subprocess, sys, os
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdDetermineBonds

def parse_residue_pdb(path):
    atoms = []  # (name, element, x, y, z)
    for ln in open(path):
        if ln.startswith(("ATOM", "HETATM")):
            name = ln[12:16].strip()
            elem = ln[76:78].strip() or name[0]
            xyz = (float(ln[30:38]), float(ln[38:46]), float(ln[46:54]))
            atoms.append((name, elem.capitalize(), xyz))
    return atoms

def xyz_block(atoms):
    s = [str(len(atoms)), "res"]
    for name, elem, (x, y, z) in atoms:
        s.append(f"{elem:2s} {x:10.4f} {y:10.4f} {z:10.4f}")
    return "\n".join(s) + "\n"

def single_bond_copy(mol):
    """RWMol copy with every bond order set to SINGLE (match on connectivity only)."""
    rw = Chem.RWMol(mol)
    for b in rw.GetBonds():
        b.SetBondType(Chem.BondType.SINGLE)
        b.SetIsAromatic(False)
    for a in rw.GetAtoms():
        a.SetIsAromatic(False)
        a.SetNoImplicit(True)
        a.SetNumExplicitHs(0)
        a.SetFormalCharge(0)
    return rw.GetMol()

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--smiles", required=True)
    ap.add_argument("--res-pdb", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--obabel", default="/public/home/genesis/miniconda3/envs/bio_env/bin/obabel")
    a = ap.parse_args()

    cryst = parse_residue_pdb(a.res_pdb)
    cryst_names = [c[0] for c in cryst]

    # crystal heavy-atom graph (connectivity only)
    q = Chem.MolFromXYZBlock(xyz_block(cryst))
    rdDetermineBonds.DetermineConnectivity(q)
    q_match = single_bond_copy(q)

    # capped molecule: correct topology + caps + bond orders
    m = Chem.AddHs(Chem.MolFromSmiles(a.smiles))
    if AllChem.EmbedMolecule(m, randomSeed=0xC0FFEE) != 0:
        AllChem.EmbedMolecule(m, randomSeed=1, useRandomCoords=True)
    AllChem.MMFFOptimizeMolecule(m, maxIters=2000)

    heavy = Chem.RemoveHs(m)
    for i, at in enumerate(heavy.GetAtoms()):
        at.SetIntProp("fullidx", i)  # RemoveHs keeps heavy order == first N atoms of m
    heavy_match = single_bond_copy(heavy)

    # Enumerate all subgraph matches and pick the one that best preserves atom
    # degree (so a terminal crystal atom like an N-methyl maps to a terminal
    # capped atom, not e.g. the ACE carbonyl carbon). Tie-break on element of
    # each matched atom's neighbours.
    matches = heavy_match.GetSubstructMatches(q_match, uniquify=False, maxMatches=100000)
    if not matches:
        sys.exit("ERROR: crystal residue graph not found in capped molecule.")
    qdeg = [a.GetDegree() for a in q.GetAtoms()]
    tdeg = [a.GetDegree() for a in heavy.GetAtoms()]
    qelem = [a.GetAtomicNum() for a in q.GetAtoms()]
    telem = [a.GetAtomicNum() for a in heavy.GetAtoms()]
    def cost(mt):
        c = 0
        for qi, ti in enumerate(mt):
            if qelem[qi] != telem[ti]:
                c += 1000                       # element mismatch is unacceptable
            c += abs(qdeg[qi] - tdeg[ti])       # prefer degree agreement
        return c
    match = min(matches, key=cost)
    if len(match) != len(cryst):
        sys.exit(f"ERROR: substructure match failed ({len(match)}/{len(cryst)} atoms).")

    # name array over full (H-bearing) molecule
    names = [None] * m.GetNumAtoms()
    for ci, hi in enumerate(match):
        names[hi] = cryst_names[ci]            # heavy idx == full idx (heavy atoms first)
    capc = 1
    for i, at in enumerate(m.GetAtoms()):
        if names[i] is not None:
            continue
        if at.GetAtomicNum() == 1:
            names[i] = None                     # fill later, after heavy
        else:
            names[i] = f"{at.GetSymbol()}P{capc}"; capc += 1   # cap heavy atom
    hc = 1
    for i, at in enumerate(m.GetAtoms()):
        if names[i] is None:
            names[i] = f"H{hc}"; hc += 1

    # ensure uniqueness
    seen = {}
    for i, nm in enumerate(names):
        if nm in seen:
            k = 1
            while f"{nm}{k}" in seen: k += 1
            names[i] = f"{nm}{k}"
        seen[names[i]] = 1

    # sdf -> obabel mol2 (SYBYL types + charges, element-only names, order preserved)
    sdf = a.out + ".tmp.sdf"
    Chem.MolToMolFile(m, sdf)
    subprocess.run([a.obabel, sdf, "-O", a.out], check=True,
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    # rewrite mol2 ATOM-name field (col 2) by atom order
    lines = open(a.out).read().splitlines()
    out, mode, ai = [], None, 0
    for ln in lines:
        if ln.startswith("@<TRIPOS>ATOM"): mode = "atom"; out.append(ln); continue
        if ln.startswith("@<TRIPOS>"): mode = None; out.append(ln); continue
        if mode == "atom" and ln.strip():
            f = ln.split()
            # fields: id name x y z type subst_id subst_name charge
            f[1] = names[ai]; ai += 1
            out.append("%7s %-8s %9s %9s %9s %-6s %3s %-8s %9s" %
                       (f[0], f[1], f[2], f[3], f[4], f[5],
                        f[6] if len(f) > 6 else "1",
                        f[7] if len(f) > 7 else "RES",
                        f[8] if len(f) > 8 else "0.0"))
        else:
            out.append(ln)
    open(a.out, "w").write("\n".join(out) + "\n")
    os.remove(sdf)
    print(f"wrote {a.out}: {ai} atoms, residue heavy names = {cryst_names}")

if __name__ == "__main__":
    main()
