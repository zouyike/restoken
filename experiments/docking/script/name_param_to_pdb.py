#!/usr/bin/env python3
"""Rename the atoms of an already-valid Rosetta polymer .params so its heavy-atom
names match a reference PDB residue (e.g. the crystal residue as it appears in an
NCAA-peptide complex). Topology/ICOOR are untouched -- only atom NAMES change --
so the param stays loadable while becoming readable directly from the complex PDB.

How: build the param's heavy-atom graph (from BOND records) and the crystal
residue's heavy-atom graph (connectivity from 3D), find the graph isomorphism
(degree- and element-aware), and apply the resulting name map across every
record that references an atom name.

Usage: name_param_to_pdb.py <in.params> <residue.pdb> <out.params>
The residue.pdb must contain ONLY the target residue's heavy ATOM/HETATM records.
"""
import sys
from rdkit import Chem
from rdkit.Chem import rdDetermineBonds

ELEM_FROM_PREFIX = {  # Rosetta atom-type -> element (fallback to name[0])
    "Nbb": "N", "CAbb": "C", "CObb": "C", "OCbb": "O",
}
PT = Chem.GetPeriodicTable()

def elem_of(name, rtype):
    if rtype in ELEM_FROM_PREFIX:
        return ELEM_FROM_PREFIX[rtype]
    n = name.lstrip("0123456789")
    if n[:2].capitalize() in ("Cl", "Br"):
        return n[:2].capitalize()
    return n[0].upper()

def parse_params(path):
    atoms = {}        # name -> (rtype, element, is_virt, is_H)
    order = []
    bonds = []
    for ln in open(path):
        if ln.startswith("ATOM "):
            f = ln.split()
            name, rtype = f[1], f[2]
            isH = name.lstrip("0123456789").startswith("H")
            isV = (rtype == "VIRT")
            atoms[name] = (rtype, elem_of(name, rtype), isV, isH)
            order.append(name)
        elif ln.startswith("BOND") and not ln.startswith("BOND_TYPE"):
            f = ln.split()
            bonds.append((f[1], f[2]))
        elif ln.startswith("BOND_TYPE"):
            f = ln.split()
            bonds.append((f[1], f[2]))
    return atoms, order, bonds

def heavy_mol(node_names, node_elems, edges):
    """RWMol with one atom per heavy node (element only, all single bonds)."""
    rw = Chem.RWMol()
    idx = {}
    for nm, el in zip(node_names, node_elems):
        a = Chem.Atom(PT.GetAtomicNumber(el))
        a.SetNoImplicit(True)
        idx[nm] = rw.AddAtom(a)
    for a, b in edges:
        if a in idx and b in idx and not rw.GetBondBetweenAtoms(idx[a], idx[b]):
            rw.AddBond(idx[a], idx[b], Chem.BondType.SINGLE)
    m = rw.GetMol()
    order = list(idx.keys())
    return m, order

def param_heavy_mol(atoms, order, bonds):
    names = [n for n in order if not (atoms[n][2] or atoms[n][3])]  # drop VIRT/H
    elems = [atoms[n][1] for n in names]
    return heavy_mol(names, elems, bonds)

def crystal_resname(path):
    for ln in open(path):
        if ln.startswith(("ATOM", "HETATM")):
            return ln[17:20].strip()
    return None

def crystal_heavy_mol(path):
    names, elems, xyz = [], [], []
    for ln in open(path):
        if ln.startswith(("ATOM", "HETATM")):
            nm = ln[12:16].strip()
            el = (ln[76:78].strip() or nm[0]).capitalize()
            if el == "H":
                continue
            names.append(nm); elems.append(el)
            xyz.append((float(ln[30:38]), float(ln[38:46]), float(ln[46:54])))
    blk = [str(len(names)), "res"]
    for el, (x, y, z) in zip(elems, xyz):
        blk.append(f"{el:2s} {x:10.4f} {y:10.4f} {z:10.4f}")
    cm = Chem.MolFromXYZBlock("\n".join(blk) + "\n")
    rdDetermineBonds.DetermineConnectivity(cm)
    edges = [(names[b.GetBeginAtomIdx()], names[b.GetEndAtomIdx()]) for b in cm.GetBonds()]
    return heavy_mol(names, elems, edges)

def main():
    inp, respdb, outp = sys.argv[1], sys.argv[2], sys.argv[3]
    atoms, order, bonds = parse_params(inp)
    resname = crystal_resname(respdb)  # name3 Rosetta matches PDB HETATM against
    pmol, pnames = param_heavy_mol(atoms, order, bonds)
    cmol, cnames = crystal_heavy_mol(respdb)
    if pmol.GetNumAtoms() != cmol.GetNumAtoms():
        sys.exit(f"ERROR: heavy-atom count differs param={pmol.GetNumAtoms()} "
                 f"crystal={cmol.GetNumAtoms()}")
    # full-graph isomorphism: same size -> substruct match of crystal onto param
    matches = pmol.GetSubstructMatches(cmol, uniquify=False, maxMatches=200000)
    if not matches:
        sys.exit("ERROR: param and crystal residue graphs are not isomorphic.")
    # mt maps crystal-atom-index -> param-atom-index; build param_name->crystal_name.
    # Prefer the iso that already agrees most (stabilises symmetric atoms).
    best, best_score = None, -1
    for mt in matches:
        mp = {pnames[ti]: cnames[ci] for ci, ti in enumerate(mt)}
        score = sum(1 for k, v in mp.items() if k == v)
        if score > best_score:
            best, best_score = mp, score
    mapping = best
    # crystal names must be unique and must not collide with kept H/VIRT names
    kept = {n for n in order if n not in mapping}
    for k, v in mapping.items():
        if v in kept:
            sys.exit(f"ERROR: crystal name {v} collides with non-mapped atom.")

    def rename_token(tok):
        return mapping.get(tok, tok)

    out = []
    for ln in open(inp):
        s = ln.rstrip("\n")
        key = s.split()[0] if s.split() else ""
        if key == "IO_STRING" and resname:
            f = s.split()
            name1 = f[2] if len(f) > 2 else "X"
            out.append("IO_STRING %s %s" % (resname, name1))
        elif key in ("ATOM",):
            f = s.split()
            f[1] = rename_token(f[1])
            out.append("ATOM %-4s %-4s %-4s %s" % (f[1], f[2], f[3], " ".join(f[4:])))
        elif key in ("BOND", "BOND_TYPE", "CUT_BOND", "CHI", "NBR_ATOM",
                     "FIRST_SIDECHAIN_ATOM", "MAINCHAIN_ATOMS", "ACT_COORD_ATOMS",
                     "METAL_BINDING_ATOMS", "VIRTUAL_SHADOW", "ICOOR_INTERNAL",
                     "LOWER_CONNECT", "UPPER_CONNECT", "CONNECT"):
            f = s.split()
            f = [f[0]] + [rename_token(t) for t in f[1:]]
            out.append(" ".join(f))
        else:
            out.append(s)
    open(outp, "w").write("\n".join(out) + "\n")
    print(f"renamed {len(mapping)} heavy atoms -> {sorted(mapping.values())}")
    print(f"identity-preserved: {best_score}/{len(mapping)}")

if __name__ == "__main__":
    main()
