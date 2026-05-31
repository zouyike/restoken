"""
Graft a new side chain / ring system onto ONE residue of crystal LUNA18, keeping
the backbone (N, CA, C, O) and the rest of the complex fixed in the bound pose.

Rationale: a single-token ResToken edit is a local substitution; the macrocycle
backbone and bound pose are conserved (that is the whole premise). So we transplant
only the new residue's side-chain atoms, rigid-aligned by the N-CA-C backbone frame
onto the crystal residue, and leave everything else at crystal coordinates.

The new residue is built from SMILES of the FREE amino acid (with -NH2 and -COOH).
Backbone atoms are detected topologically:
    Calpha = carbon bonded to both the amino N and the carboxyl C
    carboxyl C = carbon bonded to two oxygens
    amino N    = the nitrogen bonded to Calpha
Side-chain atoms = everything reachable from Calpha that is not backbone.

Usage:
    python graft_residue.py --resnum 9 --new_name A19 \
        --smiles "OC(=O)C1(N)C2CC3CC1CC(C2)C3" --tag adamantyl
"""
import argparse
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

PEP = "/scratch/genesis/NCAA_tokenization/experiments/docking/input/luna18_crystal.pdb"
KRAS = "/scratch/genesis/NCAA_tokenization/experiments/docking/input/receptor_kras.pdb"
OUT = "/scratch/genesis/NCAA_tokenization/experiments/docking/input"


def read_pdb_atoms(path):
    atoms = []
    for line in open(path):
        if line.startswith(("ATOM", "HETATM")):
            atoms.append({
                "record": line[:6], "name": line[12:16].strip(),
                "resname": line[17:20].strip(), "chain": line[21],
                "resnum": int(line[22:26]), "icode": line[26],
                "x": float(line[30:38]), "y": float(line[38:46]),
                "z": float(line[46:54]), "elem": line[76:78].strip(),
            })
    return atoms


def fmt_atom(serial, name, resname, chain, resnum, xyz, elem):
    nm = f" {name}" if len(name) < 4 else name
    return (f"HETATM{serial:>5} {nm:<4} {resname:>3} {chain}{resnum:>4}    "
            f"{xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}  1.00  0.00          {elem:>2}\n")


def build_residue(smiles):
    m = Chem.MolFromSmiles(smiles)
    m = Chem.AddHs(m)
    AllChem.EmbedMolecule(m, randomSeed=0xC0FFEE)
    AllChem.MMFFOptimizeMolecule(m)
    conf = m.GetConformer()

    # backbone detection
    carboxyl_C = None
    for a in m.GetAtoms():
        if a.GetSymbol() == "C":
            o = [n for n in a.GetNeighbors() if n.GetSymbol() == "O"]
            if len(o) == 2:
                carboxyl_C = a.GetIdx(); carboxyl_O = [x.GetIdx() for x in o]; break
    amino_N = None
    for a in m.GetAtoms():
        if a.GetSymbol() == "N":
            amino_N = a.GetIdx(); break
    # Calpha bonded to both amino_N and carboxyl_C
    ca = None
    for a in m.GetAtoms():
        nb = {n.GetIdx() for n in a.GetNeighbors()}
        if amino_N in nb and carboxyl_C in nb and a.GetSymbol() == "C":
            ca = a.GetIdx(); break
    assert ca is not None, "could not find Calpha"

    pos = {i: np.array(conf.GetAtomPosition(i)) for i in range(m.GetNumAtoms())}
    backbone = {amino_N, ca, carboxyl_C, *carboxyl_O}
    # side-chain heavy atoms = heavy atoms not in backbone
    side = [a.GetIdx() for a in m.GetAtoms()
            if a.GetIdx() not in backbone and a.GetSymbol() != "H"]
    return dict(mol=m, pos=pos, N=amino_N, CA=ca, C=carboxyl_C, side=side)


def kabsch(P, Q):
    """rotation+translation mapping P onto Q (both N x 3)."""
    Pc, Qc = P.mean(0), Q.mean(0)
    H = (P - Pc).T @ (Q - Qc)
    U, S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1, 1, d]) @ U.T
    return R, Pc, Qc


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--resnum", type=int, required=True)
    ap.add_argument("--new_name", required=True, help="3-letter PDB code for grafted residue")
    ap.add_argument("--smiles", required=True, help="free amino acid SMILES (-NH2, -COOH)")
    ap.add_argument("--tag", required=True)
    args = ap.parse_args()

    pep = read_pdb_atoms(PEP)
    target = [a for a in pep if a["resnum"] == args.resnum]
    bb_names = {"N", "CA", "C", "O"}
    crys_bb = {a["name"]: np.array([a["x"], a["y"], a["z"]]) for a in target if a["name"] in bb_names}
    assert {"N", "CA", "C"} <= set(crys_bb), f"res {args.resnum} missing backbone"

    res = build_residue(args.smiles)
    P = np.array([res["pos"][res["N"]], res["pos"][res["CA"]], res["pos"][res["C"]]])
    Q = np.array([crys_bb["N"], crys_bb["CA"], crys_bb["C"]])
    R, Pc, Qc = kabsch(P, Q)
    def xform(v):
        return R @ (v - Pc) + Qc

    # transplanted side-chain atoms for the target residue
    side_atoms = []
    elem_count = {}
    for idx in res["side"]:
        v = xform(res["pos"][idx])
        el = res["mol"].GetAtomWithIdx(idx).GetSymbol()
        elem_count[el] = elem_count.get(el, 0) + 1
        side_atoms.append({"record": "HETATM", "name": f"{el}{elem_count[el]}",
                           "resname": args.new_name, "chain": "A",
                           "resnum": args.resnum, "icode": " ",
                           "x": v[0], "y": v[1], "z": v[2], "elem": el})

    # build new peptide atom list; keep target residue's atoms CONTIGUOUS:
    # crystal backbone (N,CA,C,O) immediately followed by the new side chain.
    new_atoms = []
    side_inserted = False
    for a in pep:
        if a["resnum"] != args.resnum:
            new_atoms.append(a)
        elif a["name"] in bb_names:
            b = dict(a); b["resname"] = args.new_name
            new_atoms.append(b)
            if a["name"] == "C" and not side_inserted:  # after last backbone atom
                new_atoms.extend(side_atoms); side_inserted = True
    if not side_inserted:
        new_atoms.extend(side_atoms)

    # write grafted peptide alone
    pep_path = f"{OUT}/analogue_{args.tag}_peptide.pdb"
    with open(pep_path, "w") as f:
        for i, a in enumerate(new_atoms, 1):
            f.write(fmt_atom(i, a["name"], a["resname"], a["chain"], a["resnum"],
                             (a["x"], a["y"], a["z"]), a["elem"]))
        f.write("END\n")
    print("wrote", pep_path, f"({len(new_atoms)} atoms; +{len(res['side'])} side-chain)")

    # write full complex (peptide chain A + KRAS chain B)
    kras = read_pdb_atoms(KRAS)
    cx_path = f"{OUT}/analogue_{args.tag}_complex.pdb"
    with open(cx_path, "w") as f:
        serial = 1
        for a in new_atoms:
            f.write(fmt_atom(serial, a["name"], a["resname"], a["chain"], a["resnum"],
                             (a["x"], a["y"], a["z"]), a["elem"])); serial += 1
        f.write("TER\n")
        for a in kras:
            rec = a["record"]
            nm = f" {a['name']}" if len(a["name"]) < 4 else a["name"]
            f.write(f"{rec:<6}{serial:>5} {nm:<4} {a['resname']:>3} {a['chain']}"
                    f"{a['resnum']:>4}    {a['x']:8.3f}{a['y']:8.3f}{a['z']:8.3f}"
                    f"  1.00  0.00          {a['elem']:>2}\n"); serial += 1
        f.write("END\n")
    print("wrote", cx_path)


if __name__ == "__main__":
    main()
