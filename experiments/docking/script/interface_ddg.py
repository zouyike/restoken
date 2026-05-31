#!/usr/bin/env python3
"""Rosetta interface dG for an NCAA-peptide / protein complex whose peptide chain
uses custom polymer params (crystal-named, see name_param_to_pdb.py).

The peptide is a macrocycle with N-methyl residues, so a plain pose_from_file
crashes on auto-terminus patches. We therefore:
  - read the protein target chain with a standard pose_from_file (canonical),
  - build the peptide chain residue-by-residue (ResidueFactory + crystal xyz),
    which bypasses terminus patching,
  - join peptide to target by a jump and score the interface.

For a *fixed-backbone* interface dG, peptide-terminus and receptor-internal
terms (GDP/Mg, stripped here) cancel in the analogue-vs-parent ddG, and the
macrocyclic closure bond is not required.

Usage:
  interface_ddg.py --complex <complex.pdb> --params-dir <dir> [--params <extra.params>...]
                   [--pep-chain A] [--tgt-chain B] [--repack] [--out-pdb <p>]
"""
import argparse, glob, os, sys, tempfile
from pyrosetta import init, Pose
from pyrosetta.rosetta.core.import_pose import pose_from_file
from pyrosetta.rosetta.core.chemical import ChemicalManager, ResidueTypeFinder
from pyrosetta.rosetta.core.conformation import ResidueFactory
from pyrosetta.rosetta.numeric import xyzVector_double_t as V
from pyrosetta.rosetta.core.scoring import get_score_function
from pyrosetta.rosetta.protocols.analysis import InterfaceAnalyzerMover

AA3 = {"ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE","LEU",
       "LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL"}

def read_chain_atoms(pdb, chain):
    """ordered [(resseq, resname, {atomname: (x,y,z)})] for one chain."""
    res, order = {}, []
    for ln in open(pdb):
        if ln.startswith(("ATOM","HETATM")) and ln[21] == chain:
            nm = ln[12:16].strip(); rn = ln[17:20].strip(); sq = int(ln[22:26])
            xyz = (float(ln[30:38]), float(ln[38:46]), float(ln[46:54]))
            key = (sq, rn)
            if key not in res:
                res[key] = {}; order.append(key)
            res[key][nm] = xyz
    return [(sq, rn, res[(sq, rn)]) for (sq, rn) in order]

def write_protein_chain(pdb, chain, out):
    """canonical-AA ATOM records only (strip ligands/ions/waters/H)."""
    with open(out, "w") as fh:
        for ln in open(pdb):
            if ln.startswith(("ATOM","HETATM")) and ln[21] == chain \
               and ln[17:20].strip() in AA3:
                fh.write("ATOM  " + ln[6:])
        fh.write("TER\nEND\n")

def resolve_type(rts, rn):
    """ResidueType by full NAME, else by name3 (custom params keep a unique
    NAME like X48 while name3 matches the PDB resname A48)."""
    try:
        return rts.name_map(rn)
    except RuntimeError:
        rt = ResidueTypeFinder(rts).name3(rn).get_representative_type()
        if rt is None:
            raise KeyError(f"no ResidueType for resname/name3 {rn}")
        return rt

def build_peptide_pose(atoms_by_res, rts):
    pose = Pose()
    for k, (sq, rn, atoms) in enumerate(atoms_by_res):
        rt = resolve_type(rts, rn)
        res = ResidueFactory.create_residue(rt)
        for an, xyz in atoms.items():
            if res.has(an):
                res.set_xyz(an, V(*xyz))
        if k == 0:
            pose.append_residue_by_jump(res, 1)
        else:
            pose.append_residue_by_bond(res, False)
    return pose

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--complex", required=True)
    ap.add_argument("--params-dir", default=None)
    ap.add_argument("--params", nargs="*", default=[])
    ap.add_argument("--pep-chain", default="A")
    ap.add_argument("--tgt-chain", default="B")
    ap.add_argument("--repack", action="store_true")
    ap.add_argument("--out-pdb", default=None)
    a = ap.parse_args()

    ps = list(a.params)
    if a.params_dir:
        ps += glob.glob(os.path.join(a.params_dir, "*.params"))
    init("-beta_nov16 -mute all -load_PDB_components false -extra_res_fa " + " ".join(ps))
    rts = ChemicalManager.get_instance().residue_type_set("fa_standard")

    # receptor: standard read of the protein target chain
    tmp = tempfile.NamedTemporaryFile(suffix=".pdb", delete=False,
                                      dir="/public/home/genesis/.local/tmp").name
    write_protein_chain(a.complex, a.tgt_chain, tmp)
    pose = pose_from_file(tmp)
    n_tgt = pose.total_residue()

    # peptide: residue-by-residue build, then graft onto receptor by jump
    pep = build_peptide_pose(read_chain_atoms(a.complex, a.pep_chain), rts)
    pose.append_pose_by_jump(pep, n_tgt)   # jump connects tgt -> peptide
    pep_jump = pose.num_jump()

    sf = get_score_function()
    sf(pose)

    iam = InterfaceAnalyzerMover(pep_jump)
    iam.set_pack_separated(a.repack)
    iam.set_pack_input(a.repack)
    iam.set_compute_packstat(False)
    iam.set_calc_dSASA(True)
    iam.apply(pose)
    data = iam.get_all_data()

    print(f"complex          : {os.path.basename(a.complex)}")
    print(f"target_residues  : {n_tgt}   peptide_residues : {pep.total_residue()}")
    print(f"repack_separated : {a.repack}")
    print(f"dG_separated     : {iam.get_separated_interface_energy():.3f} REU")
    print(f"dSASA_interface  : {iam.get_interface_delta_sasa():.1f} A^2")
    print(f"dG_dSASA_ratio   : {iam.get_interface_dG()/max(iam.get_interface_delta_sasa(),1e-6)*100:.3f}")
    if a.out_pdb:
        pose.dump_pdb(a.out_pdb)
    os.remove(tmp)

if __name__ == "__main__":
    main()
