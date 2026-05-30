"""Decomposition of LUNA18 (a monocyclic, N-methyl-rich oral macrocycle).

Locks in two fixes:

1. No spurious crosslink. LUNA18 is monocyclic. A pendant dimethylamide side
   chain (a beta-amino-acid residue, atoms ...Cα-C(=O)NMe2) used to be misread
   as a second backbone carbonyl, fusing two residues and triggering a phantom
   crosslink. The backbone walk must treat an amide's own carbonyl as a hard
   barrier so the side-chain carbonyl cannot tunnel to the real backbone.

2. N-modification recognition. The backbone N-alkylation that drives LUNA18's
   permeability must be detected per residue and preserved on reconstruction,
   independent of which (mostly non-N-alkyl) library block matches the side
   chain.
"""

from rdkit import RDLogger, Chem
from rdkit.Chem import rdMolDescriptors

RDLogger.DisableLog("rdApp.*")

from restoken.src.decomposer import PeptideDecomposer
from restoken.src.library import BlockLibrary

LUNA18 = ("CC[C@H](C)[C@H]1C(=O)N([C@H](C(=O)N2CC[C@H]2C(=O)N([C@H](C(=O)N(CC(=O)N"
          "[C@H](C(=O)N3CCC[C@H]3C(=O)NC4(CCCC4)C(=O)N([C@H](C(=O)N([C@@H](CC(=O)N"
          "([C@H](C(=O)N1)CC(C)C)C)C(=O)N(C)C)C)C5CCCC5)C)CCC6=CC(=C(C(=C6)F)C(F)"
          "(F)F)F)C)CC7=CC=C(C=C7)C)CC)C)C")

_DEC = PeptideDecomposer(library=BlockLibrary())


def test_luna18_is_monocyclic_no_crosslink():
    res = _DEC.tokenize(LUNA18)
    assert res["n_crosslinks"] == 0, "LUNA18 is monocyclic; no crosslink expected"
    assert res["n_debris"] == 0, "no crosslink => no orphaned stub debris"
    assert res["n_residues"] == 11, f"expected 11 residues, got {res['n_residues']}"
    print("  OK: 11 residues, 0 crosslinks, 0 debris")


def test_luna18_nmod_recognized():
    res = _DEC.tokenize(LUNA18)
    labels = [r["nmod_label"] for r in res["residues"]]
    # 3 unmodified backbone N-H (matches LUNA18's HBD = 3), the rest N-modified
    assert labels.count("NH") == 3, f"expected 3 N-H residues, got {labels}"
    assert res["n_nmod"] == 8, f"expected 8 N-modified residues, got {res['n_nmod']}"
    assert "N-C2" in labels, "the N-ethyl residue must be recognized as N-C2"
    assert labels.count("NME") >= 4, "the N-methyl residues must be recognized"
    print(f"  OK: nmod labels = {labels}")


def test_luna18_nmod_aware_reconstruction_preserves_heavy_atoms():
    res = _DEC.tokenize(LUNA18)
    target = Chem.MolFromSmiles(LUNA18)
    t_heavy = target.GetNumHeavyAtoms()
    t_hbd = rdMolDescriptors.CalcNumHBD(target)

    plain = _DEC.reconstruct(res, nmod_correct=False)
    aware = _DEC.reconstruct(res, nmod_correct=True)
    assert plain and aware, "both reconstructions must assemble"

    pm, am = Chem.MolFromSmiles(plain), Chem.MolFromSmiles(aware)
    p_heavy, a_heavy = pm.GetNumHeavyAtoms(), am.GetNumHeavyAtoms()
    p_hbd, a_hbd = rdMolDescriptors.CalcNumHBD(pm), rdMolDescriptors.CalcNumHBD(am)

    # N-mod correction recovers the alkyl mass lost by plain block round-trip:
    # heavy-atom count becomes exact, and HBD moves toward the true 3.
    assert a_heavy == t_heavy, f"heavy atoms {a_heavy} != target {t_heavy}"
    assert a_heavy > p_heavy, "N-mod correction must add back N-alkyl atoms"
    assert a_hbd < p_hbd, "N-mod correction must remove spurious backbone N-H donors"
    # topology must stay exact
    assert (rdMolDescriptors.CalcNumRings(am)
            == rdMolDescriptors.CalcNumRings(target))
    print(f"  OK: heavy {p_heavy}->{a_heavy} (target {t_heavy}); "
          f"HBD {p_hbd}->{a_hbd} (target {t_hbd})")


if __name__ == "__main__":
    for fn in (test_luna18_is_monocyclic_no_crosslink,
               test_luna18_nmod_recognized,
               test_luna18_nmod_aware_reconstruction_preserves_heavy_atoms):
        fn()
        print("PASS", fn.__name__)
