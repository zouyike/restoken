"""
Prepare the LUNA18.KRAS-G12D docking template from PDB 7YV1.

The Fab (chains H, L) is a crystallization chaperone: LUNA18 makes 35 heavy-atom
contacts (<4 A) with KRAS and ZERO with the Fab, so H/L are stripped. Waters are
dropped. Output adopts the CycSeqSeek pyrose_fr convention:
    chain A = peptide (LUNA18)   <- the design being scored
    chain B = KRAS-G12D + GDP + MG   <- the target

Writes:
    input/parent_complex.pdb   peptide(A) + KRAS/GDP/MG(B)
    input/receptor_kras.pdb    KRAS/GDP/MG only (B)
    input/luna18_crystal.pdb   crystal LUNA18 only (A)
"""
import gemmi, sys

src = "/scratch/genesis/NCAA_tokenization/experiments/docking/structures/7yv1.cif"
out = "/scratch/genesis/NCAA_tokenization/experiments/docking/input"

st = gemmi.read_structure(src)
st.setup_entities()
m = st[0]

WATER = {"HOH", "DOD", "WAT"}

def residues(chain_name, keep_water=False):
    for ch in m:
        if ch.name != chain_name:
            continue
        for r in ch:
            if not keep_water and r.name in WATER:
                continue
            yield r

# contact sanity check: LUNA18 (I) vs KRAS (A) vs Fab (H/L)
ns = gemmi.NeighborSearch(m, st.cell, 4.5).populate()
def count_contacts(pep_chain, tgt_chains, cutoff=4.0):
    n = 0
    for ch in m:
        if ch.name != pep_chain: continue
        for r in ch:
            if r.name in WATER: continue
            for a in r:
                marks = ns.find_atoms(a.pos, '\0', radius=cutoff)
                for mk in marks:
                    cra = mk.to_cra(m)
                    if cra.chain.name in tgt_chains and cra.residue.name not in WATER:
                        n += 1
                        break
                else:
                    continue
    return n

c_kras = count_contacts("I", {"A"})
c_fab  = count_contacts("I", {"H", "L"})
print(f"LUNA18 atom-level contacts (<4A): KRAS={c_kras}  Fab={c_fab}")

def write_pdb(path, specs):
    """specs: list of (src_chain, new_chain, keep_water)"""
    new = gemmi.Structure()
    new.cell = st.cell
    new.spacegroup_hm = st.spacegroup_hm
    model = gemmi.Model("1")
    for src_chain, new_chain, kw in specs:
        ch = gemmi.Chain(new_chain)
        for r in residues(src_chain, keep_water=kw):
            ch.add_residue(r)
        model.add_chain(ch)
    new.add_model(model)
    new.setup_entities()
    new.write_pdb(path)
    print("wrote", path)

write_pdb(f"{out}/parent_complex.pdb",
          [("I", "A", False), ("A", "B", False)])
write_pdb(f"{out}/receptor_kras.pdb",
          [("A", "B", False)])
write_pdb(f"{out}/luna18_crystal.pdb",
          [("I", "A", False)])
