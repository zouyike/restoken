"""Block library loader for the ResToken 400-block NCAA library."""

import json
from pathlib import Path
from dataclasses import dataclass, field
from typing import Optional
from functools import lru_cache

_DATA_DIR = Path(__file__).resolve().parent.parent / "data"
_LLM_DICT = _DATA_DIR / "bb_dict_llm_v11.json"
_BACKEND_DICT = _DATA_DIR / "bb_dict_backend_v11.json"


@dataclass(frozen=True)
class Block:
    id: str
    helm: str
    aa_class: str
    chirality: str
    charge: int
    charge_label: str
    mc_type: str
    mc_nmod: str
    mc_len: str
    sc_bulk: str
    sc_hbd: int
    sc_hba: int
    rot_total: int
    aa_smiles: str
    sc_smiles: str
    polarity_bin: str = "?"
    flex_bin: str = "?"


class BlockLibrary:
    """Loads and indexes the 400-block NCAA library from JSON dictionaries."""

    def __init__(
        self,
        backend_path: Optional[str] = None,
        llm_path: Optional[str] = None,
    ):
        bp = Path(backend_path) if backend_path else _BACKEND_DICT
        lp = Path(llm_path) if llm_path else _LLM_DICT

        with open(bp) as f:
            self._backend = json.load(f)
        with open(lp) as f:
            self._llm = json.load(f)

        self._blocks: dict[str, Block] = {}
        self._build_index()

    def _build_index(self):
        llm_lookup = {e["id"]: e for e in self._llm["blocks"]}
        for b in self._backend["blocks"]:
            llm_entry = llm_lookup.get(b["id"], {})
            block = Block(
                id=b["id"],
                helm=b["helm"],
                aa_class=b["aa_class"],
                chirality=b["chirality"],
                charge=b["charge"],
                charge_label=b["charge_label"],
                mc_type=b["mc"]["type"],
                mc_nmod=b["mc"]["nmod"],
                mc_len=b["mc"]["len"],
                sc_bulk=b["sc"]["bulk"],
                sc_hbd=b["sc"]["hbd"],
                sc_hba=b["sc"]["hba"],
                rot_total=b["rot_total"],
                aa_smiles=b["structure"]["aa_smiles"],
                sc_smiles=b["sc"]["smiles"],
                polarity_bin=llm_entry.get("polarity_bin", "?"),
                flex_bin=llm_entry.get("flex_bin", "?"),
            )
            self._blocks[block.id] = block

    @property
    def version(self) -> str:
        return self._backend["version"]

    @property
    def size(self) -> int:
        return len(self._blocks)

    @property
    def all_ids(self) -> set[str]:
        return set(self._blocks.keys())

    def __getitem__(self, block_id: str) -> Block:
        return self._blocks[block_id]

    def __contains__(self, block_id: str) -> bool:
        return block_id in self._blocks

    def get(self, block_id: str) -> Optional[Block]:
        return self._blocks.get(block_id)

    def blocks_by_property(self, **filters) -> list[Block]:
        """Filter blocks by property values.

        Example: lib.blocks_by_property(mc_type="beta", charge_label="neu")
        """
        result = []
        for b in self._blocks.values():
            match = True
            for k, v in filters.items():
                if getattr(b, k, None) != v:
                    match = False
                    break
            if match:
                result.append(b)
        return result

    def get_llm_entry(self, block_id: str) -> Optional[dict]:
        for entry in self._llm["blocks"]:
            if entry["id"] == block_id:
                return entry
        return None

    def llm_dict_text(self) -> str:
        """Return the full LLM dictionary as formatted text for prompt injection."""
        lines = []
        for entry in self._llm["blocks"]:
            tags = ", ".join(entry.get("tags", []))
            lines.append(
                f"{entry['id']} | {entry['class']} | {entry['chirality']} | "
                f"{entry['charge_label']} | {entry['mc_type']} | {entry['mc_nmod']} | "
                f"{entry['sc_bulk']} | {entry.get('polarity_bin','?')} | "
                f"{entry.get('flex_bin','?')} | [{tags}]"
            )
        header = "ID | Class | Chiral | Charge | Backbone | N-mod | Bulk | Polarity | Flexibility | Tags"
        sep = "-" * len(header)
        return f"{header}\n{sep}\n" + "\n".join(lines)


# ── Exotic flag computation ────────────────────────────────────────────

_EXOTIC_THRESHOLDS = {
    "high_heteroatom": 0.5,
    "poly_ring": 3,
    "high_mw": 300,
    "high_flex": 8,
}

# Canonical amino acid side-chain heteroatom budgets (N, O).
# Blocks within budget+1 for N or budget+3 for O are not flagged.
_CANONICAL_SC_BUDGET = {
    'R': (3, 0),  # Arg: guanidinium
    'H': (2, 0),  # His: imidazole
    'W': (1, 0),  # Trp: indole
    'K': (1, 0),  # Lys: amine
    'N': (1, 1),  # Asn: amide
    'Q': (1, 1),  # Gln: amide
    'D': (0, 2),  # Asp: carboxyl
    'E': (0, 2),  # Glu: carboxyl
    'S': (0, 1),  # Ser: hydroxyl
    'T': (0, 1),  # Thr: hydroxyl
    'Y': (0, 1),  # Tyr: phenol
    'C': (0, 0),
    'M': (0, 0),
}

# Flags indicating synthesis difficulty (used for exclude_exotic filtering).
# high_flex is informational only — flexibility doesn't make synthesis harder.
SYNTHESIS_EXOTIC_FLAGS = frozenset({
    "halogenated", "multi_charge", "high_mw", "poly_ring",
    "high_heteroatom", "excess_oxygen", "excess_nitrogen",
    "has_phosphorus", "geminal_hetero", "strained_aminal", "aminal",
})


def _has_geminal_heteroatoms(mol) -> bool:
    """Check if any sp3 carbon has 2+ heteroatoms attached via single bonds.

    Excludes:
      - Ring-internal patterns (e.g., N-C-N in proline-like heterocycles)
      - Carbons with a double bond to a heteroatom (guanidinium, amidine,
        carbamate, urea — common pharmacophoric groups, not exotic)
    """
    ri = mol.GetRingInfo()
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 6:
            continue
        has_double_to_hetero = False
        hetero_single = []
        for nb in atom.GetNeighbors():
            if nb.GetAtomicNum() in (1, 6):
                continue
            bond = mol.GetBondBetweenAtoms(atom.GetIdx(), nb.GetIdx())
            if bond and bond.GetBondTypeAsDouble() > 1.0:
                has_double_to_hetero = True
                break
            if bond and bond.GetBondTypeAsDouble() == 1.0:
                hetero_single.append(nb.GetIdx())
        if has_double_to_hetero:
            continue
        if len(hetero_single) >= 2:
            atom_in_ring = ri.NumAtomRings(atom.GetIdx()) > 0
            all_in_same_ring = atom_in_ring and all(
                ri.NumAtomRings(h) > 0 for h in hetero_single
            )
            if not all_in_same_ring:
                return True
    return False


def _has_strained_ring_aminal(mol) -> bool:
    """Check for a carbon in a small ring (3- or 4-membered) bonded to 2+ ring
    heteroatoms via single bonds — a strained cyclic aminal / N,O-acetal.

    These motifs (e.g. 1,3-diazetidine, oxazetidine) are hydrolytically
    unstable: ring strain plus the labile N-C-N / N-C-O carbon makes them
    synthetically impractical. Distinct from geminal_hetero, which explicitly
    passes ring-internal patterns to allow stable proline-like (5+ membered)
    heterocycles.
    """
    ri = mol.GetRingInfo()
    small_rings = [set(r) for r in ri.AtomRings() if len(r) <= 4]
    if not small_rings:
        return False
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 6:
            continue
        idx = atom.GetIdx()
        if not any(idx in r for r in small_rings):
            continue
        ring_hetero = 0
        for nb in atom.GetNeighbors():
            if nb.GetAtomicNum() in (1, 6):
                continue
            bond = mol.GetBondBetweenAtoms(idx, nb.GetIdx())
            if not bond or bond.GetBondTypeAsDouble() != 1.0:
                continue
            if any(idx in r and nb.GetIdx() in r for r in small_rings):
                ring_hetero += 1
        if ring_hetero >= 2:
            return True
    return False


@lru_cache(maxsize=1)
def _aminal_pattern():
    from rdkit import Chem
    # sp3 carbon (CX4, four single bonds) single-bonded to two nitrogens.
    return Chem.MolFromSmarts("[#7]-[CX4](-[#7])")


def _has_aminal(mol) -> bool:
    """Check for any aminal carbon: an sp3 carbon single-bonded to two nitrogen
    atoms (N-C-N), regardless of ring size or ring membership.

    Aminals are the diamine analog of acetals; the N-C-N carbon is
    hydrolytically labile and these blocks are synthetically impractical. The
    CX4 constraint (four single bonds) excludes sp2 amidine / guanidinium / urea
    carbons (C=N or C=O) — stable pharmacophores, not aminals — so Arg-type
    guanidinium blocks are never flagged. Broader than strained_aminal (which is
    limited to <=4-membered rings): per directive, ALL aminal-bearing NCAAs are
    exotic.
    """
    patt = _aminal_pattern()
    return patt is not None and mol.HasSubstructMatch(patt)


def compute_exotic_flags(block: Block, raw_entry: dict) -> set[str]:
    """Compute exotic flags for a block from SMILES structure.

    Flags (any triggered = exotic):
        high_heteroatom — heteroatom/heavy ratio > 0.5
        poly_ring — more than 3 ring systems in the monomer
        high_mw — molecular weight > 300
        high_flex — rotatable bonds > 7 (informational only)
        halogenated — contains F/Cl/Br/I
        excess_oxygen — side-chain O exceeds canonical budget + 3
        excess_nitrogen — side-chain N exceeds canonical budget + 1
        has_phosphorus — any phosphorus atom
        multi_charge — |charge| >= 2
        geminal_hetero — carbon with 2+ heteroatoms via single bonds (non-ring)
        strained_aminal — aminal/acetal carbon inside a 3- or 4-membered ring
        aminal — any sp3 carbon single-bonded to two nitrogens (N-C-N)
    """
    flags = set()

    if raw_entry.get("halogen", 0) > 0:
        flags.add("halogenated")

    if block.rot_total > 7:
        flags.add("high_flex")

    if abs(block.charge) >= 2:
        flags.add("multi_charge")

    smi = block.aa_smiles
    if not smi:
        return flags

    try:
        from rdkit import Chem
        from rdkit.Chem import Descriptors, rdMolDescriptors
    except ImportError:
        return flags

    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        return flags

    mw = Descriptors.ExactMolWt(mol)
    if mw > _EXOTIC_THRESHOLDS["high_mw"]:
        flags.add("high_mw")

    n_rings = rdMolDescriptors.CalcNumRings(mol)
    if n_rings > _EXOTIC_THRESHOLDS["poly_ring"]:
        flags.add("poly_ring")

    n_heavy = mol.GetNumHeavyAtoms()
    n_hetero = sum(1 for a in mol.GetAtoms() if a.GetAtomicNum() not in (1, 6))
    if n_heavy > 0 and (n_hetero / n_heavy) > _EXOTIC_THRESHOLDS["high_heteroatom"]:
        flags.add("high_heteroatom")

    # Side-chain heteroatom counts (subtract 2N + 2O backbone/cap contribution)
    total_O = sum(1 for a in mol.GetAtoms() if a.GetAtomicNum() == 8)
    total_N = sum(1 for a in mol.GetAtoms() if a.GetAtomicNum() == 7)
    sc_O = total_O - 2
    sc_N = total_N - 2
    n_P = sum(1 for a in mol.GetAtoms() if a.GetAtomicNum() == 15)

    # Budget-aware thresholds: canonical AA analogs get more headroom
    budget_N, budget_O = _CANONICAL_SC_BUDGET.get(block.aa_class, (0, 0))
    if sc_N > budget_N + 1:
        flags.add("excess_nitrogen")
    if sc_O > budget_O + 3:
        flags.add("excess_oxygen")
    if n_P > 0:
        flags.add("has_phosphorus")

    # Geminal heteroatoms: 2+ heteroatoms single-bonded to same carbon (non-ring)
    if _has_geminal_heteroatoms(mol):
        flags.add("geminal_hetero")

    # Strained cyclic aminal/acetal: aminal carbon inside a 3- or 4-membered ring
    if _has_strained_ring_aminal(mol):
        flags.add("strained_aminal")

    # Aminal: any sp3 carbon single-bonded to two nitrogens (N-C-N), ring or not.
    if _has_aminal(mol):
        flags.add("aminal")

    return flags


def compute_all_exotic_flags(library: "BlockLibrary", raw_blocks: dict) -> dict[str, set[str]]:
    """Compute exotic flags for all blocks in the library.

    Returns:
        dict mapping block_id -> set of exotic flag strings (empty = not exotic)
    """
    result = {}
    for block_id in library.all_ids:
        block = library[block_id]
        raw = raw_blocks.get(block_id, {})
        result[block_id] = compute_exotic_flags(block, raw)
    return result
