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


def compute_exotic_flags(block: Block, raw_entry: dict) -> set[str]:
    """Compute exotic flags for a block from SMILES structure.

    Flags (any triggered = exotic):
        high_heteroatom — heteroatom/heavy ratio > 0.5
        poly_ring — 3+ ring systems in the monomer
        high_mw — molecular weight > 300
        high_flex — rotatable bonds > 7
        halogenated — contains F/Cl/Br/I
        many_oxygens — more than 3 oxygen atoms
        many_nitrogens — more than 3 nitrogen atoms
        has_phosphorus — any phosphorus atom
        multi_charge — |charge| >= 2
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
    if n_rings >= _EXOTIC_THRESHOLDS["poly_ring"]:
        flags.add("poly_ring")

    n_heavy = mol.GetNumHeavyAtoms()
    n_hetero = sum(1 for a in mol.GetAtoms() if a.GetAtomicNum() not in (1, 6))
    if n_heavy > 0 and (n_hetero / n_heavy) > _EXOTIC_THRESHOLDS["high_heteroatom"]:
        flags.add("high_heteroatom")

    # Side-chain atom counts (subtract 2N + 2O backbone contribution)
    sc_O = sum(1 for a in mol.GetAtoms() if a.GetAtomicNum() == 8) - 2
    sc_N = sum(1 for a in mol.GetAtoms() if a.GetAtomicNum() == 7) - 2
    n_P = sum(1 for a in mol.GetAtoms() if a.GetAtomicNum() == 15)
    if sc_O > 3:
        flags.add("many_oxygens")
    if sc_N > 3:
        flags.add("many_nitrogens")
    if n_P > 0:
        flags.add("has_phosphorus")

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
