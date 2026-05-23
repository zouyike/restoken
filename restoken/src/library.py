"""Block library loader for the ResToken 400-block NCAA library."""

import json
from pathlib import Path
from dataclasses import dataclass, field
from typing import Optional

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
        for b in self._backend["blocks"]:
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
