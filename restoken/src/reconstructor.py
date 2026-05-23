"""SMILES reconstructor: token sequence → full molecular SMILES.

Reconstructs peptide SMILES from ResToken sequences using the backend
dictionary. Supports both linear and cyclic peptide output.
"""

from typing import Optional

from restoken.src.library import BlockLibrary


class SMILESReconstructor:
    """Reconstruct full SMILES from ResToken ID sequences.

    Each block's AA_SMILES represents the residue as an amide-capped fragment:
        R-[C@H](NC(C)=O)C(=O)NC
    i.e., AcNH-[alpha-carbon(R)]-CONHCH3 for alpha amino acids.

    For cyclic peptide reconstruction, we:
    1. Look up each block's full AA_SMILES
    2. Return per-residue SMILES (for fingerprinting / property computation)
    3. Optionally concatenate into a linear peptide SMILES (approximate)

    Note: True macrocyclic SMILES requires ring-closure chemistry that depends
    on the specific cyclization strategy (head-to-tail, disulfide, thioether).
    This module provides the residue-level building blocks; full macrocycle
    assembly is left to dedicated cyclization tools.
    """

    def __init__(self, library: Optional[BlockLibrary] = None):
        self.lib = library or BlockLibrary()

    def residue_smiles(self, sequence: list[str]) -> list[str]:
        """Return the AA_SMILES for each block in the sequence."""
        result = []
        for tid in sequence:
            block = self.lib.get(tid)
            if block is None:
                raise KeyError(f"Unknown block ID: {tid}")
            result.append(block.aa_smiles)
        return result

    def sidechain_smiles(self, sequence: list[str]) -> list[str]:
        """Return the side-chain SMILES for each block."""
        result = []
        for tid in sequence:
            block = self.lib.get(tid)
            if block is None:
                raise KeyError(f"Unknown block ID: {tid}")
            result.append(block.sc_smiles)
        return result

    def reconstruct(self, sequence: list[str]) -> dict:
        """Full reconstruction: per-residue SMILES + aggregate info.

        Returns a dict with:
          - residue_smiles: list of AA_SMILES per position
          - sidechain_smiles: list of SC_SMILES per position
          - block_ids: the input sequence
          - properties: per-residue property summary
        """
        blocks = []
        for tid in sequence:
            block = self.lib.get(tid)
            if block is None:
                raise KeyError(f"Unknown block ID: {tid}")
            blocks.append(block)

        return {
            "block_ids": sequence,
            "residue_smiles": [b.aa_smiles for b in blocks],
            "sidechain_smiles": [b.sc_smiles for b in blocks],
            "properties": [
                {
                    "id": b.id,
                    "class": b.aa_class,
                    "chirality": b.chirality,
                    "charge": b.charge,
                    "mc_type": b.mc_type,
                    "mc_nmod": b.mc_nmod,
                    "sc_bulk": b.sc_bulk,
                    "hbd": b.sc_hbd,
                    "hba": b.sc_hba,
                    "rot": b.rot_total,
                }
                for b in blocks
            ],
        }

    def validate_smiles(self, sequence: list[str]) -> dict:
        """Reconstruct and validate each SMILES with RDKit.

        Returns dict with valid/invalid counts and any failures.
        Requires rdkit to be importable.
        """
        from rdkit import Chem

        aa_smiles = self.residue_smiles(sequence)
        results = []
        for tid, smi in zip(sequence, aa_smiles):
            mol = Chem.MolFromSmiles(smi)
            results.append({
                "id": tid,
                "smiles": smi,
                "valid": mol is not None,
                "canonical": Chem.MolToSmiles(mol) if mol else None,
            })

        n_valid = sum(1 for r in results if r["valid"])
        return {
            "total": len(results),
            "valid": n_valid,
            "invalid": len(results) - n_valid,
            "details": results,
        }
