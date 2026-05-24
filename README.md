# ResToken

**Residue-semantic tokenization for LLM-based noncanonical cyclic peptide design.**

ResToken is a curated library of 400 NCAA (noncanonical amino acid) building blocks encoded as semantic tokens, designed for use with large language models (LLMs) in cyclic peptide design.

## Key Features

- **400 building blocks** spanning alpha (250), beta (101), and gamma (49) backbone types
- **Semantic properties** — each token encodes chirality, charge, backbone type, N-modification, bulk, polarity, flexibility, and functional class
- **Dual-dictionary architecture** — LLM-safe dictionary (semantic properties only) + backend dictionary (full SMILES, SELFIES, InChIKey)
- **Constraint validator** — 7 configurable checks guarantee chemically legal sequences
- **SMILES reconstructor** — convert token sequences to full molecular structures (RDKit-validated)

## Installation

```bash
pip install -e .
```

## Quick Start

```python
from restoken.src.library import BlockLibrary
from restoken.src.validator import SequenceValidator
from restoken.src.reconstructor import SMILESReconstructor

# Load the 400-block library
lib = BlockLibrary()

# Look up a block
block = lib.get("A01")
print(f"{block.id}: {block.mc_type} backbone, {block.chirality} chirality, {block.aa_smiles}")

# Validate a sequence
validator = SequenceValidator()
result = validator.validate(["A01", "N15", "s01", "K08", "E09", "U04"])

# Reconstruct SMILES
recon = SMILESReconstructor()
info = recon.reconstruct(["A01", "N15", "s01", "K08", "E09", "U04"])
print(info["residue_smiles"])
```

## Library Statistics

| Property | Alpha | Beta | Gamma | Total |
|---|---|---|---|---|
| Count | 250 | 101 | 49 | **400** |
| L-chirality | 219 | 44 | 29 | 292 |
| D-chirality | 19 | 51 | 16 | 86 |
| Neutral | 215 | 85 | 47 | 347 |
| N-cyclic | 38 | 4 | 3 | 45 |
| Functional classes | 18 | 20 | 10 | 21 |

## Benchmark

Zero-shot generation across 5 LLMs (Gemini 2.5 Pro/Flash, GPT-4o, Qwen 3.5 9B, Gemma 3 12B):

- **Validity:** ResToken 90.4% vs SMILES 69.0% vs HELM 84.5%
- **Frozen compliance:** 100% across all models
- **Constraint satisfaction:** up to 141x enrichment over random baseline

See `experiments/` for the full benchmark suite.

## Cyclic Peptide Assembly

Convert any ResToken sequence to a macrocyclic peptide SMILES and render 2D structures:

```python
from restoken.src.cyclic_assembler import CyclicPeptideAssembler

asm = CyclicPeptideAssembler()

# Assemble to SMILES
smiles = asm.assemble("A09-a31-N03-a01-A05-A06")

# Render 2D structure
asm.render("A09-a31-N03-a01-A05-A06", "peptide.png")

# Render multiple peptides in a grid
sequences = ["A01-A09-A03-A07-A05-A06", "A09-a31-N03-a01-A05-A06"]
asm.render_grid(sequences, "grid.png", cols=2)

# Get aggregate properties
props = asm.get_properties("A09-a31-N03-a01-A05-A06")
# {'n_residues': 6, 'net_charge': 0, 'total_hbd': 0, ...}
```

3D conformer generation (ETKDG with macrocycle torsion sampling + MMFF optimization):

```python
# Generate 3D structure (SDF, PDB, or MOL)
result = asm.generate_3d("A09-a31-N03-a01-A05-A06", "peptide.sdf")
# {'smiles': '...', 'energy': 94.4, 'n_heavy_atoms': 50, ...}

result = asm.generate_3d("A09-a31-N03-a01-A05-A06", "peptide.pdb", n_confs=100)
```

CLI usage:
```bash
# Render 2D structure
python -m restoken.src.cyclic_assembler A09-a31-N03-a01-A05-A06 -o output.png

# Generate 3D conformer
python -m restoken.src.cyclic_assembler A09-a31-N03-a01-A05-A06 --mode-3d -o output.sdf

# SMILES only
python -m restoken.src.cyclic_assembler A09-a31-N03-a01-A05-A06 --smiles-only
```

Supports all 400 block types: alpha/beta/gamma backbones, L/D/achiral chirality, N-methylation (NME), and N-cyclic rings (NCY, proline-like).

## LLM-to-Molecule Design Pipeline

ResToken serves as an interface between LLMs and real chemistry. See `restoken/examples/alogp_design_demo.py` for an end-to-end example:

```bash
export GEMINI_API_KEY=your_key
python restoken/examples/alogp_design_demo.py --alogp_min 1.5 --alogp_max 3.25
```

Pipeline:
1. **Prompt** an LLM with the ResToken dictionary + property constraints (e.g., AlogP 1.5–3.25)
2. **Validate** — block IDs checked against the 400-block library (no hallucinated monomers)
3. **Assemble** — sequences converted to macrocyclic SMILES via RDKit molecular graph operations
4. **Filter** — compute AlogP/MW/TPSA from the SMILES, select hits
5. **Render** — 2D ChemDraw-style structures of the final candidates

## Citation

```
[Manuscript in preparation]
```

## License

MIT
