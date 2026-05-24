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
- **Constraint satisfaction:** up to 123x enrichment over random baseline

See `experiments/` for the full benchmark suite.

## Citation

```
[Manuscript in preparation]
```

## License

MIT
