# ResToken

**Residue-semantic tokenization for LLM-based cyclic peptide design.**

ResToken is a curated library of 400 noncanonical amino acid (NCAA) building blocks encoded as compact semantic tokens, designed for large language model (LLM) cyclic peptide design.

## Key Features

- **400 building blocks** spanning alpha (250), beta (101), and gamma (49) backbones; L (292), D (86), and achiral (22) stereochemistry; 353 unmodified (`NO`, uppercase prefix) and 47 N-modified blocks (`NCY` 38, `NXX` 6, `NME` 3; lowercase prefix).
- **8 semantic properties** exposed to the LLM: functional class, chirality, charge label, backbone type, N-modification status, side-chain bulk, polarity, and flexibility (169 unique 8-property profiles across 400 blocks).
- **Dual-dictionary architecture**: LLM-facing semantic dictionary (`restoken/data/bb_dict_llm_v11.json`) separated from the structural backend dictionary (`restoken/data/bb_dict_backend_v11.json`, containing amide-capped AA-SMILES, SELFIES, InChIKey, and computed descriptors).
- **Constraint validator and cyclic assembler**: `SequenceValidator` (7 configurable rules) and `CyclicPeptideAssembler` (RDKit molecular-graph macrocyclization and 3D conformer generation).

## Installation & Environment

Tested with Python 3.11 (compatible with Python >= 3.10). Exact pinned package versions used for all reported results are provided in `requirements.txt` and `environment.yml`.

```bash
# Option 1: Conda environment (recommended for full reproduction)
conda env create -f environment.yml
conda activate restoken

# Option 2: Pip editable install with pinned dependencies
pip install -r requirements.txt
pip install -e .
```

## Quick Start

```python
from restoken.src.library import BlockLibrary
from restoken.src.validator import SequenceValidator
from restoken.src.cyclic_assembler import CyclicPeptideAssembler

lib = BlockLibrary()
block = lib.get("A01")
print(block.id, block.mc_type, block.chirality, block.aa_smiles)

validator = SequenceValidator()
print(validator.validate(["A01", "N15", "s01", "K08", "E09", "U04"]))

asm = CyclicPeptideAssembler()
smiles = asm.assemble("A09-a31-N03-a01-A05-A06")
print(smiles)
```

## Reproducing Manuscript Results

All benchmark scripts, prompt templates, assembled macrocycle SMILES, and evaluation outputs are organized under `experiments/`:

| Manuscript Section / Figure | Script | Outputs |
|---|---|---|
| Table 1 & Fig. 2 (Library statistics & chemical space) | `experiments/script/chemical_space_analysis.py`, `experiments/script/plot_figures.py` | `restoken/data/bb_dict_*_v11.json` |
| Table 2–3 & Fig. 3 (Tasks 1–3 benchmarks & statistics) | `experiments/script/run_benchmark.py`, `experiments/script/w2_summary.py`, `experiments/script/w3_statistical_analysis.py` | `experiments/output/exp1_validity/`, `exp2_controllability/`, `exp3_property_constrained/` |
| Fig. 4 (Property-specified design demos) | `restoken/examples/alogp_design_demo.py`, `restoken/examples/cpp_design_demo.py` | `restoken/examples/` |
| Fig. 5 (Independent permeability validation on CycPeptMPDB v1.2) | `experiments/permeability_validation/script/prep_data.py`, `train_perm.py`, `assemble_smiles.py`, `score_generations.py`, `applicability_domain.py` | `experiments/permeability_validation/data/assembled_smiles.csv`, `output/generation_scores.csv`, `output/comparison_stats.json`, `output/applicability_domain.json` |
| Fig. 6 (Cross-model open-weight LLM replication) | `experiments/script/launch_open_model_bench.py`, `experiments/permeability_validation/script/open_model_score.py` | `experiments/output/open_model_bench/`, `experiments/permeability_validation/output/open_model_scores.csv`, `open_model_perm_stats.json` |
| Leakage-free splits (Butina cluster & literature source), simple baselines (1-NN, Ridge, RF), PeptideCLM-2 comparison, and library degeneracy/bias analysis | `experiments/dd_revision/perm_splits/prepare_splits_and_baselines.py`, `train_and_score_split.py` | `experiments/dd_revision/perm_splits/` |
| Fair representation controls (Monomer-SMILES, 8-property matched HELM, TxGemma-9B-Chat) | `experiments/dd_revision/rep_fairness/run_fairness_bench.py` | `experiments/dd_revision/rep_fairness/` |
| Fine-tuning representation & held-out block addressability grids | `experiments/finetune_representation/` | `experiments/finetune_representation/output/` |

### Large Model Checkpoints & External Data

- **CycPeptMPDB v1.2** (accessed 2026-06-01; 7,451 total cyclic peptides, 6,889 PAMPA-labelled unique RDKit-valid entries after filtering): download `CycPeptMPDB_Peptide_All.csv` from http://cycpeptmpdb.com/ and place under `experiments/permeability_validation/data/`.
- **Pretrained backbones**: PeptideCLM (`aaronfeller/PeptideCLM-23M-all`) and PeptideCLM-2 (`aaronfeller/peptideclm-2-hybrid-small`) are publicly hosted on Hugging Face.
- **Fine-tuned permeability predictor checkpoints** (`perm_model.pt`, 169 MB, and leakage-free split checkpoints): archived in the Zenodo release accompanying the manuscript (exceeding GitHub's 100 MB single-file limit).

## License

MIT License. See `LICENSE`.
