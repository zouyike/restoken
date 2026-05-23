# NCAA Tokenization — Lab Notebook

## Project Overview
**Goal:** JCIM Application Note — "ResToken: A Residue-Semantic Token Library Enabling LLM-Based Design of Noncanonical Cyclic Peptides"
**Target journal:** JCIM Application Note (~3000 words, 1-2 figures, 1 table)
**Started:** 2026-05-23
**Google Drive:** `D:\Gdrive\Zou\Project\邹一可_2026_LLM for NCAA tokenization\`

## Current Status
- [x] Project folder created (`/scratch/genesis/NCAA_tokenization/`)
- [x] Existing materials imported from Google Drive (v11 dictionaries, design docs, reference papers, architecture diagrams)
- [x] Outline drafted (`manuscript/drafts/JCIM_outline.md`)
- [x] 4-week experiment plan (`Lab_Notebook/experiment_plan_4week.md`)
- [x] **W1.1** ResToken validator (`restoken/src/validator.py`) — all checks implemented + tested
- [x] **W1.2** SMILES reconstructor (`restoken/src/reconstructor.py`) — 400/400 SMILES validated
- [x] **W1.3** Prompt templates (5 templates in `experiments/input/prompts/`) + builder script
- [x] **W1.4** Random baseline generator (`experiments/script/random_baseline.py`) — 4 profiles generated
- [ ] Run LLM benchmarks (ResToken vs SMILES vs HELM) — Week 2
- [ ] Prepare figures — Week 3
- [ ] Draft manuscript — Week 3-4
- [ ] Multi-agent review — Week 4
- [ ] Format for JCIM submission — Week 4

## Assets Inventory

### Data (restoken/data/)
| File | Description | Version |
|---|---|---|
| residue_tokens.csv | Full 400-block library with all properties + SMILES | v11 |
| bb_dict_llm_v11.json | LLM-safe dictionary (semantic properties only, no SMILES) | v11 |
| bb_dict_backend_v11.json | Backend dictionary (full SMILES + all computed properties) | v11 |

### Library Statistics (v11)
- Total blocks: 400
- Backbone: alpha (250), beta (101), gamma (49)
- Classes: I (88), F (46), A (46), T (42), P (38), ncaa (17), V (16), S (15), C (14), G (12), H (12), K (12), L (11), D (7), Y (7), E (6), N (4), Q (4), R (3)
- Chirality: L (292), D (86), achiral (22)
- Bulk: large (223), med (54), small (123)
- N-modification: NO (354), NCY (38), NXX (6), NME (2)

### Reference Papers (references/papers/)
1. 2025-Butcher — De novo design of all-atom biomolecules
2. 2025-Rettie — Accurate de novo design (two versions)
3. 2025-Wang — PepThink-R1: LLM for interpretable cyclic peptide
4. 2026-Liu — Geometric foundation model for enzymes

### Design Documents (references/)
- 顶刊研究计划.docx — Original Nature/Science-level research plan (includes self-critique as Reviewer #2)
- LLM规范.docx — LLM-safe block dictionary schema design rationale
- 工程步骤.docx — Engineering roadmap (Step 0-5)
- 如何让LLM生成器生成的更好.docx — Design philosophy for LLM generator improvement
- LLM for cycpep.drawio — Architecture diagram source

---

### Code (restoken/src/)
| File | Description |
|---|---|
| library.py | BlockLibrary loader — indexes 400 blocks from backend JSON, property filtering |
| validator.py | SequenceValidator — 7 checks: ID existence, length, charge, HBD, HBA, rot bonds, backbone compat |
| reconstructor.py | SMILESReconstructor — token seq → per-residue SMILES, RDKit round-trip validation |

### Prompt Templates (experiments/input/prompts/)
| File | Description | ~Tokens |
|---|---|---|
| restoken_zeroshot.txt | Full dictionary + constraints + format instructions | ~11K |
| smiles_zeroshot.txt | NCAA SMILES list + cyclization task | ~8.7K |
| helm_zeroshot.txt | HELM monomer codes + properties | ~6.8K |
| edit_frozen_restoken.txt | EDIT/FROZEN controllability prompt | varies |
| property_constrained.txt | Permeability/charge/rigidity-aware design | varies |

### Random Baselines (experiments/output/random_baseline/)
| Profile | Length | Acceptance Rate | N generated |
|---|---|---|---|
| unconstrained | 6 | 100% | 1000 |
| permeable | 6 | 2.07% | 1000 |
| charged_binder | 6 | 1.82% | 1000 |
| rigid_scaffold | 6 | 0.81% | 1000 |

---

## Experiment Log

### 2026-05-23 — W1: Build + Validate the Tool

#### W1.1 Library Loader + Validator
- Built `restoken/src/library.py`: `BlockLibrary` class loads both backend and LLM dictionaries, indexes by ID, supports property-based filtering
- Built `restoken/src/validator.py`: `SequenceValidator` with 7 configurable checks:
  1. Token ID existence (all IDs must be in 400-block library)
  2. Sequence length (configurable, default 6/8/10)
  3. Net charge = target
  4. Total HBD ≤ threshold
  5. Total HBA ≤ threshold
  6. Total rotatable bonds ≤ threshold
  7. Backbone compatibility (permissive/strict modes; strict disallows gamma-gamma and consecutive NCY)
- **Test results**: 100/100 legal sequences pass, 100/100 illegal sequences correctly rejected
- Constraint-specific tests: charge, HBD, backbone strict mode all working

#### W1.2 SMILES Reconstructor
- Built `restoken/src/reconstructor.py`: token sequence → per-residue AA_SMILES and SC_SMILES
- **RDKit round-trip validation: 400/400 blocks produce valid parseable SMILES** — zero failures
- Returns structured property dict per residue for downstream analysis

#### W1.3 Prompt Templates
- Created 5 prompt templates with placeholder variables for programmatic filling
- `build_prompts.py` fills templates with actual library data
- Token counts confirmed: all fit within 128K context windows
- Design decision: full 400-block dictionary in prompt (no chunking needed — only ~11K tokens)

#### W1.4 Random Baseline
- Rejection sampling with 4 constraint profiles
- Key finding: constrained profiles have 0.8-2.1% acceptance — validates that the design task is non-trivial
- Initial `rigid_scaffold` had max_rot=8, which was infeasible (0/5M attempts). Adjusted to max_rot=18 → 0.81% acceptance. This is a realistic constraint for the paper.
- All 4×1000 baseline sequences generated with diversity metrics (Morgan FP Tanimoto distance ~0.72-0.76)

#### Library Statistics (updated from actual data)
- mc_nmod: NO (353), NCY (45), NME (2) — note: 45 NCY, not 38 as originally reported
- charge_label: neu (347), neg (20), pos (30), zwi (3) — 3 zwitterionic blocks found
- class: includes M (8) and W (7) not listed in original inventory
- rot_total range: 1-12 per block
- HBD range: 0-10, HBA range: 0-14
