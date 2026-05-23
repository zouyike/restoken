# NCAA Tokenization — 4-Week Experiment Plan

**Paper:** JCIM Application Note — ResToken
**Start:** 2026-05-26 (Monday)
**Target completion:** 2026-06-20
**Parallel with:** RNF126 pipeline (ongoing), daily tasks

---

## Week 1 (May 26–30): Build + Validate the Tool

### W1.1 — Package the ResToken validator (Mon–Tue)
- **Goal:** A working Python module that validates LLM-output token sequences
- **Deliverable:** `restoken/src/validator.py`
- **Checks to implement:**
  - Token ID existence (is every ID in the 400-block library?)
  - Sequence length (6, 8, or 10 residues)
  - Net charge = target (sum of per-residue charges)
  - Total HBD ≤ threshold
  - Total HBA ≤ threshold
  - Total rotatable bonds ≤ threshold
  - Backbone compatibility (no illegal transitions between alpha/beta/gamma)
- **Test:** validate against 100 hand-written legal + 100 deliberately illegal sequences

### W1.2 — SMILES reconstructor (Tue–Wed)
- **Goal:** token sequence → full molecular SMILES via backend dictionary lookup
- **Deliverable:** `restoken/src/reconstructor.py`
- **Test:** reconstruct 50 random sequences, verify with RDKit `Chem.MolFromSmiles()` + canonical SMILES round-trip
- **Bonus:** compute Morgan fingerprints for downstream diversity metrics

### W1.3 — Prompt templates (Wed–Thu)
- **Goal:** Standardized prompts for each representation × model combination
- **Deliverable:** `experiments/input/prompts/` with:
  - `restoken_zeroshot.txt` — dictionary + constraints + format instructions
  - `smiles_zeroshot.txt` — equivalent task framed in SMILES (provide NCAA SMILES list)
  - `helm_zeroshot.txt` — equivalent task framed in HELM notation
  - `edit_frozen_restoken.txt` — EDIT/FROZEN controllability prompt
  - `property_constrained.txt` — permeability-aware design prompt
- **Key design decisions:**
  - ResToken prompt: full dictionary (400 blocks as text table) vs chunked (class-indexed)?
  - How much context fits in each model's window? Pre-check token counts.
  - SMILES prompt: provide all 400 NCAA SMILES or just a subset? (Full list ~30K tokens — may not fit)

### W1.4 — Baseline generator (Thu–Fri)
- **Goal:** Random baseline that samples from the 400-block library with constraint filtering
- **Deliverable:** `experiments/script/random_baseline.py`
- **Method:** rejection sampling — random sequence → validate → keep/reject
- **Output:** 1000 valid random sequences with property statistics
- **Purpose:** establishes the "how hard is it to randomly generate valid sequences?" baseline

### W1 Checkpoint
- [x] Validator passes unit tests (100 legal + 100 illegal) ✓ 2026-05-23
- [x] Reconstructor produces valid SMILES for all 400 blocks ✓ 400/400 validated
- [x] Prompt templates reviewed and token-counted ✓ ResToken ~11K, SMILES ~8.7K, HELM ~6.8K
- [x] Random baseline generates 1000 valid sequences ✓ 4 profiles × 1000 sequences
- [ ] All code committed to project folder

---

## Week 2 (Jun 2–6): LLM Benchmarks

### W2.1 — Experiment 1: Zero-Shot Generation Validity (Mon–Wed)

**Setup:**
| Model | API/Local | Context Window | Proxy needed |
|---|---|---|---|
| GPT-4o | API | 128K | Yes |
| Claude Sonnet 4.6 | API | 200K | Yes |
| Qwen-2.5-72B-Instruct | API (Dashscope) | 128K | No (domestic) |
| Qwen-2.5-7B-Instruct | Local (gpu1) | 128K | No |
| Llama-3.1-8B-Instruct | Local (gpu1) | 128K | No |

**Protocol per model × representation:**
1. Feed prompt with full dictionary + constraints
2. Generate N=200 sequences (4 batches × 50, with temperature=1.0)
3. Parse output → validate → record

**Metrics:**
| Metric | How |
|---|---|
| Format validity | % outputs that parse as token-ID sequences |
| ID validity | % where every token ID exists in library |
| Length compliance | % with correct ring size (6/8/10) |
| Constraint satisfaction | % passing all physicochemical checks |
| **Overall validity** | % passing all above |
| Chemical diversity | Mean pairwise Tanimoto distance (Morgan FP of reconstructed SMILES) |
| Uniqueness | % unique sequences out of valid ones |

**Output:** `experiments/output/exp1_validity/` — one CSV per model × representation

**Estimated cost:**
- API models: ~200 calls × 3 reps × 3 models = 1800 calls, ~$30-50
- Local models: free, ~2h GPU time

### W2.2 — Experiment 2: EDIT/FROZEN Controllability (Wed–Thu)

**Setup:**
- Select 10 diverse parent sequences from the random baseline
- For each: mark 1-2 positions as EDIT, rest as FROZEN
- Ask each model to generate 20 variants per parent (200 total per model)

**Representations:**
- ResToken: EDIT/FROZEN constraint in prompt with token IDs
- SMILES: "Modify only position 3 and 5" with SMILES notation (natural language instruction)
- HELM: analogous HELM-based instruction

**Metrics:**
| Metric | How |
|---|---|
| Edit compliance | % where only EDIT positions changed |
| Frozen compliance | % where FROZEN positions are unchanged |
| Diversity at EDIT positions | # unique substitutions per EDIT position |
| Property shift | Change in HBD/charge/bulk at EDIT vs parent |

**Output:** `experiments/output/exp2_controllability/`

### W2.3 — Experiment 3: Property-Constrained Generation (Thu–Fri)

**Setup:** 3 constraint profiles mimicking real design scenarios:

| Profile | Constraint | Rationale |
|---|---|---|
| Permeable | HBD≤1, charge=0, ≥2 NMe-backbone, bulk≥3 large | Oral bioavailability |
| Charged binder | charge=+2, HBD≥2, ≥1 aromatic | Electrostatic PPI |
| Rigid scaffold | total_rot≤8, ≥3 beta-backbone, no flexible | Conformational preorg |

**Protocol:** Each model × representation × profile: generate 100 sequences

**Metrics:**
- Constraint satisfaction rate (per-constraint and overall)
- Residue-class distribution (do models use the right types of blocks?)

**Output:** `experiments/output/exp3_property_constrained/`

### W2 Checkpoint
- [ ] All 5 models × 3 representations × 3 experiments complete
- [ ] Raw outputs + parsed results saved
- [ ] Summary statistics computed
- [ ] Any surprising findings logged in Lab Notebook

---

## Week 3 (Jun 9–13): Analysis + Figures + Start Writing

### W3.1 — Statistical Analysis (Mon–Tue)
- Aggregate results across experiments
- Compute confidence intervals (bootstrap, n=200 per condition)
- Run significance tests where needed (chi-square for validity rates)
- Identify the 2-3 key findings that tell the paper's story:
  1. ResToken validity >> SMILES/HELM (quantify the gap)
  2. EDIT/FROZEN compliance only works with structured representation
  3. Property constraints are enforceable at the token level

### W3.2 — Figure Preparation (Tue–Wed)
**Figure 1 (main):**
- **(a)** Tokenization scheme: pick 3 representative NCAAs (one alpha, one beta, one gamma) → show decomposition into [backbone | N-mod | side chain] → semantic token. Use RDKit 2D rendering for the chemical structures.
- **(b)** Representation comparison: same 6-mer peptide in SMILES / HELM / ResToken. Typeset as a comparison panel.
- **(c)** Benchmark bar chart: grouped bars for validity rate across models × representations. Use seaborn, ACS-style formatting (Arial 8pt, 300 DPI).

**Figure 2 (if space permits):**
- **(a)** EDIT/FROZEN heatmap: compliance rate per model × representation
- **(b)** Property distribution scatter: generated peptides in (HBD, rot_bonds) space, colored by representation

**Table 1:** Library statistics — backbone × class × chirality breakdown

All figures: `experiments/script/plot_figures.py` → `manuscript/figures/`

### W3.3 — Draft Methods + Results (Wed–Fri)
Following manuscript-drafting-submission skill: Methods first, then Results.

**Methods (~800 words):**
- 2.1 Building Block Library — curation source, property encoding, dual-dictionary
- 2.2 Constraint Validator — algorithm, configurable thresholds
- 2.3 Benchmark Protocol — models, representations, metrics, statistical tests

**Results (~800 words):**
- 3.1 Generation validity across representations
- 3.2 Controllability with EDIT/FROZEN
- 3.3 Property-constrained design
- Each subsection: one paragraph describing the figure/table, one paragraph interpreting

### W3 Checkpoint
- [ ] All figures at publication quality (300 DPI, ACS style)
- [ ] Methods + Results drafted
- [ ] Table 1 finalized

---

## Week 4 (Jun 16–20): Complete Draft + Review + Polish

### W4.1 — Draft Intro + Discussion + Abstract (Mon–Tue)

**Introduction (~500 words):**
- Para 1: Cyclic peptides + NCAAs as therapeutic modality
- Para 2: LLMs for molecular design — what works, what doesn't for NCAAs
- Para 3: Gap — no representation bridges NCAAs and LLMs
- Para 4: "Here we present ResToken..."

**Discussion (~300 words):**
- Interpret: why structured tokens work (closed-set prevents hallucination; semantic properties enable reasoning)
- Limitations: curated not exhaustive; binning loses precision; no property prediction
- Future: integration with scoring functions; community-contributed blocks; fine-tuning studies

**Abstract (~150 words):** Compress the whole story. Last sentence = availability.

**Title:** Finalize from the 3 candidates in outline.

### W4.2 — GitHub Repository Preparation (Tue–Wed)
- Clean `restoken/` into installable package (setup.py, __init__.py)
- Write README with install instructions + quick start
- Add example notebook
- License: MIT
- Push to GitHub (private initially, public at submission)

### W4.3 — Multi-Agent Review (Wed–Thu)
Per manuscript-drafting-submission skill:
1. Researcher agent: scientific accuracy, missing refs, method gaps
2. Assistant agent: clarity, grammar, logical flow
3. Merge feedback → prioritized revision list
- Reviews saved to `manuscript/drafts/review_*.md`

### W4.4 — Revise + Format (Thu–Fri)
- Address review comments
- Format for JCIM (ACS `achemso` LaTeX or Word template)
- Check: references (30-max), figure resolution, SI (if any)
- Prepare cover letter from skill template
- Final checklist from manuscript-drafting-submission skill

### W4 Checkpoint
- [ ] Complete manuscript draft
- [ ] GitHub repo ready (private)
- [ ] Multi-agent review complete + revisions incorporated
- [ ] Formatted per JCIM template
- [ ] Cover letter drafted
- [ ] Ready for PI review

---

## Risk Register

| Risk | Impact | Mitigation |
|---|---|---|
| Full 400-block dictionary doesn't fit in model context | Blocks Exp 1 | Pre-compute token counts; use class-indexed chunking if needed |
| SMILES baseline too hard (all models fail) | Weakens comparison | Include HELM as intermediate; document failure modes qualitatively |
| Low diversity in ResToken generation | Undermines "LLM adds value over random" | Tune temperature; report honestly if random ≈ LLM |
| API rate limits / costs exceed budget | Delays W2 | Prioritize 2 API models + 2 local; cut Qwen-72B if needed |
| RNF126 pipeline needs GPU during W2 | Delays local model benchmarks | Schedule local runs for off-peak; API models don't use GPU |
| Backbone compatibility rules unclear | Validator too permissive/strict | Start permissive, tighten based on chemical review |

---

## Compute Budget

| Resource | Estimated Use | Cost |
|---|---|---|
| GPT-4o API | ~600 calls | ~$15 |
| Claude Sonnet API | ~600 calls | ~$10 |
| Qwen-72B API (Dashscope) | ~600 calls | ~$5 |
| GPU (local models, W2) | ~4h on gpu1 (2× RTX 4090) | Free |
| Total | | ~$30 |

---

## Decision Points (flag for PI)

1. **After W1:** Is the 400-block library the right scope, or should we expand/filter?
2. **After W2:** If random baseline ≈ LLM generation quality, the story changes. Discuss.
3. **After W3:** Review figures + key results before writing Intro/Discussion.
