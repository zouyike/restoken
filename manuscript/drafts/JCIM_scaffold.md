# ResToken: A Residue-Semantic Token Library Enabling LLM-Based Design of Noncanonical Cyclic Peptides

**Journal:** JCIM Application Note (~3000 words)
**Status:** Scaffold — sections with expected content, figure/table specifications, placeholder metrics

---

## Abstract (~150 words)

Cyclic peptides incorporating noncanonical amino acids (NCAAs) are an expanding therapeutic modality, yet existing molecular representations — SMILES, SELFIES, HELM — are poorly suited for LLM-based NCAA design. SMILES strings expose structural details that cause hallucination; HELM lacks comprehensive NCAA coverage and semantic property encoding. We present ResToken, a curated library of 400 NCAA building blocks encoded as semantic tokens spanning three backbone types (alpha/beta/gamma), 20 functional classes, and L/D/achiral stereochemistry. Each token encodes eight designable properties while hiding structural details from the LLM. A constraint validator guarantees all generated sequences are chemically legal. In zero-shot benchmarks across five LLMs, ResToken achieves [XX]% overall validity versus [YY]% (SMILES) and [ZZ]% (HELM), with [XX]% constraint satisfaction under property-constrained design. The library, validator, and prompt templates are freely available as an open-source Python package.

---

## 1. Introduction (~500 words)

**Para 1 — The opportunity.** Cyclic peptides are emerging as a therapeutic modality bridging small molecules and biologics. NCAAs expand the accessible chemical space, enabling improved membrane permeability, protease resistance, and target affinity. [Cite: CycPeptMP, clinical cyclic peptides, mRNA-display libraries]

**Para 2 — LLMs for molecular design.** Recent LLM-based approaches have shown promise for peptide design (PeptideCLM, PepThink-R1, Evo-R), but these operate primarily on canonical amino acids. NCAA design with LLMs remains an open challenge because of the representation bottleneck.

**Para 3 — Why existing representations fail.**
- SMILES/SELFIES: LLMs hallucinate invalid NCAA structures; combinatorial explosion of valid SMILES makes closed-set enforcement impossible; no semantic property access at the token level
- HELM: designed for linear peptides; limited NCAA monomer coverage; no mechanism for property-aware reasoning
- One-hot/embedding: not interpretable; cannot express user constraints in natural language

**Para 4 — Our contribution.** "Here we present ResToken, a residue-semantic tokenization scheme that decouples what the LLM sees (semantic tokens with designable properties) from what the backend stores (full chemical structures). We benchmark ResToken against SMILES and HELM across five LLMs on three design tasks: unconstrained generation, controllable editing, and property-constrained design."

---

## 2. The ResToken Representation (~800 words)

### 2.1 Building Block Library

400 blocks curated from [source]. Coverage:

> **→ TABLE 1 goes here** (see below)

Each token encodes 8 semantic properties visible to the LLM:
1. `class` — functional analog (G/A/V/L/I/F/P/S/T/K/D/E/C/H/Y/N/Q/R/W/M/ncaa)
2. `chirality` — L / D / achiral
3. `charge_label` — neutral / positive / negative
4. `mc_type` — alpha / beta / gamma backbone
5. `mc_nmod` — N-modification (none / N-cyclic / N-methyl)
6. `sc_bulk` — side chain bulk (small / medium / large)
7. `polarity_bin` — low / medium / high
8. `flex_bin` — low / medium / high

Design space: 400^6 = 4.1 × 10^15 possible 6-mers.

### 2.2 Dual-Dictionary Architecture

> **→ FIGURE 1a goes here** (tokenization scheme)

- **LLM dictionary**: semantic properties only — what the LLM reads during generation. ~11K tokens for the full 400-block table. Fits within any modern LLM context window (128K+).
- **Backend dictionary**: full SMILES, SELFIES, InChIKey, all computed properties — used for structure reconstruction, scoring, synthesis planning.
- Design principle: closed-set token IDs prevent hallucination; semantic properties enable property-aware reasoning without structural detail leakage.

### 2.3 Constraint Validator

Seven configurable hard constraints checked before downstream evaluation:

| Check | Description |
|---|---|
| ID existence | All token IDs must exist in the 400-block library |
| Sequence length | Must be 6, 8, or 10 residues |
| Net charge | Sum of per-residue charges = target |
| HBD limit | Total hydrogen bond donors ≤ threshold |
| HBA limit | Total hydrogen bond acceptors ≤ threshold |
| Rotatable bonds | Total ≤ threshold |
| Backbone compatibility | No illegal backbone transitions (configurable strict/permissive) |

Rejection rate serves as a diagnostic: high rejection indicates poor prompt design or model confusion.

---

## 3. Benchmarks (~800 words)

### 3.1 Experimental Setup

> **→ FIGURE 1b goes here** (representation comparison panel)

**Models:**
| Model | Type | Context | Notes |
|---|---|---|---|
| GPT-4o | API | 128K | Via proxy |
| Claude Sonnet 4.6 | API | 200K | Via proxy |
| Qwen-2.5-72B-Instruct | API | 128K | Dashscope, domestic |
| Qwen-2.5-7B-Instruct | Local | 128K | gpu1, 1× RTX 4090 |
| Llama-3.1-8B-Instruct | Local | 128K | gpu1, 1× RTX 4090 |

**Representations:** ResToken, SMILES, HELM, Random baseline (rejection sampling)

**Task:** Generate N=200 sequences per model × representation, length 6, temperature 1.0

### 3.2 Experiment 1: Zero-Shot Generation Validity

> **→ FIGURE 2 goes here** (main benchmark bar chart)

**Metrics and expected results:**

| Metric | ResToken (expected) | SMILES (expected) | HELM (expected) | Random |
|---|---|---|---|---|
| Format validity | ~95-100% | ~60-80% | ~70-85% | 100% |
| ID / monomer validity | ~95-100% | ~20-40% | ~40-60% | 100% |
| Length compliance | ~90-95% | ~50-70% | ~60-80% | 100% |
| **Overall validity** | **~85-95%** | **~15-30%** | **~30-50%** | **100%** |
| Uniqueness (of valid) | ~99% | ~95% | ~97% | 100% |
| Chemical diversity (Tanimoto dist.) | ~0.70-0.80 | ? | ? | 0.76 |

**Key finding 1:** ResToken dramatically outperforms SMILES/HELM in generation validity because the closed-set token library eliminates hallucination of nonexistent building blocks.

**Key finding 2:** The main ResToken failure mode is format errors (wrong number of tokens, extra text) — NOT chemical errors. This is fixable with better prompting or lightweight parsing.

### 3.3 Experiment 2: EDIT/FROZEN Controllability

> **→ FIGURE 3 goes here** (controllability heatmap)

**Setup:** 10 parent sequences, 1-2 EDIT positions + rest FROZEN, 20 variants per parent per model.

**Expected results:**

| Metric | ResToken | SMILES | HELM |
|---|---|---|---|
| Edit compliance (only EDIT changed) | >90% | ~40-60% | ~50-70% |
| Frozen compliance (FROZEN unchanged) | >95% | ~30-50% | ~50-70% |
| Diversity at EDIT positions | High | Low (often invalid) | Medium |

**Key finding 3:** Structured token representation enables precise positional control. SMILES-based editing suffers from "edit bleed" — LLMs change more positions than instructed because SMILES modifications are non-local.

### 3.4 Experiment 3: Property-Constrained Design

> **→ FIGURE 4 goes here** (property distribution scatter/violin)

**Three constraint profiles:**

| Profile | Key constraints | Random baseline acceptance |
|---|---|---|
| Permeable | charge=0, HBD≤1, ≥2 NMe/NCY, ≥3 large bulk | **2.07%** |
| Charged binder | charge=+2, HBD≥2, ≥1 aromatic | **1.82%** |
| Rigid scaffold | rot≤18, ≥3 beta backbone | **0.81%** |

**Expected results:**

| Metric | ResToken | SMILES | HELM | Random |
|---|---|---|---|---|
| Constraint satisfaction (permeable) | ~50-80% | ~5-15% | ~10-20% | 2.07% |
| Constraint satisfaction (charged) | ~40-70% | ~5-10% | ~10-20% | 1.82% |
| Constraint satisfaction (rigid) | ~30-60% | ~2-5% | ~5-10% | 0.81% |

**Key finding 4:** ResToken enables LLMs to reason about constraints at the token level — each block's properties are visible in the prompt, so models can select blocks that satisfy charge/HBD/bulk requirements. SMILES provides no such property access.

**Key finding 5:** Even under tight constraints (0.81% random acceptance), LLMs with ResToken achieve [XX]× higher constraint satisfaction than random sampling, demonstrating that LLMs provide informed sampling beyond random exploration.

---

## 4. Discussion / Limitations (~300 words)

**Why it works:** Closed-set tokens prevent hallucination (the #1 failure mode for SMILES). Semantic properties enable property-level reasoning without structural detail. The dual-dictionary architecture separates LLM-facing simplicity from backend chemical completeness.

**Limitations:**
1. Curated, not exhaustive — 400 blocks cover the design-relevant space but new NCAAs require manual addition
2. Gamma backbone thin (49 blocks) — no gamma variants for 10/20 classes (C, D, E, H, K, M, N, Q, R, W, Y)
3. N-methylation nearly absent (2 NME blocks, alpha only) — a key permeability lever with minimal library coverage
4. Semantic binning (small/med/large) loses quantitative resolution — deliberate tradeoff for LLM usability
5. No property prediction — ResToken is a representation, not a scoring function
6. Zero-shot performance depends on prompt engineering; fine-tuned models expected to improve

**Future work** (brief): integration with permeability/activity predictors; expanding to stapled peptides; community-contributed NCAA blocks.

---

## 5. Software and Availability (~200 words)

Python package `restoken`:
- 400-block library (CSV + JSON, backend + LLM-safe versions)
- `SequenceValidator`: 7 configurable checks
- `SMILESReconstructor`: token sequence → full molecular SMILES (400/400 RDKit-validated)
- Prompt templates for major LLMs (ResToken / SMILES / HELM)
- Random baseline generator with 4 constraint profiles
- GitHub: [URL] | License: MIT

---

## Figures & Tables Specification

### TABLE 1: Building Block Library Statistics

| Property | Alpha | Beta | Gamma | Total |
|---|---|---|---|---|
| **Count** | 250 | 101 | 49 | **400** |
| **Chirality** | | | | |
| — L | 219 | 44 | 29 | 292 |
| — D | 19 | 51 | 16 | 86 |
| — Achiral | 12 | 6 | 4 | 22 |
| **Charge** | | | | |
| — Neutral | 215 | 85 | 47 | 347 |
| — Positive | 21 | 8 | 1 | 30 |
| — Negative | 12 | 7 | 1 | 20 |
| — Zwitterionic | 2 | 1 | 0 | 3 |
| **N-modification** | | | | |
| — None | 210 | 97 | 46 | 353 |
| — N-cyclic (NCY) | 38 | 4 | 3 | 45 |
| — N-methyl (NME) | 2 | 0 | 0 | 2 |
| **Bulk** | | | | |
| — Small | 44 | 42 | 30 | 116 |
| — Medium | 33 | 17 | 4 | 54 |
| — Large | 173 | 42 | 15 | 230 |
| **Functional classes** | 18/21 | 20/21 | 10/21 | 21 |

*400 NCAA building blocks spanning three backbone types. Beta backbone has the most balanced chirality distribution (L:D ≈ 1:1). Alpha dominates in N-cyclic variants (proline analogs). Gamma coverage is limited to 10 functional classes.*

---

### FIGURE 1: The ResToken Tokenization Scheme (full-width, main figure)

**Panel (a) — Tokenization decomposition:**
- Pick 3 structurally diverse NCAAs: one alpha-L-large (e.g., A01, isoleucine analog), one beta-D-small (e.g., N03), one gamma-NCY (e.g., S24 or similar)
- For each: show 2D chemical structure (RDKit) → arrow → decomposition into [backbone type | N-mod | side chain] → arrow → semantic token card showing all 8 properties
- Visual: chemical structure on left, property card on right, arrow showing abstraction

**Panel (b) — Representation comparison:**
- Same 6-mer peptide (e.g., A01-K03-N12-S05-E02-a07) shown three ways:

```
SMILES: CCC(C)(C)[C@H](NC(=O)[C@H](CCC(=O)[O-])CNC(=O)[C@H](CSC)CNC(=O)...
        (long, unreadable, error-prone)

HELM:   PEPTIDE1{[A01].[K03].[N12].[S05].[E02].[a07]}$PEPTIDE1,PEPTIDE1,1:R1-6:R2$$$V2.0
        (structured but no property access)

ResToken: A01-K03-N12-S05-E02-a07
          I/L/neu  E/L/neg  M/L/neu  I/L/neu  M/L/neu  P/D/neu
          alpha    beta     beta     gamma    alpha    alpha
          (short, every property readable per position)
```

- Highlight: token-level property annotation is visible in ResToken, invisible in SMILES/HELM

**Production notes:** RDKit 2D rendering for structures. Illustrator/Inkscape for layout. ACS style: Arial 8pt, 300 DPI.

---

### FIGURE 2: Zero-Shot Generation Validity Benchmark (main result)

**Type:** Grouped bar chart, 3 groups (ResToken / SMILES / HELM), 5 bars per group (one per model), y-axis = overall validity rate (%).

**Additional elements:**
- Horizontal dashed line at random baseline level (100% for unconstrained — but this is trivial; more useful: show constraint satisfaction rates)
- Error bars from 4 batches × 50 sequences
- Annotation: exact percentages on top of each bar

**Alt layout (if cleaner):** Heatmap — rows = models, columns = representations, cell color = validity rate. Simpler to read with 5 models × 3 reps = 15 cells.

**Script:** `experiments/script/plot_figures.py` → `manuscript/figures/fig2_validity_benchmark.pdf`

---

### FIGURE 3: EDIT/FROZEN Controllability (half-width)

**Type:** Grouped bar chart or heatmap.

**Data structure:**
- Rows: models (3-5)
- Column groups: Edit compliance, Frozen compliance
- Column sub-groups: ResToken, SMILES, HELM

**Expected visual story:** ResToken bars consistently high (>90%), SMILES bars low and variable (~30-60%), HELM intermediate.

**Script:** `experiments/script/plot_figures.py` → `manuscript/figures/fig3_controllability.pdf`

---

### FIGURE 4: Property-Constrained Generation (half-width)

**Type:** 2D scatter or violin plot.

**Option A — Scatter:**
- x-axis: total HBD, y-axis: total rotatable bonds
- Points colored by source: ResToken-generated (blue), SMILES-generated (red), random baseline (gray)
- Constraint boundary box drawn for "permeable" profile (HBD≤1)
- Expected: ResToken points cluster inside the constraint box; SMILES points scatter outside

**Option B — Violin/box:**
- x-axis: constraint profile (permeable / charged / rigid)
- y-axis: constraint satisfaction rate (%)
- Grouped by representation
- Shows the distribution across models

**Data:** From Experiment 3 (W2.3), using our baseline numbers as anchoring:
- Random acceptance: 2.07% (permeable), 1.82% (charged), 0.81% (rigid)
- LLM performance expected to be 10-50× higher with ResToken

**Script:** `experiments/script/plot_figures.py` → `manuscript/figures/fig4_property_constrained.pdf`

---

### FIGURE 5: Library Coverage Treemap or Sunburst (supplementary or Table 1 companion)

**Type:** Treemap or sunburst chart showing the hierarchical structure of the 400-block library.

**Hierarchy:** Backbone type → Functional class → Chirality
- Area proportional to block count
- Color: backbone type (alpha=blue, beta=green, gamma=orange)

**Purpose:** Visual complement to Table 1 — shows at a glance where the library is dense (alpha-I, alpha-A, alpha-P) and where it's sparse (gamma-anything, all-R, all-Q).

**Script:** `experiments/script/plot_figures.py` → `manuscript/figures/fig5_library_treemap.pdf`

---

### TABLE 2: Per-Model Benchmark Summary (supplementary)

| Model | Rep | Format Valid | ID Valid | Length OK | Constraint Sat | Overall Valid | Unique | Diversity |
|---|---|---|---|---|---|---|---|---|
| GPT-4o | ResToken | | | | | | | |
| GPT-4o | SMILES | | | | | | | |
| GPT-4o | HELM | | | | | | | |
| Claude Sonnet | ResToken | | | | | | | |
| ... | ... | | | | | | | |
| Random | — | 100% | 100% | 100% | varies | varies | 100% | 0.76 |

*Cells to be filled from Week 2 experiments.*

---

## Key Messages (for writing each section)

1. **The representation is the bottleneck, not the model.** Even GPT-4o fails at NCAA design with SMILES — the problem isn't model capability, it's that SMILES wasn't designed for LLM consumption.

2. **Closed-set tokens eliminate hallucination.** The single biggest win: LLMs cannot invent nonexistent NCAAs when the vocabulary is a fixed 400-block library.

3. **Semantic properties enable constraint reasoning.** Unlike SMILES where charge/bulk/flexibility are encrypted in the string, ResToken makes these properties visible and actionable.

4. **The validator is as important as the tokens.** It provides a formal guarantee that every accepted sequence is chemically legal — a safety net that SMILES/HELM approaches lack.

5. **Constrained design is genuinely hard.** Random acceptance rates of 0.8-2.1% prove that property-constrained NCAA peptide design is a non-trivial optimization problem. LLMs with ResToken provide informed sampling that significantly outperforms random.

6. **Library gaps are honest limitations, not failures.** Gamma (49 blocks), NME (2 blocks), and charged gamma (2 blocks) are sparse — reflecting the current state of synthesizable NCAAs, not a curation failure.
