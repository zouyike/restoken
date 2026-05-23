# JCIM Application Note — Outline

## Working Title (3 candidates)

1. **ResToken: A Residue-Semantic Token Library Enabling LLM-Based Design of Noncanonical Cyclic Peptides**
2. **A Semantic Building Block Representation for LLM-Guided Noncanonical Amino Acid Cyclic Peptide Design**
3. **ResToken: Bridging Noncanonical Amino Acids and Language Models through Semantic Residue Tokenization**

Recommendation: #1 — concise, names the tool, states what it does.

---

## Framing (the "so what")

**Problem:** LLMs are increasingly used for peptide design, but existing molecular representations (SMILES, SELFIES, HELM) fail for NCAAs in cyclic peptides:
- SMILES/SELFIES: LLMs hallucinate invalid chemistry, can't enforce physicochemical constraints at the token level
- HELM: designed for linear peptides; limited NCAA coverage; no semantic property encoding
- No existing representation lets an LLM reason about *why* a substitution changes permeability/rigidity/charge

**Solution:** A curated library of 400 NCAA building blocks encoded as semantic tokens — each token carries designable properties (chirality, charge, bulk, flexibility, polarity, H-bond capacity) but hides structural details (SMILES) that cause hallucination. Combined with a constraint validator that guarantees all LLM outputs are chemically legal.

**Contribution type:** Tool/resource (open-source library + validator + benchmark), not a discovery paper.

---

## Structure (JCIM App Note format, ~3000 words)

### Abstract (150 words)
- Gap: LLMs struggle with NCAA peptide design because representations expose too much or too little chemical information
- Tool: ResToken — 400 semantic tokens spanning alpha/beta/gamma backbones, 19 functional classes, L/D/achiral stereochemistry
- Key result: LLM generation validity X% vs Y% (SMILES) and Z% (HELM); constraint satisfaction rate; chemical diversity metric
- Availability: GitHub + pip install

### 1. Introduction (~500 words)
- Cyclic peptides as therapeutics (cell permeability challenge, NCAAs as solution — cite CycPeptMP, Evo-R, PepThink-R1)
- LLMs for molecular design (PeptideCLM, ChemGPT, etc.) — success on canonical AAs, failure on NCAAs
- Why existing representations fail for NCAA design:
  - SMILES: combinatorial explosion, hallucination, no property reasoning
  - HELM: limited NCAA coverage, no semantic properties
  - One-hot/embedding: not interpretable, can't constrain
- Our approach: decouple **what the LLM sees** (semantic tokens) from **what the backend stores** (full chemistry)

### 2. The ResToken Representation (~800 words)

#### 2.1 Building Block Library
- 400 blocks curated from [source — need to clarify how the 400 were selected]
- Coverage: alpha (250), beta (101), gamma (49) backbone types
- 19 functional classes: I (88), F (46), A (46), T (42), P (38), etc.
- Chirality: L (292), D (86), achiral (22)
- Each token encodes 8 semantic properties:
  - `class` (functional analog: G, A, V, L, I, F, P, S, T, K, D, C, H, Y, E, N, Q, R, ncaa)
  - `chirality` (L / D / A)
  - `charge_label` (neu / pos / neg)
  - `mc_type` (alpha / beta / gamma)
  - `mc_nmod` (NO / NCY / NME / NXX) — backbone N-modification
  - `sc_bulk` (small / med / large)
  - `polarity_bin` (low / med / high)
  - `flex_bin` (low / med / high)

**Table 1:** Library statistics — breakdown by backbone type × class × chirality. Show the design space coverage.

#### 2.2 Dual-Dictionary Architecture
- **Backend dictionary** (`bb_dict_backend.json`): full SMILES, SELFIES, InChIKey, all computed properties — used for structure reconstruction, scoring, and synthesis planning
- **LLM dictionary** (`bb_dict_llm.json`): semantic properties only — what the LLM reads during generation
- Design principle: LLM gets enough to make design decisions but not enough to hallucinate chemistry
- Bidirectional mapping: LLM output (token ID sequence) → backend → full molecular structure → RDKit/scoring

#### 2.3 Constraint Validator
- Hard constraints checked before any downstream evaluation:
  - All token IDs exist in the library
  - Sequence length within allowed range
  - Net charge = target
  - Total HBD ≤ threshold (permeability)
  - Total HBA ≤ threshold
  - Total rotatable bonds ≤ threshold
  - Backbone compatibility (no illegal beta-gamma-alpha transitions, etc.)
- Reject rate as a diagnostic: high rejection = poor prompt design or model confusion

### 3. Benchmark: LLM Generation with ResToken vs Baselines (~800 words)

#### 3.1 Experimental Setup
- Models: GPT-4o, Claude Sonnet, Qwen-2.5-72B (API-based, zero-shot)
- Also: fine-tuned Qwen-2.5-0.5B and Llama-3.1-8B (to show transferability)
- Task: generate N=1000 cyclic peptides of length 6-10, satisfying given constraints
- Representations compared:
  1. **ResToken** (this work): semantic tokens + constraint prompt
  2. **SMILES**: NCAA SMILES + cyclization instruction
  3. **HELM**: HELM notation with NCAA monomer codes
  4. **Random baseline**: uniform random sampling from the 400-block library with constraint filtering

#### 3.2 Metrics
- **Validity rate**: % of outputs that parse as legal sequences (all IDs exist, correct format)
- **Constraint satisfaction**: % that pass all physicochemical constraints (charge, HBD, HBA, rot bonds)
- **Chemical diversity**: Tanimoto distance of generated set (via ChemBERTa embeddings or Morgan FP of reconstructed SMILES)
- **Property distribution**: do generated peptides cover the target property space, or collapse to a narrow region?

#### 3.3 Results
- **Figure 1:** (a) Tokenization scheme overview — NCAA structure → [backbone | side chain] → semantic token with properties. (b) Comparison: same peptide in SMILES / HELM / ResToken. (c) Bar chart: validity × constraint satisfaction across representations and models.

- Expected findings:
  - ResToken validity ~90-95% (tokens are closed-set; main failure mode is wrong length or format)
  - SMILES validity ~20-40% (LLMs invent nonexistent NCAAs, break ring closures)
  - HELM validity ~40-60% (better than SMILES but limited NCAA coverage)
  - Constraint satisfaction: ResToken >> others (constraints are checkable at generation time)
  - Diversity: ResToken ≥ random baseline (LLM provides informed sampling, not just random)

#### 3.4 Controllability Demonstration
- **EDIT/FROZEN constraint**: give LLM a parent sequence, allow mutation at only 1-2 positions
- Show: ResToken + EDIT achieves >95% compliance; SMILES-based editing has <50% compliance (LLMs change more than instructed)
- This is the key practical advantage: controlled lead optimization, not just de novo generation

### 4. Case Study: Permeability-Aware Design (~400 words)
- Use ResToken tokens' semantic properties to guide design toward permeable peptides:
  - Prompt: "Generate a 6-mer cyclic peptide with ≤1 HBD, net charge 0, at least 2 NMe-backbone residues"
  - Show: ResToken constrains LLM to output only blocks satisfying these criteria
  - Compare: SMILES prompt with equivalent natural-language instructions → much lower compliance
- Property analysis of generated set: HBD distribution, charge, backbone type usage
- NOT claiming these are actually permeable (no PAMPA data) — showing the representation enables property-constrained generation

### 5. Software and Availability (~200 words)
- Python package: `restoken` (pip-installable)
- Contents:
  - 400-block library (CSV + JSON, backend + LLM-safe versions)
  - Validator module
  - Prompt templates for major LLMs
  - SMILES reconstruction utility (token sequence → full molecular structure)
  - Example notebooks
- GitHub: [URL]
- License: MIT

### 6. Discussion / Limitations (~300 words)
- Limitations:
  - Current library is curated, not exhaustive — new NCAAs require manual addition
  - Semantic binning (small/med/large bulk) loses quantitative resolution — deliberate tradeoff for LLM usability
  - No property prediction — ResToken is a representation, not a scoring function
  - Zero-shot LLM performance depends on prompt engineering; fine-tuned models expected to improve
- Future work (brief — don't oversell):
  - Integration with permeability/activity predictors (the active learning loop — future paper)
  - Expanding to stapled peptides, peptide-drug conjugates
  - Community-contributed NCAA blocks

---

## Figures Plan

### Figure 1 (Main — full width)
**(a)** Tokenization scheme: NCAA chemical structure → decomposition into [backbone type + N-mod + side chain properties] → semantic token with 8 property fields. Show 3 example NCAAs spanning alpha/beta/gamma.

**(b)** Representation comparison: same 6-mer NCAA cyclic peptide shown as SMILES string (long, unreadable), HELM notation (medium), and ResToken sequence (short, semantic: `A09-i01-T15-a02-K03-D01`). Highlight: token-level property readout visible in ResToken, invisible in SMILES/HELM.

**(c)** Benchmark results: grouped bar chart — validity rate and constraint satisfaction rate for ResToken / SMILES / HELM / Random, across 2-3 LLMs.

### Figure 2 (Optional — if space permits)
**(a)** EDIT/FROZEN controllability: parent sequence with marked EDIT positions → LLM output compliance rate across representations.

**(b)** Property distribution of generated peptides: 2D scatter (e.g., total HBD vs flexibility) showing ResToken-generated set vs random baseline. Should show that LLM-guided generation is more focused but still diverse.

### Table 1
Building block library statistics: rows = backbone type (alpha/beta/gamma), columns = key properties (count, class distribution, chirality, bulk, polarity, N-modification). Shows design space coverage at a glance.

---

## Estimated Timeline

| Task | Time | Notes |
|---|---|---|
| Finalize library + validator code | 1 week | v11 exists; package it |
| Benchmark: API-based LLMs (GPT-4o, Claude, Qwen) | 1 week | ~$50-100 API cost |
| Benchmark: fine-tuned models (optional) | 1-2 weeks | Use existing CycSeqSeek infra |
| Figure preparation | 1 week | RDKit for structures, matplotlib for data |
| Writing | 1-2 weeks | ~3000 words |
| Internal review (multi-agent) | 3 days | |
| **Total** | **5-7 weeks** | No wet-lab dependency |

---

## Key Positioning Decisions

1. **Do NOT claim property prediction.** The tokens encode properties but don't predict outcomes. This is a representation paper, not a model paper.

2. **Do NOT mention the active learning loop.** That's the big paper. This paper is: "here's the building block set + representation that makes NCAA design with LLMs possible."

3. **Competitor framing:** Not "we beat X." Instead: "Existing representations weren't designed for LLM-based NCAA design. We fill that gap." Cite PeptideCLM, HELM, Evo-R, PepThink-R1 as complementary, not competing.

4. **The validator is as important as the tokens.** JCIM reviewers will want to know: "how do you prevent garbage?" The validator + rejection rate analysis is the answer.

5. **Open-source is mandatory for JCIM App Note.** The library + code must be public before submission.
