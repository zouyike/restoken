# Methods and Results — Draft v1

**Status:** First draft with real W2 data (minus Sonnet). Fill Sonnet numbers after midnight run completes.

---

## 2. Methods

### 2.1 Building Block Library

The ResToken library comprises 400 NCAA building blocks curated from commercially available and synthetically accessible noncanonical amino acids. Blocks span three backbone types — alpha (250), beta (101), and gamma (49) — with L (292), D (86), and achiral (22) stereochemistry across 21 functional classes analogous to canonical amino acids (Table 1). Each block is assigned a compact alphanumeric identifier (e.g., A01, N15, s01) encoding backbone type implicitly through a prefix convention: uppercase letters denote alpha backbones, specific letter ranges denote beta and gamma types, and case encodes chirality (uppercase = L, lowercase = D).

Each token carries eight semantic properties visible to the LLM during generation: functional class, chirality, charge label, backbone type, N-modification status, side-chain bulk, polarity, and flexibility. These properties are sufficient for property-aware design decisions while abstracting away structural details (bond connectivity, stereochemistry notation, ring closures) that cause LLM hallucination in SMILES-based approaches. The full design space of 6-residue peptides is 400^6 ≈ 4.1 × 10^15.

### 2.2 Dual-Dictionary Architecture

ResToken employs a dual-dictionary architecture that decouples the LLM-facing representation from the chemical backend (Figure 1a). The **LLM dictionary** contains only semantic properties — the 400-block table fits within ~11,000 tokens, well within any modern LLM context window. The **backend dictionary** stores full AA-SMILES, side-chain SMILES, SELFIES, InChIKey, and all computed physicochemical descriptors, enabling reconstruction of complete molecular structures from any valid token sequence.

This separation enforces a closed-set vocabulary: the LLM can only output identifiers from the 400-block library, preventing hallucination of nonexistent monomers — the dominant failure mode when LLMs generate NCAA structures as SMILES strings. The backend enables bidirectional mapping: token sequences can be reconstructed to full SMILES (validated via RDKit for all 400 blocks), while existing NCAA structures can be decomposed into their nearest ResToken equivalents.

### 2.3 Sequence Validator

A constraint validator checks all LLM-generated sequences against seven configurable hard constraints before downstream evaluation: (1) token ID existence in the library, (2) sequence length compliance, (3) net charge matching the target, (4) total hydrogen bond donors ≤ threshold, (5) total hydrogen bond acceptors ≤ threshold, (6) total rotatable bonds ≤ threshold, and (7) backbone transition compatibility. The rejection rate at each checkpoint serves as a diagnostic for prompt design quality and model instruction-following capability.

### 2.4 Benchmark Design

We evaluated ResToken against two baseline representations — SMILES and HELM — across five LLMs on three tasks of increasing difficulty.

**Models.** Three API-based frontier models: Gemini 2.5 Pro, Gemini 2.5 Flash, and GPT-4o. Two locally deployed open-weight models: Qwen 3.5 9B and Gemma 3 12B (4-bit quantized), run on RTX 4090 GPUs. All models received the same system prompt containing the full LLM dictionary and task instructions. TxGemma 9B was also tested but produced zero valid output across all representations and is excluded from analysis.

**Task 1 — Zero-shot generation validity (Exp1).** Each model generated N=200 sequences of length 6 in each representation (ResToken, SMILES, HELM) with temperature 1.0 and no property constraints. Metrics: overall validity rate (format + ID + length), uniqueness among valid sequences.

**Task 2 — Positional controllability (Exp2).** Ten parent ResToken sequences were provided with 2 positions marked EDIT (mutable) and 4 marked FROZEN (immutable). Each model generated 20 variants per parent (200 total). Metrics: frozen compliance (% of outputs with all FROZEN positions unchanged), edit compliance (% with at least one EDIT position changed), and block ID validity.

**Task 3 — Property-constrained generation (Exp3).** Each model generated 100 ResToken sequences under three constraint profiles: (1) *permeable* — net charge 0, HBD ≤ 1, ≥ 2 NMe/NCY backbone residues, ≥ 3 large-bulk side chains (random baseline acceptance: 2.07%); (2) *charged binder* — net charge +2, HBD ≥ 2, ≥ 1 aromatic residue (random baseline: 1.82%); (3) *rigid scaffold* — total rotatable bonds ≤ 18, ≥ 3 beta-backbone residues (random baseline: 0.81%). Random baselines were computed by rejection sampling from 100,000 uniform-random 6-mers.

**Statistical analysis.** Validity rates are reported with bootstrap 95% confidence intervals (10,000 resamples). Representation comparisons use chi-square tests with Bonferroni correction (α = 0.05, 10 comparisons). Enrichment is reported as fold-improvement over random baseline acceptance.

---

## 3. Results

### 3.1 ResToken Achieves Superior Generation Validity

Across five LLMs, ResToken achieved a mean validity rate of 90.4 ± 9.0% (range: 77.3–98.5%), significantly outperforming SMILES (69.0 ± 27.4%) and HELM (84.5 ± 11.4%) (Figure 2a, Table 2). The best-performing model, Gemini 2.5 Pro, achieved 98.5% validity with ResToken (95% CI: 96.5–100%) versus 72.1% with SMILES (95% CI: 65.5–78.2%; χ² = 53.4, p < 10⁻¹², Bonferroni-corrected). Notably, Pro achieved 100% validity with HELM, but this came at the cost of reduced diversity compared to ResToken.

The validity advantage of ResToken was model-dependent. For frontier API models (Gemini Pro/Flash, GPT-4o), ResToken validity exceeded 92% uniformly. Open-weight local models showed greater variability: Qwen 3.5 9B achieved 77.3% validity, with failures predominantly due to format errors (wrong number of tokens) rather than invalid block IDs. Gemma 3 12B achieved 85.6% but exhibited hallucination of sequential block IDs beyond the library range (e.g., A190–A216), suggesting the model attempted to extrapolate the ID numbering scheme rather than selecting from the provided library.

SMILES validity was paradoxically higher for local models (92–98%) than for frontier models (37–72%), likely because local models produce shorter, simpler SMILES strings that happen to be valid rather than attempting complex NCAA structures. However, this high SMILES validity came with extremely low uniqueness (0.29–0.31), indicating repetitive generation of a narrow set of canonical-like structures.

### 3.2 Structured Tokens Enable Precise Positional Control

In the controllability benchmark (Exp2), all five models achieved 100% frozen-position compliance with ResToken, confirming that the discrete, position-indexed token format makes positional constraints trivial to enforce. Edit compliance — the rate at which models actually modified the designated EDIT positions — scaled with model capability: Gemini Pro and Flash both achieved 100%, GPT-4o 76%, Gemma 3 12B 61%, and Qwen 3.5 9B 19.5% (Figure 2b). The low edit compliance of local models reflects a conservative strategy: these models often returned the parent sequence unchanged rather than risk generating invalid substitutions.

Block ID validity remained high across all models (91–100%), confirming that ResToken's finite vocabulary supports reliable token-level editing even in open-weight 9–12B parameter models. This contrasts with SMILES-based editing, where modifying a single residue requires coordinating changes across multiple non-contiguous positions in the string — a task that causes "edit bleed" where LLMs inadvertently alter more positions than instructed.

### 3.3 Property Constraints Differentiate Model Reasoning Capability

Property-constrained generation (Exp3) revealed the starkest performance differences across models. Only Gemini 2.5 Pro reliably satisfied all three constraint profiles: 100% for permeable (48× enrichment over random), 80% for charged binder (44× enrichment), and 100% for rigid scaffold (123× enrichment). Gemini Flash achieved 86% on the permeable profile (42× enrichment) but failed entirely on charged binder and rigid scaffold (0%).

GPT-4o showed intermediate performance: 30% on permeable (15× enrichment), 0% on charged binder, and 70% on rigid scaffold (86× enrichment) — though the rigid scaffold result is based on only 20 valid sequences, limiting its statistical reliability. All local models achieved 0% constraint satisfaction across all profiles despite producing valid ResToken sequences (Figure 2c).

These results highlight a critical design principle: ResToken encodes the property information needed for constraint satisfaction directly in the token metadata, but exploiting this information requires strong instruction-following and in-context reasoning capabilities found only in frontier models. The mean enrichment over random baseline was 34.9× for permeable (3 models), 44.0× for charged binder (1 model), and 104.9× for rigid scaffold (2 models) — demonstrating that capable models with ResToken provide dramatically more efficient sampling than random exploration of the combinatorial space.

---

*Word count: Methods ~650, Results ~600. Total ~1250. Target was ~1600 (800+800). Can expand with Sonnet data and additional analysis.*

*TODO after Sonnet benchmark completes:*
- *Add Sonnet row to all tables*
- *Update aggregate statistics*
- *Revise Figure 2 to include 6th model*
- *Update chi-square tests with additional comparisons*
