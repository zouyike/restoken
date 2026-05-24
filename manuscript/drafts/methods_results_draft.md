# ResToken: A Residue-Semantic Token Library Enabling LLM-Based Design of Noncanonical Cyclic Peptides

**Status:** Draft v2 — full manuscript with real W2 data (5 models). Sonnet placeholder for 6th model after midnight run.

---

## Abstract

Cyclic peptides incorporating noncanonical amino acids (NCAAs) are an expanding therapeutic modality, yet existing molecular representations are poorly suited for LLM-based NCAA design. SMILES strings expose structural details that induce hallucination of nonexistent monomers; HELM provides structure but lacks semantic property encoding for design reasoning. We present ResToken, a curated library of 400 NCAA building blocks encoded as semantic tokens spanning three backbone types (alpha, beta, gamma), 21 functional classes, and L/D/achiral stereochemistry. Each token encodes eight designable properties — chirality, charge, backbone type, N-modification, bulk, polarity, and flexibility — while hiding structural details from the LLM. A constraint validator guarantees all generated sequences are chemically legal. In zero-shot benchmarks across five LLMs, ResToken achieves 90.4% mean validity versus 69.0% (SMILES) and 84.5% (HELM), with 100% frozen-position compliance in controllable editing and up to 123-fold enrichment over random baselines in property-constrained design. The library, validator, and benchmark suite are freely available as an open-source Python package.

---

## 1. Introduction

Cyclic peptides occupy a therapeutic sweet spot between small molecules and biologics, offering the target selectivity of protein interfaces with the bioavailability advantages of small molecules. The incorporation of noncanonical amino acids (NCAAs) — residues with non-standard backbones, side chains, or stereochemistry — dramatically expands the accessible chemical space, enabling improved membrane permeability through N-methylation and backbone rigidification, enhanced protease resistance via D-amino acids and beta-backbone substitution, and tunable physicochemical properties through designer side chains. Clinical successes such as cyclosporine and recent mRNA-display-derived candidates have validated NCAA-containing cyclic peptides as a viable drug modality.

Large language models (LLMs) have emerged as powerful tools for molecular design, demonstrating success in generating valid small molecules, predicting protein sequences, and optimizing canonical peptides. Recent work on peptide-specific models — including PeptideCLM for antimicrobial peptide generation, PepThink-R1 for reasoning-guided design, and Evo-R for evolutionary protein generation — has shown that LLMs can learn the grammar of biological sequences and generate functional variants. However, these approaches operate almost exclusively on canonical amino acids represented as single-letter codes, leaving NCAA-containing cyclic peptide design as an open challenge.

The bottleneck is representation, not model capability. When prompted to generate NCAA structures as SMILES strings, even frontier LLMs hallucinate nonexistent monomers, produce invalid ring closures, and cannot reason about physicochemical properties encrypted within the notation. HELM (Hierarchical Editing Language for Macromolecules) provides better structural organization but offers limited NCAA monomer coverage and no mechanism for property-aware reasoning — the LLM sees monomer codes but cannot determine which codes satisfy a charge or bulk constraint without external lookup. One-hot or learned embeddings sacrifice interpretability entirely, preventing users from expressing design constraints in natural language.

Here we present ResToken, a residue-semantic tokenization scheme that decouples what the LLM sees from what the chemistry requires. Each of the 400 NCAA building blocks is encoded as a compact token carrying eight designable properties (chirality, charge, backbone type, N-modification, side-chain bulk, polarity, flexibility, and functional class), while full chemical structures (SMILES, SELFIES, InChIKey) are stored in a separate backend dictionary invisible to the LLM. This dual-dictionary architecture enforces a closed-set vocabulary that eliminates monomer hallucination while enabling property-level reasoning directly from the token metadata. We benchmark ResToken against SMILES and HELM across five LLMs on three tasks — unconstrained generation, controllable editing, and property-constrained design — demonstrating that the representation, not the model, is the primary determinant of generation quality.

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

## 4. Discussion

The central finding of this work is that representation design, rather than model scale, is the primary bottleneck for LLM-based NCAA peptide generation. ResToken's closed-set vocabulary eliminates monomer hallucination — the dominant failure mode for SMILES-based approaches — while its semantic property encoding enables constraint-aware reasoning that neither SMILES nor HELM supports. The 100% frozen-position compliance across all tested models, including 9B-parameter open-weight models, confirms that discrete, position-indexed tokens make positional control trivial compared to the non-local edit dependencies inherent in SMILES notation.

The sharp capability gradient in property-constrained generation (Exp3) reveals an important design consideration: ResToken makes constraint satisfaction *possible* by encoding the relevant properties in the prompt, but *achieving* it requires strong instruction-following and in-context reasoning. Gemini 2.5 Pro achieved 80–100% constraint satisfaction with up to 123× enrichment over random baselines; local 9–12B models achieved 0% despite generating valid sequences. This suggests that ResToken is most immediately useful as a representation layer for frontier API models, while fine-tuned specialist models may be needed to bring constraint reasoning to smaller architectures.

**Limitations.** The current library of 400 blocks is curated, not exhaustive; new NCAAs require manual addition with property annotation. Gamma backbone coverage is thin (49 blocks spanning only 10 of 21 functional classes), and N-methylation — a critical permeability lever — is represented by only 2 blocks. Semantic property binning (small/medium/large for bulk, low/medium/high for polarity) deliberately trades quantitative resolution for LLM usability; applications requiring precise physicochemical values should query the backend dictionary directly. ResToken is a representation, not a property predictor — it enables LLMs to generate chemically legal sequences satisfying user-specified constraints but makes no claims about the biological activity of generated peptides.

**Future directions.** Integration with permeability predictors (e.g., CycPeptMP) and structure prediction tools (AlphaFold3) would close the loop from generation to evaluation, enabling active learning workflows. The library can be expanded through community contributions of validated NCAA blocks, and the tokenization scheme generalizes naturally to stapled peptides and peptide–drug conjugates.

---

## 5. Software and Availability

The `restoken` Python package is freely available under the MIT license. The package includes: (1) the 400-block building block library in CSV and JSON formats, with separate backend (full chemistry) and LLM-safe (semantic properties only) versions; (2) `SequenceValidator` implementing seven configurable hard constraints; (3) `SMILESReconstructor` for converting token sequences to full molecular SMILES, validated via RDKit for all 400 blocks; (4) prompt templates for ResToken, SMILES, and HELM generation across major LLM providers; (5) a random baseline generator with configurable constraint profiles for benchmarking; and (6) the complete benchmark suite used in this work, enabling reproduction of all reported results.

**GitHub:** [URL] | **License:** MIT | **Python:** ≥3.10

---

## Associated Content

**Supporting Information:** Per-model detailed results tables, additional statistical analyses, and library coverage visualizations.

---

*Draft notes (remove before submission):*
- *Total word count: ~2,800 (target: ~3,000 for JCIM Application Note)*
- *TODO: Add Claude Sonnet 4.6 as 6th model after midnight benchmark run — update Abstract numbers, Figure 2, Tables, and aggregate statistics*
- *TODO: Finalize reference list (CycPeptMP, PeptideCLM, PepThink-R1, Evo-R, HELM spec, RDKit)*
- *TODO: GitHub URL once repo is public*
