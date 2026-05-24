# Manuscript Review — ResToken JCIM Application Note

**Reviewer:** Assistant Agent
**Date:** 2026-05-24
**Documents reviewed:** JCIM_outline.md, JCIM_scaffold.md, methods_results_draft.md, supporting_information.md, references.md, all source code (restoken/src/), experiment scripts, raw JSON results, lab notebook.

---

## Overall Assessment

This is a well-engineered tool paper with a clearly identified problem (LLMs can't design NCAA cyclic peptides because existing molecular representations cause hallucination), a clean solution (curated closed-set tokens with semantic properties), and a systematic benchmark. The dual-dictionary architecture is genuinely clever. The code is clean and well-tested.

However, the manuscript has several significant issues that must be addressed before submission: an unfair baseline comparison that conflates representation quality with task difficulty, a missing constraint check that invalidates one experiment's results, overclaiming relative to what the data actually shows, and substantial JCIM formatting gaps.

**Recommendation:** Major revision required. The core contribution is publishable, but the framing and several experimental details need correction.

---

## Critical Issues

### C1. Rigid scaffold constraint checker is incomplete

The rigid_scaffold constraint profile communicated to the LLMs specifies three requirements:
1. Total rotatable bonds ≤ 18
2. At least 3 beta-backbone blocks
3. **Flexibility = low for at least half the residues**

However, the actual constraint checker (`run_benchmark.py:759-765`) only validates #1 and #2. The `flex_bin` property exists in `bb_dict_llm_v11.json` but is never loaded into the `Block` dataclass (`library.py`) and is never checked in `_check_profile_constraints`.

This means the reported rigid_scaffold constraint satisfaction rates (Gemini Pro 100%, GPT-4o 70%) are **overestimates** — sequences that violate the flexibility constraint are incorrectly counted as passing. The 123× enrichment figure in the abstract is derived from this incomplete check.

**Fix:** Either (a) add `flex_bin` to the Block dataclass and implement the missing check, then re-evaluate all Exp3 rigid_scaffold results, or (b) remove the flexibility constraint from the profile definition (both in the prompt and the paper) so that what the LLM is told matches what the validator actually checks.

### C2. Comparison is fundamentally unfair — must be foregrounded, not buried

The three representations differ not just in vocabulary structure but in **task difficulty**:
- **SMILES:** LLM must assemble a valid macrocyclic SMILES string character-by-character from atomic-level notation. This is a combinatorially hard generative task.
- **HELM:** LLM selects from the same 400 monomer codes but has no property metadata for constraint reasoning.
- **ResToken:** LLM selects from 400 pre-validated IDs AND gets an 8-property summary per token. This is a constrained selection task, not a generative task.

The methods_results_draft.md adds a "Comparison fairness" paragraph (commendably), but this critical caveat is buried at the end of the Methods section. When readers see "90.4% vs 69.0% vs 84.5%" in the Abstract, most will interpret this as ResToken being a superior representation, when in reality it's a fundamentally easier task disguised as the same comparison.

**Fix:** The fairness caveat must appear in the Abstract itself or the opening of the Results. Something like: "Because ResToken restricts the LLM to selecting from a finite set of pre-validated tokens rather than assembling molecular strings, the validity comparison reflects both representation design and inherent task difficulty." Without this, the paper will be criticized by reviewers #1 and #2 as misleading.

### C3. Central claim doesn't match the data

The Discussion opens with: "representation design and model capability jointly determine the quality of LLM-based NCAA peptide generation." This is better than earlier framings but still understates the data's message.

From Exp3 (the most practically relevant experiment):
- Only Gemini 2.5 Pro reliably satisfies constraints (80-100% across all profiles)
- Gemini Flash: 86% on permeable, 0% on charged/rigid
- GPT-4o: 30% on permeable, 0% on charged, 70% on n=20 (unreliable)
- Qwen 3.5 9B: 0% on everything
- Gemma 3 12B: 0% on everything

This means **4 out of 5 models fail completely on 2 out of 3 constraint profiles**. The tool's practical utility is currently limited to a single frontier model. The paper should state this directly: ResToken enables constraint-aware design but only when paired with top-tier frontier models; smaller or mid-tier models generate valid ResToken sequences but cannot reason about the semantic properties. This is actually an interesting finding that strengthens the paper — don't hide it.

### C4. GPT-4o rigid scaffold sample size (n=20) is too small for the claims built on it

GPT-4o produced only 24 parsed sequences for rigid_scaffold (vs. 100 requested), of which 20 were valid, and 14/20 passed constraints (70%). The manuscript does note "should be interpreted with caution," but this result still contributes to:
- "up to 123-fold enrichment" (abstract) — actually that's Gemini Pro, not GPT-4o
- "mean enrichment 104.9× for rigid scaffold (2 models)" — this is the mean of 123.5× (Pro) and 86.4× (GPT-4o), where the GPT-4o number has a 95% CI of roughly 46-88%

A 14/20 result is consistent with true rates anywhere from ~46% to ~88% (exact binomial CI). This is too noisy to compute meaningful enrichment. Either report the Pro result alone, or caveat the aggregate much more strongly.

---

## Major Issues

### M1. Exp2 tests only ResToken — claims about SMILES/HELM controllability are unsupported

The manuscript claims "structured tokens enable precise positional control" and the scaffold doc references "edit bleed" in SMILES. But Exp2 was run only on ResToken. The methods note this would be "a qualitatively different task" for SMILES, but that's exactly why the comparison would be informative.

Without SMILES/HELM data for Exp2, the claim is weaker than presented. The paper should either: (a) run Exp2 on SMILES/HELM to support the claim, or (b) reframe from "ResToken enables precise control vs. alternatives" to "ResToken achieves 100% frozen compliance across all tested models." The `build_exp2_smiles_prompt` and `build_exp2_helm_prompt` functions already exist in the benchmark script — the experiment is implementable.

### M2. "Mean validity" across 5 models is statistically meaningless

Averaging validity rates across 5 hand-picked models (3 API frontier, 2 local open-weight) and reporting "mean ± std" implies a sampling distribution from some population of models. But this is a convenience sample with no defined population. The variance tells you about the spread across these specific models, not about any general property of "LLMs."

**Fix:** Report individual model results prominently (Table 2 does this well). The "mean" can appear as a summary statistic but should not be the lead number in the Abstract. Instead, lead with the range: "ResToken validity ranged from 77% to 99% across five LLMs, compared to 37-98% for SMILES."

### M3. Uniqueness/diversity is underreported

GPT-4o and Gemma 3 achieve respectable ResToken validity (98.1% and 85.6%) but with uniqueness of only 0.58 and 0.56 — meaning ~40-44% of generated sequences are duplicates. This is a significant practical limitation that is mentioned only in the SI. For a design tool, diversity matters as much as validity. The main text should discuss the validity-diversity tradeoff explicitly.

Similarly, Qwen's 98.4% SMILES validity is flagged as "counterintuitive" but the explanation (shorter strings, low diversity, 0% macrocycles) should appear in the main text, not just the SI note. Without this, readers will be confused why a 9B model outperforms GPT-4o on SMILES.

### M4. Missing references for all 5 LLMs

None of the 5 benchmark models have formal citations (model cards, tech reports, or papers). For a JCIM submission, each model needs a reference. Gemini models need the Gemini technical report; GPT-4o needs the OpenAI system card; Qwen needs the Qwen 2.5 paper; Gemma 3 needs the Gemma technical report.

### M5. Scaffold vs. draft model set discrepancy (Sonnet → Gemini)

The scaffold lists Claude Sonnet 4.6 and Qwen-2.5-72B as models. The actual experiments use Gemini Pro/Flash, GPT-4o, Qwen 3.5 9B, and Gemma 3 12B. The lab notebook explains the switch (Sonnet was rate-limited, Qwen 72B was too large for available GPUs). The manuscript should explain the model selection rationale in Methods — why these 5 specific models?

### M6. Missing `polarity_bin` and `flex_bin` in Block dataclass

The paper claims each token encodes "eight semantic properties" including polarity and flexibility. These exist in `bb_dict_llm_v11.json`, but the `Block` class in `library.py` doesn't load them. The `SequenceValidator` cannot validate polarity or flexibility constraints.

This means the shipped package cannot programmatically check two of the eight advertised properties. Users who install `restoken` and try to filter blocks by `polarity_bin` or `flex_bin` will find these attributes missing.

**Fix:** Add `polarity_bin` and `flex_bin` to the `Block` dataclass and load them from the LLM dictionary (they're not in the backend dict). Or document that these are LLM-prompt-only properties not available in the Python API.

---

## Minor Issues

### m1. Functional class count inconsistency
- Scaffold abstract: "20 functional classes"
- Draft manuscript body: "21 functional classes"
- Lab notebook: 19 classes with nonzero counts listed, but Table 1 line says "21"
- The actual data has: G, A, V, L, I, F, P, S, T, K, D, E, C, H, Y, N, Q, R, W, M, ncaa = 21 classes

Pick one number and use it consistently. Verify against the actual data.

### m2. N-modification counts
Table 1 says "N-cyclic (NCY) = 38 (alpha) + 4 (beta) + 3 (gamma) = 45 total." But the lab notebook updated stats say "NCY (45)" and "NO (353), NCY (45), NXX (6), NME (2)." Note 353+45+6+2 = 406, which is 6 more than 400. This is because NXX is a separate category from NCY and NME. The manuscript Table 1 doesn't mention NXX at all — 6 blocks with NXX modification are unaccounted for in the paper's N-modification breakdown.

### m3. Rounding inconsistencies
- Abstract "123-fold" vs SI "123.5×" vs JSON 123.46×
- Manuscript "15×" for GPT-4o permeable; JSON 14.70×
- Manuscript "42×" for Flash permeable; JSON 41.61×

Pick a rounding convention (1 decimal for enrichment ≥10×, 2 for <10×) and apply consistently.

### m4. Repetitive phrasing
"Closed-set vocabulary prevents/eliminates hallucination" appears in the Abstract, Introduction, Methods §2.2, and Discussion — four times. Keep it in the Abstract and one body section. Readers get the point.

### m5. Tense inconsistency in Methods
The Methods section mixes present tense ("The ResToken library comprised…" is past, but "Each token carried…" is past while "Each block was assigned…" is past) — actually the draft is mostly consistent in past tense. But the scaffold uses present tense throughout. Whichever version becomes the submission, standardize to past tense for Methods.

### m6. Gemma rigid_scaffold anomaly unexplained
Gemma produced 431 sequences from 100 requested, with only 1 unique (extreme mode collapse). The manuscript reports 0% constraint satisfaction, which is correct. But the 431-from-100 anomaly should be noted — it suggests the model entered a degenerate repetition loop. This is interesting behavior worth a sentence in the SI.

### m7. `pyproject.toml` has no core dependencies
The `dependencies = []` field means `pip install restoken` won't pull in any packages, but the code requires at minimum `rdkit` for reconstruction validation. Consider adding runtime dependencies or at least documenting them clearly.

---

## JCIM Format Issues

These must be fixed before submission:

1. **No in-text citations.** The draft uses informal references ("[cite: ...]" or just mentions names). JCIM uses superscript numbered ACS-style citations.
2. **Abstract exceeds 150 words.** The draft abstract is ~156 words. Cut ~6 words.
3. **No TOC/abstract graphic.** Required for JCIM.
4. **Missing Keywords section.** Need 3-6 keywords after the abstract.
5. **Section numbering.** JCIM Application Notes typically don't number sections. Remove "1.", "2.", etc.
6. **No Author Information / ORCID section.**
7. **GitHub URL placeholder "[URL]" in Software section.** Must be filled.
8. **Ref 8 (PepThink-R1) is an arXiv preprint.** ACS generally doesn't accept unrefereed preprints. Check if it has been published since.
9. **Refs 9, 14 use "et al."** — ACS requires full author lists or journal-specific truncation rules.
10. **Ref 13 (RDKit) is a bare URL.** Format as software citation with version and DOI.
11. **Multiple orphan references** (Refs 4, 5, 6 never cited inline in the draft).
12. **Missing inline citations** for SELFIES (Ref 11), Bootstrap (Ref 15), and all LLMs.

---

## Reference Issues

### Missing citations for claims in the text
- All 5 benchmark LLMs (no model cards/tech reports cited)
- TxGemma (mentioned in Methods/SI, no reference)
- "Clinical successes such as cyclosporine" — Ref 5 exists but is never cited inline
- "mRNA-display-derived candidates from the RaPID platform" — Ref 4 exists but never cited inline
- SwissSidechain or similar NCAA database (relevant prior art, not cited)
- ChemBERTa (mentioned in scaffold, not in draft)
- Morgan fingerprints methodology

### Potentially wrong reference
- PeptideCLM is cited as Ref 7 (Feller & Wilke 2025, "Peptide-Aware Chemical Language Model Successfully Predicts Membrane Diffusion of Cyclic Peptides"). This is a membrane permeability prediction paper, not the AMP generation paper usually associated with "PeptideCLM." Verify this is the intended reference.

---

## Suggestions for Strengthening

1. **Add a composite metric.** A metric that combines validity × uniqueness × constraint satisfaction would give a single number per model-representation pair that captures the full picture. Currently, high validity with low uniqueness (Gemma SMILES: 92% valid, 31% unique, 0% macrocycles) looks deceptively good.

2. **Discuss why HELM performs well.** HELM achieves 100% validity for Gemini Pro — identical to ResToken. The advantage of ResToken over HELM is mainly in constraint satisfaction (which only works for Pro anyway). Be honest about when HELM is "good enough."

3. **Error analysis.** What specific errors do models make? The lab notebook has rich qualitative data (Gemma hallucinating sequential IDs like A190-A216, Qwen producing valid but non-cyclic SMILES). A brief error taxonomy in the main text would strengthen the paper significantly.

4. **Report n_parsed vs n_requested systematically.** The fact that some models produce far fewer (GPT-4o: 160/200) or far more (Gemma rigid: 431/100) parsed sequences than requested is important context that the SI tables show but the main text doesn't discuss.

5. **Consider adding Claude Sonnet.** The lab notebook shows the Sonnet benchmark was blocked by rate limiting. If this can be re-run, having 4 API models + 2 local models would strengthen the model diversity argument. If not, explain the absence.

---

## Code Quality Notes

The codebase is clean and well-organized:
- `BlockLibrary`, `SequenceValidator`, `SMILESReconstructor` are well-designed, testable, properly typed
- Test suite covers legal/illegal sequences, constraint violations, SMILES round-trip
- Benchmark scripts are parameterized and reproducible
- Data versioning (v11) is tracked

Two concerns:
1. `sys.path.insert(0, ...)` in test file — use proper package imports instead
2. The `get_llm_entry` method in `BlockLibrary` does a linear scan of the LLM dict on every call — should index by ID at load time (minor perf issue, not publication-blocking)

---

## Summary of Required Changes

| Priority | Issue | Section |
|---|---|---|
| Critical | Fix rigid_scaffold constraint checker (flex_bin missing) | C1 |
| Critical | Foreground comparison fairness caveat in Abstract/Results | C2 |
| Critical | Reframe central claim to match Exp3 data | C3 |
| Critical | Caveat or remove GPT-4o rigid scaffold from aggregates | C4 |
| Major | Run Exp2 on SMILES/HELM or reframe controllability claim | M1 |
| Major | Don't lead with "mean validity" — use range or per-model | M2 |
| Major | Discuss uniqueness/diversity tradeoff in main text | M3 |
| Major | Add LLM model references | M4 |
| Major | Add polarity_bin/flex_bin to Block dataclass | M6 |
| Minor | Fix functional class count inconsistency | m1 |
| Minor | Account for NXX blocks in Table 1 | m2 |
| Format | Fix all JCIM formatting issues (citations, abstract length, etc.) | §Format |
