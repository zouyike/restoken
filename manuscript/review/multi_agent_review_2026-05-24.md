# Multi-Agent Manuscript Review — 2026-05-24

## Review 1: Data Consistency

**CRITICAL:** None found. All numbers traceable to JSON are correct.

**WARNING (not verifiable from JSON):**
- Library composition counts (alpha=250, beta=101, gamma=49; L=292, D=86, achiral=22) — not in statistical JSON
- 21 functional classes claim — not in JSON
- "only 2 N-methylation blocks" — not in JSON
- Gamma covering "10 of 21 functional classes" — not in JSON
- SI Fig S1 per-class block counts (I=88, F=46, A=46, T=42, P=38) — not in JSON
- SI Fig S2 "neutral charge 83.5%" — not in JSON
- SI Table S2 uniqueness values — no uniqueness data in JSON

**MINOR (rounding/formatting):**
- Abstract "123-fold" vs SI "123.5x" vs JSON 123.46x — inconsistent rounding across documents
- Manuscript "15x" for GPT-4o permeable; JSON=14.70x, SI=14.7x — rounded up
- Manuscript "42x" for Flash permeable; JSON=41.61x, SI=41.6x — rounded up
- Manuscript chi-square "53.4" vs SI "53.43" — different decimal precision
- Manuscript "30%" for GPT-4o permeable; JSON=30.4% — rounded down
- Manuscript "86%" for Flash permeable; JSON=86.1% — rounded down
- Manuscript SMILES uniqueness range "0.29-0.31"; actual Gemma value is 0.314, not 0.31

---

## Review 2: Scientific Rigor

**CRITICAL:**

- **Overclaimed central finding.** Discussion states "representation design, rather than model scale, is the primary bottleneck." Exp3 shows constraint satisfaction is entirely model-dependent (only Gemini Pro succeeds broadly; 9-12B models score 0%). The paper's own results demonstrate model capability is at least co-equal. Must rewrite.

- **SMILES/HELM comparison is not fair.** SMILES prompt asks LLM to assemble full macrocyclic SMILES (harder task than picking IDs). HELM gets same 400-monomer table but no property metadata. ResToken gets both closed vocabulary AND semantic properties. Comparison conflates representation quality with task difficulty. Must discuss this confound.

- **GPT-4o rigid scaffold (n=20) reported as 70% with 86x enrichment.** 14/20 passing is statistically unreliable (95% CI ~46-88%). Abstract's "up to 123-fold enrichment" cherry-picks best cell. Both need explicit caveats.

**MAJOR:**

- No confidence intervals or significance tests for Exp3.
- "Mean validity" across 5 models is misleading (convenience sample, not random).
- Exp2 tests only ResToken — "edit bleed" claim for SMILES has no data.
- Uniqueness penalty ignored in validity comparison (Qwen/Gemma high SMILES validity but 0.29-0.31 uniqueness).
- "Mean enrichment" for charged binder based on single model (44x from Pro alone) is not meaningful aggregate.

**MINOR:**

- Abstract claims "90.4% mean validity" without noting unweighted average.
- Bonferroni for 10 comparisons appropriate but 3-way test not corrected.
- Limitations doesn't mention lack of wet-lab or computational validation (docking/MD).
- TxGemma exclusion rationale should appear in Methods, not just SI.
- Parsed counts vary widely (Gemma rigid_scaffold: 431 from 100 requested) — unexplained.

---

## Review 3: JCIM Format Compliance

**REQUIRED:**

- No in-text citations anywhere. Need superscript numbered citations (ACS style).
- Abstract is 156 words; limit is 150. Cut ~6 words.
- No keywords section. Need 3-6 keywords.
- No TOC/abstract graphic.
- Ref 8 (PepThink-R1) is arXiv preprint — ACS generally doesn't accept unrefereed preprints.
- Ref 8 missing DOI.
- Refs 9, 14 use "et al." — ACS requires full author lists or journal-specific truncation.
- Ref 13 (RDKit) is URL, not ACS format.
- GitHub URL placeholder "[URL]" must be filled.
- Draft notes (lines 111-115) must be removed.

**RECOMMENDED:**

- Ref 2 missing issue number.
- Ref 1 article number format — verify correct for Angew. Chem.
- Figure 1b never referenced in text (1a is, but not 1b).
- No Author Information / ORCID section.
- Remove section numbering (JCIM App Notes don't typically number sections).

**OPTIONAL:**

- Ref 15 (bootstrap book) uses ISBN, consider DOI.
- More specific SI cross-references from main text.
- Define all abbreviations in abstract (LLM, HELM, NCAA).

---

## Review 4: Writing Quality

**HIGH:**

- "Closed-set vocabulary eliminates hallucination" repeated 4 times (Abstract, Intro, Methods 2.2, Discussion). Keep in Abstract + one body section.
- "Representation not model" thesis stated nearly verbatim in both Intro and Discussion. Consolidate.
- Overlong sentence (62 words) about NCAA benefits — split after "chemical space."
- Overlong HELM sentence — split at em-dash.

**MEDIUM:**

- Tense inconsistency in Methods: present tense ("comprises," "span") should be past tense.
- "Confirms" overclaims from N=200 — soften to "suggests" or "indicates."
- "Edit bleed" — undefined jargon in quotes.
- Passive voice cluster in task descriptions (lines 47-55).
- "Paradoxically" — explanation removes the paradox. Use "Counterintuitively" or just state finding.

**LOW:**

- "Therapeutic sweet spot" — colloquial for JCIM.
- "Encrypted within the notation" — metaphor mismatch, use "obscured within."
- Abstract lists 7 properties for "eight designable properties" (omits functional class).
- "Thin" gamma coverage — use "limited" or "sparse."
- Final sentence joins two independent claims with "and" — split.

---

## Review 5: Reference Completeness

**CRITICAL (cited but no reference):**

- PeptideCLM — cited in Intro but Ref 7 (Feller & Wilke) may be wrong paper (membrane diffusion, not AMP generation).
- TxGemma — cited in Methods with no reference.
- All 5 LLMs (Gemini Pro/Flash, GPT-4o, Qwen 3.5, Gemma 3) — no model card/tech report references.

**CRITICAL (orphan references):**

- Ref 4 (Goto & Suga 2021, RaPID) — never cited in text.
- Ref 5 (Borel et al. 1976, cyclosporine) — only alluded to, never formally cited.
- Ref 6 (Bagal et al. 2022, MolGPT) — never cited.

**MAJOR (missing citations):**

- SELFIES (Ref 11 exists but never cited inline).
- Cyclosporine/mRNA-display claims lack inline citations.
- Bootstrap CI method (Ref 15 exists but never cited inline).
- Chi-square/Bonferroni — no methodological reference.
- No prior NCAA library/database paper (SwissSidechain, PepBDB).

**MINOR:**

- Zero numbered citations in manuscript — all references informal.
- Ref 8 is arXiv preprint.
- Ref 9 (Evo 2, Nature 2026) — verify 2026 date.
- Ref 13 (RDKit) lacks version number.
