# Independent Review of ResToken NCAA Tokenization Manuscript

Reviewer: Woody  
Date: 2026-05-24  
Scope: manuscript drafts, supporting information, benchmark summaries, scripts, and raw output summaries under `/scratch/genesis/NCAA_tokenization/`; I did **not** read any files under `manuscript/review/`.

## Overall recommendation

**Major revision before submission.**

The paper has a clear and useful idea: represent NCAA cyclic peptides as closed-set, residue-semantic tokens so LLMs can generate syntactically legal residue sequences and reason over residue-level properties. This is well aligned with JCIM Application Note scope if the package is public and reproducible. However, the current manuscript overclaims in several places, and the benchmark definition needs tightening. The largest concern is that the reported SMILES validity appears to count RDKit-parseable molecules as valid even when they are not macrocyclic peptides, which makes several comparisons hard to interpret.

## Strengths

1. **Good representation-level contribution.** The dual-dictionary design is practical: LLM-facing semantic tokens plus backend chemical dictionaries is exactly the right engineering abstraction for NCAA design.

2. **Closed-set validity is a real advantage.** ResToken directly prevents nonexistent monomer hallucination and is especially useful for lead optimization where substitutions must remain inside a known building-block inventory.

3. **The EDIT/FROZEN benchmark is compelling.** The updated raw summaries show 100% frozen-position compliance across models and variable edit compliance, which supports the claim that discrete residue tokens make positional control easier.

4. **Property-constrained generation is the strongest scientific story.** Gemini Pro achieving 80–100% constraint satisfaction with 44–123x enrichment over random baselines is a strong result, provided the constraints and baselines are clearly defined.

5. **The limitations section is honest in several places.** The manuscript correctly states that ResToken is not a property predictor, has limited gamma/N-methyl coverage, and lacks experimental or computational validation of generated peptides.

## Major concerns

### 1. SMILES validity metric is not equivalent to cyclic peptide validity

The SMILES summaries include separate fields for `chemical_valid`, `has_macrocycle`, and `avg_amide_bonds`, but the manuscript reports `overall_valid` values that appear to equal chemical parse validity rather than cyclic-peptide validity. For example:

- Qwen 3.5 9B SMILES: 98.4% `overall_valid`, but `has_macrocycle = 0.0`.
- Gemma 3 12B SMILES: 92.2% `overall_valid`, but `has_macrocycle = 0.0`.
- Gemini Flash SMILES: 37.0% `overall_valid`, but `has_macrocycle = 0.76%`.
- GPT-4o SMILES: 45.2% `overall_valid`, but `has_macrocycle = 17.7%`.

For a cyclic peptide paper, RDKit-parseable small molecules should not be counted as valid cyclic peptide outputs. The current metric risks making SMILES look both better and worse in confusing ways: local models look highly valid because they emit simple parseable non-macrocycles, while the real failure is that they did not generate the requested modality.

**Required fix:** split the metric into at least:

1. parse validity;
2. macrocycle/ring validity;
3. peptide-like amide count or residue-count compliance;
4. overall task validity = parseable + macrocyclic + peptide-like + requested length/monomer constraints.

Then update Figure 2, Table 2, SI tables, and all textual claims. The key result may become even stronger for ResToken, but it must be stated correctly.

### 2. The comparison to HELM is underdeveloped and partly contradicts the narrative

The manuscript argues that HELM lacks semantic property encoding and therefore should underperform. However, the reported Exp1 results show HELM is often close to or better than ResToken for validity: Gemini Pro reaches 100% HELM validity, Qwen reaches 92.5%, and the mean is 84.5% versus 90.4% for ResToken. This is not a large universal gap.

Also, the prompt setup should be clarified. If the HELM prompt includes a table of monomer properties, then HELM is being augmented with semantic metadata and the claim that HELM cannot support property-aware reasoning is weakened. If the HELM prompt does not include properties, then the baseline may be artificially disadvantaged for property constraints.

**Required fix:** state exactly what information HELM receives. If HELM monomers are provided with metadata, call it `HELM+property table`; if not, acknowledge that the comparison tests a non-semantic HELM use case. Avoid broad claims that HELM cannot be made property-aware with an external monomer table.

### 3. Exp3 only benchmarks ResToken, so cross-representation claims about constraint reasoning are not fully supported

The manuscript concludes that ResToken enables property-aware reasoning that SMILES and HELM do not support, but SI Table S4 reports Exp3 only for ResToken. The notebook suggests some SMILES/HELM constrained runs were performed for permeable constraints, but these are not incorporated systematically into the manuscript.

**Required fix:** either:

1. add full SMILES and HELM constrained-generation results for all three profiles and all models; or
2. explicitly narrow the claim: “Within ResToken, capable models can exploit token metadata for constraint satisfaction.”

Without this, do not claim benchmarked superiority over HELM/SMILES for property-constrained design.

### 4. The manuscript overstates “ResToken achieves superior generation validity”

The mean validity is higher for ResToken, but the model-specific story is more nuanced:

- Gemini Pro: HELM 100% vs ResToken 98.5%.
- Qwen: SMILES 98.4% and HELM 92.5% vs ResToken 77.3%, although SMILES is likely invalid as cyclic peptide output.
- Gemma: SMILES 92.2% vs ResToken 85.6%, again likely non-macrocyclic.

The conclusion should be reframed around **task-valid cyclic peptide generation** and **closed-set monomer validity**, not generic validity.

**Required fix:** after fixing the SMILES task-validity metric, rewrite the Results heading and claims to reflect the corrected metrics. Use “ResToken improves task-valid NCAA sequence generation and positional control” rather than broad “superior generation validity” unless the corrected data support it.

### 5. Library curation and synthetic accessibility are insufficiently documented

The manuscript says the 400 blocks are “commercially available and synthetically accessible,” but `residue_tokens.csv` has no obvious source/vendor/provenance column. Reviewers will ask how these 400 were selected, whether stereochemistry was manually verified, and whether the set is reproducible.

**Required fix:** add a curation subsection and SI table/metadata fields covering:

1. source/provenance for each block;
2. inclusion/exclusion rules;
3. stereochemistry assignment method;
4. whether each block is commercial, literature-known, or designed;
5. any manual correction steps.

For JCIM, this is essential because the library itself is a main contribution.

### 6. “Constraint validator guarantees all generated sequences are chemically legal” is too strong

The validator checks token existence, length, charge, HBD/HBA, rotatable bonds, and backbone compatibility. It does not guarantee synthetic feasibility, macrocyclization feasibility, conformational viability, permeability, or biological activity. “Chemically legal” is acceptable only if narrowly defined as “composed of known library block IDs and passing configured rule checks.”

**Required fix:** replace “guarantees all generated sequences are chemically legal” with “guarantees all accepted sequences are composed of known library blocks and satisfy predefined rule-based constraints.”

### 7. Generated peptide diversity needs clearer chemical definition

The manuscript mentions diversity and SI reports uniqueness, but chemical diversity is not consistently reported for all representations. For SMILES/HELM, diversity is partly missing or not directly comparable. Sequence uniqueness alone is not enough, especially if local models generate repetitive simple molecules or non-macrocycles.

**Required fix:** define diversity metrics consistently:

1. sequence uniqueness for token outputs;
2. Morgan/Tanimoto diversity on reconstructed full structures for valid ResToken outputs;
3. analogous RDKit fingerprints for task-valid SMILES/HELM structures where possible;
4. report failures explicitly instead of leaving blanks.

### 8. The manuscript should distinguish representation contribution from model capability more sharply

Exp3 shows local models get 0% constraint satisfaction despite valid ResToken sequences. This is important: ResToken makes constraints visible, but reasoning remains model-limited. The Discussion says this, but the Abstract and Results still sound too representation-deterministic.

**Required fix:** revise the Abstract sentence “with up to 123-fold enrichment” to specify that this is achieved by Gemini 2.5 Pro, not by ResToken universally. Add model-capability caveat in the main claim.

## Minor comments

1. The Introduction cites Evo 2 as 2026 Nature and PepThink-R1 as 2025 NeurIPS/arXiv. Verify all bibliographic details before submission; some may be future/preprint entries.

2. The library statistics differ across files: the older outline says 19 functional classes and 38 NCY blocks, while the manuscript says 21 classes and 45 NCY. Ensure all drafts, figures, and tables use the same final v11 values.

3. “HELM lacks comprehensive NCAA coverage” should be phrased carefully. HELM is a notation; coverage depends on the monomer database provided. Better: “standard HELM deployments require explicit monomer dictionaries and do not intrinsically expose design-oriented semantic properties.”

4. The package name, GitHub link, and license should be verified as live before submission. JCIM reviewers may check immediately.

5. Define token naming conventions completely. The manuscript says uppercase/lowercase encodes chirality, but examples include mixed prefixes and IDs such as `S21`, `a27`, `A143`; make sure this rule is correct for all blocks.

6. The Application Note should not overpromise “bidirectional mapping” unless structure-to-token decomposition is implemented and validated, not just conceptually possible.

7. The Abstract says “while hiding structural details from the LLM.” This is accurate for ResToken, but if SMILES/HELM prompts include full dictionaries, explain that the representation comparison tests different information exposures.

8. Consider renaming “validity” to “acceptance rate” for ResToken validator outputs, because validity includes arbitrary user constraints in some contexts.

## Suggested revised central claim

A safer and stronger framing:

> ResToken converts NCAA cyclic peptide design into closed-set residue-token selection with explicit residue-level properties. This improves task-valid sequence generation and positional editing relative to atom-level SMILES generation, and enables capable frontier LLMs to perform enriched property-constrained sampling. The representation does not by itself solve property prediction or biological validation, and constraint satisfaction remains strongly model-dependent.

## Suggested additional analyses before submission

1. Recompute SMILES task-validity using parseable + macrocyclic + peptide-like criteria.
2. Add full constrained-generation baselines for HELM and SMILES, or remove cross-representation Exp3 superiority claims.
3. Add a provenance table for the 400 NCAA blocks.
4. Report confidence intervals for Exp2 and Exp3, not only Exp1.
5. Add a failure-mode figure/table: invalid IDs, wrong length, non-macrocycle, no edit, mode collapse.
6. Include a minimal reproducibility command in SI: how to regenerate tables/figures from raw outputs.

## Bottom line

The idea is publishable and useful, but the manuscript needs a stricter benchmark definition and more conservative claims. The highest-priority fix is the SMILES validity/task-validity issue; without it, reviewers may challenge the central comparison. Once corrected, the paper should emphasize ResToken’s true strengths: closed-set NCAA sequence generation, interpretable property metadata, precise positional control, and frontier-model-enriched constrained sampling.
