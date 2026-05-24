# Review of ResToken Manuscript

Reviewer: Engineer / Codex  
Date: 2026-05-24  
Scope: Manuscript drafts, SI, figures, tables, benchmark scripts, result summaries, and package code under `/scratch/genesis/NCAA_tokenization/`. I did not read existing files in `manuscript/review/`.

## Overall Assessment

This is a strong and timely application-note concept. The core idea is valuable: represent noncanonical amino-acid cyclic peptide designs as closed-set, property-annotated residue tokens rather than forcing LLMs to generate atom-level strings. The package already has a usable 400-block library, validators, prompt templates, benchmark outputs, and manuscript-ready figures. The practical positioning as a representation/tool paper is appropriate for JCIM.

My recommendation is **major revision before submission**. The paper's direction is sound, but several headline claims are currently vulnerable because benchmark denominators, validity definitions, and baseline framing are inconsistent. These are fixable, but they should be fixed before external review.

## Major Issues

### 1. Validity denominators are not stable enough for headline claims

The Methods state that each model generated `N=200` sequences per representation, but the reported validity rates are calculated over parsed sequences, and parsed counts vary widely: for example, Gemini Flash SMILES has 395 parsed outputs from 200 requested, Qwen HELM has only 53, and Qwen ResToken has 256. This makes cross-model and cross-representation validity comparisons difficult to interpret.

Report three separate quantities everywhere:

- requested outputs
- parsed candidate sequences
- valid candidates

For headline validity, prefer a per-requested-output success rate, or at least make the parsed-denominator rate explicit in every table, figure, and abstract-level statement.

### 2. The phrase "seven hard constraints" overstates what Exp1 validity measures

The manuscript and SI describe validity as passing seven checks, including charge, HBD, HBA, rotatable bonds, and backbone compatibility. In the benchmark outputs, Exp1 ResToken validity appears to be format + ID + length validity. SMILES validity is RDKit parseability, and HELM cyclization validity appears to be computed separately rather than always included in `overall_valid`.

This is a serious wording issue. The paper should distinguish:

- representation syntax validity
- token/monomer identity validity
- task length validity
- cyclic peptide/macrocycle validity
- property-constraint validity

Avoid saying Exp1 sequences passed all seven rule checks unless the code actually applies all seven checks in that experiment.

### 3. SMILES validity is misleading unless macrocycle validity is included

For local models, SMILES `overall_valid` can be high while `has_macrocycle = 0`. Qwen 3.5 9B and Gemma 3 12B generate RDKit-parseable molecules, but not cyclic peptides. Since the task is cyclic peptide generation, a chemically valid non-macrocycle should not be counted as task-valid in the main comparison.

The paper should either redefine SMILES validity as "RDKit-parseable molecule validity" or count only valid macrocyclic peptide outputs for the main task-validity comparison. Otherwise the table gives a false impression of what the SMILES baseline achieved.

### 4. Property-constrained superiority over SMILES/HELM is under-supported

The draft says ResToken is benchmarked against SMILES and HELM across three tasks, but property-constrained constraint satisfaction is reported mainly for ResToken. Some SMILES/HELM constrained runs exist only for selected profiles. This supports "ResToken enables property-constrained generation," but not a full claim that ResToken outperforms SMILES and HELM across all constrained profiles.

Either run the full constrained baseline set for SMILES and HELM, or narrow the claim. A safe version is: "We evaluate constrained generation with ResToken and include limited baseline comparisons; complete constrained cross-representation benchmarking remains future work."

### 5. HELM baseline framing is too broad

The manuscript says HELM lacks semantic property encoding, but the benchmark prompt supplies HELM monomer properties to the model. A reviewer may object that the baseline is not plain HELM; it is HELM plus an external property table.

The fairer claim is that ResToken integrates semantic metadata into the tokenization scheme, while HELM requires auxiliary monomer metadata to expose comparable information to an LLM. Also revise Figure 1b, which currently implies HELM is "monomer IDs only."

### 6. Rigid scaffold prompt and evaluator are inconsistent

The rigid scaffold prompt includes a low-flexibility requirement, but the evaluator and random baseline appear to use rotatable-bond and beta-backbone criteria without explicitly enforcing "low flexibility for at least half the residues." This means the reported rigid-scaffold enrichment may correspond to a different task than the one given to the models.

Make the prompt, evaluator, manuscript profile definition, and random baseline use exactly the same criteria. Then regenerate the rigid-scaffold numbers.

### 7. Enrichment claims need denominator and model caveats

"Up to 123-fold enrichment" is true for the best Gemini 2.5 Pro rigid-scaffold cell, but it is not representative across models. Local models achieve 0% constraint satisfaction, and GPT-4o's rigid-scaffold result is based on only 20 valid sequences from 100 requested.

In the abstract, write this as a model-specific best-case result, e.g. "Gemini 2.5 Pro achieved up to 123-fold enrichment..." The Results should also report pass/requested alongside pass/valid.

### 8. Library provenance is not reproducible enough

"Curated from commercially available and synthetically accessible noncanonical amino acids" is not sufficient for a resource paper. Reviewers will ask:

- which vendors, catalogs, or databases were used
- inclusion and exclusion rules
- duplicate and stereoisomer handling
- how functional class labels were assigned
- how bulk, polarity, flexibility, HBD/HBA, and charge bins were computed
- what manual curation was performed after automated annotation

This can be solved in Methods or SI. The paper does not need to over-explain every block, but the curation pipeline must be reproducible.

### 9. "Chemically legal" and "guarantees" language is too strong

The validator guarantees that outputs use known IDs and satisfy predefined rule checks. It does not guarantee synthesis feasibility, macrocyclization feasibility, conformational viability, permeability, or activity. Replace broad wording like "chemically legal sequences" with "library-valid sequences passing predefined rule checks."

This matters because the paper's own limitations state that no docking, MD, predictor, wet-lab, permeability, or activity validation has been performed.

## Minor Issues

- The draft should avoid saying "HELM lacks property-aware reasoning" without clarifying that the benchmark HELM prompt was given external property metadata.
- The statistical methods say Pearson chi-square tests, but the code uses `chi2_contingency(..., correction=True)`, i.e. Yates-corrected chi-square. Report the correction or change the code.
- Random baseline methods say 100,000 uniform-random 6-mers, while the result files show rejection sampling until 1,000 accepted sequences with profile-specific attempt counts. Align the text with the implementation.
- TxGemma is excluded from most tables after producing zero parseable output. That is acceptable, but the exclusion rule should be stated consistently.
- The permeability profile combines NME and NCY, but the library contains only 2 NME blocks. Call this a heuristic low-HBD/N-modified profile rather than true permeability design.
- Figure 2 needs denominator labels or a caption explaining parsed versus requested rates. Error bars should be defined.
- Figure 1a is conceptually clear but visually large for an application note; consider simplifying if the final article has strict figure limits.
- The GitHub URL and package availability should be verified before submission.

## Suggested Revision Plan

1. Recompute or re-tabulate all benchmark metrics with requested, parsed, valid, and constraint-passing counts.
2. Define validity separately for ResToken, SMILES, and HELM, and make task-validity include cyclic peptide/macrocycle validity where relevant.
3. Decide whether Exp3 will be ResToken-only or a full cross-representation benchmark. Then align the claims accordingly.
4. Fix the rigid-scaffold profile mismatch and regenerate that profile.
5. Add a reproducible library-curation subsection to Methods or SI.
6. Replace guarantee/chemical-legality language with precise rule-check language.
7. Update Figure 1b and Figure 2 captions to match the actual benchmark setup.

## Bottom Line

The work is worth pursuing and likely publishable as a JCIM application note, but the current manuscript should not be submitted yet. The central idea is strong; the weak point is benchmark hygiene. If the denominator, validity-definition, and baseline-framing issues are cleaned up, the paper will become much more defensible.
