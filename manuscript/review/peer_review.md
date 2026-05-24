# Peer Review Report

**Manuscript:** ResToken: A Residue-Semantic Token Library Enabling LLM-Based Design of Noncanonical Cyclic Peptides
**Journal:** Journal of Chemical Information and Modeling (Application Note)

## Summary

The authors present "ResToken", a novel residue-semantic token library designed to facilitate the generation of noncanonical amino acid (NCAA) cyclic peptides by Large Language Models (LLMs). The core problem addressed is that existing molecular representations either cause LLMs to hallucinate invalid chemistry (SMILES) or fail to provide the physicochemical context necessary for property-aware design (HELM). By decoupling the LLM-facing semantic properties (via a token dictionary) from the chemical backend (full SMILES/SELFIES), ResToken ensures a closed-set vocabulary and constraint-aware generation. 

The manuscript is well-written, the problem is highly relevant to modern computational chemistry and AI-driven drug discovery, and the dual-dictionary architectural approach is elegant and practical. The benchmarking across frontier and local open-weights models is rigorous, offering a candid look at both the successes of ResToken and the current reasoning limitations of smaller LLMs. I recommend the manuscript for publication as an Application Note after the following major and minor comments are addressed.

## Major Comments

1. **Library Composition and Functional Bias:** The library currently consists of 400 building blocks, which is a solid foundation. However, the authors note that N-methylation is represented by only 2 blocks. Given that N-methylation is one of the most critical structural modifications for improving the membrane permeability of cyclic peptides, having only 2 blocks severely limits the practical utility of the library for permeability-constrained design tasks (which is the exact focus of Experiment 3 and Section 4). The authors should either expand the library to include a more representative set of N-methylated blocks or, at minimum, provide a clear programmatic mechanism (e.g., a script) in the `restoken` package for users to dynamically generate and automatically annotate new N-methylated blocks.

2. **HELM Baseline Fairness in Experiment 3:** In the property-constrained generation task, the authors state that HELM fails because it lacks semantic properties. It is necessary to clarify whether the LLM prompt for the HELM baseline included a property lookup table equivalent to the one provided for ResToken. If no property metadata was provided in the HELM prompt, the comparison conflates the *representation format* with the *information provided in the prompt*. Discussing this inherent difference and acknowledging whether the LLM could realistically succeed zero-shot without an explicit lookup table would strengthen the fairness and context of the benchmark.

3. **Macrocyclization Feasibility:** The sequence validator successfully checks for backbone transition compatibility (e.g., preventing illegal beta-gamma transitions). However, it is unclear if the validator or the SMILES reconstruction step guarantees that the resulting cyclic peptide is geometrically and physically feasible (e.g., assessing ring strain or steric clashes during macrocyclization). The authors should briefly discuss how users might assess the 3D cyclization feasibility post-generation to ensure the validity metric translates to physically realizable peptides.

4. **Context Window Impact for Local Models:** The manuscript states that the LLM dictionary consumes approximately 11,000 tokens. While this fits comfortably within the context windows of frontier models (GPT-4o, Gemini), this is a demanding prompt length for local 9B–12B models, potentially exceeding or pushing the limits of their effective context (e.g., standard Qwen 3.5 has an 8k default context, though extended variants exist). The authors should briefly discuss whether this large prompt length contributed to the failure of local models in Experiment 3, potentially due to "lost in the middle" context phenomena or simply exceeding their trained reasoning capacity.

## Minor Comments

1. **Inconsistency in Functional Classes:** There is a minor numerical inconsistency regarding the functional classes. The Abstract and Methods sections mention "21 functional classes," but the Outline (Section 2.1) states "19 functional classes." The authors should reconcile this number throughout the text and Supporting Information.
2. **Collapsed SMILES for Local Models:** The observation that local models achieved artificially high SMILES validity (92–98%) by collapsing to repetitive, simple strings (uniqueness 0.29–0.31) is an excellent finding. Providing a specific example of these repetitive, collapsed SMILES strings in the SI would be highly illustrative of this LLM failure mode.
3. **Diversity Metric with HELM:** In Experiment 1, it is noted that Gemini Pro achieved 100% validity with HELM but exhibited lower sequence diversity compared to ResToken (SI Table S2). Adding a brief sentence analyzing why the HELM representation might restrict the generative diversity of the model (e.g., limited monomer vocabulary known to the LLM) would enrich the discussion.

## Recommendation
Minor Revision
