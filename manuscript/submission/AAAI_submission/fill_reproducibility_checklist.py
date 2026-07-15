#!/usr/bin/env python3
"""Fill the AAAI reproducibility checklist from the ResToken manuscript evidence."""

from pathlib import Path
import argparse

from docx import Document


RESPONSES = {
    2: "partial — Figure 1 and the Methods describe the dual-dictionary architecture and rule-based validator, but no pseudocode is included.",
    3: "yes — Results are separated from interpretation, and qualifications and unvalidated claims are explicitly identified in the Discussion and Limitations.",
    4: "yes — Background references are provided for cyclic peptides, molecular representations, LLM-based molecular design, permeability prediction, and statistical methods.",
    6: "no — The paper introduces and evaluates a representation and software framework; it does not claim formal theoretical contributions.",
    8: "NA — The paper does not make theoretical contributions.",
    9: "NA — The paper does not make theorem-style claims.",
    10: "NA — No formal proofs are required or claimed.",
    11: "NA — The paper does not present theoretical results requiring proof sketches.",
    12: "NA — No theoretical tools are used to establish formal claims.",
    13: "NA — The paper makes empirical rather than theoretical claims.",
    14: "NA — The paper does not use experiments to eliminate or disprove theoretical claims.",
    16: "yes — The work uses CycPeptMPDB, Protein Data Bank structures, a curated 400-building-block library, and generated benchmark datasets.",
    18: "yes — The manuscript explains why CycPeptMPDB is used for independent permeability evaluation and why the selected PDB complexes test structure-guided design.",
    19: "partial — The novel building-block library and generated benchmark outputs are provided in the anonymized repository and Supplementary Information rather than a dedicated data appendix.",
    20: "yes — The novel library and benchmark outputs will be publicly released; the software is under the MIT license and will be archived with a citable Zenodo DOI.",
    21: "yes — CycPeptMPDB and the Protein Data Bank structures are cited with a literature reference or accession code.",
    22: "yes — CycPeptMPDB and the PDB structures used in the study are publicly available.",
    23: "NA — The study does not rely on a third-party dataset that is unavailable to the public.",
    25: "yes — The paper reports LLM generation benchmarks, rule-based validation, permeability prediction, and structure-based scoring experiments.",
    27: "partial — Final task settings and model configurations are reported, but the number and range of every development-stage hyperparameter trial are not fully enumerated.",
    28: "partial — Data-processing and reconstruction code is available in the anonymized repository, but it is not reproduced in a code appendix.",
    29: "partial — The anonymized repository contains the benchmark and analysis code, but the source code is not reproduced in a code appendix.",
    30: "yes — The ResToken package and experiment code will be publicly released under the MIT license and archived with a Zenodo DOI.",
    31: "partial — The released package documents the principal components and constraints, but not every implementation step is cross-referenced to the manuscript.",
    32: "no — The manuscript does not currently provide a complete seed-setting procedure for all stochastic generation and sampling steps.",
    33: "partial — Model names, quantization, Python requirement, and RTX 4090 use are stated, but CPU, memory, operating system, and all library versions are not fully specified.",
    34: "yes — Generation validity, constraint enrichment, positional compliance, and permeability shift are formally defined and motivated by the study objectives.",
    35: "yes — The manuscript reports N=200 per representation/model for Task 1, 200 variants per model for Task 2, N=100 per constraint profile for Task 3, and 100,000 random baseline samples.",
    36: "yes — Results include ranges, standard deviations, bootstrap 95% confidence intervals, medians, Cliff's delta, and distributional comparisons.",
    37: "yes — Representation differences use Pearson chi-square tests with Bonferroni correction, and permeability contrasts use Mann–Whitney tests with effect sizes.",
    38: "partial — Final generation settings and model configurations are reported in the Methods and Supplementary Information, but not every predictor-training and scoring parameter is listed in the main paper.",
}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("input", type=Path)
    ap.add_argument("output", type=Path)
    args = ap.parse_args()

    doc = Document(args.input)
    for idx, response in RESPONSES.items():
        paragraph = doc.paragraphs[idx]
        marker = "Type your response here"
        if marker not in paragraph.text:
            raise ValueError(f"Expected placeholder missing in paragraph {idx}")
        question = paragraph.text.replace(marker, "").rstrip()
        paragraph.clear()
        paragraph.add_run(question + " ")
        answer_run = paragraph.add_run("Response: " + response)
        answer_run.bold = True

    doc.core_properties.last_modified_by = "Yike"
    doc.save(args.output)

    check = Document(args.output)
    text = "\n".join(p.text for p in check.paragraphs)
    assert "Type your response here" not in text
    assert len(RESPONSES) == 31


if __name__ == "__main__":
    main()
