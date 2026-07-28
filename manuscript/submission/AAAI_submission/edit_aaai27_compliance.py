#!/usr/bin/env python3
"""Create AAAI-27-compliant ResToken manuscript, SI, and checklist DOCX files.

The edit removes author-owned web references, moves detailed material to the SI,
and shortens the main paper without changing margins, type size, or line spacing.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
from pathlib import Path

from docx import Document
from docx.text.paragraph import Paragraph


def remove_paragraph(paragraph: Paragraph) -> None:
    element = paragraph._element
    element.getparent().remove(element)
    paragraph._p = paragraph._element = None


def replace_text_preserve_drawings(paragraph: Paragraph, text: str) -> None:
    """Replace paragraph text while retaining floating/inline drawing anchors."""
    nodes = paragraph._p.xpath(".//w:t")
    if not nodes:
        paragraph.add_run(text)
        return
    nodes[0].text = text
    for node in nodes[1:]:
        node.text = ""


def set_neutral_metadata(doc: Document) -> None:
    props = doc.core_properties
    props.author = "1-office"
    props.last_modified_by = "1-office"
    props.title = ""
    props.subject = ""
    props.keywords = ""
    props.category = ""
    props.comments = ""
    props.content_status = ""
    props.identifier = ""
    props.version = ""
    props.language = ""
    props.revision = 1
    epoch = datetime(2000, 1, 1, tzinfo=timezone.utc)
    props.created = epoch
    props.modified = epoch
    props.last_printed = epoch


def edit_manuscript(src: Path, dst: Path) -> None:
    doc = Document(src)
    p = doc.paragraphs

    replace_text_preserve_drawings(
        p[39],
        "Beyond the benchmark, ResToken converted a natural-language prompt into "
        "position-9 analogues of the clinical-stage LUNA18-KRAS-G12D complex. Three "
        "analogues improved both Rosetta interface delta-delta-G and MOE interaction "
        "energy (best: -2.4 REU and -2.6 kcal/mol), and none appeared among the 45 "
        "compounds in the published campaign. A second analysis of cyclosporin A bound "
        "to cyclophilin A identified norleucine as the only dual-engine improver. These "
        "in silico candidates require permeability testing and synthesis; structures, "
        "scoring protocols, full results, and caveats are reported in SI Sections S7-S8 "
        "and Tables S9-S13.",
    )

    p[43].text = (
        "ResToken separates semantic instructions from chemical reconstruction. Across "
        "five models, its closed vocabulary prevented monomer hallucination, improved "
        "task-valid generation, and gave 100% frozen-position compliance. HELM-BERT "
        "learns predictive representations from explicit chemistry and topology, whereas "
        "ResToken provides an interpretable design interface for general-purpose LLMs; "
        "HELMify instead addresses monomer registration.17,18"
    )
    p[43].style = "Text-Indent"

    # Replacing the whole paragraph intentionally removes the two drawings that form
    # Figure 7. The complete case-study structures and tables remain in SI Sections S7-S8.
    p[44].text = (
        "Representation alone was insufficient for strict multi-constraint generation. "
        "Gemini 2.5 Pro achieved 80-100% constraint satisfaction with up to 141-fold "
        "enrichment, whereas local 9-12B models failed the conjunctive hard gate despite "
        "producing valid sequences. Nevertheless, an independent permeability predictor "
        "showed directional steering across four open-model families, indicating that "
        "soft property control generalizes more broadly than literal constraint satisfaction."
    )
    p[44].style = "Text-Indent"

    p[45].text = (
        "Limitations. The curated 400-block library is not exhaustive, gamma-backbone "
        "and N-methylated coverage remain limited, and categorical property bins trade "
        "resolution for usability. The representation comparison reflects inherent "
        "differences in vocabulary and metadata. ResToken is not a property predictor: "
        "it generates library-valid sequences under rule-based constraints but does not "
        "establish activity, permeability, or synthetic feasibility. The case-study "
        "analogues are computational hypotheses without wet-lab validation."
    )
    p[45].style = "Text-Indent"

    p[46].text = (
        "Future work will couple generation to independent permeability and structure "
        "predictors, expand validated NCAA coverage, and test the proposed analogues "
        "experimentally. The same tokenization principle may extend to stapled peptides "
        "and peptide-drug conjugates."
    )
    p[46].style = "Text-Indent"

    p[48].text = (
        "We compared ResToken, SMILES, and HELM across Gemini 2.5 Pro, Gemini 2.5 Flash, "
        "GPT-4o, Qwen 3.5 9B, and Gemma 3 12B. Tasks measured zero-shot validity (200 "
        "six-residue sequences per representation and model), positional control (10 "
        "parents, 20 variants each), and property-constrained generation (100 sequences "
        "per profile). A seven-rule validator checked token identity, length, charge, "
        "hydrogen-bond donors and acceptors, rotatable bonds, and backbone compatibility. "
        "Validity confidence intervals used 10,000 bootstrap resamples; representation "
        "comparisons used Pearson chi-square tests with Bonferroni correction. Random "
        "baselines used 100,000 uniformly sampled six-mers. Full prompts, constraint "
        "definitions, model settings, parsing rules, equations, statistical procedures, "
        "and comparison-fairness details are provided in SI Sections S1-S6."
    )
    p[48].style = "Text-Indent"

    p[70].text = (
        "The building-block library, dictionaries, benchmark outputs, assembled structures, "
        "and predictor outputs are included in the Supplementary Information and will be "
        "released publicly upon publication. The permeability model used CycPeptMPDB, and "
        "case-study structures are available from the Protein Data Bank under 7YV1 and 1CWA."
    )
    p[70].style = "Text-Indent"
    p[72].text = (
        "The ResToken package and experiment code will be released under the MIT license "
        "and archived with a citable DOI upon publication. No author-owned web resource is "
        "linked in this anonymous submission."
    )
    p[72].style = "Text-Indent"
    p[74].text = (
        "The Supplementary Information contains complete methods, prompts, statistical "
        "analyses, benchmark tables, case-study structures and scores, permeability "
        "validation, and cross-model generalization results."
    )
    p[74].style = "Text-Indent"

    # Remove the TOC graphic so every page after page 7 contains references only,
    # plus the second long case-study paragraph and the expanded Methods body.
    for idx in [99, 98, 68, 67, 66, 65, 64, 63, 62, 61, 60, 59, 58, 57, 56, 55,
                54, 53, 52, 51, 50, 49, 41, 40]:
        remove_paragraph(p[idx])

    set_neutral_metadata(doc)
    dst.parent.mkdir(parents=True, exist_ok=True)
    doc.save(dst)


def edit_si(src: Path, dst: Path) -> None:
    doc = Document(src)
    target = doc.paragraphs[86]
    target.text = (
        "All LLM evaluations used standardized prompt templates with the full building-"
        "block library injected at the {DICTIONARY_TABLE} placeholder. Complete templates "
        "for all five conditions (ResToken, SMILES, and HELM zero-shot generation; EDIT "
        "and FROZEN controllability; and property-constrained generation) are provided "
        "in Section S1 of this Supplementary Information."
    )
    set_neutral_metadata(doc)
    dst.parent.mkdir(parents=True, exist_ok=True)
    doc.save(dst)


def edit_checklist(src: Path, dst: Path) -> None:
    doc = Document(src)
    replacements = {
        19: "All novel datasets introduced in this paper are included in a data appendix "
            "(yes/partial/no/NA) Response: partial. The novel building-block library and "
            "benchmark summaries are provided in the Supplementary Information; complete "
            "machine-readable outputs will be released upon publication.",
        28: "Any code required for pre-processing data is included in the appendix "
            "(yes/partial/no) Response: partial. The preprocessing and reconstruction "
            "procedures are specified in the Supplementary Information; executable code "
            "will be released upon publication.",
        29: "All source code required for conducting and analyzing the experiments is "
            "included in a code appendix (yes/partial/no) Response: no. Complete methods, "
            "prompts, and statistical procedures are provided in the Supplementary "
            "Information, and the source code will be released upon publication.",
    }
    for idx, text in replacements.items():
        doc.paragraphs[idx].text = text
        doc.paragraphs[idx].style = "Checklist item"
    set_neutral_metadata(doc)
    dst.parent.mkdir(parents=True, exist_ok=True)
    doc.save(dst)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--source-dir", type=Path, required=True)
    ap.add_argument("--output-dir", type=Path, required=True)
    args = ap.parse_args()
    edit_manuscript(
        args.source_dir / "ResToken_manuscript_AAAI_v1.docx",
        args.output_dir / "ResToken_manuscript_AAAI_v2_compliant.docx",
    )
    edit_si(
        args.source_dir / "ResToken_SI.docx",
        args.output_dir / "ResToken_SI_v2_compliant.docx",
    )
    edit_checklist(
        args.source_dir / "ReproducibilityChecklist.docx",
        args.output_dir / "ReproducibilityChecklist_v2_compliant.docx",
    )


if __name__ == "__main__":
    main()
