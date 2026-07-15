#!/usr/bin/env python3
"""Convert the ResToken manuscript into the AAAI anonymous Word layout."""

from __future__ import annotations

import argparse
import copy
import re
import zipfile
from pathlib import Path
from typing import cast

from docx import Document
from docx.document import Document as DocumentType
from docx.enum.section import WD_SECTION
from docx.enum.text import WD_ALIGN_PARAGRAPH, WD_BREAK, WD_LINE_SPACING
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
from docx.shared import Inches, Pt
from docx.styles.style import _ParagraphStyle
from docxcompose.composer import Composer


W = "http://schemas.openxmlformats.org/wordprocessingml/2006/main"


def clear_body(doc: DocumentType) -> None:
    body = doc._element.body
    for child in list(body):
        if child.tag != qn("w:sectPr"):
            body.remove(child)


def set_cols(sectpr, count: int) -> None:
    for old in sectpr.xpath("./w:cols"):
        sectpr.remove(old)
    cols = OxmlElement("w:cols")
    cols.set(qn("w:num"), str(count))
    cols.set(qn("w:space"), "540")  # 0.375 inch
    sectpr.append(cols)


def set_page_geometry(section, first: bool) -> None:
    section.page_width = Inches(8.5)
    section.page_height = Inches(11)
    section.left_margin = Inches(0.75)
    section.right_margin = Inches(0.75)
    section.top_margin = Inches(1.25 if first else 0.75)
    section.bottom_margin = Inches(1.25)
    section.header_distance = Inches(0.3)
    section.footer_distance = Inches(0.3)


def set_style(doc: DocumentType, name: str, size: float = 10, bold=False,
              italic=False, align=None, before=0, after=0, line=12,
              first_indent=None, keep_with_next=False) -> None:
    style = cast(_ParagraphStyle, doc.styles[name])
    style.font.name = "Times New Roman"
    getattr(style, "_element").get_or_add_rPr().get_or_add_rFonts().set(
        qn("w:eastAsia"), "Times New Roman"
    )
    style.font.size = Pt(size)
    style.font.bold = bold
    style.font.italic = italic
    pf = style.paragraph_format
    pf.space_before = Pt(before)
    pf.space_after = Pt(after)
    pf.line_spacing_rule = WD_LINE_SPACING.EXACTLY
    pf.line_spacing = Pt(line)
    pf.keep_with_next = keep_with_next
    if align is not None:
        pf.alignment = align
    if first_indent is not None:
        pf.first_line_indent = Inches(first_indent)


def configure_styles(doc: DocumentType) -> None:
    set_style(doc, "Paper Title", 15, bold=True, align=WD_ALIGN_PARAGRAPH.CENTER,
              after=10, line=17, keep_with_next=True)
    set_style(doc, "Author Name", 12, align=WD_ALIGN_PARAGRAPH.CENTER,
              after=4, line=14, keep_with_next=True)
    set_style(doc, "Affiliation and Address", 10, align=WD_ALIGN_PARAGRAPH.CENTER,
              after=10, line=12, keep_with_next=True)
    set_style(doc, "Abstract Head", 10, bold=True, align=WD_ALIGN_PARAGRAPH.CENTER,
              before=3, after=3, keep_with_next=True)
    set_style(doc, "Abstract Text", 9, line=11, after=4)
    set_style(doc, "Text", 10, line=12, after=0)
    set_style(doc, "Text-Indent", 10, line=12, after=0, first_indent=0.15)
    set_style(doc, "Section Heading", 11, bold=True,
              align=WD_ALIGN_PARAGRAPH.CENTER, before=9, after=4,
              line=13, keep_with_next=True)
    set_style(doc, "Subsection Heading", 10, bold=True, before=6, after=0,
              line=12, keep_with_next=True)
    set_style(doc, "Figure Caption", 9, line=10, after=5)
    set_style(doc, "References", 9, line=10, after=1)
    set_style(doc, "External Paper Links", 9, line=10, after=0)
    set_style(doc, "Normal", 10, line=12, after=0)


def strip_headers_footers(doc: DocumentType) -> None:
    for section in doc.sections:
        section.header.is_linked_to_previous = False
        section.footer.is_linked_to_previous = False
        for container in (section.header, section.footer):
            for p in container.paragraphs:
                p._element.getparent().remove(p._element)
        for ref in section._sectPr.xpath("./w:headerReference|./w:footerReference"):
            section._sectPr.remove(ref)


def flatten_links_and_bookmarks(doc: DocumentType) -> None:
    """AAAI forbids embedded links/bookmarks; retain visible black URL text."""
    root = doc._element
    for hyperlink in list(root.xpath(".//w:hyperlink")):
        parent = hyperlink.getparent()
        pos = parent.index(hyperlink)
        for child in list(hyperlink):
            hyperlink.remove(child)
            parent.insert(pos, child)
            pos += 1
        parent.remove(hyperlink)
    for bookmark in list(root.xpath(".//w:bookmarkStart|.//w:bookmarkEnd")):
        bookmark.getparent().remove(bookmark)
    for rid, rel in list(doc.part.rels.items()):
        if rel.reltype.endswith("/hyperlink"):
            doc.part.rels.pop(rid)


def normalize_runs(paragraph) -> None:
    for run in paragraph.runs:
        run.font.name = "Times New Roman"
        run._element.get_or_add_rPr().rFonts.set(qn("w:eastAsia"), "Times New Roman")
        run.font.color.rgb = None


def classify_paragraph(p, idx: int) -> str:
    text = p.text.strip()
    old = p.style.name
    if idx == 0:
        return "Paper Title"
    if idx == 1:
        return "Author Name"
    if idx == 2:
        return "Affiliation and Address"
    if text == "ABSTRACT":
        return "Abstract Head"
    if old == "BD_Abstract":
        return "Abstract Text"
    if old in {"VA_Figure_Caption", "VD_Table_Title"} or re.match(r"^(Figure|Table)\s+\d+\.", text):
        return "Figure Caption"
    if old == "TF_References_Section":
        return "Section Heading" if text.upper() == "REFERENCES" else "References"
    if text in {"METHODS", "Results", "Discussion", "Methods", "Data availability",
                "Code availability", "Supplementary Information"}:
        return "Section Heading"
    if old == "TE_Supporting_Information":
        return "Text"
    # Source uses bold leading runs for subsection heads.
    if text and len(text) < 130 and (text.endswith(".") or text.endswith(":")):
        if p.runs and p.runs[0].bold:
            return "Subsection Heading"
    return "Text-Indent"


def add_first_section_break(after_p, doc: DocumentType) -> None:
    p = OxmlElement("w:p")
    ppr = OxmlElement("w:pPr")
    sectpr = copy.deepcopy(doc._element.body.sectPr)
    set_cols(sectpr, 1)
    pgmar = sectpr.find(qn("w:pgMar"))
    pgmar.set(qn("w:top"), "1800")
    ppr.append(sectpr)
    p.append(ppr)
    after_p._p.addnext(p)


def section_break_paragraph(reference_sectpr, columns: int):
    p = OxmlElement("w:p")
    ppr = OxmlElement("w:pPr")
    sectpr = copy.deepcopy(reference_sectpr)
    set_cols(sectpr, columns)
    pgmar = sectpr.find(qn("w:pgMar"))
    pgmar.set(qn("w:top"), "1080")
    ppr.append(sectpr)
    p.append(ppr)
    return p


def add_full_width_object_sections(doc: DocumentType) -> None:
    """Give figures and tables full text width inside the two-column paper."""
    body = doc._element.body
    reference = body.sectPr
    children = list(body)
    i = 0
    while i < len(children):
        child = children[i]
        is_drawing = child.tag == qn("w:p") and bool(child.xpath(".//w:drawing"))
        is_table_caption = False
        if child.tag == qn("w:p"):
            text = "".join(child.itertext()).strip()
            is_table_caption = bool(re.match(r"^Table\s+\d+\.", text))
        if is_drawing:
            start = child
            end = child
            if i + 1 < len(children) and children[i + 1].tag == qn("w:p"):
                nxt_text = "".join(children[i + 1].itertext()).strip()
                if re.match(r"^Figure\s+\d+\.", nxt_text):
                    end = children[i + 1]
                    i += 1
            start.addprevious(section_break_paragraph(reference, 2))
            end.addnext(section_break_paragraph(reference, 1))
        elif is_table_caption and i + 1 < len(children) and children[i + 1].tag == qn("w:tbl"):
            child.addprevious(section_break_paragraph(reference, 2))
            children[i + 1].addnext(section_break_paragraph(reference, 1))
            i += 1
        i += 1


def scrub_package_metadata(path: Path) -> None:
    """Force anonymous OOXML metadata after libraries finish saving."""
    with zipfile.ZipFile(path, "r") as zin:
        items = {name: zin.read(name) for name in zin.namelist()}
    core = items.get("docProps/core.xml", b"")
    core = re.sub(rb"(<dc:creator[^>]*>).*?(</dc:creator>)", rb"\1Anonymous\2", core)
    core = re.sub(
        rb"(<cp:lastModifiedBy[^>]*>).*?(</cp:lastModifiedBy>)",
        rb"\1Anonymous\2",
        core,
    )
    items["docProps/core.xml"] = core
    tmp = path.with_suffix(".metadata.docx")
    with zipfile.ZipFile(tmp, "w", zipfile.ZIP_DEFLATED) as zout:
        for name, payload in items.items():
            zout.writestr(name, payload)
    tmp.replace(path)


def remove_metadata(doc: DocumentType) -> None:
    cp = doc.core_properties
    cp.author = "Anonymous"
    cp.last_modified_by = "Anonymous"
    cp.title = "ResToken: Semantic Tokens for LLM-Based Cyclic Peptide Design"
    cp.subject = "Anonymous submission"
    cp.comments = ""
    cp.keywords = ""
    cp.category = ""


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("source", type=Path)
    ap.add_argument("template", type=Path)
    ap.add_argument("output", type=Path)
    args = ap.parse_args()

    master = Document(args.template)
    clear_body(master)
    Composer(master).append(Document(args.source))
    master.save(args.output)

    doc = Document(args.output)
    configure_styles(doc)
    remove_metadata(doc)

    # Anonymous title block required by the official template.
    doc.paragraphs[1].text = "Anonymous Submission"
    doc.paragraphs[2].text = "Anonymous affiliation"

    for idx, p in enumerate(doc.paragraphs):
        if not p.text.strip() and not p._p.xpath(".//w:drawing|.//m:oMath|.//m:oMathPara"):
            continue
        p.style = classify_paragraph(p, idx)
        normalize_runs(p)
        p.paragraph_format.widow_control = True

    # Title block is full width; abstract and body are two columns.
    set_cols(doc._element.body.sectPr, 2)
    add_first_section_break(doc.paragraphs[2], doc)
    add_full_width_object_sections(doc)
    for i, section in enumerate(doc.sections):
        set_page_geometry(section, i == 0)

    # Remove page numbers, headers, footers, embedded links/bookmarks, and personal metadata.
    strip_headers_footers(doc)
    flatten_links_and_bookmarks(doc)
    for rpr in doc._element.xpath(".//w:rPr"):
        for child in list(rpr):
            if child.tag in {qn("w:color"), qn("w:u")}: rpr.remove(child)

    # Prevent figure/table rows from splitting and keep captions with their objects.
    for table in doc.tables:
        table.autofit = True
        for row in table.rows:
            trpr = row._tr.get_or_add_trPr()
            trpr.append(OxmlElement("w:cantSplit"))
            for cell in row.cells:
                for p in cell.paragraphs:
                    normalize_runs(p)
                    p.paragraph_format.space_after = Pt(0)
                    p.paragraph_format.line_spacing_rule = WD_LINE_SPACING.EXACTLY
                    p.paragraph_format.line_spacing = Pt(10)

    doc.settings.element.append(OxmlElement("w:doNotTrackMoves"))
    doc.save(args.output)
    scrub_package_metadata(args.output)

    # Parseability gate.
    check = Document(args.output)
    assert check.paragraphs[0].text.startswith("ResToken:")
    assert check.paragraphs[1].text == "Anonymous Submission"
    assert "Yike Zou" not in "\n".join(p.text for p in check.paragraphs)
    assert len(check.inline_shapes) >= 1


if __name__ == "__main__":
    main()
