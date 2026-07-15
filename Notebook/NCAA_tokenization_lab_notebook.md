# NCAA Tokenization Lab Notebook

### 2026-07-15 | aaai_word_anonymous_v2

- Objective: convert the live ResToken manuscript to the official AAAI 2027 anonymous Word layout.
- Parameters: US Letter; 0.75-inch side margins; 1.25/0.75-inch first/subsequent top margins; 2 columns with 0.375-inch gap; 10-point Times New Roman body with exact 12-point leading.
- Result: preserved 9 drawings, 3 tables, and 9 OMML equation nodes; removed author/affiliation/ORCID, personal metadata, headers, footers, page numbers, embedded hyperlinks, and bookmarks.
- Conclusion: OOXML package is parseable and ZIP-valid; wide figures/tables are placed in full-width continuous sections.

### 2026-07-15 | aaai_word_v3_figures_fixed

- Objective: correct the v2 display failure in which only one figure appeared in Word/Drive preview.
- Diagnosis: all 9 media hashes and relationships were preserved, but 22 alternating continuous section breaks made rendering unstable.
- Parameters: removed figure-specific section breaks; retained two document sections; scaled all 9 inline figures proportionally to 3.12-inch column width.
- Result: 9/9 distinct image relationships resolve to nonempty media; 9/9 figure widths are 3.12 inches; 3 tables and 9 OMML nodes retained.
- Conclusion: v3 uses stable AAAI two-column flow; relationship-level image validation is now a mandatory build gate.
