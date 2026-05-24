# Design Guidelines for Cell-Penetrating Cyclic Peptides

**Purpose:** Provide residue-composition and structural-feature rules for a cyclic peptide generation model with a ~400 non-canonical amino acid (NCAA) building-block library.

**Scope:** These rules are optimized for generating **cell-penetrating cyclic peptides**. They distinguish two permeability mechanisms:

1. **CPP-like endocytic uptake + endosomal escape**: cationic, amphipathic, Arg-rich cyclic peptides.
2. **Passive membrane diffusion**: more hydrophobic/chameleonic cyclic peptides with reduced exposed polarity.

For intracellular target binders, use CPP-like rules when uptake is the primary goal, and passive-diffusion rules when drug-like membrane permeability/oral-like behavior is the primary goal.

---

## 1. Core Design Decision

Before generation, assign each candidate to one permeability mode.

| Mode | Best use case | Ring size | Charge | Key structural logic |
|---|---:|---:|---:|---|
| **CPP-like cyclic peptide** | Cytosolic delivery, endosomal escape, intracellular cargo delivery | **6–9 residues** | **+3 to +5**, usually Arg-rich | Constrained amphipathic ring with cationic face + hydrophobic/aromatic face |
| **Passive permeable cyclic peptide** | Drug-like permeability, lower charge, oral-like optimization | **5–10 residues**, occasionally 11 | **−1 to +1**, ideally neutral | Low exposed HBD/HBA, N-methylation, intramolecular H-bonds, steric shielding |
| **Bicyclic CPP-binder** | Target binding needs larger surface area | CPP ring **6–8** + binding ring variable | CPP ring +3 to +5 | One ring for uptake, second ring for target recognition |

**Default recommendation for the generation model:** start from **CPP-like 7–8 residue monocyclic peptides** and **bicyclic designs with a 6–8 residue CPP ring**.

---

## 2. CPP-like Cell-Penetrating Cyclic Peptide Rules

### 2.1 Ring size

Use the following ring-size priority:

```text
Priority 1: 6–8 residues
Priority 2: 9 residues
Priority 3: 10–11 residues only if target-binding motif requires it
Avoid: >11 residues for pure CPP function
```

**Reasoning:** cyclohexa- to cyclooctapeptides are the best-established cyclic CPP size window. Increasing ring size progressively reduces cellular entry efficiency. Rings ≤9 residues are repeatedly associated with high cyclic CPP activity.

---

### 2.2 Positive charged residues

Use **Arg or Arg-like guanidinium residues** as the main cationic elements.

| Residue class | Use | Recommendation |
|---|---|---|
| **L-Arg** | Primary cationic residue | Use frequently |
| **D-Arg** | Protease-resistant cationic residue | Strongly recommended |
| **homoArg / guanidinylated Lys analogs** | Arg mimics | Good NCAA options |
| **Lys / Orn / Dab / Dap** | Cationic but weaker membrane interaction than Arg | Use sparingly unless target binding requires |
| **His / methyl-His / aromatic cationic residues** | pH-sensitive/endosomal effect | Useful as minor support, not Arg replacement |

**Recommended number of positive residues:**

```text
6-mer: 3 Arg-like residues
7-mer: 3–4 Arg-like residues
8-mer: 3–4 Arg-like residues
9-mer: 4 Arg-like residues
>9-mer: 4–5 Arg-like residues only if solubility/endocytosis is desired
```

**Do not use cyclic R6 or all-cationic rings as default.** Cyclic R6 has weak CPP activity compared with amphipathic rings containing both Arg and hydrophobic residues.

---

### 2.3 Where to place the positive residues

The positively charged residues should form a **cationic face or cationic arc**, not a uniformly distributed charge cloud.

Preferred 3D arrangement:

```text
One face: Arg / D-Arg / guanidinium residues
Opposite or adjacent face: aromatic hydrophobic residues
Remaining positions: turn, linker, target-binding, or polarity-tuning residues
```

Sequence-level patterns that often help generate this 3D distribution:

```text
Pattern A: cyclo[Ar-Hyd-Arg-Arg-Arg-Arg-Polar]
Pattern B: cyclo[Ar-Ar-Arg-Hyd-Arg-Hyd-Polar]
Pattern C: cyclo[Ar-DAr-Arg-DArg-Arg-DArg-Polar]
Pattern D: cyclo[Ar-Hyd-Arg-X-Arg-Hyd-Arg-X]
```

Where:

```text
Ar = aromatic hydrophobe, preferably Nal / Trp / Phe / halogenated Phe / biphenyl-Ala
Hyd = hydrophobic residue, e.g., Leu / Ile / Val / Cha / tBuGly / cycloalkyl-Ala
Arg = Arg-like guanidinium residue
DArg = D-Arg or D-guanidinium analog
Polar = Gln / Ser / Thr / Tyr / weakly polar NCAA / linker residue
X = target-binding residue or conformational control residue
```

**Practical rule:** In generated 3D structures, at least **3 guanidinium groups** should be solvent-exposed and spatially close enough to create a membrane-binding patch. They should not be buried inside the ring.

---

### 2.4 Hydrophobic/aromatic residues

CPP-like cyclic peptides need **1–2 hydrophobic residues**, preferably **aromatic hydrophobes**.

| Class | Examples | Role |
|---|---|---|
| Aromatic hydrophobes | Phe, D-Phe, Trp, Nal, Bip, halogenated Phe | Membrane insertion, lipid interaction, endosomal escape |
| Bulky aliphatic hydrophobes | Leu, Ile, Val, Cha, tBuGly, cyclohexyl-Ala | Hydrophobic face, steric shielding |
| Too hydrophobic | Multiple large aromatics + long alkyl chains | Risk of aggregation, membrane toxicity, poor solubility |

**Recommended count:**

```text
6–7 residues: 1–2 hydrophobic/aromatic residues
8–9 residues: 2 hydrophobic/aromatic residues
>9 residues: 2–3 hydrophobic/aromatic residues, but monitor aggregation
```

**Strong starting motif:**

```text
cyclo[Phe/Nal - Arg - Arg - Arg - Arg - Gln/Ser]
```

or for better stability:

```text
cyclo[D-Phe/Nal - Arg - D-Arg - Arg - D-Arg - Gln/Ser]
```

---

### 2.5 Amphipathicity requirement

High-priority candidates should have:

```text
Cationic face: 3–4 Arg-like side chains
Hydrophobic face: 1–2 aromatic or bulky hydrophobic side chains
Neutral/polar tuning site: 1 residue for solubility or target binding
```

Avoid these patterns:

```text
All positive: Arg-Arg-Arg-Arg-Arg-Arg
All hydrophobic: Leu/Phe/Nal-rich with no cationic face
Random charge distribution with no amphipathic topology
Too many acidic residues near Arg residues
```

---

## 3. Passive-Permeable Cyclic Peptide Rules

These rules are different from CPP-like uptake. Passive permeability does **not** mainly come from high positive charge.

### 3.1 Ring size

```text
Best: 5–7 residues
Useful: 8–10 residues
Special cases: 11 residues, e.g., cyclosporine-like chameleonic peptides
Avoid: >11 residues unless strong conformational shielding is designed
```

---

### 3.2 Composition

Passive-permeable candidates should favor:

```text
Hydrophobic residues: 40–70%
Aromatic residues: 1–2
Backbone N-methylated residues: 1–4 depending on ring size
D-residues / Pro / Pip / Aib / beta-turn residues: 1–3
Net charge: −1 to +1, ideally 0
Acidic/basic residues: minimized unless target binding requires
```

Useful NCAA classes:

| NCAA class | Why useful |
|---|---|
| **N-methyl amino acids** | Remove backbone HBD, reduce desolvation penalty |
| **D-amino acids** | Protease resistance, turn induction, conformational control |
| **Pro / D-Pro / Pip / pipecolic acid analogs** | Turn control, conformational restriction |
| **Aib / tBuGly / cycloalkyl Gly** | Steric shielding of backbone amides |
| **β- or γ-amino acid units / statine-like residues** | Can promote intramolecular H-bonding |
| **Peptoid-like residues** | Remove backbone HBD and add side-chain diversity |

---

### 3.3 Polarity shielding

Passive-permeable designs should minimize exposed polar surface in membrane-like conformations.

Prioritize candidates that satisfy:

```text
Backbone HBD count exposed to solvent: low
Intramolecular H-bonds: ≥2 for 6–10 residue rings
Sterically shielded backbone amides: present
Conformational switch possible: open in water, closed in membrane-like environment
```

Good features:

```text
N-methylation at solvent-exposed amides
D/L alternation to create turns
Bulky hydrophobes projecting over backbone amides
Intramolecular side-chain-to-backbone H-bonds
Transannular backbone H-bonds
```

Avoid:

```text
Many exposed amide NH groups
Multiple unpaired charged side chains
High flexibility without stable low-polar conformer
Large ring size without intramolecular H-bond network
```

---

## 4. How to Use the 400 NCAA Building Blocks

Classify every NCAA into functional bins before model generation.

### 4.1 Required building-block bins

Each NCAA should be annotated with:

```yaml
charge_class:
  - guanidinium_cation
  - primary_amine_cation
  - imidazole_weak_cation
  - acidic_anion
  - neutral_polar
  - neutral_hydrophobic

hydrophobicity_class:
  - aromatic_hydrophobe
  - bulky_aliphatic
  - small_aliphatic
  - polar
  - charged

backbone_feature:
  - normal_alpha
  - D_alpha
  - N_methyl_alpha
  - beta_amino_acid
  - gamma_amino_acid
  - peptoid_like
  - proline_like
  - turn_inducer

permeability_role:
  - cationic_membrane_binding
  - aromatic_membrane_insertion
  - hydrophobic_shielding
  - intramolecular_Hbond_promoter
  - turn_or_rigidity
  - solubility_tuning
  - target_binding
  - avoid_or_penalize
```

---

### 4.2 NCAA selection logic for CPP-like cyclic peptides

For each generated CPP-like cyclic peptide:

```text
Choose 3–4 from: Arg / D-Arg / guanidinium NCAA
Choose 1–2 from: Nal / Trp / Phe / D-Phe / aromatic hydrophobic NCAA
Choose 0–2 from: Leu / Cha / tBuGly / bulky hydrophobic NCAA
Choose 1 from: Gln / Ser / Thr / Tyr / neutral polar / target-binding residue
Choose 1–2 from: D-residue / Pro-like / Pip / Aib if conformational control is needed
Avoid >1 acidic residue unless required for target binding
Avoid >3 large aromatic residues unless aggregation filters pass
```

**Default generator recipe:**

```text
Ring size: 7 or 8
Arg-like residues: 4
Aromatic hydrophobes: 2
Neutral/tuning residue: 1–2
D-residue fraction: 25–50%
Net charge: +3 to +5
```

---

### 4.3 NCAA selection logic for passive-permeable cyclic peptides

For each generated passive-permeable cyclic peptide:

```text
Choose 2–4 from: N-methyl AA / peptoid-like AA / ester-like surrogate
Choose 2–5 from: hydrophobic AA / aromatic hydrophobic AA
Choose 1–3 from: D-AA / Pro-like / Aib / turn-inducing NCAA
Choose 0–2 from: neutral polar AA capable of intramolecular H-bonding
Avoid strong cationic or acidic residues unless target binding requires
Net charge: −1 to +1
```

**Default generator recipe:**

```text
Ring size: 6–9
N-methyl or peptoid-like positions: 2–3
Hydrophobic/aromatic residues: 3–5
D/Pro-like turn residues: 1–2
Net charge: 0 or +1
Predicted exposed HBD/HBA: low
```

---

## 5. Positive Residue Placement Rules

### 5.1 CPP-like mode

Positive residues should be **surface-exposed** and arranged as a **patch**.

Use these layouts:

```text
Layout 1: contiguous cationic arc
cyclo[Ar-Hyd-Arg-Arg-Arg-Arg-Polar]

Layout 2: alternating cationic-hydrophobic amphipathic ring
cyclo[Trp-Arg-Trp-Arg-Trp-Arg-Trp-Arg]

Layout 3: mixed stereochemistry cationic face
cyclo[D-Phe-Nal-Arg-DArg-Arg-DArg-Gln]

Layout 4: target-binding + CPP hybrid
cyclo[Target1-Target2-Nal-Arg-Arg-DArg-DArg]
```

Avoid:

```text
Buried Arg side chains
Arg paired tightly with acidic residues inside the ring
Uniformly distributed positive residues without hydrophobic face
Long flexible cationic side chains all pointing in random directions
```

---

### 5.2 Passive-permeable mode

Positive residues should be minimized.

Use:

```text
0 Arg/Lys: preferred for passive diffusion
1 Arg/Lys: acceptable if target binding requires it
2 Arg/Lys: only if compensated by N-methylation, hydrophobic shielding, and intramolecular H-bonds
≥3 Arg/Lys: switch design label to CPP-like uptake, not passive diffusion
```

---

## 6. Scoring Rules for Model Filtering

### 6.1 CPP-like permeability score

Use a simple first-pass score:

```text
CPP_score =
  +2.0 × number_of_Arg_like_residues_clipped_to_4
  +2.0 × number_of_aromatic_hydrophobes_clipped_to_2
  +1.5 × amphipathic_face_score
  +1.0 × ring_constraint_score
  +1.0 × D_residue_or_NCAA_stability_score
  −2.0 × number_of_acidic_residues
  −2.0 × aggregation_penalty
  −2.0 × excessive_flexibility_penalty
```

Recommended acceptance:

```text
CPP_score ≥ 10: strong CPP-like candidate
CPP_score 7–9: acceptable, keep for diversity
CPP_score <7: deprioritize
```

---

### 6.2 Passive permeability score

```text
Passive_score =
  +2.0 × intramolecular_Hbond_count
  +2.0 × N_methyl_or_peptoid_count
  +1.5 × steric_backbone_shielding_score
  +1.0 × hydrophobic_balance_score
  +1.0 × conformational_preorganization_score
  −2.0 × exposed_backbone_HBD_penalty
  −2.0 × net_charge_penalty
  −2.0 × large_ring_penalty
  −2.0 × aggregation_penalty
```

Recommended acceptance:

```text
Passive_score ≥ 9: strong passive-permeability candidate
Passive_score 6–8: acceptable if target binding is strong
Passive_score <6: deprioritize
```

---

## 7. Hard Filters

### 7.1 CPP-like cyclic peptides

Reject or strongly penalize:

```text
Ring size >11 for pure CPP design
Arg-like count <3
No aromatic hydrophobic residue
Net charge <+2
>2 acidic residues
Hydrophobic/aromatic count >4 in ≤8 residue ring
Predicted severe aggregation
Predicted unstable/unstructured ensemble with no amphipathic face
```

### 7.2 Passive-permeable cyclic peptides

Reject or strongly penalize:

```text
Ring size >11
Net charge >+2 or <−2
Exposed backbone NH count high
No N-methyl / no intramolecular H-bond / no steric shielding
Too many rotatable side chains
High polar surface without chameleonic shielding
```

---

## 8. Recommended Starting Libraries

### Library A: Minimal CPP-like rings

```text
Ring size: 6–7
Composition: 3–4 Arg-like + 1–2 aromatic hydrophobe + 1 neutral polar/turn residue
Goal: high uptake, compact scaffold
```

Example templates:

```text
cyclo[Nal-Arg-DArg-Arg-DArg-Gln]
cyclo[Phe-Nal-Arg-DArg-Arg-DArg-Gln]
cyclo[Trp-Arg-Trp-Arg-Arg-DArg-Ser]
```

---

### Library B: Balanced CPP-binder rings

```text
Ring size: 8–10
Composition: 3–4 Arg-like + 2 hydrophobes + 2–3 target-binding residues
Goal: combine cell penetration and target recognition
```

Example templates:

```text
cyclo[Target1-Nal-Arg-DArg-Target2-Arg-Hyd-DArg]
cyclo[Phe-Target1-Arg-DArg-Nal-Target2-Arg-Gln]
```

---

### Library C: Passive-permeable rings

```text
Ring size: 6–9
Composition: 2–3 N-methyl/peptoid-like residues + 3–5 hydrophobes + 1–2 turn residues
Goal: low exposed polarity and passive diffusion
```

Example templates:

```text
cyclo[NMeLeu-DPro-Leu-NMePhe-Val-Aib]
cyclo[NMeVal-Leu-DPhe-NMeLeu-Pro-Tyr]
cyclo[Leu-NMeLeu-DLeu-Pro-NMePhe-tBuGly]
```

---

### Library D: Bicyclic CPP-delivery designs

```text
CPP ring: 6–8 residues
Binding ring: target-dependent
Composition of CPP ring: 3–4 Arg-like + 1–2 aromatic hydrophobes
Goal: preserve CPP ring while allowing larger target-binding surface
```

Example architecture:

```text
Ring 1, CPP ring:
cyclo[Phe-Nal-Arg-DArg-Arg-DArg]

Ring 2, binder ring:
target-specific residues selected by docking/design model
```

---

## 9. Model Output Annotation

For every generated peptide, report:

```yaml
sequence:
ring_size:
cyclization_type:
net_charge_pH7_4:
Arg_like_count:
D_Arg_like_count:
aromatic_hydrophobe_count:
hydrophobic_fraction:
N_methyl_count:
D_residue_count:
Pro_or_turn_residue_count:
acidic_residue_count:
predicted_amphipathic_face_score:
predicted_intramolecular_Hbond_count:
predicted_exposed_backbone_HBD_count:
predicted_SASA_polar:
predicted_aggregation_risk:
mode_label:
  - CPP_like_endocytic
  - passive_diffusion
  - bicyclic_CPP_binder
priority:
  - high
  - medium
  - low
reason_for_priority:
```

---

## 10. Final Practical Rules

### For CPP-like cyclic peptide generation

Use this as the main rule:

```text
Small constrained ring + 3–4 Arg-like residues + 1–2 aromatic hydrophobes + amphipathic 3D display.
```

Best starting design:

```text
7–8 residues
+3 to +5 charge
3–4 Arg/D-Arg/guanidinium NCAA
1–2 Nal/Trp/Phe-like NCAA
1–2 D/Pro/Aib/Pip-like conformational control residues
0–1 neutral polar residue
```

---

### For passive-permeable cyclic peptide generation

Use this as the main rule:

```text
Small hydrophobic/chameleonic ring + low exposed polarity + N-methylation/peptoid residues + intramolecular H-bonding.
```

Best starting design:

```text
6–9 residues
net charge −1 to +1
2–3 N-methyl/peptoid-like residues
3–5 hydrophobic residues
1–2 D/Pro/Aib-like turn residues
≥2 intramolecular H-bonds or clear backbone shielding
```

---

## 11. Design Prioritization

For your generation model, rank candidates in this order:

1. **8-residue CPP-like amphipathic rings** with 4 Arg-like residues and 2 aromatic hydrophobes.
2. **7-residue CPP-like rings** with 3–4 Arg-like residues, 1–2 aromatic hydrophobes, and mixed D/L stereochemistry.
3. **Bicyclic CPP-binder designs** where the CPP ring is preserved and target-binding residues are placed in the second ring.
4. **Passive-permeable 6–9 residue rings** with N-methylation, hydrophobic shielding, and low exposed polarity.
5. Larger macrocycles only if target binding requires them and permeability is handled by a CPP ring or strong chameleonic shielding.

---

## 12. Key References Used

- Dougherty, P. G.; Sahni, A.; Pei, D. **Understanding Cell Penetration of Cyclic Peptides.** *Chemical Reviews* 2019, 119, 10241–10287.
- Rettie, S. A. et al. **Accurate de novo design of high-affinity protein-binding macrocycles using deep learning.** *Nature Chemical Biology* 2025.
- Xie, X. et al. **CyclicBoltz1: fast and accurately predicting structures of cyclic peptides and complexes containing non-canonical amino acids using AlphaFold 3 framework.** bioRxiv 2025.
- Cao, S. et al. **Accurate structure prediction of cyclic peptides containing unnatural amino acids using HighFold3.** *Briefings in Bioinformatics* 2025.
