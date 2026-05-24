# Antibiotic Cyclic Peptide

```json
{
  "modes": {
    "cationic_amphipathic": {
      "display_name": "Cationic Amphipathic Antibiotic",
      "ring_sizes": [7, 8, 9, 10],
      "charge_range": [2, 6],
      "recipes": {
        "7": {"primary_amine_cation": [2, 3], "bulky_aliphatic": [2, 3]},
        "8": {"primary_amine_cation": [3, 4], "bulky_aliphatic": [2, 3]},
        "9": {"primary_amine_cation": [3, 5], "bulky_aliphatic": [2, 4]},
        "10": {"primary_amine_cation": [3, 5], "bulky_aliphatic": [3, 4]}
      },
      "filler_bins": ["D_residue", "turn_inducer", "aromatic_hydrophobe"],
      "hard_filters": {
        "max_ring_size": 12,
        "min_cationic": 2,
        "min_charge": 2,
        "max_acidic": 1,
        "max_aromatic": 3
      },
      "priority_thresholds": [10, 7]
    },
    "broad_spectrum": {
      "display_name": "Broad-Spectrum Amphipathic",
      "ring_sizes": [8, 9, 10],
      "charge_range": [2, 5],
      "recipes": {
        "8": {"all_cationic": [2, 3], "aromatic_hydrophobe": [1, 2], "bulky_aliphatic": [1, 2]},
        "9": {"all_cationic": [2, 4], "aromatic_hydrophobe": [1, 2], "bulky_aliphatic": [2, 3]},
        "10": {"all_cationic": [3, 4], "aromatic_hydrophobe": [1, 2], "bulky_aliphatic": [2, 3]}
      },
      "filler_bins": ["D_residue", "turn_inducer", "neutral_polar"],
      "hard_filters": {
        "max_ring_size": 12,
        "min_cationic": 2,
        "min_charge": 2,
        "max_acidic": 2
      },
      "priority_thresholds": [10, 7]
    }
  }
}
```

## Design Context

### Cationic amphipathic mode
Polymyxin/battacin-inspired: Dab/Orn/Lys-rich cationic ring + hydrophobic residues for membrane insertion.
Best starting design: 8-10 residues, +3 to +5 charge, 2-4 bulky aliphatic hydrophobes, 1-3 D-amino acids.
Prioritize primary amine cations (Lys/Orn/Dab-like) over Arg for antibiotic selectivity.
Note: lipid tails (C8-C12) are not available in the 400-block library; hydrophobic residues substitute.

### Broad-spectrum mode
General amphipathic antimicrobial peptide: mixed cationic + aromatic + aliphatic hydrophobes.
Best starting design: 8-10 residues, +2 to +5 charge, 35-60% hydrophobic fraction, 1-2 aromatics, 1-3 D/Pro turn residues.
Allows Arg alongside Lys/Orn for broader cationic diversity.

### Block bins reference
- **primary_amine_cation** (18): Lys/Orn/Dab-like, charge +1, no guanidinium
- **all_cationic** (21): all charge +1 blocks including Arg
- **guanidinium_cation** (3): Arg-like with guanidinium group
- **bulky_aliphatic** (106): large neutral aliphatics (Leu/Ile/Val/Cha-like)
- **aromatic_hydrophobe** (78): Phe/Trp/Tyr-like, neutral
- **D_residue** (86): D-chirality blocks
- **turn_inducer** (125): D-residues + proline-like + glycine
- **acidic_anion** (10): charge -1 (Asp/Glu-like)
- **halogenated** (5): halogen-containing blocks
- **neutral_polar** (155): neutral, has HBD or HBA
