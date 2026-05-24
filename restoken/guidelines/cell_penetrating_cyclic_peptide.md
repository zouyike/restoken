# Cell-Penetrating Cyclic Peptide

```json
{
  "modes": {
    "cpp_like": {
      "display_name": "CPP-like Cell-Penetrating",
      "ring_sizes": [6, 7, 8, 9],
      "charge_range": [2, 5],
      "recipes": {
        "6": {"all_cationic": [3, 3], "aromatic_hydrophobe": [1, 2]},
        "7": {"all_cationic": [3, 4], "aromatic_hydrophobe": [1, 2]},
        "8": {"all_cationic": [3, 4], "aromatic_hydrophobe": [2, 2]},
        "9": {"all_cationic": [4, 4], "aromatic_hydrophobe": [2, 2]}
      },
      "filler_bins": ["neutral_polar", "turn_inducer", "neutral_hydrophobic"],
      "hard_filters": {
        "max_ring_size": 11,
        "min_cationic": 3,
        "min_aromatic": 1,
        "min_charge": 2,
        "max_acidic": 2
      },
      "priority_thresholds": [10, 7]
    },
    "passive_permeable": {
      "display_name": "Passive Permeable",
      "ring_sizes": [5, 6, 7, 8, 9],
      "charge_range": [-1, 1],
      "recipes": {
        "5": {"all_hydrophobic": [2, 3], "backbone_modifier": [1, 2]},
        "6": {"all_hydrophobic": [3, 4], "backbone_modifier": [1, 2]},
        "7": {"all_hydrophobic": [3, 4], "backbone_modifier": [2, 3]},
        "8": {"all_hydrophobic": [3, 5], "backbone_modifier": [2, 3]},
        "9": {"all_hydrophobic": [3, 5], "backbone_modifier": [2, 3]}
      },
      "filler_bins": ["turn_inducer", "neutral_hydrophobic"],
      "hard_filters": {
        "max_ring_size": 11,
        "max_abs_charge": 2
      },
      "priority_thresholds": [9, 6]
    }
  }
}
```

## Design Context

### CPP-like mode
Small constrained ring + 3-4 Arg-like residues + 1-2 aromatic hydrophobes + amphipathic 3D display.
Best starting design: 7-8 residues, +3 to +5 charge, mixed D/L stereochemistry.
Cationic face for membrane binding, aromatic residues for membrane insertion and endosomal escape.

### Passive permeable mode
Small hydrophobic/chameleonic ring + low exposed polarity + N-methylation + intramolecular H-bonding.
Best starting design: 6-9 residues, net charge -1 to +1, 2-3 backbone modifiers (NME/NCY).
Reduce exposed backbone HBD, promote intramolecular H-bonds, steric backbone shielding.

### Block bins reference
- **all_cationic** (21): Arg/Lys/amine blocks, charge +1
- **guanidinium_cation** (3): Arg-like with guanidinium group
- **aromatic_hydrophobe** (78): Phe/Trp/Tyr-like, neutral
- **all_hydrophobic** (214): neutral, no HBD/HBA
- **bulky_aliphatic** (106): large neutral aliphatics
- **backbone_modifier** (47): NME (2) + NCY/proline-like (45)
- **turn_inducer** (125): D-residues + proline-like + glycine
- **neutral_polar** (155): neutral, has HBD or HBA
- **acidic_anion** (10): charge -1
