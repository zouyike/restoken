# Cell-Permeable Cyclic Peptide

```json
{
  "modes": {
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

### Passive permeable mode
Small hydrophobic/chameleonic ring + low exposed polarity + N-methylation + intramolecular H-bonding.
Best starting design: 6-9 residues, net charge -1 to +1, 2-3 backbone modifiers (NME/NCY).
Reduce exposed backbone HBD, promote intramolecular H-bonds, steric backbone shielding.

### Block bins reference
- **all_hydrophobic** (214): neutral, no HBD/HBA
- **bulky_aliphatic** (106): large neutral aliphatics
- **backbone_modifier** (47): NME (2) + NCY/proline-like (45)
- **turn_inducer** (125): D-residues + proline-like + glycine
- **neutral_hydrophobic** (136): neutral, no polar groups
- **neutral_polar** (155): neutral, has HBD or HBA
