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
    }
  }
}
```

## Design Context

### CPP-like mode
Small constrained ring + 3-4 Arg-like residues + 1-2 aromatic hydrophobes + amphipathic 3D display.
Best starting design: 7-8 residues, +3 to +5 charge, mixed D/L stereochemistry.
Cationic face for membrane binding, aromatic residues for membrane insertion and endosomal escape.

### Block bins reference
- **all_cationic** (21): Arg/Lys/amine blocks, charge +1
- **guanidinium_cation** (3): Arg-like with guanidinium group
- **aromatic_hydrophobe** (78): Phe/Trp/Tyr-like, neutral
- **bulky_aliphatic** (106): large neutral aliphatics
- **turn_inducer** (125): D-residues + proline-like + glycine
- **neutral_polar** (155): neutral, has HBD or HBA
- **neutral_hydrophobic** (136): neutral, no polar groups
