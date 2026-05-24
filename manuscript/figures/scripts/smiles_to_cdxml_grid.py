#!/usr/bin/env python3
"""
Convert SMILES/SDF to CDXML grid figure with ACS Document 1996 style.

Usage:
    # From JSON (sequence + smiles + labels):
    python smiles_to_cdxml_grid.py --input peptides.json --output fig.cdxml

    # From SMILES directly:
    python smiles_to_cdxml_grid.py --smiles "CCO,CCCO" --names "Ethanol,Propanol" --output fig.cdxml

    # From SDF files:
    python smiles_to_cdxml_grid.py --sdf mol_0.sdf,mol_1.sdf --names "Mol A,Mol B" --output fig.cdxml

    # ResToken sequences (auto-assembles SMILES):
    python smiles_to_cdxml_grid.py --restoken "A11-A02-i01-a31-A05-A03,A01-A04-A03-i01-i01-A09" \
        --output fig.cdxml

    # Custom grid layout:
    python smiles_to_cdxml_grid.py --input peptides.json --cols 3 --output fig.cdxml

JSON format:
    [{"smiles": "CCO", "label": "Ethanol  AlogP=..."}, ...]
    or
    [{"sequence": "A11-A02-...", "alogp": 2.35, "mw": 649}, ...]  (ResToken mode)

Requires: rdkit
Environment: /public/home/genesis/miniconda3/envs/bio_env/bin/python
"""

import argparse
import json
import math
import os
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import Crippen, Descriptors, rdDepictor, AllChem

rdDepictor.SetPreferCoordGen(True)

ACS_BOND_LENGTH = 14.40


def make_root():
    root = ET.Element('CDXML')
    attrs = {
        'CreationProgram': 'RDKit+Python',
        'BondLength': '14.40', 'LineWidth': '0.60', 'BoldWidth': '2',
        'BondSpacing': '18', 'HashSpacing': '2.50', 'MarginWidth': '1.60',
        'ChainAngle': '120', 'LabelFont': '3', 'LabelSize': '10',
        'LabelFace': '96', 'CaptionFont': '3', 'CaptionSize': '10',
        'ShowNonTerminalCarbonLabels': 'no', 'ShowTerminalCarbonLabels': 'no',
        'HideImplicitHydrogens': 'no', 'LabelJustification': 'Auto',
        'CaptionJustification': 'Left', 'FractionalWidths': 'yes',
        'InterpretChemically': 'yes', 'ShowAtomQuery': 'yes',
        'ShowAtomStereo': 'no', 'ShowAtomEnhancedStereo': 'yes',
        'ShowBondQuery': 'yes', 'ShowBondRxn': 'yes',
        'ShowBondStereo': 'no', 'color': '0', 'bgcolor': '1',
    }
    for k, v in attrs.items():
        root.set(k, v)

    ct = ET.SubElement(root, 'colortable')
    for r, g, b in [(1,1,1), (0,0,0), (1,0,0), (1,1,0), (0,1,0), (0,1,1), (0,0,1), (1,0,1)]:
        ET.SubElement(ct, 'color', r=str(r), g=str(g), b=str(b))

    ft = ET.SubElement(root, 'fonttable')
    ET.SubElement(ft, 'font', id="3", charset="iso-8859-1", name="Arial")
    return root


def prep_mol(mol):
    mol = Chem.RemoveHs(mol)
    Chem.Kekulize(mol, clearAromaticFlags=True)
    rdDepictor.Compute2DCoords(mol)
    Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
    try:
        Chem.WedgeMolBonds(mol, mol.GetConformer())
    except Exception:
        pass
    return mol


def measure_avg_bond_length(mol):
    conf = mol.GetConformer()
    lengths = []
    for bond in mol.GetBonds():
        p1 = conf.GetAtomPosition(bond.GetBeginAtomIdx())
        p2 = conf.GetAtomPosition(bond.GetEndAtomIdx())
        d = math.sqrt((p1.x - p2.x)**2 + (p1.y - p2.y)**2)
        if d > 0.01:
            lengths.append(d)
    return sum(lengths) / len(lengths) if lengths else 1.5


def mol_bbox_scaled(mol, scale):
    conf = mol.GetConformer()
    xs = [conf.GetAtomPosition(i).x * scale for i in range(mol.GetNumAtoms())]
    ys = [conf.GetAtomPosition(i).y * scale for i in range(mol.GetNumAtoms())]
    return min(xs), max(xs), min(ys), max(ys)


def add_fragment(page, mol, cx, cy, scale, nid):
    conf = mol.GetConformer()
    xs = [conf.GetAtomPosition(i).x for i in range(mol.GetNumAtoms())]
    ys = [conf.GetAtomPosition(i).y for i in range(mol.GetNumAtoms())]
    mcx = (min(xs) + max(xs)) / 2
    mcy = (min(ys) + max(ys)) / 2

    frag = ET.SubElement(page, 'fragment', id=str(nid))
    nid += 1
    atom_map = {}

    for atom in mol.GetAtoms():
        ai = atom.GetIdx()
        pos = conf.GetAtomPosition(ai)
        px = cx + (pos.x - mcx) * scale
        py = cy - (pos.y - mcy) * scale

        aid = nid
        nid += 1
        atom_map[ai] = aid

        n = ET.SubElement(frag, 'n', id=str(aid), p=f"{px:.2f} {py:.2f}")

        elem = atom.GetAtomicNum()
        if elem != 6:
            n.set('Element', str(elem))
            nh = atom.GetTotalNumHs()
            n.set('NumHydrogens', str(nh))

        charge = atom.GetFormalCharge()
        if charge != 0:
            n.set('Charge', str(charge))

    for bond in mol.GetBonds():
        bid = nid
        nid += 1
        b_elem = ET.SubElement(frag, 'b', id=str(bid),
                               B=str(atom_map[bond.GetBeginAtomIdx()]),
                               E=str(atom_map[bond.GetEndAtomIdx()]))

        bt = bond.GetBondType()
        if bt == Chem.BondType.DOUBLE:
            b_elem.set('Order', '2')
        elif bt == Chem.BondType.TRIPLE:
            b_elem.set('Order', '3')
        elif bt == Chem.BondType.AROMATIC:
            b_elem.set('Order', '2')

        bd = bond.GetBondDir()
        if bd == Chem.BondDir.BEGINWEDGE:
            b_elem.set('Display', 'WedgeBegin')
        elif bd == Chem.BondDir.BEGINDASH:
            b_elem.set('Display', 'WedgedHashBegin')

    return nid


def add_label(page, text, cx, y, nid, size="10", bold=False):
    t = ET.SubElement(page, 't', id=str(nid), p=f"{cx:.2f} {y:.2f}")
    t.set('Justification', 'Center')
    t.set('CaptionJustification', 'Center')
    t.set('Warning', 'Chemical Interpretation is not possible for this label')
    nid += 1
    face = "1" if bold else "1"
    s = ET.SubElement(t, 's', font="3", size=size, face=face)
    s.text = text
    return nid


def write_cdxml(root, out_path):
    ET.indent(root, space="  ")
    xml_str = ET.tostring(root, encoding='unicode')
    header = ('<?xml version="1.0" encoding="UTF-8" ?>\n'
              '<!DOCTYPE CDXML SYSTEM "http://www.cambridgesoft.com/xml/cdxml.dtd">\n')
    with open(out_path, 'w', encoding='utf-8') as f:
        f.write(header)
        f.write(xml_str)
        f.write('\n')


def load_mol_from_sdf(sdf_path):
    suppl = Chem.SDMolSupplier(str(sdf_path), removeHs=True)
    for mol in suppl:
        if mol is not None:
            return mol
    return None


def build_grid_cdxml(entries, output_path, cols=2, cell_w=200, cell_h=180,
                     label_offset=12, scale_factor=None):
    """
    Build a CDXML grid from a list of (mol, label) pairs.

    entries: list of dicts with keys 'mol' (RDKit mol) and 'label' (str)
    """
    n = len(entries)
    rows = (n + cols - 1) // cols

    page_w = cols * cell_w + 40
    page_h = rows * cell_h + 40

    root = make_root()
    root.set('BoundingBox', f"0 0 {page_w} {page_h}")
    page = ET.SubElement(root, 'page', id="1",
                         BoundingBox=f"0 0 {page_w} {page_h}",
                         HeaderPosition="36", FooterPosition="36",
                         HeightPages="1", WidthPages="1")

    nid = 100

    for idx, entry in enumerate(entries):
        mol = entry['mol']
        label = entry.get('label', '')
        row = idx // cols
        col = idx % cols

        cx = 20 + col * cell_w + cell_w / 2
        cy = 20 + row * cell_h + cell_h / 2 + label_offset

        avg_bl = measure_avg_bond_length(mol)
        if scale_factor:
            scale = scale_factor
        else:
            scale = ACS_BOND_LENGTH / avg_bl

        # Check if molecule fits in cell, scale down if needed
        xmin, xmax, ymin, ymax = mol_bbox_scaled(mol, scale)
        mol_w = xmax - xmin
        mol_h = ymax - ymin
        max_w = cell_w * 0.85
        max_h = (cell_h - label_offset - 10) * 0.85
        if mol_w > max_w or mol_h > max_h:
            shrink = min(max_w / mol_w, max_h / mol_h)
            scale *= shrink

        nid = add_fragment(page, mol, cx, cy, scale, nid)

        if label:
            label_y = 20 + row * cell_h + 8
            nid = add_label(page, label, cx, label_y, nid, size="8")

    write_cdxml(root, output_path)
    print(f"Saved CDXML: {output_path} ({n} structures, {rows}x{cols} grid)")


def main():
    parser = argparse.ArgumentParser(description="Convert SMILES/SDF to CDXML grid")
    parser.add_argument("--input", help="JSON file with entries")
    parser.add_argument("--smiles", help="Comma-separated SMILES strings")
    parser.add_argument("--sdf", help="Comma-separated SDF file paths")
    parser.add_argument("--restoken", help="Comma-separated ResToken sequences")
    parser.add_argument("--names", help="Comma-separated labels (for --smiles/--sdf mode)")
    parser.add_argument("--output", required=True, help="Output CDXML path")
    parser.add_argument("--cols", type=int, default=2, help="Columns in grid (default: 2)")
    parser.add_argument("--cell-width", type=int, default=200, help="Cell width in page units")
    parser.add_argument("--cell-height", type=int, default=180, help="Cell height in page units")
    parser.add_argument("--backend", default=None, help="ResToken backend dict path")
    args = parser.parse_args()

    entries = []

    if args.input:
        with open(args.input) as f:
            data = json.load(f)
        for item in data:
            if 'smiles' in item:
                mol = Chem.MolFromSmiles(item['smiles'])
            elif 'sequence' in item:
                sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
                from restoken.src.cyclic_assembler import CyclicPeptideAssembler
                asm = CyclicPeptideAssembler(backend_path=args.backend)
                smi = asm.assemble(item['sequence'])
                mol = Chem.MolFromSmiles(smi) if smi else None
                if mol and 'alogp' not in item:
                    item['alogp'] = Crippen.MolLogP(mol)
                if mol and 'mw' not in item:
                    item['mw'] = Descriptors.MolWt(mol)
            elif 'sdf' in item:
                mol = load_mol_from_sdf(item['sdf'])
            else:
                continue

            if mol is None:
                print(f"  SKIP: {item}", file=sys.stderr)
                continue

            mol = prep_mol(mol)

            label = item.get('label', '')
            if not label:
                parts = []
                if 'sequence' in item:
                    parts.append(item['sequence'])
                if 'alogp' in item:
                    parts.append(f"AlogP={item['alogp']:.2f}")
                if 'mw' in item:
                    parts.append(f"MW={item['mw']:.0f}")
                if 'name' in item:
                    parts.append(item['name'])
                label = "  ".join(parts)

            entries.append({'mol': mol, 'label': label})

    elif args.smiles:
        smiles_list = [s.strip() for s in args.smiles.split(",")]
        names_list = [n.strip() for n in args.names.split(",")] if args.names else [""] * len(smiles_list)
        for smi, name in zip(smiles_list, names_list):
            mol = Chem.MolFromSmiles(smi)
            if mol:
                mol = prep_mol(mol)
                entries.append({'mol': mol, 'label': name})

    elif args.sdf:
        sdf_list = [s.strip() for s in args.sdf.split(",")]
        names_list = [n.strip() for n in args.names.split(",")] if args.names else [""] * len(sdf_list)
        for sdf_path, name in zip(sdf_list, names_list):
            mol = load_mol_from_sdf(sdf_path)
            if mol:
                mol = prep_mol(mol)
                entries.append({'mol': mol, 'label': name})

    elif args.restoken:
        sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
        from restoken.src.cyclic_assembler import CyclicPeptideAssembler
        asm = CyclicPeptideAssembler(backend_path=args.backend)
        seq_list = [s.strip() for s in args.restoken.split(",")]
        for seq in seq_list:
            smi = asm.assemble(seq)
            if smi is None:
                continue
            mol = Chem.MolFromSmiles(smi)
            if mol is None:
                continue
            mol = prep_mol(mol)
            alogp = Crippen.MolLogP(mol)
            mw = Descriptors.MolWt(mol)
            label = f"{seq}  AlogP={alogp:.2f}  MW={mw:.0f}"
            entries.append({'mol': mol, 'label': label})

    else:
        parser.error("Provide --input, --smiles, --sdf, or --restoken")

    if not entries:
        print("No valid molecules to render.", file=sys.stderr)
        sys.exit(1)

    build_grid_cdxml(entries, args.output, cols=args.cols,
                     cell_w=args.cell_width, cell_h=args.cell_height)


if __name__ == "__main__":
    main()
