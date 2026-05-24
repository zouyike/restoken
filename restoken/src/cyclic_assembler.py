"""
Assemble ResToken block IDs into a cyclic peptide SMILES and render 2D structure.

Usage:
    from restoken.src.cyclic_assembler import CyclicPeptideAssembler
    asm = CyclicPeptideAssembler()
    smiles = asm.assemble("A01-K03-N15-S14-E02-B07")
    asm.render("A01-K03-N15-S14-E02-B07", "output.png")
"""

import json, io
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import AllChem, rdDepictor, Draw, Descriptors
from rdkit.Chem.Draw import rdMolDraw2D

rdDepictor.SetPreferCoordGen(True)

_DATA_DIR = Path(__file__).resolve().parent.parent / "data"

_AC_SMARTS = Chem.MolFromSmarts("[CH3]C(=O)[N]")
_NME_SMARTS = Chem.MolFromSmarts("C(=O)[NH,N][CH3]")


class CyclicPeptideAssembler:

    def __init__(self, backend_path=None):
        bp = Path(backend_path) if backend_path else _DATA_DIR / "bb_dict_backend_v11.json"
        with open(bp) as f:
            data = json.load(f)
        self._blocks = {b["id"]: b for b in data["blocks"]}

    def assemble(self, sequence_str, validate=True):
        """
        Assemble block IDs into cyclic peptide SMILES.

        Operates at the RDKit molecular graph level: parses each block's
        aa_smiles, strips Ac/NMe caps, forms inter-residue amide bonds,
        and closes the macrocycle.

        Args:
            sequence_str: dash-separated block IDs, e.g. "A01-K03-N15-S14-E02-B07"
            validate: if True, verify round-trip SMILES with RDKit

        Returns:
            str: cyclic peptide SMILES, or None if assembly fails
        """
        ids = [s.strip() for s in sequence_str.split("-")]

        for bid in ids:
            if bid not in self._blocks:
                raise ValueError(f"Unknown block ID: {bid}")

        mols = []
        cap_infos = []
        for bid in ids:
            block = self._blocks[bid]
            mol = Chem.MolFromSmiles(block["structure"]["aa_smiles"])
            if mol is None:
                return None

            ac = mol.GetSubstructMatch(_AC_SMARTS)
            nme = mol.GetSubstructMatch(_NME_SMARTS)
            if not ac or not nme:
                return None

            mols.append(mol)
            cap_infos.append({
                "ac_ch3": ac[0], "ac_co": ac[1], "ac_o": ac[2], "backbone_n": ac[3],
                "nme_co": nme[0], "nme_o": nme[1], "nme_n": nme[2], "nme_ch3": nme[3],
            })

        combined = mols[0]
        offsets = [0]
        for m in mols[1:]:
            offsets.append(combined.GetNumAtoms())
            combined = Chem.CombineMols(combined, m)

        rw = Chem.RWMol(combined)
        n = len(ids)

        for i in range(n):
            o = offsets[i]
            ci = cap_infos[i]
            rw.RemoveBond(o + ci["nme_co"], o + ci["nme_n"])
            rw.RemoveBond(o + ci["ac_co"], o + ci["backbone_n"])

        for i in range(n):
            next_i = (i + 1) % n
            co_idx = offsets[i] + cap_infos[i]["nme_co"]
            n_idx = offsets[next_i] + cap_infos[next_i]["backbone_n"]
            rw.AddBond(co_idx, n_idx, Chem.BondType.SINGLE)

        atoms_to_remove = set()
        for i in range(n):
            o = offsets[i]
            ci = cap_infos[i]
            atoms_to_remove.update([
                o + ci["ac_ch3"], o + ci["ac_co"], o + ci["ac_o"],
                o + ci["nme_n"], o + ci["nme_ch3"],
            ])

        for idx in sorted(atoms_to_remove, reverse=True):
            rw.RemoveAtom(idx)

        try:
            Chem.SanitizeMol(rw)
            smi = Chem.MolToSmiles(rw)
        except Exception:
            return None

        if validate:
            mol = Chem.MolFromSmiles(smi)
            if mol is None:
                return None
            if len(Chem.GetMolFrags(mol)) != 1:
                return None

        return smi

    def render(self, sequence_str, output_path, width=1600, height=1200,
               label=True, bond_width=2.0):
        """
        Assemble and render cyclic peptide as 2D PNG.

        Args:
            sequence_str: dash-separated block IDs
            output_path: path to save PNG
            width, height: image dimensions
            label: if True, add sequence label
            bond_width: line width for bonds
        """
        smi = self.assemble(sequence_str)
        if smi is None:
            raise ValueError(f"Failed to assemble: {sequence_str}")

        mol = Chem.MolFromSmiles(smi)
        rdDepictor.Compute2DCoords(mol)

        drawer = rdMolDraw2D.MolDraw2DCairo(width, height)
        opts = drawer.drawOptions()
        opts.addStereoAnnotation = True
        opts.bondLineWidth = bond_width
        opts.additionalAtomLabelPadding = 0.15

        drawer.DrawMolecule(mol)
        drawer.FinishDrawing()

        png_data = drawer.GetDrawingText()

        if label:
            try:
                from PIL import Image, ImageDraw, ImageFont
                img = Image.open(io.BytesIO(png_data))
                draw_ctx = ImageDraw.Draw(img)
                try:
                    font = ImageFont.truetype("/usr/share/fonts/dejavu/DejaVuSans.ttf", 28)
                except Exception:
                    font = ImageFont.load_default()
                draw_ctx.text((20, height - 50), sequence_str, fill=(0, 0, 0), font=font)
                img.save(output_path)
            except ImportError:
                with open(output_path, "wb") as f:
                    f.write(png_data)
        else:
            with open(output_path, "wb") as f:
                f.write(png_data)

        return smi

    def render_grid(self, sequences, output_path, cols=3, mol_size=(600, 500)):
        """Render multiple cyclic peptides in a grid."""
        mols = []
        legends = []
        for seq in sequences:
            smi = self.assemble(seq)
            if smi:
                mol = Chem.MolFromSmiles(smi)
                if mol:
                    rdDepictor.Compute2DCoords(mol)
                    mols.append(mol)
                    legends.append(seq)

        if not mols:
            raise ValueError("No valid molecules to render")

        img = Draw.MolsToGridImage(
            mols, molsPerRow=cols, subImgSize=mol_size,
            legends=legends, useSVG=False
        )
        img.save(output_path)
        return [Chem.MolToSmiles(m) for m in mols]

    def get_properties(self, sequence_str):
        """Get aggregate properties for a sequence."""
        ids = [s.strip() for s in sequence_str.split("-")]
        props = {
            "n_residues": len(ids),
            "net_charge": 0,
            "total_hbd": 0,
            "total_hba": 0,
            "total_rot": 0,
            "backbone_types": [],
            "n_nme": 0,
            "n_ncy": 0,
            "n_beta": 0,
        }
        for bid in ids:
            b = self._blocks[bid]
            props["net_charge"] += b["charge"]
            props["total_hbd"] += b["sc"]["hbd"]
            props["total_hba"] += b["sc"]["hba"]
            props["total_rot"] += b["rot_total"]
            props["backbone_types"].append(b["mc"]["type"])
            if b["mc"]["nmod"] == "NME":
                props["n_nme"] += 1
            if b["mc"]["nmod"] == "NCY":
                props["n_ncy"] += 1
            if b["mc"]["type"] == "beta":
                props["n_beta"] += 1
        return props

    def generate_3d(self, sequence_str, output_path, n_confs=50,
                    optimize=True, energy_window=50.0):
        """
        Generate 3D conformer for a cyclic peptide and save to file.

        Uses RDKit's ETKDG with macrocycle sampling to produce a
        low-energy 3D structure. Supports SDF, MOL, PDB output formats.

        Args:
            sequence_str: dash-separated block IDs
            output_path: path to save (extension determines format: .sdf, .mol, .pdb)
            n_confs: number of conformers to generate (best is kept)
            optimize: if True, run MMFF force field minimization
            energy_window: kcal/mol window for conformer pruning

        Returns:
            dict with 'smiles', 'mol', 'energy', 'n_atoms', or None on failure
        """
        smi = self.assemble(sequence_str)
        if smi is None:
            raise ValueError(f"Failed to assemble: {sequence_str}")

        mol = Chem.MolFromSmiles(smi)
        mol = Chem.AddHs(mol)

        params = AllChem.ETKDGv3()
        params.useSmallRingTorsions = True
        params.useMacrocycleTorsions = True
        params.numThreads = 0
        params.randomSeed = 42
        params.pruneRmsThresh = 0.5

        cids = AllChem.EmbedMultipleConfs(mol, numConfs=n_confs, params=params)
        if not cids:
            params.useRandomCoords = True
            cids = AllChem.EmbedMultipleConfs(mol, numConfs=n_confs, params=params)
        if not cids:
            return None

        best_cid = 0
        best_energy = float("inf")

        if optimize:
            results = AllChem.MMFFOptimizeMoleculeConfs(
                mol, numThreads=0, maxIters=2000)
            for cid, (converged, energy) in enumerate(results):
                if energy < best_energy:
                    best_energy = energy
                    best_cid = cid
        else:
            ff_props = AllChem.MMFFGetMoleculeProperties(mol)
            if ff_props:
                for cid in cids:
                    ff = AllChem.MMFFGetMoleculeForceField(mol, ff_props, confId=cid)
                    if ff:
                        energy = ff.CalcEnergy()
                        if energy < best_energy:
                            best_energy = energy
                            best_cid = cid

        output_path = str(output_path)
        ext = output_path.rsplit(".", 1)[-1].lower()

        if ext == "pdb":
            Chem.MolToPDBFile(mol, output_path, confId=best_cid)
        elif ext in ("sdf", "mol"):
            writer = Chem.SDWriter(output_path)
            writer.write(mol, confId=best_cid)
            writer.close()
        else:
            Chem.MolToMolFile(mol, output_path, confId=best_cid)

        mol_noh = Chem.RemoveHs(mol)
        return {
            "smiles": Chem.MolToSmiles(mol_noh),
            "mol": mol,
            "conf_id": best_cid,
            "energy": best_energy,
            "n_atoms": mol.GetNumAtoms(),
            "n_heavy_atoms": mol_noh.GetNumAtoms(),
            "output_path": output_path,
        }


def main():
    import argparse
    parser = argparse.ArgumentParser(
        description="Assemble ResToken block IDs into cyclic peptide structure")
    parser.add_argument("sequence", help="Dash-separated block IDs, e.g. A01-K03-N15-S14-E02-B07")
    parser.add_argument("-o", "--output", default="cyclic_peptide.png",
                        help="Output path (.png for 2D, .sdf/.pdb/.mol for 3D)")
    parser.add_argument("--width", type=int, default=1600)
    parser.add_argument("--height", type=int, default=1200)
    parser.add_argument("--smiles-only", action="store_true", help="Print SMILES only, no image")
    parser.add_argument("--mode-3d", action="store_true",
                        help="Generate 3D conformer (output as .sdf/.pdb/.mol)")
    parser.add_argument("--n-confs", type=int, default=50,
                        help="Number of conformers for 3D generation")
    parser.add_argument("--backend", default=None, help="Path to backend dict JSON")
    args = parser.parse_args()

    asm = CyclicPeptideAssembler(backend_path=args.backend)

    if args.smiles_only:
        smi = asm.assemble(args.sequence)
        if smi:
            print(smi)
        else:
            print("ERROR: Assembly failed", file=__import__("sys").stderr)
            raise SystemExit(1)
    elif args.mode_3d:
        ext = args.output.rsplit(".", 1)[-1].lower()
        if ext == "png":
            args.output = args.output.replace(".png", ".sdf")
        result = asm.generate_3d(args.sequence, args.output, n_confs=args.n_confs)
        if result:
            print(f"SMILES: {result['smiles']}")
            print(f"Energy: {result['energy']:.1f} kcal/mol")
            print(f"Atoms: {result['n_heavy_atoms']} heavy, {result['n_atoms']} total")
            print(f"Saved: {result['output_path']}")
        else:
            print("ERROR: 3D generation failed", file=__import__("sys").stderr)
            raise SystemExit(1)
    else:
        smi = asm.render(args.sequence, args.output, width=args.width, height=args.height)
        props = asm.get_properties(args.sequence)
        print(f"SMILES: {smi}")
        print(f"Properties: charge={props['net_charge']}, HBD={props['total_hbd']}, "
              f"HBA={props['total_hba']}, rot={props['total_rot']}, "
              f"beta={props['n_beta']}, NMe={props['n_nme']}, NCY={props['n_ncy']}")
        print(f"Saved: {args.output}")


if __name__ == "__main__":
    main()
