#!/usr/bin/env python3
"""
Render a demo grid figure: each peptide shown as 2D + 3D side by side.
Layout: N rows × 2 columns, each cell = [label, 2D, 3D].

Usage:
    # From JSON candidates file:
    python render_demo_grid.py --input candidates.json --output fig_demo.png

    # From sequence list directly:
    python render_demo_grid.py --sequences "A11-A02-i01-a31-A05-A03,A01-A04-A03-i01-i01-A09" \
        --output fig_demo.png

    # Custom layout:
    python render_demo_grid.py --input candidates.json --cols 3 --output fig_demo.png

JSON format:
    [{"sequence": "A11-A02-i01-a31-A05-A03", "alogp": 2.35, "mw": 649}, ...]
    Fields besides "sequence" are optional — will be computed if missing.

Requires: rdkit, pymol (headless), PIL/Pillow
Environment: /public/home/genesis/miniconda3/envs/bio_env/bin/python
"""

import argparse
import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import Crippen, Descriptors, rdDepictor, AllChem
from rdkit.Chem.Draw import rdMolDraw2D

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from restoken.src.cyclic_assembler import CyclicPeptideAssembler

rdDepictor.SetPreferCoordGen(True)

PYMOL_BIN = "/public/home/genesis/miniconda3/envs/bio_env/bin/pymol"


def render_2d(smiles, output_path, width=800, height=600):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False
    rdDepictor.Compute2DCoords(mol)
    drawer = rdMolDraw2D.MolDraw2DCairo(width, height)
    opts = drawer.drawOptions()
    opts.addStereoAnnotation = True
    opts.bondLineWidth = 2.0
    opts.padding = 0.15
    drawer.DrawMolecule(mol)
    drawer.FinishDrawing()
    with open(output_path, "wb") as f:
        f.write(drawer.GetDrawingText())
    return True


def write_pymol_script(workdir, n_mols, ray_width=800, ray_height=600):
    script = f"""
import pymol
from pymol import cmd
import os
import numpy as np
from PIL import Image

tmpdir = {workdir!r}

for i in range({n_mols}):
    cmd.delete("all")
    cmd.load(os.path.join(tmpdir, f"mol_{{i}}.sdf"), f"mol_{{i}}", format="sdf")
    cmd.set("connect_mode", 1)
    cmd.set("connect_cutoff", 0.1)
    cmd.hide("everything")
    cmd.show("sticks", f"mol_{{i}}")
    cmd.set("stick_radius", 0.15)
    cmd.color("gray50", f"mol_{{i}} and elem C")
    cmd.color("red", f"mol_{{i}} and elem O")
    cmd.color("blue", f"mol_{{i}} and elem N")
    cmd.color("tv_yellow", f"mol_{{i}} and elem S")
    cmd.set("stick_color", "gray20")
    cmd.bg_color("white")
    cmd.set("ray_opaque_background", 1)
    cmd.set("antialias", 2)
    cmd.set("ray_trace_mode", 0)

    best_img = None
    best_score = -1
    for rot_x in [0, 45, 90, 135]:
        for rot_y in [0, 90]:
            cmd.orient(f"mol_{{i}}")
            cmd.zoom(f"mol_{{i}}", buffer=3.0)
            cmd.rotate("x", rot_x)
            cmd.rotate("y", rot_y)
            out = os.path.join(tmpdir, f"mol_{{i}}_3d_test.png")
            cmd.ray({ray_width}, {ray_height})
            cmd.png(out, dpi=150)
            img = Image.open(out)
            arr = np.array(img)
            non_white = np.any(arr < 250, axis=-1)
            if non_white.any():
                rows = np.where(non_white.any(axis=1))[0]
                cols = np.where(non_white.any(axis=0))[0]
                spread = (rows[-1] - rows[0]) * (cols[-1] - cols[0])
                if spread > best_score:
                    best_score = spread
                    best_img = (rot_x, rot_y)

    cmd.orient(f"mol_{{i}}")
    cmd.zoom(f"mol_{{i}}", buffer=3.0)
    if best_img:
        cmd.rotate("x", best_img[0])
        cmd.rotate("y", best_img[1])
    cmd.ray({ray_width}, {ray_height})
    cmd.png(os.path.join(tmpdir, f"mol_{{i}}_3d.png"), dpi=150)
    print(f"  3D rendered: mol_{{i}}")

cmd.quit()
"""
    script_path = os.path.join(workdir, "render_3d.py")
    with open(script_path, "w") as f:
        f.write(script)
    return script_path


def render_3d_all(workdir, n_mols):
    script_path = write_pymol_script(workdir, n_mols)
    env = os.environ.copy()
    env["DISPLAY"] = ""
    cmd = [PYMOL_BIN, "-cq", script_path]
    result = subprocess.run(cmd, capture_output=True, text=True, env=env, timeout=600)
    if result.returncode != 0:
        print(f"PyMOL stderr: {result.stderr[:500]}", file=sys.stderr)
    return result.returncode == 0


def autocrop(img, padding=10):
    """Crop whitespace around the molecule."""
    import numpy as np
    arr = np.array(img.convert("RGB"))
    non_white = np.any(arr < 250, axis=-1)
    if not non_white.any():
        return img
    rows = np.where(non_white.any(axis=1))[0]
    cols = np.where(non_white.any(axis=0))[0]
    top = max(0, rows[0] - padding)
    bottom = min(arr.shape[0], rows[-1] + padding)
    left = max(0, cols[0] - padding)
    right = min(arr.shape[1], cols[-1] + padding)
    return img.crop((left, top, right, bottom))


def composite_grid(workdir, peptides, output_path, cols=2,
                   cell_width=600, cell_height=500, label_height=30):
    from PIL import Image, ImageDraw, ImageFont

    n = len(peptides)
    rows = (n + cols - 1) // cols

    fig_w = cols * cell_width * 2  # 2D + 3D per peptide
    fig_h = rows * (cell_height + label_height)
    canvas = Image.new("RGB", (fig_w, fig_h), "white")
    draw = ImageDraw.Draw(canvas)

    try:
        font = ImageFont.truetype("/usr/share/fonts/dejavu/DejaVuSans.ttf", 16)
    except (OSError, IOError):
        try:
            font = ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf", 16)
        except (OSError, IOError):
            font = ImageFont.load_default()

    for idx, pep in enumerate(peptides):
        row = idx // cols
        col = idx % cols

        x_base = col * cell_width * 2
        y_base = row * (cell_height + label_height)

        seq = pep["sequence"].replace("-", "-")
        label_parts = [seq]
        if "alogp" in pep:
            label_parts.append(f"AlogP={pep['alogp']:.2f}")
        if "mw" in pep:
            label_parts.append(f"MW={pep['mw']:.0f}")
        label = "    ".join(label_parts)

        bbox = draw.textbbox((0, 0), label, font=font)
        tw = bbox[2] - bbox[0]
        tx = x_base + (cell_width * 2 - tw) // 2
        draw.text((tx, y_base + 4), label, fill="black", font=font)

        img_y = y_base + label_height

        # 2D image
        path_2d = os.path.join(workdir, f"mol_{idx}_2d.png")
        if os.path.exists(path_2d):
            img_2d = autocrop(Image.open(path_2d))
            scale = min(cell_width / img_2d.width, cell_height / img_2d.height) * 0.9
            new_w = int(img_2d.width * scale)
            new_h = int(img_2d.height * scale)
            img_2d = img_2d.resize((new_w, new_h), Image.LANCZOS)
            paste_x = x_base + (cell_width - new_w) // 2
            paste_y = img_y + (cell_height - new_h) // 2
            canvas.paste(img_2d, (paste_x, paste_y))

        # 3D image
        path_3d = os.path.join(workdir, f"mol_{idx}_3d.png")
        if os.path.exists(path_3d):
            img_3d = autocrop(Image.open(path_3d))
            scale = min(cell_width / img_3d.width, cell_height / img_3d.height) * 0.9
            new_w = int(img_3d.width * scale)
            new_h = int(img_3d.height * scale)
            img_3d = img_3d.resize((new_w, new_h), Image.LANCZOS)
            paste_x = x_base + cell_width + (cell_width - new_w) // 2
            paste_y = img_y + (cell_height - new_h) // 2
            canvas.paste(img_3d, (paste_x, paste_y))

    canvas.save(output_path, dpi=(300, 300))
    print(f"Saved composite: {output_path} ({canvas.width}x{canvas.height})")


def main():
    parser = argparse.ArgumentParser(description="Render cyclic peptide demo grid (2D+3D)")
    parser.add_argument("--input", help="JSON file with peptide candidates")
    parser.add_argument("--sequences", help="Comma-separated sequences (alternative to --input)")
    parser.add_argument("--output", required=True, help="Output PNG path")
    parser.add_argument("--cols", type=int, default=2, help="Peptides per row (default: 2)")
    parser.add_argument("--cell-width", type=int, default=600, help="Width per 2D or 3D panel")
    parser.add_argument("--cell-height", type=int, default=500, help="Height per panel")
    parser.add_argument("--backend", default=None, help="Path to backend dictionary JSON")
    parser.add_argument("--workdir", default=None, help="Working directory (default: auto temp)")
    parser.add_argument("--keep-workdir", action="store_true", help="Don't clean up workdir")
    parser.add_argument("--skip-3d", action="store_true", help="Skip 3D rendering (2D only)")
    args = parser.parse_args()

    if args.input:
        with open(args.input) as f:
            peptides = json.load(f)
    elif args.sequences:
        peptides = [{"sequence": s.strip()} for s in args.sequences.split(",")]
    else:
        parser.error("Provide --input or --sequences")

    asm = CyclicPeptideAssembler(backend_path=args.backend)

    workdir = args.workdir or tempfile.mkdtemp(prefix="demo_grid_",
                                                dir="/public/home/genesis/.local/tmp")
    os.makedirs(workdir, exist_ok=True)
    print(f"Working directory: {workdir}")

    for idx, pep in enumerate(peptides):
        seq = pep["sequence"]
        print(f"[{idx+1}/{len(peptides)}] {seq}")

        smi = asm.assemble(seq)
        if smi is None:
            print(f"  SKIP: assembly failed")
            continue
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            print(f"  SKIP: invalid SMILES")
            continue

        if "alogp" not in pep:
            pep["alogp"] = Crippen.MolLogP(mol)
        if "mw" not in pep:
            pep["mw"] = Descriptors.MolWt(mol)

        render_2d(smi, os.path.join(workdir, f"mol_{idx}_2d.png"))
        print(f"  2D rendered")

        if not args.skip_3d:
            result = asm.generate_3d(seq, os.path.join(workdir, f"mol_{idx}.sdf"))
            if result:
                print(f"  3D conformer: E={result['energy']:.1f} kcal/mol")
            else:
                print(f"  WARN: 3D conformer generation failed")

    if not args.skip_3d:
        print("Rendering 3D with PyMOL...")
        render_3d_all(workdir, len(peptides))

    print("Compositing grid...")
    composite_grid(workdir, peptides, args.output,
                   cols=args.cols, cell_width=args.cell_width,
                   cell_height=args.cell_height)

    if not args.keep_workdir and not args.workdir:
        import shutil
        shutil.rmtree(workdir, ignore_errors=True)


if __name__ == "__main__":
    main()
