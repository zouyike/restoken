"""
End-to-end demo: LLM-guided cyclic peptide design with AlogP filtering.

Pipeline:
    1. Load ResToken LLM dictionary (semantic properties)
    2. Prompt an LLM (Gemini) to generate 6-residue cyclic peptides targeting AlogP 1.5-3.25
    3. Validate block IDs against the library
    4. Assemble valid sequences into macrocyclic SMILES (RDKit)
    5. Compute AlogP and filter by target range
    6. Render 2D structures of hits

Requirements:
    pip install rdkit requests
    export GEMINI_API_KEY=your_key

Usage:
    python alogp_design_demo.py
    python alogp_design_demo.py --n_generate 20 --alogp_min 2.0 --alogp_max 4.0
    python alogp_design_demo.py --output my_peptides.png --smiles_out results.csv
"""

import argparse
import csv
import json
import os
import re
import sys
from pathlib import Path

import requests
from rdkit import Chem
from rdkit.Chem import Crippen, Descriptors, rdDepictor, Draw

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent))
from restoken.src.cyclic_assembler import CyclicPeptideAssembler

rdDepictor.SetPreferCoordGen(True)

_DATA_DIR = Path(__file__).resolve().parent.parent / "data"


def load_llm_dict(path=None):
    p = Path(path) if path else _DATA_DIR / "bb_dict_llm_v11.json"
    with open(p) as f:
        return json.load(f)


def build_prompt(llm_dict, n_generate, alogp_min, alogp_max):
    good_blocks = [
        b for b in llm_dict["blocks"]
        if b["charge_label"] == "neu" and b["polarity_bin"] in ("low", "med")
    ]

    block_lines = [
        f"{b['id']}: class={b['class']}, chir={b['chirality']}, "
        f"bb={b['mc_type']}, nmod={b['mc_nmod']}, bulk={b['sc_bulk']}"
        for b in good_blocks
    ]
    block_table = "\n".join(block_lines)

    target_center = (alogp_min + alogp_max) / 2
    return f"""Design cyclic peptides with moderate lipophilicity (AlogP ~{alogp_min}-{alogp_max}).

Below are pre-filtered building blocks (neutral charge, low/medium polarity).
Use ONLY these block IDs.

AVAILABLE BLOCKS:
{block_table}

RULES:
1. Each peptide has exactly 6 residues
2. Use 3-4 blocks with pol=low and 2-3 with pol=med
3. Include 1-2 blocks with nmod=NME or NCY (N-methylation increases lipophilicity)
4. Use mostly bulk=medium or bulk=large
5. Mix L and D chirality for diversity

Here are example sequences that achieve AlogP ~{target_center:.1f} (VERIFIED):
- A146-A66-A87-a12-e01-A79 (AlogP=2.36)
- A86-a31-a16-N18-N23-A42 (AlogP=2.34)
- A68-A13-a07-a02-N18-N07 (AlogP=1.75)

Generate {n_generate} NEW sequences (different from examples).
Each line: six block IDs separated by dashes. No other text."""


def call_gemini(prompt, api_key, model="gemini-2.5-flash"):
    url = (
        f"https://generativelanguage.googleapis.com/v1beta/"
        f"models/{model}:generateContent?key={api_key}"
    )
    payload = {
        "contents": [{"parts": [{"text": prompt}]}],
        "generationConfig": {"temperature": 0.9, "maxOutputTokens": 8192},
    }
    proxy_url = os.environ.get("HTTPS_PROXY") or os.environ.get("https_proxy")
    proxies = {"https": proxy_url, "http": proxy_url} if proxy_url else None
    resp = requests.post(url, json=payload, proxies=proxies, timeout=120)
    resp.raise_for_status()
    result = resp.json()
    for part in result["candidates"][0]["content"]["parts"]:
        if "text" in part:
            return part["text"]
    return ""


def parse_sequences(raw_text, valid_ids):
    lines = [l.strip() for l in raw_text.strip().split("\n") if l.strip()]
    lines = [re.sub(r"^[\d\.\)\-\*\s`]+", "", l).strip().strip("`") for l in lines]
    parsed = []
    for l in lines:
        ids = l.split("-")
        if len(ids) == 6 and all(bid in valid_ids for bid in ids):
            parsed.append(l)
    return parsed


def main():
    parser = argparse.ArgumentParser(description="LLM-guided cyclic peptide design demo")
    parser.add_argument("--n_generate", type=int, default=12)
    parser.add_argument("--alogp_min", type=float, default=1.5)
    parser.add_argument("--alogp_max", type=float, default=3.25)
    parser.add_argument("--output", default="alogp_design_results.png")
    parser.add_argument("--smiles_out", default=None, help="CSV output path")
    parser.add_argument("--model", default="gemini-2.5-flash")
    parser.add_argument("--backend", default=None)
    args = parser.parse_args()

    api_key = os.environ.get("GEMINI_API_KEY")
    if not api_key:
        print("ERROR: Set GEMINI_API_KEY environment variable", file=sys.stderr)
        sys.exit(1)

    llm_dict = load_llm_dict()
    valid_ids = {b["id"] for b in llm_dict["blocks"]}
    print(f"Loaded {len(valid_ids)} block IDs")

    prompt = build_prompt(llm_dict, args.n_generate, args.alogp_min, args.alogp_max)
    print(f"Calling {args.model} to generate {args.n_generate} candidates...")
    raw = call_gemini(prompt, api_key, model=args.model)

    sequences = parse_sequences(raw, valid_ids)
    print(f"Parsed {len(sequences)} valid sequences from LLM output")

    asm = CyclicPeptideAssembler(backend_path=args.backend)
    results = []
    for seq in sequences:
        try:
            smi = asm.assemble(seq)
            if smi is None:
                continue
            mol = Chem.MolFromSmiles(smi)
            if mol is None:
                continue
            alogp = Crippen.MolLogP(mol)
            mw = Descriptors.MolWt(mol)
            hbd = Descriptors.NumHDonors(mol)
            hba = Descriptors.NumHAcceptors(mol)
            tpsa = Descriptors.TPSA(mol)
            in_range = args.alogp_min <= alogp <= args.alogp_max
            tag = "HIT" if in_range else ""
            print(f"  {seq}: AlogP={alogp:.2f} {tag} | MW={mw:.0f} HBD={hbd} HBA={hba} TPSA={tpsa:.0f}")
            results.append({
                "seq": seq, "smi": smi, "mol": mol,
                "alogp": alogp, "mw": mw, "hbd": hbd, "hba": hba,
                "tpsa": tpsa, "in_range": in_range,
            })
        except Exception as e:
            print(f"  {seq}: error — {e}")

    if not results:
        print("No valid molecules produced.")
        sys.exit(1)

    target_center = (args.alogp_min + args.alogp_max) / 2
    best = sorted(results, key=lambda v: abs(v["alogp"] - target_center))[:6]

    mols, legends = [], []
    for v in best:
        rdDepictor.Compute2DCoords(v["mol"])
        mols.append(v["mol"])
        mark = "HIT" if v["in_range"] else "near"
        legends.append(f"{v['seq']}  |  AlogP={v['alogp']:.2f} ({mark})  MW={v['mw']:.0f}")

    img = Draw.MolsToGridImage(mols, molsPerRow=3, subImgSize=(1000, 800),
                                legends=legends, useSVG=False)
    img.save(args.output)

    n_hit = sum(1 for v in results if v["in_range"])
    print(f"\nResults: {n_hit}/{len(results)} hit AlogP {args.alogp_min}-{args.alogp_max} "
          f"({100 * n_hit / len(results):.0f}%)")
    print(f"AlogP range: {min(v['alogp'] for v in results):.2f} to "
          f"{max(v['alogp'] for v in results):.2f}")
    print(f"Saved: {args.output}")

    if args.smiles_out:
        with open(args.smiles_out, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["sequence", "smiles", "alogp", "mw", "hbd", "hba", "tpsa", "in_range"])
            for v in results:
                w.writerow([v["seq"], v["smi"], f"{v['alogp']:.2f}", f"{v['mw']:.0f}",
                            v["hbd"], v["hba"], f"{v['tpsa']:.0f}", v["in_range"]])
        print(f"Saved: {args.smiles_out}")


if __name__ == "__main__":
    main()
