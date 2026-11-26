#!/usr/bin/env python3
"""Digital Discovery Major Revision — Fair Representation Comparison Benchmark.

Implements:
1. Monomer-SMILES sequence generation (R3-1): LLM outputs 6 semicolon-separated
   amide-capped monomer SMILES instead of full macrocycle SMILES; assembled
   programmatically with CyclicPeptideAssembler.assemble_advanced().
2. Matched-Metadata HELM (R1-2): HELM prompt given the exact same 8 semantic
   properties + tags as ResToken, isolating notation syntax from metadata volume.
3. Compact 100-block 4-notation benchmark for chemistry-pretrained TxGemma-9B-Chat
   (8,192-token context limit, R2-4b).
"""
import argparse
import json
import os
import re
import sys
import time
from pathlib import Path

from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")

ROOT = Path("/scratch/genesis/NCAA_tokenization")
sys.path.insert(0, str(ROOT))
from restoken.src.cyclic_assembler import CyclicPeptideAssembler
from restoken.src.library import BlockLibrary

OUT_ROOT = ROOT / "experiments" / "dd_revision" / "rep_fairness"
OUT_ROOT.mkdir(parents=True, exist_ok=True)


def get_subset_ids(lib: BlockLibrary, k: int = 100):
    """Stratified 100-block subset preserving alpha/beta/gamma proportions."""
    import random
    rng = random.Random(42)
    by_mc = {"alpha": [], "beta": [], "gamma": []}
    for bid in sorted(lib.all_ids):
        by_mc[lib[bid].mc_type].append(bid)
    # 62 alpha, 25 beta, 13 gamma = 100
    chosen = (
        sorted(rng.sample(by_mc["alpha"], 62))
        + sorted(rng.sample(by_mc["beta"], 25))
        + sorted(rng.sample(by_mc["gamma"], 13))
    )
    return chosen


def build_prompt(lib: BlockLibrary, representation: str, n: int, length: int = 6, subset_ids=None):
    ids = subset_ids if subset_ids is not None else sorted(lib.all_ids)
    n_blocks = len(ids)

    if representation == "monomer_smiles":
        lines = []
        for bid in ids:
            b = lib[bid]
            lines.append(
                f"{bid} | {b.aa_smiles} | class={b.aa_class}, chiral={b.chirality}, "
                f"charge={b.charge}, backbone={b.mc_type}"
            )
        table = "\n".join(lines)
        ex1 = "; ".join(lib[ids[i]].aa_smiles for i in [0, 1, 2, 3, 4, 5])
        return (
            f"You are a cyclic peptide designer specializing in noncanonical amino acid (NCAA) peptides.\n\n"
            f"## Task\n\n"
            f"Design {n} unique cyclic peptides, each consisting of {length} noncanonical amino acid residues. "
            f"Instead of writing the cyclized macrocycle SMILES, output each peptide as an ordered sequence of "
            f"{length} monomer SMILES strings (amide-capped fragment form copied from the table below) separated "
            f"by semicolons (;). An external assembler will close the head-to-tail macrocycle.\n\n"
            f"## Available NCAA Residues\n\n"
            f"Below is a list of {n_blocks} NCAA residues available for design. Each entry shows the residue code "
            f"and its monomer SMILES (amide-capped fragment form). You may ONLY use residues from this list.\n\n"
            f"{table}\n\n"
            f"## Constraints\n\n"
            f"- All {length} monomer SMILES must come from the library above\n"
            f"- Sequence length = {length} monomers per line\n\n"
            f"## Output Format\n\n"
            f"Output each peptide as {length} monomer SMILES separated by '; ' on a single line, one peptide per line. "
            f"Do NOT include block IDs, numbering, explanations, or markdown formatting.\n\n"
            f"Example format:\n{ex1}\n\nBegin:"
        )

    elif representation == "helm_8prop":
        lines = [
            "HELM | Class | Chiral | Charge | Backbone | N-mod | Bulk | Polarity | Flexibility | Tags",
            "-" * 90,
        ]
        for bid in ids:
            b = lib[bid]
            tags = ", ".join(b.tags) if hasattr(b, "tags") and b.tags else ""
            lines.append(
                f"{b.helm} | {b.aa_class} | {b.chirality} | {b.charge_label} | {b.mc_type} | "
                f"{b.mc_nmod} | {b.sc_bulk} | {b.polarity_bin} | {b.flex_bin} | [{tags}]"
            )
        table = "\n".join(lines)
        return (
            f"You are a cyclic peptide designer specializing in noncanonical amino acid (NCAA) peptides, "
            f"using HELM (Hierarchical Editing Language for Macromolecules) notation.\n\n"
            f"## Task\n\n"
            f"Design {n} unique cyclic peptides, each consisting of {length} noncanonical amino acid residues. "
            f"Output each peptide in HELM notation.\n\n"
            f"## Available NCAA Monomers\n\n"
            f"Below is the complete list of {n_blocks} NCAA monomers available for design with eight semantic properties. "
            f"You must ONLY use monomer codes from this list.\n\n"
            f"{table}\n\n"
            f"## Constraints\n\n"
            f"- All monomer codes must exist in the library above\n"
            f"- Sequence length = {length}\n\n"
            f"## Output Format\n\n"
            f"Output each cyclic peptide in HELM notation, one per line. Head-to-tail cyclization is encoded "
            f"by connecting residue {length}:R2 to residue 1:R1. Do NOT include explanations or numbering.\n\n"
            f"Example format:\n"
            f"PEPTIDE1{{[A01].[A02].[A03].[A04].[A05].[A06]}}$PEPTIDE1,PEPTIDE1,{length}:R2-1:R1$$$V2.0\n\nBegin:"
        )

    elif representation == "restoken":
        lines = [
            "ID | Class | Chiral | Charge | Backbone | N-mod | Bulk | Polarity | Flexibility",
            "-" * 80,
        ]
        for bid in ids:
            b = lib[bid]
            lines.append(
                f"{bid} | {b.aa_class} | {b.chirality} | {b.charge_label} | {b.mc_type} | "
                f"{b.mc_nmod} | {b.sc_bulk} | {b.polarity_bin} | {b.flex_bin}"
            )
        table = "\n".join(lines)
        return (
            f"You are a cyclic peptide designer using ResToken semantic tokens.\n\n"
            f"## Building Block Library ({n_blocks} NCAA blocks)\n\n"
            f"{table}\n\n"
            f"## Task\n"
            f"Design {n} unique cyclic peptides of length {length} using ONLY block IDs from the table above.\n"
            f"Output one dash-separated sequence per line (e.g. {ids[0]}-{ids[1]}-{ids[2]}-{ids[3]}-{ids[4]}-{ids[5]}). "
            f"No explanations or numbering.\n\nBegin:"
        )

    elif representation == "smiles":
        lines = []
        for bid in ids:
            b = lib[bid]
            lines.append(f"{bid} | {b.aa_smiles} | class={b.aa_class}, chiral={b.chirality}, charge={b.charge}, backbone={b.mc_type}")
        table = "\n".join(lines)
        return (
            f"You are a cyclic peptide designer specializing in noncanonical amino acid (NCAA) peptides.\n\n"
            f"## Task\n"
            f"Design {n} unique cyclic peptides of {length} residues from the {n_blocks} residues below. "
            f"Output each peptide as a single macrocyclic SMILES string (head-to-tail cyclized), one per line. "
            f"No explanations or numbering.\n\n"
            f"## Available NCAA Residues\n\n{table}\n\nBegin:"
        )

    else:
        raise ValueError(representation)


def parse_monomer_smiles(text: str):
    seqs = []
    for line in text.strip().splitlines():
        line = line.strip()
        if not line or line.startswith("#") or line.startswith("```"):
            continue
        line = re.sub(r"^\d+[\.\)]\s*", "", line)
        line = re.sub(r"[`*]", "", line).strip()
        if ";" in line:
            parts = [p.strip() for p in line.split(";") if p.strip()]
        elif "|" in line:
            parts = [p.strip() for p in line.split("|") if p.strip()]
        else:
            continue
        if len(parts) >= 3:
            # Strip accidental leading block IDs like "A01: CNC..."
            cleaned = []
            for p in parts:
                p = re.sub(r"^[A-Za-z]\d{2,3}\s*[:=]\s*", "", p).strip()
                cleaned.append(p)
            seqs.append(cleaned)
    return seqs


def validate_monomer_smiles(seqs, lib: BlockLibrary, allowed_ids, length: int = 6):
    asm = CyclicPeptideAssembler()
    lib_can_iso = {}
    lib_can_flat = {}
    for bid in allowed_ids:
        m = Chem.MolFromSmiles(lib[bid].aa_smiles)
        if m is not None:
            lib_can_iso[Chem.MolToSmiles(m, isomericSmiles=True)] = bid
            lib_can_flat[Chem.MolToSmiles(m, isomericSmiles=False)] = bid

    results = []
    for parts in seqs:
        r = {
            "n_monomers": len(parts),
            "length_valid": len(parts) == length,
            "monomer_parse_valid": False,
            "n_monomers_parsed": 0,
            "library_match_iso": False,
            "library_match_flat": False,
            "n_monomers_in_lib_iso": 0,
            "assembly_valid_any": False,
            "task_valid_closed_set": False,
            "assembled_smiles": None,
        }
        mols = [Chem.MolFromSmiles(p) for p in parts]
        r["n_monomers_parsed"] = sum(m is not None for m in mols)
        r["monomer_parse_valid"] = r["length_valid"] and all(m is not None for m in mols)

        if r["monomer_parse_valid"]:
            iso_s = [Chem.MolToSmiles(m, isomericSmiles=True) for m in mols]
            flat_s = [Chem.MolToSmiles(m, isomericSmiles=False) for m in mols]
            r["n_monomers_in_lib_iso"] = sum(s in lib_can_iso for s in iso_s)
            r["library_match_iso"] = all(s in lib_can_iso for s in iso_s)
            r["library_match_flat"] = all(s in lib_can_flat for s in flat_s)

            # Try programmatic assembly directly from the 6 monomer SMILES
            try:
                cyc = asm.assemble_advanced(parts, crosslinks=None, validate=True)
            except Exception:
                cyc = None
            if cyc is not None:
                cmol = Chem.MolFromSmiles(cyc)
                if cmol is not None:
                    ri = cmol.GetRingInfo()
                    if any(len(ring) >= 12 for ring in ri.AtomRings()):
                        r["assembly_valid_any"] = True
                        r["assembled_smiles"] = cyc
            r["task_valid_closed_set"] = r["library_match_iso"] and r["assembly_valid_any"]
        results.append(r)
    return results


def parse_helm(text: str):
    seqs = []
    for line in text.strip().splitlines():
        line = line.strip()
        if not line:
            continue
        line = re.sub(r"^\d+[\.\)]\s*", "", line)
        line = re.sub(r"[`*]", "", line).strip()
        if re.search(r"PEPTIDE\d+\{", line):
            seqs.append(line)
    return seqs


def validate_helm(seqs, allowed_ids, length: int = 6):
    allowed = set(allowed_ids)
    results = []
    for h in seqs:
        m = re.search(r"(PEPTIDE\d+)\{(.*?)\}", h)
        r = {
            "raw": h,
            "format_valid": bool(m),
            "strict_peptide1": bool(m and m.group(1) == "PEPTIDE1"),
            "length_valid": False,
            "monomers_valid": False,
            "cyclization_valid": False,
            "overall_valid": False,
        }
        if m:
            mons = re.findall(r"\[([^\]]+)\]", m.group(2))
            r["n_monomers"] = len(mons)
            r["length_valid"] = len(mons) == length
            r["monomers_valid"] = len(mons) > 0 and all(x in allowed for x in mons)
            r["cyclization_valid"] = bool(re.search(r"1:R1.*R2|R2.*1:R1", h))
            r["overall_valid"] = r["format_valid"] and r["length_valid"] and r["monomers_valid"]
        results.append(r)
    return results


def parse_restoken(text: str):
    seqs = []
    for line in text.strip().splitlines():
        line = line.strip()
        if not line or line.startswith("#") or line.startswith("```"):
            continue
        line = re.sub(r"^\d+[\.\)]\s*", "", line)
        line = re.sub(r"[`*]", "", line).strip()
        ids = re.findall(r"[A-Za-z]\d{2,3}", line)
        if len(ids) >= 3:
            seqs.append(ids)
    return seqs


def validate_restoken(seqs, allowed_ids, length: int = 6):
    allowed = set(allowed_ids)
    results = []
    for s in seqs:
        ok = (len(s) == length) and all(x in allowed for x in s)
        results.append({"raw": "-".join(s), "length_valid": len(s) == length, "id_valid": all(x in allowed for x in s), "overall_valid": ok})
    return results


def parse_smiles(text: str):
    seqs = []
    for line in text.strip().splitlines():
        line = line.strip()
        if not line or line.startswith("#") or line.startswith("```"):
            continue
        line = re.sub(r"^\d+[\.\)]\s*", "", line)
        line = re.sub(r"[`*]", "", line).strip()
        if len(line) < 15:
            continue
        smiles_chars = set("CNOSPFIBrcnospfib[]()=#+-./\\@%0123456789H")
        ratio = sum(1 for c in line if c in smiles_chars) / max(len(line), 1)
        if ratio > 0.7:
            seqs.append(line)
    return seqs


def validate_smiles(seqs):
    results = []
    for smi in seqs:
        mol = Chem.MolFromSmiles(smi)
        parse_ok = mol is not None
        macro_ok = False
        if mol is not None:
            macro_ok = any(len(ring) >= 12 for ring in mol.GetRingInfo().AtomRings())
        results.append({"raw": smi, "chemical_valid": parse_ok, "has_macrocycle": macro_ok, "task_valid": parse_ok and macro_ok})
    return results


class OpenAIBackend:
    def __init__(self, model="gpt-4o"):
        self.model = model
        self.api_key = os.environ.get("OPENAI_API_KEY", "")

    def generate(self, prompt, temperature=1.0, max_tokens=16384):
        import urllib.request
        url = "https://api.openai.com/v1/chat/completions"
        body = json.dumps({
            "model": self.model,
            "messages": [{"role": "user", "content": prompt}],
            "temperature": temperature,
            "max_completion_tokens": min(max_tokens, 16384),
        }).encode()
        req = urllib.request.Request(
            url, data=body,
            headers={"Content-Type": "application/json", "Authorization": f"Bearer {self.api_key}"},
            method="POST",
        )
        proxy_url = os.environ.get("HTTPS_PROXY", os.environ.get("https_proxy", "http://127.0.0.1:48890"))
        opener = urllib.request.build_opener(urllib.request.ProxyHandler({"https": proxy_url, "http": proxy_url}))
        with opener.open(req, timeout=240) as resp:
            data = json.loads(resp.read())
        return data["choices"][0]["message"]["content"]


class GeminiBackend:
    def __init__(self, model="gemini-2.5-flash"):
        self.model = model
        self.api_key = os.environ.get("GEMINI_API_KEY", "")

    def generate(self, prompt, temperature=1.0, max_tokens=32768):
        import urllib.request
        url = f"https://generativelanguage.googleapis.com/v1beta/models/{self.model}:generateContent?key={self.api_key}"
        body = json.dumps({
            "contents": [{"parts": [{"text": prompt}]}],
            "generationConfig": {"temperature": temperature, "maxOutputTokens": max_tokens},
        }).encode()
        req = urllib.request.Request(url, data=body, headers={"Content-Type": "application/json"}, method="POST")
        proxy_url = os.environ.get("HTTPS_PROXY", os.environ.get("https_proxy", "http://127.0.0.1:48890"))
        opener = urllib.request.build_opener(urllib.request.ProxyHandler({"https": proxy_url, "http": proxy_url}))
        with opener.open(req, timeout=240) as resp:
            data = json.loads(resp.read())
        return data["candidates"][0]["content"]["parts"][0]["text"]


class LocalHFBackend:
    def __init__(self, model_path, load_mode="prequant", max_new=4096):
        import torch
        from transformers import AutoModelForCausalLM, AutoTokenizer
        self.max_new = max_new
        self.tok = AutoTokenizer.from_pretrained(model_path)
        kwargs = {"device_map": "auto"}
        if load_mode == "bf16":
            kwargs["torch_dtype"] = torch.bfloat16
        elif load_mode == "8bit":
            from transformers import BitsAndBytesConfig
            kwargs["quantization_config"] = BitsAndBytesConfig(load_in_8bit=True)
        # prequant models (unsloth bnb-4bit) load their built-in quant config
        self.model = AutoModelForCausalLM.from_pretrained(model_path, **kwargs)
        self.model.eval()
        print(f"[LocalHF] loaded {model_path} mode={load_mode} ctx={getattr(self.model.config, 'max_position_embeddings', '?')}", flush=True)

    def generate(self, prompt, temperature=1.0, max_tokens=4096):
        import torch
        msgs = [{"role": "user", "content": prompt}]
        try:
            text = self.tok.apply_chat_template(msgs, tokenize=False, add_generation_prompt=True, enable_thinking=False)
        except TypeError:
            text = self.tok.apply_chat_template(msgs, tokenize=False, add_generation_prompt=True)
        inp = self.tok(text, return_tensors="pt").to(self.model.device)
        print(f"  prompt tokens: {inp['input_ids'].shape[1]}", flush=True)
        with torch.no_grad():
            out = self.model.generate(
                **inp,
                max_new_tokens=min(max_tokens, self.max_new),
                temperature=max(temperature, 0.01),
                do_sample=True,
                top_p=0.95,
            )
        return self.tok.decode(out[0][inp["input_ids"].shape[1]:], skip_special_tokens=True)


def summarize(rep, results, n_requested):
    n = len(results)
    if rep == "monomer_smiles":
        return {
            "n_requested": n_requested,
            "n_parsed_lines": n,
            "length_valid": sum(r["length_valid"] for r in results),
            "monomer_parse_valid": sum(r["monomer_parse_valid"] for r in results),
            "library_match_iso": sum(r["library_match_iso"] for r in results),
            "library_match_flat": sum(r["library_match_flat"] for r in results),
            "assembly_valid_any": sum(r["assembly_valid_any"] for r in results),
            "task_valid_closed_set": sum(r["task_valid_closed_set"] for r in results),
            "rate_monomer_parse": sum(r["monomer_parse_valid"] for r in results) / max(1, n),
            "rate_library_match_iso": sum(r["library_match_iso"] for r in results) / max(1, n),
            "rate_assembly_any": sum(r["assembly_valid_any"] for r in results) / max(1, n),
            "rate_task_closed_set": sum(r["task_valid_closed_set"] for r in results) / max(1, n),
            "rate_task_closed_set_over_req": sum(r["task_valid_closed_set"] for r in results) / max(1, n_requested),
        }
    elif rep == "helm_8prop":
        v = sum(r["overall_valid"] for r in results)
        c = sum(r["overall_valid"] and r["cyclization_valid"] for r in results)
        return {
            "n_requested": n_requested,
            "n_parsed_lines": n,
            "overall_valid": v,
            "overall_and_cyc_valid": c,
            "validity_rate": v / max(1, n),
            "validity_and_cyc_rate": c / max(1, n),
        }
    elif rep == "restoken":
        v = sum(r["overall_valid"] for r in results)
        return {"n_requested": n_requested, "n_parsed_lines": n, "overall_valid": v, "validity_rate": v / max(1, n)}
    else:
        p = sum(r["chemical_valid"] for r in results)
        t = sum(r["task_valid"] for r in results)
        return {
            "n_requested": n_requested,
            "n_parsed_lines": n,
            "parse_valid": p,
            "task_valid": t,
            "parse_rate": p / max(1, n),
            "task_rate": t / max(1, n),
        }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--backend", required=True, choices=["openai", "gemini", "hf"])
    ap.add_argument("--model", required=True)
    ap.add_argument("--model_path", default=None)
    ap.add_argument("--load_mode", default="prequant")
    ap.add_argument("--reps", nargs="+", default=["monomer_smiles", "helm_8prop"])
    ap.add_argument("--n", type=int, default=200)
    ap.add_argument("--batch_size", type=int, default=50)
    ap.add_argument("--subset_k", type=int, default=0, help="0 = full 400 blocks; >0 = stratified k blocks (e.g. 100 for 8k-ctx models)")
    ap.add_argument("--max_new", type=int, default=4096)
    ap.add_argument("--tag_prefix", default="")
    args = ap.parse_args()

    lib = BlockLibrary()
    subset_ids = get_subset_ids(lib, args.subset_k) if args.subset_k > 0 else None
    allowed_ids = subset_ids if subset_ids is not None else sorted(lib.all_ids)

    if args.backend == "openai":
        be = OpenAIBackend(args.model)
    elif args.backend == "gemini":
        be = GeminiBackend(args.model)
    else:
        be = LocalHFBackend(args.model_path or args.model, load_mode=args.load_mode, max_new=args.max_new)

    for rep in args.reps:
        tag = f"{args.tag_prefix}{args.model}_{rep}_k{len(allowed_ids)}_n{args.n}"
        sum_path = OUT_ROOT / f"{tag}_summary.json"
        if sum_path.exists():
            print(f"SKIP (exists): {tag}", flush=True)
            continue

        print(f"\n=== Running {tag} ===", flush=True)
        n_batches = (args.n + args.batch_size - 1) // args.batch_size
        raws = []
        parsed_all = []
        for bi in range(n_batches):
            bn = min(args.batch_size, args.n - bi * args.batch_size)
            prompt = build_prompt(lib, rep, bn, length=6, subset_ids=subset_ids)
            t0 = time.time()
            raw = be.generate(prompt, temperature=1.0, max_tokens=args.max_new)
            dt = time.time() - t0
            raws.append(raw)
            if rep == "monomer_smiles":
                p = parse_monomer_smiles(raw)
            elif rep == "helm_8prop":
                p = parse_helm(raw)
            elif rep == "restoken":
                p = parse_restoken(raw)
            else:
                p = parse_smiles(raw)
            parsed_all.extend(p)
            print(f"  batch {bi+1}/{n_batches}: parsed={len(p)} (total={len(parsed_all)}) in {dt:.1f}s", flush=True)

        if rep == "monomer_smiles":
            res = validate_monomer_smiles(parsed_all[: args.n], lib, allowed_ids, length=6)
        elif rep == "helm_8prop":
            res = validate_helm(parsed_all[: args.n], allowed_ids, length=6)
        elif rep == "restoken":
            res = validate_restoken(parsed_all[: args.n], allowed_ids, length=6)
        else:
            res = validate_smiles(parsed_all[: args.n])

        summ = summarize(rep, res, args.n)
        summ["model"] = args.model
        summ["representation"] = rep
        summ["n_blocks"] = len(allowed_ids)
        print("  SUMMARY:", json.dumps(summ, indent=2), flush=True)

        (OUT_ROOT / f"{tag}_raw.json").write_text(json.dumps(raws, indent=2))
        (OUT_ROOT / f"{tag}_results.json").write_text(json.dumps(res, indent=2))
        sum_path.write_text(json.dumps(summ, indent=2))


if __name__ == "__main__":
    main()
