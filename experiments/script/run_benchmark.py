#!/usr/bin/env python3
"""ResToken W2 benchmark runner — LLM generation + validation + metrics.

Supports 3 experiments × 3 representations × multiple model backends.

Usage:
    # Exp1: Zero-shot validity (API model)
    python run_benchmark.py --experiment exp1 --backend gemini \
        --model gemini-2.5-pro --representation restoken --n 200

    # Exp1: Zero-shot validity (local model on GPU)
    python run_benchmark.py --experiment exp1 --backend hf \
        --model Qwen/Qwen2.5-7B-Instruct --representation restoken --n 200

    # Exp2: EDIT/FROZEN controllability (ResToken only)
    python run_benchmark.py --experiment exp2 --backend gemini \
        --model gemini-2.5-flash --n_variants 20

    # Exp3: Property-constrained (all profiles)
    python run_benchmark.py --experiment exp3 --backend gemini \
        --model gemini-2.5-pro --representation restoken \
        --profile permeable --n 100
"""

import argparse
import csv
import json
import os
import re
import sys
import time
import random
from pathlib import Path
from collections import Counter

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from restoken.src.library import BlockLibrary
from restoken.src.validator import SequenceValidator

PROJECT_ROOT = Path(__file__).resolve().parents[2]
PROMPT_DIR = PROJECT_ROOT / "experiments" / "input" / "prompts"
OUTPUT_ROOT = PROJECT_ROOT / "experiments" / "output"

CONSTRAINT_PROFILES = {
    "unconstrained": "- All block IDs must exist in the library\n- Sequence length = {length}",
    "permeable": (
        "- Net charge = 0 (select only neutral blocks, or balance positive + negative)\n"
        "- Total HBD ≤ 1 (minimize hydrogen bond donors for membrane permeability)\n"
        "- At least 2 blocks with N-mod = NME or NCY (N-methylation improves permeability)\n"
        "- At least 3 blocks with Bulk = large (steric shielding of backbone)\n"
        "- Sequence length = {length}"
    ),
    "charged_binder": (
        "- Net charge = +2 (select blocks so total positive charges minus negatives = +2)\n"
        "- Total HBD ≥ 2\n"
        "- At least 1 block with aromatic character (Class = F, Y, W, or H)\n"
        "- Sequence length = {length}"
    ),
    "rigid_scaffold": (
        "- Total rotatable bonds ≤ 18 (select low-flexibility blocks)\n"
        "- At least 3 blocks with Backbone = beta\n"
        "- Flexibility = low for at least half the residues\n"
        "- Sequence length = {length}"
    ),
}


# ─── Model Backends ───────────────────────────────────────────────────────────

class GeminiBackend:
    """Google Gemini API via REST."""

    def __init__(self, model="gemini-2.5-flash"):
        import urllib.request
        self.model = model
        self.api_key = os.environ.get("GEMINI_API_KEY", "")
        self.name = model

    def generate(self, prompt, temperature=1.0, max_tokens=32768):
        import urllib.request, urllib.error
        url = (f"https://generativelanguage.googleapis.com/v1beta/"
               f"models/{self.model}:generateContent?key={self.api_key}")
        body = json.dumps({
            "contents": [{"parts": [{"text": prompt}]}],
            "generationConfig": {
                "temperature": temperature,
                "maxOutputTokens": max_tokens,
                "thinkingConfig": {"thinkingBudget": 4096},
            },
        }).encode()
        req = urllib.request.Request(url, data=body,
                                     headers={"Content-Type": "application/json"},
                                     method="POST")
        proxy_url = os.environ.get("HTTPS_PROXY",
                                    os.environ.get("https_proxy", "http://127.0.0.1:48890"))
        opener = urllib.request.build_opener(
            urllib.request.ProxyHandler({"https": proxy_url, "http": proxy_url}))
        with opener.open(req, timeout=180) as resp:
            result = json.loads(resp.read())
        try:
            return result["candidates"][0]["content"]["parts"][0]["text"]
        except (KeyError, IndexError):
            return f"[ERROR] {json.dumps(result)[:500]}"


class AnthropicBackend:
    """Anthropic Messages API via Claude.ai OAuth."""

    def __init__(self, model="claude-sonnet-4-6"):
        import urllib.request
        self.model = model
        creds_path = os.path.expanduser("~/.claude/.credentials.json")
        with open(creds_path) as f:
            creds = json.load(f)
        self.api_key = creds["claudeAiOauth"]["accessToken"]
        self.name = model

    def generate(self, prompt, temperature=1.0, max_tokens=16384):
        import urllib.request, urllib.error, time
        url = "https://api.anthropic.com/v1/messages"
        body = json.dumps({
            "model": self.model,
            "messages": [{"role": "user", "content": prompt}],
            "temperature": temperature,
            "max_tokens": max_tokens,
        }).encode()
        headers = {
            "Content-Type": "application/json",
            "Authorization": f"Bearer {self.api_key}",
            "anthropic-version": "2023-06-01",
        }
        proxy_url = os.environ.get("HTTPS_PROXY",
                                    os.environ.get("https_proxy", "http://127.0.0.1:48890"))
        opener = urllib.request.build_opener(
            urllib.request.ProxyHandler({"https": proxy_url, "http": proxy_url}))
        for attempt in range(8):
            req = urllib.request.Request(url, data=body, headers=headers, method="POST")
            try:
                with opener.open(req, timeout=180) as resp:
                    result = json.loads(resp.read())
                return result["content"][0]["text"]
            except urllib.error.HTTPError as e:
                if e.code == 429:
                    wait = min(60 * (attempt + 1), 300)
                    print(f"    Rate limited, waiting {wait}s (attempt {attempt+1}/8)...")
                    time.sleep(wait)
                    continue
                if e.code == 401:
                    self._refresh_token()
                    headers["Authorization"] = f"Bearer {self.api_key}"
                    continue
                body_text = e.read().decode()[:300]
                return f"[ERROR] HTTP {e.code}: {body_text}"
            except (KeyError, IndexError):
                return f"[ERROR] {json.dumps(result)[:500]}"
        return "[ERROR] Rate limited after 8 retries"

    def _refresh_token(self):
        creds_path = os.path.expanduser("~/.claude/.credentials.json")
        with open(creds_path) as f:
            creds = json.load(f)
        self.api_key = creds["claudeAiOauth"]["accessToken"]


class OpenAIBackend:
    """OpenAI Chat Completions API."""

    def __init__(self, model="gpt-4o"):
        import urllib.request
        self.model = model
        self.api_key = os.environ.get("OPENAI_API_KEY", "")
        self.name = model

    def generate(self, prompt, temperature=1.0, max_tokens=16384):
        import urllib.request, urllib.error
        url = "https://api.openai.com/v1/chat/completions"
        body = json.dumps({
            "model": self.model,
            "messages": [{"role": "user", "content": prompt}],
            "temperature": temperature,
            "max_completion_tokens": min(max_tokens, 16384),
        }).encode()
        req = urllib.request.Request(
            url, data=body,
            headers={
                "Content-Type": "application/json",
                "Authorization": f"Bearer {self.api_key}",
            },
            method="POST")
        proxy_url = os.environ.get("HTTPS_PROXY",
                                    os.environ.get("https_proxy", "http://127.0.0.1:48890"))
        opener = urllib.request.build_opener(
            urllib.request.ProxyHandler({"https": proxy_url, "http": proxy_url}))
        with opener.open(req, timeout=180) as resp:
            result = json.loads(resp.read())
        try:
            return result["choices"][0]["message"]["content"]
        except (KeyError, IndexError):
            return f"[ERROR] {json.dumps(result)[:500]}"


class HFLocalBackend:
    """Hugging Face transformers local inference."""

    def __init__(self, model_name, load_in_4bit=False):
        from transformers import AutoModelForCausalLM, AutoTokenizer
        import torch

        self.name = model_name.split("/")[-1]
        kwargs = {"torch_dtype": torch.float16, "device_map": "auto"}
        if load_in_4bit:
            from transformers import BitsAndBytesConfig
            kwargs["quantization_config"] = BitsAndBytesConfig(
                load_in_4bit=True, bnb_4bit_compute_dtype=torch.float16)

        self.tokenizer = AutoTokenizer.from_pretrained(model_name)
        self.model = AutoModelForCausalLM.from_pretrained(model_name, **kwargs)
        self.model.eval()

    def generate(self, prompt, temperature=1.0, max_tokens=8192):
        import torch

        messages = [{"role": "user", "content": prompt}]
        text = self.tokenizer.apply_chat_template(
            messages, tokenize=False, add_generation_prompt=True)
        inputs = self.tokenizer(text, return_tensors="pt").to(self.model.device)

        with torch.no_grad():
            out = self.model.generate(
                **inputs,
                max_new_tokens=max_tokens,
                temperature=max(temperature, 0.01),
                do_sample=True,
                top_p=0.95,
                repetition_penalty=1.05,
            )
        return self.tokenizer.decode(out[0][inputs["input_ids"].shape[1]:],
                                      skip_special_tokens=True)


# ─── Prompt Construction ──────────────────────────────────────────────────────

def build_prompt(lib, representation, n, length, profile="unconstrained"):
    """Build a filled prompt for the given representation + profile."""
    constraints = CONSTRAINT_PROFILES[profile].format(length=length)

    if representation == "restoken":
        tmpl = (PROMPT_DIR / "restoken_zeroshot.txt").read_text()
        return tmpl.format(
            N_BLOCKS=lib.size, DICTIONARY_TABLE=lib.llm_dict_text(),
            N_SEQUENCES=n, SEQ_LENGTH=length, CONSTRAINTS=constraints)
    elif representation == "smiles":
        tmpl = (PROMPT_DIR / "smiles_zeroshot.txt").read_text()
        lines = []
        for bid in sorted(lib.all_ids):
            b = lib[bid]
            lines.append(f"{bid} | {b.aa_smiles} | class={b.aa_class}, "
                         f"chiral={b.chirality}, charge={b.charge}, backbone={b.mc_type}")
        return tmpl.format(
            N_BLOCKS=lib.size, SMILES_TABLE="\n".join(lines),
            N_SEQUENCES=n, SEQ_LENGTH=length, CONSTRAINTS=constraints)
    elif representation == "helm":
        tmpl = (PROMPT_DIR / "helm_zeroshot.txt").read_text()
        lines = []
        for bid in sorted(lib.all_ids):
            b = lib[bid]
            lines.append(f"{b.helm} | class={b.aa_class}, chiral={b.chirality}, "
                         f"charge={b.charge_label}, backbone={b.mc_type}, bulk={b.sc_bulk}")
        return tmpl.format(
            N_BLOCKS=lib.size, HELM_TABLE="\n".join(lines),
            N_SEQUENCES=n, SEQ_LENGTH=length, CONSTRAINTS=constraints)
    else:
        raise ValueError(f"Unknown representation: {representation}")


def build_exp2_prompt(lib, parent_ids, edit_positions, n_variants, constraints=""):
    """Build EDIT/FROZEN prompt for experiment 2."""
    tmpl = (PROMPT_DIR / "edit_frozen_restoken.txt").read_text()
    parent_str = "-".join(parent_ids)
    annotations = []
    for i, bid in enumerate(parent_ids):
        tag = "EDIT" if i in edit_positions else "FROZEN"
        annotations.append(f"Position {i+1}: {bid} → {tag}")
    if not constraints:
        constraints = "- All block IDs must exist in the library"
    return tmpl.format(
        DICTIONARY_TABLE=lib.llm_dict_text(),
        PARENT_SEQUENCE=parent_str,
        POSITION_ANNOTATIONS="\n".join(annotations),
        N_VARIANTS=n_variants,
        CONSTRAINTS=constraints)


def build_exp3_prompt(lib, n, length, profile):
    """Build property-constrained prompt for experiment 3."""
    tmpl = (PROMPT_DIR / "property_constrained.txt").read_text()
    constraints = CONSTRAINT_PROFILES[profile].format(length=length)
    return tmpl.format(
        DICTIONARY_TABLE=lib.llm_dict_text(),
        N_SEQUENCES=n, SEQ_LENGTH=length,
        PROPERTY_CONSTRAINTS=constraints)


# ─── Output Parsers ───────────────────────────────────────────────────────────

def parse_restoken_output(text):
    """Extract ResToken ID sequences from free-form LLM output.

    Returns list of lists of token ID strings.
    """
    sequences = []
    for line in text.strip().split("\n"):
        line = line.strip()
        if not line or line.startswith("#") or line.startswith("```"):
            continue
        line = re.sub(r"^\d+[\.\)]\s*", "", line)
        line = re.sub(r"[`*]", "", line).strip()
        ids = re.findall(r"[A-Za-z]\d{2,3}", line)
        if len(ids) >= 3:
            sequences.append(ids)
    return sequences


def parse_smiles_output(text):
    """Extract SMILES strings from LLM output."""
    sequences = []
    for line in text.strip().split("\n"):
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
            sequences.append(line)
        elif ratio > 0.4:
            match = re.search(r"[CNO\[\(][A-Za-z0-9@#\(\)\[\]=+\-\./\\%\{\}]{15,}", line)
            if match:
                sequences.append(match.group())
    return sequences


def parse_helm_output(text):
    """Extract HELM notations from LLM output."""
    sequences = []
    for line in text.strip().split("\n"):
        line = line.strip()
        if not line:
            continue
        line = re.sub(r"^\d+[\.\)]\s*", "", line)
        line = re.sub(r"[`*]", "", line).strip()
        if "PEPTIDE1{" in line or "PEPTIDE1{" in line.replace(" ", ""):
            sequences.append(line)
    return sequences


# ─── Validation ───────────────────────────────────────────────────────────────

def validate_restoken(sequences, lib, length):
    """Validate ResToken sequences. Returns list of result dicts."""
    results = []
    for seq in sequences:
        r = {"raw": "-".join(seq), "n_tokens": len(seq)}
        r["format_valid"] = True
        r["id_valid"] = all(sid in lib.all_ids for sid in seq)
        r["length_valid"] = len(seq) == length
        if r["id_valid"]:
            charge = sum(lib[s].charge for s in seq)
            hbd = sum(lib[s].sc_hbd for s in seq)
            hba = sum(lib[s].sc_hba for s in seq)
            rot = sum(lib[s].rot_total for s in seq)
            r["properties"] = {"charge": charge, "hbd": hbd, "hba": hba, "rot": rot}
        else:
            r["properties"] = {}
            r["bad_ids"] = [s for s in seq if s not in lib.all_ids]
        r["overall_valid"] = r["format_valid"] and r["id_valid"] and r["length_valid"]
        results.append(r)
    return results


def validate_smiles(sequences):
    """Validate SMILES strings with RDKit. Returns list of result dicts."""
    try:
        from rdkit import Chem
        from rdkit.Chem import Descriptors, rdMolDescriptors
    except ImportError:
        return [{"raw": s, "chemical_valid": False, "error": "no rdkit"} for s in sequences]

    results = []
    for smi in sequences:
        r = {"raw": smi, "format_valid": True}
        mol = Chem.MolFromSmiles(smi)
        r["chemical_valid"] = mol is not None
        if mol:
            ri = mol.GetRingInfo()
            r["has_macrocycle"] = any(len(ring) >= 12 for ring in ri.AtomRings())
            r["n_rings"] = ri.NumRings()
            amide = Chem.MolFromSmarts("[C](=O)[NH]")
            n_amide = Chem.MolFromSmarts("[C](=O)[N]")
            r["n_amide"] = len(mol.GetSubstructMatches(n_amide)) if n_amide else 0
            r["mw"] = Descriptors.MolWt(mol)
            r["n_atoms"] = mol.GetNumHeavyAtoms()
        else:
            r["has_macrocycle"] = False
            r["n_rings"] = 0
            r["n_amide"] = 0
            r["mw"] = 0
            r["n_atoms"] = 0
        r["overall_valid"] = r["chemical_valid"]
        results.append(r)
    return results


def validate_helm(sequences, lib, length):
    """Validate HELM sequences. Returns list of result dicts."""
    results = []
    for helm in sequences:
        r = {"raw": helm, "format_valid": False}
        m = re.search(r"PEPTIDE1\{(.*?)\}", helm)
        if m:
            r["format_valid"] = True
            monomer_str = m.group(1)
            monomers = re.findall(r"\[([^\]]+)\]", monomer_str)
            r["monomers"] = monomers
            r["n_monomers"] = len(monomers)
            r["length_valid"] = len(monomers) == length
            r["monomers_valid"] = all(mid in lib.all_ids for mid in monomers)
            if not r["monomers_valid"]:
                r["bad_monomers"] = [mid for mid in monomers if mid not in lib.all_ids]
            has_cycle = bool(re.search(r"1:R1.*R2|R1.*1:R1", helm))
            r["cyclization_valid"] = has_cycle
        else:
            r["length_valid"] = False
            r["monomers_valid"] = False
            r["cyclization_valid"] = False
            r["monomers"] = []
            r["n_monomers"] = 0
        r["overall_valid"] = (r["format_valid"] and r.get("monomers_valid", False)
                              and r.get("length_valid", False))
        results.append(r)
    return results


# ─── Metrics ──────────────────────────────────────────────────────────────────

def compute_diversity_restoken(valid_results, lib):
    """Compute Tanimoto diversity for valid ResToken sequences."""
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
        from rdkit import DataStructs
    except ImportError:
        return {"mean_tanimoto": -1, "n_fps": 0}

    fps = []
    for r in valid_results:
        if not r["overall_valid"]:
            continue
        ids = r["raw"].split("-")
        combined = ".".join(lib[s].aa_smiles for s in ids)
        mol = Chem.MolFromSmiles(combined)
        if mol:
            fps.append(AllChem.GetMorganFingerprintAsBitVect(mol, 2, 2048))

    if len(fps) < 2:
        return {"mean_tanimoto": 0, "n_fps": len(fps)}

    sims = []
    n = min(len(fps), 200)
    sample = random.sample(range(len(fps)), n) if len(fps) > n else range(len(fps))
    for i in sample:
        for j in sample:
            if i < j:
                sims.append(DataStructs.TanimotoSimilarity(fps[i], fps[j]))

    return {
        "mean_tanimoto": sum(sims) / len(sims) if sims else 0,
        "mean_distance": 1 - sum(sims) / len(sims) if sims else 0,
        "n_fps": len(fps),
    }


def summarize_results(results, representation):
    """Compute aggregate metrics from validation results."""
    n = len(results)
    if n == 0:
        return {"n_parsed": 0}
    s = {"n_parsed": n}
    s["format_valid"] = sum(1 for r in results if r.get("format_valid")) / n
    if representation == "restoken":
        s["id_valid"] = sum(1 for r in results if r.get("id_valid")) / n
        s["length_valid"] = sum(1 for r in results if r.get("length_valid")) / n
    elif representation == "smiles":
        s["chemical_valid"] = sum(1 for r in results if r.get("chemical_valid")) / n
        s["has_macrocycle"] = sum(1 for r in results if r.get("has_macrocycle")) / n
        avg_amide = sum(r.get("n_amide", 0) for r in results) / n
        s["avg_amide_bonds"] = round(avg_amide, 1)
    elif representation == "helm":
        s["monomers_valid"] = sum(1 for r in results if r.get("monomers_valid")) / n
        s["length_valid"] = sum(1 for r in results if r.get("length_valid")) / n
        s["cyclization_valid"] = sum(1 for r in results if r.get("cyclization_valid")) / n
    s["overall_valid"] = sum(1 for r in results if r.get("overall_valid")) / n
    valid_seqs = [r["raw"] for r in results if r.get("overall_valid")]
    s["n_valid"] = len(valid_seqs)
    s["n_unique"] = len(set(valid_seqs))
    s["uniqueness"] = s["n_unique"] / max(s["n_valid"], 1)
    return s


# ─── Experiment Runners ───────────────────────────────────────────────────────

def run_exp1(backend, lib, representation, n, length, batch_size, profile):
    """Experiment 1: Zero-shot generation validity."""
    out_dir = OUTPUT_ROOT / "exp1_validity"
    out_dir.mkdir(parents=True, exist_ok=True)
    model_tag = backend.name.replace("/", "_")
    tag = f"{model_tag}_{representation}_{profile}_len{length}"

    all_raw = []
    all_parsed = []
    n_remaining = n
    batch_num = 0

    while n_remaining > 0:
        batch_n = min(batch_size, n_remaining)
        batch_num += 1
        print(f"  Batch {batch_num}: requesting {batch_n} sequences...")

        prompt = build_prompt(lib, representation, batch_n, length, profile)
        t0 = time.time()
        raw = backend.generate(prompt, temperature=1.0)
        elapsed = time.time() - t0
        print(f"    Response received ({elapsed:.1f}s, {len(raw)} chars)")

        all_raw.append({"batch": batch_num, "raw": raw, "elapsed": elapsed})

        if representation == "restoken":
            parsed = parse_restoken_output(raw)
        elif representation == "smiles":
            parsed = parse_smiles_output(raw)
        elif representation == "helm":
            parsed = parse_helm_output(raw)
        else:
            parsed = []

        all_parsed.extend(parsed)
        n_remaining -= batch_n
        print(f"    Parsed {len(parsed)} sequences (total: {len(all_parsed)})")

        if batch_num < (n // batch_size):
            time.sleep(15)

    if representation == "restoken":
        results = validate_restoken(all_parsed, lib, length)
        diversity = compute_diversity_restoken(results, lib)
    elif representation == "smiles":
        results = validate_smiles(all_parsed)
        diversity = {}
    elif representation == "helm":
        results = validate_helm(all_parsed, lib, length)
        diversity = {}
    else:
        results = []
        diversity = {}

    summary = summarize_results(results, representation)
    summary["diversity"] = diversity
    summary["model"] = backend.name
    summary["representation"] = representation
    summary["profile"] = profile
    summary["n_requested"] = n
    summary["length"] = length

    (out_dir / f"{tag}_raw.json").write_text(json.dumps(all_raw, indent=2))
    (out_dir / f"{tag}_results.json").write_text(
        json.dumps(results, indent=2, default=str))
    (out_dir / f"{tag}_summary.json").write_text(json.dumps(summary, indent=2))

    print(f"\n  === {tag} Summary ===")
    for k, v in summary.items():
        if k not in ("diversity", "model", "representation"):
            print(f"    {k}: {v}")
    if diversity:
        print(f"    diversity: {diversity}")
    print(f"  Saved to: {out_dir / tag}_*.json")
    return summary


def run_exp2(backend, lib, n_parents, n_variants, length):
    """Experiment 2: EDIT/FROZEN controllability."""
    out_dir = OUTPUT_ROOT / "exp2_controllability"
    out_dir.mkdir(parents=True, exist_ok=True)
    model_tag = backend.name.replace("/", "_")

    random.seed(42)
    all_ids = sorted(lib.all_ids)
    parents = []
    for _ in range(n_parents):
        seq = random.sample(all_ids, length)
        edit_pos = sorted(random.sample(range(length), 2))
        parents.append({"ids": seq, "edit_positions": edit_pos})

    all_results = []
    all_raw = []
    for pi, parent in enumerate(parents):
        print(f"  Parent {pi+1}/{n_parents}: {'-'.join(parent['ids'])} "
              f"(EDIT positions: {parent['edit_positions']})")

        prompt = build_exp2_prompt(lib, parent["ids"], parent["edit_positions"], n_variants)
        t0 = time.time()
        raw = backend.generate(prompt, temperature=1.0)
        elapsed = time.time() - t0
        all_raw.append({"parent_idx": pi, "parent": "-".join(parent["ids"]),
                        "edit_positions": parent["edit_positions"],
                        "raw": raw, "elapsed": elapsed})

        parsed = parse_restoken_output(raw)
        print(f"    Got {len(parsed)} variants ({elapsed:.1f}s)")

        for seq in parsed:
            r = {"parent": "-".join(parent["ids"]),
                 "edit_positions": parent["edit_positions"],
                 "variant": "-".join(seq),
                 "n_tokens": len(seq)}

            if len(seq) != length:
                r["length_match"] = False
                r["edit_compliance"] = False
                r["frozen_compliance"] = False
            else:
                r["length_match"] = True
                frozen_ok = all(seq[i] == parent["ids"][i]
                                for i in range(length)
                                if i not in parent["edit_positions"])
                edit_ok = all(seq[i] != parent["ids"][i]
                              for i in parent["edit_positions"]
                              if i < len(seq))
                r["frozen_compliance"] = frozen_ok
                r["edit_compliance"] = edit_ok
                r["id_valid"] = all(s in lib.all_ids for s in seq)
                r["positions_changed"] = [i for i in range(length)
                                          if seq[i] != parent["ids"][i]]
            all_results.append(r)

        if pi < n_parents - 1:
            time.sleep(15)

    summary = {
        "model": backend.name,
        "n_parents": n_parents,
        "n_variants_requested": n_variants,
        "n_total_variants": len(all_results),
        "length_match": sum(1 for r in all_results if r.get("length_match", False)) / max(len(all_results), 1),
        "frozen_compliance": sum(1 for r in all_results if r.get("frozen_compliance", False)) / max(len(all_results), 1),
        "edit_compliance": sum(1 for r in all_results if r.get("edit_compliance", False)) / max(len(all_results), 1),
        "id_valid": sum(1 for r in all_results if r.get("id_valid", False)) / max(len(all_results), 1),
    }

    tag = f"{model_tag}_exp2_controllability"
    (out_dir / f"{tag}_raw.json").write_text(json.dumps(all_raw, indent=2))
    (out_dir / f"{tag}_results.json").write_text(
        json.dumps(all_results, indent=2, default=str))
    (out_dir / f"{tag}_summary.json").write_text(json.dumps(summary, indent=2))

    print(f"\n  === Exp2 {model_tag} Summary ===")
    for k, v in summary.items():
        print(f"    {k}: {v}")
    return summary


def run_exp3(backend, lib, representation, n, length, profile):
    """Experiment 3: Property-constrained generation."""
    out_dir = OUTPUT_ROOT / "exp3_property_constrained"
    out_dir.mkdir(parents=True, exist_ok=True)
    model_tag = backend.name.replace("/", "_")
    tag = f"{model_tag}_{representation}_{profile}_len{length}"

    if representation == "restoken":
        prompt = build_exp3_prompt(lib, n, length, profile)
    else:
        prompt = build_prompt(lib, representation, n, length, profile)

    t0 = time.time()
    raw = backend.generate(prompt, temperature=1.0)
    elapsed = time.time() - t0
    print(f"  Response: {len(raw)} chars ({elapsed:.1f}s)")

    if representation == "restoken":
        parsed = parse_restoken_output(raw)
        results = validate_restoken(parsed, lib, length)
    elif representation == "smiles":
        parsed = parse_smiles_output(raw)
        results = validate_smiles(parsed)
    elif representation == "helm":
        parsed = parse_helm_output(raw)
        results = validate_helm(parsed, lib, length)
    else:
        parsed, results = [], []

    valid_results = [r for r in results if r.get("overall_valid")]
    constraint_pass = 0

    if representation == "restoken":
        for r in valid_results:
            ids = r["raw"].split("-")
            props = r.get("properties", {})
            passes = _check_profile_constraints(ids, props, lib, profile, length)
            r["constraint_satisfied"] = passes
            if passes:
                constraint_pass += 1

    summary = summarize_results(results, representation)
    summary["model"] = backend.name
    summary["representation"] = representation
    summary["profile"] = profile
    summary["n_requested"] = n
    if representation == "restoken":
        summary["constraint_satisfaction"] = constraint_pass / max(len(valid_results), 1)
        summary["n_constraint_pass"] = constraint_pass

    (out_dir / f"{tag}_raw.json").write_text(json.dumps([{"raw": raw}], indent=2))
    (out_dir / f"{tag}_results.json").write_text(
        json.dumps(results, indent=2, default=str))
    (out_dir / f"{tag}_summary.json").write_text(json.dumps(summary, indent=2))

    print(f"\n  === Exp3 {tag} Summary ===")
    for k, v in summary.items():
        print(f"    {k}: {v}")
    return summary


def _check_profile_constraints(ids, props, lib, profile, length):
    """Check if a valid sequence satisfies the given constraint profile."""
    if profile == "unconstrained":
        return True
    elif profile == "permeable":
        if props.get("charge", 999) != 0:
            return False
        if props.get("hbd", 999) > 1:
            return False
        nmod_count = sum(1 for s in ids if lib[s].mc_nmod in ("NCY", "NME"))
        if nmod_count < 2:
            return False
        large_count = sum(1 for s in ids if lib[s].sc_bulk == "large")
        if large_count < 3:
            return False
        return True
    elif profile == "charged_binder":
        if props.get("charge", 0) != 2:
            return False
        if props.get("hbd", 0) < 2:
            return False
        aromatic_classes = {"F", "Y", "W", "H"}
        aro_count = sum(1 for s in ids if lib[s].aa_class in aromatic_classes)
        return aro_count >= 1
    elif profile == "rigid_scaffold":
        if props.get("rot", 999) > 18:
            return False
        beta_count = sum(1 for s in ids if lib[s].mc_type == "beta")
        if beta_count < 3:
            return False
        return True
    return True


# ─── SMILES/HELM EDIT/FROZEN Prompts ─────────────────────────────────────────

def build_exp2_smiles_prompt(lib, parent_ids, edit_positions, n_variants):
    """Build EDIT/FROZEN prompt using SMILES representation."""
    parent_smiles = [lib[bid].aa_smiles for bid in parent_ids]
    annotations = []
    for i, (bid, smi) in enumerate(zip(parent_ids, parent_smiles)):
        tag = "EDIT" if i in edit_positions else "FROZEN"
        annotations.append(f"Position {i+1}: {smi} ({bid}) → {tag}")

    lines = []
    for bid in sorted(lib.all_ids):
        b = lib[bid]
        lines.append(f"{bid} | {b.aa_smiles}")

    return (
        f"You are a cyclic peptide designer using NCAA amino acids.\n\n"
        f"## Available Residues (SMILES)\n\n" + "\n".join(lines) +
        f"\n\n## Task: Controlled Editing\n\n"
        f"Parent peptide (6 residues):\n" +
        "\n".join(f"  Pos {i+1}: {smi}" for i, smi in enumerate(parent_smiles)) +
        f"\n\nPosition annotations:\n" + "\n".join(annotations) +
        f"\n\nGenerate {n_variants} variants. For EDIT positions, substitute a "
        f"DIFFERENT residue SMILES from the list. For FROZEN positions, keep the "
        f"exact same SMILES.\n\n"
        f"Output each variant as 6 SMILES strings separated by ' | ', one variant "
        f"per line. No other text.\n\nBegin:\n"
    )


def build_exp2_helm_prompt(lib, parent_ids, edit_positions, n_variants):
    """Build EDIT/FROZEN prompt using HELM representation."""
    annotations = []
    for i, bid in enumerate(parent_ids):
        tag = "EDIT" if i in edit_positions else "FROZEN"
        annotations.append(f"Position {i+1}: [{bid}] → {tag}")

    lines = []
    for bid in sorted(lib.all_ids):
        b = lib[bid]
        lines.append(f"[{bid}] | class={b.aa_class}, backbone={b.mc_type}")

    n = len(parent_ids)
    parent_helm = ("PEPTIDE1{" + ".".join(f"[{bid}]" for bid in parent_ids) +
                   f"}}$PEPTIDE1,PEPTIDE1,1:R1-{n}:R2$$$V2.0")
    return (
        f"You are a cyclic peptide designer using HELM notation.\n\n"
        f"## Available Monomers\n\n" + "\n".join(lines) +
        f"\n\n## Task: Controlled Editing\n\n"
        f"Parent: {parent_helm}\n\n"
        f"Annotations:\n" + "\n".join(annotations) +
        f"\n\nGenerate {n_variants} variants in HELM notation. For EDIT positions, "
        f"substitute a DIFFERENT monomer. For FROZEN, keep the same.\n\n"
        f"Output one HELM string per line, no other text.\n\nBegin:\n"
    )


def run_exp2_all_reps(backend, lib, n_parents, n_variants, length):
    """Run Exp2 for all 3 representations."""
    out_dir = OUTPUT_ROOT / "exp2_controllability"
    out_dir.mkdir(parents=True, exist_ok=True)
    model_tag = backend.name.replace("/", "_")

    random.seed(42)
    all_ids = sorted(lib.all_ids)
    parents = []
    for _ in range(n_parents):
        seq = random.sample(all_ids, length)
        edit_pos = sorted(random.sample(range(length), 2))
        parents.append({"ids": seq, "edit_positions": edit_pos})

    for rep in ["restoken", "smiles", "helm"]:
        print(f"\n--- Exp2: {model_tag} × {rep} ---")
        all_results = []
        for pi, parent in enumerate(parents):
            print(f"  Parent {pi+1}/{n_parents}: {'-'.join(parent['ids'])}")

            if rep == "restoken":
                prompt = build_exp2_prompt(lib, parent["ids"], parent["edit_positions"],
                                           n_variants)
            elif rep == "smiles":
                prompt = build_exp2_smiles_prompt(lib, parent["ids"],
                                                  parent["edit_positions"], n_variants)
            elif rep == "helm":
                prompt = build_exp2_helm_prompt(lib, parent["ids"],
                                                parent["edit_positions"], n_variants)

            t0 = time.time()
            raw = backend.generate(prompt, temperature=1.0)
            elapsed = time.time() - t0

            if rep == "restoken":
                parsed = parse_restoken_output(raw)
            elif rep == "smiles":
                parsed_lines = raw.strip().split("\n")
                parsed = []
                for line in parsed_lines:
                    line = re.sub(r"^\d+[\.\)]\s*", "", line.strip())
                    parts = [p.strip() for p in line.split("|")]
                    if len(parts) == length:
                        parsed.append(parts)
            elif rep == "helm":
                parsed = parse_helm_output(raw)
                parsed_monomers = []
                for h in parsed:
                    m = re.search(r"PEPTIDE1\{(.*?)\}", h)
                    if m:
                        monomers = re.findall(r"\[([^\]]+)\]", m.group(1))
                        parsed_monomers.append(monomers)
                parsed = parsed_monomers

            for seq_data in parsed:
                if rep == "restoken":
                    seq = seq_data
                elif rep == "smiles":
                    seq = seq_data
                elif rep == "helm":
                    seq = seq_data

                r = {"parent": "-".join(parent["ids"]),
                     "edit_positions": parent["edit_positions"],
                     "representation": rep}

                if rep == "restoken":
                    r["variant"] = "-".join(seq)
                    if len(seq) != length:
                        r["length_match"] = False
                        r["frozen_compliance"] = False
                        r["edit_compliance"] = False
                    else:
                        r["length_match"] = True
                        r["id_valid"] = all(s in lib.all_ids for s in seq)
                        r["frozen_compliance"] = all(
                            seq[i] == parent["ids"][i]
                            for i in range(length)
                            if i not in parent["edit_positions"])
                        r["edit_compliance"] = all(
                            seq[i] != parent["ids"][i]
                            for i in parent["edit_positions"])
                elif rep == "smiles":
                    r["variant"] = " | ".join(seq)
                    if len(seq) != length:
                        r["length_match"] = False
                        r["frozen_compliance"] = False
                        r["edit_compliance"] = False
                    else:
                        r["length_match"] = True
                        parent_smiles = [lib[bid].aa_smiles for bid in parent["ids"]]
                        r["frozen_compliance"] = all(
                            seq[i] == parent_smiles[i]
                            for i in range(length)
                            if i not in parent["edit_positions"])
                        r["edit_compliance"] = all(
                            seq[i] != parent_smiles[i]
                            for i in parent["edit_positions"])
                elif rep == "helm":
                    r["variant"] = "-".join(seq)
                    if len(seq) != length:
                        r["length_match"] = False
                        r["frozen_compliance"] = False
                        r["edit_compliance"] = False
                    else:
                        r["length_match"] = True
                        r["monomers_valid"] = all(s in lib.all_ids for s in seq)
                        r["frozen_compliance"] = all(
                            seq[i] == parent["ids"][i]
                            for i in range(length)
                            if i not in parent["edit_positions"])
                        r["edit_compliance"] = all(
                            seq[i] != parent["ids"][i]
                            for i in parent["edit_positions"])

                all_results.append(r)

            print(f"    {len(parsed)} variants parsed ({elapsed:.1f}s)")
            time.sleep(2)

        n_r = len(all_results)
        summary = {
            "model": backend.name,
            "representation": rep,
            "n_parents": n_parents,
            "n_total": n_r,
            "length_match": sum(1 for r in all_results if r.get("length_match")) / max(n_r, 1),
            "frozen_compliance": sum(1 for r in all_results if r.get("frozen_compliance")) / max(n_r, 1),
            "edit_compliance": sum(1 for r in all_results if r.get("edit_compliance")) / max(n_r, 1),
        }

        tag = f"{model_tag}_{rep}_exp2"
        (out_dir / f"{tag}_results.json").write_text(
            json.dumps(all_results, indent=2, default=str))
        (out_dir / f"{tag}_summary.json").write_text(json.dumps(summary, indent=2))
        print(f"  Summary: frozen={summary['frozen_compliance']:.1%}, "
              f"edit={summary['edit_compliance']:.1%}")


# ─── Main ─────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description="ResToken W2 Benchmark Runner")
    parser.add_argument("--experiment", required=True,
                        choices=["exp1", "exp2", "exp3"])
    parser.add_argument("--backend", required=True, choices=["gemini", "hf", "openai", "anthropic"])
    parser.add_argument("--model", required=True,
                        help="Model name (e.g., gemini-2.5-pro, Qwen/Qwen2.5-7B-Instruct)")
    parser.add_argument("--representation", default="restoken",
                        choices=["restoken", "smiles", "helm"])
    parser.add_argument("--n", type=int, default=200,
                        help="Sequences to generate (exp1/exp3)")
    parser.add_argument("--length", type=int, default=6)
    parser.add_argument("--batch_size", type=int, default=50,
                        help="Sequences per API call (exp1)")
    parser.add_argument("--profile", default="unconstrained",
                        choices=list(CONSTRAINT_PROFILES.keys()))
    parser.add_argument("--n_parents", type=int, default=10,
                        help="Number of parent sequences (exp2)")
    parser.add_argument("--n_variants", type=int, default=20,
                        help="Variants per parent (exp2)")
    parser.add_argument("--load_in_4bit", action="store_true",
                        help="Load HF model in 4-bit quantization")
    args = parser.parse_args()

    lib = BlockLibrary()
    print(f"Library loaded: {lib.size} blocks")

    if args.backend == "gemini":
        backend = GeminiBackend(args.model)
    elif args.backend == "openai":
        backend = OpenAIBackend(args.model)
    elif args.backend == "anthropic":
        backend = AnthropicBackend(args.model)
    elif args.backend == "hf":
        backend = HFLocalBackend(args.model, load_in_4bit=args.load_in_4bit)

    print(f"Backend: {backend.name}")
    print(f"Experiment: {args.experiment}")

    if args.experiment == "exp1":
        print(f"Running Exp1: {args.representation}, n={args.n}, "
              f"length={args.length}, profile={args.profile}")
        run_exp1(backend, lib, args.representation, args.n, args.length,
                 args.batch_size, args.profile)

    elif args.experiment == "exp2":
        print(f"Running Exp2: n_parents={args.n_parents}, "
              f"n_variants={args.n_variants}")
        if args.representation == "restoken":
            run_exp2(backend, lib, args.n_parents, args.n_variants, args.length)
        else:
            run_exp2_all_reps(backend, lib, args.n_parents, args.n_variants,
                              args.length)

    elif args.experiment == "exp3":
        print(f"Running Exp3: {args.representation}, profile={args.profile}, n={args.n}")
        run_exp3(backend, lib, args.representation, args.n, args.length,
                 args.profile)

    print("\nDone.")


if __name__ == "__main__":
    main()
