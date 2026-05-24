"""
Guideline-based cyclic peptide generator.

Reads design guidelines from Markdown files and generates cyclic peptide
sequences using constrained sampling or LLM-assisted generation.

Usage:
    from restoken.src.guideline_generator import GuidelineGenerator
    gen = GuidelineGenerator("restoken/guidelines/cell_penetrating_cyclic_peptide.md")

    # Constrained sampling
    candidates = gen.generate(n=50, mode="cpp_like")
    for c in candidates[:5]:
        print(f"{c['sequence']}  score={c['score']:.1f}  charge={c['net_charge']}")

    # Score an existing sequence
    ann = gen.annotate("K06-A09-K06-A09-K06-A09-A05")

    # Build LLM prompt (for use with any LLM API)
    prompt = gen.build_llm_prompt(mode="cpp_like", n=20)
"""

import copy
import csv
import json
import re
import random
from pathlib import Path
from collections import defaultdict
from typing import Optional

from restoken.src.library import BlockLibrary, Block

_DATA_DIR = Path(__file__).resolve().parent.parent / "data"
_GUIDELINES_DIR = Path(__file__).resolve().parent.parent / "guidelines"


# ── Default mode configurations ──────────────────────────────────────

_DEFAULT_CPP = {
    "display_name": "CPP-like Cell-Penetrating",
    "ring_sizes": [6, 7, 8, 9],
    "charge_range": (2, 5),
    "recipes": {
        6: {"all_cationic": (3, 3), "aromatic_hydrophobe": (1, 2)},
        7: {"all_cationic": (3, 4), "aromatic_hydrophobe": (1, 2)},
        8: {"all_cationic": (3, 4), "aromatic_hydrophobe": (2, 2)},
        9: {"all_cationic": (4, 4), "aromatic_hydrophobe": (2, 2)},
    },
    "filler_bins": ["neutral_polar", "turn_inducer", "neutral_hydrophobic"],
}

_DEFAULT_PASSIVE = {
    "display_name": "Passive Permeable",
    "ring_sizes": [5, 6, 7, 8, 9],
    "charge_range": (-1, 1),
    "recipes": {
        5: {"all_hydrophobic": (2, 3), "backbone_modifier": (1, 2)},
        6: {"all_hydrophobic": (3, 4), "backbone_modifier": (1, 2)},
        7: {"all_hydrophobic": (3, 4), "backbone_modifier": (2, 3)},
        8: {"all_hydrophobic": (3, 5), "backbone_modifier": (2, 3)},
        9: {"all_hydrophobic": (3, 5), "backbone_modifier": (2, 3)},
    },
    "filler_bins": ["turn_inducer", "neutral_hydrophobic"],
}


# ── Block classification ─────────────────────────────────────────────

def classify_blocks(library, raw_blocks):
    """Classify all blocks into functional bins for guideline-based generation.

    Bin categories:
      Charge:        guanidinium_cation, primary_amine_cation, all_cationic,
                     acidic_anion, neutral_polar, neutral_hydrophobic
      Hydrophobicity: aromatic_hydrophobe, bulky_aliphatic, all_hydrophobic
      Backbone:      D_residue, N_methyl, proline_like, backbone_modifier,
                     beta_amino_acid, gamma_amino_acid, turn_inducer
      Other:         halogenated

    Args:
        library: BlockLibrary instance
        raw_blocks: dict of block_id -> raw backend dict entry

    Returns:
        (bins, block_to_bins): bins maps bin_name -> [block_ids],
        block_to_bins maps block_id -> {bin_names}
    """
    bins = defaultdict(list)
    block_to_bins = defaultdict(set)

    for block_id in sorted(library.all_ids):
        b = library[block_id]
        raw = raw_blocks.get(block_id, {})
        aro = raw.get("aromatic_level", 0)

        # Charge
        if b.charge > 0:
            if b.aa_class == "R":
                bins["guanidinium_cation"].append(b.id)
                block_to_bins[b.id].add("guanidinium_cation")
            else:
                bins["primary_amine_cation"].append(b.id)
                block_to_bins[b.id].add("primary_amine_cation")
            bins["all_cationic"].append(b.id)
            block_to_bins[b.id].add("all_cationic")
        elif b.charge < 0:
            bins["acidic_anion"].append(b.id)
            block_to_bins[b.id].add("acidic_anion")
        elif b.sc_hbd > 0 or b.sc_hba > 0:
            bins["neutral_polar"].append(b.id)
            block_to_bins[b.id].add("neutral_polar")
        else:
            bins["neutral_hydrophobic"].append(b.id)
            block_to_bins[b.id].add("neutral_hydrophobic")

        # Hydrophobicity
        if aro > 0 and b.charge == 0:
            bins["aromatic_hydrophobe"].append(b.id)
            block_to_bins[b.id].add("aromatic_hydrophobe")

        if b.charge == 0 and b.sc_hbd == 0 and b.sc_hba == 0:
            if b.sc_bulk == "large":
                bins["bulky_aliphatic"].append(b.id)
                block_to_bins[b.id].add("bulky_aliphatic")
            bins["all_hydrophobic"].append(b.id)
            block_to_bins[b.id].add("all_hydrophobic")

        # Backbone features
        if b.chirality == "D":
            bins["D_residue"].append(b.id)
            block_to_bins[b.id].add("D_residue")

        if b.mc_nmod == "NME":
            bins["N_methyl"].append(b.id)
            block_to_bins[b.id].add("N_methyl")

        if b.mc_nmod == "NCY":
            bins["proline_like"].append(b.id)
            block_to_bins[b.id].add("proline_like")

        if b.mc_nmod in ("NME", "NCY"):
            bins["backbone_modifier"].append(b.id)
            block_to_bins[b.id].add("backbone_modifier")

        if b.mc_type == "beta":
            bins["beta_amino_acid"].append(b.id)
            block_to_bins[b.id].add("beta_amino_acid")

        if b.mc_type == "gamma":
            bins["gamma_amino_acid"].append(b.id)
            block_to_bins[b.id].add("gamma_amino_acid")

        # Turn inducers: D-residues, proline-like, and glycine
        if b.chirality == "D" or b.mc_nmod == "NCY" or b.aa_class == "G":
            bins["turn_inducer"].append(b.id)
            block_to_bins[b.id].add("turn_inducer")

        if raw.get("halogen", 0) > 0:
            bins["halogenated"].append(b.id)
            block_to_bins[b.id].add("halogenated")

    return dict(bins), dict(block_to_bins)


# ── Generator ────────────────────────────────────────────────────────

class GuidelineGenerator:
    """Generate cyclic peptides following design guideline rules.

    Reads a Markdown guideline file, classifies the 400 building blocks
    into functional bins, and generates candidate sequences via constrained
    random sampling. Also provides scoring, annotation, filtering, and
    LLM prompt building.
    """

    def __init__(self, guideline_path, backend_path=None):
        self.guideline_path = Path(guideline_path)
        self.guideline_text = self.guideline_path.read_text()
        self.library = BlockLibrary(backend_path=backend_path)

        bp = Path(backend_path) if backend_path else _DATA_DIR / "bb_dict_backend_v11.json"
        with open(bp) as f:
            raw = json.load(f)
        self._raw_blocks = {b["id"]: b for b in raw["blocks"]}

        self.bins, self.block_to_bins = classify_blocks(self.library, self._raw_blocks)
        self.modes = self._parse_guideline()

    # ── Guideline parsing ────────────────────────────────────────────

    def _parse_guideline(self):
        """Detect design modes from guideline and build configurations.

        Looks for CPP-like and passive-permeable keywords. Extracts ring
        size ranges where possible, falls back to defaults.
        """
        text = self.guideline_text
        modes = {}

        if re.search(r"(?i)CPP|cell[- ]penetrat", text):
            cfg = copy.deepcopy(_DEFAULT_CPP)
            sizes = self._extract_ring_sizes(text, r"(?i)CPP.*ring\s*size|2\.1.*[Rr]ing")
            if sizes:
                cfg["ring_sizes"] = sizes
            modes["cpp_like"] = cfg

        if re.search(r"(?i)passive[- ](?:permea|diffus)", text):
            cfg = copy.deepcopy(_DEFAULT_PASSIVE)
            sizes = self._extract_ring_sizes(text, r"(?i)passive.*ring\s*size|3\.1.*[Rr]ing")
            if sizes:
                cfg["ring_sizes"] = sizes
            modes["passive_permeable"] = cfg

        if not modes:
            modes["default"] = copy.deepcopy(_DEFAULT_CPP)

        return modes

    def _extract_ring_sizes(self, text, section_pattern):
        """Extract ring sizes from the section matching the pattern."""
        lines = text.split("\n")
        section_lines = []
        capturing = False
        for line in lines:
            if re.search(section_pattern, line):
                capturing = True
                continue
            elif capturing:
                if re.match(r"^#{1,3}\s", line):
                    break
                section_lines.append(line)

        if not section_lines:
            return None

        section_text = "\n".join(section_lines)
        sizes = set()
        for m in re.finditer(r"(\d+)[–\-](\d+)\s*residues?", section_text):
            lo, hi = int(m.group(1)), int(m.group(2))
            for s in range(lo, min(hi + 1, 16)):
                sizes.add(s)
        return sorted(sizes) if sizes else None

    # ── Public API ───────────────────────────────────────────────────

    def list_modes(self):
        """Return available design mode names."""
        return list(self.modes.keys())

    def block_summary(self):
        """Return a summary of block classification bins."""
        lines = ["Block Classification Summary:", ""]
        for bin_name in sorted(self.bins.keys()):
            ids = self.bins[bin_name]
            lines.append(f"  {bin_name}: {len(ids)} blocks")
            if len(ids) <= 15:
                lines.append(f"    [{', '.join(ids)}]")
        return "\n".join(lines)

    def generate(self, n=50, mode="cpp_like", ring_size=None, seed=None):
        """Generate candidate sequences via constrained sampling.

        Args:
            n: number of candidates to produce
            mode: design mode ("cpp_like", "passive_permeable")
            ring_size: fixed ring size, or None to sample from mode range
            seed: random seed for reproducibility

        Returns:
            list of dicts sorted by score (descending), each containing:
            sequence, score, and annotation fields
        """
        if seed is not None:
            random.seed(seed)

        if mode not in self.modes:
            raise ValueError(f"Unknown mode '{mode}'. Available: {list(self.modes.keys())}")

        config = self.modes[mode]
        candidates = []
        seen = set()

        for _ in range(n * 100):
            seq = self._generate_one(config, ring_size)
            if seq is None:
                continue

            key = tuple(sorted(seq))
            if key in seen:
                continue
            seen.add(key)

            seq_str = "-".join(seq)
            score = self._score_sequence(seq, mode)
            ann = self._annotate_sequence(seq, mode)
            candidates.append({"sequence": seq_str, "blocks": seq, "score": score, **ann})

            if len(candidates) >= n:
                break

        candidates.sort(key=lambda x: x["score"], reverse=True)
        return candidates

    def score(self, sequence_str, mode="cpp_like"):
        """Score a dash-separated sequence string."""
        ids = [s.strip() for s in sequence_str.split("-")]
        for bid in ids:
            if bid not in self.library:
                raise ValueError(f"Unknown block ID: {bid}")
        return self._score_sequence(ids, mode)

    def annotate(self, sequence_str, mode="cpp_like"):
        """Produce full annotation for a sequence.

        Returns a dict matching the output format in guideline section 9.
        """
        ids = [s.strip() for s in sequence_str.split("-")]
        for bid in ids:
            if bid not in self.library:
                raise ValueError(f"Unknown block ID: {bid}")

        score = self._score_sequence(ids, mode)
        ann = self._annotate_sequence(ids, mode)
        return {"sequence": sequence_str, "score": score, **ann}

    def filter(self, candidates, mode="cpp_like", min_score=None):
        """Apply hard filters from the guideline.

        Args:
            candidates: list of annotation dicts (from generate or annotate)
            mode: design mode
            min_score: score cutoff (default: 7 for CPP, 6 for passive)

        Returns:
            filtered list
        """
        if min_score is None:
            min_score = 7.0 if mode in ("cpp_like", "default") else 6.0

        out = []
        for c in candidates:
            if c["score"] < min_score:
                continue
            rs = c["ring_size"]

            if mode in ("cpp_like", "default"):
                if rs > 11:
                    continue
                if c["cationic_count"] < 3:
                    continue
                if c["aromatic_hydrophobe_count"] < 1:
                    continue
                if c["net_charge"] < 2:
                    continue
                if c["acidic_count"] > 2:
                    continue
                if c["aromatic_hydrophobe_count"] > 4 and rs <= 8:
                    continue
            elif mode == "passive_permeable":
                if rs > 11:
                    continue
                if abs(c["net_charge"]) > 2:
                    continue

            out.append(c)
        return out

    def build_llm_prompt(self, mode="cpp_like", n=20, extra_constraints=None):
        """Build an LLM prompt combining the guideline with classified blocks.

        The returned prompt can be sent to any LLM API (Gemini, GPT, Claude).

        Args:
            mode: design mode
            n: number of sequences to request
            extra_constraints: optional additional instruction string

        Returns:
            str: ready-to-use prompt
        """
        if mode not in self.modes:
            raise ValueError(f"Unknown mode: {mode}")

        config = self.modes[mode]

        relevant_bins = set()
        for recipe in config["recipes"].values():
            relevant_bins.update(recipe.keys())
        relevant_bins.update(config["filler_bins"])

        seen_ids = set()
        block_lines = []
        for bin_name in sorted(relevant_bins):
            for bid in self.bins.get(bin_name, []):
                if bid in seen_ids:
                    continue
                seen_ids.add(bid)
                b = self.library[bid]
                block_lines.append(
                    f"{bid}: role={bin_name} | chiral={b.chirality} | charge={b.charge} | "
                    f"bb={b.mc_type} | nmod={b.mc_nmod} | bulk={b.sc_bulk} | "
                    f"hbd={b.sc_hbd} | hba={b.sc_hba}"
                )

        block_table = "\n".join(block_lines)
        ring_str = ", ".join(str(s) for s in config["ring_sizes"])
        lo, hi = config["charge_range"]
        charge_str = f"{lo:+d} to {hi:+d}" if lo != hi else f"{lo:+d}"
        extra = f"\nAdditional constraints: {extra_constraints}" if extra_constraints else ""

        return (
            f"Design cyclic peptides following these guidelines:\n\n"
            f"{self.guideline_text}\n\n---\n\n"
            f"AVAILABLE BUILDING BLOCKS (classified by functional role):\n\n"
            f"{block_table}\n\n"
            f"GENERATION TASK:\n"
            f"- Mode: {config['display_name']}\n"
            f"- Ring size: {ring_str} residues\n"
            f"- Net charge: {charge_str}\n"
            f"- Use ONLY the block IDs listed above\n"
            f"{extra}\n\n"
            f"Generate {n} unique cyclic peptide sequences.\n"
            f"Format: one sequence per line, block IDs separated by dashes.\n"
            f"No other text or explanations."
        )

    # ── Internal: generation ─────────────────────────────────────────

    def _generate_one(self, config, target_ring_size=None):
        rs = target_ring_size or random.choice(config["ring_sizes"])

        recipe = config["recipes"].get(rs)
        if recipe is None:
            closest = min(config["recipes"].keys(), key=lambda k: abs(k - rs))
            recipe = config["recipes"][closest]

        comp = self._sample_composition(rs, recipe, config["filler_bins"])
        if comp is None:
            return None

        blocks = []
        used = set()
        for bin_name, count in comp.items():
            pool = [bid for bid in self.bins.get(bin_name, []) if bid not in used]
            if len(pool) < count:
                return None
            selected = random.sample(pool, count)
            blocks.extend(selected)
            used.update(selected)

        total_charge = sum(self.library[bid].charge for bid in blocks)
        lo, hi = config["charge_range"]
        if not (lo <= total_charge <= hi):
            return None

        random.shuffle(blocks)
        return blocks

    def _sample_composition(self, ring_size, recipe, filler_bins):
        comp = {}
        total = 0
        for bin_name, (mn, _mx) in recipe.items():
            comp[bin_name] = mn
            total += mn

        if total > ring_size:
            return None

        remaining = ring_size - total

        expandable = [(bn, mx - comp[bn]) for bn, (_mn, mx) in recipe.items() if mx > comp[bn]]
        random.shuffle(expandable)
        for bn, headroom in expandable:
            if remaining <= 0:
                break
            add = random.randint(0, min(headroom, remaining))
            comp[bn] += add
            remaining -= add

        if remaining > 0 and filler_bins:
            filler = random.choice(filler_bins)
            comp[filler] = comp.get(filler, 0) + remaining

        return comp

    # ── Internal: scoring ────────────────────────────────────────────

    def _score_sequence(self, block_ids, mode):
        blocks = [self.library[bid] for bid in block_ids]
        if mode in ("cpp_like", "default"):
            return self._score_cpp(block_ids, blocks)
        elif mode == "passive_permeable":
            return self._score_passive(block_ids, blocks)
        return 0.0

    def _score_cpp(self, block_ids, blocks):
        """CPP-like permeability score (guideline section 6.1)."""
        score = 0.0
        rs = len(blocks)

        n_cationic = sum(1 for b in blocks if b.charge > 0)
        score += 2.0 * min(n_cationic, 4)

        n_aro = sum(1 for bid in block_ids if "aromatic_hydrophobe" in self.block_to_bins.get(bid, set()))
        score += 2.0 * min(n_aro, 2)

        if n_cationic >= 3 and n_aro >= 1:
            score += 1.5

        if 6 <= rs <= 8:
            score += 1.0
        elif rs == 9:
            score += 0.5

        n_d = sum(1 for b in blocks if b.chirality == "D")
        n_ncy = sum(1 for b in blocks if b.mc_nmod == "NCY")
        score += 1.0 * min(n_d + n_ncy, 2)

        n_acidic = sum(1 for b in blocks if b.charge < 0)
        score -= 2.0 * n_acidic

        if n_aro > 4:
            score -= 2.0

        total_rot = sum(b.rot_total for b in blocks)
        if total_rot > rs * 5:
            score -= 2.0

        return score

    def _score_passive(self, block_ids, blocks):
        """Passive permeability score (guideline section 6.2)."""
        score = 0.0
        rs = len(blocks)

        n_mod = sum(1 for bid in block_ids if "backbone_modifier" in self.block_to_bins.get(bid, set()))
        score += 2.0 * min(n_mod, 3)

        n_hydro = sum(1 for bid in block_ids if "all_hydrophobic" in self.block_to_bins.get(bid, set()))
        frac = n_hydro / rs if rs else 0
        if 0.4 <= frac <= 0.7:
            score += 2.0
        elif frac > 0.3:
            score += 1.0

        n_bulky = sum(1 for bid in block_ids if "bulky_aliphatic" in self.block_to_bins.get(bid, set()))
        score += 1.5 * min(n_bulky / max(rs, 1), 1.0)

        n_turn = sum(1 for bid in block_ids if "turn_inducer" in self.block_to_bins.get(bid, set()))
        score += 1.0 * min(n_turn, 2)

        n_hbond = sum(1 for b in blocks if b.sc_hbd > 0 and b.sc_hba > 0)
        score += 1.0 * min(n_hbond, 2)

        n_exposed = sum(1 for b in blocks if b.mc_nmod == "NO")
        score -= 0.5 * max(n_exposed - 4, 0)

        net_charge = sum(b.charge for b in blocks)
        score -= 2.0 * abs(net_charge)

        if rs > 9:
            score -= 2.0

        return score

    # ── Internal: annotation ─────────────────────────────────────────

    def _annotate_sequence(self, block_ids, mode="cpp_like"):
        blocks = [self.library[bid] for bid in block_ids]
        rs = len(blocks)

        n_cationic = sum(1 for b in blocks if b.charge > 0)
        n_arg = sum(1 for bid in block_ids if "guanidinium_cation" in self.block_to_bins.get(bid, set()))
        n_d_arg = sum(
            1 for bid in block_ids
            if "guanidinium_cation" in self.block_to_bins.get(bid, set())
            and self.library[bid].chirality == "D"
        )
        n_aro = sum(1 for bid in block_ids if "aromatic_hydrophobe" in self.block_to_bins.get(bid, set()))
        n_hydro = sum(1 for bid in block_ids if "all_hydrophobic" in self.block_to_bins.get(bid, set()))
        n_nme = sum(1 for b in blocks if b.mc_nmod == "NME")
        n_d = sum(1 for b in blocks if b.chirality == "D")
        n_turn = sum(1 for bid in block_ids if "turn_inducer" in self.block_to_bins.get(bid, set()))
        n_acidic = sum(1 for b in blocks if b.charge < 0)
        net_charge = sum(b.charge for b in blocks)
        total_hbd = sum(b.sc_hbd for b in blocks)
        total_hba = sum(b.sc_hba for b in blocks)
        total_rot = sum(b.rot_total for b in blocks)

        if net_charge >= 2 and n_cationic >= 3:
            mode_label = "CPP_like_endocytic"
        elif abs(net_charge) <= 1:
            mode_label = "passive_diffusion"
        else:
            mode_label = mode

        score = self._score_sequence(block_ids, mode)
        if mode in ("cpp_like", "default"):
            priority = "high" if score >= 10 else ("medium" if score >= 7 else "low")
        else:
            priority = "high" if score >= 9 else ("medium" if score >= 6 else "low")

        return {
            "ring_size": rs,
            "cyclization_type": "head_to_tail",
            "net_charge": net_charge,
            "cationic_count": n_cationic,
            "arg_like_count": n_arg,
            "D_arg_like_count": n_d_arg,
            "aromatic_hydrophobe_count": n_aro,
            "hydrophobic_fraction": round(n_hydro / rs, 2) if rs else 0,
            "N_methyl_count": n_nme,
            "D_residue_count": n_d,
            "turn_residue_count": n_turn,
            "acidic_count": n_acidic,
            "total_hbd": total_hbd,
            "total_hba": total_hba,
            "total_rotatable_bonds": total_rot,
            "backbone_types": [b.mc_type for b in blocks],
            "mode_label": mode_label,
            "priority": priority,
        }


# ── CLI ──────────────────────────────────────────────────────────────

def main():
    import argparse
    import sys

    parser = argparse.ArgumentParser(
        description="Generate cyclic peptides from design guidelines")
    sub = parser.add_subparsers(dest="command")

    p_gen = sub.add_parser("generate", help="Generate candidates via constrained sampling")
    p_gen.add_argument("guideline", help="Path to guideline MD file")
    p_gen.add_argument("-n", type=int, default=50, help="Number of candidates")
    p_gen.add_argument("--mode", default="cpp_like", help="Design mode")
    p_gen.add_argument("--ring-size", type=int, default=None)
    p_gen.add_argument("--seed", type=int, default=None)
    p_gen.add_argument("--min-score", type=float, default=None)
    p_gen.add_argument("-o", "--output", default=None, help="CSV output path")
    p_gen.add_argument("--backend", default=None)

    p_cls = sub.add_parser("classify", help="Show block classification summary")
    p_cls.add_argument("guideline", help="Path to guideline MD file")
    p_cls.add_argument("--backend", default=None)

    p_sc = sub.add_parser("score", help="Score a sequence")
    p_sc.add_argument("guideline", help="Path to guideline MD file")
    p_sc.add_argument("sequence", help="Dash-separated block IDs")
    p_sc.add_argument("--mode", default="cpp_like")
    p_sc.add_argument("--backend", default=None)

    p_pr = sub.add_parser("prompt", help="Build LLM prompt from guideline")
    p_pr.add_argument("guideline", help="Path to guideline MD file")
    p_pr.add_argument("--mode", default="cpp_like")
    p_pr.add_argument("-n", type=int, default=20)
    p_pr.add_argument("--backend", default=None)

    args = parser.parse_args()

    if not args.command:
        parser.print_help()
        sys.exit(1)

    gen = GuidelineGenerator(args.guideline, backend_path=getattr(args, "backend", None))

    if args.command == "classify":
        print(gen.block_summary())
        print(f"\nAvailable modes: {gen.list_modes()}")

    elif args.command == "generate":
        candidates = gen.generate(
            n=args.n, mode=args.mode, ring_size=args.ring_size, seed=args.seed
        )
        filtered = gen.filter(candidates, mode=args.mode, min_score=args.min_score)

        print(f"Generated {len(candidates)} candidates, {len(filtered)} passed filters\n")
        for c in filtered[:20]:
            print(
                f"  {c['sequence']}  score={c['score']:.1f}  "
                f"charge={c['net_charge']:+d}  priority={c['priority']}"
            )

        if args.output and filtered:
            fields = ["sequence", "score", "ring_size", "net_charge", "cationic_count",
                       "arg_like_count", "aromatic_hydrophobe_count", "hydrophobic_fraction",
                       "N_methyl_count", "D_residue_count", "turn_residue_count",
                       "acidic_count", "total_hbd", "total_hba", "mode_label", "priority"]
            with open(args.output, "w", newline="") as f:
                w = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
                w.writeheader()
                for c in filtered:
                    w.writerow(c)
            print(f"\nSaved {len(filtered)} candidates to {args.output}")

    elif args.command == "score":
        ann = gen.annotate(args.sequence, mode=args.mode)
        print(f"Sequence: {ann['sequence']}")
        print(f"Score: {ann['score']:.1f} ({ann['priority']})")
        for k in ["ring_size", "net_charge", "cationic_count", "arg_like_count",
                   "aromatic_hydrophobe_count", "hydrophobic_fraction", "N_methyl_count",
                   "D_residue_count", "turn_residue_count", "acidic_count",
                   "total_hbd", "total_hba", "mode_label"]:
            print(f"  {k}: {ann[k]}")

    elif args.command == "prompt":
        print(gen.build_llm_prompt(mode=args.mode, n=args.n))


if __name__ == "__main__":
    main()
