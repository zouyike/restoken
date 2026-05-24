"""
Guideline-based cyclic peptide generator.

Reads design guidelines from Markdown files (with JSON config blocks)
or accepts parameters directly. Generates cyclic peptide sequences via
constrained sampling with scoring, filtering, and annotation.

Built-in guidelines:
    cell_penetrating_cyclic_peptide.md  — CPP-like (Arg-rich, cationic)
    cell_permeable_cyclic_peptide.md    — Passive permeable (hydrophobic, low polarity)
    antibiotic_cyclic_peptide.md        — Antimicrobial (cationic amphipathic)

Usage:
    from restoken.src.guideline_generator import GuidelineGenerator

    gen = GuidelineGenerator("restoken/guidelines/cell_penetrating_cyclic_peptide.md")
    candidates = gen.generate(n=50, mode="cpp_like")
    ann = gen.annotate("K06-A09-N20-A10-K08-A05-A97")
    prompt = gen.build_llm_prompt(mode="cpp_like", n=20)

    # Or from explicit parameters (no file needed)
    gen = GuidelineGenerator.from_params({
        "my_mode": {
            "ring_sizes": [6, 7, 8],
            "charge_range": [0, 2],
            "recipes": {7: {"all_cationic": [2, 3], "aromatic_hydrophobe": [2, 3]}},
            "filler_bins": ["neutral_polar"],
        }
    })
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
# These serve as fallbacks when the guideline file omits fields.

_MODE_DEFAULTS = {
    "cpp_like": {
        "display_name": "CPP-like Cell-Penetrating",
        "ring_sizes": [6, 7, 8, 9],
        "charge_range": [2, 5],
        "recipes": {
            6: {"all_cationic": [3, 3], "aromatic_hydrophobe": [1, 2]},
            7: {"all_cationic": [3, 4], "aromatic_hydrophobe": [1, 2]},
            8: {"all_cationic": [3, 4], "aromatic_hydrophobe": [2, 2]},
            9: {"all_cationic": [4, 4], "aromatic_hydrophobe": [2, 2]},
        },
        "filler_bins": ["neutral_polar", "turn_inducer", "neutral_hydrophobic"],
        "hard_filters": {
            "max_ring_size": 11,
            "min_cationic": 3,
            "min_aromatic": 1,
            "min_charge": 2,
            "max_acidic": 2,
        },
        "priority_thresholds": [10, 7],
    },
    "passive_permeable": {
        "display_name": "Passive Permeable",
        "ring_sizes": [5, 6, 7, 8, 9],
        "charge_range": [-1, 1],
        "recipes": {
            5: {"all_hydrophobic": [2, 3], "backbone_modifier": [1, 2]},
            6: {"all_hydrophobic": [3, 4], "backbone_modifier": [1, 2]},
            7: {"all_hydrophobic": [3, 4], "backbone_modifier": [2, 3]},
            8: {"all_hydrophobic": [3, 5], "backbone_modifier": [2, 3]},
            9: {"all_hydrophobic": [3, 5], "backbone_modifier": [2, 3]},
        },
        "filler_bins": ["turn_inducer", "neutral_hydrophobic"],
        "hard_filters": {
            "max_ring_size": 11,
            "max_abs_charge": 2,
        },
        "priority_thresholds": [9, 6],
    },
}


def _deep_merge(base, override):
    """Deep merge override dict into base, returning a new dict."""
    result = copy.deepcopy(base)
    for k, v in override.items():
        if k in result and isinstance(result[k], dict) and isinstance(v, dict):
            result[k] = _deep_merge(result[k], v)
        else:
            result[k] = copy.deepcopy(v)
    return result


def _normalize_config(cfg):
    """Convert JSON-parsed config to internal format (int keys, tuples)."""
    if "recipes" in cfg:
        cfg["recipes"] = {
            int(k): {bn: tuple(v) for bn, v in recipe.items()}
            for k, recipe in cfg["recipes"].items()
        }
    if "charge_range" in cfg and isinstance(cfg["charge_range"], list):
        cfg["charge_range"] = tuple(cfg["charge_range"])
    if "priority_thresholds" in cfg and isinstance(cfg["priority_thresholds"], list):
        cfg["priority_thresholds"] = tuple(cfg["priority_thresholds"])
    return cfg


# ── Block classification ─────────────────────────────────────────────

def classify_blocks(library, raw_blocks):
    """Classify all blocks into functional bins.

    Bins: guanidinium_cation, primary_amine_cation, all_cationic,
    acidic_anion, neutral_polar, neutral_hydrophobic, aromatic_hydrophobe,
    bulky_aliphatic, all_hydrophobic, D_residue, N_methyl, proline_like,
    backbone_modifier, beta_amino_acid, gamma_amino_acid, turn_inducer,
    halogenated.
    """
    bins = defaultdict(list)
    block_to_bins = defaultdict(set)

    for block_id in sorted(library.all_ids):
        b = library[block_id]
        raw = raw_blocks.get(block_id, {})
        aro = raw.get("aromatic_level", 0)

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

        if aro > 0 and b.charge == 0:
            bins["aromatic_hydrophobe"].append(b.id)
            block_to_bins[b.id].add("aromatic_hydrophobe")

        if b.charge == 0 and b.sc_hbd == 0 and b.sc_hba == 0:
            if b.sc_bulk == "large":
                bins["bulky_aliphatic"].append(b.id)
                block_to_bins[b.id].add("bulky_aliphatic")
            bins["all_hydrophobic"].append(b.id)
            block_to_bins[b.id].add("all_hydrophobic")

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
        if b.chirality == "D" or b.mc_nmod == "NCY" or b.aa_class == "G":
            bins["turn_inducer"].append(b.id)
            block_to_bins[b.id].add("turn_inducer")
        if raw.get("halogen", 0) > 0:
            bins["halogenated"].append(b.id)
            block_to_bins[b.id].add("halogenated")

    return dict(bins), dict(block_to_bins)


# ── Generator ────────────────────────────────────────────────────────

class GuidelineGenerator:
    """Generate cyclic peptides from guideline configs or direct parameters."""

    def __init__(self, guideline_path, backend_path=None):
        """Load from a guideline markdown file with a JSON config block.

        The file should contain a fenced ```json ... ``` block defining
        modes with ring_sizes, charge_range, recipes, filler_bins,
        hard_filters, and priority_thresholds. Any prose after the JSON
        block is stored as context for LLM prompt building.
        """
        self.guideline_path = Path(guideline_path)
        self.guideline_text = self.guideline_path.read_text()
        self._init_library(backend_path)
        self.modes = self._parse_guideline()

    @classmethod
    def from_params(cls, modes=None, backend_path=None):
        """Create generator from explicit parameters, no guideline file.

        Args:
            modes: dict of mode_name -> config dict. Each config can have:
                display_name, ring_sizes, charge_range, recipes, filler_bins,
                hard_filters, priority_thresholds. Missing fields use defaults.
                If None, uses cpp_like + passive_permeable defaults.
            backend_path: path to backend dict JSON

        Example:
            gen = GuidelineGenerator.from_params({
                "custom": {
                    "ring_sizes": [6, 7, 8],
                    "charge_range": [0, 3],
                    "recipes": {7: {"all_cationic": [2,3], "aromatic_hydrophobe": [1,2]}},
                    "filler_bins": ["neutral_polar"],
                }
            })
        """
        gen = cls.__new__(cls)
        gen.guideline_path = None
        gen.guideline_text = ""
        gen._init_library(backend_path)

        if modes is None:
            gen.modes = {
                k: _normalize_config(copy.deepcopy(v))
                for k, v in _MODE_DEFAULTS.items()
            }
        else:
            gen.modes = {}
            for name, user_cfg in modes.items():
                defaults = _MODE_DEFAULTS.get(name, _MODE_DEFAULTS["cpp_like"])
                merged = _deep_merge(defaults, user_cfg)
                gen.modes[name] = _normalize_config(merged)

        return gen

    def _init_library(self, backend_path):
        self.library = BlockLibrary(backend_path=backend_path)
        bp = Path(backend_path) if backend_path else _DATA_DIR / "bb_dict_backend_v11.json"
        with open(bp) as f:
            raw = json.load(f)
        self._raw_blocks = {b["id"]: b for b in raw["blocks"]}
        self.bins, self.block_to_bins = classify_blocks(self.library, self._raw_blocks)

    def _parse_guideline(self):
        """Parse JSON config block from guideline markdown.

        Looks for a fenced ```json ... ``` block containing a "modes" dict.
        Merges with defaults so users only need to specify overrides.
        Falls back to default modes if no JSON block found.
        """
        m = re.search(r"```json\s*\n(.*?)\n```", self.guideline_text, re.DOTALL)
        if not m:
            return {k: _normalize_config(copy.deepcopy(v))
                    for k, v in _MODE_DEFAULTS.items()}

        try:
            config = json.loads(m.group(1))
        except json.JSONDecodeError as e:
            raise ValueError(f"Invalid JSON in guideline: {e}")

        modes = {}
        for name, user_cfg in config.get("modes", {}).items():
            defaults = _MODE_DEFAULTS.get(name, _MODE_DEFAULTS["cpp_like"])
            merged = _deep_merge(defaults, user_cfg)
            modes[name] = _normalize_config(merged)

        return modes if modes else {
            k: _normalize_config(copy.deepcopy(v))
            for k, v in _MODE_DEFAULTS.items()
        }

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
            mode: design mode name
            ring_size: fixed ring size, or None to sample from mode range
            seed: random seed for reproducibility

        Returns:
            list of dicts sorted by score (descending)
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
        """Produce full annotation for a sequence."""
        ids = [s.strip() for s in sequence_str.split("-")]
        for bid in ids:
            if bid not in self.library:
                raise ValueError(f"Unknown block ID: {bid}")
        score = self._score_sequence(ids, mode)
        ann = self._annotate_sequence(ids, mode)
        return {"sequence": sequence_str, "score": score, **ann}

    def filter(self, candidates, mode="cpp_like", min_score=None):
        """Apply hard filters from the mode config.

        Filter thresholds are read from the config's hard_filters dict.
        Score cutoff defaults to the mode's medium priority threshold.
        """
        config = self.modes.get(mode, {})
        hf = config.get("hard_filters", {})
        thresholds = config.get("priority_thresholds", (10, 7))

        if min_score is None:
            min_score = thresholds[1]

        out = []
        for c in candidates:
            if c["score"] < min_score:
                continue
            rs = c["ring_size"]
            if hf.get("max_ring_size") and rs > hf["max_ring_size"]:
                continue
            if hf.get("min_cationic") and c["cationic_count"] < hf["min_cationic"]:
                continue
            if hf.get("min_aromatic") and c["aromatic_hydrophobe_count"] < hf["min_aromatic"]:
                continue
            if hf.get("min_charge") and c["net_charge"] < hf["min_charge"]:
                continue
            if hf.get("max_acidic") and c["acidic_count"] > hf["max_acidic"]:
                continue
            if hf.get("max_abs_charge") and abs(c["net_charge"]) > hf["max_abs_charge"]:
                continue
            if hf.get("max_aromatic") and c["aromatic_hydrophobe_count"] > hf["max_aromatic"]:
                continue
            out.append(c)
        return out

    def build_llm_prompt(self, mode="cpp_like", n=20, extra_constraints=None):
        """Build an LLM prompt combining guideline context with classified blocks."""
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

        context = self.guideline_text or f"Mode: {config.get('display_name', mode)}"

        return (
            f"Design cyclic peptides following these guidelines:\n\n"
            f"{context}\n\n---\n\n"
            f"AVAILABLE BUILDING BLOCKS (classified by functional role):\n\n"
            f"{block_table}\n\n"
            f"GENERATION TASK:\n"
            f"- Mode: {config.get('display_name', mode)}\n"
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
        if "passive" in mode:
            return self._score_passive(block_ids, blocks)
        if "antibiotic" in mode or "amphipathic" in mode or "broad" in mode or "cationic_amphipathic" in mode:
            return self._score_antibiotic(block_ids, blocks)
        return self._score_cpp(block_ids, blocks)

    def _score_cpp(self, block_ids, blocks):
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
        score -= 2.0 * sum(1 for b in blocks if b.charge < 0)
        if n_aro > 4:
            score -= 2.0
        if sum(b.rot_total for b in blocks) > rs * 5:
            score -= 2.0
        return score

    def _score_antibiotic(self, block_ids, blocks):
        score = 0.0
        rs = len(blocks)
        n_primary = sum(1 for bid in block_ids if "primary_amine_cation" in self.block_to_bins.get(bid, set()))
        n_cationic = sum(1 for b in blocks if b.charge > 0)
        n_arg = sum(1 for bid in block_ids if "guanidinium_cation" in self.block_to_bins.get(bid, set()))
        score += 2.0 * min(n_primary, 4)
        score += 1.0 * min(n_arg, 1)
        if n_arg > 2:
            score -= 1.5
        n_bulky = sum(1 for bid in block_ids if "bulky_aliphatic" in self.block_to_bins.get(bid, set()))
        n_hydro = sum(1 for bid in block_ids if "all_hydrophobic" in self.block_to_bins.get(bid, set()))
        frac = n_hydro / rs if rs else 0
        if 0.35 <= frac <= 0.6:
            score += 2.0
        elif 0.25 <= frac < 0.35:
            score += 1.0
        if frac > 0.7:
            score -= 2.0
        score += 1.5 * min(n_bulky, 3)
        n_aro = sum(1 for bid in block_ids if "aromatic_hydrophobe" in self.block_to_bins.get(bid, set()))
        score += 1.0 * min(n_aro, 2)
        n_d = sum(1 for b in blocks if b.chirality == "D")
        score += 1.0 * min(n_d, 3)
        n_turn = sum(1 for bid in block_ids if "turn_inducer" in self.block_to_bins.get(bid, set()))
        score += 0.5 * min(n_turn, 2)
        if n_cationic >= 3 and n_bulky >= 2:
            score += 1.5
        score -= 2.0 * sum(1 for b in blocks if b.charge < 0)
        if rs > 10:
            score -= 1.0
        return score

    def _score_passive(self, block_ids, blocks):
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
        score -= 2.0 * abs(sum(b.charge for b in blocks))
        if rs > 9:
            score -= 2.0
        return score

    # ── Internal: annotation ─────────────────────────────────────────

    def _annotate_sequence(self, block_ids, mode="cpp_like"):
        blocks = [self.library[bid] for bid in block_ids]
        rs = len(blocks)
        config = self.modes.get(mode, {})
        thresholds = config.get("priority_thresholds", (10, 7))

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

        n_primary = sum(1 for bid in block_ids if "primary_amine_cation" in self.block_to_bins.get(bid, set()))
        n_bulky = sum(1 for bid in block_ids if "bulky_aliphatic" in self.block_to_bins.get(bid, set()))

        if "antibiotic" in mode or "amphipathic" in mode or "broad" in mode or "cationic_amphipathic" in mode:
            if n_primary >= 3 and n_bulky >= 2:
                mode_label = "cationic_membrane_disruptor"
            elif net_charge >= 2 and n_cationic >= 2:
                mode_label = "amphipathic_antimicrobial"
            else:
                mode_label = mode
        elif net_charge >= 2 and n_cationic >= 3:
            mode_label = "CPP_like_endocytic"
        elif abs(net_charge) <= 1:
            mode_label = "passive_diffusion"
        else:
            mode_label = mode

        score = self._score_sequence(block_ids, mode)
        if score >= thresholds[0]:
            priority = "high"
        elif score >= thresholds[1]:
            priority = "medium"
        else:
            priority = "low"

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
            "total_hbd": sum(b.sc_hbd for b in blocks),
            "total_hba": sum(b.sc_hba for b in blocks),
            "total_rotatable_bonds": sum(b.rot_total for b in blocks),
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

    p_gen = sub.add_parser("generate", help="Generate candidates")
    p_gen.add_argument("guideline", nargs="?", default=None,
                       help="Path to guideline MD file (omit for defaults)")
    p_gen.add_argument("-n", type=int, default=50, help="Number of candidates")
    p_gen.add_argument("--mode", default="cpp_like", help="Design mode")
    p_gen.add_argument("--ring-size", type=int, default=None)
    p_gen.add_argument("--seed", type=int, default=None)
    p_gen.add_argument("--min-score", type=float, default=None)
    p_gen.add_argument("-o", "--output", default=None, help="CSV output path")
    p_gen.add_argument("--backend", default=None)

    p_cls = sub.add_parser("classify", help="Show block classification")
    p_cls.add_argument("--backend", default=None)

    p_sc = sub.add_parser("score", help="Score a sequence")
    p_sc.add_argument("sequence", help="Dash-separated block IDs")
    p_sc.add_argument("--guideline", default=None)
    p_sc.add_argument("--mode", default="cpp_like")
    p_sc.add_argument("--backend", default=None)

    p_pr = sub.add_parser("prompt", help="Build LLM prompt")
    p_pr.add_argument("guideline", help="Path to guideline MD file")
    p_pr.add_argument("--mode", default="cpp_like")
    p_pr.add_argument("-n", type=int, default=20)
    p_pr.add_argument("--backend", default=None)

    args = parser.parse_args()
    if not args.command:
        parser.print_help()
        sys.exit(1)

    backend = getattr(args, "backend", None)

    if args.command == "classify":
        gen = GuidelineGenerator.from_params(backend_path=backend)
        print(gen.block_summary())

    elif args.command == "generate":
        if args.guideline:
            gen = GuidelineGenerator(args.guideline, backend_path=backend)
        else:
            gen = GuidelineGenerator.from_params(backend_path=backend)

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
        if args.guideline:
            gen = GuidelineGenerator(args.guideline, backend_path=backend)
        else:
            gen = GuidelineGenerator.from_params(backend_path=backend)
        ann = gen.annotate(args.sequence, mode=args.mode)
        print(f"Sequence: {ann['sequence']}")
        print(f"Score: {ann['score']:.1f} ({ann['priority']})")
        for k in ["ring_size", "net_charge", "cationic_count", "arg_like_count",
                   "aromatic_hydrophobe_count", "hydrophobic_fraction", "N_methyl_count",
                   "D_residue_count", "turn_residue_count", "acidic_count",
                   "total_hbd", "total_hba", "mode_label"]:
            print(f"  {k}: {ann[k]}")

    elif args.command == "prompt":
        gen = GuidelineGenerator(args.guideline, backend_path=backend)
        print(gen.build_llm_prompt(mode=args.mode, n=args.n))


if __name__ == "__main__":
    main()
