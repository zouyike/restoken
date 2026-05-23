"""Build filled prompts from templates + library data.

Usage:
    python build_prompts.py --template restoken_zeroshot --length 6 --n 50
    python build_prompts.py --template smiles_zeroshot --length 6 --n 50
    python build_prompts.py --template helm_zeroshot --length 6 --n 50
    python build_prompts.py --template property_constrained --length 6 --n 50 --profile permeable
"""

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from restoken.src.library import BlockLibrary

PROMPT_DIR = Path(__file__).resolve().parent.parent / "input" / "prompts"
OUTPUT_DIR = Path(__file__).resolve().parent.parent / "output" / "filled_prompts"

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


def build_restoken_prompt(lib: BlockLibrary, n: int, length: int, profile: str) -> str:
    template = (PROMPT_DIR / "restoken_zeroshot.txt").read_text()
    constraints = CONSTRAINT_PROFILES[profile].format(length=length)
    return template.format(
        N_BLOCKS=lib.size,
        DICTIONARY_TABLE=lib.llm_dict_text(),
        N_SEQUENCES=n,
        SEQ_LENGTH=length,
        CONSTRAINTS=constraints,
    )


def build_smiles_prompt(lib: BlockLibrary, n: int, length: int, profile: str) -> str:
    template = (PROMPT_DIR / "smiles_zeroshot.txt").read_text()

    lines = []
    for bid in sorted(lib.all_ids):
        b = lib[bid]
        lines.append(f"{bid} | {b.aa_smiles} | class={b.aa_class}, chiral={b.chirality}, "
                      f"charge={b.charge}, backbone={b.mc_type}")

    smiles_table = "\n".join(lines)
    constraints = CONSTRAINT_PROFILES[profile].format(length=length)

    return template.format(
        N_BLOCKS=lib.size,
        SMILES_TABLE=smiles_table,
        N_SEQUENCES=n,
        SEQ_LENGTH=length,
        CONSTRAINTS=constraints,
    )


def build_helm_prompt(lib: BlockLibrary, n: int, length: int, profile: str) -> str:
    template = (PROMPT_DIR / "helm_zeroshot.txt").read_text()

    lines = []
    for bid in sorted(lib.all_ids):
        b = lib[bid]
        lines.append(f"{b.helm} | class={b.aa_class}, chiral={b.chirality}, "
                      f"charge={b.charge_label}, backbone={b.mc_type}, "
                      f"bulk={b.sc_bulk}")

    helm_table = "\n".join(lines)
    constraints = CONSTRAINT_PROFILES[profile].format(length=length)

    return template.format(
        N_BLOCKS=lib.size,
        HELM_TABLE=helm_table,
        N_SEQUENCES=n,
        SEQ_LENGTH=length,
        CONSTRAINTS=constraints,
    )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--template", required=True,
                        choices=["restoken_zeroshot", "smiles_zeroshot", "helm_zeroshot"])
    parser.add_argument("--length", type=int, default=6)
    parser.add_argument("--n", type=int, default=50)
    parser.add_argument("--profile", default="unconstrained",
                        choices=list(CONSTRAINT_PROFILES.keys()))
    parser.add_argument("--output", default=None)
    args = parser.parse_args()

    lib = BlockLibrary()

    builders = {
        "restoken_zeroshot": build_restoken_prompt,
        "smiles_zeroshot": build_smiles_prompt,
        "helm_zeroshot": build_helm_prompt,
    }

    prompt = builders[args.template](lib, args.n, args.length, args.profile)

    if args.output:
        out = Path(args.output)
    else:
        OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
        out = OUTPUT_DIR / f"{args.template}_{args.profile}_len{args.length}_n{args.n}.txt"

    out.write_text(prompt)
    print(f"Saved: {out} ({len(prompt)} chars, ~{len(prompt)//4} tokens)")


if __name__ == "__main__":
    main()
