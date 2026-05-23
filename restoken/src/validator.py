"""ResToken sequence validator.

Validates LLM-output token sequences against the 400-block library
and user-specified physicochemical constraints.
"""

from dataclasses import dataclass, field
from typing import Optional

from restoken.src.library import BlockLibrary, Block


@dataclass
class ValidationResult:
    valid: bool
    sequence: list[str]
    checks: dict[str, bool] = field(default_factory=dict)
    errors: list[str] = field(default_factory=list)
    properties: dict = field(default_factory=dict)

    def summary(self) -> str:
        status = "PASS" if self.valid else "FAIL"
        lines = [f"[{status}] {'-'.join(self.sequence)}"]
        for check, passed in self.checks.items():
            mark = "OK" if passed else "FAIL"
            lines.append(f"  {mark}: {check}")
        if self.errors:
            for e in self.errors:
                lines.append(f"  ! {e}")
        if self.properties:
            props = ", ".join(f"{k}={v}" for k, v in self.properties.items())
            lines.append(f"  props: {props}")
        return "\n".join(lines)


# Backbone type ordering for adjacency compatibility
_BACKBONE_ORDER = {"alpha": 0, "beta": 1, "gamma": 2}


def parse_sequence(raw: str) -> list[str]:
    """Parse a token sequence string into a list of block IDs.

    Accepts formats:
      - Dash-separated: "A01-K03-N12-S05-E02-a07"
      - Space-separated: "A01 K03 N12 S05 E02 a07"
      - Comma-separated: "A01, K03, N12, S05, E02, a07"
      - Bracket-wrapped HELM-like: "[A01]-[K03]-[N12]"
    """
    raw = raw.strip()
    raw = raw.replace("[", "").replace("]", "")

    if "," in raw:
        tokens = [t.strip() for t in raw.split(",")]
    elif "-" in raw:
        tokens = [t.strip() for t in raw.split("-")]
    else:
        tokens = raw.split()

    return [t for t in tokens if t]


class SequenceValidator:
    """Validates ResToken sequences against library and constraints.

    Parameters
    ----------
    library : BlockLibrary, optional
        Block library instance. If None, loads the default v11 library.
    allowed_lengths : list[int], optional
        Allowed sequence lengths. Default: [6, 8, 10].
    max_hbd : int, optional
        Maximum total hydrogen bond donors. Default: no limit.
    max_hba : int, optional
        Maximum total hydrogen bond acceptors. Default: no limit.
    max_rot : int, optional
        Maximum total rotatable bonds. Default: no limit.
    target_charge : int or None, optional
        Required net charge. Default: no constraint.
    check_backbone_compat : bool, optional
        Whether to check backbone adjacency compatibility. Default: True.
    backbone_rules : str, optional
        Backbone compatibility rule set. "permissive" allows all transitions,
        "strict" disallows gamma→gamma adjacency. Default: "permissive".
    """

    def __init__(
        self,
        library: Optional[BlockLibrary] = None,
        allowed_lengths: Optional[list[int]] = None,
        max_hbd: Optional[int] = None,
        max_hba: Optional[int] = None,
        max_rot: Optional[int] = None,
        target_charge: Optional[int] = None,
        check_backbone_compat: bool = True,
        backbone_rules: str = "permissive",
    ):
        self.lib = library or BlockLibrary()
        self.allowed_lengths = allowed_lengths or [6, 8, 10]
        self.max_hbd = max_hbd
        self.max_hba = max_hba
        self.max_rot = max_rot
        self.target_charge = target_charge
        self.check_backbone_compat = check_backbone_compat
        self.backbone_rules = backbone_rules

    def validate(self, sequence: list[str] | str) -> ValidationResult:
        """Validate a single token sequence.

        Parameters
        ----------
        sequence : list[str] or str
            Either a list of block IDs or a raw string to parse.

        Returns
        -------
        ValidationResult
        """
        if isinstance(sequence, str):
            sequence = parse_sequence(sequence)

        result = ValidationResult(valid=True, sequence=sequence)
        blocks: list[Block] = []

        # 1. Token ID existence
        unknown = []
        for tid in sequence:
            if tid not in self.lib:
                unknown.append(tid)
            else:
                blocks.append(self.lib[tid])

        id_ok = len(unknown) == 0
        result.checks["id_existence"] = id_ok
        if not id_ok:
            result.errors.append(f"Unknown token IDs: {unknown}")
            result.valid = False

        # 2. Sequence length
        len_ok = len(sequence) in self.allowed_lengths
        result.checks["length"] = len_ok
        if not len_ok:
            result.errors.append(
                f"Length {len(sequence)} not in allowed {self.allowed_lengths}"
            )
            result.valid = False

        if not blocks:
            return result

        # Compute aggregate properties
        total_charge = sum(b.charge for b in blocks)
        total_hbd = sum(b.sc_hbd for b in blocks)
        total_hba = sum(b.sc_hba for b in blocks)
        total_rot = sum(b.rot_total for b in blocks)
        backbone_types = [b.mc_type for b in blocks]
        n_mods = [b.mc_nmod for b in blocks]

        result.properties = {
            "net_charge": total_charge,
            "total_hbd": total_hbd,
            "total_hba": total_hba,
            "total_rot": total_rot,
            "backbone_types": backbone_types,
            "n_nme": sum(1 for m in n_mods if m == "NME"),
            "n_ncy": sum(1 for m in n_mods if m == "NCY"),
        }

        # 3. Net charge
        if self.target_charge is not None:
            charge_ok = total_charge == self.target_charge
            result.checks["net_charge"] = charge_ok
            if not charge_ok:
                result.errors.append(
                    f"Net charge {total_charge} != target {self.target_charge}"
                )
                result.valid = False

        # 4. HBD limit
        if self.max_hbd is not None:
            hbd_ok = total_hbd <= self.max_hbd
            result.checks["max_hbd"] = hbd_ok
            if not hbd_ok:
                result.errors.append(
                    f"Total HBD {total_hbd} > max {self.max_hbd}"
                )
                result.valid = False

        # 5. HBA limit
        if self.max_hba is not None:
            hba_ok = total_hba <= self.max_hba
            result.checks["max_hba"] = hba_ok
            if not hba_ok:
                result.errors.append(
                    f"Total HBA {total_hba} > max {self.max_hba}"
                )
                result.valid = False

        # 6. Rotatable bonds limit
        if self.max_rot is not None:
            rot_ok = total_rot <= self.max_rot
            result.checks["max_rot"] = rot_ok
            if not rot_ok:
                result.errors.append(
                    f"Total rotatable bonds {total_rot} > max {self.max_rot}"
                )
                result.valid = False

        # 7. Backbone compatibility (cyclic — last→first also checked)
        if self.check_backbone_compat and len(blocks) >= 2:
            compat_ok, compat_msg = self._check_backbone_compat(blocks)
            result.checks["backbone_compat"] = compat_ok
            if not compat_ok:
                result.errors.append(compat_msg)
                result.valid = False

        return result

    def _check_backbone_compat(self, blocks: list[Block]) -> tuple[bool, str]:
        """Check backbone adjacency compatibility for a cyclic sequence."""
        n = len(blocks)
        for i in range(n):
            j = (i + 1) % n
            bi, bj = blocks[i], blocks[j]

            if self.backbone_rules == "strict":
                # Disallow gamma→gamma adjacency (ring strain)
                if bi.mc_type == "gamma" and bj.mc_type == "gamma":
                    return False, (
                        f"Pos {i}→{j}: gamma-gamma adjacency disallowed "
                        f"({bi.id}→{bj.id})"
                    )
                # Disallow consecutive NCY (N-cyclic, proline-like)
                if bi.mc_nmod == "NCY" and bj.mc_nmod == "NCY":
                    return False, (
                        f"Pos {i}→{j}: consecutive NCY blocks disallowed "
                        f"({bi.id}→{bj.id})"
                    )

        return True, ""

    def validate_batch(
        self, sequences: list[list[str] | str]
    ) -> list[ValidationResult]:
        """Validate a batch of sequences."""
        return [self.validate(s) for s in sequences]

    def batch_stats(self, results: list[ValidationResult]) -> dict:
        """Compute aggregate statistics from a batch of validation results."""
        n = len(results)
        if n == 0:
            return {}

        n_valid = sum(1 for r in results if r.valid)
        check_names = set()
        for r in results:
            check_names.update(r.checks.keys())

        per_check = {}
        for cn in sorted(check_names):
            passed = sum(1 for r in results if r.checks.get(cn, True))
            per_check[cn] = {"passed": passed, "total": n, "rate": passed / n}

        return {
            "total": n,
            "valid": n_valid,
            "validity_rate": n_valid / n,
            "per_check": per_check,
        }
