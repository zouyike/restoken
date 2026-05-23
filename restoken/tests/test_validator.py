"""Unit tests for ResToken validator and reconstructor.

Generates 100 legal + 100 illegal sequences and validates each.
Also validates all 400 block SMILES via RDKit round-trip.
"""

import random
import sys
sys.path.insert(0, "/scratch/genesis/NCAA_tokenization")

from restoken.src.library import BlockLibrary
from restoken.src.validator import SequenceValidator, parse_sequence
from restoken.src.reconstructor import SMILESReconstructor


def make_legal_sequence(lib: BlockLibrary, length: int = 6) -> list[str]:
    """Generate a random sequence guaranteed to be format-legal (valid IDs, correct length)."""
    all_ids = list(lib.all_ids)
    return random.sample(all_ids, length)


def make_illegal_sequence_bad_id(lib: BlockLibrary, length: int = 6) -> list[str]:
    """Sequence with at least one fake ID."""
    all_ids = list(lib.all_ids)
    seq = random.sample(all_ids, length - 1)
    fake_ids = ["FAKE1", "ZZ99", "X00", "BOGUS", "AA999", "??"]
    seq.insert(random.randint(0, len(seq)), random.choice(fake_ids))
    return seq


def make_illegal_sequence_bad_length(lib: BlockLibrary) -> list[str]:
    """Sequence with wrong length (not 6, 8, or 10)."""
    all_ids = list(lib.all_ids)
    bad_lengths = [1, 2, 3, 4, 5, 7, 9, 11, 12, 15]
    length = random.choice(bad_lengths)
    return [random.choice(all_ids) for _ in range(length)]


def test_parse_formats():
    """Test all supported input formats."""
    print("=== Parse Format Tests ===")

    cases = [
        ("A01-K03-N12", ["A01", "K03", "N12"]),
        ("A01, K03, N12", ["A01", "K03", "N12"]),
        ("A01 K03 N12", ["A01", "K03", "N12"]),
        ("[A01]-[K03]-[N12]", ["A01", "K03", "N12"]),
        ("  A01 - K03 - N12  ", ["A01", "K03", "N12"]),
    ]

    for raw, expected in cases:
        result = parse_sequence(raw)
        assert result == expected, f"FAIL: parse('{raw}') = {result}, expected {expected}"
        print(f"  OK: '{raw}' -> {result}")

    print(f"  All {len(cases)} parse tests passed.\n")


def test_100_legal(lib: BlockLibrary, validator: SequenceValidator):
    """Generate and validate 100 legal sequences."""
    print("=== 100 Legal Sequences ===")
    passed = 0
    failed_details = []

    for i in range(100):
        length = random.choice([6, 8, 10])
        seq = make_legal_sequence(lib, length)
        result = validator.validate(seq)

        # For unconstrained validator, all should pass id_existence + length
        id_ok = result.checks.get("id_existence", False)
        len_ok = result.checks.get("length", False)

        if id_ok and len_ok:
            passed += 1
        else:
            failed_details.append((i, seq, result.errors))

    print(f"  Passed (id + length): {passed}/100")
    if failed_details:
        for idx, seq, errors in failed_details[:5]:
            print(f"  FAIL #{idx}: {'-'.join(seq)} -> {errors}")

    assert passed == 100, f"Expected 100/100 legal sequences to pass, got {passed}"
    print("  All 100 legal sequences passed.\n")


def test_100_illegal(lib: BlockLibrary, validator: SequenceValidator):
    """Generate and validate 100 illegal sequences (mix of failure types)."""
    print("=== 100 Illegal Sequences ===")
    correctly_rejected = 0

    for i in range(100):
        if i < 50:
            # Bad ID
            seq = make_illegal_sequence_bad_id(lib)
            result = validator.validate(seq)
            if not result.checks.get("id_existence", True):
                correctly_rejected += 1
        else:
            # Bad length
            seq = make_illegal_sequence_bad_length(lib)
            result = validator.validate(seq)
            if not result.checks.get("length", True):
                correctly_rejected += 1

    print(f"  Correctly rejected: {correctly_rejected}/100")
    assert correctly_rejected == 100, f"Expected 100/100 rejections, got {correctly_rejected}"
    print("  All 100 illegal sequences correctly rejected.\n")


def test_constraint_checks(lib: BlockLibrary):
    """Test specific constraint violations."""
    print("=== Constraint Violation Tests ===")

    # Charge constraint
    sv = SequenceValidator(library=lib, target_charge=0)
    neg_blocks = [b.id for b in lib.blocks_by_property(charge_label="neg")]
    if len(neg_blocks) >= 6:
        seq = neg_blocks[:6]
        r = sv.validate(seq)
        assert not r.checks.get("net_charge", True), "Should fail charge check"
        print(f"  OK: all-negative sequence correctly fails charge=0 check")

    # HBD constraint
    sv2 = SequenceValidator(library=lib, max_hbd=0)
    hbd_blocks = [b for b in lib.blocks_by_property() if b.sc_hbd > 0]
    if len(hbd_blocks) >= 6:
        seq = [b.id for b in hbd_blocks[:6]]
        r = sv2.validate(seq)
        assert not r.checks.get("max_hbd", True), "Should fail HBD check"
        print(f"  OK: high-HBD sequence correctly fails max_hbd=0 check")

    # Backbone strict mode: gamma-gamma
    sv3 = SequenceValidator(library=lib, backbone_rules="strict")
    gamma_blocks = [b.id for b in lib.blocks_by_property(mc_type="gamma")]
    if len(gamma_blocks) >= 6:
        seq = gamma_blocks[:6]
        r = sv3.validate(seq)
        assert not r.checks.get("backbone_compat", True), "Should fail backbone check"
        print(f"  OK: all-gamma sequence correctly fails strict backbone check")

    print("  All constraint tests passed.\n")


def test_smiles_roundtrip(lib: BlockLibrary):
    """Validate all 400 block SMILES via RDKit."""
    print("=== SMILES Round-Trip (400 blocks) ===")
    from rdkit import Chem

    recon = SMILESReconstructor(library=lib)
    all_ids = sorted(lib.all_ids)
    result = recon.validate_smiles(all_ids)

    print(f"  Total: {result['total']}, Valid: {result['valid']}, Invalid: {result['invalid']}")

    if result["invalid"] > 0:
        for d in result["details"]:
            if not d["valid"]:
                print(f"  INVALID: {d['id']} -> {d['smiles']}")

    assert result["invalid"] == 0, f"{result['invalid']} blocks have invalid SMILES"
    print("  All 400 blocks validated.\n")


def test_reconstruction_properties(lib: BlockLibrary):
    """Test that reconstruction returns correct properties."""
    print("=== Reconstruction Property Tests ===")
    recon = SMILESReconstructor(library=lib)

    seq = ["A01", "K03", "N12"]
    result = recon.reconstruct(seq)

    assert len(result["residue_smiles"]) == 3
    assert len(result["sidechain_smiles"]) == 3
    assert len(result["properties"]) == 3
    assert result["block_ids"] == seq
    assert result["properties"][0]["id"] == "A01"
    assert result["properties"][0]["class"] == "I"
    print("  OK: reconstruction returns correct structure and properties")

    # Test error on unknown ID
    try:
        recon.reconstruct(["FAKE"])
        assert False, "Should have raised KeyError"
    except KeyError:
        print("  OK: KeyError raised for unknown block ID")

    print("  All reconstruction tests passed.\n")


if __name__ == "__main__":
    random.seed(42)
    lib = BlockLibrary()

    # Unconstrained validator (only checks ID + length)
    sv_basic = SequenceValidator(library=lib)

    test_parse_formats()
    test_100_legal(lib, sv_basic)
    test_100_illegal(lib, sv_basic)
    test_constraint_checks(lib)
    test_smiles_roundtrip(lib)
    test_reconstruction_properties(lib)

    print("=" * 50)
    print("ALL TESTS PASSED")
    print("=" * 50)
