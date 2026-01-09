# -*- coding: utf-8 -*-
"""
test_niche_cases.py - HELIX DNA Storage System
Comprehensive tests for edge cases and niche scenarios
"""

from helix import HelixCodec
from analyzer import SequenceAnalyzer
import mapping
import differential
from rll_constraint import RLLCodec
from gc_balance import GCBalancer


def test_edge_cases():
    """Test edge cases and boundary conditions."""
    print("=" * 70)
    print("TESTING EDGE CASES")
    print("=" * 70)

    codec = HelixCodec(ell=3, epsilon=0.05)
    test_cases = []

    # 1. Empty and minimal inputs
    print("\n1. Minimal Inputs")
    print("-" * 70)

    cases = [
        ("00", "Empty-like (2 bits)"),
        ("0000", "Four zeros"),
        ("11", "Two ones"),
        ("01", "Single byte minimum"),
        ("10", "Alternating minimum"),
    ]

    for binary, desc in cases:
        try:
            dna = codec.encode(binary, verbose=False)
            decoded = codec.decode(dna, verbose=False)
            match = binary == decoded
            test_cases.append((desc, match))
            status = "PASS" if match else "FAIL"
            print(f"  {desc:30s} {status}")
            print(f"    Binary: {binary} -> DNA: {dna} -> Decoded: {decoded}")
        except Exception as e:
            test_cases.append((desc, False))
            print(f"  {desc:30s} X ERROR: {str(e)[:50]}")

    # 2. All zeros and all ones
    print("\n2. Homogeneous Patterns")
    print("-" * 70)

    for length in [8, 16, 24]:
        # All zeros
        binary = "0" * length
        try:
            dna = codec.encode(binary, verbose=False)
            decoded = codec.decode(dna, verbose=False)
            match = binary == decoded
            test_cases.append((f"All zeros ({length} bits)", match))
            status = "PASS" if match else "FAIL"
            print(f"  All zeros ({length} bits): {status}")
        except Exception as e:
            test_cases.append((f"All zeros ({length} bits)", False))
            print(f"  All zeros ({length} bits): X ERROR: {str(e)[:40]}")

        # All ones
        binary = "1" * length
        try:
            dna = codec.encode(binary, verbose=False)
            decoded = codec.decode(dna, verbose=False)
            match = binary == decoded
            test_cases.append((f"All ones ({length} bits)", match))
            status = "PASS" if match else "FAIL"
            print(f"  All ones ({length} bits): {status}")
        except Exception as e:
            test_cases.append((f"All ones ({length} bits)", False))
            print(f"  All ones ({length} bits): X ERROR: {str(e)[:40]}")

    # 3. Sequences that naturally contain escape patterns
    print("\n3. Natural Escape Patterns")
    print("-" * 70)

    # Binary that produces [3,2] in quaternary
    escape_patterns = [
        ("1110", "[3,2] in quaternary"),
        ("111010", "Multiple [3,2] patterns"),
        ("11101110", "Adjacent [3,2] patterns"),
    ]

    for binary, desc in escape_patterns:
        quat = mapping.binary_to_quaternary(binary)
        print(f"  {desc}: Binary {binary} -> Quaternary {quat}")
        try:
            dna = codec.encode(binary, verbose=False)
            decoded = codec.decode(dna, verbose=False)
            match = binary == decoded
            test_cases.append((desc, match))
            status = "PASS" if match else "FAIL"
            print(f"    Result: {status}")
        except Exception as e:
            test_cases.append((desc, False))
            print(f"    Result: X ERROR: {str(e)[:50]}")

    # 4. Sequences with maximum runlength
    print("\n4. Maximum Runlength Boundaries")
    print("-" * 70)

    # Test sequences with exactly ell consecutive symbols
    runlength_cases = [
        ("000000", "6 zeros (ell=3, should have 2 forbidden)"),
        ("001100110011", "Exactly at boundary"),
        ("00110011001100", "Multiple boundary cases"),
    ]

    for binary, desc in runlength_cases:
        try:
            dna = codec.encode(binary, verbose=False)
            decoded = codec.decode(dna, verbose=False)
            match = binary == decoded
            analyzer = SequenceAnalyzer(ell=3)
            max_run = analyzer.compute_max_runlength(dna)
            test_cases.append((desc, match))
            print(f"  {desc}")
            print(f"    DNA max runlength: {max_run} (limit: 3)")
            status = "PASS" if match else "FAIL"
            print(f"    Result: {status}")
        except Exception as e:
            test_cases.append((desc, False))
            print(f"  {desc}: X ERROR: {str(e)[:50]}")

    # 5. Long sequences
    print("\n5. Long Sequences")
    print("-" * 70)

    for length in [64, 128, 256]:
        # Alternating pattern
        binary = "01" * (length // 2)
        try:
            dna = codec.encode(binary, verbose=False)
            decoded = codec.decode(dna, verbose=False)
            match = binary == decoded
            test_cases.append((f"Alternating {length} bits", match))
            efficiency = len(binary) / (len(dna) * 2) * 100
            status = "PASS" if match else "FAIL"
            print(f"  {length} bits alternating: {status}")
            print(f"    DNA length: {len(dna)} bp, Efficiency: {efficiency:.1f}%")
        except Exception as e:
            test_cases.append((f"Alternating {length} bits", False))
            print(f"  {length} bits alternating: X ERROR: {str(e)[:40]}")

    # 6. Sequences that stress differential encoding
    print("\n6. Differential Encoding Stress Tests")
    print("-" * 70)

    diff_stress = [
        ("00000000", "All same quaternary"),
        ("01010101", "Alternating in binary"),
        ("00110011", "Alternating pairs"),
        ("11100100", "Random mix"),
    ]

    for binary, desc in diff_stress:
        quat = mapping.binary_to_quaternary(binary)
        diff = differential.differential_encode(quat)
        print(f"  {desc}: Quat {quat} -> Diff {diff}")
        try:
            dna = codec.encode(binary, verbose=False)
            decoded = codec.decode(dna, verbose=False)
            match = binary == decoded
            test_cases.append((desc, match))
            status = "PASS" if match else "FAIL"
            print(f"    Result: {status}")
        except Exception as e:
            test_cases.append((desc, False))
            print(f"    Result: X ERROR: {str(e)[:50]}")

    return test_cases


def test_constraint_boundaries():
    """Test GC-balance and runlength at boundary conditions."""
    print("\n" + "=" * 70)
    print("TESTING CONSTRAINT BOUNDARIES")
    print("=" * 70)

    test_cases = []

    # Test with different constraint parameters
    print("\n1. Strict Constraints (ell=2, epsilon=0.03)")
    print("-" * 70)

    codec_strict = HelixCodec(ell=2, epsilon=0.03)
    test_data = ["11010011", "10101010", "11110000"]

    for binary in test_data:
        try:
            dna = codec_strict.encode(binary, verbose=False)
            decoded = codec_strict.decode(dna, verbose=False)
            match = binary == decoded

            analyzer = SequenceAnalyzer(ell=2, epsilon=0.03)
            analysis = analyzer.analyze_dna(dna)

            test_cases.append((f"Strict: {binary}", match))
            status = "PASS" if match else "FAIL"
            print(f"  Binary: {binary}")
            print(f"    GC-content: {analysis['gc_content']:.2%} (target: 47-53%)")
            print(f"    Max runlength: {analysis['max_runlength']} (limit: 2)")
            print(f"    Roundtrip: {status}")
        except Exception as e:
            test_cases.append((f"Strict: {binary}", False))
            print(f"  Binary: {binary}: X ERROR: {str(e)[:50]}")

    # Test with relaxed constraints
    print("\n2. Relaxed Constraints (ell=4, epsilon=0.1)")
    print("-" * 70)

    codec_relaxed = HelixCodec(ell=4, epsilon=0.1)

    for binary in test_data:
        try:
            dna = codec_relaxed.encode(binary, verbose=False)
            decoded = codec_relaxed.decode(dna, verbose=False)
            match = binary == decoded

            analyzer = SequenceAnalyzer(ell=4, epsilon=0.1)
            analysis = analyzer.analyze_dna(dna)

            test_cases.append((f"Relaxed: {binary}", match))
            status = "PASS" if match else "FAIL"
            print(f"  Binary: {binary}")
            print(f"    GC-content: {analysis['gc_content']:.2%} (target: 40-60%)")
            print(f"    Max runlength: {analysis['max_runlength']} (limit: 4)")
            print(f"    Roundtrip: {status}")
        except Exception as e:
            test_cases.append((f"Relaxed: {binary}", False))
            print(f"  Binary: {binary}: X ERROR: {str(e)[:50]}")

    return test_cases


def test_special_sequences():
    """Test special bit patterns and sequences."""
    print("\n" + "=" * 70)
    print("TESTING SPECIAL SEQUENCES")
    print("=" * 70)

    codec = HelixCodec(ell=3, epsilon=0.05)
    test_cases = []

    # 1. Sequences that produce specific quaternary patterns
    print("\n1. Specific Quaternary Patterns")
    print("-" * 70)

    special = [
        ("00001111", "Step pattern"),
        ("11110000", "Reverse step"),
        ("01100110", "Symmetric pattern"),
        ("10011001", "Inverse symmetric"),
    ]

    for binary, desc in special:
        quat = mapping.binary_to_quaternary(binary)
        print(f"  {desc}: {binary} -> {quat}")
        try:
            dna = codec.encode(binary, verbose=False)
            decoded = codec.decode(dna, verbose=False)
            match = binary == decoded
            test_cases.append((desc, match))
            print(f"    DNA: {dna}")
            status = "PASS" if match else "FAIL"
            print(f"    Result: {status}")
        except Exception as e:
            test_cases.append((desc, False))
            print(f"    Result: X ERROR: {str(e)[:50]}")

    # 2. Power of 2 lengths
    print("\n2. Power of 2 Bit Lengths")
    print("-" * 70)

    for exp in [2, 3, 4, 5, 6]:
        length = 2 ** exp
        binary = "10" * (length // 2)
        try:
            dna = codec.encode(binary, verbose=False)
            decoded = codec.decode(dna, verbose=False)
            match = binary == decoded
            test_cases.append((f"{length} bits", match))
            status = "PASS" if match else "FAIL"
            print(f"  {length} bits: {status}")
        except Exception as e:
            test_cases.append((f"{length} bits", False))
            print(f"  {length} bits: X ERROR: {str(e)[:40]}")

    # 3. Sequences with many forbidden substrings after differential encoding
    print("\n3. Multiple Forbidden Substrings")
    print("-" * 70)

    # These should produce many zeros after differential encoding
    multi_forbidden = [
        ("00000000", "8 zeros"),
        ("0000000000000000", "16 zeros"),
        ("00001010000010100000", "Multiple zero runs"),
    ]

    for binary, desc in multi_forbidden:
        quat = mapping.binary_to_quaternary(binary)
        diff = differential.differential_encode(quat)
        rll_codec = RLLCodec(ell=3)
        forbidden_count = sum(1 for i in range(len(diff) - 2)
                             if all(diff[i+j] == 0 for j in range(3)))

        print(f"  {desc}: {forbidden_count} forbidden substrings in differential")
        try:
            dna = codec.encode(binary, verbose=False)
            decoded = codec.decode(dna, verbose=False)
            match = binary == decoded
            test_cases.append((desc, match))
            status = "PASS" if match else "FAIL"
            print(f"    Result: {status}")
        except Exception as e:
            test_cases.append((desc, False))
            print(f"    Result: X ERROR: {str(e)[:50]}")

    return test_cases


def test_codec_combinations():
    """Test different codec parameter combinations."""
    print("\n" + "=" * 70)
    print("TESTING CODEC PARAMETER COMBINATIONS")
    print("=" * 70)

    test_cases = []
    binary = "1101001110101010"

    # Test various parameter combinations
    params = [
        (2, 0.05, True, "Strict RLL, standard GC, EC on"),
        (3, 0.05, True, "Standard RLL, standard GC, EC on"),
        (4, 0.05, True, "Relaxed RLL, standard GC, EC on"),
        (3, 0.03, True, "Standard RLL, strict GC, EC on"),
        (3, 0.10, True, "Standard RLL, relaxed GC, EC on"),
        (3, 0.05, False, "Standard RLL, standard GC, EC off"),
        (2, 0.03, True, "Strict RLL, strict GC, EC on"),
        (4, 0.10, False, "Relaxed RLL, relaxed GC, EC off"),
    ]

    for ell, epsilon, use_ec, desc in params:
        print(f"\n{desc}")
        print(f"  ell={ell}, epsilon={epsilon}, EC={use_ec}")

        try:
            codec = HelixCodec(ell=ell, epsilon=epsilon, use_error_correction=use_ec)
            dna = codec.encode(binary, verbose=False)
            decoded = codec.decode(dna, verbose=False)
            match = binary == decoded

            analyzer = SequenceAnalyzer(ell=ell, epsilon=epsilon)
            analysis = analyzer.analyze_dna(dna)

            test_cases.append((desc, match))
            print(f"  DNA: {dna} ({len(dna)} bp)")
            print(f"  GC: {analysis['gc_content']:.2%}, Max run: {analysis['max_runlength']}")
            status = "PASS" if match else "FAIL"
            print(f"  Result: {status}")
        except Exception as e:
            test_cases.append((desc, False))
            print(f"  Result: X ERROR: {str(e)[:50]}")

    return test_cases


def test_error_conditions():
    """Test error handling and invalid inputs."""
    print("\n" + "=" * 70)
    print("TESTING ERROR CONDITIONS")
    print("=" * 70)

    codec = HelixCodec(ell=3, epsilon=0.05)
    test_cases = []

    # Test invalid inputs
    print("\n1. Invalid Binary Strings")
    print("-" * 70)

    invalid_inputs = [
        ("", "Empty string"),
        ("2", "Invalid digit (2)"),
        ("abc", "Non-numeric"),
        ("101 01", "Contains space"),
        ("1.01", "Contains decimal"),
    ]

    for invalid, desc in invalid_inputs:
        try:
            dna = codec.encode(invalid, verbose=False)
            test_cases.append((desc, False))
            print(f"  {desc:30s} X FAIL (should have raised error)")
        except Exception as e:
            test_cases.append((desc, True))
            print(f"  {desc:30s} OK PASS (caught error: {type(e).__name__})")

    # Test odd-length binary (should be padded)
    print("\n2. Odd-Length Binary Strings")
    print("-" * 70)

    odd_lengths = ["1", "101", "10101", "1010101"]

    for binary in odd_lengths:
        try:
            dna = codec.encode(binary, verbose=False)
            decoded = codec.decode(dna, verbose=False)
            # Odd-length gets padded with leading zero
            expected = "0" + binary if len(binary) % 2 == 1 else binary
            match = decoded == binary or decoded == expected
            test_cases.append((f"Odd length {len(binary)}", match))
            print(f"  Length {len(binary)}: {binary} -> decoded: {decoded}")
            status = "PASS" if match else "FAIL"
            print(f"    Result: {status}")
        except Exception as e:
            test_cases.append((f"Odd length {len(binary)}", False))
            print(f"  Length {len(binary)}: X ERROR: {str(e)[:50]}")

    return test_cases


def print_summary(all_test_cases):
    """Print summary of all tests."""
    print("\n" + "=" * 70)
    print("TEST SUMMARY")
    print("=" * 70)

    total = len(all_test_cases)
    passed = sum(1 for _, result in all_test_cases if result)
    failed = total - passed
    pass_rate = (passed / total * 100) if total > 0 else 0

    print(f"\nTotal tests:  {total}")
    print(f"Passed:       {passed} OK")
    print(f"Failed:       {failed} X")
    print(f"Pass rate:    {pass_rate:.1f}%")

    if failed > 0:
        print(f"\nFailed tests:")
        for name, result in all_test_cases:
            if not result:
                print(f"  X {name}")

    print("\n" + "=" * 70)


def main():
    """Run all niche case tests."""
    print("=" * 70)
    print("HELIX DNA Storage System - Niche Case Testing")
    print("=" * 70)
    print("\nTesting edge cases, boundary conditions, and special scenarios...")

    all_test_cases = []

    # Run all test suites
    all_test_cases.extend(test_edge_cases())
    all_test_cases.extend(test_constraint_boundaries())
    all_test_cases.extend(test_special_sequences())
    all_test_cases.extend(test_codec_combinations())
    all_test_cases.extend(test_error_conditions())

    # Print summary
    print_summary(all_test_cases)


if __name__ == "__main__":
    main()
