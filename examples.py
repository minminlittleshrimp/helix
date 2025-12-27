"""
examples.py - HELIX DNA Storage System
Demonstration scripts and usage examples
"""

import helix
import analyzer
import mapping

def example_basic_encoding():
    """Basic encoding and decoding example."""
    print("\n" + "=" * 70)
    print("Example 1: Basic Encoding and Decoding")
    print("=" * 70)

    codec = helix.HelixCodec(ell=3, epsilon=0.05)

    # Encode a simple message
    binary_data = "11010011"
    print(f"\nOriginal binary: {binary_data}")

    dna = codec.encode(binary_data, verbose=True)
    print(f"\nEncoded DNA: {dna}")

    decoded = codec.decode(dna, verbose=True)
    print(f"\nDecoded binary: {decoded}")
    print(f"Match: {binary_data == decoded}")


def example_text_encoding():
    """Encode text message as DNA."""
    print("\n" + "=" * 70)
    print("Example 2: Text Message Encoding")
    print("=" * 70)

    codec = helix.HelixCodec(ell=3, epsilon=0.05)

    # Convert text to binary
    message = "HELIX"
    binary_data = ''.join(format(ord(char), '08b') for char in message)

    print(f"\nOriginal message: {message}")
    print(f"Binary representation: {binary_data}")
    print(f"Length: {len(binary_data)} bits")

    # Encode to DNA
    dna = codec.encode(binary_data, verbose=False)
    print(f"\nEncoded DNA: {dna}")
    print(f"DNA length: {len(dna)} bp")

    # Analyze DNA
    seq_analyzer = analyzer.SequenceAnalyzer(ell=3, epsilon=0.05)
    analysis = seq_analyzer.analyze_dna(dna)
    seq_analyzer.print_analysis(analysis)

    # Decode back
    decoded_binary = codec.decode(dna, verbose=False)
    decoded_message = ''.join(
        chr(int(decoded_binary[i:i+8], 2))
        for i in range(0, len(decoded_binary), 8)
    )

    print(f"\nDecoded message: {decoded_message}")
    print(f"Match: {message == decoded_message}")


def example_constraint_comparison():
    """Compare sequences with different constraints."""
    print("\n" + "=" * 70)
    print("Example 3: Constraint Comparison")
    print("=" * 70)

    binary_data = "1101001110101010"

    # Test with different parameters
    configs = [
        (2, 0.05, "Strict runlength"),
        (3, 0.05, "Standard"),
        (4, 0.05, "Relaxed runlength"),
        (3, 0.10, "Relaxed GC-content"),
    ]

    print(f"\nInput binary: {binary_data} ({len(binary_data)} bits)")
    print("\nComparing different configurations:")
    print("-" * 70)

    for ell, epsilon, description in configs:
        codec = helix.HelixCodec(ell=ell, epsilon=epsilon, use_error_correction=False)
        dna = codec.encode(binary_data, verbose=False)

        seq_analyzer = analyzer.SequenceAnalyzer(ell=ell, epsilon=epsilon)
        analysis = seq_analyzer.analyze_dna(dna)

        print(f"\n{description} (ell={ell}, epsilon={epsilon}):")
        print(f"  DNA: {dna}")
        print(f"  Length: {len(dna)} bp")
        print(f"  GC-content: {analysis['gc_content']:.2%}")
        print(f"  Max runlength: {analysis['max_runlength']}")
        print(f"  Efficiency: {len(binary_data) / (len(dna) * 2):.2%}")

        constraints = seq_analyzer.validate_constraints(dna)
        all_pass = all(constraints.values())
        print(f"  Constraints: {'PASS' if all_pass else 'FAIL'}")


def example_batch_processing():
    """Batch encode multiple sequences."""
    print("\n" + "=" * 70)
    print("Example 4: Batch Processing")
    print("=" * 70)

    codec = helix.HelixCodec(ell=3, epsilon=0.05)

    # Multiple binary sequences
    sequences = [
        "1010",
        "11001100",
        "111000111",
        "10101010101010",
    ]

    print("\nEncoding batch of sequences:")
    print("-" * 70)

    results = []
    for i, binary in enumerate(sequences, 1):
        dna = codec.encode(binary, verbose=False)
        decoded = codec.decode(dna, verbose=False)
        match = (binary == decoded)

        results.append({
            'index': i,
            'input': binary,
            'dna': dna,
            'decoded': decoded,
            'match': match,
            'compression_ratio': len(binary) / (len(dna) * 2)
        })

        print(f"\nSequence {i}:")
        print(f"  Input:  {binary} ({len(binary)} bits)")
        print(f"  DNA:    {dna} ({len(dna)} bp)")
        print(f"  Verify: {'PASS' if match else 'FAIL'}")
        print(f"  Efficiency: {results[-1]['compression_ratio']:.2%}")

    # Summary statistics
    print("\n" + "-" * 70)
    print("Batch Summary:")
    total_sequences = len(results)
    successful = sum(1 for r in results if r['match'])
    avg_efficiency = sum(r['compression_ratio'] for r in results) / len(results)

    print(f"  Total sequences: {total_sequences}")
    print(f"  Successful: {successful}")
    print(f"  Success rate: {successful/total_sequences*100:.1f}%")
    print(f"  Average efficiency: {avg_efficiency:.2%}")


def example_error_correction():
    """Demonstrate error detection."""
    print("\n" + "=" * 70)
    print("Example 5: Error Detection")
    print("=" * 70)

    codec = helix.HelixCodec(ell=3, epsilon=0.05, use_error_correction=True)

    binary_data = "11010011"
    print(f"\nOriginal binary: {binary_data}")

    # Encode
    dna = codec.encode(binary_data, verbose=False)
    print(f"Encoded DNA: {dna}")

    # Simulate error by mutating one nucleotide
    dna_list = list(dna)
    if len(dna_list) > 5:
        original_nucleotide = dna_list[5]
        # Mutate to different nucleotide
        mutations = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
        dna_list[5] = mutations[original_nucleotide]
        corrupted_dna = ''.join(dna_list)

        print(f"\nSimulated error at position 5:")
        print(f"  Original:  {dna}")
        print(f"  Corrupted: {corrupted_dna}")
        print(f"  Changed:   {original_nucleotide} -> {dna_list[5]}")

        # Convert to quaternary for error detection
        quaternary = mapping.dna_to_quaternary(corrupted_dna)

        # Check syndrome (simplified)
        ec = codec.error_corrector
        syndrome = ec.compute_syndrome(quaternary)
        checksum = ec.compute_checksum(quaternary)

        print(f"\nError detection:")
        print(f"  Syndrome: {syndrome}")
        print(f"  Checksum: {checksum}")
        print(f"  Note: Full error correction requires additional implementation")


def example_direct_conversion():
    """Direct conversions between representations."""
    print("\n" + "=" * 70)
    print("Example 6: Direct Format Conversions")
    print("=" * 70)

    # Binary -> DNA
    binary = "11010011"
    print(f"\nBinary: {binary}")

    quaternary = mapping.binary_to_quaternary(binary)
    print(f"Quaternary: {quaternary}")

    dna = mapping.quaternary_to_dna(quaternary)
    print(f"DNA: {dna}")

    # Reverse
    print("\nReverse conversion:")
    quaternary_back = mapping.dna_to_quaternary(dna)
    print(f"Quaternary: {quaternary_back}")

    binary_back = mapping.quaternary_to_binary(quaternary_back)
    print(f"Binary: {binary_back}")

    print(f"\nRoundtrip match: {binary == binary_back}")


def example_analysis_tools():
    """Demonstrate sequence analysis capabilities."""
    print("\n" + "=" * 70)
    print("Example 7: Sequence Analysis Tools")
    print("=" * 70)

    seq_analyzer = analyzer.SequenceAnalyzer(ell=3, epsilon=0.05)

    test_sequences = [
        ("Balanced", "ATCGATCG"),
        ("GC-rich", "GCGCGCGC"),
        ("AT-rich", "ATATATAT"),
        ("Long run", "AAAATTTCCCGGG"),
    ]

    for name, dna in test_sequences:
        print(f"\n{name}: {dna}")
        analysis = seq_analyzer.analyze_dna(dna)

        print(f"  Length: {analysis['length']}")
        print(f"  GC-content: {analysis['gc_content']:.2%}")
        print(f"  Max runlength: {analysis['max_runlength']}")
        print(f"  GC-balanced: {analysis['gc_balanced']}")
        print(f"  Runlength OK: {analysis['runlength_ok']}")

        if analysis['homopolymer_runs']:
            print(f"  Homopolymer runs:")
            for run in analysis['homopolymer_runs']:
                print(f"    {run['nucleotide']} x {run['length']}")


def main():
    """Run all examples."""
    print("=" * 70)
    print("HELIX - DNA Storage System Examples")
    print("=" * 70)

    examples = [
        ("Basic Encoding", example_basic_encoding),
        ("Text Encoding", example_text_encoding),
        ("Constraint Comparison", example_constraint_comparison),
        ("Batch Processing", example_batch_processing),
        ("Error Detection", example_error_correction),
        ("Direct Conversion", example_direct_conversion),
        ("Analysis Tools", example_analysis_tools),
    ]

    print("\nAvailable examples:")
    for i, (name, _) in enumerate(examples, 1):
        print(f"  {i}. {name}")

    print("\nRunning all examples...")

    for name, example_func in examples:
        try:
            example_func()
        except Exception as e:
            print(f"\nError in {name}: {e}")
            import traceback
            traceback.print_exc()

    print("\n" + "=" * 70)
    print("All examples completed")
    print("=" * 70)


if __name__ == "__main__":
    main()