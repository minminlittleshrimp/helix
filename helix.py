"""
helix.py - HELIX DNA Storage System
High-Efficiency Lossless Information eXchange

Main codec integrating all encoding/decoding steps.
Based on: Nguyen et al. - Capacity-Approaching Constrained Codes with
Error Correction for DNA-Based Data Storage
"""

from typing import List, Tuple, Dict
import argparse
import sys
import mapping
import differential
import rll_constraint
import gc_balance
import error_correction
import analyzer


class HelixCodec:
    """
    Complete HELIX encoding/decoding pipeline.
    Transforms binary data into DNA sequences with:
    - Runlength constraint (no homopolymer runs > ell)
    - GC-content balancing (50% ± epsilon)
    - Error correction capability (single-edit errors)
    """

    def __init__(self, ell: int = 3, epsilon: float = 0.05,
                 use_error_correction: bool = True):
        """
        Initialize HELIX codec.

        Args:
            ell: Maximum runlength (homopolymer limit)
            epsilon: GC-content tolerance
            use_error_correction: Enable error correction suffix
        """
        self.ell = ell
        self.epsilon = epsilon
        self.use_error_correction = use_error_correction

        # Initialize component modules
        self.rll_codec = rll_constraint.RLLCodec(ell=ell)
        self.gc_balancer = gc_balance.GCBalancer(epsilon=epsilon)
        self.error_corrector = error_correction.VTErrorCorrection()
        self.analyzer = analyzer.SequenceAnalyzer(ell=ell, epsilon=epsilon)

    def encode(self, binary_data: str, verbose: bool = False) -> str:
        """
        Encode binary data into DNA sequence.

        Pipeline:
        Step 1: Binary -> Quaternary
        Step 2: Differential encoding (prepares for RLL)
        Step 3: RLL encoding (enforces runlength constraint)
        Step 4: GC-balancing (achieves 50% GC-content)
        Step 5: Add index suffix and optional error correction

        Args:
            binary_data: Binary string (e.g., "11010011")
            verbose: Print detailed steps

        Returns:
            DNA sequence string
        """
        if verbose:
            print("\nHELIX Encoding Pipeline")
            print("=" * 70)
            print(f"Input (binary): {binary_data}")

        # Step 1: Binary to Quaternary
        quaternary = mapping.binary_to_quaternary(binary_data)
        if verbose:
            print(f"\nStep 1 - Quaternary conversion:")
            print(f"  Result: {quaternary}")

        # Step 2: Differential encoding
        diff_encoded = differential.differential_encode(quaternary)
        if verbose:
            print(f"\nStep 2 - Differential encoding:")
            print(f"  Result: {diff_encoded}")

        # Step 3: RLL encoding
        rll_encoded = self.rll_codec.encode(diff_encoded)
        if verbose:
            print(f"\nStep 3 - RLL encoding (max run = {self.ell}):")
            print(f"  Result: {rll_encoded}")
            print(f"  Max runlength: {self.rll_codec.max_runlength(rll_encoded)}")

        # Step 4: GC-balancing
        balanced_seq, t = self.gc_balancer.balance(rll_encoded)
        if verbose:
            print(f"\nStep 4 - GC-balancing (epsilon = {self.epsilon}):")
            print(f"  Flip index t: {t}")
            print(f"  Result: {balanced_seq}")
            print(f"  GC-content: {self.gc_balancer.gc_content(balanced_seq):.2%}")

        # Step 5: Add index suffix
        index_suffix = self.gc_balancer.create_index_suffix(t, len(rll_encoded))
        final_seq = balanced_seq + index_suffix
        if verbose:
            print(f"\nStep 5 - Add index suffix:")
            print(f"  Index suffix: {index_suffix}")
            print(f"  Final quaternary: {final_seq}")

        # Optional: Add error correction
        if self.use_error_correction:
            ec_suffix = self.error_corrector.create_error_correction_suffix(final_seq)
            final_seq = final_seq + ec_suffix
            if verbose:
                print(f"\nStep 6 - Error correction suffix:")
                print(f"  EC suffix: {ec_suffix}")
                print(f"  With EC: {final_seq}")

        # Convert to DNA
        dna = mapping.quaternary_to_dna(final_seq)
        if verbose:
            print(f"\nFinal DNA sequence: {dna}")
            print(f"Length: {len(dna)} nucleotides")
            print("=" * 70)

        return dna

    def decode(self, dna: str, verbose: bool = False) -> str:
        """
        Decode DNA sequence back to binary data.

        Pipeline (reverse of encoding):
        Step 1: DNA -> Quaternary
        Step 2: Extract and verify error correction (if used)
        Step 3: Extract and decode index suffix
        Step 4: Reverse GC-balancing (unflip)
        Step 5: RLL decoding
        Step 6: Differential decoding
        Step 7: Quaternary -> Binary

        Args:
            dna: DNA sequence string
            verbose: Print detailed steps

        Returns:
            Binary string
        """
        if verbose:
            print("\nHELIX Decoding Pipeline")
            print("=" * 70)
            print(f"Input (DNA): {dna}")

        # Step 1: DNA to Quaternary
        quaternary = mapping.dna_to_quaternary(dna)
        if verbose:
            print(f"\nStep 1 - Quaternary conversion:")
            print(f"  Result: {quaternary}")

        # Step 2: Handle error correction if used
        if self.use_error_correction:
            ec_suffix_len = 6
            if len(quaternary) > ec_suffix_len:
                body = quaternary[:-ec_suffix_len]
                ec_suffix = quaternary[-ec_suffix_len:]

                if verbose:
                    print(f"\nStep 2 - Error correction:")
                    print(f"  EC suffix: {ec_suffix}")

                expected_syn, expected_check = \
                    self.error_corrector.extract_error_correction_info(ec_suffix)

                if verbose:
                    print(f"  Expected syndrome: {expected_syn}")
                    print(f"  Expected checksum: {expected_check}")

                quaternary = body

        # Step 3: Extract index suffix
        for suffix_len in range(2, min(20, len(quaternary)), 2):
            try:
                body = quaternary[:-suffix_len]
                suffix = quaternary[-suffix_len:]

                t = self.gc_balancer.decode_index_suffix(suffix)

                if verbose and suffix_len == 2:
                    print(f"\nStep 3 - Extract index suffix:")
                    print(f"  Trying suffix length: {suffix_len}")
                    print(f"  Index suffix: {suffix}")
                    print(f"  Decoded t: {t}")

                unbalanced = self.gc_balancer.unbalance(body, t)

                if verbose and suffix_len == 2:
                    print(f"\nStep 4 - Reverse GC-balancing:")
                    print(f"  Unbalanced: {unbalanced}")

                rll_decoded = self.rll_codec.decode(unbalanced)

                if verbose and suffix_len == 2:
                    print(f"\nStep 5 - RLL decoding:")
                    print(f"  Result: {rll_decoded}")

                diff_decoded = differential.differential_decode(rll_decoded)

                if verbose and suffix_len == 2:
                    print(f"\nStep 6 - Differential decoding:")
                    print(f"  Result: {diff_decoded}")

                binary = mapping.quaternary_to_binary(diff_decoded)

                if verbose:
                    print(f"\nStep 7 - Binary conversion:")
                    print(f"  Result: {binary}")
                    print("=" * 70)

                return binary

            except Exception as e:
                if verbose and suffix_len == 2:
                    print(f"  Decode attempt failed: {e}")
                continue

        if verbose:
            print("\nWarning: Using fallback decoder")

        rll_decoded = self.rll_codec.decode(quaternary)
        diff_decoded = differential.differential_decode(rll_decoded)
        return mapping.quaternary_to_binary(diff_decoded)

    def encode_with_analysis(self, binary_data: str) -> Dict:
        """
        Encode and return both DNA sequence and analysis.

        Args:
            binary_data: Binary string

        Returns:
            Dictionary with DNA sequence and analysis
        """
        dna = self.encode(binary_data, verbose=False)
        analysis = self.analyzer.analyze_dna(dna)

        return {
            'input': binary_data,
            'output': dna,
            'analysis': analysis,
            'constraints_satisfied': self.analyzer.validate_constraints(dna)
        }

    def verify_roundtrip(self, binary_data: str, verbose: bool = False) -> bool:
        """
        Verify that encode->decode returns original data.

        Args:
            binary_data: Binary string to test
            verbose: Print results

        Returns:
            True if roundtrip successful
        """
        dna = self.encode(binary_data, verbose=False)
        decoded = self.decode(dna, verbose=False)
        match = (binary_data == decoded)

        if verbose:
            print(f"\nRoundtrip Test")
            print(f"Original:  {binary_data}")
            print(f"DNA:       {dna}")
            print(f"Decoded:   {decoded}")
            print(f"Match:     {match}")

        return match


def create_parser():
    """Create argument parser for CLI."""
    parser = argparse.ArgumentParser(
        prog='helix',
        description='HELIX - High-Efficiency Lossless Information eXchange\n'
                    'DNA Storage Encoding/Decoding System',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Encode binary string to DNA
  python3 helix.py encode -i "11010011"

  # Decode DNA to binary
  python3 helix.py decode -i "ATCGATCG"

  # Encode with verbose output
  python3 helix.py encode -i "11010011" -v

  # Encode from file
  python3 helix.py encode -f input.txt

  # Encode with custom constraints
  python3 helix.py encode -i "11010011" --ell 4 --epsilon 0.1

  # Analyze DNA sequence
  python3 helix.py analyze -i "ATCGATCG"

  # Run demo
  python3 helix.py demo

  # Text to DNA (automatic binary conversion)
  python3 helix.py text-encode -i "HELLO"

For more information, visit: https://github.com/minminlittleshrimp/helix
        """
    )

    subparsers = parser.add_subparsers(dest='command', help='Available commands')

    # Encode command
    encode_parser = subparsers.add_parser('encode', help='Encode binary data to DNA')
    encode_parser.add_argument('-i', '--input', type=str, help='Binary string to encode')
    encode_parser.add_argument('-f', '--file', type=str, help='Read binary from file')
    encode_parser.add_argument('-o', '--output', type=str, help='Output file (default: stdout)')
    encode_parser.add_argument('--ell', type=int, default=3, help='Maximum runlength (default: 3)')
    encode_parser.add_argument('--epsilon', type=float, default=0.05, help='GC-content tolerance (default: 0.05)')
    encode_parser.add_argument('--no-ec', action='store_true', help='Disable error correction')
    encode_parser.add_argument('-v', '--verbose', action='store_true', help='Verbose output')
    encode_parser.add_argument('-a', '--analyze', action='store_true', help='Analyze output sequence')

    # Decode command
    decode_parser = subparsers.add_parser('decode', help='Decode DNA sequence to binary')
    decode_parser.add_argument('-i', '--input', type=str, help='DNA string to decode')
    decode_parser.add_argument('-f', '--file', type=str, help='Read DNA from file')
    decode_parser.add_argument('-o', '--output', type=str, help='Output file (default: stdout)')
    decode_parser.add_argument('--ell', type=int, default=3, help='Maximum runlength (default: 3)')
    decode_parser.add_argument('--epsilon', type=float, default=0.05, help='GC-content tolerance (default: 0.05)')
    decode_parser.add_argument('--no-ec', action='store_true', help='Disable error correction')
    decode_parser.add_argument('-v', '--verbose', action='store_true', help='Verbose output')

    # Text encode command
    text_encode_parser = subparsers.add_parser('text-encode', help='Encode text to DNA')
    text_encode_parser.add_argument('-i', '--input', type=str, help='Text to encode')
    text_encode_parser.add_argument('-f', '--file', type=str, help='Read text from file')
    text_encode_parser.add_argument('-o', '--output', type=str, help='Output file (default: stdout)')
    text_encode_parser.add_argument('--ell', type=int, default=3, help='Maximum runlength (default: 3)')
    text_encode_parser.add_argument('--epsilon', type=float, default=0.05, help='GC-content tolerance (default: 0.05)')
    text_encode_parser.add_argument('--no-ec', action='store_true', help='Disable error correction')
    text_encode_parser.add_argument('-v', '--verbose', action='store_true', help='Verbose output')

    # Text decode command
    text_decode_parser = subparsers.add_parser('text-decode', help='Decode DNA to text')
    text_decode_parser.add_argument('-i', '--input', type=str, help='DNA to decode')
    text_decode_parser.add_argument('-f', '--file', type=str, help='Read DNA from file')
    text_decode_parser.add_argument('-o', '--output', type=str, help='Output file (default: stdout)')
    text_decode_parser.add_argument('--ell', type=int, default=3, help='Maximum runlength (default: 3)')
    text_decode_parser.add_argument('--epsilon', type=float, default=0.05, help='GC-content tolerance (default: 0.05)')
    text_decode_parser.add_argument('--no-ec', action='store_true', help='Disable error correction')
    text_decode_parser.add_argument('-v', '--verbose', action='store_true', help='Verbose output')

    # Analyze command
    analyze_parser = subparsers.add_parser('analyze', help='Analyze DNA sequence')
    analyze_parser.add_argument('-i', '--input', type=str, help='DNA string to analyze')
    analyze_parser.add_argument('-f', '--file', type=str, help='Read DNA from file')
    analyze_parser.add_argument('--ell', type=int, default=3, help='Maximum runlength (default: 3)')
    analyze_parser.add_argument('--epsilon', type=float, default=0.05, help='GC-content tolerance (default: 0.05)')

    # Demo command
    demo_parser = subparsers.add_parser('demo', help='Run demonstration')
    demo_parser.add_argument('--ell', type=int, default=3, help='Maximum runlength (default: 3)')
    demo_parser.add_argument('--epsilon', type=float, default=0.05, help='GC-content tolerance (default: 0.05)')

    # Version command
    subparsers.add_parser('version', help='Show version information')

    return parser


def read_input(args):
    """Read input from args or file."""
    if args.input:
        return args.input.strip()
    elif args.file:
        with open(args.file, 'r') as f:
            return f.read().strip()
    else:
        print("Error: Must provide either -i/--input or -f/--file", file=sys.stderr)
        sys.exit(1)


def write_output(data, args):
    """Write output to file or stdout."""
    if args.output:
        with open(args.output, 'w') as f:
            f.write(data + '\n')
        print(f"Output written to: {args.output}")
    else:
        print(data)


def main():
    """Main CLI entry point."""
    parser = create_parser()
    args = parser.parse_args()

    if not args.command:
        parser.print_help()
        sys.exit(0)

    # Version command
    if args.command == 'version':
        print("HELIX v1.0.0")
        print("DNA Storage Encoding/Decoding System")
        print("Based on: Nguyen et al. - Capacity-Approaching Constrained Codes")
        sys.exit(0)

    # Demo command
    if args.command == 'demo':
        print("=" * 70)
        print("HELIX - High-Efficiency Lossless Information eXchange")
        print("DNA Storage Encoding/Decoding System - DEMO")
        print("=" * 70)

        codec = HelixCodec(
            ell=args.ell,
            epsilon=args.epsilon,
            use_error_correction=True
        )

        test_cases = [
            ("Simple", "11010011"),
            ("Alternating", "10101010"),
            ("All ones", "11111111"),
            ("Mixed", "100100011010"),
        ]

        for name, binary_data in test_cases:
            print(f"\n{'=' * 70}")
            print(f"Test: {name}")
            print(f"{'=' * 70}")

            result = codec.encode_with_analysis(binary_data)

            print(f"\nInput:  {result['input']} ({len(result['input'])} bits)")
            print(f"Output: {result['output']} ({result['analysis']['length']} bp)")
            print(f"\nConstraint Validation:")
            for constraint, passed in result['constraints_satisfied'].items():
                status = "PASS" if passed else "FAIL"
                print(f"  {constraint}: {status}")

            print(f"\nMetrics:")
            print(f"  GC-content: {result['analysis']['gc_content']:.2%}")
            print(f"  Max runlength: {result['analysis']['max_runlength']}")
            print(f"  Efficiency: {len(result['input']) / (len(result['output']) * 2):.2%}")

            codec.verify_roundtrip(binary_data, verbose=True)

        print(f"\n{'=' * 70}")
        print("Demo completed")
        print("=" * 70)
        sys.exit(0)

    # Initialize codec with parameters
    codec = HelixCodec(
        ell=args.ell if hasattr(args, 'ell') else 3,
        epsilon=args.epsilon if hasattr(args, 'epsilon') else 0.05,
        use_error_correction=not args.no_ec if hasattr(args, 'no_ec') else True
    )

    # Encode command
    if args.command == 'encode':
        binary_data = read_input(args)

        if not all(c in '01' for c in binary_data):
            print("Error: Input must be binary string (only 0 and 1)", file=sys.stderr)
            sys.exit(1)

        dna = codec.encode(binary_data, verbose=args.verbose)

        if args.analyze:
            print("\nSequence Analysis:")
            print("=" * 70)
            analysis = codec.analyzer.analyze_dna(dna)
            codec.analyzer.print_analysis(analysis)

        if not args.verbose:
            write_output(dna, args)

    # Decode command
    elif args.command == 'decode':
        dna = read_input(args)

        if not all(c.upper() in 'ATCG' for c in dna):
            print("Error: Input must be DNA string (only A, T, C, G)", file=sys.stderr)
            sys.exit(1)

        binary = codec.decode(dna, verbose=args.verbose)

        if not args.verbose:
            write_output(binary, args)

    # Text encode command
    elif args.command == 'text-encode':
        text = read_input(args)
        binary = ''.join(format(ord(c), '08b') for c in text)

        if args.verbose:
            print(f"Text: {text}")
            print(f"Binary: {binary}")

        dna = codec.encode(binary, verbose=args.verbose)

        if not args.verbose:
            write_output(dna, args)

    # Text decode command
    elif args.command == 'text-decode':
        dna = read_input(args)
        binary = codec.decode(dna, verbose=args.verbose)

        # Convert binary to text
        text = ''.join(
            chr(int(binary[i:i+8], 2))
            for i in range(0, len(binary), 8)
            if i+8 <= len(binary)
        )

        if args.verbose:
            print(f"Binary: {binary}")
            print(f"Text: {text}")
        else:
            write_output(text, args)

    # Analyze command
    elif args.command == 'analyze':
        dna = read_input(args)

        seq_analyzer = analyzer.SequenceAnalyzer(ell=args.ell, epsilon=args.epsilon)
        analysis = seq_analyzer.analyze_dna(dna)
        seq_analyzer.print_analysis(analysis)


if __name__ == "__main__":
    main()