"""
analyzer.py - HELIX DNA Storage System
Sequence analysis and constraint validation tools
"""

from typing import List, Dict
import mapping


class SequenceAnalyzer:
    """
    Analyzes DNA sequences for constraint satisfaction and quality metrics.
    """

    def __init__(self, ell: int = 3, epsilon: float = 0.05):
        """
        Initialize sequence analyzer.

        Args:
            ell: Maximum allowed runlength
            epsilon: GC-content tolerance
        """
        self.ell = ell
        self.epsilon = epsilon

    def analyze_dna(self, dna: str) -> Dict:
        """
        Complete analysis of DNA sequence.

        Args:
            dna: DNA string

        Returns:
            Dictionary containing all analysis metrics
        """
        quaternary = mapping.dna_to_quaternary(dna)

        return {
            'sequence': dna,
            'length': len(dna),
            'quaternary': quaternary,
            'gc_content': self.compute_gc_content(dna),
            'gc_balanced': self.is_gc_balanced(dna),
            'gc_target_range': (0.5 - self.epsilon, 0.5 + self.epsilon),
            'max_runlength': self.compute_max_runlength(dna),
            'runlength_ok': self.check_runlength_constraint(dna),
            'runlength_limit': self.ell,
            'nucleotide_counts': self.count_nucleotides(dna),
            'homopolymer_runs': self.find_homopolymer_runs(dna),
        }

    def analyze_quaternary(self, quaternary: List[int]) -> Dict:
        """
        Analyze quaternary sequence.

        Args:
            quaternary: Quaternary sequence

        Returns:
            Dictionary containing analysis metrics
        """
        dna = mapping.quaternary_to_dna(quaternary)
        return self.analyze_dna(dna)

    def compute_gc_content(self, dna: str) -> float:
        """
        Calculate GC-content ratio.

        Args:
            dna: DNA string

        Returns:
            GC-content as ratio [0.0, 1.0]
        """
        if not dna:
            return 0.0
        gc_count = sum(1 for nucleotide in dna.upper() if nucleotide in ['C', 'G'])
        return gc_count / len(dna)

    def is_gc_balanced(self, dna: str) -> bool:
        """
        Check if DNA sequence meets GC-balance constraint.

        Args:
            dna: DNA string

        Returns:
            True if |GC_content - 0.5| <= epsilon
        """
        gc_content = self.compute_gc_content(dna)
        return abs(gc_content - 0.5) <= self.epsilon

    def compute_max_runlength(self, dna: str) -> int:
        """
        Find maximum homopolymer runlength in DNA sequence.

        Args:
            dna: DNA string

        Returns:
            Length of longest homopolymer run
        """
        if not dna:
            return 0

        max_run = 1
        current_run = 1

        for i in range(1, len(dna)):
            if dna[i].upper() == dna[i-1].upper():
                current_run += 1
                max_run = max(max_run, current_run)
            else:
                current_run = 1

        return max_run

    def check_runlength_constraint(self, dna: str) -> bool:
        """
        Check if DNA sequence meets runlength constraint.

        Args:
            dna: DNA string

        Returns:
            True if max_runlength <= ell
        """
        return self.compute_max_runlength(dna) <= self.ell

    def count_nucleotides(self, dna: str) -> Dict[str, int]:
        """
        Count occurrences of each nucleotide.

        Args:
            dna: DNA string

        Returns:
            Dictionary mapping nucleotide to count
        """
        counts = {'A': 0, 'T': 0, 'C': 0, 'G': 0}
        for nucleotide in dna.upper():
            if nucleotide in counts:
                counts[nucleotide] += 1
        return counts

    def find_homopolymer_runs(self, dna: str) -> List[Dict]:
        """
        Find all homopolymer runs in sequence.

        Args:
            dna: DNA string

        Returns:
            List of dictionaries describing each run:
            {'nucleotide': str, 'start': int, 'length': int}
        """
        if not dna:
            return []

        runs = []
        current_nucleotide = dna[0].upper()
        current_start = 0
        current_length = 1

        for i in range(1, len(dna)):
            nucleotide = dna[i].upper()
            if nucleotide == current_nucleotide:
                current_length += 1
            else:
                if current_length > 1:
                    runs.append({
                        'nucleotide': current_nucleotide,
                        'start': current_start,
                        'length': current_length
                    })
                current_nucleotide = nucleotide
                current_start = i
                current_length = 1

        # Don't forget the last run
        if current_length > 1:
            runs.append({
                'nucleotide': current_nucleotide,
                'start': current_start,
                'length': current_length
            })

        return runs

    def print_analysis(self, analysis: Dict) -> None:
        """
        Pretty-print analysis results.

        Args:
            analysis: Analysis dictionary from analyze_dna()
        """
        print("\nSequence Analysis")
        print("=" * 70)
        print(f"DNA Sequence:     {analysis['sequence']}")
        print(f"Length:           {analysis['length']} bp")
        print(f"Quaternary:       {analysis['quaternary']}")
        print()
        print(f"GC-Content:       {analysis['gc_content']:.2%}")
        print(f"GC-Balanced:      {analysis['gc_balanced']}")
        print(f"Target Range:     {analysis['gc_target_range'][0]:.2%} - "
              f"{analysis['gc_target_range'][1]:.2%}")
        print()
        print(f"Max Runlength:    {analysis['max_runlength']}")
        print(f"Runlength OK:     {analysis['runlength_ok']}")
        print(f"Runlength Limit:  {analysis['runlength_limit']}")
        print()
        print("Nucleotide Counts:")
        for nucleotide in ['A', 'T', 'C', 'G']:
            count = analysis['nucleotide_counts'][nucleotide]
            pct = count / analysis['length'] * 100 if analysis['length'] > 0 else 0
            print(f"  {nucleotide}: {count:3d} ({pct:5.1f}%)")

        if analysis['homopolymer_runs']:
            print()
            print("Homopolymer Runs:")
            for run in analysis['homopolymer_runs']:
                print(f"  {run['nucleotide']} x {run['length']} at position {run['start']}")

        print("=" * 70)

    def compare_sequences(self, dna1: str, dna2: str) -> Dict:
        """
        Compare two DNA sequences.

        Args:
            dna1: First DNA string
            dna2: Second DNA string

        Returns:
            Dictionary with comparison metrics
        """
        analysis1 = self.analyze_dna(dna1)
        analysis2 = self.analyze_dna(dna2)

        return {
            'length_diff': analysis2['length'] - analysis1['length'],
            'gc_content_diff': analysis2['gc_content'] - analysis1['gc_content'],
            'runlength_diff': analysis2['max_runlength'] - analysis1['max_runlength'],
            'both_gc_balanced': analysis1['gc_balanced'] and analysis2['gc_balanced'],
            'both_runlength_ok': analysis1['runlength_ok'] and analysis2['runlength_ok'],
        }

    def validate_constraints(self, dna: str) -> Dict[str, bool]:
        """
        Check all constraints for a DNA sequence.

        Args:
            dna: DNA string

        Returns:
            Dictionary mapping constraint name to pass/fail
        """
        return {
            'gc_balanced': self.is_gc_balanced(dna),
            'runlength_ok': self.check_runlength_constraint(dna),
            'valid_nucleotides': all(n.upper() in ['A', 'T', 'C', 'G'] for n in dna),
        }


if __name__ == "__main__":
    # Test sequence analyzer
    print("Testing analyzer.py")
    print("-" * 70)

    analyzer = SequenceAnalyzer(ell=3, epsilon=0.05)

    test_sequences = [
        "ATCGATCG",           # Balanced, no long runs
        "AAAATTTCCCGGG",      # Long runs, balanced
        "ATATATATATAT",       # No long runs, AT-rich
        "CGCGCGCGCGCG",       # No long runs, GC-rich
    ]

    for dna in test_sequences:
        analysis = analyzer.analyze_dna(dna)
        analyzer.print_analysis(analysis)
        print()

        constraints = analyzer.validate_constraints(dna)
        print("Constraint Validation:")
        for constraint, passed in constraints.items():
            status = "PASS" if passed else "FAIL"
            print(f"  {constraint}: {status}")
        print()

    # Test comparison
    print("\n" + "-" * 70)
    print("Sequence Comparison")
    print("-" * 70)

    dna1 = "ATCGATCG"
    dna2 = "ATCGGGGATCG"

    print(f"Sequence 1: {dna1}")
    print(f"Sequence 2: {dna2}")

    comparison = analyzer.compare_sequences(dna1, dna2)
    print(f"\nComparison:")
    print(f"  Length difference:     {comparison['length_diff']}")
    print(f"  GC-content difference: {comparison['gc_content_diff']:.2%}")
    print(f"  Runlength difference:  {comparison['runlength_diff']}")
    print(f"  Both GC-balanced:      {comparison['both_gc_balanced']}")
    print(f"  Both runlength OK:     {comparison['both_runlength_ok']}")

    print("\n" + "-" * 70)
    print("All tests completed")