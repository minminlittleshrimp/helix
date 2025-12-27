"""
error_correction.py - HELIX DNA Storage System
Varshamov-Tenengolts (VT) syndrome for single-edit error correction
"""

from typing import List, Optional, Tuple


class VTErrorCorrection:
    """
    Varshamov-Tenengolts error correction for single insertion, deletion, or substitution.
    """

    def __init__(self):
        """Initialize VT error correction system."""
        pass

    def compute_syndrome(self, sequence: List[int]) -> int:
        """
        Compute Varshamov-Tenengolts syndrome.

        Formula: Syn(x) = sum(i * x[i]) mod (2n)
        where i is 1-indexed

        Args:
            sequence: Quaternary sequence

        Returns:
            VT syndrome value
        """
        n = len(sequence)
        if n == 0:
            return 0

        syndrome = sum((i + 1) * sequence[i] for i in range(n))
        return syndrome % (2 * n)

    def compute_checksum(self, sequence: List[int]) -> int:
        """
        Compute simple checksum: sum of all symbols mod 4.

        Args:
            sequence: Quaternary sequence

        Returns:
            Checksum value [0-3]
        """
        return sum(sequence) % 4

    def create_error_correction_suffix(self, sequence: List[int]) -> List[int]:
        """
        Create error correction suffix containing syndrome and checksum.
        Uses interleaved encoding to maintain balance.

        Args:
            sequence: Quaternary sequence

        Returns:
            Error correction suffix (interleaved)
        """
        syndrome = self.compute_syndrome(sequence)
        checksum = self.compute_checksum(sequence)

        # Encode syndrome and checksum as balanced suffix
        # Convert to quaternary digits
        syndrome_quat = self._int_to_quaternary(syndrome, min_length=2)
        checksum_quat = [checksum]

        # Interleave with flipped versions
        ec_suffix = []
        flip_map = {0: 2, 2: 0, 1: 3, 3: 1}

        for symbol in syndrome_quat + checksum_quat:
            ec_suffix.append(symbol)
            ec_suffix.append(flip_map[symbol])

        return ec_suffix

    def extract_error_correction_info(self, suffix: List[int]) -> Tuple[int, int]:
        """
        Extract syndrome and checksum from error correction suffix.

        Args:
            suffix: Error correction suffix

        Returns:
            Tuple of (syndrome, checksum)
        """
        # Extract every other symbol (the original values)
        original = [suffix[i] for i in range(0, len(suffix), 2)]

        # Last symbol is checksum, rest is syndrome
        if len(original) < 1:
            return 0, 0

        checksum = original[-1]
        syndrome_quat = original[:-1]

        # Convert syndrome from quaternary to int
        syndrome = 0
        for digit in syndrome_quat:
            syndrome = syndrome * 4 + digit

        return syndrome, checksum

    def verify_sequence(self, sequence: List[int], expected_syndrome: int,
                       expected_checksum: int) -> bool:
        """
        Verify if sequence matches expected syndrome and checksum.

        Args:
            sequence: Quaternary sequence to verify
            expected_syndrome: Expected VT syndrome
            expected_checksum: Expected checksum

        Returns:
            True if sequence is valid
        """
        actual_syndrome = self.compute_syndrome(sequence)
        actual_checksum = self.compute_checksum(sequence)

        return (actual_syndrome == expected_syndrome and
                actual_checksum == expected_checksum)

    def detect_error(self, sequence: List[int], expected_syndrome: int,
                    expected_checksum: int) -> Optional[str]:
        """
        Detect type of single-edit error if present.

        Args:
            sequence: Possibly corrupted sequence
            expected_syndrome: Expected VT syndrome
            expected_checksum: Expected checksum

        Returns:
            Error type: None, 'insertion', 'deletion', or 'substitution'
        """
        if self.verify_sequence(sequence, expected_syndrome, expected_checksum):
            return None

        # Simplified error detection based on syndrome difference
        actual_syndrome = self.compute_syndrome(sequence)
        actual_checksum = self.compute_checksum(sequence)

        syndrome_diff = (actual_syndrome - expected_syndrome) % (2 * len(sequence))
        checksum_diff = (actual_checksum - expected_checksum) % 4

        # These are heuristics; full correction requires more sophisticated analysis
        if checksum_diff == 0:
            return 'substitution'  # Same length, different values
        elif syndrome_diff != 0:
            if len(sequence) > 0:
                # Could be insertion or deletion
                return 'insertion_or_deletion'

        return 'unknown'

    def _int_to_quaternary(self, value: int, min_length: int = 1) -> List[int]:
        """
        Convert integer to quaternary representation.

        Args:
            value: Integer to convert
            min_length: Minimum length of output

        Returns:
            List of quaternary digits
        """
        if value == 0:
            result = [0]
        else:
            result = []
            temp = value
            while temp > 0:
                result.append(temp % 4)
                temp //= 4
            result.reverse()

        # Pad to minimum length
        while len(result) < min_length:
            result.insert(0, 0)

        return result


class ExtendedVT:
    """
    Extended VT codes for dual-strand DNA storage.
    Each DNA strand has upper and lower sequences.
    """

    def __init__(self):
        """Initialize extended VT system."""
        self.vt = VTErrorCorrection()

    def compute_dual_syndromes(self, upper: List[int],
                               lower: List[int]) -> Tuple[int, int]:
        """
        Compute VT syndromes for both upper and lower strands.

        Args:
            upper: Upper strand sequence
            lower: Lower strand sequence

        Returns:
            Tuple of (upper_syndrome, lower_syndrome)
        """
        upper_syn = self.vt.compute_syndrome(upper)
        lower_syn = self.vt.compute_syndrome(lower)
        return upper_syn, lower_syn

    def create_dual_error_correction(self, upper: List[int],
                                     lower: List[int]) -> List[int]:
        """
        Create combined error correction suffix for dual strands.

        Args:
            upper: Upper strand sequence
            lower: Lower strand sequence

        Returns:
            Combined error correction suffix
        """
        upper_suffix = self.vt.create_error_correction_suffix(upper)
        lower_suffix = self.vt.create_error_correction_suffix(lower)

        # Combine suffixes
        return upper_suffix + lower_suffix


if __name__ == "__main__":
    # Test VT error correction
    print("Testing error_correction.py")
    print("-" * 50)

    ec = VTErrorCorrection()

    test_sequences = [
        [1, 2, 3, 0, 1, 2],
        [0, 0, 1, 1, 2, 2, 3, 3],
        [3, 2, 1, 0],
    ]

    for i, seq in enumerate(test_sequences, 1):
        print(f"\nTest Case {i}")
        print(f"Sequence:  {seq}")

        syndrome = ec.compute_syndrome(seq)
        print(f"Syndrome:  {syndrome}")

        checksum = ec.compute_checksum(seq)
        print(f"Checksum:  {checksum}")

        ec_suffix = ec.create_error_correction_suffix(seq)
        print(f"EC Suffix: {ec_suffix}")

        # Decode suffix
        syn_decoded, check_decoded = ec.extract_error_correction_info(ec_suffix)
        print(f"Decoded syndrome: {syn_decoded}, checksum: {check_decoded}")
        print(f"Match: {syndrome == syn_decoded and checksum == check_decoded}")

        # Verify
        is_valid = ec.verify_sequence(seq, syndrome, checksum)
        print(f"Valid: {is_valid}")

        # Test error detection with corrupted sequence
        corrupted = seq.copy()
        if len(corrupted) > 0:
            corrupted[0] = (corrupted[0] + 1) % 4  # Introduce error
            error_type = ec.detect_error(corrupted, syndrome, checksum)
            print(f"Corrupted sequence: {corrupted}")
            print(f"Error detected: {error_type}")

    print("\n" + "-" * 50)

    # Test dual-strand
    print("\nTesting dual-strand VT")
    ext_vt = ExtendedVT()

    upper = [1, 2, 3, 0]
    lower = [0, 3, 2, 1]

    print(f"Upper: {upper}")
    print(f"Lower: {lower}")

    upper_syn, lower_syn = ext_vt.compute_dual_syndromes(upper, lower)
    print(f"Upper syndrome: {upper_syn}, Lower syndrome: {lower_syn}")

    dual_ec = ext_vt.create_dual_error_correction(upper, lower)
    print(f"Dual EC suffix: {dual_ec}")

    print("\n" + "-" * 50)
    print("All tests completed")