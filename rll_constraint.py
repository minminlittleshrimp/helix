"""
rll_constraint.py - HELIX DNA Storage System
Method B: Runlength constraint encoder/decoder
Prevents homopolymer runs longer than ell
"""

from typing import List


class RLLCodec:
    """
    Runlength-Limited codec implementing Method B from the paper.
    Prevents sequences from having runs of identical symbols longer than ell.
    """

    def __init__(self, ell: int = 3):
        """
        Initialize RLL codec.

        Args:
            ell: Maximum allowed runlength (homopolymer limit)
        """
        self.ell = ell

    def encode(self, data: List[int]) -> List[int]:
        """
        Encode quaternary sequence to satisfy runlength constraint.

        Algorithm:
        1. Append termination symbol '0'
        2. Scan for forbidden substrings (ell consecutive zeros)
        3. Replace each forbidden substring with pointer 'Re' where e != 0
        4. Repeat until no forbidden substrings remain

        Args:
            data: Input quaternary sequence

        Returns:
            RLL-encoded sequence (no runs of 0s longer than ell)
        """
        # Append termination symbol
        x = data + [0]

        # Iteratively remove forbidden substrings
        max_iterations = len(x) * 2  # Safety limit
        iteration = 0

        while iteration < max_iterations:
            forbidden_pos = self._find_forbidden_substring(x)

            if forbidden_pos is None:
                # No more forbidden substrings
                break

            # Replace forbidden substring with pointer
            # Use pointer pattern: [1, 1] (R=1, e=1)
            x = x[:forbidden_pos] + [1, 1] + x[forbidden_pos + self.ell:]
            iteration += 1

        return x

    def decode(self, encoded: List[int]) -> List[int]:
        """
        Decode RLL-encoded sequence back to original.

        Algorithm:
        1. Scan from right to left
        2. Identify pointer patterns [1, 1]
        3. Replace each pointer with forbidden substring (ell zeros)
        4. Remove termination symbol

        Args:
            encoded: RLL-encoded sequence

        Returns:
            Original quaternary sequence
        """
        x = encoded.copy()

        # Scan from right to left
        i = len(x) - 1
        while i >= 1:
            # Check for pointer pattern [1, 1]
            if x[i] == 1 and x[i-1] == 1:
                # Replace pointer with forbidden substring
                x = x[:i-1] + [0] * self.ell + x[i+1:]
                i -= 1
            else:
                i -= 1

        # Remove termination symbol (last 0)
        if x and x[-1] == 0:
            x = x[:-1]

        return x

    def _find_forbidden_substring(self, sequence: List[int]) -> int:
        """
        Find the first occurrence of a forbidden substring (ell consecutive zeros).

        Args:
            sequence: Quaternary sequence to search

        Returns:
            Index of first forbidden substring, or None if not found
        """
        for i in range(len(sequence) - self.ell + 1):
            if all(sequence[i+j] == 0 for j in range(self.ell)):
                return i
        return None

    def has_forbidden_substring(self, sequence: List[int]) -> bool:
        """
        Check if sequence contains any forbidden substrings.

        Args:
            sequence: Quaternary sequence

        Returns:
            True if forbidden substring exists
        """
        return self._find_forbidden_substring(sequence) is not None

    def max_runlength(self, sequence: List[int]) -> int:
        """
        Calculate maximum runlength in sequence.

        Args:
            sequence: Quaternary sequence

        Returns:
            Length of longest run of identical symbols
        """
        if not sequence:
            return 0

        max_run = 1
        current_run = 1

        for i in range(1, len(sequence)):
            if sequence[i] == sequence[i-1]:
                current_run += 1
                max_run = max(max_run, current_run)
            else:
                current_run = 1

        return max_run


if __name__ == "__main__":
    # Test RLL codec
    print("Testing rll_constraint.py")
    print("-" * 50)

    codec = RLLCodec(ell=3)

    test_cases = [
        [0, 0, 0, 1, 2],           # Contains forbidden substring
        [1, 0, 0, 0, 0, 2],        # Multiple zeros
        [0, 1, 0, 1, 0],           # No forbidden substring
        [0, 0, 0, 0, 0, 0],        # All zeros
    ]

    for i, original in enumerate(test_cases, 1):
        print(f"\nTest Case {i}")
        print(f"Original:     {original}")
        print(f"Has forbidden: {codec.has_forbidden_substring(original)}")
        print(f"Max runlength: {codec.max_runlength(original)}")

        encoded = codec.encode(original)
        print(f"Encoded:      {encoded}")
        print(f"Max runlength (encoded): {codec.max_runlength(encoded)}")

        decoded = codec.decode(encoded)
        print(f"Decoded:      {decoded}")
        print(f"Match: {original == decoded}")

    print("\n" + "-" * 50)
    print("All tests completed")