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
        1. Escape any existing [3,2] patterns to [3,1,2] to avoid conflicts
        2. Append termination symbol '0'
        3. Scan for forbidden substrings (ell consecutive zeros)
        4. Replace each forbidden substring with pointer [3,2]
        5. Track number of pointers for proper decoding
        6. Repeat until no forbidden substrings remain

        Args:
            data: Input quaternary sequence

        Returns:
            RLL-encoded sequence (no runs of 0s longer than ell)
        """
        # Step 1: Escape existing [3,2] patterns
        x = self._escape_pointer_pattern(data)

        # Step 2: Append termination symbol
        x = x + [0]

        # Step 3-6: Iteratively remove forbidden substrings
        pointer_count = 0
        max_iterations = len(x) * 2  # Safety limit
        iteration = 0

        while iteration < max_iterations:
            forbidden_pos = self._find_forbidden_substring(x)

            if forbidden_pos is None:
                # No more forbidden substrings
                break

            # Replace forbidden substring with pointer pattern [3, 2]
            x = x[:forbidden_pos] + [3, 2] + x[forbidden_pos + self.ell:]
            pointer_count += 1
            iteration += 1

        # Always encode pointer count into the sequence (even if 0)
        # Use marker [2, 2] followed by count as quaternary digits (LSB first)
        # Encode count as base-4 with exactly 4 digits (supports 0-255 pointers)
        # Note: Method B from paper supports arbitrary lengths via linear O(n) processing
        count_quat = [pointer_count % 4, (pointer_count // 4) % 4,
                      (pointer_count // 16) % 4, (pointer_count // 64) % 4]

        # Add glue symbol before [2,2] marker to prevent runlength violations
        # Junction Rule (Corollary 24): Insert symbol that differs from both neighbors
        last_symbol = x[-1] if x else 0
        first_marker = 2
        if last_symbol == first_marker:
            # Need glue: last symbol same as first marker symbol
            glue_options = [s for s in [0, 1, 2, 3] if s != last_symbol and s != first_marker]
            if not glue_options:
                glue_options = [0, 1, 3]  # Fallback (anything != 2)
            glue1 = glue_options[0]
            x = x + [glue1]

        # Add marker
        x = x + [2, 2]

        # Add glue between marker and count if needed
        # (if first count digit is also 2)
        first_count = count_quat[0]
        if first_count == 2:
            glue2 = 0  # Any symbol != 2
            x = x + [glue2]

        x = x + count_quat
        return x

    def decode(self, encoded: List[int]) -> List[int]:
        """
        Decode RLL-encoded sequence back to original.

        Algorithm:
        1. Check for pointer count marker [2, 2, count]
        2. Scan from LEFT to RIGHT exactly 'count' times
        3. Replace each pointer [3,2] with forbidden substring (ell zeros)
        4. Remove termination symbol
        5. Un-escape [3,1,2] patterns back to [3,2]

        Args:
            encoded: RLL-encoded sequence

        Returns:
            Original quaternary sequence

        Raises:
            ValueError: If RLL marker [2, 2, *, *] not found at end
        """
        x = encoded.copy()

        # Decode RLL marker and count from end
        # Format: [data] + [opt glue1] + [2, 2] + [opt glue2] + [d0, d1, d2, d3]
        #
        # Rules:
        # - glue1 exists if last data symbol == 2
        # - glue2 exists if d0 == 2
        # - Count is encoded as 4 quaternary digits (LSB first, supports 0-255)

        if len(x) < 6:  # Minimum: marker(2) + count(4)
            raise ValueError(f"Sequence too short to contain RLL marker")

        # Step 1: Extract count digits from end (4 digits now)
        count_d0 = x[-4]
        count_d1 = x[-3]
        count_d2 = x[-2]
        count_d3 = x[-1]
        pointer_count = count_d0 + count_d1 * 4 + count_d2 * 16 + count_d3 * 64

        # Step 2: Find marker position based on whether glue2 exists
        if count_d0 == 2:
            # glue2 exists at position -5
            # marker should be at positions -7 and -6
            if len(x) < 7 or x[-7] != 2 or x[-6] != 2:
                raise ValueError(f"RLL marker [2, 2] not found at expected position")
            marker_end = len(x) - 7  # Position before marker
        else:
            # No glue2, marker should be at positions -6 and -5
            if x[-6] != 2 or x[-5] != 2:
                raise ValueError(f"RLL marker [2, 2] not found at expected position")
            marker_end = len(x) - 6  # Position before marker

        # Step 3: Remove everything from marker onwards
        x = x[:marker_end]

        # Step 3b: Check for glue1
        # Glue1 is added if body ends with 2 after all pointer replacements
        # To detect this, we simulate replacing the last pointer and check if
        # the resulting body would end with 2

        if len(x) >= 3 and x[-3] == 3 and x[-2] == 2 and x[-1] in {0, 1, 3}:
            # Find all [3,2] patterns
            patterns = [i for i in range(len(x) - 1) if x[i] == 3 and x[i+1] == 2]
            end_pattern_pos = len(x) - 3

            if end_pattern_pos in patterns and pointer_count > 0:
                pattern_index = patterns.index(end_pattern_pos)
                # Check if this is the LAST pointer (0-indexed: pointer_count - 1)
                if pattern_index == pointer_count - 1:
                    # Simulate replacing this last [3,2] with the forbidden substring
                    # If body would end with 2 after replacement, glue was added
                    simulated = x[:end_pattern_pos] + [0] * self.ell + x[end_pattern_pos + 2:-1]
                    # After replacement, does body end with 2?
                    # If yes, then glue was added to prevent [2,2,2]
                    if simulated and simulated[-1] == 2:
                        x = x[:-1]  # Remove glue

        # Step 4: Replace exactly 'pointer_count' [3,2] patterns from LEFT to RIGHT
        # At this point, ALL [3,2] patterns in the sequence are pointers
        # (Original [3,2] patterns in data were escaped to [3,1,2] by encoder)
        replacements_made = 0
        i = 0
        while i < len(x) - 1 and replacements_made < pointer_count:
            # Check for pointer pattern [3, 2]
            if x[i] == 3 and x[i+1] == 2:
                # Replace pointer with forbidden substring
                x = x[:i] + [0] * self.ell + x[i+2:]
                replacements_made += 1
                # After replacement, skip past the inserted zeros
                i += self.ell
            else:
                i += 1

        # Step 5: Remove termination symbol (last 0)
        # The encoder always appends a 0 at step 2, so we must remove it
        if x and x[-1] == 0:
            x = x[:-1]

        # Step 6: Un-escape [3,1,2] patterns back to [3,2]
        # This restores original [3,2] patterns that were in the input data
        x = self._unescape_pointer_pattern(x)

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

    def _escape_pointer_pattern(self, sequence: List[int]) -> List[int]:
        """
        Escape existing [3,2] patterns to [3,1,2] to avoid conflicts with pointers.

        Args:
            sequence: Quaternary sequence

        Returns:
            Sequence with [3,2] patterns escaped
        """
        result = []
        i = 0
        while i < len(sequence):
            if i < len(sequence) - 1 and sequence[i] == 3 and sequence[i+1] == 2:
                # Found [3,2], escape it to [3,1,2]
                result.extend([3, 1, 2])
                i += 2
            else:
                result.append(sequence[i])
                i += 1
        return result

    def _unescape_pointer_pattern(self, sequence: List[int]) -> List[int]:
        """
        Un-escape [3,1,2] patterns back to [3,2].

        Args:
            sequence: Quaternary sequence with escaped patterns

        Returns:
            Sequence with [3,1,2] patterns un-escaped to [3,2]
        """
        result = []
        i = 0
        while i < len(sequence):
            if (i < len(sequence) - 2 and
                sequence[i] == 3 and sequence[i+1] == 1 and sequence[i+2] == 2):
                # Found escaped pattern [3,1,2], restore to [3,2]
                result.extend([3, 2])
                i += 3
            else:
                result.append(sequence[i])
                i += 1
        return result


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