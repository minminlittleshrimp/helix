"""
gc_balance.py - HELIX DNA Storage System
Method D: GC-content balancing through prefix flipping
"""

from typing import List, Tuple


class GCBalancer:
    """
    GC-content balancing using Method D from the paper.
    Ensures GC-content (percentage of C and G) is within 50% ± epsilon.
    """

    def __init__(self, epsilon: float = 0.05):
        """
        Initialize GC balancer.

        Args:
            epsilon: GC-content tolerance (target: 0.5 ± epsilon)
        """
        self.epsilon = epsilon
        # Flipping rule: f(0)=2, f(2)=0, f(1)=3, f(3)=1
        self.flip_map = {0: 2, 2: 0, 1: 3, 3: 1}

    def flip_symbol(self, symbol: int) -> int:
        """
        Apply flipping rule to a single symbol.
        Rule: f(0)=2, f(2)=0, f(1)=3, f(3)=1
        Swaps non-GC (A,T) with GC (C,G)

        Args:
            symbol: Quaternary symbol [0-3]

        Returns:
            Flipped symbol
        """
        return self.flip_map[symbol]

    def flip_sequence(self, sequence: List[int], length: int) -> List[int]:
        """
        Flip the first 'length' symbols of the sequence.

        Args:
            sequence: Quaternary sequence
            length: Number of symbols to flip from the start

        Returns:
            Sequence with first 'length' symbols flipped
        """
        result = sequence.copy()
        for i in range(min(length, len(sequence))):
            result[i] = self.flip_symbol(result[i])
        return result

    def gc_weight(self, sequence: List[int]) -> int:
        """
        Count GC symbols (2=C and 3=G) in sequence.

        Args:
            sequence: Quaternary sequence

        Returns:
            Number of GC symbols
        """
        return sum(1 for s in sequence if s in [2, 3])

    def gc_content(self, sequence: List[int]) -> float:
        """
        Calculate GC-content as a ratio.

        Args:
            sequence: Quaternary sequence

        Returns:
            GC-content ratio [0.0, 1.0]
        """
        if not sequence:
            return 0.0
        return self.gc_weight(sequence) / len(sequence)

    def is_balanced(self, sequence: List[int]) -> bool:
        """
        Check if sequence is epsilon-balanced.
        Balanced means: |GC_content - 0.5| <= epsilon

        Args:
            sequence: Quaternary sequence

        Returns:
            True if sequence is balanced
        """
        if not sequence:
            return True
        gc_ratio = self.gc_content(sequence)
        return abs(gc_ratio - 0.5) <= self.epsilon

    def generate_search_set(self, n: int) -> List[int]:
        """
        Generate search set S_{epsilon,n} = {0, n} ∪ {2⌊εn⌋, 4⌊εn⌋, ...}

        Enhanced to include more candidates for better coverage.

        Args:
            n: Length of sequence

        Returns:
            Sorted list of candidate positions
        """
        S = set([0, n])

        # Add the standard search points
        step = 2 * int(self.epsilon * n)
        if step > 0:
            i = step
            while i < n:
                S.add(i)
                i += step

        # For short sequences or when step is 0, add more granular search points
        if n <= 20 or step == 0:
            # Add every position for short sequences
            for i in range(n + 1):
                S.add(i)
        else:
            # Add quarter points for better coverage
            S.add(n // 4)
            S.add(n // 2)
            S.add(3 * n // 4)

        return sorted(S)

    def balance(self, sequence: List[int]) -> Tuple[List[int], int]:
        """
        Find balancing index t and flip prefix to achieve epsilon-balance.

        Algorithm:
        1. Generate search set of candidate flip positions
        2. For each position t, flip first t symbols
        3. Check if resulting sequence is balanced
        4. Return first balanced sequence found

        Args:
            sequence: Input quaternary sequence

        Returns:
            Tuple of (balanced_sequence, balancing_index)
        """
        n = len(sequence)
        if n == 0:
            return sequence, 0

        search_set = self.generate_search_set(n)

        # Try each candidate position
        for t in search_set:
            test_seq = self.flip_sequence(sequence, t)
            if self.is_balanced(test_seq):
                return test_seq, t

        # If no exact balance found, find the closest one
        best_t = 0
        best_diff = float('inf')
        best_seq = sequence

        for t in search_set:
            test_seq = self.flip_sequence(sequence, t)
            gc = self.gc_content(test_seq)
            diff = abs(gc - 0.5)
            if diff < best_diff:
                best_diff = diff
                best_t = t
                best_seq = test_seq

        return best_seq, best_t

    def create_index_suffix(self, t: int, n: int) -> List[int]:
        """
        Create balanced suffix p by interleaving quaternary representation
        of t with its flipped version.

        This suffix encodes the balancing index so it can be recovered during decoding.

        Args:
            t: Balancing index
            n: Original sequence length (for validation)

        Returns:
            Interleaved index suffix
        """
        # Convert t to quaternary representation
        if t == 0:
            tau = [0]
        else:
            tau = []
            temp = t
            while temp > 0:
                tau.append(temp % 4)
                temp //= 4
            tau.reverse()

        # Interleave tau with f(tau)
        p = []
        for symbol in tau:
            p.append(symbol)
            p.append(self.flip_symbol(symbol))

        return p

    def decode_index_suffix(self, suffix: List[int]) -> int:
        """
        Decode the balancing index from interleaved suffix.
        Extract every other symbol starting from index 0.

        Args:
            suffix: Interleaved index suffix

        Returns:
            Decoded balancing index t

        Raises:
            ValueError: If suffix is not valid interleaved format
        """
        # Validate interleaved format: suffix[i+1] should be flip of suffix[i]
        flip_map = {0: 2, 2: 0, 1: 3, 3: 1}
        if len(suffix) % 2 != 0:
            raise ValueError("Index suffix must have even length")

        for i in range(0, len(suffix), 2):
            if suffix[i+1] != flip_map[suffix[i]]:
                raise ValueError(f"Index suffix not properly interleaved at position {i}")

        # Extract every other symbol (the original tau)
        tau = [suffix[i] for i in range(0, len(suffix), 2)]

        # Convert from quaternary to decimal
        t = 0
        for symbol in tau:
            t = t * 4 + symbol
        return t

    def unbalance(self, sequence: List[int], t: int) -> List[int]:
        """
        Reverse the balancing operation by flipping first t symbols again.

        Args:
            sequence: Balanced sequence
            t: Balancing index

        Returns:
            Original unbalanced sequence
        """
        return self.flip_sequence(sequence, t)


if __name__ == "__main__":
    # Test GC balancing
    print("Testing gc_balance.py")
    print("-" * 50)

    balancer = GCBalancer(epsilon=0.05)

    test_cases = [
        [0, 0, 0, 0, 1, 1, 1, 1],  # All non-GC
        [2, 2, 2, 2, 3, 3, 3, 3],  # All GC
        [0, 1, 2, 3, 0, 1, 2, 3],  # Balanced
        [0, 0, 0, 1, 1, 1, 2, 3],  # Slightly unbalanced
    ]

    for i, original in enumerate(test_cases, 1):
        print(f"\nTest Case {i}")
        print(f"Original:     {original}")
        print(f"GC-content:   {balancer.gc_content(original):.2%}")
        print(f"Is balanced:  {balancer.is_balanced(original)}")

        balanced, t = balancer.balance(original)
        print(f"Balanced:     {balanced}")
        print(f"Flip index t: {t}")
        print(f"GC-content:   {balancer.gc_content(balanced):.2%}")
        print(f"Is balanced:  {balancer.is_balanced(balanced)}")

        # Test index suffix encoding/decoding
        suffix = balancer.create_index_suffix(t, len(original))
        print(f"Index suffix: {suffix}")
        t_decoded = balancer.decode_index_suffix(suffix)
        print(f"Decoded t:    {t_decoded}")
        print(f"Match:        {t == t_decoded}")

        # Test unbalancing
        unbalanced = balancer.unbalance(balanced, t)
        print(f"Unbalanced:   {unbalanced}")
        print(f"Recovered:    {original == unbalanced}")

    print("\n" + "-" * 50)
    print("All tests completed")