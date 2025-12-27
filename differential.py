"""
differential.py - HELIX DNA Storage System
Differential encoding/decoding for quaternary sequences
"""

from typing import List


def differential_encode(sequence: List[int]) -> List[int]:
    """
    Apply differential encoding to quaternary sequence.

    Formula:
        y[0] = x[0]
        y[i] = (x[i] - x[i-1]) mod 4  for i > 0

    This encoding helps with runlength constraint by converting
    repeating symbols into zeros.

    Args:
        sequence: Input quaternary sequence

    Returns:
        Differentially encoded sequence

    Example:
        [2, 2, 2, 3] -> [2, 0, 0, 1]
        (2, 2-2=0, 2-2=0, 3-2=1)
    """
    if not sequence:
        return []

    encoded = [sequence[0]]
    for i in range(1, len(sequence)):
        diff = (sequence[i] - sequence[i-1]) % 4
        encoded.append(diff)

    return encoded


def differential_decode(encoded: List[int]) -> List[int]:
    """
    Reverse differential encoding to recover original sequence.

    Formula:
        x[0] = y[0]
        x[i] = (x[i-1] + y[i]) mod 4  for i > 0

    Args:
        encoded: Differentially encoded sequence

    Returns:
        Original quaternary sequence

    Example:
        [2, 0, 0, 1] -> [2, 2, 2, 3]
    """
    if not encoded:
        return []

    decoded = [encoded[0]]
    for i in range(1, len(encoded)):
        symbol = (decoded[i-1] + encoded[i]) % 4
        decoded.append(symbol)

    return decoded


def verify_differential_pair(original: List[int], encoded: List[int]) -> bool:
    """
    Verify that encoding/decoding are inverses of each other.

    Args:
        original: Original sequence
        encoded: Encoded sequence

    Returns:
        True if decode(encode(original)) == original
    """
    if len(original) != len(encoded):
        return False

    decoded = differential_decode(encoded)
    return original == decoded


if __name__ == "__main__":
    # Test differential encoding
    print("Testing differential.py")
    print("-" * 50)

    test_cases = [
        [2, 2, 2, 3],
        [0, 1, 2, 3],
        [3, 3, 3, 3],
        [1, 0, 3, 2, 1],
    ]

    for original in test_cases:
        print(f"\nOriginal:  {original}")
        encoded = differential_encode(original)
        print(f"Encoded:   {encoded}")
        decoded = differential_decode(encoded)
        print(f"Decoded:   {decoded}")
        print(f"Match: {original == decoded}")

    print("\n" + "-" * 50)
    print("All tests completed")