"""
mapping.py - HELIX DNA Storage System
Binary <-> Quaternary <-> DNA conversions
"""

from typing import List

# Nucleotide mapping: 0->A, 1->T, 2->C, 3->G
NUCLEOTIDE_MAP = {0: 'A', 1: 'T', 2: 'C', 3: 'G'}
REVERSE_MAP = {'A': 0, 'T': 1, 'C': 2, 'G': 3}


def binary_to_quaternary(binary_str: str) -> List[int]:
    """
    Convert binary string to quaternary (base-4) sequence.
    Each pair of bits becomes one quaternary symbol.

    Args:
        binary_str: Binary string (e.g., "11010011")

    Returns:
        List of quaternary symbols [0-3]

    Example:
        "1101" -> [3, 1] (binary 11=3, 01=1)
    """
    # Pad to even length if needed
    if len(binary_str) % 2 != 0:
        binary_str = '0' + binary_str

    quaternary = []
    for i in range(0, len(binary_str), 2):
        two_bits = binary_str[i:i+2]
        quaternary.append(int(two_bits, 2))

    return quaternary


def quaternary_to_binary(quaternary: List[int]) -> str:
    """
    Convert quaternary sequence back to binary string.

    Args:
        quaternary: List of quaternary symbols [0-3]

    Returns:
        Binary string with leading zeros removed

    Example:
        [3, 1] -> "1101"
    """
    binary = ''.join(format(q, '02b') for q in quaternary)
    # Remove leading zeros, but keep at least one '0'
    return binary.lstrip('0') or '0'


def quaternary_to_dna(quaternary: List[int]) -> str:
    """
    Convert quaternary sequence to DNA string.
    Mapping: 0->A, 1->T, 2->C, 3->G

    Args:
        quaternary: List of quaternary symbols [0-3]

    Returns:
        DNA string using alphabet {A, T, C, G}

    Example:
        [0, 1, 2, 3] -> "ATCG"
    """
    return ''.join(NUCLEOTIDE_MAP[q] for q in quaternary)


def dna_to_quaternary(dna: str) -> List[int]:
    """
    Convert DNA string to quaternary sequence.
    Mapping: A->0, T->1, C->2, G->3

    Args:
        dna: DNA string using alphabet {A, T, C, G}

    Returns:
        List of quaternary symbols [0-3]

    Example:
        "ATCG" -> [0, 1, 2, 3]
    """
    return [REVERSE_MAP[nucleotide] for nucleotide in dna.upper()]


def binary_to_dna(binary_str: str) -> str:
    """
    Direct conversion: Binary -> DNA

    Args:
        binary_str: Binary string

    Returns:
        DNA string
    """
    quaternary = binary_to_quaternary(binary_str)
    return quaternary_to_dna(quaternary)


def dna_to_binary(dna: str) -> str:
    """
    Direct conversion: DNA -> Binary

    Args:
        dna: DNA string

    Returns:
        Binary string
    """
    quaternary = dna_to_quaternary(dna)
    return quaternary_to_binary(quaternary)


if __name__ == "__main__":
    # Test the mapping functions
    print("Testing mapping.py")
    print("-" * 50)

    # Test case 1
    binary = "11010011"
    print(f"Binary: {binary}")
    quat = binary_to_quaternary(binary)
    print(f"Quaternary: {quat}")
    dna = quaternary_to_dna(quat)
    print(f"DNA: {dna}")

    # Reverse
    quat_back = dna_to_quaternary(dna)
    print(f"Quaternary (decoded): {quat_back}")
    binary_back = quaternary_to_binary(quat_back)
    print(f"Binary (decoded): {binary_back}")
    print(f"Match: {binary == binary_back}")

    print("\n" + "-" * 50)

    # Test case 2
    dna2 = "ATCGATCG"
    print(f"DNA: {dna2}")
    binary2 = dna_to_binary(dna2)
    print(f"Binary: {binary2}")
    dna2_back = binary_to_dna(binary2)
    print(f"DNA (re-encoded): {dna2_back}")
    print(f"Match: {dna2 == dna2_back}")