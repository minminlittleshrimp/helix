# HELIX - High-Efficiency Lossless Information eXchange

A Python implementation of capacity-approaching constrained codes for DNA-based data storage, based on the paper:

**"Capacity-Approaching Constrained Codes with Error Correction for DNA-Based Data Storage"**
by Tuan Thanh Nguyen, Kui Cai, Kees A. Schouhamer Immink, and Han Mao Kiah

## Overview

HELIX transforms binary data into DNA sequences while satisfying critical biochemical constraints:

- **Runlength Constraint (RLL)**: Prevents homopolymer runs longer than `ell` nucleotides, reducing sequencing errors
- **GC-Content Balancing**: Maintains GC-content near 50% +- `epsilon` for DNA stability
- **Error Correction**: Incorporates Varshamov-Tenengolts (VT) syndromes for single-edit error detection

## Architecture

The system consists of modular components, each handling one step of the encoding/decoding pipeline:

```
mapping.py          - Binary <-> Quaternary <-> DNA conversions
differential.py     - Differential encoding/decoding
rll_constraint.py   - Method B: Runlength limiting
gc_balance.py       - Method D: GC-content balancing
error_correction.py - VT syndrome & error detection
analyzer.py         - Sequence analysis & validation
helix.py           - Main codec integration
examples.py        - Demo scripts
```

## Installation

```bash
# Clone or download the repository
# No external dependencies required - uses only Python standard library

# Verify installation
python helix.py
```

## Quick Start

### Basic Usage

```python
from helix import HelixCodec

# Initialize codec with constraints
codec = HelixCodec(ell=3, epsilon=0.05)

# Encode binary data to DNA
binary_data = "11010011"
dna_sequence = codec.encode(binary_data)
print(f"DNA: {dna_sequence}")

# Decode DNA back to binary
decoded = codec.decode(dna_sequence)
print(f"Decoded: {decoded}")
```

### Encoding Text

```python
from helix import HelixCodec

codec = HelixCodec()

# Convert text to binary
message = "HELIX"
binary = ''.join(format(ord(c), '08b') for c in message)

# Encode to DNA
dna = codec.encode(binary)
print(f"Text '{message}' encoded as DNA: {dna}")
```

### Analyzing Sequences

```python
from analyzer import SequenceAnalyzer

analyzer = SequenceAnalyzer(ell=3, epsilon=0.05)

dna = "ATCGATCG"
analysis = analyzer.analyze_dna(dna)
analyzer.print_analysis(analysis)
```

## Encoding Pipeline

The HELIX encoding process follows these steps:

1. **Binary -> Quaternary**: Convert binary string to base-4 representation (2 bits -> 1 quaternary symbol)
2. **Differential Encoding**: Transform to facilitate runlength constraint enforcement
3. **RLL Encoding**: Apply Method B to prevent long homopolymer runs
4. **GC-Balancing**: Use Method D with prefix flipping to achieve 50% GC-content
5. **Index Suffix**: Append interleaved suffix encoding the balancing index
6. **Error Correction**: Add VT syndrome and checksum (optional)
7. **DNA Conversion**: Map quaternary to nucleotides (0->A, 1->T, 2->C, 3->G)

Decoding reverses these steps to recover the original binary data.

## Module Documentation

### mapping.py
Handles conversions between binary, quaternary, and DNA representations.

**Key functions:**
- `binary_to_quaternary(binary_str)` - Convert binary to quaternary
- `quaternary_to_dna(quaternary)` - Convert quaternary to DNA
- `dna_to_quaternary(dna)` - Convert DNA to quaternary
- `quaternary_to_binary(quaternary)` - Convert quaternary to binary

### differential.py
Implements differential encoding that transforms repeated symbols into zeros, facilitating RLL constraint enforcement.

**Key functions:**
- `differential_encode(sequence)` - Apply differential encoding
- `differential_decode(encoded)` - Reverse differential encoding

### rll_constraint.py
Enforces runlength-limited (RLL) constraint using Method B from the paper.

**Key class: `RLLCodec`**
- `encode(data)` - Encode to satisfy runlength constraint
- `decode(encoded)` - Decode RLL-encoded sequence
- `max_runlength(sequence)` - Calculate maximum runlength

### gc_balance.py
Achieves GC-content balancing through prefix flipping (Method D).

**Key class: `GCBalancer`**
- `balance(sequence)` - Find optimal flip index and balance sequence
- `create_index_suffix(t, n)` - Create suffix encoding flip index
- `decode_index_suffix(suffix)` - Decode flip index from suffix
- `flip_symbol(symbol)` - Apply flipping rule: f(0)=2, f(1)=3, f(2)=0, f(3)=1

### error_correction.py
Implements Varshamov-Tenengolts (VT) error correction framework.

**Key class: `VTErrorCorrection`**
- `compute_syndrome(sequence)` - Calculate VT syndrome
- `compute_checksum(sequence)` - Calculate simple checksum
- `detect_error(sequence, expected_syndrome, expected_checksum)` - Detect errors

### analyzer.py
Provides tools for analyzing DNA sequences and validating constraints.

**Key class: `SequenceAnalyzer`**
- `analyze_dna(dna)` - Complete sequence analysis
- `validate_constraints(dna)` - Check constraint satisfaction
- `compute_gc_content(dna)` - Calculate GC-content ratio
- `compute_max_runlength(dna)` - Find longest homopolymer run

### helix.py
Main codec integrating all components into complete pipeline.

**Key class: `HelixCodec`**
- `encode(binary_data, verbose)` - Complete encoding pipeline
- `decode(dna, verbose)` - Complete decoding pipeline
- `encode_with_analysis(binary_data)` - Encode and analyze
- `verify_roundtrip(binary_data)` - Test encode/decode correctness

## Examples

Run the examples script to see demonstrations of all features:

```bash
python examples.py
```

Examples include:
1. Basic encoding and decoding
2. Text message encoding
3. Constraint comparison with different parameters
4. Batch processing of multiple sequences
5. Error detection demonstration
6. Direct format conversions
7. Sequence analysis tools

## Parameters

### Codec Parameters
- **`ell`** (default: 3): Maximum allowed runlength (homopolymer limit)
  - Smaller values = stricter constraint, more robust to sequencing errors
  - Larger values = more efficient encoding, but longer homopolymer runs

- **`epsilon`** (default: 0.05): GC-content tolerance
  - Target: 50% +- epsilon
  - Smaller values = tighter GC-balance, better DNA stability
  - Larger values = more relaxed constraint, higher efficiency

- **`use_error_correction`** (default: True): Enable error correction suffix
  - Adds VT syndrome and checksum to sequence
  - Enables single-edit error detection

## Testing

Each module includes built-in tests. Run individual modules:

```bash
python mapping.py
python differential.py
python rll_constraint.py
python gc_balance.py
python error_correction.py
python analyzer.py
```

## Performance Characteristics

- **Efficiency**: Typically 40-60% depending on constraints
  - Measured as: (input bits) / (output DNA length × 2)
- **Constraint Satisfaction**: 100% for valid parameter choices
- **Error Tolerance**: Single insertion, deletion, or substitution detectable

## Limitations

Current implementation:
- Error **detection** implemented, full error **correction** requires extension
- Simplified suffix length handling in decoder (fixed assumptions)
- Not optimized for large-scale data (designed for demonstration)

## Paper Reference

This implementation follows the methods described in:

```
Nguyen, T. T., Cai, K., Immink, K. A. S., & Kiah, H. M.
"Capacity-Approaching Constrained Codes with Error Correction
for DNA-Based Data Storage"
```

Key concepts implemented:
- Method B: RLL constraint with pointer replacement
- Method D: GC-balancing through prefix flipping
- Corollary 24: Symbol selection for concatenation
- VT codes: Single-edit error syndrome calculation

## License

MIT

This is a research implementation for educational and academic purposes.

## Contributing

This is a demonstration implementation. For production use, consider:
- Optimizing suffix length encoding/decoding
- Implementing full VT error correction (not just detection)
- Adding support for longer sequences with segmentation
- Parallel processing for batch operations
- More sophisticated gamma selection algorithm

## Contact

For questions about the implementation or the underlying paper, please refer to the original publication.