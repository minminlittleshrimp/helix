# HELIX: Soft Explanation of DNA Storage Methods

## Why Do We Need These Methods?

Imagine you want to store digital data (like a photo or document) in DNA molecules.
But DNA isn't like a hard drive - it has **biological constraints** that affect how
reliably it can be synthesized, stored, and read back

### The Three Main Problems:

1. **Homopolymer Runs** (like AAAAAA or GGGGGG)
   - DNA sequencing machines get confused when the same base repeats many times
   - Like trying to count "how many A's?" in AAAAAAAAAA - you might lose track!
   - **Solution**: Limit runs to maximum length `ell` (typically 3-4)

2. **GC-Content Imbalance**
   - DNA has 4 bases: A, T, C, G
   - C and G bond together more strongly (3 hydrogen bonds vs 2)
   - If you have too many C's and G's, the DNA is too stable
   - If you have too many A's and T's, it's too weak
   - **Ideal**: ~50% GC content (balanced)
   - **Solution**: Keep GC-content within 50% ± epsilon (e.g., 45%-55%)

3. **Sequencing Errors**
   - DNA synthesis/sequencing isn't perfect
   - Bases can be inserted, deleted, or substituted by mistake
   - **Solution**: Add error correction capability

---

## The Encoding Pipeline

Think of it like packing a suitcase for a flight:

```
Your Data → Pack it → Check weight → Add tags → Ready to fly!
Binary    → RLL     → GC-balance → Error Code → DNA
```

---

## Method B: Runlength-Limited (RLL) Encoding

**Problem**: After converting binary to DNA, you might accidentally create long
runs like "00000" (which becomes AAAAA in DNA).

**Method B's Solution - "Pointer Replacement"**:

### The Idea:
1. First, append a termination marker '0' to your data
2. Look for forbidden patterns (e.g., "000" = three zeros in a row)
3. When you find one, **replace it with a pointer** "Re" (like a shortcut)
4. Keep doing this until no forbidden patterns remain

### Example:
```
Original:    [1, 0, 0, 0, 2]  ← Has "000" (forbidden if ell=3)
Add marker:  [1, 0, 0, 0, 2, 0]
Replace:     [1, 1, 1, 2, 0]  ← Replaced "000" with pointer "11"
```

### Why It Works:
- The pointer "11" tells the decoder: "This used to be 000"
- During decoding, scan from right to left
- When you see "11", expand it back to "000"
- Remove the termination marker

**Real-World Analogy**: Like compressing a file - replace repeated patterns with
shorter codes, but make sure you can uncompress it perfectly!

---

## Method D: GC-Content Balancing

**Problem**: Your DNA sequence might have too many A's and T's (low GC%) or too
many C's and G's (high GC%).

**Method D's Solution - "Prefix Flipping"**:

### The Flipping Rule:
There's a special rule that swaps non-GC bases with GC bases:
```
f(0) = 2    →    A <-> C
f(1) = 3    →    T <-> G
f(2) = 0    →    C <-> A
f(3) = 1    →    G <-> T
```

### The Algorithm:
1. **Generate search set**: Try different "flip positions" t = {0, 2εn, 4εn, ..., n}
2. **For each position t**:
   - Flip the first t symbols using the rule above
   - Check if GC-content is now balanced (between 45%-55%)
3. **When balanced**: Stop and remember the flip position t
4. **Create suffix**: Encode t as a small suffix so decoder knows where you flipped

### Example:
```
Original: [0, 0, 0, 0]  -> AAAA (0% GC - BAD!)
Try t=2:  [2, 2, 0, 0]  -> CCAA (50% GC - GOOD!)

Encode t=2 as suffix: [2] -> [1, 0] in quaternary -> [1, 3, 0, 2] interleaved
Final: [2, 2, 0, 0] + [1, 3, 0, 2] → CCAATGAC
```

### Why Interleaved Suffix?
- We encode t and f(t) interleaved: [t[0], f(t[0]), t[1], f(t[1]), ...]
- This suffix itself is GC-balanced!
- Like adding a balanced "address label" to tell decoder where you flipped

**Real-World Analogy**: Like balancing a see-saw - you add weights to one side
(flip some bases) until it's level, then mark where you added them!

---

## Varshamov-Tenengolts (VT) Error Correction

**Problem**: What if one base gets deleted, inserted, or changed during sequencing?

**VT Solution - "Syndrome Fingerprint"**:

### The Syndrome Formula:
```
Syn(x) = (1.x1 + 2.x2 + 3.x3 + ... + n.xn) mod 2n
```

Think of it as a **weighted checksum**:
- Multiply each symbol by its position
- Add them all up
- Take remainder when divided by 2n

### Why This Works:
Each type of error creates a **unique signature**:

1. **Substitution** (wrong base): Changes the value but not the length
2. **Deletion** (missing base): Shifts all positions to the left
3. **Insertion** (extra base): Shifts all positions to the right

The syndrome can detect which type and approximately where it happened!

### Additional Checksum:
We also compute: `sum(all symbols) mod 4`

This helps distinguish between error types.

### Example:
```
Original:  [1, 2, 3, 0]
Syndrome:  (1×1 + 2×2 + 3×3 + 4×0) mod 8 = 14 mod 8 = 6
Checksum:  (1+2+3+0) mod 4 = 2

If corrupted: [1, 2, 0]  ← deleted the 3
New Syndrome: (1×1 + 2×2 + 3×0) mod 6 = 5
New Checksum: (1+2+0) mod 4 = 3

Different syndrome + different checksum = Error detected!
```

**Real-World Analogy**: Like a shipping barcode - if the package is damaged, the
barcode won't scan correctly, alerting you to the problem!

---

## Putting It All Together: The Complete Pipeline

### Encoding (Binary → DNA):

```
Step 1: Binary → Quaternary
  "11010011" → [3, 1, 0, 3]
  (Group every 2 bits: 11=3, 01=1, 00=0, 11=3)

Step 2: Differential Encoding
  [3, 1, 0, 3] → [3, 2, 3, 3]
  (First stays, rest are differences mod 4)

Step 3: RLL Encoding (Method B)
  Replace forbidden "000" patterns with pointers

Step 4: GC-Balancing (Method D)
  Find flip position t that balances GC-content

Step 5: Add Index Suffix
  Encode t in interleaved format

Step 6: Add Error Correction
  Compute VT syndrome and checksum, append as suffix

Step 7: Convert to DNA
  [3, 1, 0, 3] → GTAG
```

### Decoding (DNA → Binary):

```
Step 1: DNA → Quaternary
  GTAG -> [3, 1, 0, 3]

Step 2: Extract & Verify Error Correction
  Check syndrome matches (if error, can locate it)

Step 3: Extract Index Suffix
  Decode t from interleaved suffix

Step 4: Reverse GC-Balancing
  Unflip first t symbols

Step 5: RLL Decoding
  Replace pointers back to "000" patterns

Step 6: Differential Decoding
  Reverse the difference operation

Step 7: Quaternary -> Binary
  [3, 1, 0, 3] -> "11010011"
```

---

## Key Advantages of This Approach

1. **High Efficiency**: Approaches theoretical capacity (~1.92 bits per base)
2. **Low Complexity**: Simple operations (no complex math)
3. **Error Resilience**: Can detect and correct single errors
4. **Practical**: Works with real DNA synthesis/sequencing constraints

---

## Analogy Summary

Think of the entire system as **preparing a message for international shipping**:

1. **Method B (RLL)**: Pack items so nothing is too bulky (no long runs)
2. **Method D (GC-balance)**: Balance the weight distribution (50% GC)
3. **VT Error Correction**: Add tracking barcode (detect if damaged)
4. **Index Suffix**: Add return address label (how to unpack)

The goal: Get your data safely through the "DNA shipping system"!

---

## Why "Capacity-Approaching"?

The paper proves these methods achieve rates very close to the **theoretical maximum**
possible given the constraints. It's like getting 95% of the absolute best possible efficiency
- much better than previous methods that only got 70-80%!

---

## Practical Example from Your Demo:

```
Input:  "11010011" (8 bits)
Output: "GCGGAACACACTG" (13 bp = 26 bits capacity)
Efficiency: 8/26 = 30.77%

Why not 50%? Because we added:
- RLL redundancy (pointers)
- GC-balance suffix (index t)
- Error correction suffix (syndrome + checksum)

These are the "packing materials" to protect your data!
```

The actual usable rate approaches ~45-48% for longer sequences, which is excellent
given all three constraints!