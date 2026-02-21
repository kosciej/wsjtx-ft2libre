# Differences: Our Reverse-Engineering vs Known Protocol

Comparison between our blind signal analysis (`protocol.md`) and the confirmed
FT2 protocol specification (`protocol_from_decodium.md`).

## What We Got Right

### 4-GFSK Modulation
- **Our finding**: 4 tones, NOT 8 as claimed on ft2.it
- **Confirmed**: 4 tones {0,1,2,3}, Gray-coded
- This was our most important breakthrough — the ft2.it website is wrong about 8-GFSK

### Tone Spacing ~41.3 Hz
- **Our measurement**: 41.3 Hz (from instantaneous frequency histogram)
- **Actual**: h=1.0, baud=41.667 Hz, so tone spacing = 41.667 Hz
- Our measurement was within 0.9% — excellent agreement

### NSPS = 288
- **Our finding**: NSPS=288 was listed as our top candidate ("most likely on theoretical grounds")
- **Confirmed**: NSPS=288 exactly
- We correctly identified h≈1.0 as the likely modulation index

### Sample Rate 12000 Hz
- **Both agree**: 12000 Hz RX sample rate

### LDPC(174,91)
- **Our speculation**: "If LDPC(174,91) is used (as in FT8/FT4)..."
- **Confirmed**: LDPC(174,91) with KK=91, same as FT8

### 77-bit Message Format
- **Our question**: "What message packing format is used? (77-bit like FT8, or different?)"
- **Confirmed**: 77-bit packing, identical to FT8

## What We Got Wrong

### T/R Period
- **Our finding**: 3.8 seconds
- **Actual**: 3.75 seconds (NMAX=45000 = 3.75*12000)
- Close but not exact — our burst timing measurement had ~1.3% error

### TX Duration
- **Our finding**: ~2.16-2.32 seconds
- **Actual**: 2.52 seconds (105*288/12000)
- We significantly underestimated — likely due to burst detection thresholds
  cutting off the ramp-up/ramp-down portions

### "3 identical bursts per T/R period"
- **Our finding**: We observed 3 bursts in each capture
- **Actual**: Only 1 TX burst per T/R period (105 symbols = 2.52s within 3.75s cycle)
- What we saw as "3 bursts" was actually 3 separate T/R periods in our ~13s captures
  (13.47s / 3.75s ≈ 3.6 periods). We misinterpreted consecutive T/R cycles as
  repeated bursts within a single cycle.

### Symbol Count
- **Our finding**: ~72-96 symbols (uncertain)
- **Actual**: 103 sync+data symbols + 2 ramp = 105 total channel symbols

### Sync Pattern — Completely Wrong Approach
- **Our finding**: "No Costas arrays found", tried 7-symbol Costas arrays
- **Actual**: Uses four **4-symbol** Costas arrays (not 7-symbol!)
  - A=[0,1,3,2], B=[1,0,2,3], C=[2,3,1,0], D=[3,2,0,1]
  - Placed at positions 0, 33, 66, 99 in the 103-symbol frame
- We searched for 7-element Costas arrays (FT8 pattern) — wrong size!
  FT2 uses 4x4 permutations of {0,1,2,3}, matching the 4-tone alphabet

### Tone Frequencies
- **Our measurement**: ~1500.6, ~1541.9, ~1583.2, ~1624.5 Hz
- **Actual**: Tones are relative to the TX frequency, not absolute.
  The base frequency is user-configurable. Tone offsets
  are 0, 41.667, 83.333, 125.0 Hz from f0.
- Our measured frequencies were the actual tones in the specific captures,
  not a protocol property. The ~1500 Hz base was just the chosen TX frequency.

## What We Didn't Discover

### Frame Structure
- `r1 + s4 + d29 + s4 + d29 + s4 + d29 + s4 + r1`
- 4 sync blocks of 4 symbols each, evenly spaced every 33 symbols
- 3 data blocks of 29 symbols each
- 2 ramp symbols (cosine-tapered) at start and end

### Gray Code Mapping
- Bits 00→tone 0, 01→tone 1, 11→tone 2, 10→tone 3

### Scrambling Vector (rvec)
- 77-bit pseudorandom vector XORed with message bits before LDPC encoding

### BT Product
- BT = 1.0 for the Gaussian filter (we listed this as an open question)

### CRC
- 14-bit CRC (77 message bits + 14 CRC = 91 LDPC information bits)

### Downsample Factor
- Factor 9 (12000→1333.33 Hz) in the decoder

### Fox Mode
- Multi-slot Fox mode with 200 Hz slot spacing (vs 60 Hz for FT8)

## Assessment

Our blind analysis correctly identified the most important parameter — **4-GFSK
modulation with h≈1.0 and NSPS=288** — despite the ft2.it website claiming 8-GFSK.
The instantaneous frequency histogram technique was decisive for this.

The main failures were:
1. **Sync pattern size**: We assumed FT8-style 7-symbol Costas arrays instead of
   4-symbol arrays. With 4 tones, 4-symbol permutation arrays are the natural choice.
2. **Timing measurements**: Our burst detection was too aggressive, leading to
   underestimated TX duration and wrong T/R period.
3. **Burst interpretation**: Mistaking 3 consecutive T/R periods for 3 bursts
   within one period.

The cross-reference technique (comparing two captures with different messages)
was sound in principle but was hampered by incorrect NSPS candidates and
alignment optimization bias.

## Impact on ft2libre

The existing `wsjtx/lib/ft2libre/` code needs to be updated:

| Parameter | ft2libre (wrong) | Correct |
|-----------|-------------------|---------|
| Modulation | 8-GFSK | 4-GFSK |
| h | 0.75 | 1.0 |
| NSPS | 360 | 288 |
| Symbols | 79 | 103 + 2 ramp = 105 |
| Sync | 7-symbol Costas x3 | 4-symbol Costas x4 |
| Costas arrays | FT8 [3,1,4,0,6,5,2] | A=[0,1,3,2] B=[1,0,2,3] C=[2,3,1,0] D=[3,2,0,1] |
| Tones | 0-7 | 0-3 |
| BT | (FT8 default 2.0) | 1.0 |
| Data symbols | 58 | 87 |
| Sync symbols | 21 | 16 |
| Bits/symbol | 3 | 2 |
