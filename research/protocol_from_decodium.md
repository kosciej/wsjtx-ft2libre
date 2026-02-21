# FT2 Protocol Specification

Derived from publicly available source code analysis.

## Core Parameters

| Parameter | Value | Description |
|-----------|-------|-------------|
| KK        | 91    | Information bits (77 message + 14 CRC) |
| ND        | 87    | Data symbols |
| NS        | 16    | Sync symbols (4 x 4-symbol Costas arrays) |
| NN        | 103   | Total sync + data symbols (NS+ND) |
| NN2       | 105   | Total channel symbols including 2 ramp-up/ramp-down |
| NSPS      | 288   | Samples per symbol at 12000 S/s |
| NZ        | 29664 | Sync+data samples (NSPS*NN = 288*103) |
| NZ2       | 30240 | Total samples including ramps (NSPS*NN2 = 288*105) |
| NMAX      | 45000 | Max samples in iwave (3.75*12000) |
| NFFT1     | 1152  | FFT length for symbol spectra |
| NSTEP     | 288   | Coarse time-sync step (= NSPS) |
| NDOWN     | 9     | Downsample factor (12000/9 = 1333.33 S/s) |

## Modulation

- **4-GFSK** (Gaussian Frequency Shift Keying with 4 tones)
- **Tones**: {0, 1, 2, 3}
- **Modulation index h = 1.0**
- **BT = 1.0** (Gaussian filter bandwidth-time product)
- **Baud rate**: 12000/288 = **41.667 Hz**
- **Tone spacing**: h * baud = 1.0 * 41.667 = **41.667 Hz**
- **Signal bandwidth**: 3 * 41.667 = 125 Hz center-to-center, ~167 Hz effective with GFSK spreading

## Gray Code Mapping

Bit pairs are mapped to tones using Gray code:

| Bits | Tone |
|------|------|
| 00   | 0    |
| 01   | 1    |
| 11   | 2    |
| 10   | 3    |

## FEC / Error Correction

- **LDPC(174,91)** — same code as FT8
- 77 message bits + 14 CRC bits = 91 information bits
- Encoded to 174 coded bits = 87 data symbols (2 bits/symbol)

## Message Format

- **77-bit message** — same packing as FT8
- **Scrambling vector (rvec)**: applied to message bits before encoding:
  ```
  0,1,0,0,1,0,1,0,0,1,0,1,1,1,1,0,1,0,0,0,1,0,0,1,1,0,1,1,0,
  1,0,0,1,0,1,1,0,0,0,0,1,0,0,0,1,0,1,0,0,1,1,1,1,0,0,1,0,1,
  0,1,0,1,0,1,1,0,1,1,1,1,1,0,0,0,1,0,1
  ```

## Frame Structure

Total: **105 channel symbols** (103 sync+data + 2 ramp symbols)

Layout: `r1 + s4 + d29 + s4 + d29 + s4 + d29 + s4 + r1`

Where:
- `r1` = 1 ramp-up/ramp-down symbol
- `s4` = 4-symbol Costas sync array
- `d29` = 29 data symbols

Detailed symbol positions (0-indexed into the 103 sync+data array):

| Position | Type | Content |
|----------|------|---------|
| 0-3      | Sync | Costas A: [0, 1, 3, 2] |
| 4-32     | Data | 29 data symbols |
| 33-36    | Sync | Costas B: [1, 0, 2, 3] |
| 37-65    | Data | 29 data symbols |
| 66-69    | Sync | Costas C: [2, 3, 1, 0] |
| 70-98    | Data | 29 data symbols |
| 99-102   | Sync | Costas D: [3, 2, 0, 1] |

Plus ramp-up symbol at the start and ramp-down symbol at the end (total 105).

### Costas Arrays (4x4)

Four distinct 4-element Costas-like arrays using tones {0,1,2,3}:

| Array | Tones      |
|-------|------------|
| A     | [0,1,3,2]  |
| B     | [1,0,2,3]  |
| C     | [2,3,1,0]  |
| D     | [3,2,0,1]  |

Note: These are each permutations of {0,1,2,3}. Each array is 4 symbols long
(not 7 like FT8 Costas arrays).

## Timing

- **T/R period**: 3.75 seconds
- **TX duration**: 105 * 288/12000 = **2.52 seconds**
- **Guard time**: 3.75 - 2.52 = 1.23 seconds
- **TX delay**: 500 ms
- **Max samples**: 45000 = 3.75 * 12000

## Waveform Generation

- GFSK pulse shaping with BT=1.0
- Pulse extends over 3 symbol periods (3*NSPS samples)
- Phase is continuous (accumulated dphi)
- Ramp-up: first NSPS samples shaped with cosine taper `(1-cos)/2`
- Ramp-down: last NSPS samples shaped with cosine taper `(1+cos)/2`
- At 48000 Hz TX sample rate: NSPS=4*288=1152 samples/symbol

## Decoding Pipeline

### 1. Candidate Detection
- Nuttall-windowed FFT (NFFT1=1152) stepping by NSPS
- Spectral averaging and baseline estimation
- Peak detection in smoothed spectrum
- Frequency range: 200-4910 Hz
- Frequency offset correction: `-1.5 * 12000/NSPS`

### 2. Downsampling
- Factor 9 downsampling: 12000 Hz -> 1333.33 Hz
- Bandpass filter: flat bandwidth = 4*baud, transition = 0.5*baud
- Frequency-shifted to baseband

### 3. Sync Detection
- Correlates against all 4 Costas arrays simultaneously
- Expected positions: i0, i0+33*NSS, i0+66*NSS, i0+99*NSS (where NSS=NSPS/NDOWN=32)
- Coarse search: step=4, range=-688 to 2024
- Fine search: step=1, +/-5 around best
- Frequency search: +/-12 Hz coarse (step 3), +/-4 Hz fine (step 1)

### 4. Bit Metrics
- FFT each symbol (NSS=32 samples) to extract tone powers
- 3 coherence modes: 1-symbol, 2-symbol, 4-symbol integration
- Gray-code-aware soft bit computation
- Sync quality check: requires >= 4/16 correct hard sync symbols

### 5. LDPC Decoding
- Up to 5 metric passes (3 coherence modes + best-of + average)
- Additional AP (a priori) passes depending on QSO progress
- LDPC decode with OSD (ordered statistics decoding)
- Max iterations: 40
- Signal subtraction after successful decode for multi-decode

## Fox/Hound Mode

- Multi-slot transmission: up to 3 slots
- Slot spacing: **200 Hz** (vs 60 Hz for FT8)
- Signal BW per slot: ~167 Hz
- Fox waveform generated at 48000 Hz, NSPS=4*288=1152
- Spectral filtering with 50 Hz transition width

## Sample Rates

| Context | Sample Rate | NSPS |
|---------|-------------|------|
| RX (decoder input) | 12000 Hz | 288 |
| RX (after downsample) | 1333.33 Hz | 32 |
| TX (waveform) | 48000 Hz | 1152 |

## Sync Quality Thresholds

- Minimum sync power: 0.90
- Minimum correct hard sync symbols: 4 out of 16
- Refined sync check threshold: 15 out of 32 hard-decoded sync bits
