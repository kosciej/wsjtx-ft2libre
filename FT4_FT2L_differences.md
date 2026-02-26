# FT4 vs FT2L: Implementation Differences & Protocol Specification

This document details the specifications of the **FT2L** protocol (a clean-room implementation compatible with "FT2") and provides a line-by-line technical comparison with **FT4**.

---

## 1. FT2L Protocol Specification

### Physical Layer
| Parameter | Value | Description |
|-----------|-------|-------------|
| **Modulation** | 4-GFSK | Gaussian Frequency Shift Keying |
| **Tone Spacing** | 20.833 Hz | Frequency difference between adjacent tones |
| **Symbol Rate** | 41.667 Baud | 2x the speed of FT4 |
| **Bandwidth** | ~150 Hz | Occupied spectral width |
| **T/R Cycle** | 3.75 seconds | Total duration of one transmission/reception period |
| **Symbol Duration** | 24 ms | Time duration of a single symbol (288 samples @ 12kHz) |
| **NSPS** | 288 | Samples per symbol at 12000 S/s |

### Frame Structure
| Component | Symbols | Description |
|-----------|---------|-------------|
| **Sync (Costas)** | 16 | Four 4x4 Costas arrays (Indices: 1-4, 34-37, 67-70, 100-103) |
| **Data Symbols** | 87 | Carries the 174-bit encoded message (2 bits/symbol) |
| **Total Channel Symbols** | 103 | 16 Sync + 87 Data |
| **Total with Ramps** | 105 | Includes power ramp-up and ramp-down symbols |
| **Effective Duration** | ~2.52s | Total time on air per burst |

### Coding & Scrambling
- **Information Bits**: 77 bits (standard WSJT-X message).
- **CRC**: 14 bits.
- **FEC**: LDPC (174, 91).
- **Scrambling**: Synchronous XOR with a 77-bit pseudorandom vector (`rvec`).
- **Tone Mapping**: Gray-coded (00=0, 01=1, 11=2, 10=3).

---

## 2. Technical Differences: FT4 vs. FT2L

### Core Parameters (`ft2l_params.f90`)
| Feature | FT4 | FT2L |
|---------|-----|------|
| **NSPS (12kHz)** | 576 | **288** |
| **NDOWN (Downsample)**| 18 (666.67 Hz) | **9 (1333.33 Hz)** |
| **NMAX (Samples)** | 72576 (6.0s) | **45000 (3.75s)** |
| **NFFT1 (Spectral)** | 2304 | **1152** |
| **NHSYM (Spectra)** | 125 | **152** |

### Implementation Logic Differences

#### 1. Timing & Modulator (`Modulator.cpp` / `mainwindow.cpp`)
- **Cycle Timing**: FT4 uses a 7.5s cycle; FT2L uses **3.75s**.
- **Modulator Start Delay**: FT4 starts at 300ms; FT2L was increased to **500ms** to ensure TCI/PTT stability during the rapid cycle.
- **Symbol Scaling**: FT4 `framesPerSymbol` is 576.0; FT2L is **288.0**.

#### 2. Decoding Strategy (`ft2l_decode.f90`)
- **Sampling Rate**: FT2L operates at **1333.33 S/s** after downsampling (twice the rate of FT4).
- **Synchronization Search**:
    - **ibstp (Coarse step)**: FT4 uses 4; FT2L uses **1** (finer search required due to shorter cycle).
    - **Search Ranges**: Ibmin/ibmax indices were doubled/halved appropriately to maintain time-domain coverage within the 3.75s window.
- **Thresholds (Sensitivity)**:
    - **syncmin**: 1.18 (FT4) → **1.05** (FT2L).
    - **smax check**: 1.2 (FT4) → **1.0** (FT2L).
    - **nsync_qual**: 20 (FT4) → **15** (FT2L).
- **AP Decoding**:
    - **Passes**: Increased from 3 to **5** base passes.
    - **Pass 4 & 5**: Added new passes using `llrd` and `llre` metrics derived from bitmetrics 4 and 5.
    - **AP Masking**: Fixed a critical bug in the baseline where masks were limited to 77 bits; FT2L now correctly masks up to **116** or **135** bits depending on `iaptype`.

#### 3. Memory & Stability
- **Common Blocks**: FT2L uses `common/heap2l/` to isolate its subtraction memory from FT4's `common/heap4/`.
- **Array Bounds**:
    - `ft2l_downsample.f90`: Padded `x(NMAX+2)` to prevent overflow during complex-to-real equivalence.
    - `getcandidates2l.f90`: Padded `x(NFFT1+2)` for the same reason.
- **Bit Metrics**: Added logic to explicitly handle symbol indices up to 206 (2 * NN) to prevent OOB access.

#### 4. GUI & Orchestration
- **Early Decoding**:
    - FT4: Typically one pass at symbol 50.
    - FT2L: **Two early passes** added at symbols **12** and **17** to provide sub-second decoding results.
- **Frequency Plan**: FT2L default frequencies are shifted (e.g., 14.084 MHz vs 14.080 MHz) to prevent interference with established FT4/FT8 segments.
