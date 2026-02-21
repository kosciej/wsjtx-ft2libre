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

## Impact on ft2libre (original errors, now fixed)

The existing `wsjtx/lib/ft2libre/` code was updated from wrong 8-GFSK to correct
4-GFSK parameters:

| Parameter | ft2libre (was wrong) | Correct |
|-----------|----------------------|---------|
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

These parameters were corrected in a first implementation pass. However, direct
comparison with the decodium3 reference implementation revealed further
incompatibilities documented below.

---

## ft2libre vs Decodium3: Remaining Incompatibilities

Direct code comparison between our `wsjtx/lib/ft2libre/` and the reference
`decodium3-build/lib/ft2/` implementation. These are the differences that prevent
our decoder from being compatible with the reference FT2 encoder/decoder.

### 1. NHSYM Calculation (Minor)

| | ft2libre | decodium3 |
|---|---|---|
| Formula | `NMAX/NSTEP - 3` = 153 | `(NMAX-NFFT1)/NSTEP` = 152 |

Decodium3 accounts for the FFT window length when computing the number of valid
symbol spectra. Our formula yields one extra spectrum that may read past the
valid audio data. **Fix**: use `(NMAX-NFFT1)/NSTEP`.

### 2. Waveform Generation — Dummy Symbols (TX Incompatibility)

**Our code** (`gen_ft2libre_wave.f90:46-47`) adds "dummy symbols" extending the
first and last tone values into the ramp regions:

```fortran
! ft2libre EXTRA code — NOT in reference:
dphi(0:2*nsps-1) = dphi(0:2*nsps-1) + dphi_peak*itone(1)*pulse(nsps+1:3*nsps)
dphi(nsym*nsps:...) = dphi(nsym*nsps:...) + dphi_peak*itone(nsym)*pulse(1:2*nsps)
```

**Decodium3** does NOT add any dummy symbols. The ramp regions have zero frequency
deviation (carrier only). This causes our TX waveform to have a different phase
trajectory. While the ramp amplitude is tapered by cosine envelope, the phase
difference propagates into the first and last data symbols.

**Fix**: Remove the two dummy-symbol dphi additions entirely.

### 3. Waveform Generation — BT Parameter (Minor)

Our code passes `bt` as a parameter to `gen_ft2libre_wave`. Decodium3 hardcodes
`hmod=1.0` and uses `gfsk_pulse(1.0, tt)` — BT is always 1.0.

**Fix**: Hardcode BT=1.0, remove as parameter.

### 4. Sync Correlation — Fundamentally Different Approach (RX Incompatibility)

**Decodium3** (`sync2d.f90`):
- Builds 4 continuous-phase reference arrays (`csynca/b/c/d`), each 2*NSS=64 samples
- Phase accumulates continuously across all 4 symbols within each Costas array
- Uses NSS/2=16 samples per symbol (stride-2 access: `cd0(i1:i1+4*NSS-1:2)`)
- Applies `fac=1/(2*NSS)` scaling in the power function
- Returns `p(z1)+p(z2)+p(z3)+p(z4)` where p = abs(z*fac) (not squared)

**ft2libre** (`sync8d_ft2libre.f90`):
- Builds per-tone reference `csync(0:3, NSPSD)` with independent phase per tone
- Correlates each of 16 sync symbols independently (4 per Costas array)
- Uses every sample (no stride): `cd0(i1:i1+NSPSD-1)`
- No scaling factor in power function
- Accumulates 16 separate p(z) values where p = |z|² (squared, not abs)

The decodium3 approach gives **coherent gain** across the 4 symbols of each Costas
array — the complex correlation sums coherently before taking magnitude. Our approach
sums powers non-coherently (16 separate correlations). This means:
- Different sync threshold magnitudes (not directly comparable)
- ~6 dB worse sensitivity in our approach at weak signals
- Different power function: decodium3 uses `abs()`, we use `|z|²`

**Fix**: Rewrite sync8d_ft2libre to match sync2d exactly:
1. Build 4 continuous-phase reference arrays of length 2*NSS
2. Correlate with stride-2 (`cd0(i1:i1+4*NSS-1:2)`)
3. Use `fac=1/(2*NSS)` and `p(z)=abs(z*fac)` (not squared)
4. Sum 4 correlations (one per Costas array), not 16

### 5. Candidate Search Strategy (Different Approach)

**Decodium3**: Two-stage approach:
1. `getcandidates2()` — spectral peak detection with Nuttall window, returns (freq, peak_height)
2. For each candidate, 3-segment DT search using `sync2d()` with coarse/fine passes

**ft2libre**: Single-stage `sync_ft2libre()` — combined Costas-correlation spectral
search that returns (freq, dt, sync_quality) directly.

This is an algorithmic choice that affects sensitivity but not protocol compatibility.
The decodium3 approach has wider DT search range and more thorough frequency search.

**Fix**: Replace `sync_ft2libre()` with `getcandidates2()`-style spectral peak search,
then use `sync2d()` for fine time/frequency refinement. Or keep our approach but
ensure the DT search range covers at least -688 to 2024 (downsampled samples).

### 6. DT Search Range (RX Sensitivity)

| | ft2libre | decodium3 |
|---|---|---|
| Range | i0-10 to i0+10 | -688 to 2024 (3 segments) |
| Strategy | Single narrow search around initial estimate | 3-segment coarse/fine with segmented exit |

**Fix**: Implement 3-segment search matching decodium3 ranges:
- Segment 1: ibmin=216, ibmax=1120
- Segment 2: ibmin=1120, ibmax=2024
- Segment 3: ibmin=-688, ibmax=216
With early exit: skip segment 2/3 if sync < 0.9 or sync < segment 1's best.

### 7. Bit Metrics — Scope and Extraction (RX Incompatibility)

**Decodium3** (`get_ft2_bitmetrics.f90`):
- Computes metrics over ALL 103 symbols (including sync): `ks=1, NN-nsym+1, nsym`
- Produces `bitmetrics(2*NN, 3)` = 206 bit positions
- Then extracts data-only bits in `ft2_decode.f90`:
  ```
  llra(1:58) = bitmetrics(9:66, 1)       ! Skip sync bits 1-8
  llra(59:116) = bitmetrics(75:132, 1)    ! Skip sync bits 67-74
  llra(117:174) = bitmetrics(141:198, 1)  ! Skip sync bits 133-140
  ```
- Edge patching: `bitmetrics(205:206,2) = bitmetrics(205:206,1)` etc.
- Additional hard sync quality check: 15/32 bits must match Gray-coded Costas

**ft2libre** (`ft2libreb.f90`):
- Computes metrics over data symbols only (3 blocks of 29)
- Produces 174 bit metrics directly: `i32=1+(k-1)*2+(iblock-1)*58`
- No sync quality check on hard-decoded bit metrics
- Has 5 metric variants (a,b,c,d,e) where d=confidence-ratio, e=best-of

The multi-symbol coherent integration boundaries differ. In decodium3, the 4-symbol
integration window can span across the sync/data boundary (e.g. symbols 1-4 include
3 sync + 1 data), giving slightly different LLR values near boundaries.

**Fix**: Rewrite to match decodium3 exactly:
1. Compute bitmetrics over all 103 symbols → 206 bit positions
2. Apply edge patching for multi-symbol modes
3. Extract data bits by skipping sync positions (9:66, 75:132, 141:198)
4. Add hard sync quality check (15/32 threshold)
5. Change metric combination to match decodium3's best-of and average

### 8. Signal Subtraction (RX Incompatibility)

| | ft2libre | decodium3 |
|---|---|---|
| nstart | `dt*12000+1` | `dt*12000+1-NSPS` |
| NFILT | 600 | 700 |
| End correction | Yes (cos² edge correction) | No |
| Reference gen | `gen_ft2libre_wave(itone,NN,NSPS,1.0,12000.0,...)` | `gen_ft2wave(itone,103,NSPS,12000.0,...)` |
| Common block | `heap2libre` | `heap2` |

**nstart offset**: Decodium3 starts one full symbol (NSPS samples) earlier than we do.
This aligns the subtraction window differently relative to the detected signal.

**Fix**:
1. Change nstart to `dt*12000+1-NSPS`
2. Change NFILT from 600 to 700
3. Remove end correction
4. Match the reference waveform generation call signature

### 9. AP (A Priori) Decoding — Completely Missing (RX Sensitivity)

Decodium3 has full AP decoding with 6 AP types based on QSO progress state,
contest modes (NA_VHF, EU_VHF, FIELD DAY, RTTY, WW_DIGI, FOX, HOUND), and known
callsigns. This provides ~3-6 dB effective sensitivity gain for known QSO partners.

Our ft2libre has NO AP decoding (all passes use iaptype=0, apmask=0).

**Fix**: Port the full AP decoding logic from `ft2_decode.f90` lines 96-183
(setup) and 353-426 (AP pass logic). This requires:
- Precomputed AP bit patterns for CQ variants and known callsigns
- rvec-scrambled AP patterns for mcq/mcqru/mcqfd/mcqtest/mcqww/mrrr/m73/mrr73
- QSO progress state tracking (nQSOProgress 0-5)
- nappasses and naptypes configuration tables

### 10. Decoder Parameters

| | ft2libre | decodium3 |
|---|---|---|
| max_iterations | 30 | 40 |
| maxosd | 2 (or -1 for shallow) | 3 (or 4 near nfqso) |
| norder/ndeep | 2 | 3 |
| Subtraction passes | depends on lsubtract flag | 3 passes (ndepth≥2) |
| scalefac | 2.83 | 2.83 |

**Fix**: Change max_iterations to 40, maxosd to 3 (4 near nfqso), norder/ndeep to 3.

### 11. DT Reporting

| | ft2libre | decodium3 |
|---|---|---|
| Formula | `(ibest-1)*dt2` | `ibest/1333.33 - 0.5` |

The -0.5 offset in decodium3 accounts for the TX delay convention.

**Fix**: Use `ibest/1333.33 - 0.5`.

### 12. Normalization of Downsampled Signal

Decodium3 normalizes `cd2` after downsampling:
```fortran
sum2 = sum(cd2*conjg(cd2)) / (real(NMAX)/real(NDOWN))
if(sum2.gt.0.0) cd2 = cd2/sqrt(sum2)
```

And again after extracting the signal region:
```fortran
sum2 = sum(abs(cb)**2) / (real(NSS)*NN)
if(sum2.gt.0.0) cb = cb/sqrt(sum2)
```

Our ft2libre does not normalize after downsampling. This affects the absolute
magnitude of sync and bit metric values.

**Fix**: Add normalization after ft2libre_downsample.

### Summary: Priority Order for Fixes

1. **Waveform generation** (remove dummy symbols) — TX compatibility
2. **Sync correlation** (continuous-phase, stride-2, abs not squared) — RX detection
3. **Bit metrics** (all 103 symbols, extract data bits, edge patching) — RX decoding
4. **Signal subtraction** (nstart-NSPS, NFILT=700, no end correction) — multi-decode
5. **Normalization** of downsampled signal — magnitude calibration
6. **DT search range** — wider coverage
7. **Decoder parameters** (max_iterations=40, maxosd=3) — sensitivity
8. **AP decoding** — sensitivity for known QSO partners
9. **DT reporting** formula — UI display
10. **NHSYM** formula — minor correctness
