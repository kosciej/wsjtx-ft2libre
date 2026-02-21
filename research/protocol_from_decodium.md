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

- GFSK pulse shaping with BT=1.0, hmod=1.0 (hardcoded, NOT a parameter)
- Pulse extends over 3 symbol periods (3*NSPS samples)
- Phase is continuous (accumulated dphi)
- Ramp-up: first NSPS samples shaped with cosine taper `(1-cos)/2`
- Ramp-down: last NSPS samples shaped with cosine taper `(1+cos)/2`
- At 48000 Hz TX sample rate: NSPS=4*288=1152 samples/symbol

### Critical: No dummy symbols in ramp regions

The dphi computation iterates ONLY over the nsym data+sync symbols. The ramp
regions (first and last NSPS samples of the (nsym+2)*nsps output) receive NO
additional frequency deviation — they carry only the carrier frequency f0.
Specifically:

```fortran
! Reference gen_ft2wave.f90 — THIS IS THE COMPLETE dphi loop:
dphi=0.0
do j=1,nsym
   ib=(j-1)*nsps
   ie=ib+3*nsps-1
   dphi(ib:ie) = dphi(ib:ie) + dphi_peak*pulse(1:3*nsps)*itone(j)
enddo
! NO further dphi additions. Ramp samples have dphi = 0 (carrier only).
dphi = dphi + twopi*f0*dt    ! Then shift everything up by f0
```

The ramp cosine taper is applied to the waveform amplitude, not frequency:

```fortran
! Ramp-up: first nsps samples
cwave(1:nsps) = cwave(1:nsps) * (1.0-cos(twopi*(/(i,i=0,nsps-1)/)/(2.0*nsps)))/2.0
! Ramp-down: last nsps samples
k1=(nsym+1)*nsps+1
cwave(k1:k1+nsps-1) = cwave(k1:k1+nsps-1) * (1.0+cos(twopi*(/(i,i=0,nsps-1)/)/(2.0*nsps)))/2.0
```

## Decoding Pipeline

### 1. Candidate Detection (`getcandidates2.f90`)
- Nuttall-windowed FFT (NFFT1=1152) stepping by NSPS
- NHSYM = (NMAX-NFFT1)/NSTEP = (45000-1152)/288 = **152** symbol spectra
- Power spectrum: `s(i,j) = abs(cx(i))**2`
- Average spectrum computed, then 15-bin smoothed: `savsm(i)=sum(savg(i-7:i+7))/15`
- Baseline estimated via `ft2_baseline()`, spectrum normalized by baseline
- Peak detection: local maxima in smoothed spectrum with parabolic interpolation
- Frequency offset correction: `f_offset = -1.5 * 12000/NSPS`
- Frequency range: 200-4910 Hz
- Candidates returned as (frequency, spectral peak height) — NO time estimate
- Candidates near nfqso (within 20 Hz) placed first in list

### 2. Downsampling
- Factor 9 downsampling: 12000 Hz -> 1333.33 Hz
- Bandpass filter: flat bandwidth = 4*baud, transition = 0.5*baud
- Frequency-shifted to baseband

### 3. Sync Detection (`sync2d.f90`)

Sync detection uses **continuous-phase coherent correlation** across all 4 symbols
of each Costas array simultaneously.

#### Reference array construction (computed once, saved):

Four complex reference arrays `csynca/b/c/d`, each of length `2*NSS` (=64 samples),
are built with **continuous phase** across the 4 Costas symbols:

```fortran
! For each Costas array (a,b,c,d):
k=1; phi=0.0
do i=0,3                            ! 4 symbols in the Costas array
   dphi = 2*twopi*icos(i)/real(NSS) ! Phase increment for this tone
   do j=1,NSS/2                     ! Only NSS/2 = 16 samples per symbol
      csync(k) = cmplx(cos(phi), sin(phi))
      phi = mod(phi+dphi, twopi)    ! Phase is CONTINUOUS across symbols
      k = k+1
   enddo
enddo
! Result: 4*(NSS/2) = 2*NSS = 64 samples, phase continuous across all 4 symbols
```

Key details:
- Each reference array spans 4 symbols but uses only **every 2nd sample** (NSS/2
  per symbol = 16 samples), yielding 2*NSS=64 total reference samples
- Phase accumulates continuously from symbol to symbol within each Costas array
- A frequency tweak `ctwk(2*NSS)` is applied: `csync2 = ctwk * csynca`

#### Correlation:

The correlation uses **stride-2 indexing** into the downsampled signal:

```fortran
i1=i0                      ! Costas A at position 0
i2=i0+33*NSS               ! Costas B at position 33
i3=i0+66*NSS               ! Costas C at position 66
i4=i0+99*NSS               ! Costas D at position 99

! Each correlation spans 4*NSS samples but takes every 2nd sample:
z1 = sum(cd0(i1:i1+4*NSS-1:2) * conjg(csync2))   ! stride-2
```

Power function includes a scaling factor:

```fortran
fac = 1.0/(2.0*NSS)
p(z1) = (real(z1*fac)**2 + aimag(z1*fac)**2)**0.5  ! = abs(z1*fac)
sync = p(z1) + p(z2) + p(z3) + p(z4)
```

Edge handling: partial correlations when signal extends past buffer boundaries.

#### Search strategy (in `ft2_decode.f90`):

Three DT segments searched sequentially, each with 2-pass (coarse then fine):

```
Segment 1: ibmin=216,  ibmax=1120  (early-normal DT range)
Segment 2: ibmin=1120, ibmax=2024  (late DT range)
Segment 3: ibmin=-688, ibmax=216   (very early DT range)
```

Coarse pass: freq -12 to +12 Hz (step 3), time step 4
Fine pass: freq +/-4 Hz (step 1), time +/-5 (step 1)

Frequency tweaks precomputed: `ctwk2(2*NSS, -16:16)` via `twkfreq1()`

### 4. Bit Metrics (`get_ft2_bitmetrics.f90`)

#### Symbol extraction:

FFT each of the NN=103 symbols (NSS=32 samples each):

```fortran
do k=1,NN
   i1=(k-1)*NSS
   csymb = cd(i1:i1+NSS-1)
   call four2a(csymb, NSS, 1, -1, 1)   ! c2c FFT
   cs(0:3,k) = csymb(1:4)              ! Complex amplitudes for tones 0-3
   s4(0:3,k) = abs(csymb(1:4))         ! Power for tones 0-3
enddo
```

#### Sync quality check:

Hard-decode all 16 sync symbols, count correct:

```fortran
do k=1,4
   ip=maxloc(s4(:,k))         ! Costas A at positions 1-4
   if(icos4a(k-1).eq.(ip(1)-1)) is1=is1+1
   ip=maxloc(s4(:,k+33))      ! Costas B at positions 34-37
   if(icos4b(k-1).eq.(ip(1)-1)) is2=is2+1
   ip=maxloc(s4(:,k+66))      ! Costas C at positions 67-70
   if(icos4c(k-1).eq.(ip(1)-1)) is3=is3+1
   ip=maxloc(s4(:,k+99))      ! Costas D at positions 100-103
   if(icos4d(k-1).eq.(ip(1)-1)) is4=is4+1
enddo
nsync=is1+is2+is3+is4   ! 0-16, reject if < 4
```

#### Coherent integration over ALL 103 symbols:

The bit metric loop iterates over **ALL NN=103 symbols** (including sync), NOT just
data symbols. This produces `bitmetrics(2*NN, 3)` = 206 bit positions x 3 modes.

```fortran
do nseq=1,3                       ! 3 modes: 1-sym, 2-sym, 4-sym
   nsym = [1, 2, 4](nseq)
   nt = 4**nsym                    ! Hypothesis count
   do ks=1,NN-nsym+1,nsym          ! Step through ALL symbols
      ! For each hypothesis, compute coherent sum of cs() values
      ! Then extract soft bits using max-of-matching vs max-of-non-matching
      ipt = 1+(ks-1)*2             ! Bit position (1-indexed, 2 bits per symbol)
      do ib=0,ibmax
         bm = maxval(s2, one(0:nt-1,ibmax-ib)) - maxval(s2, .not.one(0:nt-1,ibmax-ib))
         bitmetrics(ipt+ib, nseq) = bm
      enddo
   enddo
enddo
```

Edge patching for incomplete multi-symbol groups:

```fortran
bitmetrics(205:206, 2) = bitmetrics(205:206, 1)   ! 2-sym fallback
bitmetrics(201:204, 3) = bitmetrics(201:204, 2)   ! 4-sym fallback
bitmetrics(205:206, 3) = bitmetrics(205:206, 1)   ! 4-sym fallback
```

Each mode is normalized independently via `normalizebmet()`.

#### Data bit extraction (in `ft2_decode.f90`):

Sync bit positions are **stripped out** when building the 174-bit LLR arrays:

```fortran
! Sync occupies bit positions 1-8, 67-74, 133-140, 199-206
! Data occupies bit positions 9-66, 75-132, 141-198
llra(  1: 58) = scalefac * bitmetrics(  9: 66, 1)   ! Block 1 data
llra( 59:116) = scalefac * bitmetrics( 75:132, 1)   ! Block 2 data
llra(117:174) = scalefac * bitmetrics(141:198, 1)   ! Block 3 data
```

#### Additional sync quality check on hard-decoded bit metrics:

```fortran
hbits=0
where(bitmetrics(:,1).ge.0) hbits=1
ns1=count(hbits(  1:  8).eq.(/0,0,0,1,1,0,1,1/))   ! Costas A Gray-coded
ns2=count(hbits( 67: 74).eq.(/0,1,0,0,1,1,1,0/))   ! Costas B Gray-coded
ns3=count(hbits(133:140).eq.(/1,1,1,0,0,1,0,0/))   ! Costas C Gray-coded
ns4=count(hbits(199:206).eq.(/1,0,1,1,0,0,0,1/))   ! Costas D Gray-coded
nsync_qual=ns1+ns2+ns3+ns4                          ! Reject if < 15 out of 32
```

This maps each Costas tone to its 2-bit Gray code representation:
- Costas A [0,1,3,2] → bits [00,01,11,10] → [0,0,0,1,1,0,1,1] → hbits check pattern
- etc.

#### Multi-metric combination:

```fortran
! llrd = best-of (max absolute value among a,b,c)
! llre = average of a,b,c
```

### 5. LDPC Decoding (`ft2_decode.f90`)

#### Metric passes:

5 base passes + AP passes:

| Pass | LLR source | Description |
|------|------------|-------------|
| 1 | llra | 1-symbol coherent |
| 2 | llrb | 2-symbol coherent |
| 3 | llrc | 4-symbol coherent |
| 4 | llrd | Best-of (max abs value among a,b,c) |
| 5 | llre | Average of a,b,c |

Scale factor: `scalefac = 2.83`

#### LDPC decoder call:

```fortran
Keff=91
maxosd=3                    ! OSD order (4 if near nfqso)
ndeep=3
call decode174_91(llr, Keff, maxosd, ndeep, apmask, message91, cw, ntype, nharderror, dmin)
```

Max iterations: 40

#### Descrambling after decode:

```fortran
message77 = mod(message77 + rvec, 2)   ! Remove rvec scrambling
```

Note: decodium3 uses `mod(+,2)` for descrambling (equivalent to `ieor` for binary).

#### AP (a priori) Decoding:

Additional passes 6+ use a priori information based on QSO progress state:

| nQSOProgress | Extra AP passes | AP types |
|--------------|-----------------|----------|
| 0 (CQ) | 3 | 1 (CQ), 2 (MyCall) |
| 1 (Tx1) | 3 | 2 (MyCall), 3 (MyCall+DxCall) |
| 2 (Tx2) | 3 | 2, 3 |
| 3 (Tx3) | 3 | 3 (MyCall+DxCall), 6 (RR73) |
| 4 (Tx4) | 3 | 3, 6 |
| 5 (Tx5) | 4 | 3 (MyCall+DxCall), 1 (CQ), 2 (MyCall) |

AP types inject known bit patterns (scrambled by rvec) into the LLR array at
high magnitude (`apmag = maxval(abs(llra)) * 1.1`) with corresponding apmask.

Contest-specific AP patterns supported: NA_VHF, EU_VHF, FIELD DAY, RTTY, WW_DIGI,
FOX, HOUND.

AP bits are precomputed from mycall/hiscall:
```fortran
message = trim(mycall)//' '//trim(hiscall0)//' RR73'
call pack77(message, ...)
call encode174_91(message77, cw)
apbits = 2*cw - 1          ! Convert {0,1} to {-1,+1}
```

#### Signal subtraction after successful decode (`subtractft2.f90`):

```fortran
nstart = dt*12000 + 1 - NSPS        ! NOTE: starts 1 symbol BEFORE dt
nsym = 103
NFRAME = (103+2)*NSPS = 30240
NFILT = 700                          ! Filter half-width
```

Steps:
1. Generate complex reference: `gen_ft2wave(itone, 103, NSPS, 12000.0, f0, cref, ...)`
   - Generated at **12000 Hz** (not 48000), no BT parameter
2. Compute complex amplitude: `camp(i) = dd(nstart-1+i) * conjg(cref(i))`
3. LPF via FFT: `cfilt = FFT(camp) * CW; cfilt = IFFT(cfilt)`
   - Filter: cos²(pi*j/NFILT) window, NO end correction applied
4. Subtract: `dd(j) = dd(j) - 2*REAL(cfilt(i)*cref(i))`

Key differences from a naive implementation:
- `nstart` uses `dt*12000+1-NSPS` (one full symbol before the nominal start time)
- NFILT=700 (not 600)
- No end-correction on the cos² filter edges
- Reference waveform generated at 12000 Hz sample rate

#### Subtraction loop structure:

The outer decode loop runs up to 3 passes (nsp=3 for ndepth>=2):
- Pass 1: initial decode
- Pass 2: re-search after subtracting decoded signals from pass 1
- Pass 3: re-search after subtracting decoded signals from pass 2
Each pass exits early if no new decodes were found in the previous pass.

#### DT reporting:

```fortran
xdt = ibest/1333.33 - 0.5    ! DT in seconds, offset by -0.5s
```

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

- Minimum sync power: 0.90 (used in `getcandidates2` spectral peak threshold)
- Minimum correct hard sync symbols: 4 out of 16 (in `get_ft2_bitmetrics`)
- Refined sync check: 15 out of 32 hard-decoded sync bits (in `ft2_decode.f90`)

## Encoder Details (`genft2.f90`)

### Scrambling:

```fortran
msgbits = mod(msgbits + rvec, 2)      ! XOR with scrambling vector
call encode174_91(msgbits, codeword)   ! LDPC encode the scrambled bits
```

### Gray code tone mapping:

```fortran
do i=1,ND
   is = codeword(2*i) + 2*codeword(2*i-1)   ! 2-bit index from consecutive coded bits
   if(is.le.1) itmp(i) = is                  ! 0→0, 1→1
   if(is.eq.2) itmp(i) = 3                   ! 2→3
   if(is.eq.3) itmp(i) = 2                   ! 3→2
enddo
```

Note the bit ordering: `codeword(2*i-1)` is the MSB, `codeword(2*i)` is the LSB.

### Frame assembly (1-indexed):

```fortran
i4tone(1:4)    = icos4a              ! [0,1,3,2]
i4tone(5:33)   = itmp(1:29)          ! Data block 1
i4tone(34:37)  = icos4b              ! [1,0,2,3]
i4tone(38:66)  = itmp(30:58)         ! Data block 2
i4tone(67:70)  = icos4c              ! [2,3,1,0]
i4tone(71:99)  = itmp(59:87)         ! Data block 3
i4tone(100:103)= icos4d              ! [3,2,0,1]
```

### Entry point for decoder:

```fortran
entry get_ft2_tones_from_77bits(msgbits, i4tone)
```

Takes 77-bit message (NOT scrambled), applies scrambling internally, encodes,
and returns the 103-element tone array. Used by the subtraction routine after
successful decode.
