# FT2 Protocol — Reverse-Engineered Findings

Status: Work in progress. Based on analysis of two clean audio captures from Decodium.

## Confirmed Findings

### Modulation: 4-GFSK (NOT 8-GFSK)

The ft2.it website states 8-GFSK. Our analysis conclusively shows **4 tones only**.

Evidence: Instantaneous frequency histograms of both captures show exactly 4 sharp peaks
with no intermediate tone levels. The ~3-5% occupancy between peaks is consistent with
GFSK transition smoothing, not additional tones.

### Tone Grid

| Tone | Frequency (Hz) |
|------|----------------|
| 0    | ~1500.6        |
| 1    | ~1541.9        |
| 2    | ~1583.2        |
| 3    | ~1624.5        |

- **Tone spacing**: 41.3 Hz (extremely consistent, std dev ~0 across both captures)
- **Bits per symbol**: 2 (log2(4) = 2)
- **Bandwidth**: ~124 Hz (tone 0 to tone 3) plus GFSK rolloff

### Timing

- **Sample rate**: 12000 Hz
- **T/R cycle**: 3.8 seconds (TX occurs within each cycle)
- **TX burst duration**: ~2.16-2.32 seconds
- **3 identical bursts per T/R period** (first may be incomplete depending on timing)

## Tentative / Under Investigation

### Baud Rate and NSPS

The exact baud rate has not been pinned down. Candidates:

| NSPS | Baud (Hz) | h = spacing/baud | Notes |
|------|-----------|-------------------|-------|
| 288  | 41.67     | ~0.99 (≈1.0)     | Clean h value; 72 symbols in 2.16s |
| 270  | 44.44     | ~0.93             | Median transition interval |
| 260  | 46.15     | ~0.89             | Best cross-reference match count |

The modulation index h = tone_spacing / baud_rate. An integer or simple-fraction h
(like 1.0 or 0.5) would be typical for designed modulation schemes.

NSPS=288 with h≈1.0 is the most likely candidate on theoretical grounds.

### Symbol Count

With NSPS=288 and burst duration ~2.16s: 72000/288 = ~90 symbols possible in 2.5s,
but actual TX content appears to be ~72-96 symbols depending on NSPS.

### Sync Pattern

Cross-referencing two captures with different messages (same stations, different content)
shows matching tone positions (= sync) and differing positions (= data).

With the 4-tone model at NSPS=260, we found ~54/96 matching positions, which is much
higher than the 25% random-match rate for 4 tones. This suggests either:
- A very large sync overhead (>50% of symbols), or
- The alignment optimization is biasing toward false matches

The sync structure does NOT appear to use Costas arrays (exhaustively tested for
7-symbol Costas patterns in the 8-tone model, none found; 4-tone model has not been
tested for Costas arrays yet).

Contiguous sync blocks observed at NSPS=260: [0:12] and [20:46], but these need
independent verification.

### FEC / Error Correction

Unknown. If LDPC(174,91) is used (as in FT8/FT4):
- With 4-GFSK (2 bits/symbol): need 87 data symbols for 174 coded bits
- This leaves very few symbols for sync in a 72-symbol frame
- A different FEC code or different message format may be in use

### What the Current ft2libre Code Gets Wrong

The existing `wsjtx/lib/ft2libre/` implementation uses:
- **h = 0.75** — wrong (actual tone spacing is 41.3 Hz, not consistent with h=0.75)
- **8-GFSK** — wrong (only 4 tones observed)
- **NSPS = 360** — likely wrong (baud ~33 Hz gives h≈1.24 which doesn't match)
- **79 symbols** — likely wrong (derived from FT8 frame structure)
- **Costas sync** — no Costas arrays found in the signal

## Open Questions

1. What is the exact NSPS / baud rate?
2. What FEC code is used? Is it LDPC(174,91) or something else?
3. What is the sync pattern structure?
4. What message packing format is used? (77-bit like FT8, or different?)
5. Is there a CRC? How many bits?
6. What is the exact Gaussian BT product for the GFSK filter?

## Source Data

- `ft2_capture.wav` — 12kHz mono, ~13.47s, unknown message, captured from Decodium via BlackHole
- `ft2_capture2.wav` — 12kHz mono, ~13.44s, message "CQ CQ SP6TLS KO02", captured from Decodium via BlackHole
- `AUD-20260220-WA0029.mp3` — original off-air recording (lower quality)
