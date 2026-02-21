# FT2 Research Methodology

How we reverse-engineered the FT2 protocol parameters from audio captures.

## Signal Acquisition

1. Installed Decodium (macOS) — an FT2-capable SDR application
2. Used BlackHole virtual audio driver to route Decodium's audio output directly
   to a recording application, producing clean 12kHz mono WAV files with no
   analog noise or acoustic path degradation
3. Captured two transmissions with different messages to enable cross-referencing:
   - `ft2_capture.wav` — unknown message content
   - `ft2_capture2.wav` — "CQ CQ SP6TLS KO02"

## Phase 1: Initial Analysis (8-GFSK Assumption)

Based on the ft2.it website claiming "8-GFSK", we initially assumed 8 tones.

### Scripts and what they found

- `analyze_ft2_signal.py`, `analyze_ft2_detailed.py`, `analyze_ft2_zoom.py`,
  `analyze_ft2_ultrazoom.py` — Spectrogram analysis at various resolutions.
  Confirmed signal is FSK with ~41 Hz tone spacing.
- `ft2_measure_tones.py` — Attempted to measure 8 tone levels. Found tones
  but couldn't reliably identify 8 distinct levels.
- `ft2_sync_detect.py` — Looked for repeating sync patterns.
- `ft2_find_sync.py` — Systematic sync pattern search assuming 8-tone Costas arrays.
- `ft2_costas_all.py` — Tested all 200 order-7 Costas arrays. None matched.
- `ft2_frame_detect.py` — TX boundary detection, measured burst duration (~2.16s),
  confirmed 3 bursts per T/R period.
- `ft2_middle_burst.py` — Focused on cleanest (middle) burst.
- `ft2_test_nsps288.py` — Systematic NSPS scan from 240-400 at multiple h values.
  Found NSPS=360 gave best quantization RMS.
- `ft2_nsps360_focused.py` — Deep dive at NSPS=360. Still no Costas sync found.
- `ft2_verify_costas.py` — Final exhaustive Costas verification. Negative.

### Cross-reference approach (8-tone)

After obtaining the second capture, we compared tone sequences between captures.
Positions where both captures produce the same tone = sync (fixed pattern).
Positions where tones differ = data (message-dependent).

- `ft2_cross_reference.py` — Initial cross-reference. Failed because capture 2
  was initially silent (re-recorded).
- `ft2_xref_robust.py` — Robust extraction. Failed for capture 2: frequency span
  was 148.4 Hz, too wide for 8 tones at h=0.5 (expected 116.7 Hz). This was the
  first hint something was wrong with the 8-tone assumption.
- `ft2_xref_v2.py` — Found 2 persistent outlier symbols (positions 17, 32) that
  couldn't be quantized to 0-7 regardless of h value.
- `ft2_xref_v3.py` — Outlier-tolerant version. Got scattered sync pattern (16/72).
- `ft2_xref_joint.py` — Joint alignment optimization. Got 30/72 matches but no
  clear block structure.
- `ft2_xref_fine.py` — Fine alignment scan at NSPS=359 and 360. Max 24 matches.

All 8-tone cross-reference results were noisy and inconsistent.

## Phase 2: The Breakthrough — Fresh Analysis Without Assumptions

**Key insight**: The user pointed out we were assuming FT2 is a derivative of FT8.
We dropped all assumptions and analyzed the raw signal.

### The critical script: `ft2_raw_analysis.py`

This script applies no prior assumptions. It uses:

1. **Hilbert transform** to compute the analytic signal
2. **Phase unwrapping** of the analytic signal phase
3. **Numerical differentiation** of the unwrapped phase to get instantaneous frequency
4. **Histogram of instantaneous frequency** during the TX burst (excluding transitions)

The instantaneous frequency histogram showed **exactly 4 sharp peaks**, not 8:

```
Tone 0: ~1500.6 Hz
Tone 1: ~1541.9 Hz
Tone 2: ~1583.2 Hz
Tone 3: ~1624.5 Hz
Spacing: ~41.3 Hz (consistent across both captures)
```

We verified there were no intermediate tones by:
- Widening the bandpass to 1100-1900 Hz — still only 4 peaks
- Checking the histogram with finer bins — the ~3-5% occupancy between peaks is
  exactly what GFSK transition smoothing produces (continuous frequency sweep
  between tone centers during symbol transitions)

This was the decisive finding: **FT2 uses 4-GFSK, not 8-GFSK**.

### Baud rate measurement attempts

Multiple methods were tried but gave ambiguous results:

- **Autocorrelation of instantaneous frequency**: Peak too broad, multiple candidates
- **|df/dt| spectrum**: Noisy, multiple peaks
- **Tone transition interval detection**: Measured intervals between frequency
  threshold crossings. Median ~261 samples, mean ~274-300. Most consistent with
  NSPS around 270-290 but not conclusive.
- **Quantization RMS scan**: With 41 Hz spacing and 4 tones, quantization is
  trivially easy at any NSPS, so RMS doesn't discriminate well.

## Phase 3: 4-Tone Cross-Reference

- `ft2_xref_4tone.py` — Cross-reference with 4-tone model, scanned NSPS 260-398.
  Best at NSPS=260: 54/96 matches (56%). Shows contiguous sync blocks but the
  high match rate (vs 25% random for 4 tones) suggests possible alignment bias.

## Key Lessons

1. **Don't trust the spec**: The ft2.it website says "8-GFSK" but the signal clearly
   uses only 4 tones. The website may be describing a planned/future version, or the
   specification may simply be wrong.

2. **Instantaneous frequency histograms are definitive**: While spectrograms and FFT-based
   tone detection can be ambiguous (especially with GFSK smoothing), the instantaneous
   frequency histogram unambiguously shows the number and spacing of tones.

3. **Cross-referencing captures with different messages is powerful**: It cleanly separates
   sync (constant) from data (variable) without needing to know the sync pattern a priori.
   But the alignment optimization can introduce bias — the best approach would be to
   fix alignment by an independent method (e.g., burst onset detection) rather than
   optimizing for maximum match count.

## Tools Used

- Python 3.x with NumPy, SciPy, Matplotlib, soundfile
- All scripts in `research/` directory
- Output plots in `output_raw/` (key results) and `output_xref/` (cross-reference)
