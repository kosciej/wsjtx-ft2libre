#!/usr/bin/env python3
"""
Robust cross-reference of two FT2 captures to identify sync vs data symbols.

Step 1: Examine both captures to find signal frequency bands
Step 2: Extract tones with parameters adapted to each capture
Step 3: Cross-reference to identify sync positions

Parameters: NSPS=360, h=0.5, 72 symbols, 12000 Hz sample rate
"""

import numpy as np
import soundfile as sf
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.signal import butter, filtfilt

OUTPUT_DIR = Path("output_xref")
OUTPUT_DIR.mkdir(exist_ok=True)

FS = 12000
NSPS = 360
BAUD = FS / NSPS  # 33.333 Hz
H = 0.5
TONE_SPACING = H * BAUD  # 16.667 Hz
NSYMBOLS = 72
NFFT = NSPS * 4  # 1440-point FFT

# First capture tone sequence (from previous analysis, middle burst)
CAPTURE1_TONES = np.array([
    0, 0, 2, 7, 4, 0, 7, 4, 3, 7, 7, 0, 3, 7, 2, 2,
    0, 2, 7, 5, 6, 6, 5, 5, 3, 2, 3, 7, 0, 3, 3, 1,
    2, 6, 7, 0, 2, 4, 3, 5, 7, 3, 7, 7, 4, 7, 0, 0,
    2, 6, 3, 1, 1, 0, 0, 2, 0, 5, 0, 7, 0, 5, 5, 5,
    7, 5, 3, 7, 7, 6, 5, 5
])


def find_burst_region(signal, fs=FS):
    """Find the middle complete burst by power envelope."""
    window_len = int(0.05 * fs)
    power = np.convolve(signal**2, np.ones(window_len)/window_len, mode='same')
    threshold = np.max(power) * 0.01

    above = power > threshold
    transitions = np.diff(above.astype(int))
    starts = np.where(transitions == 1)[0]
    ends = np.where(transitions == -1)[0]

    if above[0]:
        starts = np.concatenate([[0], starts])
    if above[-1]:
        ends = np.concatenate([ends, [len(signal)-1]])

    n = min(len(starts), len(ends))
    bursts = []
    for i in range(n):
        duration = (ends[i] - starts[i]) / fs
        if duration > 1.5:
            bursts.append((starts[i], ends[i], duration))

    print(f"  Found {len(bursts)} bursts:")
    for i, (s, e, d) in enumerate(bursts):
        print(f"    Burst {i}: {s/fs:.3f}s - {e/fs:.3f}s ({d:.3f}s)")

    # Skip first burst, use first complete one after that
    if len(bursts) >= 2:
        print(f"  -> Using burst 1 (skipping first)")
        return bursts[1][:2]
    elif len(bursts) == 1:
        print(f"  -> Only one burst found, using it")
        return bursts[0][:2]
    else:
        return None


def extract_peak_freqs(signal, start_sample, nsymbols=NSYMBOLS, nsps=NSPS,
                       fs=FS, nfft=NFFT, freq_lo=200, freq_hi=2500):
    """Extract peak frequency for each symbol."""
    freqs = []
    snrs = []
    freq_bins = np.fft.rfftfreq(nfft, 1/fs)
    band_mask = (freq_bins >= freq_lo) & (freq_bins <= freq_hi)

    for i in range(nsymbols):
        sym_start = start_sample + i * nsps
        sym_end = sym_start + nsps
        if sym_end > len(signal):
            freqs.append(np.nan)
            snrs.append(0)
            continue

        segment = signal[sym_start:sym_end]
        windowed = segment * np.hanning(len(segment))
        spectrum = np.abs(np.fft.rfft(windowed, n=nfft))

        band_spectrum = spectrum.copy()
        band_spectrum[~band_mask] = 0

        peak_bin = np.argmax(band_spectrum)

        # Parabolic interpolation
        if 0 < peak_bin < len(spectrum) - 1:
            alpha = np.log(spectrum[peak_bin-1] + 1e-30)
            beta = np.log(spectrum[peak_bin] + 1e-30)
            gamma = np.log(spectrum[peak_bin+1] + 1e-30)
            denom = alpha - 2*beta + gamma
            if abs(denom) > 1e-10:
                delta = 0.5 * (alpha - gamma) / denom
            else:
                delta = 0
            peak_freq = (peak_bin + delta) * fs / nfft
        else:
            peak_freq = peak_bin * fs / nfft

        # SNR: peak power vs median
        peak_power = spectrum[peak_bin]**2
        noise = np.median(band_spectrum[band_mask])**2 + 1e-30
        snr = peak_power / noise

        freqs.append(peak_freq)
        snrs.append(snr)

    return np.array(freqs), np.array(snrs)


def quantize_tones(freqs, tone_spacing=TONE_SPACING):
    """Find best base frequency and quantize to tone indices 0-7."""
    valid = ~np.isnan(freqs)
    f = freqs[valid]

    best_rms = 999
    best_tones = None
    best_base = 0

    # Scan base frequencies
    f_min = np.min(f)
    f_max = np.max(f)
    for base in np.arange(f_min - 2*tone_spacing, f_min + tone_spacing, 0.25):
        tones_float = (f - base) / tone_spacing
        tones_int = np.round(tones_float).astype(int)
        residuals = tones_float - tones_int

        if np.any(tones_int < 0) or np.any(tones_int > 7):
            continue

        rms = np.sqrt(np.mean(residuals**2))
        unique = len(np.unique(tones_int))

        if rms < best_rms and unique >= 5:
            best_rms = rms
            best_base = base
            best_tones_valid = tones_int.copy()

    if best_tones is None and best_rms < 999:
        # Map back to full array
        best_tones = np.full(len(freqs), -1, dtype=int)
        best_tones[valid] = best_tones_valid

    return best_tones, best_rms, best_base


def main():
    # ================================================================
    # STEP 1: Load both captures
    # ================================================================
    print("="*80)
    print("LOADING CAPTURES")
    print("="*80)

    sig1, _ = sf.read("ft2_capture.wav")
    if len(sig1.shape) > 1:
        sig1 = sig1[:, 0]
    print(f"\nCapture 1: {len(sig1)/FS:.2f}s, max={np.max(np.abs(sig1)):.4f}")

    sig2, _ = sf.read("ft2_capture2.wav")
    if len(sig2.shape) > 1:
        sig2 = sig2[:, 0]
    print(f"Capture 2: {len(sig2)/FS:.2f}s, max={np.max(np.abs(sig2)):.4f}")

    # ================================================================
    # STEP 2: Find middle burst in each capture
    # ================================================================
    print("\n" + "="*80)
    print("BURST DETECTION")
    print("="*80)

    print("\nCapture 1:")
    burst1 = find_burst_region(sig1)
    print("\nCapture 2:")
    burst2 = find_burst_region(sig2)

    if burst1 is None or burst2 is None:
        print("ERROR: Could not find bursts in one or both captures")
        return

    # ================================================================
    # STEP 3: Extract frequencies with fine alignment scan
    # ================================================================
    print("\n" + "="*80)
    print("FREQUENCY EXTRACTION WITH FINE ALIGNMENT")
    print("="*80)

    for label, signal, (burst_start, burst_end) in [
        ("Capture 1", sig1, burst1),
        ("Capture 2", sig2, burst2)
    ]:
        margin = int(0.05 * FS)
        nominal_start = burst_start + margin

        print(f"\n{label}: burst at {burst_start/FS:.3f}s - {burst_end/FS:.3f}s")

        # Quick frequency overview at nominal alignment
        freqs, snrs = extract_peak_freqs(signal, nominal_start)
        print(f"  Freq range: {np.nanmin(freqs):.1f} - {np.nanmax(freqs):.1f} Hz")
        print(f"  Freq span: {np.nanmax(freqs) - np.nanmin(freqs):.1f} Hz")
        print(f"  Expected span for 8 tones: {7 * TONE_SPACING:.1f} Hz")
        print(f"  Mean SNR: {np.nanmean(snrs):.1f}")

    # ================================================================
    # STEP 4: Fine alignment scan for both captures
    # ================================================================
    print("\n" + "="*80)
    print("FINE ALIGNMENT SCAN")
    print("="*80)

    results = {}
    for label, signal, (burst_start, burst_end) in [
        ("Capture 1", sig1, burst1),
        ("Capture 2", sig2, burst2)
    ]:
        margin = int(0.05 * FS)
        nominal_start = burst_start + margin

        best_rms = 999
        best_result = None

        # Scan offsets
        step = NSPS // 36  # 10-sample steps
        for offset in range(-NSPS, NSPS, step):
            start = nominal_start + offset
            if start < 0:
                continue
            if start + NSYMBOLS * NSPS > len(signal):
                continue

            freqs, snrs = extract_peak_freqs(signal, start)
            tones, rms, base = quantize_tones(freqs)

            if tones is not None and rms < best_rms:
                best_rms = rms
                best_result = {
                    'tones': tones,
                    'rms': rms,
                    'base': base,
                    'offset': offset,
                    'freqs': freqs,
                    'snrs': snrs,
                    'start': start
                }

        if best_result:
            t = best_result['tones']
            unique = len(np.unique(t[t >= 0]))
            print(f"\n{label}:")
            print(f"  Best RMS: {best_result['rms']:.4f}")
            print(f"  Base freq: {best_result['base']:.1f} Hz")
            print(f"  Offset: {best_result['offset']} samples ({best_result['offset']/FS*1000:.1f}ms)")
            print(f"  Unique tones: {unique}")
            print(f"  Tones: {best_result['tones'].tolist()}")
            results[label] = best_result
        else:
            print(f"\n{label}: FAILED to extract tones!")

    if len(results) < 2:
        print("\nCannot cross-reference without both captures.")
        return

    # ================================================================
    # STEP 5: CROSS-REFERENCE
    # ================================================================
    print("\n" + "="*80)
    print("CROSS-REFERENCE ANALYSIS")
    print("="*80)

    t1 = results["Capture 1"]["tones"]
    t2 = results["Capture 2"]["tones"]

    # Also compare with the pre-computed capture 1 tones
    print(f"\nCapture 1 (fresh extraction): {t1.tolist()}")
    print(f"Capture 1 (previous result):  {CAPTURE1_TONES.tolist()}")
    prev_match = np.sum(t1 == CAPTURE1_TONES)
    print(f"  Consistency check: {prev_match}/{NSYMBOLS} match with previous extraction")

    print(f"\nCapture 2 (fresh extraction): {t2.tolist()}")

    matches = (t1 == t2)
    match_count = np.sum(matches)
    print(f"\nMatches: {match_count}/{NSYMBOLS} ({match_count/NSYMBOLS*100:.1f}%)")

    sync_pos = np.where(matches)[0]
    data_pos = np.where(~matches)[0]

    print(f"\nSYNC positions ({len(sync_pos)} symbols — same in both captures):")
    print(f"  {sync_pos.tolist()}")
    print(f"  Tone values at sync: {t1[sync_pos].tolist()}")

    print(f"\nDATA positions ({len(data_pos)} symbols — different between captures):")
    print(f"  {data_pos.tolist()}")

    # ================================================================
    # STEP 6: Look for sync structure
    # ================================================================
    print("\n" + "="*80)
    print("SYNC STRUCTURE ANALYSIS")
    print("="*80)

    # Mark each position
    frame_map = ['S' if matches[i] else 'D' for i in range(NSYMBOLS)]
    print(f"\nFrame map (S=sync, D=data):")
    # Print in rows of 36
    for row_start in range(0, NSYMBOLS, 36):
        row_end = min(row_start + 36, NSYMBOLS)
        positions = ''.join(f'{i:3d}' for i in range(row_start, row_end))
        symbols = ''.join(f'  {frame_map[i]}' for i in range(row_start, row_end))
        t1_vals = ''.join(f'{t1[i]:3d}' for i in range(row_start, row_end))
        t2_vals = ''.join(f'{t2[i]:3d}' for i in range(row_start, row_end))
        print(f"  Pos:  {positions}")
        print(f"  Map:  {symbols}")
        print(f"  Cap1: {t1_vals}")
        print(f"  Cap2: {t2_vals}")
        print()

    # Look for contiguous sync blocks
    print("Contiguous sync blocks:")
    in_block = False
    block_start = 0
    for i in range(NSYMBOLS + 1):
        if i < NSYMBOLS and matches[i]:
            if not in_block:
                block_start = i
                in_block = True
        else:
            if in_block:
                block_len = i - block_start
                block_tones = t1[block_start:i]
                norm = block_tones - np.min(block_tones)
                print(f"  [{block_start}:{i}] length={block_len} "
                      f"tones={block_tones.tolist()} normalized={norm.tolist()}")
                in_block = False

    # Check: are sync positions at known patterns like FT8's [0:7, 36:43, 72:79]?
    print(f"\nChecking known sync position patterns:")
    patterns = {
        "FT8-like [0:7, 36:43, 65:72]": list(range(0,7)) + list(range(36,43)) + list(range(65,72)),
        "FT8-like [0:7, 33:40, 65:72]": list(range(0,7)) + list(range(33,40)) + list(range(65,72)),
        "Start+End [0:7, 65:72]": list(range(0,7)) + list(range(65,72)),
        "Start+Mid [0:7, 36:43]": list(range(0,7)) + list(range(36,43)),
        "FT4-like [0:4, 17:21, 34:38, 51:55, 68:72]": (
            list(range(0,4)) + list(range(17,21)) + list(range(34,38)) +
            list(range(51,55)) + list(range(68,72))),
        "Evenly3 [0:7, 22:29, 44:51]": list(range(0,7)) + list(range(22,29)) + list(range(44,51)),
    }

    for name, positions in patterns.items():
        positions = [p for p in positions if p < NSYMBOLS]
        sync_at = sum(1 for p in positions if matches[p])
        data_at = sum(1 for p in range(NSYMBOLS) if p not in positions and not matches[p])
        total_data = NSYMBOLS - len(positions)
        print(f"  {name}: sync={sync_at}/{len(positions)}, "
              f"data_correct={data_at}/{total_data}")

    # Count data symbols for each pattern
    print(f"\nData capacity check (need 58 data symbols for 174 bits):")
    for name, positions in patterns.items():
        positions = [p for p in positions if p < NSYMBOLS]
        n_data = NSYMBOLS - len(positions)
        n_bits = n_data * 3
        print(f"  {name}: {len(positions)} sync + {n_data} data = {n_bits} bits "
              f"({'OK' if n_bits >= 174 else 'TOO FEW'})")

    # ================================================================
    # STEP 7: Visualization
    # ================================================================
    fig, axes = plt.subplots(3, 1, figsize=(20, 12))

    # Tone comparison
    x = np.arange(NSYMBOLS)
    width = 0.35
    colors1 = ['green' if matches[i] else 'red' for i in range(NSYMBOLS)]
    colors2 = ['green' if matches[i] else 'blue' for i in range(NSYMBOLS)]

    axes[0].bar(x - width/2, t1, width, color=colors1, alpha=0.7, label='Capture 1')
    axes[0].bar(x + width/2, t2, width, color=colors2, alpha=0.7, label='Capture 2')
    axes[0].set_ylabel("Tone (0-7)")
    axes[0].set_title("Cross-Reference: Capture 1 vs Capture 2\nGreen=SYNC (same), Red/Blue=DATA (different)")
    axes[0].legend()
    axes[0].set_ylim(-0.5, 8)

    # Frame structure map
    frame_colors = ['green' if matches[i] else 'red' for i in range(NSYMBOLS)]
    axes[1].bar(range(NSYMBOLS), [1]*NSYMBOLS, color=frame_colors, edgecolor='black', linewidth=0.5)
    axes[1].set_title("Frame Structure: Green=SYNC, Red=DATA")
    axes[1].set_yticks([])
    for i in range(NSYMBOLS):
        if matches[i]:
            axes[1].text(i, 0.5, str(t1[i]), ha='center', va='center', fontsize=6, fontweight='bold')

    # Frequency comparison
    f1 = results["Capture 1"]["freqs"]
    f2 = results["Capture 2"]["freqs"]
    axes[2].plot(x, f1, 'ro-', markersize=3, linewidth=0.5, label='Capture 1')
    axes[2].plot(x, f2, 'bs-', markersize=3, linewidth=0.5, label='Capture 2')
    axes[2].set_ylabel("Frequency (Hz)")
    axes[2].set_xlabel("Symbol Position")
    axes[2].set_title("Raw Peak Frequencies")
    axes[2].legend()
    axes[2].grid(True, alpha=0.3)

    for ax in axes:
        ax.set_xlim(-0.5, NSYMBOLS - 0.5)
        ax.set_xticks(range(0, NSYMBOLS, 5))

    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / "cross_reference.png", dpi=150)
    plt.close()
    print(f"\nSaved cross_reference.png")

    # Symbol-aligned waterfall for capture 2
    print("\nGenerating symbol-aligned waterfall for capture 2...")
    start2 = results["Capture 2"]["start"]
    nfft_wf = NSPS * 2
    waterfall = np.zeros((NSYMBOLS, nfft_wf // 2 + 1))
    freq_axis = np.fft.rfftfreq(nfft_wf, 1/FS)

    for i in range(NSYMBOLS):
        sym_start = start2 + i * NSPS
        sym_end = sym_start + NSPS
        if sym_end <= len(sig2):
            segment = sig2[sym_start:sym_end]
            windowed = segment * np.hanning(len(segment))
            spectrum = np.abs(np.fft.rfft(windowed, n=nfft_wf))
            waterfall[i] = spectrum

    fig, ax = plt.subplots(figsize=(16, 10))
    f_lo, f_hi = 500, 2000
    f_mask = (freq_axis >= f_lo) & (freq_axis <= f_hi)
    extent = [freq_axis[f_mask][0], freq_axis[f_mask][-1], NSYMBOLS-0.5, -0.5]
    ax.imshow(waterfall[:, f_mask], aspect='auto', extent=extent, cmap='viridis',
              interpolation='nearest')

    # Mark sync positions
    for pos in sync_pos:
        ax.plot([freq_axis[f_mask][0], freq_axis[f_mask][-1]], [pos, pos],
                'r-', alpha=0.3, linewidth=0.5)

    ax.set_xlabel("Frequency (Hz)")
    ax.set_ylabel("Symbol index")
    ax.set_title("Capture 2 - Symbol-aligned waterfall\n(red lines = sync positions)")
    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / "capture2_waterfall.png", dpi=150)
    plt.close()
    print("Saved capture2_waterfall.png")

    print("\n" + "="*80)
    print("SUMMARY")
    print("="*80)
    print(f"SYNC symbols: {len(sync_pos)}/{NSYMBOLS}")
    print(f"DATA symbols: {len(data_pos)}/{NSYMBOLS}")
    if len(data_pos) > 0:
        n_data_bits = len(data_pos) * 3
        print(f"Data capacity: {len(data_pos)} × 3 = {n_data_bits} bits")
        print(f"LDPC(174,91) needs: 58 data symbols = 174 bits")
        if n_data_bits == 174:
            print("  -> PERFECT MATCH!")
        elif n_data_bits > 174:
            print(f"  -> {n_data_bits - 174} excess bits (some 'sync' positions may be data)")
        else:
            print(f"  -> {174 - n_data_bits} bits short (some 'data' positions may actually be sync)")


if __name__ == "__main__":
    main()
