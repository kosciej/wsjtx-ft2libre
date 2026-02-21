#!/usr/bin/env python3
"""
Cross-reference two FT2 captures to identify sync vs data symbols.

Capture 1: ft2_capture.wav (13.47s, unknown message) - tone sequence already known
Capture 2: ft2_capture2.wav (13.44s, message "CQ CQ SP6TLS KO02")

Parameters: NSPS=360, h=0.5, 72 symbols, 12000 Hz sample rate
"""

import numpy as np
import soundfile as sf
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.signal import butter, filtfilt

OUTPUT_DIR = Path("output_xref")
OUTPUT_DIR.mkdir(exist_ok=True)

# Known parameters
FS = 12000
NSPS = 360
BAUD = FS / NSPS  # 33.333 Hz
H = 0.5
TONE_SPACING = H * BAUD  # 16.667 Hz
NSYMBOLS = 72
TX_DURATION = NSYMBOLS * NSPS / FS  # 2.16s

# First capture tone sequence (from previous analysis)
CAPTURE1_TONES = np.array([
    0, 0, 2, 7, 4, 0, 7, 4, 3, 7, 7, 0, 3, 7, 2, 2,
    0, 2, 7, 5, 6, 6, 5, 5, 3, 2, 3, 7, 0, 3, 3, 1,
    2, 6, 7, 0, 2, 4, 3, 5, 7, 3, 7, 7, 4, 7, 0, 0,
    2, 6, 3, 1, 1, 0, 0, 2, 0, 5, 0, 7, 0, 5, 5, 5,
    7, 5, 3, 7, 7, 6, 5, 5
])


def detect_bursts(signal, fs=FS, min_gap=1.0, min_duration=1.5):
    """Detect TX bursts using power envelope."""
    # Bandpass filter 200-2000 Hz (FT2 signal band)
    nyq = fs / 2
    b, a = butter(4, [200/nyq, 2000/nyq], btype='band')
    filtered = filtfilt(b, a, signal)

    # Power envelope with smoothing
    window_len = int(0.05 * fs)  # 50ms window
    power = np.convolve(filtered**2, np.ones(window_len)/window_len, mode='same')

    # Adaptive threshold - use absolute power threshold since silence is true zero
    peak_power = np.percentile(power, 98)
    threshold = peak_power * 0.01  # 1% of peak power

    # Find regions above threshold
    above = power > threshold
    transitions = np.diff(above.astype(int))
    starts = np.where(transitions == 1)[0]
    ends = np.where(transitions == -1)[0]

    if above[0]:
        starts = np.concatenate([[0], starts])
    if above[-1]:
        ends = np.concatenate([ends, [len(signal)-1]])

    # Ensure equal length
    n = min(len(starts), len(ends))
    starts = starts[:n]
    ends = ends[:n]

    # Filter by minimum duration and merge close bursts
    bursts = []
    for s, e in zip(starts, ends):
        duration = (e - s) / fs
        if duration >= min_duration:
            bursts.append((s, e))

    # Merge bursts separated by small gaps
    if len(bursts) > 1:
        merged = [bursts[0]]
        for s, e in bursts[1:]:
            prev_s, prev_e = merged[-1]
            gap = (s - prev_e) / fs
            if gap < min_gap:
                merged[-1] = (prev_s, e)
            else:
                merged.append((s, e))
        bursts = merged

    return bursts, power, threshold


def extract_tones(signal, start_sample, fs=FS, nsps=NSPS, nsymbols=NSYMBOLS,
                  tone_spacing=TONE_SPACING, nfft=None):
    """Extract tone sequence from a burst with fine alignment."""
    if nfft is None:
        nfft = nsps * 4  # Zero-padded FFT for better frequency resolution

    best_rms = 999
    best_tones = None
    best_offset = 0
    best_base_freq = 0

    # Fine alignment scan
    for offset in range(-nsps//2, nsps//2, nsps//20):
        start = start_sample + offset
        if start < 0 or start + nsymbols * nsps > len(signal):
            continue

        # Extract peak frequencies for each symbol
        freqs = []
        for i in range(nsymbols):
            sym_start = start + i * nsps
            sym_end = sym_start + nsps
            segment = signal[sym_start:sym_end]

            # Apply Hann window
            windowed = segment * np.hanning(len(segment))

            # Zero-padded FFT
            spectrum = np.abs(np.fft.rfft(windowed, n=nfft))

            # Find peak in signal band (200-2000 Hz)
            freq_bins = np.fft.rfftfreq(nfft, 1/fs)
            band_mask = (freq_bins >= 200) & (freq_bins <= 2000)
            band_spectrum = spectrum.copy()
            band_spectrum[~band_mask] = 0

            peak_bin = np.argmax(band_spectrum)

            # Parabolic interpolation
            if 0 < peak_bin < len(spectrum) - 1:
                alpha = np.log(spectrum[peak_bin-1] + 1e-30)
                beta = np.log(spectrum[peak_bin] + 1e-30)
                gamma = np.log(spectrum[peak_bin+1] + 1e-30)
                delta = 0.5 * (alpha - gamma) / (alpha - 2*beta + gamma + 1e-30)
                peak_freq = (peak_bin + delta) * fs / nfft
            else:
                peak_freq = peak_bin * fs / nfft

            freqs.append(peak_freq)

        freqs = np.array(freqs)

        # Try different base frequencies
        for base_try in np.arange(np.min(freqs) - tone_spacing, np.min(freqs) + tone_spacing, 0.5):
            tones_float = (freqs - base_try) / tone_spacing
            tones_int = np.round(tones_float).astype(int)
            residuals = tones_float - tones_int
            rms = np.sqrt(np.mean(residuals**2))

            # Check if all tones are in valid range (0-7)
            if np.all(tones_int >= 0) and np.all(tones_int <= 7):
                unique_tones = len(np.unique(tones_int))
                if rms < best_rms and unique_tones >= 6:
                    best_rms = rms
                    best_tones = tones_int.copy()
                    best_offset = offset
                    best_base_freq = base_try

    return best_tones, best_rms, best_offset, best_base_freq


def main():
    # Load capture 2
    print("Loading ft2_capture2.wav...")
    signal2, fs2 = sf.read("ft2_capture2.wav")
    if len(signal2.shape) > 1:
        signal2 = signal2[:, 0]
    print(f"  Duration: {len(signal2)/fs2:.2f}s, Sample rate: {fs2} Hz")

    if fs2 != FS:
        from scipy.signal import resample
        signal2 = resample(signal2, int(len(signal2) * FS / fs2))
        print(f"  Resampled to {FS} Hz")

    # Detect bursts
    print("\nDetecting TX bursts...")
    bursts, power, threshold = detect_bursts(signal2)
    print(f"  Found {len(bursts)} bursts")

    # Plot burst detection
    fig, axes = plt.subplots(2, 1, figsize=(16, 8))
    t = np.arange(len(signal2)) / FS

    axes[0].plot(t, signal2, linewidth=0.3)
    axes[0].set_title("ft2_capture2.wav - Waveform")
    axes[0].set_ylabel("Amplitude")
    for i, (s, e) in enumerate(bursts):
        axes[0].axvspan(s/FS, e/FS, alpha=0.2, color='red')
        axes[0].text((s+e)/2/FS, axes[0].get_ylim()[1]*0.9, f'B{i}',
                     ha='center', fontsize=10, fontweight='bold')

    t_pow = np.arange(len(power)) / FS
    axes[1].plot(t_pow, 10*np.log10(power + 1e-30), linewidth=0.5)
    axes[1].axhline(10*np.log10(threshold), color='r', linestyle='--', label='Threshold')
    axes[1].set_title("Power Envelope")
    axes[1].set_ylabel("Power (dB)")
    axes[1].set_xlabel("Time (s)")
    axes[1].legend()

    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / "burst_detection.png", dpi=150)
    plt.close()

    # Filter: skip first burst (likely incomplete) and any too-short bursts
    complete_bursts = []
    for i, (s, e) in enumerate(bursts):
        duration = (e - s) / FS
        print(f"  Burst {i}: {s/FS:.3f}s - {e/FS:.3f}s (duration: {duration:.3f}s)", end="")
        if i == 0:
            print(" [SKIPPED - first burst, likely incomplete]")
        elif duration < 2.0:
            print(" [SKIPPED - too short]")
        else:
            print(" [OK]")
            complete_bursts.append((i, s, e))

    print(f"\n  Using {len(complete_bursts)} complete burst(s)")

    # Extract tones from complete bursts
    print("\nExtracting tone sequences...")
    burst_tones = []
    for i, s, e in complete_bursts:
        duration = (e - s) / FS
        print(f"\n  Burst {i}: {s/FS:.3f}s - {e/FS:.3f}s (duration: {duration:.3f}s)")

        # Start from slightly after the burst start to avoid ramp-up
        margin = int(0.05 * FS)  # 50ms margin
        start_sample = s + margin

        tones, rms, offset, base_freq = extract_tones(signal2, start_sample)

        if tones is not None:
            unique = len(np.unique(tones))
            print(f"    Tones: {tones.tolist()}")
            print(f"    RMS: {rms:.4f}, Unique tones: {unique}, Base freq: {base_freq:.1f} Hz")
            print(f"    Offset: {offset} samples ({offset/FS*1000:.1f} ms)")
            burst_tones.append((i, tones, rms, offset, base_freq))
        else:
            print(f"    FAILED to extract tones")
            burst_tones.append((i, None, None, None, None))

    # Cross-reference with capture 1
    print("\n" + "="*80)
    print("CROSS-REFERENCE ANALYSIS")
    print("="*80)
    print(f"\nCapture 1 tones: {CAPTURE1_TONES.tolist()}")

    for burst_idx, tones, rms, offset, base_freq in burst_tones:
        if tones is None:
            print(f"\nBurst {burst_idx}: No tones extracted")
            continue

        print(f"\nBurst {burst_idx} tones: {tones.tolist()}")
        print(f"  RMS error: {rms:.4f}")

        # Compare with capture 1
        matches = (tones == CAPTURE1_TONES)
        match_count = np.sum(matches)
        print(f"  Matches with capture 1: {match_count}/{NSYMBOLS} ({match_count/NSYMBOLS*100:.1f}%)")

        # Show position-by-position comparison
        print(f"\n  Position-by-position comparison:")
        print(f"  {'Pos':>4} {'Cap1':>4} {'Cap2':>4} {'Match':>6}")
        print(f"  " + "-"*22)
        for j in range(NSYMBOLS):
            match_str = " SAME" if matches[j] else " DIFF"
            print(f"  {j:4d} {CAPTURE1_TONES[j]:4d} {tones[j]:4d} {match_str}")

        # Identify potential sync positions (matching symbols)
        sync_positions = np.where(matches)[0]
        data_positions = np.where(~matches)[0]

        print(f"\n  Potential SYNC positions ({len(sync_positions)}): {sync_positions.tolist()}")
        print(f"  Potential DATA positions ({len(data_positions)}): {data_positions.tolist()}")

        # Check if sync positions form a recognizable pattern
        if len(sync_positions) >= 7:
            # Check for 7-symbol blocks
            print(f"\n  Looking for 7-symbol sync blocks...")
            for start in range(NSYMBOLS - 6):
                if all(matches[start:start+7]):
                    block_tones = tones[start:start+7]
                    # Normalize to start from 0
                    block_norm = block_tones - np.min(block_tones)
                    print(f"    Block at [{start}:{start+7}]: "
                          f"raw={block_tones.tolist()}, "
                          f"normalized={block_norm.tolist()}")

            # Check for scattered sync pattern (e.g., every Nth symbol)
            print(f"\n  Looking for periodic sync patterns...")
            for period in range(2, NSYMBOLS//7 + 1):
                for phase in range(period):
                    positions = list(range(phase, NSYMBOLS, period))
                    if len(positions) >= 7:
                        matches_at_pos = sum(1 for p in positions if matches[p])
                        if matches_at_pos >= len(positions) * 0.8:
                            print(f"    Period={period}, Phase={phase}: "
                                  f"{matches_at_pos}/{len(positions)} match")

    # Cross-reference between all capture 2 bursts
    if len(burst_tones) >= 2:
        print("\n" + "="*80)
        print("INTER-BURST COMPARISON (within capture 2)")
        print("="*80)
        valid_bursts = [(idx, t) for idx, t, r, o, b in burst_tones if t is not None]

        for i in range(len(valid_bursts)):
            for j in range(i+1, len(valid_bursts)):
                idx_i, tones_i = valid_bursts[i]
                idx_j, tones_j = valid_bursts[j]
                matches_ij = (tones_i == tones_j)
                match_count = np.sum(matches_ij)
                print(f"\n  Burst {idx_i} vs Burst {idx_j}: "
                      f"{match_count}/{NSYMBOLS} match ({match_count/NSYMBOLS*100:.1f}%)")

                if match_count < NSYMBOLS:  # Not identical
                    diff_positions = np.where(~matches_ij)[0]
                    print(f"    Different at positions: {diff_positions.tolist()}")

        # If all capture 2 bursts are identical (same message repeated),
        # find consensus tone sequence
        if len(valid_bursts) >= 2:
            all_tones = np.array([t for _, t in valid_bursts])
            # Use mode (most common value) at each position
            from scipy.stats import mode
            consensus = mode(all_tones, axis=0, keepdims=False)
            consensus_tones = consensus.mode
            consensus_count = consensus.count

            print(f"\n  Consensus tone sequence (mode of all bursts):")
            print(f"    {consensus_tones.tolist()}")
            print(f"    Agreement counts: {consensus_count.tolist()}")

            # Final cross-reference: consensus of capture 2 vs capture 1
            print(f"\n" + "="*80)
            print("FINAL CROSS-REFERENCE: Capture 2 Consensus vs Capture 1")
            print("="*80)

            matches_final = (consensus_tones == CAPTURE1_TONES)
            match_count = np.sum(matches_final)
            print(f"  Matches: {match_count}/{NSYMBOLS} ({match_count/NSYMBOLS*100:.1f}%)")

            sync_pos = np.where(matches_final)[0]
            data_pos = np.where(~matches_final)[0]

            print(f"\n  SYNC positions (same in both captures): {sync_pos.tolist()}")
            print(f"  DATA positions (different): {data_pos.tolist()}")

            # Extract sync tones at sync positions
            sync_tones = consensus_tones[sync_pos]
            print(f"  Sync tone values: {sync_tones.tolist()}")

            # Check for Costas-like patterns in sync
            if len(sync_pos) >= 7:
                print(f"\n  Checking sync for 7-symbol contiguous blocks:")
                for k in range(len(sync_pos) - 6):
                    block_pos = sync_pos[k:k+7]
                    if np.all(np.diff(block_pos) == 1):  # Contiguous
                        block_vals = consensus_tones[block_pos]
                        block_norm = block_vals - np.min(block_vals)
                        print(f"    Block at [{block_pos[0]}:{block_pos[-1]+1}]: "
                              f"values={block_vals.tolist()}, "
                              f"normalized={block_norm.tolist()}")

                        # Check if it's a Costas array
                        if len(set(block_norm)) == 7 and np.max(block_norm) == 6:
                            print(f"      -> This IS a valid permutation of 0-6!")
                            # Check Costas property
                            diffs = []
                            is_costas = True
                            for a in range(7):
                                for b in range(a+1, 7):
                                    d = (block_norm[b] - block_norm[a], b - a)
                                    if d in diffs:
                                        is_costas = False
                                    diffs.append(d)
                            print(f"      -> Costas property: {is_costas}")

            # Visualization
            fig, ax = plt.subplots(figsize=(20, 6))
            x = np.arange(NSYMBOLS)
            width = 0.35

            bars1 = ax.bar(x - width/2, CAPTURE1_TONES, width, label='Capture 1 (unknown msg)',
                          color=['green' if m else 'lightcoral' for m in matches_final], alpha=0.7)
            bars2 = ax.bar(x + width/2, consensus_tones, width, label='Capture 2 (CQ CQ CQ CQ CQ)',
                          color=['green' if m else 'salmon' for m in matches_final], alpha=0.7)

            # Mark sync positions
            for pos in sync_pos:
                ax.axvline(pos, color='green', alpha=0.15, linewidth=3)

            ax.set_xlabel("Symbol Position")
            ax.set_ylabel("Tone Value (0-7)")
            ax.set_title(f"Cross-Reference: Capture 1 vs Capture 2\n"
                        f"Green = SYNC (same), Red = DATA (different)")
            ax.set_xticks(range(0, NSYMBOLS, 5))
            ax.legend()
            ax.set_ylim(-0.5, 8)
            ax.grid(axis='y', alpha=0.3)

            plt.tight_layout()
            plt.savefig(OUTPUT_DIR / "cross_reference.png", dpi=150)
            plt.close()
            print(f"\n  Saved cross_reference.png")

            # Summary figure showing sync structure
            fig, ax = plt.subplots(figsize=(20, 3))
            colors = ['green' if m else 'red' for m in matches_final]
            ax.bar(range(NSYMBOLS), [1]*NSYMBOLS, color=colors, edgecolor='black', linewidth=0.5)
            ax.set_xlabel("Symbol Position")
            ax.set_title("Frame Structure: Green=SYNC, Red=DATA")
            ax.set_yticks([])
            ax.set_xticks(range(0, NSYMBOLS, 1))
            ax.set_xticklabels([str(i) for i in range(NSYMBOLS)], fontsize=5, rotation=90)
            plt.tight_layout()
            plt.savefig(OUTPUT_DIR / "frame_structure.png", dpi=150)
            plt.close()
            print(f"  Saved frame_structure.png")

    print("\nDone!")


if __name__ == "__main__":
    main()
