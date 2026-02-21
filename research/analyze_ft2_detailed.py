#!/usr/bin/env python3
"""
FT2 Detailed Signal Analyzer — high-resolution spectrogram and targeted analysis.

Produces fine-grained spectrograms to visually identify tone spacing, baud rate,
and sync structure. Tests candidate baud rates against the signal.

Usage:
    uv run python analyze_ft2_detailed.py <wav_file> [--fmin 600] [--fmax 1600]
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import soundfile as sf
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from scipy import signal as sig
from scipy.fft import fft, fftfreq, ifft


COSTAS_7x7 = np.array([3, 1, 4, 0, 6, 5, 2])
FT8_SYNC_POS = [0, 36, 72]


def load_wav(path):
    data, sr = sf.read(path, dtype="float64")
    if data.ndim > 1:
        data = data.mean(axis=1)
    return data, sr


def high_res_spectrogram(data, sr, nfft, hop, fmin, fmax, title, outpath,
                          candidate_bauds=None, t_offset=0):
    """Very high resolution spectrogram with optional symbol grid overlay."""
    freqs, times, Sxx = sig.spectrogram(
        data, fs=sr, nperseg=nfft, noverlap=nfft - hop,
        window="hann", scaling="spectrum"
    )
    times = times + t_offset

    mask = (freqs >= fmin) & (freqs <= fmax)
    Sxx_band = Sxx[mask]
    freqs_band = freqs[mask]

    Sxx_db = 10 * np.log10(Sxx_band + 1e-20)
    vmin = np.percentile(Sxx_db, 20)
    vmax = np.percentile(Sxx_db, 99.5)

    fig, ax = plt.subplots(1, 1, figsize=(20, 10))
    ax.pcolormesh(times, freqs_band, Sxx_db, shading="gouraud",
                  cmap="inferno", vmin=vmin, vmax=vmax)
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Frequency (Hz)")
    ax.set_title(title)

    # Overlay candidate baud rate grids
    if candidate_bauds:
        colors = ["cyan", "lime", "magenta", "yellow"]
        for idx, (baud, label) in enumerate(candidate_bauds):
            color = colors[idx % len(colors)]
            sym_dur = 1.0 / baud
            t = times[0]
            first = True
            while t < times[-1]:
                ax.axvline(t, color=color, alpha=0.3, linewidth=0.5,
                          label=label if first else None)
                first = False
                t += sym_dur
            # Also show tone spacing
            tone_sp = baud  # assuming h=1
            for i in range(8):
                pass  # tone lines added separately if needed

    if candidate_bauds:
        ax.legend(loc="upper right", fontsize=8)

    plt.tight_layout()
    plt.savefig(outpath, dpi=200, bbox_inches="tight")
    plt.close()
    print(f"Saved: {outpath}")
    return freqs_band, times, Sxx_band


def find_signals(data, sr, fmin, fmax, nfft=2048):
    """Find individual signal centers by looking for spectral peaks."""
    freqs, times, Sxx = sig.spectrogram(
        data, fs=sr, nperseg=nfft, noverlap=nfft - nfft // 4,
        window="hann", scaling="spectrum"
    )
    mask = (freqs >= fmin) & (freqs <= fmax)
    avg_psd = Sxx[mask].mean(axis=1)
    freqs_sub = freqs[mask]

    # Smooth
    kernel = np.ones(5) / 5
    smoothed = np.convolve(avg_psd, kernel, mode="same")

    # Find peaks — expect narrow signals (~150 Hz wide)
    min_distance = int(100 / (freqs_sub[1] - freqs_sub[0]))  # at least 100 Hz apart
    peaks, props = sig.find_peaks(smoothed, height=np.median(smoothed) * 2,
                                   distance=min_distance, prominence=np.median(smoothed))

    signals = []
    for p in peaks:
        freq = freqs_sub[p]
        strength = smoothed[p]
        signals.append((freq, strength))
        print(f"  Signal at {freq:.1f} Hz (strength: {strength:.2e})")

    # Plot
    fig, ax = plt.subplots(figsize=(14, 4))
    ax.semilogy(freqs_sub, avg_psd, "b-", alpha=0.5, label="Raw")
    ax.semilogy(freqs_sub, smoothed, "r-", label="Smoothed")
    for freq, _ in signals:
        ax.axvline(freq, color="green", alpha=0.5, linestyle="--")
    ax.set_xlabel("Frequency (Hz)")
    ax.set_ylabel("Power")
    ax.set_title("Signal Detection")
    ax.legend()
    plt.tight_layout()
    plt.savefig("output/signal_detection.png", dpi=150)
    plt.close()

    return signals


def extract_single_signal(data, sr, center_freq, bw=300):
    """Bandpass filter and mix a single signal to baseband."""
    nyq = sr / 2
    low = (center_freq - bw / 2) / nyq
    high = (center_freq + bw / 2) / nyq
    low = max(low, 0.001)
    high = min(high, 0.999)
    b, a = sig.butter(5, [low, high], btype="band")
    filtered = sig.filtfilt(b, a, data)
    return filtered


def measure_baud_from_spectrogram(data, sr, center_freq, bw=300,
                                    candidate_nsps_list=None):
    """
    Measure baud rate by computing symbol-rate-matched spectrograms
    and checking how well symbols align.
    """
    print(f"\n--- Symbol Rate Detection ---")

    # Bandpass around the signal
    filtered = extract_single_signal(data, sr, center_freq, bw)

    # Analytic signal
    analytic = sig.hilbert(filtered)

    # Instantaneous frequency
    phase = np.unwrap(np.angle(analytic))
    inst_freq = np.diff(phase) * sr / (2 * np.pi)

    # Try candidate NSPS values
    if candidate_nsps_list is None:
        candidate_nsps_list = [
            (1920, "FT8 NSPS=1920 (6.25 Bd)"),
            (576, "FT4/FT2? NSPS=576 (20.83 Bd)"),
            (480, "NSPS=480 (25.0 Bd)"),
            (360, "FT2Libre NSPS=360 (33.33 Bd)"),
            (288, "NSPS=288 (41.67 Bd)"),
            (192, "NSPS=192 (62.5 Bd)"),
        ]

    results = []
    for nsps, label in candidate_nsps_list:
        baud = sr / nsps
        # Compute symbol-aligned energy variance
        # For each candidate, segment the IF into symbols and measure
        # how consistent the frequency is within each symbol vs between symbols
        n_symbols = len(inst_freq) // nsps
        if n_symbols < 10:
            continue

        within_var = 0
        between_vals = []
        for s in range(n_symbols):
            start = s * nsps
            end = start + nsps
            segment = inst_freq[start:end]
            mean_freq = np.mean(segment)
            within_var += np.var(segment)
            between_vals.append(mean_freq)

        within_var /= n_symbols
        between_var = np.var(between_vals)

        # High ratio = good symbol alignment (frequency stable within, varying between)
        ratio = between_var / (within_var + 1e-20)
        results.append((nsps, baud, label, ratio, within_var, between_var))
        print(f"  {label}: within_var={within_var:.1f}, between_var={between_var:.1f}, ratio={ratio:.3f}")

    # Sort by ratio (higher = better alignment)
    results.sort(key=lambda x: x[3], reverse=True)
    print(f"\n  Best candidate: {results[0][2]} (ratio={results[0][3]:.3f})")
    return results


def extract_tone_sequence(data, sr, center_freq, nsps, bw=300):
    """Extract the tone index for each symbol period."""
    filtered = extract_single_signal(data, sr, center_freq, bw)
    analytic = sig.hilbert(filtered)
    phase = np.unwrap(np.angle(analytic))

    baud = sr / nsps
    n_symbols = len(phase) // nsps

    tone_freqs = []
    for s in range(n_symbols):
        start = s * nsps
        end = start + nsps
        seg_phase = phase[start:end]
        # Average instantaneous frequency relative to center
        avg_freq = np.mean(np.diff(seg_phase)) * sr / (2 * np.pi)
        tone_freqs.append(avg_freq)

    tone_freqs = np.array(tone_freqs)
    return tone_freqs


def analyze_tone_spacing(tone_freqs, baud):
    """Analyze tone frequency histogram to find spacing and tone count."""
    print(f"\n--- Tone Spacing Analysis ---")

    # Remove outliers
    median = np.median(tone_freqs)
    mad = np.median(np.abs(tone_freqs - median))
    mask = np.abs(tone_freqs - median) < 5 * mad
    clean = tone_freqs[mask]

    # Histogram with fine bins
    nbins = 200
    hist, edges = np.histogram(clean, bins=nbins)
    centers = (edges[:-1] + edges[1:]) / 2

    # Smooth histogram
    kernel = np.ones(3) / 3
    smoothed = np.convolve(hist, kernel, mode="same")

    # Find peaks
    peaks, props = sig.find_peaks(smoothed, height=np.max(smoothed) * 0.05,
                                   distance=3)

    if len(peaks) < 2:
        print("  Not enough tone peaks found")
        return 0, 0, []

    peak_freqs = sorted(centers[peaks])
    spacings = np.diff(peak_freqs)

    print(f"  Tone center frequencies:")
    for i, f in enumerate(peak_freqs):
        print(f"    Tone {i}: {f:.2f} Hz")
    print(f"  Spacings: {[f'{s:.2f}' for s in spacings]}")
    print(f"  Median spacing: {np.median(spacings):.2f} Hz")
    print(f"  Mean spacing: {np.mean(spacings):.2f} Hz")
    print(f"  Number of tones: {len(peak_freqs)}")
    print(f"  Modulation index h = spacing/baud = {np.median(spacings)/baud:.3f}")

    # Plot
    fig, ax = plt.subplots(figsize=(12, 4))
    ax.bar(centers, hist, width=centers[1] - centers[0], alpha=0.5, label="Histogram")
    ax.plot(centers, smoothed, "r-", label="Smoothed")
    for f in peak_freqs:
        ax.axvline(f, color="green", linestyle="--", alpha=0.7)
    ax.set_xlabel("Frequency (Hz)")
    ax.set_ylabel("Count")
    ax.set_title(f"Tone Frequency Distribution (baud={baud:.2f})")
    ax.legend()
    plt.tight_layout()
    plt.savefig("output/tone_histogram.png", dpi=150)
    plt.close()

    return len(peak_freqs), np.median(spacings), peak_freqs


def scan_for_costas(tone_freqs, tone_spacing, n_tones, baud):
    """Scan for Costas 7x7 sync pattern in the tone sequence."""
    print(f"\n--- Costas 7x7 Scan ---")

    if tone_spacing <= 0:
        print("  No valid tone spacing")
        return

    # Quantize to tone indices
    tone_min = np.min(tone_freqs)
    tone_indices = np.round((tone_freqs - tone_min) / tone_spacing).astype(int)

    # Clip to 0-7
    tone_indices = np.clip(tone_indices, 0, 7)

    print(f"  Quantized tone range: {tone_indices.min()} to {tone_indices.max()}")
    print(f"  First 100 tones: {tone_indices[:100].tolist()}")

    # Scan all offsets for Costas match
    costas = COSTAS_7x7
    n = len(tone_indices)
    max_score = 0
    best_offset = 0
    scores = []

    for offset in range(n - 79):
        score = 0
        for sync_start in FT8_SYNC_POS:
            seg = tone_indices[offset + sync_start:offset + sync_start + 7]
            if len(seg) < 7:
                continue
            # Try all possible base tone offsets
            for base in range(max(0, seg.min() - 7), seg.max() + 1):
                matches = np.sum((seg - base) == costas)
                score = max(score, matches)
        scores.append(score)
        if score > max_score:
            max_score = score
            best_offset = offset

    scores = np.array(scores)

    print(f"  Best Costas score: {max_score}/21 at offset {best_offset}")
    print(f"  Score distribution: max={scores.max()}, mean={scores.mean():.1f}, "
          f"scores>=14: {np.sum(scores >= 14)}, scores>=17: {np.sum(scores >= 17)}")

    # Show top matches
    top_indices = np.argsort(scores)[::-1][:10]
    print(f"\n  Top 10 matches:")
    for i, idx in enumerate(top_indices):
        t_sec = idx / baud
        print(f"    #{i+1}: offset={idx} (t={t_sec:.3f}s), score={scores[idx]}/21")
        if scores[idx] >= 14:
            for j, pos in enumerate(FT8_SYNC_POS):
                seg = tone_indices[idx + pos:idx + pos + 7]
                print(f"      Costas block {j} (pos {pos}): {seg.tolist()}")

    # Plot
    fig, axes = plt.subplots(2, 1, figsize=(16, 8))

    ax = axes[0]
    ax.plot(tone_indices, ".-", markersize=1, linewidth=0.3)
    ax.set_xlabel("Symbol index")
    ax.set_ylabel("Tone index")
    ax.set_title("Quantized Tone Sequence")
    # Highlight best match region
    if max_score >= 10:
        ax.axvspan(best_offset, best_offset + 79, alpha=0.15, color="red",
                  label=f"Best 79-sym window (score={max_score})")
        for pos in FT8_SYNC_POS:
            ax.axvspan(best_offset + pos, best_offset + pos + 7,
                      alpha=0.3, color="orange")
        ax.legend()
    ax.grid(True, alpha=0.3)

    ax = axes[1]
    ax.plot(scores)
    ax.set_xlabel("Symbol offset")
    ax.set_ylabel("Costas score (/21)")
    ax.set_title("Costas 7x7 Correlation Scan")
    ax.axhline(14, color="red", linestyle="--", alpha=0.5, label="Strong (14)")
    ax.axhline(10, color="orange", linestyle="--", alpha=0.5, label="Moderate (10)")
    ax.grid(True, alpha=0.3)
    ax.legend()

    plt.tight_layout()
    plt.savefig("output/costas_scan.png", dpi=150)
    plt.close()
    print(f"  Saved: output/costas_scan.png")

    return scores, best_offset, tone_indices


def plot_fine_spectrogram_with_tones(data, sr, center_freq, nsps, tone_spacing,
                                      tone_indices, best_offset, fmin, fmax):
    """Very fine spectrogram with tone grid and Costas highlights."""
    baud = sr / nsps
    nfft = nsps * 4
    hop = nsps // 4

    freqs, times, Sxx = sig.spectrogram(
        data, fs=sr, nperseg=nfft, noverlap=nfft - hop,
        window="hann", scaling="spectrum"
    )
    mask = (freqs >= fmin) & (freqs <= fmax)

    Sxx_db = 10 * np.log10(Sxx[mask] + 1e-20)
    vmin = np.percentile(Sxx_db, 30)
    vmax = np.percentile(Sxx_db, 99.5)

    fig, ax = plt.subplots(figsize=(24, 8))
    ax.pcolormesh(times, freqs[mask], Sxx_db, shading="gouraud",
                  cmap="inferno", vmin=vmin, vmax=vmax)

    # Symbol boundaries
    sym_dur = 1.0 / baud
    t = 0
    while t < times[-1]:
        ax.axvline(t, color="white", alpha=0.15, linewidth=0.3)
        t += sym_dur

    # Tone grid
    if tone_spacing > 0:
        base = center_freq - 3.5 * tone_spacing
        for i in range(8):
            f = base + i * tone_spacing
            if fmin <= f <= fmax:
                ax.axhline(f, color="cyan", alpha=0.2, linewidth=0.5)

    # Costas regions
    if best_offset is not None:
        for pos in FT8_SYNC_POS:
            t_start = (best_offset + pos) * sym_dur
            t_end = (best_offset + pos + 7) * sym_dur
            ax.axvspan(t_start, t_end, alpha=0.1, color="cyan")

    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Frequency (Hz)")
    ax.set_title(f"Fine Spectrogram (NSPS={nsps}, baud={baud:.2f}, tone_spacing={tone_spacing:.2f} Hz)")
    plt.tight_layout()
    plt.savefig("output/fine_spectrogram.png", dpi=200, bbox_inches="tight")
    plt.close()
    print(f"  Saved: output/fine_spectrogram.png")


def main():
    parser = argparse.ArgumentParser(description="FT2 Detailed Analyzer")
    parser.add_argument("wav_file")
    parser.add_argument("--fmin", type=float, default=500)
    parser.add_argument("--fmax", type=float, default=1800)
    args = parser.parse_args()

    Path("output").mkdir(exist_ok=True)

    data, sr = load_wav(args.wav_file)
    print(f"Loaded {args.wav_file}: {sr} Hz, {len(data)/sr:.2f}s")

    # Step 1: Find individual signals
    print(f"\n=== Step 1: Signal Detection ===")
    signals = find_signals(data, sr, args.fmin, args.fmax)

    if not signals:
        print("No signals found!")
        return

    # Step 2: High-res overview spectrogram
    print(f"\n=== Step 2: Overview Spectrogram ===")
    candidate_bauds = [
        (12000 / 576, "NSPS=576 (20.83 Bd)"),
        (12000 / 360, "NSPS=360 (33.33 Bd)"),
    ]
    high_res_spectrogram(data, sr, nfft=2048, hop=64,
                          fmin=args.fmin, fmax=args.fmax,
                          title="FT2 Signal Overview",
                          outpath="output/overview_hires.png",
                          candidate_bauds=candidate_bauds)

    # Step 3: For each detected signal, try baud rate detection
    # Use the strongest signal
    strongest = max(signals, key=lambda x: x[1])
    center = strongest[0]
    print(f"\n=== Step 3: Analyzing strongest signal at {center:.1f} Hz ===")

    # Measure baud rate
    baud_results = measure_baud_from_spectrogram(data, sr, center, bw=300)

    # Step 4: Extract tone sequences for top candidates
    print(f"\n=== Step 4: Tone Analysis for Top Candidates ===")
    for nsps, baud, label, ratio, wv, bv in baud_results[:3]:
        print(f"\n  --- {label} ---")
        tone_freqs = extract_tone_sequence(data, sr, center, nsps, bw=300)

        n_tones, spacing, peak_freqs = analyze_tone_spacing(tone_freqs, baud)

        if n_tones >= 2 and spacing > 0:
            scan_result = scan_for_costas(tone_freqs, spacing, n_tones, baud)

            if scan_result:
                scores, best_offset, tone_indices = scan_result
                plot_fine_spectrogram_with_tones(
                    data, sr, center, nsps, spacing,
                    tone_indices, best_offset,
                    center - 150, center + 150
                )

    # Step 5: Summary
    print(f"\n{'='*60}")
    print(f"  ANALYSIS COMPLETE")
    print(f"  Check output/ directory for plots")
    print(f"{'='*60}")


if __name__ == "__main__":
    main()
