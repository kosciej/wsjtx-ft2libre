#!/usr/bin/env python3
"""
FT2 Tone Spacing Measurement — precise tone spacing from the NSPS=576 trace.
Uses parabolic interpolation for sub-bin peak frequency accuracy.

Usage:
    uv run python ft2_measure_tones.py ft2_signal.wav
"""

import sys
from pathlib import Path

import numpy as np
import soundfile as sf
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import signal as sig
from scipy.fft import fft


def load_wav(path):
    data, sr = sf.read(path, dtype="float64")
    if data.ndim > 1:
        data = data.mean(axis=1)
    return data, sr


def parabolic_peak(spectrum, peak_bin, df):
    """Parabolic interpolation for sub-bin frequency accuracy."""
    if peak_bin <= 0 or peak_bin >= len(spectrum) - 1:
        return peak_bin * df

    alpha = np.log(spectrum[peak_bin - 1] + 1e-20)
    beta = np.log(spectrum[peak_bin] + 1e-20)
    gamma = np.log(spectrum[peak_bin + 1] + 1e-20)

    p = 0.5 * (alpha - gamma) / (alpha - 2 * beta + gamma)
    return (peak_bin + p) * df


def extract_precise_peak_freqs(data, sr, nsps, fmin, fmax, t_start, t_end):
    """Extract peak frequency per symbol with parabolic interpolation."""
    i_start = int(t_start * sr)
    i_end = int(t_end * sr)
    segment = data[i_start:i_end]

    nfft = nsps * 8  # high oversampling for precise peaks
    df = sr / nfft
    n_symbols = len(segment) // nsps

    peak_freqs = []
    peak_powers = []

    for s in range(n_symbols):
        start = s * nsps
        end = start + nsps
        sym = segment[start:end]

        window = np.hanning(len(sym))
        spectrum = np.abs(fft(sym * window, n=nfft))[:nfft // 2]
        freqs = np.arange(nfft // 2) * df

        band = (freqs >= fmin) & (freqs <= fmax)
        if not np.any(band):
            peak_freqs.append(0)
            peak_powers.append(0)
            continue

        spec_band = spectrum[band]
        freqs_band = freqs[band]
        peak_bin_relative = np.argmax(spec_band)
        peak_bin_absolute = np.where(band)[0][peak_bin_relative]

        # Parabolic interpolation
        freq = parabolic_peak(spectrum, peak_bin_absolute, df)
        peak_freqs.append(freq)
        peak_powers.append(spec_band[peak_bin_relative])

    return np.array(peak_freqs), np.array(peak_powers)


def analyze_steps(peak_freqs, peak_powers, baud, label=""):
    """Analyze frequency step sizes to determine tone spacing."""
    print(f"\n{'='*60}")
    print(f"  TONE SPACING ANALYSIS — {label}")
    print(f"{'='*60}")

    # Only use symbols with good power (signal present)
    power_threshold = np.percentile(peak_powers, 50)
    good = peak_powers > power_threshold
    good_freqs = peak_freqs[good]
    good_indices = np.where(good)[0]

    print(f"  Total symbols: {len(peak_freqs)}")
    print(f"  Good-power symbols: {len(good_freqs)}")

    if len(good_freqs) < 5:
        print("  Not enough good symbols")
        return

    # Compute all consecutive frequency differences
    all_diffs = []
    for i in range(len(good_indices) - 1):
        if good_indices[i + 1] == good_indices[i] + 1:  # adjacent symbols
            diff = good_freqs[i + 1] - good_freqs[i]  # signed
            all_diffs.append(diff)

    all_diffs = np.array(all_diffs)
    abs_diffs = np.abs(all_diffs)

    print(f"\n  Consecutive frequency differences (signed):")
    print(f"    {[f'{d:.2f}' for d in all_diffs[:30]]}")

    # Filter out near-zero diffs (same tone) and analyze non-zero steps
    min_step = 3.0  # Hz — minimum meaningful step
    nonzero_abs = abs_diffs[abs_diffs > min_step]

    if len(nonzero_abs) == 0:
        print("  No significant frequency steps found")
        return

    print(f"\n  Non-zero step sizes (absolute):")
    print(f"    Min: {nonzero_abs.min():.2f} Hz")
    print(f"    Max: {nonzero_abs.max():.2f} Hz")
    print(f"    Mean: {nonzero_abs.mean():.2f} Hz")
    print(f"    Median: {np.median(nonzero_abs):.2f} Hz")

    # Histogram of absolute step sizes — fine bins
    fig, axes = plt.subplots(3, 1, figsize=(16, 12))

    ax = axes[0]
    ax.hist(nonzero_abs, bins=80, range=(0, 200), alpha=0.7, edgecolor="black")
    ax.set_xlabel("Step size (Hz)")
    ax.set_ylabel("Count")
    ax.set_title(f"Frequency Step Size Histogram — {label}")
    # Mark candidate tone spacings
    for h, color, lbl in [(0.5, "red", "h=0.5"), (0.75, "orange", "h=0.75"), (1.0, "green", "h=1.0")]:
        spacing = h * baud
        for mult in range(1, 8):
            ax.axvline(mult * spacing, color=color, linestyle="--", alpha=0.4,
                      label=f"{lbl} ×{mult}={mult*spacing:.1f}" if mult <= 3 else None)
    ax.legend(fontsize=7)
    ax.grid(True, alpha=0.3)

    # Peak frequencies scatter
    ax = axes[1]
    ax.plot(good_freqs, ".-", markersize=4)
    ax.set_xlabel("Symbol index (good only)")
    ax.set_ylabel("Peak frequency (Hz)")
    ax.set_title("Peak Frequency Trace (good-power symbols)")
    ax.grid(True, alpha=0.3)

    # Try to quantize to a tone grid
    ax = axes[2]
    # Test multiple tone spacings
    for h, color in [(0.5, "red"), (0.75, "orange"), (1.0, "green")]:
        spacing = h * baud
        # Find best base frequency
        test_bases = np.arange(good_freqs.min(), good_freqs.min() + spacing, 0.5)
        best_residual = 1e20
        best_base = 0
        for base in test_bases:
            quantized = np.round((good_freqs - base) / spacing) * spacing + base
            residual = np.mean((good_freqs - quantized) ** 2)
            if residual < best_residual:
                best_residual = residual
                best_base = base

        rms_error = np.sqrt(best_residual)
        quantized = np.round((good_freqs - best_base) / spacing)
        tone_range = int(quantized.max() - quantized.min())

        ax.plot(quantized - quantized.min(), ".-", markersize=4,
               label=f"h={h} (spacing={spacing:.1f}Hz, RMS err={rms_error:.2f}Hz, "
                     f"tones={tone_range+1})", color=color)

        print(f"\n  Quantization with h={h} (spacing={spacing:.2f} Hz):")
        print(f"    Best base: {best_base:.2f} Hz")
        print(f"    RMS quantization error: {rms_error:.2f} Hz")
        print(f"    Tone range: {int(quantized.min())} to {int(quantized.max())} "
              f"({tone_range + 1} tones)")
        print(f"    Tone indices: {(quantized - quantized.min()).astype(int).tolist()[:40]}")

    ax.set_xlabel("Symbol index")
    ax.set_ylabel("Tone index")
    ax.set_title("Quantized Tone Sequence (testing different h values)")
    ax.legend(fontsize=7)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(f"output/tone_spacing_{label}.png", dpi=150)
    plt.close()
    print(f"\n  Saved: output/tone_spacing_{label}.png")


def main():
    wav_path = sys.argv[1] if len(sys.argv) > 1 else "ft2_signal.wav"
    Path("output").mkdir(exist_ok=True)

    data, sr = load_wav(wav_path)
    print(f"Loaded: {sr} Hz, {len(data)/sr:.2f}s")

    nsps = 576  # Best candidate from visual analysis
    baud = sr / nsps
    print(f"Using NSPS={nsps}, baud={baud:.2f} Bd")

    # Analyze each burst separately
    bursts = [
        (0.5, 2.8, 1150, 1400, "Burst0_SigB"),
        (0.5, 2.8, 680, 920, "Burst0_SigA"),
        (4.5, 8.0, 1150, 1400, "Burst1_SigB"),
        (4.5, 8.0, 680, 920, "Burst1_SigA"),
    ]

    for t0, t1, fmin, fmax, label in bursts:
        print(f"\n  Extracting: {label} ({t0}-{t1}s, {fmin}-{fmax} Hz)")
        peak_freqs, peak_powers = extract_precise_peak_freqs(
            data, sr, nsps, fmin, fmax, t0, t1)
        analyze_steps(peak_freqs, peak_powers, baud, label)


if __name__ == "__main__":
    main()
