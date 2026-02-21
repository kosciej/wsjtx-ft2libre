#!/usr/bin/env python3
"""
FT2 Zoomed Analysis — high-resolution spectrograms of individual bursts
to visually measure symbol duration and tone spacing.

Usage:
    uv run python analyze_ft2_zoom.py ft2_signal.wav
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


def find_bursts(data, sr, min_gap_s=0.5):
    """Find signal bursts by energy envelope."""
    # Short-time energy
    frame = int(sr * 0.01)  # 10ms frames
    n_frames = len(data) // frame
    energy = np.array([
        np.sum(data[i * frame:(i + 1) * frame] ** 2) / frame
        for i in range(n_frames)
    ])

    # Smooth
    kernel = np.ones(10) / 10
    smooth_energy = np.convolve(energy, kernel, mode="same")

    threshold = np.percentile(smooth_energy, 70)
    active = smooth_energy > threshold

    # Find burst boundaries
    transitions = np.diff(active.astype(int))
    starts = np.where(transitions == 1)[0]
    ends = np.where(transitions == -1)[0]

    if active[0]:
        starts = np.concatenate([[0], starts])
    if active[-1]:
        ends = np.concatenate([ends, [len(active) - 1]])

    # Merge close bursts
    min_gap = int(min_gap_s / 0.01)
    bursts = []
    for s, e in zip(starts, ends):
        t_start = s * frame / sr
        t_end = e * frame / sr
        duration = t_end - t_start
        if duration < 0.5:  # skip very short
            continue
        if bursts and (t_start - bursts[-1][1]) < min_gap_s:
            bursts[-1] = (bursts[-1][0], t_end)
        else:
            bursts.append((t_start, t_end))

    return bursts


def spectrogram_at_nsps(data, sr, nsps, fmin, fmax, t_start, t_end, title, outpath):
    """
    Generate spectrogram using FFT length matched to candidate NSPS.
    For 8-GFSK with h=1, using nfft=nsps gives frequency bins exactly at tone centers.
    """
    baud = sr / nsps
    nfft = nsps * 2  # 2x oversampling for cleaner display
    hop = nsps // 8  # fine time resolution

    # Extract the time window
    i_start = int(t_start * sr)
    i_end = int(t_end * sr)
    segment = data[i_start:i_end]

    freqs, times, Sxx = sig.spectrogram(
        segment, fs=sr, nperseg=nfft, noverlap=nfft - hop,
        window="hann", scaling="spectrum"
    )
    times = times + t_start

    mask = (freqs >= fmin) & (freqs <= fmax)
    Sxx_db = 10 * np.log10(Sxx[mask] + 1e-20)
    vmin = np.percentile(Sxx_db, 20)
    vmax = np.percentile(Sxx_db, 99.5)

    fig, ax = plt.subplots(figsize=(24, 8))
    im = ax.pcolormesh(times, freqs[mask], Sxx_db, shading="gouraud",
                        cmap="inferno", vmin=vmin, vmax=vmax)

    # Symbol boundary grid
    sym_dur = 1.0 / baud
    t = t_start
    while t < t_end:
        ax.axvline(t, color="white", alpha=0.25, linewidth=0.3)
        t += sym_dur

    # Tone frequency grid (assuming h=1, spacing = baud)
    tone_spacing = baud  # h=1
    center = (fmin + fmax) / 2
    for i in range(-5, 13):
        f = center + i * tone_spacing
        if fmin <= f <= fmax:
            ax.axhline(f, color="cyan", alpha=0.15, linewidth=0.5)

    df = sr / nfft
    ax.set_xlabel(f"Time (s)")
    ax.set_ylabel("Frequency (Hz)")
    ax.set_title(f"{title}\nNSPS={nsps}, baud={baud:.2f} Bd, "
                 f"tone_spacing={tone_spacing:.2f} Hz (h=1), "
                 f"df={df:.2f} Hz, sym_dur={1000/baud:.1f} ms")

    plt.colorbar(im, ax=ax, label="Power (dB)")
    plt.tight_layout()
    plt.savefig(outpath, dpi=200, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {outpath}")


def peak_frequency_per_frame(data, sr, nsps, fmin, fmax):
    """Extract peak frequency in each symbol-length frame."""
    nfft = nsps * 4  # oversampled FFT
    n_symbols = len(data) // nsps
    baud = sr / nsps

    peak_freqs = []
    for s in range(n_symbols):
        start = s * nsps
        end = start + nsps
        segment = data[start:end]

        # Apply window
        window = np.hanning(len(segment))
        windowed = segment * window

        # FFT
        spectrum = np.abs(fft(windowed, n=nfft))
        freqs = np.fft.fftfreq(nfft, 1 / sr)

        # Only positive frequencies in range
        pos_mask = (freqs >= fmin) & (freqs <= fmax)
        if not np.any(pos_mask):
            peak_freqs.append(0)
            continue

        peak_idx = np.argmax(spectrum[pos_mask])
        peak_freq = freqs[pos_mask][peak_idx]
        peak_freqs.append(peak_freq)

    return np.array(peak_freqs)


def analyze_burst_tones(data, sr, t_start, t_end, fmin, fmax, nsps, outpath):
    """Extract and analyze tone sequence from a single burst."""
    baud = sr / nsps
    i_start = int(t_start * sr)
    i_end = int(t_end * sr)
    segment = data[i_start:i_end]

    peak_freqs = peak_frequency_per_frame(segment, sr, nsps, fmin, fmax)

    # Remove frames where no signal (freq = 0 or very low power)
    valid = peak_freqs > 0
    if np.sum(valid) < 10:
        print(f"  Not enough valid frames for tone analysis")
        return None

    # Histogram of peak frequencies
    freq_range = (fmin, fmax)
    hist, edges = np.histogram(peak_freqs[valid], bins=300, range=freq_range)
    centers = (edges[:-1] + edges[1:]) / 2

    # Find tone clusters
    kernel = np.ones(3) / 3
    smooth_hist = np.convolve(hist, kernel, mode="same")
    peaks, props = sig.find_peaks(smooth_hist, height=np.max(smooth_hist) * 0.08,
                                   distance=int(baud * 0.5 / ((fmax - fmin) / 300)))

    tone_centers = sorted(centers[peaks])
    if len(tone_centers) >= 2:
        spacings = np.diff(tone_centers)
        print(f"  Tone centers: {[f'{f:.1f}' for f in tone_centers]}")
        print(f"  Spacings: {[f'{s:.1f}' for s in spacings]}")
        print(f"  Median spacing: {np.median(spacings):.2f} Hz")
        print(f"  h = {np.median(spacings)/baud:.3f}")

    # Plot
    fig, axes = plt.subplots(3, 1, figsize=(16, 10))

    ax = axes[0]
    t_axis = np.arange(len(peak_freqs)) / baud + t_start
    ax.plot(t_axis[valid], peak_freqs[valid], ".", markersize=2)
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Peak Frequency (Hz)")
    ax.set_title(f"Peak Frequency per Symbol (NSPS={nsps}, baud={baud:.2f})")
    ax.set_ylim(fmin, fmax)
    ax.grid(True, alpha=0.3)
    # Add tone grid
    if len(tone_centers) >= 2:
        for f in tone_centers:
            ax.axhline(f, color="red", alpha=0.3, linewidth=0.5)

    ax = axes[1]
    ax.bar(centers, hist, width=centers[1] - centers[0], alpha=0.6)
    ax.plot(centers, smooth_hist, "r-", linewidth=1)
    for f in tone_centers:
        ax.axvline(f, color="green", linestyle="--", alpha=0.7)
    ax.set_xlabel("Frequency (Hz)")
    ax.set_ylabel("Count")
    ax.set_title("Peak Frequency Histogram (tone detection)")

    # Quantized tone sequence
    ax = axes[2]
    if len(tone_centers) >= 2:
        spacing = np.median(spacings)
        base_freq = tone_centers[0]
        tone_indices = np.round((peak_freqs - base_freq) / spacing).astype(int)
        tone_indices = np.clip(tone_indices, -1, 9)
        ax.plot(tone_indices[valid], ".-", markersize=2, linewidth=0.5)
        ax.set_xlabel("Symbol index")
        ax.set_ylabel("Tone index")
        ax.set_title("Quantized Tone Sequence")
        ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(outpath, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {outpath}")

    return peak_freqs, tone_centers


def main():
    wav_path = sys.argv[1] if len(sys.argv) > 1 else "ft2_signal.wav"
    Path("output").mkdir(exist_ok=True)

    data, sr = load_wav(wav_path)
    print(f"Loaded: {sr} Hz, {len(data)/sr:.2f}s, {len(data)} samples")

    # Find bursts — use visible windows from overview spectrogram
    # The energy-based detector misses weak signals, so use manually identified windows
    print(f"\n=== Burst Detection (from spectrogram) ===")
    total_dur = len(data) / sr
    bursts_auto = find_bursts(data, sr, min_gap_s=0.3)
    print(f"  Auto-detected: {[(f'{t0:.2f}-{t1:.2f}') for t0,t1 in bursts_auto]}")

    # Use wider windows covering visible activity in the spectrogram
    # Multiple FT2 cycles visible across the ~14s recording
    bursts = [
        (0.5, 4.0),    # First cycle
        (4.5, 8.0),    # Second cycle
        (8.5, 13.5),   # Third cycle (possibly two overlapping)
    ]
    print(f"  Using manual windows: {[(f'{t0:.1f}-{t1:.1f}') for t0,t1 in bursts]}")

    # Two signal groups visible in spectrogram:
    # Group A: ~700-850 Hz
    # Group B: ~1150-1400 Hz
    signal_bands = [
        (700, 900, "Signal_A_700-900Hz"),
        (1150, 1400, "Signal_B_1150-1400Hz"),
    ]

    # Candidate NSPS values to test
    candidate_nsps = [576, 480, 360]

    # For each burst + signal band + candidate NSPS, generate spectrogram
    for bi, (t0, t1) in enumerate(bursts[:4]):  # first 4 bursts
        for fmin, fmax, band_name in signal_bands:
            for nsps in candidate_nsps:
                baud = sr / nsps
                outname = f"zoom_burst{bi}_{band_name}_nsps{nsps}.png"
                print(f"\n  Burst {bi} ({t0:.1f}-{t1:.1f}s), {band_name}, NSPS={nsps}")
                spectrogram_at_nsps(data, sr, nsps, fmin, fmax, t0, t1,
                                     f"Burst {bi} - {band_name}",
                                     f"output/{outname}")

    # Detailed tone analysis for the cleanest-looking burst (typically the first or strongest)
    print(f"\n=== Tone Analysis ===")
    if len(bursts) >= 1:
        for bi in range(min(3, len(bursts))):
            t0, t1 = bursts[bi]
            for fmin, fmax, band_name in signal_bands:
                for nsps in candidate_nsps:
                    baud = sr / nsps
                    print(f"\n--- Burst {bi}, {band_name}, NSPS={nsps} ({baud:.2f} Bd) ---")
                    result = analyze_burst_tones(data, sr, t0, t1, fmin, fmax, nsps,
                                                  f"output/tones_burst{bi}_{band_name}_nsps{nsps}.png")

    print(f"\n=== Done. Check output/ for spectrograms. ===")


if __name__ == "__main__":
    main()
