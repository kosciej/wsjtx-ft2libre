#!/usr/bin/env python3
"""
FT2 Ultra-Zoom — sub-second spectrogram windows to count individual symbols
and precisely measure baud rate and tone spacing.

Usage:
    uv run python analyze_ft2_ultrazoom.py ft2_signal.wav
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


def ultrazoom_spectrogram(data, sr, t_start, t_end, fmin, fmax, nfft, title, outpath):
    """Ultra-zoomed spectrogram with fine time/frequency grid."""
    i_start = int(t_start * sr)
    i_end = int(t_end * sr)
    segment = data[i_start:i_end]

    hop = nfft // 16  # very fine hop for smooth display
    freqs, times, Sxx = sig.spectrogram(
        segment, fs=sr, nperseg=nfft, noverlap=nfft - hop,
        window="hann", scaling="spectrum"
    )
    times = times + t_start

    mask = (freqs >= fmin) & (freqs <= fmax)
    Sxx_db = 10 * np.log10(Sxx[mask] + 1e-20)
    vmin = np.percentile(Sxx_db, 10)
    vmax = np.percentile(Sxx_db, 99.5)

    fig, ax = plt.subplots(figsize=(28, 10))
    im = ax.pcolormesh(times, freqs[mask], Sxx_db, shading="gouraud",
                        cmap="inferno", vmin=vmin, vmax=vmax)

    df = sr / nfft
    ax.set_xlabel(f"Time (s)", fontsize=12)
    ax.set_ylabel("Frequency (Hz)", fontsize=12)
    ax.set_title(f"{title}\nnfft={nfft}, df={df:.2f} Hz, hop={hop}", fontsize=11)

    # Fine frequency grid every df Hz
    for f in np.arange(fmin, fmax, df):
        ax.axhline(f, color="cyan", alpha=0.08, linewidth=0.3)

    plt.colorbar(im, ax=ax, label="Power (dB)")
    plt.tight_layout()
    plt.savefig(outpath, dpi=250, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {outpath}")


def symbol_peak_analysis(data, sr, t_start, t_end, fmin, fmax, nsps_list):
    """For each candidate NSPS, measure peak sharpness per symbol."""
    i_start = int(t_start * sr)
    i_end = int(t_end * sr)
    segment = data[i_start:i_end]

    print(f"\n=== Symbol Peak Sharpness Analysis ({t_start:.2f}-{t_end:.2f}s) ===")

    results = []
    for nsps in nsps_list:
        baud = sr / nsps
        nfft = nsps * 8  # high oversampling
        n_symbols = len(segment) // nsps
        if n_symbols < 3:
            continue

        peak_to_avg = []
        peak_freqs_list = []
        for s in range(n_symbols):
            start = s * nsps
            end = start + nsps
            sym = segment[start:end]

            window = np.hanning(len(sym))
            spectrum = np.abs(fft(sym * window, n=nfft))[:nfft // 2]
            freqs = np.arange(nfft // 2) * sr / nfft

            band = (freqs >= fmin) & (freqs <= fmax)
            if not np.any(band):
                continue
            spec_band = spectrum[band]
            freqs_band = freqs[band]

            if spec_band.max() > 0:
                peak_to_avg.append(spec_band.max() / (spec_band.mean() + 1e-20))
                peak_freqs_list.append(freqs_band[np.argmax(spec_band)])

        if peak_to_avg:
            avg_pta = np.mean(peak_to_avg)
            med_pta = np.median(peak_to_avg)
            results.append((nsps, baud, avg_pta, med_pta, n_symbols, peak_freqs_list))
            print(f"  NSPS={nsps:4d} ({baud:6.2f} Bd): "
                  f"peak/avg={avg_pta:6.2f} (median={med_pta:6.2f}), "
                  f"n_sym={n_symbols}")

    # Sort by median peak-to-average ratio
    results.sort(key=lambda x: x[3], reverse=True)
    print(f"\n  Best: NSPS={results[0][0]} ({results[0][1]:.2f} Bd)")
    return results


def plot_symbol_aligned_waterfall(data, sr, t_start, t_end, fmin, fmax,
                                    nsps, title, outpath):
    """
    Each row = one symbol, columns = frequency bins.
    This makes the tone sequence very clear.
    """
    i_start = int(t_start * sr)
    i_end = int(t_end * sr)
    segment = data[i_start:i_end]

    baud = sr / nsps
    nfft = nsps * 4
    n_symbols = len(segment) // nsps

    waterfall = np.zeros((n_symbols, nfft // 2))
    freqs = np.arange(nfft // 2) * sr / nfft

    for s in range(n_symbols):
        start = s * nsps
        end = start + nsps
        sym = segment[start:end]
        window = np.hanning(len(sym))
        spectrum = np.abs(fft(sym * window, n=nfft))[:nfft // 2]
        waterfall[s, :] = spectrum

    mask = (freqs >= fmin) & (freqs <= fmax)
    wf_band = waterfall[:, mask]
    freqs_band = freqs[mask]

    wf_db = 10 * np.log10(wf_band + 1e-20)
    vmin = np.percentile(wf_db, 20)
    vmax = np.percentile(wf_db, 99)

    fig, ax = plt.subplots(figsize=(16, max(6, n_symbols * 0.15)))
    im = ax.imshow(wf_db, aspect="auto", origin="lower",
                    extent=[freqs_band[0], freqs_band[-1], 0, n_symbols],
                    cmap="inferno", vmin=vmin, vmax=vmax)
    ax.set_xlabel("Frequency (Hz)")
    ax.set_ylabel("Symbol index")
    ax.set_title(f"{title}\nNSPS={nsps}, baud={baud:.2f} Bd, df={sr/nfft:.2f} Hz")

    # Mark symbol boundaries
    for i in range(n_symbols):
        ax.axhline(i, color="white", alpha=0.1, linewidth=0.3)

    plt.colorbar(im, ax=ax, label="Power (dB)")
    plt.tight_layout()
    plt.savefig(outpath, dpi=200, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {outpath}")

    # Return peak frequencies per symbol
    peak_freqs = freqs_band[np.argmax(wf_band, axis=1)]
    return peak_freqs


def plot_peak_freq_trace(peak_freqs_dict, t_start, outpath):
    """Compare peak frequency traces for different NSPS values."""
    fig, axes = plt.subplots(len(peak_freqs_dict), 1,
                              figsize=(18, 3 * len(peak_freqs_dict)), sharex=False)
    if len(peak_freqs_dict) == 1:
        axes = [axes]

    for ax, (nsps, (baud, peak_freqs)) in zip(axes, peak_freqs_dict.items()):
        times = np.arange(len(peak_freqs)) / baud + t_start
        ax.plot(times, peak_freqs, ".-", markersize=3, linewidth=0.8)
        ax.set_ylabel("Freq (Hz)")
        ax.set_title(f"NSPS={nsps} ({baud:.2f} Bd) — {len(peak_freqs)} symbols")
        ax.grid(True, alpha=0.3)

    axes[-1].set_xlabel("Time (s)")
    plt.tight_layout()
    plt.savefig(outpath, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {outpath}")


def measure_tone_spacing_from_peaks(peak_freqs, baud, label=""):
    """Precise tone spacing measurement from peak frequency differences."""
    # Compute all pairwise differences
    diffs = []
    for i in range(1, len(peak_freqs)):
        d = abs(peak_freqs[i] - peak_freqs[i - 1])
        if d > 1:  # ignore zero-change
            diffs.append(d)

    if not diffs:
        print(f"  {label}: No frequency transitions found")
        return 0

    diffs = np.array(diffs)

    # The smallest common difference should be the tone spacing
    # Use histogram to find it
    max_diff = np.percentile(diffs, 95)
    hist, edges = np.histogram(diffs, bins=100, range=(1, max_diff))
    centers = (edges[:-1] + edges[1:]) / 2

    # Find peaks in the histogram
    peaks, props = sig.find_peaks(hist, height=max(hist) * 0.1)

    if len(peaks) > 0:
        # The first peak should be the fundamental tone spacing
        tone_spacing = centers[peaks[0]]
        print(f"  {label}: Smallest common step = {tone_spacing:.2f} Hz, "
              f"h = {tone_spacing / baud:.3f}")

        # Check if higher peaks are multiples
        for p in peaks[1:4]:
            ratio = centers[p] / tone_spacing
            print(f"    Step {centers[p]:.2f} Hz = {ratio:.2f}x base")

        return tone_spacing

    return 0


def main():
    wav_path = sys.argv[1] if len(sys.argv) > 1 else "ft2_signal.wav"
    Path("output").mkdir(exist_ok=True)

    data, sr = load_wav(wav_path)
    print(f"Loaded: {sr} Hz, {len(data) / sr:.2f}s")

    # Focus on the strongest signal area from previous analysis
    # Burst 0: Signal B (1150-1400 Hz) has clearest FT2 around t=0.6-2.5s
    # Burst 1: Signal B around t=5.5-8.0s also strong

    windows = [
        (0.6, 2.8, 1150, 1400, "Burst0_SigB"),
        (0.6, 2.8, 680, 920, "Burst0_SigA"),
        (5.0, 8.0, 1150, 1400, "Burst1_SigB"),
        (9.0, 13.0, 1150, 1400, "Burst2_SigB"),
    ]

    candidate_nsps = [192, 240, 288, 320, 360, 384, 420, 448, 480, 512, 540, 576, 640, 720, 768, 960]

    # Step 1: Peak sharpness analysis for each window
    print(f"\n{'='*60}")
    print(f"  STEP 1: SYMBOL RATE SCAN")
    print(f"{'='*60}")

    best_nsps_per_window = {}
    for t0, t1, fmin, fmax, label in windows:
        print(f"\n  Window: {label} ({t0:.1f}-{t1:.1f}s, {fmin}-{fmax} Hz)")
        results = symbol_peak_analysis(data, sr, t0, t1, fmin, fmax, candidate_nsps)
        if results:
            best_nsps_per_window[label] = results[0][0]

    # Step 2: Symbol-aligned waterfalls for top candidates
    print(f"\n{'='*60}")
    print(f"  STEP 2: SYMBOL-ALIGNED WATERFALLS")
    print(f"{'='*60}")

    # Use the first clean window for detailed analysis
    t0, t1, fmin, fmax, label = windows[0]
    top_nsps = [576, 480, 360]  # Always test these

    # Add any unique best candidates
    for lbl, nsps in best_nsps_per_window.items():
        if nsps not in top_nsps:
            top_nsps.append(nsps)

    peak_freqs_dict = {}
    for nsps in sorted(top_nsps):
        baud = sr / nsps
        print(f"\n  NSPS={nsps} ({baud:.2f} Bd)")
        pf = plot_symbol_aligned_waterfall(
            data, sr, t0, t1, fmin, fmax, nsps,
            f"{label} Symbol Waterfall",
            f"output/waterfall_{label}_nsps{nsps}.png"
        )
        peak_freqs_dict[nsps] = (baud, pf)
        measure_tone_spacing_from_peaks(pf, baud, f"NSPS={nsps}")

    # Step 3: Compare peak frequency traces
    plot_peak_freq_trace(peak_freqs_dict, t0, f"output/peak_traces_{label}.png")

    # Step 4: Ultra-zoomed spectrograms (0.5s windows)
    print(f"\n{'='*60}")
    print(f"  STEP 3: ULTRA-ZOOM SPECTROGRAMS")
    print(f"{'='*60}")

    for nsps in [576, 480, 360]:
        baud = sr / nsps
        # Zoom into 0.5s of the cleanest part
        t_zoom_start = 1.0
        t_zoom_end = 1.5
        ultrazoom_spectrogram(
            data, sr, t_zoom_start, t_zoom_end, fmin, fmax,
            nfft=nsps * 2,
            title=f"Ultra-Zoom NSPS={nsps} ({baud:.2f} Bd)",
            outpath=f"output/ultrazoom_nsps{nsps}.png"
        )

    print(f"\n{'='*60}")
    print(f"  DONE — check output/ directory")
    print(f"{'='*60}")


if __name__ == "__main__":
    main()
