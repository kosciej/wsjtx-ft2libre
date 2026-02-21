#!/usr/bin/env python3
"""
FT2 Signal Analyzer — Reverse-engineer FT2 protocol parameters from audio samples.

Measures: baud rate, tone count, tone spacing, modulation index, bandwidth,
sync pattern, and symbol structure. Compares against FT8/FT4 known parameters.

Usage:
    uv run python analyze_ft2_signal.py <wav_file> [--freq <center_freq_hz>] [--output <dir>]
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import soundfile as sf
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import signal as sig
from scipy.fft import fft, fftfreq, ifft


# Known FT8 Costas 7x7 sync pattern
COSTAS_7x7 = np.array([3, 1, 4, 0, 6, 5, 2])
# FT8 sync symbol positions within 79 symbols
FT8_SYNC_POS = [0, 36, 72]  # three Costas arrays


def load_audio(wav_path: str) -> tuple[np.ndarray, int]:
    """Load WAV file, convert to mono float64, return (samples, sample_rate)."""
    data, sr = sf.read(wav_path, dtype="float64")
    if data.ndim > 1:
        data = data.mean(axis=1)
    print(f"Loaded: {wav_path}")
    print(f"  Sample rate: {sr} Hz")
    print(f"  Duration: {len(data)/sr:.3f} s")
    print(f"  Samples: {len(data)}")
    return data, sr


def resample_to_12000(data: np.ndarray, sr: int) -> np.ndarray:
    """Resample to 12000 Hz if needed (WSJT-X standard rate)."""
    if sr == 12000:
        return data
    target_sr = 12000
    num_samples = int(len(data) * target_sr / sr)
    resampled = sig.resample(data, num_samples)
    print(f"  Resampled: {sr} -> {target_sr} Hz ({num_samples} samples)")
    return resampled


def compute_spectrogram(data: np.ndarray, sr: int, nfft: int = 1024,
                         overlap_frac: float = 0.875) -> tuple:
    """Compute spectrogram, return (Sxx, freqs, times)."""
    noverlap = int(nfft * overlap_frac)
    freqs, times, Sxx = sig.spectrogram(
        data, fs=sr, nperseg=nfft, noverlap=noverlap,
        window="hann", scaling="spectrum"
    )
    return Sxx, freqs, times


def find_signal_band(Sxx: np.ndarray, freqs: np.ndarray,
                      min_freq: float = 200, max_freq: float = 4000) -> tuple:
    """Find the frequency band containing the FT2 signal."""
    mask = (freqs >= min_freq) & (freqs <= max_freq)
    avg_power = Sxx[mask].mean(axis=1)
    freqs_sub = freqs[mask]

    # Smooth and find peak
    kernel_size = max(5, len(avg_power) // 50)
    if kernel_size % 2 == 0:
        kernel_size += 1
    smoothed = np.convolve(avg_power, np.ones(kernel_size) / kernel_size, mode="same")

    peak_idx = np.argmax(smoothed)
    peak_freq = freqs_sub[peak_idx]

    # Find -10dB bandwidth
    threshold = smoothed[peak_idx] * 0.1
    above = smoothed > threshold
    indices = np.where(above)[0]
    if len(indices) > 0:
        f_low = freqs_sub[indices[0]]
        f_high = freqs_sub[indices[-1]]
    else:
        f_low = peak_freq - 100
        f_high = peak_freq + 100

    center = (f_low + f_high) / 2
    bandwidth = f_high - f_low

    print(f"\nSignal band detection:")
    print(f"  Peak frequency: {peak_freq:.1f} Hz")
    print(f"  Band: {f_low:.1f} - {f_high:.1f} Hz")
    print(f"  Center: {center:.1f} Hz")
    print(f"  Bandwidth (-10dB): {bandwidth:.1f} Hz")
    return center, f_low, f_high, bandwidth


def extract_baseband(data: np.ndarray, sr: int, center_freq: float,
                      bandwidth: float = 300) -> np.ndarray:
    """Mix signal to baseband and low-pass filter."""
    t = np.arange(len(data)) / sr
    # Mix down
    analytic = sig.hilbert(data)
    baseband = analytic * np.exp(-2j * np.pi * center_freq * t)
    # Low-pass filter
    cutoff = bandwidth / 2
    nyq = sr / 2
    b, a = sig.butter(6, cutoff / nyq, btype="low")
    baseband = sig.filtfilt(b, a, baseband)
    return baseband


def measure_baud_rate(data: np.ndarray, sr: int, center_freq: float,
                       bw: float, nfft: int = 1024) -> float:
    """Measure symbol rate from autocorrelation of spectral energy envelope."""
    print(f"\n--- Baud Rate Measurement ---")

    # Extract baseband signal
    baseband = extract_baseband(data, sr, center_freq, bw * 1.5)

    # Compute instantaneous frequency changes (tone transitions)
    phase = np.unwrap(np.angle(baseband))
    inst_freq = np.diff(phase) * sr / (2 * np.pi)

    # Smooth slightly
    win = max(3, sr // 2000)
    inst_freq_smooth = np.convolve(inst_freq, np.ones(win) / win, mode="same")

    # Autocorrelation of the instantaneous frequency derivative (detects transitions)
    freq_diff = np.diff(inst_freq_smooth)
    freq_diff = freq_diff - freq_diff.mean()

    # Use power of freq_diff as transition detector
    power = freq_diff ** 2
    power = power - power.mean()

    # Autocorrelation via FFT
    n = len(power)
    fft_power = fft(power, n=2 * n)
    acf = np.real(ifft(fft_power * np.conj(fft_power)))[:n]
    acf = acf / acf[0]

    # Find first significant peak after zero
    min_lag = int(sr * 0.005)  # minimum 5ms symbol (200 Bd max)
    max_lag = int(sr * 0.2)    # maximum 200ms symbol (5 Bd min)
    max_lag = min(max_lag, len(acf) - 1)

    search = acf[min_lag:max_lag]
    peaks, properties = sig.find_peaks(search, height=0.05, distance=int(sr * 0.003))

    if len(peaks) == 0:
        print("  WARNING: No clear periodicity found in autocorrelation")
        return 0.0

    # The first strong peak gives the symbol period
    peak_heights = properties["peak_heights"]
    best_peak = peaks[np.argmax(peak_heights)]
    symbol_samples = best_peak + min_lag
    symbol_duration = symbol_samples / sr
    baud = 1.0 / symbol_duration

    print(f"  Symbol period: {symbol_samples} samples = {symbol_duration*1000:.2f} ms")
    print(f"  Baud rate: {baud:.2f} Bd")
    print(f"  NSPS (at 12000): {12000/baud:.1f}")

    # Check top peaks
    print(f"\n  Top autocorrelation peaks:")
    sorted_idx = np.argsort(peak_heights)[::-1][:5]
    for i, idx in enumerate(sorted_idx):
        lag = peaks[idx] + min_lag
        dur = lag / sr * 1000
        rate = sr / lag
        print(f"    #{i+1}: lag={lag} samples, {dur:.2f} ms, {rate:.2f} Bd (height={peak_heights[idx]:.3f})")

    return baud


def measure_tones(data: np.ndarray, sr: int, center_freq: float,
                   baud: float, bw: float) -> tuple:
    """Analyze tone structure — count tones, measure spacing."""
    print(f"\n--- Tone Analysis ---")

    if baud <= 0:
        print("  Cannot analyze tones without valid baud rate")
        return 0, 0.0

    nsps = int(round(sr / baud))
    print(f"  Using NSPS={nsps} (baud={baud:.2f})")

    # Use a longer FFT for better frequency resolution
    nfft = nsps * 4  # 4x oversampling in frequency

    # Compute spectrogram at symbol rate
    noverlap = nsps * 3  # 75% overlap
    freqs, times, Sxx = sig.spectrogram(
        data, fs=sr, nperseg=nfft, noverlap=noverlap,
        window="hann", scaling="spectrum"
    )

    # Focus on signal band
    f_margin = 50
    mask = (freqs >= center_freq - bw / 2 - f_margin) & (freqs <= center_freq + bw / 2 + f_margin)
    Sxx_band = Sxx[mask]
    freqs_band = freqs[mask]

    # Average power spectrum across time
    avg_spectrum = Sxx_band.mean(axis=1)

    # Also look at peak frequencies per time slice
    peak_freqs = freqs_band[np.argmax(Sxx_band, axis=0)]

    # Histogram of peak frequencies to find tone centers
    freq_range = (center_freq - bw / 2 - f_margin, center_freq + bw / 2 + f_margin)
    hist, bin_edges = np.histogram(peak_freqs, bins=200, range=freq_range)
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2

    # Find clusters (tones)
    hist_smooth = np.convolve(hist, np.ones(3) / 3, mode="same")
    peaks, props = sig.find_peaks(hist_smooth, height=np.max(hist_smooth) * 0.1,
                                   distance=5)

    if len(peaks) >= 2:
        tone_freqs = sorted(bin_centers[peaks])
        spacings = np.diff(tone_freqs)
        avg_spacing = np.median(spacings)
        n_tones = len(tone_freqs)

        print(f"  Detected {n_tones} tones")
        print(f"  Tone frequencies: {[f'{f:.1f}' for f in tone_freqs]}")
        print(f"  Spacings: {[f'{s:.1f}' for s in spacings]}")
        print(f"  Median tone spacing: {avg_spacing:.2f} Hz")
        print(f"  Modulation index h = spacing/baud = {avg_spacing/baud:.3f}")

        return n_tones, avg_spacing
    else:
        print("  Could not reliably detect individual tones")
        return 0, 0.0


def detect_costas_sync(data: np.ndarray, sr: int, center_freq: float,
                        baud: float, tone_spacing: float) -> dict:
    """Try to detect Costas 7x7 sync arrays in the signal."""
    print(f"\n--- Costas Sync Detection ---")

    if baud <= 0 or tone_spacing <= 0:
        print("  Cannot detect sync without valid baud and tone spacing")
        return {}

    nsps = int(round(sr / baud))
    n_symbols = int(len(data) / nsps)
    print(f"  Total symbols in recording: ~{n_symbols}")

    # Extract tone sequence by finding peak frequency per symbol
    baseband = extract_baseband(data, sr, center_freq, tone_spacing * 10)
    phase = np.unwrap(np.angle(baseband))

    tone_sequence = []
    for i in range(n_symbols):
        start = i * nsps
        end = start + nsps
        if end >= len(phase):
            break
        # Average instantaneous frequency in this symbol
        segment_phase = phase[start:end]
        avg_freq = np.mean(np.diff(segment_phase)) * sr / (2 * np.pi)
        tone_idx = round(avg_freq / tone_spacing)
        tone_sequence.append(tone_idx)

    tone_sequence = np.array(tone_sequence)
    # Normalize to 0-7 range
    tone_min = tone_sequence.min()
    tone_norm = tone_sequence - tone_min

    print(f"  Extracted {len(tone_norm)} symbol tones")
    print(f"  Tone range: {tone_norm.min()} to {tone_norm.max()}")

    # Try to find Costas 7x7 pattern
    costas = COSTAS_7x7
    best_score = 0
    best_offset = 0

    scores = []
    for offset in range(len(tone_norm) - 79):
        score = 0
        for sync_start in FT8_SYNC_POS:
            segment = tone_norm[offset + sync_start:offset + sync_start + 7]
            if len(segment) == 7:
                # Correlation with Costas pattern (try all offsets)
                for tone_off in range(-3, 11):
                    matches = np.sum((segment - tone_off) == costas)
                    score = max(score, matches)
        scores.append(score)
        if score > best_score:
            best_score = score
            best_offset = offset

    scores = np.array(scores)
    print(f"  Best Costas correlation: {best_score}/21 at symbol offset {best_offset}")

    if best_score >= 14:
        print(f"  STRONG Costas 7x7 sync detected!")
        # Extract the actual sync tones
        for i, pos in enumerate(FT8_SYNC_POS):
            seg = tone_norm[best_offset + pos:best_offset + pos + 7]
            print(f"    Costas block {i} (pos {pos}): {seg.tolist()}")
    elif best_score >= 10:
        print(f"  Possible Costas sync (moderate correlation)")
    else:
        print(f"  Costas 7x7 sync NOT clearly detected")

    return {
        "best_score": best_score,
        "best_offset": best_offset,
        "n_symbols": len(tone_norm),
        "tone_sequence": tone_norm,
        "scores": scores,
    }


def measure_timing(data: np.ndarray, sr: int, baud: float) -> dict:
    """Measure total transmission duration and symbol count."""
    print(f"\n--- Timing Analysis ---")

    if baud <= 0:
        print("  Cannot measure timing without valid baud rate")
        return {}

    # Find signal start/end using energy envelope
    frame_size = int(sr / baud)  # one symbol
    n_frames = len(data) // frame_size
    energy = np.array([
        np.sum(data[i * frame_size:(i + 1) * frame_size] ** 2)
        for i in range(n_frames)
    ])

    threshold = np.max(energy) * 0.01
    active = energy > threshold
    if not np.any(active):
        print("  No signal detected")
        return {}

    first_symbol = np.argmax(active)
    last_symbol = len(active) - 1 - np.argmax(active[::-1])
    n_active = last_symbol - first_symbol + 1
    tx_duration = n_active / baud

    print(f"  Signal start: symbol {first_symbol} ({first_symbol/baud:.3f} s)")
    print(f"  Signal end: symbol {last_symbol} ({last_symbol/baud:.3f} s)")
    print(f"  Active symbols: {n_active}")
    print(f"  TX duration: {tx_duration:.3f} s")
    print(f"  Expected for 79 symbols: {79/baud:.3f} s")

    return {
        "first_symbol": first_symbol,
        "last_symbol": last_symbol,
        "n_active": n_active,
        "tx_duration": tx_duration,
    }


def plot_spectrogram(data: np.ndarray, sr: int, center_freq: float,
                      bw: float, baud: float, output_path: str):
    """Generate detailed spectrogram plot."""
    nsps = max(int(round(sr / baud)), 256) if baud > 0 else 512
    nfft = nsps * 2

    fig, axes = plt.subplots(3, 1, figsize=(16, 14))

    # Full spectrogram
    ax = axes[0]
    freqs, times, Sxx = sig.spectrogram(
        data, fs=sr, nperseg=nfft, noverlap=int(nfft * 0.875),
        window="hann", scaling="spectrum"
    )
    Sxx_db = 10 * np.log10(Sxx + 1e-20)
    ax.pcolormesh(times, freqs, Sxx_db, shading="gouraud", cmap="viridis")
    ax.set_ylabel("Frequency (Hz)")
    ax.set_title("Full Spectrogram")
    ax.set_ylim(0, sr / 2)

    # Zoomed spectrogram around signal
    ax = axes[1]
    margin = max(bw * 0.5, 100)
    f_lo = center_freq - bw / 2 - margin
    f_hi = center_freq + bw / 2 + margin
    mask = (freqs >= f_lo) & (freqs <= f_hi)
    ax.pcolormesh(times, freqs[mask], Sxx_db[mask], shading="gouraud", cmap="viridis")
    ax.set_ylabel("Frequency (Hz)")
    ax.set_title(f"Zoomed Spectrogram (center={center_freq:.0f} Hz, BW={bw:.0f} Hz)")

    # Add symbol grid if baud is known
    if baud > 0:
        symbol_dur = 1.0 / baud
        t = 0
        while t < times[-1]:
            ax.axvline(t, color="white", alpha=0.15, linewidth=0.5)
            t += symbol_dur

    # Average power spectrum
    ax = axes[2]
    avg_psd = Sxx.mean(axis=1)
    ax.semilogy(freqs, avg_psd)
    ax.set_xlabel("Frequency (Hz)")
    ax.set_ylabel("Power")
    ax.set_title("Average Power Spectrum")
    ax.set_xlim(f_lo, f_hi)
    ax.axvline(center_freq, color="r", alpha=0.5, linestyle="--", label="Center")
    ax.legend()

    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches="tight")
    print(f"\nSpectrogram saved to: {output_path}")
    plt.close()


def plot_tone_analysis(data: np.ndarray, sr: int, center_freq: float,
                        baud: float, tone_spacing: float, sync_info: dict,
                        output_path: str):
    """Plot tone sequence and sync correlation."""
    fig, axes = plt.subplots(2, 1, figsize=(16, 8))

    if "tone_sequence" in sync_info:
        tones = sync_info["tone_sequence"]
        ax = axes[0]
        ax.plot(tones, ".-", markersize=2, linewidth=0.5)
        ax.set_xlabel("Symbol index")
        ax.set_ylabel("Tone index")
        ax.set_title("Extracted Tone Sequence")
        ax.grid(True, alpha=0.3)

        # Mark Costas sync positions if detected
        if sync_info.get("best_score", 0) >= 10:
            off = sync_info["best_offset"]
            for pos in FT8_SYNC_POS:
                ax.axvspan(off + pos, off + pos + 7, alpha=0.2, color="red",
                          label="Costas" if pos == 0 else None)
            ax.axvline(off, color="green", linestyle="--", alpha=0.5, label="Message start")
            ax.axvline(off + 79, color="green", linestyle="--", alpha=0.5, label="Message end")
            ax.legend()

    if "scores" in sync_info:
        ax = axes[1]
        ax.plot(sync_info["scores"])
        ax.set_xlabel("Symbol offset")
        ax.set_ylabel("Costas correlation score (/21)")
        ax.set_title("Costas 7x7 Sync Correlation")
        ax.axhline(14, color="red", linestyle="--", alpha=0.5, label="Strong threshold")
        ax.grid(True, alpha=0.3)
        ax.legend()

    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches="tight")
    print(f"Tone analysis saved to: {output_path}")
    plt.close()


def print_comparison(baud: float, n_tones: int, tone_spacing: float, bw: float):
    """Compare measured parameters with known modes."""
    print(f"\n{'='*60}")
    print(f"  PARAMETER COMPARISON")
    print(f"{'='*60}")
    print(f"{'Parameter':<20} {'Measured':<15} {'FT8':<12} {'FT4':<12} {'FT2Libre':<12}")
    print(f"{'-'*60}")
    print(f"{'Baud rate (Bd)':<20} {baud:<15.2f} {'6.25':<12} {'20.83':<12} {'33.33':<12}")

    nsps_meas = 12000 / baud if baud > 0 else 0
    print(f"{'NSPS @12kHz':<20} {nsps_meas:<15.1f} {'1920':<12} {'576':<12} {'360':<12}")
    print(f"{'Tones':<20} {n_tones:<15d} {'8':<12} {'4':<12} {'8':<12}")
    print(f"{'Tone spacing (Hz)':<20} {tone_spacing:<15.2f} {'6.25':<12} {'20.83':<12} {'25.00':<12}")

    h_meas = tone_spacing / baud if baud > 0 else 0
    print(f"{'Mod index h':<20} {h_meas:<15.3f} {'1.0':<12} {'1.0':<12} {'0.75':<12}")
    print(f"{'Bandwidth (Hz)':<20} {bw:<15.1f} {'~50':<12} {'~83':<12} {'~175':<12}")

    sym_dur = 1000 / baud if baud > 0 else 0
    print(f"{'Symbol dur (ms)':<20} {sym_dur:<15.2f} {'160':<12} {'48':<12} {'30':<12}")

    tx_79 = 79 / baud if baud > 0 else 0
    print(f"{'79-sym TX (s)':<20} {tx_79:<15.3f} {'12.64':<12} {'3.79':<12} {'2.37':<12}")
    print(f"{'='*60}")


def main():
    parser = argparse.ArgumentParser(description="FT2 Signal Analyzer")
    parser.add_argument("wav_file", help="Path to WAV file containing FT2 signal")
    parser.add_argument("--freq", type=float, default=0,
                        help="Approximate center frequency of signal (Hz). Auto-detected if 0.")
    parser.add_argument("--output", default="./output",
                        help="Output directory for plots (default: ./output)")
    args = parser.parse_args()

    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Load audio
    data, sr = load_audio(args.wav_file)

    # Resample to 12000 Hz
    data = resample_to_12000(data, sr)
    sr = 12000

    # Compute spectrogram and find signal
    Sxx, freqs, times = compute_spectrogram(data, sr)
    center, f_low, f_high, bw = find_signal_band(Sxx, freqs)

    if args.freq > 0:
        center = args.freq
        print(f"  Using user-specified center: {center:.1f} Hz")

    # Measure baud rate
    baud = measure_baud_rate(data, sr, center, bw)

    # Measure tone structure
    n_tones, tone_spacing = measure_tones(data, sr, center, baud, bw)

    # Detect Costas sync
    sync_info = detect_costas_sync(data, sr, center, baud, tone_spacing)

    # Measure timing
    timing = measure_timing(data, sr, baud)

    # Generate plots
    plot_spectrogram(data, sr, center, bw, baud,
                     str(output_dir / "spectrogram.png"))

    if sync_info:
        plot_tone_analysis(data, sr, center, baud, tone_spacing, sync_info,
                          str(output_dir / "tone_analysis.png"))

    # Print comparison
    print_comparison(baud, n_tones, tone_spacing, bw)

    print(f"\nAll outputs saved to: {output_dir}/")


if __name__ == "__main__":
    main()
