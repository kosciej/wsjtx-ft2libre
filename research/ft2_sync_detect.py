#!/usr/bin/env python3
"""
FT2 Sync Detector — uses the FT8 Costas 7x7 sync pattern to detect
FT2 signals and determine the correct NSPS/baud rate.

This mimics what sync8.f90 does in WSJT-X: compute symbol spectra,
then correlate with Costas arrays at expected positions.

Usage:
    uv run python ft2_sync_detect.py ft2_signal.wav
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


# FT8 Costas 7x7 sync pattern and positions in 79-symbol frame
COSTAS = np.array([3, 1, 4, 0, 6, 5, 2])
SYNC_POS = [0, 36, 72]  # Start of each Costas block within 79 symbols


def load_wav(path):
    data, sr = sf.read(path, dtype="float64")
    if data.ndim > 1:
        data = data.mean(axis=1)
    return data, sr


def compute_symbol_spectra(data, sr, nsps):
    """
    Compute power spectrum for each symbol-length window.
    Uses nfft = 2*nsps for 2x frequency oversampling (like FT8).
    Step = nsps/4 for fine time resolution (like FT8's NSTEP).

    Returns: spectra array [n_time_steps x n_freq_bins], freq array
    """
    nfft = 2 * nsps
    step = nsps // 4
    n_steps = (len(data) - nfft) // step + 1

    if n_steps <= 0:
        return np.array([]), np.array([])

    spectra = np.zeros((n_steps, nfft // 2))
    window = np.hanning(nfft)

    for i in range(n_steps):
        start = i * step
        segment = data[start:start + nfft] * window
        spec = np.abs(fft(segment, n=nfft))[:nfft // 2]
        spectra[i] = spec ** 2  # power spectrum

    freqs = np.arange(nfft // 2) * sr / nfft
    return spectra, freqs


def sync_search(spectra, freqs, nsps, sr):
    """
    Search for FT8-style Costas sync in the symbol spectra.
    Scans all time offsets and frequency offsets.

    For h=1: tone spacing = 1 bin in nfft = 2*nsps
    For other h: tone spacing = nfos*h bins where nfos = nfft/nsps = 2

    Returns: array of sync scores indexed by (time_step, freq_bin)
    """
    nfft = 2 * nsps
    nfos = nfft // nsps  # frequency oversampling factor = 2
    step = nsps // 4
    nssy = nsps // step  # time steps per symbol = 4

    n_steps = spectra.shape[0]
    n_freq = spectra.shape[1]

    # For the 79-symbol frame, need nssy*79 time steps
    # Total symbols needed: 79
    n_sym_steps = 79 * nssy  # time steps for full frame

    # Try different modulation indices
    h_candidates = [0.5, 0.75, 1.0]

    best_overall = {
        "score": 0, "t_step": 0, "f_bin": 0, "h": 0, "t_sec": 0, "f_hz": 0
    }

    results_per_h = {}

    for h in h_candidates:
        tone_bins = nfos * h  # tone spacing in FFT bins
        if tone_bins < 0.5:
            continue

        # For integer tone_bins, use direct bin lookup
        # For fractional, interpolate

        max_freq_bin = int(n_freq - 8 * tone_bins)
        if max_freq_bin <= 0:
            continue

        # Scan all time offsets and frequency offsets
        scores = np.zeros((n_steps - n_sym_steps, max_freq_bin))

        for jt in range(scores.shape[0]):
            for jf in range(scores.shape[1]):
                score = 0
                # Sum power at Costas tone positions for each sync block
                for sync_start in SYNC_POS:
                    for k in range(7):
                        # Symbol time step
                        t_idx = jt + (sync_start + k) * nssy
                        if t_idx >= n_steps:
                            continue
                        # Expected frequency bin for this Costas tone
                        tone = COSTAS[k]
                        f_idx = jf + int(round(tone * tone_bins))
                        if 0 <= f_idx < n_freq:
                            score += spectra[t_idx, f_idx]

                # Subtract average power in non-sync positions for normalization
                # (simplified — just use the sync sum)
                scores[jt, jf] = score

        # Find peak
        if scores.size > 0:
            peak_idx = np.unravel_index(np.argmax(scores), scores.shape)
            peak_score = scores[peak_idx]

            # Normalize by noise floor
            noise_floor = np.median(scores)
            snr = peak_score / (noise_floor + 1e-20)

            t_sec = peak_idx[0] * step / sr
            f_hz = freqs[peak_idx[1]]

            results_per_h[h] = {
                "score": peak_score,
                "snr": snr,
                "t_step": peak_idx[0],
                "f_bin": peak_idx[1],
                "t_sec": t_sec,
                "f_hz": f_hz,
                "scores": scores,
            }

            if snr > best_overall.get("snr", 0):
                best_overall = {
                    "score": peak_score,
                    "snr": snr,
                    "t_step": peak_idx[0],
                    "f_bin": peak_idx[1],
                    "h": h,
                    "t_sec": t_sec,
                    "f_hz": f_hz,
                }

    return best_overall, results_per_h


def main():
    wav_path = sys.argv[1] if len(sys.argv) > 1 else "ft2_signal.wav"
    Path("output").mkdir(exist_ok=True)

    data, sr = load_wav(wav_path)
    print(f"Loaded: {sr} Hz, {len(data) / sr:.2f}s")

    # Candidate NSPS values
    candidate_nsps = [360, 384, 420, 448, 480, 512, 540, 576, 640, 720]

    print(f"\n{'='*70}")
    print(f"  FT2 SYNC DETECTION — Costas 7x7 correlation")
    print(f"{'='*70}")

    all_results = []

    for nsps in candidate_nsps:
        baud = sr / nsps
        print(f"\n  NSPS={nsps} ({baud:.2f} Bd)")

        spectra, freqs = compute_symbol_spectra(data, sr, nsps)
        if spectra.size == 0:
            print(f"    Signal too short")
            continue

        best, per_h = sync_search(spectra, freqs, nsps, sr)

        for h, res in sorted(per_h.items()):
            snr = res["snr"]
            marker = " ***" if snr > 3.0 else ""
            print(f"    h={h:.2f}: SNR={snr:.2f}, "
                  f"t={res['t_sec']:.3f}s, f={res['f_hz']:.1f}Hz{marker}")
            all_results.append((nsps, baud, h, snr, res['t_sec'], res['f_hz'], res))

    # Sort all results by SNR
    all_results.sort(key=lambda x: x[3], reverse=True)

    print(f"\n{'='*70}")
    print(f"  TOP RESULTS (sorted by sync SNR)")
    print(f"{'='*70}")
    print(f"  {'NSPS':<6} {'Baud':<8} {'h':<6} {'SNR':<8} {'Time':<8} {'Freq':<8}")
    print(f"  {'-'*50}")
    for nsps, baud, h, snr, t, f, _ in all_results[:15]:
        print(f"  {nsps:<6} {baud:<8.2f} {h:<6.2f} {snr:<8.2f} {t:<8.3f} {f:<8.1f}")

    # Plot sync scores for the best result
    if all_results:
        best = all_results[0]
        nsps, baud, h, snr, t_best, f_best, res = best
        scores = res["scores"]

        fig, axes = plt.subplots(2, 1, figsize=(16, 10))

        # Sync score vs time (max over frequency)
        ax = axes[0]
        step = nsps // 4
        max_over_freq = scores.max(axis=1)
        t_axis = np.arange(len(max_over_freq)) * step / sr
        ax.plot(t_axis, max_over_freq)
        ax.axvline(t_best, color="red", linestyle="--",
                  label=f"Best: t={t_best:.3f}s")
        ax.set_xlabel("Time (s)")
        ax.set_ylabel("Sync Score")
        ax.set_title(f"Costas Sync Score vs Time\n"
                     f"Best: NSPS={nsps} ({baud:.2f} Bd), h={h}, SNR={snr:.2f}")
        ax.legend()
        ax.grid(True, alpha=0.3)

        # Sync score vs frequency at best time
        ax = axes[1]
        best_t_idx = res["t_step"]
        freqs_axis = np.arange(scores.shape[1]) * sr / (2 * nsps)
        ax.plot(freqs_axis, scores[best_t_idx])
        ax.axvline(f_best, color="red", linestyle="--",
                  label=f"Best: f={f_best:.1f} Hz")
        ax.set_xlabel("Frequency (Hz)")
        ax.set_ylabel("Sync Score")
        ax.set_title("Sync Score vs Frequency (at best time)")
        ax.legend()
        ax.grid(True, alpha=0.3)

        plt.tight_layout()
        plt.savefig("output/sync_detection.png", dpi=150)
        plt.close()
        print(f"\n  Saved: output/sync_detection.png")


if __name__ == "__main__":
    main()
