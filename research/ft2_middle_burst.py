#!/usr/bin/env python3
"""
FT2 Middle Burst Analysis — focused on the cleanest signal.

The middle TX in the capture should be the most complete.
Based on 3.8s T/R cycle timing, it should span roughly 5.0-9.0s.

This script:
1. Uses a wide window (4.5-9.5s) to capture full TX including ramp up/down
2. Low power threshold to catch signal edges
3. Precisely measures TX boundaries
4. Extracts complete tone sequence
5. Tries multiple NSPS values to find exact symbol alignment
6. Displays the complete decoded tone sequence for manual inspection

Usage:
    uv run python ft2_middle_burst.py ft2_capture.wav
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
    if peak_bin <= 0 or peak_bin >= len(spectrum) - 1:
        return peak_bin * df
    alpha = np.log(spectrum[peak_bin - 1] + 1e-30)
    beta = np.log(spectrum[peak_bin] + 1e-30)
    gamma = np.log(spectrum[peak_bin + 1] + 1e-30)
    denom = alpha - 2 * beta + gamma
    if abs(denom) < 1e-20:
        return peak_bin * df
    p = 0.5 * (alpha - gamma) / denom
    return (peak_bin + p) * df


def find_tx_boundaries(data, sr, fmin, fmax, t_start, t_end, threshold_db_below=-15):
    """Find precise TX boundaries using rolling power in signal band."""
    from scipy.signal import butter, filtfilt

    i0 = int(t_start * sr)
    i1 = int(t_end * sr)
    seg = data[i0:i1]

    # Bandpass filter
    nyq = sr / 2
    b, a = butter(5, [fmin / nyq, fmax / nyq], btype='band')
    filtered = filtfilt(b, a, seg)

    # Very fine power measurement (1ms frames)
    frame_ms = 1
    frame_samples = max(int(sr * frame_ms / 1000), 1)
    n_frames = len(filtered) // frame_samples
    power = np.array([
        np.mean(filtered[i * frame_samples:(i + 1) * frame_samples] ** 2)
        for i in range(n_frames)
    ])
    times = np.arange(n_frames) * frame_ms / 1000 + t_start

    # Smooth with 20ms window
    kernel_size = max(int(20 / frame_ms), 1)
    kernel = np.ones(kernel_size) / kernel_size
    smooth_power = np.convolve(power, kernel, mode='same')
    power_db = 10 * np.log10(smooth_power + 1e-20)

    peak_power = np.max(power_db)
    threshold = peak_power + threshold_db_below

    active = power_db > threshold

    # Find first and last active frame
    active_idx = np.where(active)[0]
    if len(active_idx) == 0:
        return t_start, t_end, times, power_db

    tx_on = times[active_idx[0]]
    tx_off = times[active_idx[-1]]

    return tx_on, tx_off, times, power_db


def extract_symbols(data, sr, nsps, fmin, fmax, t_start, t_end, h):
    """Extract symbol information for the given parameters."""
    i0 = int(t_start * sr)
    i1 = int(t_end * sr)
    seg = data[i0:i1]
    nfft = nsps * 8
    df = sr / nfft
    baud = sr / nsps
    spacing = h * baud
    n_sym = len(seg) // nsps

    freqs = np.zeros(n_sym)
    powers = np.zeros(n_sym)
    snrs = np.zeros(n_sym)

    for s in range(n_sym):
        chunk = seg[s * nsps:(s + 1) * nsps]
        w = np.hanning(len(chunk))
        spec = np.abs(fft(chunk * w, n=nfft))[:nfft // 2]
        f = np.arange(nfft // 2) * df
        band = (f >= fmin) & (f <= fmax)
        if not np.any(band):
            continue
        sb = spec[band]
        pk = np.argmax(sb)
        pk_abs = np.where(band)[0][pk]
        freqs[s] = parabolic_peak(spec, pk_abs, df)
        powers[s] = sb[pk]
        med = np.median(sb)
        snrs[s] = sb[pk] / (med + 1e-20)

    return freqs, powers, snrs, n_sym


def quantize_tones(freqs, powers, snrs, spacing, snr_threshold=2.0):
    """Quantize frequencies to tone indices."""
    good = snrs > snr_threshold
    if np.sum(good) < 5:
        return np.zeros(len(freqs), dtype=int), good, 1e20

    gf = freqs[good]
    best_rms = 1e20
    best_base = gf.min()
    for base in np.arange(gf.min() - spacing, gf.min() + spacing, 0.3):
        q = np.round((gf - base) / spacing) * spacing + base
        rms = np.sqrt(np.mean((gf - q) ** 2))
        if rms < best_rms:
            best_rms = rms
            best_base = base

    tone_idx = np.round((freqs - best_base) / spacing).astype(int)
    return tone_idx, good, best_rms


def scan_start_offset(data, sr, nsps, fmin, fmax, h, t_center, scan_range_ms=200):
    """
    Fine-scan the start time offset to find optimal symbol alignment.
    The idea: correct alignment maximizes the average peak-to-noise ratio
    across all symbols.
    """
    baud = sr / nsps
    spacing = h * baud
    sym_dur = nsps / sr

    offsets = np.arange(-scan_range_ms, scan_range_ms + 1, 1) / 1000.0
    scores = []

    for offset in offsets:
        t0 = t_center + offset
        t1 = t0 + 79 * sym_dur
        if t0 < 0 or t1 > len(data) / sr:
            scores.append(0)
            continue

        freqs, powers, snrs, n_sym = extract_symbols(
            data, sr, nsps, fmin, fmax, t0, t1, h)
        if n_sym < 50:
            scores.append(0)
            continue

        # Score: sum of per-symbol SNR for symbols with SNR > 2
        good_snr = snrs[snrs > 2.0]
        scores.append(np.sum(good_snr) if len(good_snr) > 0 else 0)

    scores = np.array(scores)
    best_idx = np.argmax(scores)
    best_offset = offsets[best_idx]

    return best_offset, offsets, scores


def main():
    wav_path = sys.argv[1] if len(sys.argv) > 1 else "ft2_capture.wav"
    Path("output_middle").mkdir(exist_ok=True)

    data, sr = load_wav(wav_path)
    total_dur = len(data) / sr
    print(f"Loaded: {sr} Hz, {total_dur:.2f}s")

    fmin, fmax = 1440, 1720

    # === STEP 1: Find precise TX boundaries of middle burst ===
    print(f"\n{'='*60}")
    print(f"  STEP 1: PRECISE TX BOUNDARIES (middle burst)")
    print(f"{'='*60}")

    # Scan a wide window around the expected middle burst
    tx_on, tx_off, times, power_db = find_tx_boundaries(
        data, sr, fmin, fmax, 4.0, 10.0, threshold_db_below=-20)

    tx_dur = tx_off - tx_on
    print(f"  TX on:  {tx_on:.4f}s")
    print(f"  TX off: {tx_off:.4f}s")
    print(f"  Duration: {tx_dur:.4f}s")

    for nsps in [360, 480, 576]:
        baud = sr / nsps
        n_sym = tx_dur * baud
        print(f"  NSPS={nsps} ({baud:.2f} Bd): {n_sym:.2f} symbols")

    # Also try different thresholds
    for thresh in [-10, -15, -20, -25, -30]:
        on, off, _, _ = find_tx_boundaries(
            data, sr, fmin, fmax, 4.0, 10.0, threshold_db_below=thresh)
        dur = off - on
        print(f"  Threshold {thresh} dB: {on:.3f}-{off:.3f}s ({dur:.3f}s) "
              f"= {dur*20.83:.1f} sym@576, {dur*25.0:.1f} sym@480, "
              f"{dur*33.33:.1f} sym@360")

    # Plot power envelope
    fig, ax = plt.subplots(figsize=(18, 5))
    ax.plot(times, power_db, linewidth=0.5)
    ax.axvline(tx_on, color='green', linewidth=2, label=f'TX on: {tx_on:.3f}s')
    ax.axvline(tx_off, color='red', linewidth=2, label=f'TX off: {tx_off:.3f}s')
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Power (dB)")
    ax.set_title(f"Middle Burst Power Envelope (threshold: -20dB)")
    ax.legend()
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig("output_middle/power_envelope_middle.png", dpi=150)
    plt.close()
    print(f"\n  Saved: output_middle/power_envelope_middle.png")

    # === STEP 2: Fine-tune symbol alignment ===
    print(f"\n{'='*60}")
    print(f"  STEP 2: SYMBOL ALIGNMENT OPTIMIZATION")
    print(f"{'='*60}")

    fig, axes = plt.subplots(3, 1, figsize=(16, 12))
    alignment_results = {}

    for ax, nsps in zip(axes, [360, 480, 576]):
        baud = sr / nsps
        h = 1.0
        spacing = h * baud

        best_offset, offsets, scores = scan_start_offset(
            data, sr, nsps, fmin, fmax, h, tx_on, scan_range_ms=int(1000 * nsps / sr))

        best_t_start = tx_on + best_offset
        alignment_results[nsps] = best_t_start

        print(f"\n  NSPS={nsps}: best offset = {best_offset*1000:.1f}ms "
              f"(t_start = {best_t_start:.4f}s)")

        ax.plot(offsets * 1000, scores)
        ax.axvline(best_offset * 1000, color='red', linewidth=1, linestyle='--')
        ax.set_xlabel("Offset (ms)")
        ax.set_ylabel("Sum SNR")
        ax.set_title(f"Symbol Alignment — NSPS={nsps} ({baud:.2f} Bd)")
        ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig("output_middle/alignment_scan.png", dpi=150)
    plt.close()
    print(f"  Saved: output_middle/alignment_scan.png")

    # === STEP 3: Extract complete tone sequences ===
    print(f"\n{'='*60}")
    print(f"  STEP 3: COMPLETE TONE SEQUENCES")
    print(f"{'='*60}")

    # Use wider window to capture full frame
    # At 3.8s cycle, TX could be up to 3.8s
    t_wide_start = tx_on - 0.2  # small margin
    t_wide_end = min(tx_on + 4.0, total_dur - 0.1)

    for nsps in [576, 480, 360]:
        baud = sr / nsps
        h = 1.0
        spacing = h * baud

        t_start = alignment_results.get(nsps, tx_on)
        sym_dur = nsps / sr
        # Extend to capture max possible symbols
        t_end = min(t_start + 100 * sym_dur, total_dur - 0.1)

        freqs, powers, snrs, n_sym = extract_symbols(
            data, sr, nsps, fmin, fmax, t_start, t_end, h)
        tones, good, rms = quantize_tones(freqs, powers, snrs, spacing, snr_threshold=2.0)

        n_good = int(np.sum(good))
        good_tones = tones[good]
        t_range = int(good_tones.max() - good_tones.min()) + 1 if n_good > 0 else 0

        print(f"\n  NSPS={nsps} ({baud:.2f} Bd):")
        print(f"    Start: {t_start:.4f}s, {n_sym} symbols total")
        print(f"    Good symbols: {n_good}, RMS: {rms:.2f} Hz")
        print(f"    Tone range: {t_range} tones")

        if n_good < 10:
            continue

        # Print tone sequence with quality indicators
        print(f"\n    TONE SEQUENCE (SNR>2 marked with *):")
        line = "    "
        for s in range(min(n_sym, 120)):
            if good[s]:
                line += f"{tones[s]:2d}* "
            else:
                line += f"{tones[s]:2d}  "
            if (s + 1) % 20 == 0:
                print(line)
                line = "    "
        if line.strip():
            print(line)

        # Show power/SNR per symbol
        print(f"\n    Per-symbol detail:")
        for s in range(min(n_sym, 100)):
            marker = "***" if snrs[s] > 10 else "** " if snrs[s] > 3 else "*  " if snrs[s] > 2 else "   "
            if powers[s] > 0:
                print(f"      [{s:3d}] tone={tones[s]:2d}  "
                      f"f={freqs[s]:7.1f} Hz  "
                      f"snr={snrs[s]:6.1f} {marker}")

    # === STEP 4: Try h=0.5 specifically for NSPS=576 ===
    print(f"\n{'='*60}")
    print(f"  STEP 4: NSPS=576 WITH h=0.5")
    print(f"{'='*60}")

    nsps = 576
    h = 0.5
    baud = sr / nsps
    spacing = h * baud

    t_start = alignment_results.get(nsps, tx_on)
    sym_dur = nsps / sr
    t_end = min(t_start + 100 * sym_dur, total_dur - 0.1)

    freqs, powers, snrs, n_sym = extract_symbols(
        data, sr, nsps, fmin, fmax, t_start, t_end, h)
    tones, good, rms = quantize_tones(freqs, powers, snrs, spacing, snr_threshold=2.0)

    n_good = int(np.sum(good))
    print(f"  NSPS=576, h=0.5 (spacing={spacing:.2f} Hz):")
    print(f"  {n_sym} symbols, {n_good} good, RMS={rms:.2f} Hz")
    if n_good > 0:
        good_tones = tones[good]
        print(f"  Tone range: {good_tones.min()}-{good_tones.max()} "
              f"({good_tones.max()-good_tones.min()+1} tones)")

    # === STEP 5: Symbol-aligned waterfall of middle burst ===
    print(f"\n{'='*60}")
    print(f"  STEP 5: SYMBOL-ALIGNED WATERFALL")
    print(f"{'='*60}")

    for nsps in [576, 480, 360]:
        baud = sr / nsps
        t_start = alignment_results.get(nsps, tx_on)
        sym_dur = nsps / sr
        t_end = min(t_start + 85 * sym_dur, total_dur - 0.1)

        i0 = int(t_start * sr)
        i1 = int(t_end * sr)
        seg = data[i0:i1]
        nfft = nsps * 4
        df = sr / nfft
        n_sym = len(seg) // nsps

        waterfall = np.zeros((n_sym, nfft // 2))
        wf_freqs = np.arange(nfft // 2) * df

        for s in range(n_sym):
            chunk = seg[s * nsps:(s + 1) * nsps]
            w = np.hanning(len(chunk))
            waterfall[s] = np.abs(fft(chunk * w, n=nfft))[:nfft // 2]

        mask = (wf_freqs >= fmin) & (wf_freqs <= fmax)
        wf_db = 10 * np.log10(waterfall[:, mask] + 1e-20)
        vmin, vmax = np.percentile(wf_db, [10, 99])

        fig, ax = plt.subplots(figsize=(16, max(8, n_sym * 0.15)))
        ax.imshow(wf_db, aspect='auto', origin='lower',
                  extent=[wf_freqs[mask][0], wf_freqs[mask][-1], 0, n_sym],
                  cmap='inferno', vmin=vmin, vmax=vmax)

        # Mark symbol boundaries
        for i in range(n_sym):
            ax.axhline(i, color='white', alpha=0.08, linewidth=0.3)

        # Mark 79-symbol frame boundary
        ax.axhline(79, color='cyan', alpha=0.6, linewidth=1.5,
                  label='79 sym (FT8 frame)')

        ax.set_xlabel("Frequency (Hz)")
        ax.set_ylabel("Symbol index")
        ax.set_title(f"Middle Burst Waterfall — NSPS={nsps} ({baud:.2f} Bd), "
                     f"df={df:.2f} Hz\n"
                     f"Start: {t_start:.3f}s, {n_sym} symbols")
        ax.legend(loc='upper right')
        plt.colorbar(ax.images[0], ax=ax, label='dB')
        plt.tight_layout()
        plt.savefig(f"output_middle/waterfall_middle_nsps{nsps}.png", dpi=200)
        plt.close()
        print(f"  Saved: output_middle/waterfall_middle_nsps{nsps}.png")

    # === STEP 6: Also look at first burst with wider window ===
    print(f"\n{'='*60}")
    print(f"  STEP 6: FIRST BURST (TAIL END)")
    print(f"{'='*60}")

    # First burst ends around 1.4s — check how far back signal goes
    tx_on1, tx_off1, _, _ = find_tx_boundaries(
        data, sr, fmin, fmax, 0.0, 3.0, threshold_db_below=-20)
    print(f"  First burst: {tx_on1:.3f}-{tx_off1:.3f}s ({tx_off1-tx_on1:.3f}s)")

    for nsps in [576, 480, 360]:
        baud = sr / nsps
        n_sym = (tx_off1 - tx_on1) * baud
        print(f"    NSPS={nsps}: {n_sym:.1f} symbols")

    # === STEP 7: Third burst ===
    print(f"\n{'='*60}")
    print(f"  STEP 7: THIRD BURST")
    print(f"{'='*60}")

    tx_on3, tx_off3, _, _ = find_tx_boundaries(
        data, sr, fmin, fmax, 11.0, total_dur, threshold_db_below=-20)
    print(f"  Third burst: {tx_on3:.3f}-{tx_off3:.3f}s ({tx_off3-tx_on3:.3f}s)")

    for nsps in [576, 480, 360]:
        baud = sr / nsps
        n_sym = (tx_off3 - tx_on3) * baud
        print(f"    NSPS={nsps}: {n_sym:.1f} symbols")

    print(f"\n{'='*60}")
    print(f"  DONE — check output_middle/")
    print(f"{'='*60}")


if __name__ == "__main__":
    main()
