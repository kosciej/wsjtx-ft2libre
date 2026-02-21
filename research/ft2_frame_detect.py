#!/usr/bin/env python3
"""
FT2 Frame Structure Detection — find frame boundaries without assuming
a specific sync pattern.

Approaches:
1. Spectrogram autocorrelation — find periodicity matching frame length
2. Power envelope analysis — TX on/off transitions mark frame boundaries
3. Instantaneous frequency analysis — measure exact baud rate
4. Cross-burst frequency alignment — same base freq in different bursts?
5. Direct symbol count — count symbols from TX start to TX end

Usage:
    uv run python ft2_frame_detect.py ft2_capture.wav
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


def bandpass_filter(data, sr, fmin, fmax, order=5):
    """Apply bandpass filter."""
    nyq = sr / 2
    b, a = sig.butter(order, [fmin / nyq, fmax / nyq], btype='band')
    return sig.filtfilt(b, a, data)


def measure_baud_rate(data, sr, fmin, fmax, t_start, t_end):
    """Measure exact baud rate using autocorrelation of instantaneous frequency."""
    i0 = int(t_start * sr)
    i1 = int(t_end * sr)
    segment = data[i0:i1]

    # Bandpass filter
    filtered = bandpass_filter(segment, sr, fmin, fmax)

    # Analytic signal
    analytic = sig.hilbert(filtered)
    inst_phase = np.unwrap(np.angle(analytic))
    inst_freq = np.diff(inst_phase) * sr / (2 * np.pi)

    # Autocorrelation of instantaneous frequency
    # Look for the period of tone changes = symbol duration
    if_centered = inst_freq - inst_freq.mean()
    max_lag = int(sr * 0.1)  # up to 100ms
    min_lag = int(sr * 0.02)  # from 20ms

    autocorr = np.correlate(if_centered[:min(len(if_centered), sr * 2)],
                            if_centered[:min(len(if_centered), sr * 2)],
                            mode='full')
    autocorr = autocorr[len(autocorr) // 2:]  # positive lags
    autocorr = autocorr / autocorr[0]

    # Find first trough then first peak (symbol period)
    search = autocorr[min_lag:max_lag]
    lags = np.arange(min_lag, min_lag + len(search))

    # Find peaks
    peaks, props = sig.find_peaks(search, height=0.05, distance=50)
    troughs, _ = sig.find_peaks(-search, height=0.05, distance=50)

    print(f"\n  Instantaneous frequency autocorrelation:")
    if len(peaks) > 0:
        for p in peaks[:5]:
            lag_samples = lags[p]
            lag_ms = lag_samples / sr * 1000
            baud = sr / lag_samples
            nsps = lag_samples
            print(f"    Peak at lag={lag_samples} samples ({lag_ms:.1f} ms) "
                  f"→ baud={baud:.2f}, NSPS={nsps}")
    else:
        print("    No clear peaks found")

    return autocorr, lags if len(peaks) > 0 else (autocorr, None)


def measure_tx_duration(data, sr, fmin, fmax):
    """Precisely measure TX on/off times from power envelope."""
    filtered = bandpass_filter(data, sr, fmin, fmax)

    # Short-time power
    frame_ms = 5
    frame_samples = int(sr * frame_ms / 1000)
    n_frames = len(filtered) // frame_samples
    power = np.array([
        np.mean(filtered[i * frame_samples:(i + 1) * frame_samples] ** 2)
        for i in range(n_frames)
    ])
    times = np.arange(n_frames) * frame_ms / 1000

    # Smooth
    kernel = np.ones(5) / 5
    smooth_power = np.convolve(power, kernel, mode='same')
    power_db = 10 * np.log10(smooth_power + 1e-20)

    # Find TX transitions (sharp rise/fall in power)
    threshold_db = np.max(power_db) - 10  # 10 dB below peak
    active = power_db > threshold_db

    # Find edges
    transitions = np.diff(active.astype(int))
    rises = np.where(transitions == 1)[0]
    falls = np.where(transitions == -1)[0]

    print(f"\n  TX Power Envelope Analysis:")
    print(f"    Peak power: {np.max(power_db):.1f} dB")
    print(f"    Threshold: {threshold_db:.1f} dB")

    tx_periods = []
    if len(rises) > 0:
        print(f"    TX on times:  {[f'{times[r]:.3f}' for r in rises]}")
    if len(falls) > 0:
        print(f"    TX off times: {[f'{times[f]:.3f}' for f in falls]}")

    # Match rises to falls
    for r in rises:
        next_falls = falls[falls > r]
        if len(next_falls) > 0:
            f = next_falls[0]
            duration = times[f] - times[r]
            n_sym_576 = duration * 20.83
            n_sym_480 = duration * 25.0
            n_sym_360 = duration * 33.33
            tx_periods.append((times[r], times[f], duration))
            print(f"    TX period: {times[r]:.3f}-{times[f]:.3f}s "
                  f"({duration:.3f}s = {n_sym_576:.1f} sym@576, "
                  f"{n_sym_480:.1f} sym@480, {n_sym_360:.1f} sym@360)")

    return times, power_db, tx_periods


def spectrogram_autocorrelation(data, sr, nsps, fmin, fmax, t_start, t_end):
    """
    Compute autocorrelation of the spectrogram along time axis.
    Periodicity at frame_length × symbol_duration indicates frame structure.
    """
    i0 = int(t_start * sr)
    i1 = int(t_end * sr)
    seg = data[i0:i1]

    nfft = nsps * 2
    df = sr / nfft
    n_sym = len(seg) // nsps

    # Symbol-aligned spectrogram
    spectra = np.zeros((n_sym, nfft // 2))
    for s in range(n_sym):
        chunk = seg[s * nsps:(s + 1) * nsps]
        w = np.hanning(nsps)
        spectra[s] = np.abs(fft(chunk * w, n=nfft))[:nfft // 2] ** 2

    # Restrict to signal band
    freqs = np.arange(nfft // 2) * df
    band = (freqs >= fmin) & (freqs <= fmax)
    spec_band = spectra[:, band]

    # Autocorrelation along time for each frequency bin, then average
    max_lag = min(n_sym - 1, 120)
    autocorr = np.zeros(max_lag)

    for fb in range(spec_band.shape[1]):
        col = spec_band[:, fb]
        col_centered = col - col.mean()
        for lag in range(max_lag):
            autocorr[lag] += np.sum(col_centered[:n_sym - lag] * col_centered[lag:n_sym])

    autocorr /= autocorr[0] + 1e-20

    return autocorr


def detailed_tone_trace(data, sr, nsps, fmin, fmax, t_start, t_end, h):
    """
    Extract detailed per-symbol information: frequency, power, SNR.
    """
    i0 = int(t_start * sr)
    i1 = int(t_end * sr)
    seg = data[i0:i1]
    nfft = nsps * 8
    df = sr / nfft
    baud = sr / nsps
    spacing = h * baud
    n_sym = len(seg) // nsps

    result = {
        'freqs': np.zeros(n_sym),
        'powers': np.zeros(n_sym),
        'snrs': np.zeros(n_sym),
        'n_sym': n_sym,
    }

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

        # Parabolic interpolation
        if 0 < pk_abs < len(spec) - 1:
            alpha = np.log(spec[pk_abs - 1] + 1e-30)
            beta = np.log(spec[pk_abs] + 1e-30)
            gamma = np.log(spec[pk_abs + 1] + 1e-30)
            denom = alpha - 2 * beta + gamma
            if abs(denom) > 1e-20:
                p = 0.5 * (alpha - gamma) / denom
                result['freqs'][s] = (pk_abs + p) * df
            else:
                result['freqs'][s] = pk_abs * df
        else:
            result['freqs'][s] = pk_abs * df

        result['powers'][s] = sb[pk]
        # SNR = peak / median of band
        result['snrs'][s] = sb[pk] / (np.median(sb) + 1e-20)

    return result


def main():
    wav_path = sys.argv[1] if len(sys.argv) > 1 else "ft2_capture.wav"
    Path("output_frame").mkdir(exist_ok=True)

    data, sr = load_wav(wav_path)
    print(f"Loaded: {sr} Hz, {len(data)/sr:.2f}s")

    fmin, fmax = 1440, 1720

    # === STEP 1: Precise TX duration measurement ===
    print(f"\n{'='*60}")
    print(f"  STEP 1: TX DURATION MEASUREMENT")
    print(f"{'='*60}")

    times, power_db, tx_periods = measure_tx_duration(data, sr, fmin, fmax)

    # Plot power envelope
    fig, ax = plt.subplots(figsize=(18, 5))
    ax.plot(times, power_db)
    for t0, t1, dur in tx_periods:
        ax.axvspan(t0, t1, alpha=0.2, color='green')
        ax.text((t0 + t1) / 2, np.max(power_db) - 3, f"{dur:.3f}s",
                ha='center', fontsize=9)
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Power (dB)")
    ax.set_title("TX Power Envelope (bandpass filtered)")
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig("output_frame/power_envelope.png", dpi=150)
    plt.close()
    print(f"  Saved: output_frame/power_envelope.png")

    # === STEP 2: Baud rate from instantaneous frequency ===
    print(f"\n{'='*60}")
    print(f"  STEP 2: BAUD RATE FROM INST. FREQ AUTOCORRELATION")
    print(f"{'='*60}")

    # Use the strongest TX period
    if tx_periods:
        strongest = max(tx_periods, key=lambda x: x[2])
        t0, t1, dur = strongest
        print(f"  Using TX period: {t0:.3f}-{t1:.3f}s ({dur:.3f}s)")
    else:
        t0, t1 = 0.5, 5.0
        print(f"  No TX periods detected, using {t0}-{t1}s")

    autocorr, lags = measure_baud_rate(data, sr, fmin, fmax, t0, t1)

    # Plot
    fig, ax = plt.subplots(figsize=(16, 5))
    ax.plot(autocorr[:int(sr * 0.08)])
    for nsps_candidate in [360, 480, 576]:
        ax.axvline(nsps_candidate, color='red', linestyle='--', alpha=0.5,
                  label=f'NSPS={nsps_candidate} ({sr/nsps_candidate:.2f} Bd)')
    ax.set_xlabel("Lag (samples)")
    ax.set_ylabel("Autocorrelation")
    ax.set_title("Instantaneous Frequency Autocorrelation")
    ax.legend()
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig("output_frame/baud_autocorrelation.png", dpi=150)
    plt.close()
    print(f"  Saved: output_frame/baud_autocorrelation.png")

    # === STEP 3: Symbol count from TX duration ===
    print(f"\n{'='*60}")
    print(f"  STEP 3: SYMBOL COUNT FROM TX DURATION")
    print(f"{'='*60}")

    if tx_periods:
        for t0, t1, dur in tx_periods:
            print(f"\n  TX: {t0:.3f}-{t1:.3f}s ({dur:.3f}s)")
            for nsps in [360, 480, 576]:
                baud = sr / nsps
                n_sym = dur * baud
                print(f"    NSPS={nsps} ({baud:.2f} Bd): {n_sym:.1f} symbols")
                # Check if close to known frame lengths
                for frame_len in [77, 79, 85, 87, 103, 105]:
                    n_frames = n_sym / frame_len
                    if abs(n_frames - round(n_frames)) < 0.15:
                        print(f"      → {round(n_frames)} × {frame_len} symbols "
                              f"(residual: {n_sym - round(n_frames) * frame_len:.1f})")

    # === STEP 4: Spectrogram autocorrelation ===
    print(f"\n{'='*60}")
    print(f"  STEP 4: SPECTROGRAM AUTOCORRELATION")
    print(f"{'='*60}")

    if tx_periods:
        strongest = max(tx_periods, key=lambda x: x[2])
        t0, t1, dur = strongest
    else:
        t0, t1 = 0.2, 5.5

    fig, axes = plt.subplots(3, 1, figsize=(16, 12))
    for ax, nsps in zip(axes, [360, 480, 576]):
        baud = sr / nsps
        autocorr = spectrogram_autocorrelation(data, sr, nsps, fmin, fmax, t0, t1)

        ax.plot(autocorr)
        # Mark known frame lengths
        for fl, color in [(77, 'orange'), (79, 'red'), (85, 'green'),
                          (87, 'blue'), (103, 'purple')]:
            if fl < len(autocorr):
                ax.axvline(fl, color=color, linestyle='--', alpha=0.5,
                          label=f'{fl} sym')
        ax.set_xlabel("Lag (symbols)")
        ax.set_ylabel("Autocorrelation")
        ax.set_title(f"Spectrogram Autocorrelation — NSPS={nsps} ({baud:.2f} Bd)")
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig("output_frame/spec_autocorrelation.png", dpi=150)
    plt.close()
    print(f"  Saved: output_frame/spec_autocorrelation.png")

    # === STEP 5: Detailed tone trace with SNR ===
    print(f"\n{'='*60}")
    print(f"  STEP 5: DETAILED TONE TRACE (NSPS=576, h=1.0)")
    print(f"{'='*60}")

    nsps = 576
    h = 1.0
    baud = sr / nsps
    spacing = h * baud

    if tx_periods:
        for ti, (t0, t1, dur) in enumerate(tx_periods):
            result = detailed_tone_trace(data, sr, nsps, fmin, fmax, t0, t1, h)

            # Find high-SNR symbols
            high_snr = result['snrs'] > np.percentile(result['snrs'], 40)
            good_freqs = result['freqs'][high_snr]

            if len(good_freqs) > 5:
                # Quantize
                best_rms = 1e20
                best_base = good_freqs.min()
                for base in np.arange(good_freqs.min() - spacing,
                                      good_freqs.min() + spacing, 0.5):
                    q = np.round((good_freqs - base) / spacing) * spacing + base
                    rms = np.sqrt(np.mean((good_freqs - q) ** 2))
                    if rms < best_rms:
                        best_rms = rms
                        best_base = base

                tone_idx = np.round((result['freqs'] - best_base) / spacing).astype(int)
                good_tones = tone_idx[high_snr]

                print(f"\n  TX{ti} ({t0:.3f}-{t1:.3f}s, {dur:.3f}s):")
                print(f"    {result['n_sym']} symbols, {int(np.sum(high_snr))} high-SNR")
                print(f"    RMS quantization error: {best_rms:.2f} Hz")
                print(f"    Base frequency: {best_base:.2f} Hz")
                print(f"    Tone range: {int(good_tones.min())}-{int(good_tones.max())} "
                      f"({int(good_tones.max() - good_tones.min()) + 1} tones)")

                # Print full tone sequence with power/SNR annotation
                print(f"\n    Symbol-by-symbol (tone, power, snr):")
                for s in range(result['n_sym']):
                    marker = "*" if high_snr[s] else " "
                    if result['powers'][s] > 0:
                        print(f"      sym {s:3d}: tone={tone_idx[s]:2d}  "
                              f"f={result['freqs'][s]:.1f}  "
                              f"pwr={result['powers'][s]:.1f}  "
                              f"snr={result['snrs'][s]:.1f} {marker}")

    # === STEP 6: Visual comparison of tone sequences ===
    print(f"\n{'='*60}")
    print(f"  STEP 6: VISUAL TONE SEQUENCE COMPARISON")
    print(f"{'='*60}")

    if tx_periods and len(tx_periods) >= 2:
        fig, axes = plt.subplots(len(tx_periods), 1,
                                  figsize=(20, 4 * len(tx_periods)))
        if len(tx_periods) == 1:
            axes = [axes]

        for ax, (t0, t1, dur) in zip(axes, tx_periods):
            result = detailed_tone_trace(data, sr, nsps, fmin, fmax, t0, t1, h)
            high_snr = result['snrs'] > np.percentile(result['snrs'], 40)

            if np.sum(high_snr) > 5:
                good_freqs = result['freqs'][high_snr]
                best_base = good_freqs.min()
                for base in np.arange(good_freqs.min() - spacing,
                                      good_freqs.min() + spacing, 0.5):
                    q = np.round((good_freqs - base) / spacing) * spacing + base
                    rms = np.sqrt(np.mean((good_freqs - q) ** 2))
                    if rms < 1e20:
                        best_base = base

                tone_idx = np.round((result['freqs'] - best_base) / spacing).astype(int)

                x = np.arange(result['n_sym'])
                colors = np.where(high_snr, result['snrs'], 0)

                scatter = ax.scatter(x[high_snr], tone_idx[high_snr],
                                     c=result['snrs'][high_snr],
                                     cmap='viridis', s=20, edgecolors='none')
                ax.scatter(x[~high_snr], tone_idx[~high_snr],
                          c='lightgray', s=5, alpha=0.3)

                # Mark frame boundaries at 79 symbols
                for i in range(0, result['n_sym'], 79):
                    ax.axvline(i, color='red', alpha=0.4, linewidth=0.8)

                ax.set_ylabel("Tone index")
                ax.set_title(f"TX {t0:.3f}-{t1:.3f}s ({dur:.3f}s, "
                            f"{result['n_sym']} sym)")
                ax.grid(True, alpha=0.2)
                ax.set_ylim(-2, 9)

                plt.colorbar(scatter, ax=ax, label="SNR")

        axes[-1].set_xlabel("Symbol index")
        plt.suptitle(f"Tone Sequences — NSPS={nsps}, h={h}", fontsize=13)
        plt.tight_layout()
        plt.savefig("output_frame/tone_sequences_detailed.png", dpi=150)
        plt.close()
        print(f"  Saved: output_frame/tone_sequences_detailed.png")

    print(f"\n{'='*60}")
    print(f"  DONE — check output_frame/")
    print(f"{'='*60}")


if __name__ == "__main__":
    main()
