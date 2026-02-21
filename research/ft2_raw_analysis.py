#!/usr/bin/env python3
"""
Fresh analysis of FT2 signal from raw spectrogram.
No assumptions about NSPS, h, symbol count, or frame structure.
Measure everything directly from the signal.
"""

import numpy as np
import soundfile as sf
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.signal import spectrogram, butter, filtfilt

OUTPUT_DIR = Path("output_raw")
OUTPUT_DIR.mkdir(exist_ok=True)

FS = 12000


def load_middle_burst(filename, skip_first=False):
    """Load signal and extract middle burst."""
    sig, _ = sf.read(filename)
    if len(sig.shape) > 1:
        sig = sig[:, 0]

    # Power envelope
    wl = int(0.05 * FS)
    pw = np.convolve(sig**2, np.ones(wl)/wl, mode='same')
    th = np.max(pw) * 0.01
    ab = pw > th
    tr = np.diff(ab.astype(int))
    starts = np.where(tr == 1)[0]
    ends = np.where(tr == -1)[0]
    if ab[0]: starts = np.concatenate([[0], starts])
    if ab[-1]: ends = np.concatenate([ends, [len(sig)-1]])
    n = min(len(starts), len(ends))
    bursts = [(starts[i], ends[i]) for i in range(n)
              if (ends[i]-starts[i])/FS > 1.5]

    if skip_first and len(bursts) > 1:
        s, e = bursts[1]
    elif bursts:
        s, e = max(bursts, key=lambda b: b[1]-b[0])
    else:
        return None, 0, 0

    return sig[s:e], s, e


def main():
    # ================================================================
    # Load both bursts
    # ================================================================
    burst1, s1, e1 = load_middle_burst("ft2_capture.wav")
    burst2, s2, e2 = load_middle_burst("ft2_capture2.wav", skip_first=True)

    for label, burst, s, e in [("Capture 1", burst1, s1, e1),
                                 ("Capture 2", burst2, s2, e2)]:
        dur = len(burst) / FS
        print(f"{label}: {s/FS:.3f}-{e/FS:.3f}s, duration={dur:.3f}s, "
              f"samples={len(burst)}")

    # ================================================================
    # STEP 1: High-resolution spectrogram
    # ================================================================
    print("\n=== STEP 1: HIGH-RESOLUTION SPECTROGRAM ===")

    for label, burst in [("Capture 1", burst1), ("Capture 2", burst2)]:
        fig, axes = plt.subplots(3, 1, figsize=(20, 16))

        # Ultra-high frequency resolution (long window)
        nfft_hires = 4096  # df = 12000/4096 = 2.93 Hz
        hop = 32  # ~2.67ms time resolution
        f, t, Sxx = spectrogram(burst, fs=FS, nperseg=nfft_hires,
                                noverlap=nfft_hires-hop, nfft=nfft_hires)
        f_mask = (f >= 1000) & (f <= 2000)
        Sxx_db = 10 * np.log10(Sxx[f_mask] + 1e-30)

        axes[0].pcolormesh(t, f[f_mask], Sxx_db, shading='gouraud', cmap='viridis')
        axes[0].set_ylabel("Frequency (Hz)")
        axes[0].set_title(f"{label} - High freq resolution (df={FS/nfft_hires:.1f}Hz)")

        # High time resolution (short window)
        nfft_time = 512  # df = 23.4 Hz, dt = 42.7ms
        hop2 = 16
        f2, t2, Sxx2 = spectrogram(burst, fs=FS, nperseg=nfft_time,
                                     noverlap=nfft_time-hop2, nfft=nfft_time*2)
        f2_mask = (f2 >= 1000) & (f2 <= 2000)
        Sxx2_db = 10 * np.log10(Sxx2[f2_mask] + 1e-30)

        axes[1].pcolormesh(t2, f2[f2_mask], Sxx2_db, shading='gouraud', cmap='viridis')
        axes[1].set_ylabel("Frequency (Hz)")
        axes[1].set_title(f"{label} - High time resolution (dt={nfft_time/FS*1000:.1f}ms)")

        # Balanced resolution
        nfft_bal = 1024
        hop3 = 64
        f3, t3, Sxx3 = spectrogram(burst, fs=FS, nperseg=nfft_bal,
                                     noverlap=nfft_bal-hop3, nfft=nfft_bal*2)
        f3_mask = (f3 >= 1000) & (f3 <= 2000)
        Sxx3_db = 10 * np.log10(Sxx3[f3_mask] + 1e-30)

        axes[2].pcolormesh(t3, f3[f3_mask], Sxx3_db, shading='gouraud', cmap='viridis')
        axes[2].set_ylabel("Frequency (Hz)")
        axes[2].set_xlabel("Time (s)")
        axes[2].set_title(f"{label} - Balanced resolution")

        plt.tight_layout()
        plt.savefig(OUTPUT_DIR / f"spectrogram_{label.replace(' ', '_').lower()}.png", dpi=200)
        plt.close()

    # ================================================================
    # STEP 2: Instantaneous frequency tracking
    # ================================================================
    print("\n=== STEP 2: INSTANTANEOUS FREQUENCY ===")

    for label, burst in [("Capture 1", burst1), ("Capture 2", burst2)]:
        # Bandpass filter around signal
        nyq = FS / 2
        b, a = butter(4, [1300/nyq, 1700/nyq], btype='band')
        filtered = filtfilt(b, a, burst)

        # Analytic signal via Hilbert transform
        from scipy.signal import hilbert
        analytic = hilbert(filtered)
        inst_phase = np.unwrap(np.angle(analytic))
        inst_freq = np.diff(inst_phase) / (2 * np.pi) * FS

        # Smooth instantaneous frequency
        smooth_win = 64
        inst_freq_smooth = np.convolve(inst_freq,
                                        np.ones(smooth_win)/smooth_win,
                                        mode='same')

        t_if = np.arange(len(inst_freq)) / FS

        fig, axes = plt.subplots(2, 1, figsize=(20, 10))
        axes[0].plot(t_if, inst_freq, linewidth=0.2, alpha=0.3, label='Raw')
        axes[0].plot(t_if, inst_freq_smooth, linewidth=0.8, color='red',
                     label=f'Smoothed ({smooth_win} samples)')
        axes[0].set_ylabel("Frequency (Hz)")
        axes[0].set_title(f"{label} - Instantaneous Frequency")
        axes[0].legend()
        axes[0].set_ylim(1300, 1700)
        axes[0].grid(True, alpha=0.3)

        # Histogram of instantaneous frequency -> reveals tone grid
        valid_if = inst_freq_smooth[(inst_freq_smooth > 1400) &
                                     (inst_freq_smooth < 1700)]
        axes[1].hist(valid_if, bins=500, density=True, alpha=0.7)
        axes[1].set_xlabel("Frequency (Hz)")
        axes[1].set_ylabel("Density")
        axes[1].set_title(f"{label} - Frequency Histogram (reveals tone grid)")
        axes[1].grid(True, alpha=0.3)

        plt.tight_layout()
        plt.savefig(OUTPUT_DIR / f"inst_freq_{label.replace(' ', '_').lower()}.png",
                    dpi=150)
        plt.close()

        # Find peaks in histogram -> tone frequencies
        from scipy.signal import find_peaks
        hist_vals, hist_bins = np.histogram(valid_if, bins=1000)
        hist_centers = (hist_bins[:-1] + hist_bins[1:]) / 2
        # Smooth histogram
        hist_smooth = np.convolve(hist_vals, np.ones(5)/5, mode='same')
        peaks, props = find_peaks(hist_smooth, height=np.max(hist_smooth)*0.1,
                                   distance=20)
        tone_freqs = hist_centers[peaks]
        tone_heights = hist_smooth[peaks]

        # Sort by frequency
        order = np.argsort(tone_freqs)
        tone_freqs = tone_freqs[order]
        tone_heights = tone_heights[order]

        print(f"\n{label} - Detected tone frequencies:")
        for i, (f, h) in enumerate(zip(tone_freqs, tone_heights)):
            print(f"  Tone {i}: {f:.1f} Hz (height={h:.0f})")

        if len(tone_freqs) >= 2:
            spacings = np.diff(tone_freqs)
            print(f"  Spacings: {np.round(spacings, 1).tolist()}")
            print(f"  Mean spacing: {np.mean(spacings):.2f} Hz")
            print(f"  Std spacing: {np.std(spacings):.2f} Hz")
            if np.mean(spacings) > 0:
                h_est = np.mean(spacings) / (FS / 360)
                print(f"  Implied h (if NSPS=360): {h_est:.4f}")
                baud_est = np.mean(spacings) / 0.5
                nsps_est = FS / baud_est
                print(f"  Implied baud (if h=0.5): {baud_est:.2f}, NSPS={nsps_est:.1f}")

    # ================================================================
    # STEP 3: Measure baud rate from frequency transitions
    # ================================================================
    print("\n=== STEP 3: BAUD RATE FROM TRANSITIONS ===")

    for label, burst in [("Capture 1", burst1), ("Capture 2", burst2)]:
        nyq = FS / 2
        b, a = butter(4, [1300/nyq, 1700/nyq], btype='band')
        filtered = filtfilt(b, a, burst)

        from scipy.signal import hilbert
        analytic = hilbert(filtered)
        inst_phase = np.unwrap(np.angle(analytic))
        inst_freq = np.diff(inst_phase) / (2 * np.pi) * FS

        smooth_win = 32
        inst_freq_smooth = np.convolve(inst_freq,
                                        np.ones(smooth_win)/smooth_win,
                                        mode='same')

        # Derivative of frequency -> large spikes at symbol transitions
        freq_deriv = np.abs(np.diff(inst_freq_smooth))

        # Autocorrelation of frequency derivative -> peaks at symbol period
        max_lag = int(0.1 * FS)  # up to 100ms
        min_lag = int(0.02 * FS)  # from 20ms
        acf = np.correlate(freq_deriv[:int(1.5*FS)],
                           freq_deriv[:int(1.5*FS)], mode='full')
        acf = acf[len(acf)//2:]  # positive lags only
        acf = acf / acf[0]  # normalize

        # Find first peak after min_lag
        from scipy.signal import find_peaks
        acf_peaks, _ = find_peaks(acf[min_lag:max_lag], height=0.1)
        acf_peaks += min_lag

        lags = np.arange(len(acf)) / FS * 1000  # in ms

        fig, axes = plt.subplots(3, 1, figsize=(20, 12))

        t_d = np.arange(len(freq_deriv)) / FS
        axes[0].plot(t_d[:int(0.5*FS)], freq_deriv[:int(0.5*FS)],
                     linewidth=0.3)
        axes[0].set_ylabel("|df/dt|")
        axes[0].set_title(f"{label} - Frequency derivative (first 0.5s)")
        axes[0].set_xlabel("Time (s)")

        axes[1].plot(lags[min_lag:max_lag], acf[min_lag:max_lag], 'b-')
        for pk in acf_peaks:
            axes[1].axvline(pk/FS*1000, color='r', alpha=0.5)
            axes[1].text(pk/FS*1000, acf[pk], f'{pk/FS*1000:.1f}ms',
                        fontsize=8, ha='center', va='bottom')
        axes[1].set_xlabel("Lag (ms)")
        axes[1].set_ylabel("Autocorrelation")
        axes[1].set_title(f"{label} - ACF of freq derivative (symbol period)")
        axes[1].grid(True, alpha=0.3)

        # Also try direct frequency ACF
        valid_range = inst_freq_smooth[int(0.1*FS):int(1.5*FS)]
        acf_freq = np.correlate(valid_range - np.mean(valid_range),
                                valid_range - np.mean(valid_range), mode='full')
        acf_freq = acf_freq[len(acf_freq)//2:]
        acf_freq = acf_freq / acf_freq[0]

        axes[2].plot(lags[min_lag:max_lag], acf_freq[min_lag:max_lag], 'g-')
        axes[2].set_xlabel("Lag (ms)")
        axes[2].set_ylabel("Autocorrelation")
        axes[2].set_title(f"{label} - ACF of instantaneous frequency")
        axes[2].grid(True, alpha=0.3)

        plt.tight_layout()
        plt.savefig(OUTPUT_DIR / f"baud_rate_{label.replace(' ', '_').lower()}.png",
                    dpi=150)
        plt.close()

        if len(acf_peaks) > 0:
            first_peak = acf_peaks[0]
            period_ms = first_peak / FS * 1000
            baud = FS / first_peak
            nsps = first_peak
            print(f"\n{label}:")
            print(f"  First ACF peak at lag {first_peak} samples = {period_ms:.2f} ms")
            print(f"  Implied baud rate: {baud:.2f} Hz")
            print(f"  Implied NSPS: {nsps}")
        else:
            print(f"\n{label}: No clear ACF peak found")

    # ================================================================
    # STEP 4: Average spectrum -> frequency grid
    # ================================================================
    print("\n=== STEP 4: AVERAGE SPECTRUM ===")

    for label, burst in [("Capture 1", burst1), ("Capture 2", burst2)]:
        # Compute average power spectrum using many overlapping windows
        window_len = 2048
        hop = 128
        n_windows = (len(burst) - window_len) // hop
        avg_spec = np.zeros(window_len // 2 + 1)

        for i in range(n_windows):
            start = i * hop
            seg = burst[start:start+window_len] * np.hanning(window_len)
            spec = np.abs(np.fft.rfft(seg))**2
            avg_spec += spec

        avg_spec /= n_windows
        freq_axis = np.fft.rfftfreq(window_len, 1/FS)

        fig, ax = plt.subplots(figsize=(16, 6))
        f_mask = (freq_axis >= 1400) & (freq_axis <= 1700)
        ax.plot(freq_axis[f_mask], 10*np.log10(avg_spec[f_mask] + 1e-30),
                linewidth=0.8)
        ax.set_xlabel("Frequency (Hz)")
        ax.set_ylabel("Power (dB)")
        ax.set_title(f"{label} - Average spectrum (reveals tone grid)")
        ax.grid(True, alpha=0.3)

        plt.tight_layout()
        plt.savefig(OUTPUT_DIR / f"avg_spectrum_{label.replace(' ', '_').lower()}.png",
                    dpi=150)
        plt.close()

        # Find peaks in average spectrum
        from scipy.signal import find_peaks
        band = avg_spec[f_mask]
        band_f = freq_axis[f_mask]
        peaks, props = find_peaks(band, height=np.max(band)*0.05, distance=10)

        print(f"\n{label} - Average spectrum peaks:")
        peak_freqs = band_f[peaks]
        peak_powers = 10*np.log10(band[peaks] + 1e-30)
        for f, p in sorted(zip(peak_freqs, peak_powers), key=lambda x: -x[1]):
            print(f"  {f:.1f} Hz ({p:.1f} dB)")

    # ================================================================
    # STEP 5: Zoomed spectrogram of first ~10 symbols
    # ================================================================
    print("\n=== STEP 5: ZOOMED SPECTROGRAM (first symbols) ===")

    for label, burst in [("Capture 1", burst1), ("Capture 2", burst2)]:
        # Use very high resolution spectrogram for first 0.5s
        first_part = burst[:int(0.5*FS)]

        nfft = 2048
        hop = 8  # very fine time steps
        f, t, Sxx = spectrogram(first_part, fs=FS, nperseg=nfft,
                                noverlap=nfft-hop, nfft=nfft)
        f_mask = (f >= 1400) & (f <= 1700)

        fig, ax = plt.subplots(figsize=(20, 8))
        Sxx_db = 10*np.log10(Sxx[f_mask] + 1e-30)
        ax.pcolormesh(t*1000, f[f_mask], Sxx_db, shading='gouraud', cmap='viridis')
        ax.set_xlabel("Time (ms)")
        ax.set_ylabel("Frequency (Hz)")
        ax.set_title(f"{label} - Zoomed spectrogram (first 500ms)\n"
                     f"Count horizontal lines = tones, measure time between transitions = symbol period")
        ax.grid(True, alpha=0.2, color='white')

        plt.tight_layout()
        plt.savefig(OUTPUT_DIR / f"zoom_{label.replace(' ', '_').lower()}.png",
                    dpi=200)
        plt.close()

    print("\nDone! Check output_raw/ for all plots.")


if __name__ == "__main__":
    main()
