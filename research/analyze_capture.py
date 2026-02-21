#!/usr/bin/env python3
"""
Targeted FT2 analysis for the clean BlackHole capture.
Signal at ~1450-1700 Hz, three bursts visible.
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

COSTAS = np.array([3, 1, 4, 0, 6, 5, 2])
SYNC_POS = [0, 36, 72]


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


def extract_peaks(data, sr, nsps, fmin, fmax, t_start, t_end):
    """Extract peak frequency per symbol with parabolic interpolation."""
    i0 = int(t_start * sr)
    i1 = int(t_end * sr)
    seg = data[i0:i1]
    nfft = nsps * 8
    df = sr / nfft
    n_sym = len(seg) // nsps

    freqs_out = np.zeros(n_sym)
    powers_out = np.zeros(n_sym)

    for s in range(n_sym):
        chunk = seg[s * nsps:(s + 1) * nsps]
        w = np.hanning(len(chunk))
        spec = np.abs(fft(chunk * w, n=nfft))[:nfft // 2]
        f = np.arange(nfft // 2) * df
        band = (f >= fmin) & (f <= fmax)
        if not np.any(band):
            continue
        sb = spec[band]
        fb = f[band]
        pk = np.argmax(sb)
        pk_abs = np.where(band)[0][pk]
        freqs_out[s] = parabolic_peak(spec, pk_abs, df)
        powers_out[s] = sb[pk]

    return freqs_out, powers_out


def waterfall(data, sr, nsps, fmin, fmax, t_start, t_end, title, outpath):
    """Symbol-aligned waterfall."""
    i0 = int(t_start * sr)
    i1 = int(t_end * sr)
    seg = data[i0:i1]
    baud = sr / nsps
    nfft = nsps * 4
    df = sr / nfft
    n_sym = len(seg) // nsps

    wf = np.zeros((n_sym, nfft // 2))
    freqs = np.arange(nfft // 2) * df

    for s in range(n_sym):
        chunk = seg[s * nsps:(s + 1) * nsps]
        w = np.hanning(len(chunk))
        wf[s] = np.abs(fft(chunk * w, n=nfft))[:nfft // 2]

    mask = (freqs >= fmin) & (freqs <= fmax)
    wf_db = 10 * np.log10(wf[:, mask] + 1e-20)
    vmin, vmax = np.percentile(wf_db, [20, 99])

    fig, ax = plt.subplots(figsize=(16, max(6, n_sym * 0.12)))
    ax.imshow(wf_db, aspect="auto", origin="lower",
              extent=[freqs[mask][0], freqs[mask][-1], 0, n_sym],
              cmap="inferno", vmin=vmin, vmax=vmax)
    ax.set_xlabel("Frequency (Hz)")
    ax.set_ylabel("Symbol index")
    ax.set_title(f"{title}\nNSPS={nsps}, baud={baud:.2f}, df={df:.2f} Hz")
    for i in range(n_sym):
        ax.axhline(i, color="white", alpha=0.08, linewidth=0.3)
    plt.colorbar(ax.images[0], ax=ax, label="dB")
    plt.tight_layout()
    plt.savefig(outpath, dpi=200, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {outpath}")


def costas_search(peak_freqs, peak_powers, baud, tone_spacing, label):
    """Search for Costas 7x7 in extracted tone sequence."""
    print(f"\n--- Costas 7x7 Search ({label}) ---")
    good = peak_powers > np.percentile(peak_powers, 30)
    # Quantize to tones
    base = np.min(peak_freqs[good])
    tone_idx = np.round((peak_freqs - base) / tone_spacing).astype(int)
    tone_idx = np.clip(tone_idx, -2, 12)

    n = len(tone_idx)
    best_score = 0
    best_off = 0
    scores = []
    for off in range(max(1, n - 79)):
        score = 0
        for sp in SYNC_POS:
            seg = tone_idx[off + sp:off + sp + 7]
            if len(seg) < 7:
                continue
            for tb in range(max(0, int(seg.min()) - 7), int(seg.max()) + 1):
                m = int(np.sum((seg - tb) == COSTAS))
                score = max(score, m)
        scores.append(score)
        if score > best_score:
            best_score = score
            best_off = off

    scores = np.array(scores)
    print(f"  Best: {best_score}/21 at offset {best_off} (t={best_off/baud:.3f}s)")
    print(f"  Scores >=14: {np.sum(scores >= 14)}, >=17: {np.sum(scores >= 17)}")

    if best_score >= 10:
        for i, sp in enumerate(SYNC_POS):
            seg = tone_idx[best_off + sp:best_off + sp + 7]
            if len(seg) >= 7:
                print(f"  Costas block {i} (pos {sp}): {seg.tolist()}")

    return scores, best_off, tone_idx


def main():
    wav_path = sys.argv[1] if len(sys.argv) > 1 else "ft2_capture.wav"
    Path("output_capture").mkdir(exist_ok=True)

    data, sr = load_wav(wav_path)
    print(f"Loaded: {sr} Hz, {len(data)/sr:.2f}s")

    # Signal is at ~1450-1700 Hz based on spectrogram
    fmin, fmax = 1440, 1720

    # Three visible bursts
    bursts = [
        (0.2, 5.5, "burst0"),
        (5.8, 8.5, "burst1"),
        (8.5, 13.0, "burst2"),
    ]

    # Test candidate NSPS values
    candidate_nsps = [288, 320, 360, 384, 420, 448, 480, 512, 540, 576, 640, 720]

    # === STEP 1: Waterfalls for top candidates on cleanest burst ===
    print(f"\n{'='*60}")
    print(f"  STEP 1: SYMBOL-ALIGNED WATERFALLS (burst0, 1440-1720 Hz)")
    print(f"{'='*60}")
    t0, t1, label = bursts[0]
    for nsps in [360, 480, 576]:
        waterfall(data, sr, nsps, fmin, fmax, t0, t1,
                  f"Clean Capture — {label}",
                  f"output_capture/waterfall_{label}_nsps{nsps}.png")

    # === STEP 2: Peak frequency traces ===
    print(f"\n{'='*60}")
    print(f"  STEP 2: PEAK FREQUENCY TRACES")
    print(f"{'='*60}")

    fig, axes = plt.subplots(len(candidate_nsps), 1,
                              figsize=(18, 2.5 * len(candidate_nsps)))
    for ax, nsps in zip(axes, candidate_nsps):
        baud = sr / nsps
        pf, pp = extract_peaks(data, sr, nsps, fmin, fmax, t0, t1)
        good = pp > np.percentile(pp, 30)
        ax.plot(pf, ".-", markersize=2, linewidth=0.5, alpha=0.3)
        ax.plot(np.where(good, pf, np.nan), ".-", markersize=3, linewidth=0.8)
        ax.set_ylabel("Hz")
        ax.set_title(f"NSPS={nsps} ({baud:.2f} Bd) — {int(np.sum(good))} good symbols",
                     fontsize=9)
        ax.set_ylim(fmin, fmax)
        ax.grid(True, alpha=0.3)
    axes[-1].set_xlabel("Symbol index")
    plt.tight_layout()
    plt.savefig("output_capture/peak_traces_all.png", dpi=150)
    plt.close()
    print(f"  Saved: output_capture/peak_traces_all.png")

    # === STEP 3: Tone spacing measurement ===
    print(f"\n{'='*60}")
    print(f"  STEP 3: TONE SPACING & QUANTIZATION")
    print(f"{'='*60}")

    for nsps in [360, 480, 576]:
        baud = sr / nsps
        print(f"\n--- NSPS={nsps} ({baud:.2f} Bd) ---")

        for bi, (bt0, bt1, blabel) in enumerate(bursts):
            pf, pp = extract_peaks(data, sr, nsps, fmin, fmax, bt0, bt1)
            good = pp > np.percentile(pp, 50)
            gf = pf[good]

            if len(gf) < 5:
                continue

            # Step histogram
            diffs = np.abs(np.diff(gf))
            diffs = diffs[diffs > 2]
            if len(diffs) == 0:
                continue

            # Quantize with different h values
            for h in [0.5, 0.75, 1.0]:
                spacing = h * baud
                best_rms = 1e20
                for base in np.arange(gf.min(), gf.min() + spacing, 0.5):
                    q = np.round((gf - base) / spacing) * spacing + base
                    rms = np.sqrt(np.mean((gf - q) ** 2))
                    if rms < best_rms:
                        best_rms = rms
                        best_base = base

                q_idx = np.round((gf - best_base) / spacing).astype(int)
                n_tones = int(q_idx.max() - q_idx.min()) + 1
                print(f"  {blabel} h={h:.2f}: RMS={best_rms:.2f} Hz, "
                      f"tones={n_tones}, range={int(q_idx.min())}-{int(q_idx.max())}")

    # === STEP 4: Costas sync search ===
    print(f"\n{'='*60}")
    print(f"  STEP 4: COSTAS 7x7 SYNC DETECTION")
    print(f"{'='*60}")

    for nsps in [360, 480, 576]:
        baud = sr / nsps
        for h in [0.5, 0.75, 1.0]:
            spacing = h * baud
            for bt0, bt1, blabel in bursts:
                pf, pp = extract_peaks(data, sr, nsps, fmin, fmax, bt0, bt1)
                scores, best_off, tone_idx = costas_search(
                    pf, pp, baud, spacing,
                    f"NSPS={nsps}, h={h}, {blabel}")

    print(f"\n  Done. Check output_capture/ for plots.")


if __name__ == "__main__":
    main()
