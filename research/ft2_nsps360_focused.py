#!/usr/bin/env python3
"""
FT2 Focused Analysis — NSPS=360, baud≈33.33.

Key constraints:
- TX duration: 2.16s = 72 symbols at 33.33 Bd
- 8-GFSK modulation (8 tones)
- Need enough bits for LDPC(174,91) → 58 data symbols minimum
- With 2 Costas blocks (14 sync): 72 - 14 = 58 data → 174 bits = EXACTLY enough
- With 1 Costas block (7 sync): 72 - 7 = 65 data → 195 bits (room for CRC, etc.)

Tests:
1. h=0.5 and h=0.75 quantization quality
2. All Costas arrays at various 2-block and 1-block positions
3. Cross-burst consistency
4. Symbol-aligned waterfalls

Usage:
    uv run python ft2_nsps360_focused.py ft2_capture.wav
"""

import sys
from pathlib import Path
from itertools import permutations

import numpy as np
import soundfile as sf
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
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


def is_costas(perm):
    n = len(perm)
    for d in range(1, n):
        diffs = set()
        for i in range(n - d):
            diff = perm[i + d] - perm[i]
            if diff in diffs:
                return False
            diffs.add(diff)
    return True


def generate_all_costas(n):
    costas_arrays = []
    for perm in permutations(range(n)):
        if is_costas(perm):
            costas_arrays.append(np.array(perm))
    return costas_arrays


def extract_symbols(data, sr, nsps, fmin, fmax, t_start, n_sym_max, h):
    """Extract per-symbol frequencies."""
    i0 = int(t_start * sr)
    seg_data = data[i0:]
    nfft = nsps * 8
    df = sr / nfft
    baud = sr / nsps
    spacing = h * baud
    n_sym = min(len(seg_data) // nsps, n_sym_max)

    freqs = np.zeros(n_sym)
    powers = np.zeros(n_sym)
    snrs = np.zeros(n_sym)

    for s in range(n_sym):
        chunk = seg_data[s * nsps:(s + 1) * nsps]
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


def quantize_tones(freqs, snrs, spacing, snr_threshold=2.0):
    good = snrs > snr_threshold
    if np.sum(good) < 5:
        return np.zeros(len(freqs), dtype=int), good, 1e20, 0

    gf = freqs[good]
    best_rms = 1e20
    best_base = gf.min()
    for base in np.arange(gf.min() - spacing, gf.min() + spacing, 0.2):
        q = np.round((gf - base) / spacing) * spacing + base
        rms = np.sqrt(np.mean((gf - q) ** 2))
        if rms < best_rms:
            best_rms = rms
            best_base = base

    tone_idx = np.round((freqs - best_base) / spacing).astype(int)
    return tone_idx, good, best_rms, best_base


def scan_alignment(data, sr, nsps, fmin, fmax, h, t_center, n_target_sym):
    """Fine-scan start time for optimal alignment."""
    sym_dur = nsps / sr
    offsets = np.arange(-50, 51, 0.5) / 1000.0  # ±50ms in 0.5ms steps
    scores = []

    for offset in offsets:
        t0 = t_center + offset
        if t0 < 0:
            scores.append(0)
            continue
        freqs, powers, snrs, n_sym = extract_symbols(
            data, sr, nsps, fmin, fmax, t0, n_target_sym, h)
        scores.append(np.mean(snrs))

    scores = np.array(scores)
    best_idx = np.argmax(scores)
    return offsets[best_idx], offsets, scores


def costas_search(tones, costas_arrays, positions, n_frame, tolerance=1):
    """Search for best Costas match."""
    n = len(tones)
    best = (0, -1, 0, 0, [])

    for ci, costas in enumerate(costas_arrays):
        clen = len(costas)
        for off in range(max(1, n - n_frame + 1)):
            for tb in range(-2, 10):
                score = 0
                for sp in positions:
                    seg_slice = tones[off + sp:off + sp + clen]
                    if len(seg_slice) < clen:
                        continue
                    for k in range(clen):
                        if abs(seg_slice[k] - (tb + costas[k])) <= tolerance:
                            score += 1
                if score > best[0]:
                    best = (score, ci, off, tb, costas.tolist())

    return best


def main():
    wav_path = sys.argv[1] if len(sys.argv) > 1 else "ft2_capture.wav"
    Path("output_360").mkdir(exist_ok=True)

    data, sr = load_wav(wav_path)
    total_dur = len(data) / sr
    print(f"Loaded: {sr} Hz, {total_dur:.2f}s")

    fmin, fmax = 1440, 1720
    nsps = 360
    baud = sr / nsps
    print(f"NSPS={nsps}, baud={baud:.2f} Bd")
    print(f"72 symbols in 2.16s")

    # Middle burst TX boundary
    tx_on = 5.879

    # Generate all Costas arrays
    all_costas_7 = generate_all_costas(7)
    all_costas_4 = generate_all_costas(4)
    print(f"Costas arrays: {len(all_costas_7)} order-7, {len(all_costas_4)} order-4")

    # ================================================================
    # Test h values
    # ================================================================
    print(f"\n{'='*60}")
    print(f"  STEP 1: h VALUE DETERMINATION")
    print(f"{'='*60}")

    h_results = {}

    for h in [0.25, 0.5, 0.75, 1.0]:
        spacing = h * baud
        bw = 7 * spacing

        # Align
        best_off, _, _ = scan_alignment(data, sr, nsps, fmin, fmax, h, tx_on, 72)
        t_start = tx_on + best_off

        freqs, powers, snrs, n_sym = extract_symbols(
            data, sr, nsps, fmin, fmax, t_start, 72, h)
        tones, good, rms, base = quantize_tones(freqs, snrs, spacing)

        n_good = int(np.sum(good))
        gt = tones[good]
        t_range = int(gt.max() - gt.min()) + 1 if n_good > 0 else 0

        h_results[h] = {
            'tones': tones, 'good': good, 'freqs': freqs, 'snrs': snrs,
            'rms': rms, 'base': base, 't_start': t_start, 'n_good': n_good,
            'spacing': spacing,
        }

        print(f"\n  h={h:.2f} (spacing={spacing:.2f} Hz, BW≈{bw:.1f} Hz):")
        print(f"    Good: {n_good}/{n_sym}, RMS={rms:.2f} Hz")
        print(f"    Tones: {t_range} ({gt.min()}-{gt.max()})")
        print(f"    Tone sequence: {gt.tolist()[:40]}")

    # ================================================================
    # Costas search for h=0.5 (best RMS match)
    # ================================================================
    print(f"\n{'='*60}")
    print(f"  STEP 2: COSTAS SEARCH (h=0.5, 72 symbols)")
    print(f"{'='*60}")

    h = 0.5
    r = h_results[h]
    tones = r['tones']

    # With 2 Costas blocks (14 sync + 58 data):
    # Possible positions for 2 blocks of 7 in 72-symbol frame
    two_block_positions = [
        ([0, 65], 72, "start+end"),
        ([0, 36], 72, "start+mid"),
        ([0, 44], 72, "start+44"),
        ([0, 58], 72, "start+58"),
        ([33, 65], 72, "33+end"),
        ([0, 32], 72, "start+32"),
    ]

    print("\n  Two Costas blocks (order 7):")
    for positions, n_frame, label in two_block_positions:
        max_poss = 7 * len(positions)
        best = costas_search(tones, all_costas_7, positions, n_frame)
        score, ci, off, tb, pat = best
        frac = score / max_poss
        marker = " ***" if frac > 0.7 else " **" if frac > 0.6 else ""
        print(f"    {label:15s}: {score:2d}/{max_poss} ({frac:.0%}) "
              f"array #{ci:3d}: {pat} off={off} base={tb}{marker}")

    # With 1 Costas block (7 sync + 65 data):
    one_block_positions = [
        ([0], 72, "start"),
        ([65], 72, "end"),
        ([33], 72, "middle"),
    ]

    print("\n  One Costas block (order 7):")
    for positions, n_frame, label in one_block_positions:
        max_poss = 7
        best = costas_search(tones, all_costas_7, positions, n_frame)
        score, ci, off, tb, pat = best
        frac = score / max_poss
        marker = " ***" if frac > 0.85 else " **" if frac > 0.7 else ""
        print(f"    {label:15s}: {score:2d}/{max_poss} ({frac:.0%}) "
              f"array #{ci:3d}: {pat} off={off} base={tb}{marker}")

    # With 3 Costas blocks (21 sync + 51 data):
    three_block_positions = [
        ([0, 33, 65], 72, "evenly_3x"),
        ([0, 26, 52], 72, "[0,26,52]"),
        ([0, 24, 48], 72, "[0,24,48]"),
    ]

    print("\n  Three Costas blocks (order 7):")
    for positions, n_frame, label in three_block_positions:
        max_poss = 7 * 3
        best = costas_search(tones, all_costas_7, positions, n_frame)
        score, ci, off, tb, pat = best
        frac = score / max_poss
        marker = " ***" if frac > 0.7 else " **" if frac > 0.6 else ""
        print(f"    {label:15s}: {score:2d}/{max_poss} ({frac:.0%}) "
              f"array #{ci:3d}: {pat} off={off} base={tb}{marker}")

    # Order-4 Costas with 4 blocks (FT4-like):
    four_block_positions_4 = [
        ([0, 18, 36, 54], 72, "evenly_4x"),
        ([0, 15, 37, 54], 72, "[0,15,37,54]"),
        ([0, 18, 40, 58], 72, "[0,18,40,58]"),
        ([0, 18, 54, 68], 72, "[0,18,54,68]"),
    ]

    print("\n  Four Costas blocks (order 4):")
    for positions, n_frame, label in four_block_positions_4:
        max_poss = 4 * 4
        best = costas_search(tones, all_costas_4, positions, n_frame)
        score, ci, off, tb, pat = best
        frac = score / max_poss
        marker = " ***" if frac > 0.7 else " **" if frac > 0.6 else ""
        print(f"    {label:15s}: {score:2d}/{max_poss} ({frac:.0%}) "
              f"array #{ci:3d}: {pat} off={off} base={tb}{marker}")

    # ================================================================
    # Same search for h=0.75
    # ================================================================
    print(f"\n{'='*60}")
    print(f"  STEP 3: COSTAS SEARCH (h=0.75, 72 symbols)")
    print(f"{'='*60}")

    h = 0.75
    r = h_results[h]
    tones_075 = r['tones']

    print("\n  Two Costas blocks (order 7):")
    for positions, n_frame, label in two_block_positions:
        max_poss = 14
        best = costas_search(tones_075, all_costas_7, positions, n_frame)
        score, ci, off, tb, pat = best
        frac = score / max_poss
        marker = " ***" if frac > 0.7 else " **" if frac > 0.6 else ""
        print(f"    {label:15s}: {score:2d}/{max_poss} ({frac:.0%}) "
              f"array #{ci:3d}: {pat} off={off} base={tb}{marker}")

    print("\n  Three Costas blocks (order 7):")
    for positions, n_frame, label in three_block_positions:
        max_poss = 21
        best = costas_search(tones_075, all_costas_7, positions, n_frame)
        score, ci, off, tb, pat = best
        frac = score / max_poss
        marker = " ***" if frac > 0.7 else " **" if frac > 0.6 else ""
        print(f"    {label:15s}: {score:2d}/{max_poss} ({frac:.0%}) "
              f"array #{ci:3d}: {pat} off={off} base={tb}{marker}")

    # ================================================================
    # Cross-burst validation
    # ================================================================
    print(f"\n{'='*60}")
    print(f"  STEP 4: CROSS-BURST VALIDATION")
    print(f"{'='*60}")

    # First burst (tail end visible)
    # Third burst (start visible)
    burst_windows = [
        (0.0, "burst0_tail"),
        (5.879, "burst1_full"),
        (12.68, "burst2_start"),
    ]

    for h in [0.5, 0.75]:
        spacing = h * baud
        print(f"\n  h={h}:")

        for t_start, label in burst_windows:
            freqs, powers, snrs, n_sym = extract_symbols(
                data, sr, nsps, fmin, fmax, t_start, 80, h)
            tones_b, good_b, rms_b, base_b = quantize_tones(freqs, snrs, spacing)
            n_good_b = int(np.sum(good_b))
            if n_good_b < 5:
                print(f"    {label}: too few good symbols ({n_good_b})")
                continue

            gt = tones_b[good_b]
            print(f"    {label}: {n_good_b} good, {gt.max()-gt.min()+1} tones, "
                  f"RMS={rms_b:.2f}")
            print(f"      Tones: {gt[:30].tolist()}")

    # ================================================================
    # Waterfall comparison
    # ================================================================
    print(f"\n{'='*60}")
    print(f"  STEP 5: WATERFALLS")
    print(f"{'='*60}")

    for h in [0.5, 0.75]:
        r = h_results[h]
        t_start = r['t_start']
        spacing = h * baud
        sym_dur = nsps / sr

        i0 = int(t_start * sr)
        seg = data[i0:i0 + 72 * nsps]
        nfft = nsps * 4
        df = sr / nfft

        waterfall = np.zeros((72, nfft // 2))
        wf_freqs = np.arange(nfft // 2) * df

        for s in range(72):
            chunk = seg[s * nsps:(s + 1) * nsps]
            if len(chunk) < nsps:
                break
            w = np.hanning(nsps)
            waterfall[s] = np.abs(fft(chunk * w, n=nfft))[:nfft // 2]

        mask = (wf_freqs >= fmin) & (wf_freqs <= fmax)
        wf_db = 10 * np.log10(waterfall[:, mask] + 1e-20)
        vmin, vmax = np.percentile(wf_db, [10, 99])

        fig, ax = plt.subplots(figsize=(16, 12))
        ax.imshow(wf_db, aspect='auto', origin='lower',
                  extent=[wf_freqs[mask][0], wf_freqs[mask][-1], 0, 72],
                  cmap='inferno', vmin=vmin, vmax=vmax)

        for i in range(72):
            ax.axhline(i, color='white', alpha=0.06, linewidth=0.3)

        # Mark tone grid
        base_freq = r['base']
        for t in range(-1, 10):
            f = base_freq + t * spacing
            if fmin <= f <= fmax:
                ax.axvline(f, color='cyan', alpha=0.2, linewidth=0.8,
                          label=f'tone {t}' if t <= 2 else None)

        ax.set_xlabel("Frequency (Hz)")
        ax.set_ylabel("Symbol index")
        ax.set_title(f"NSPS=360, h={h} — Middle Burst (72 symbols)\n"
                     f"baud={baud:.2f}, spacing={spacing:.2f} Hz, "
                     f"RMS={r['rms']:.2f} Hz")
        ax.legend(loc='upper right', fontsize=8)
        plt.colorbar(ax.images[0], ax=ax, label='dB')
        plt.tight_layout()
        plt.savefig(f"output_360/waterfall_h{h}.png", dpi=200)
        plt.close()
        print(f"  Saved: output_360/waterfall_h{h}.png")

    # ================================================================
    # Detailed tone sequence for h=0.5
    # ================================================================
    print(f"\n{'='*60}")
    print(f"  STEP 6: COMPLETE TONE SEQUENCE (h=0.5)")
    print(f"{'='*60}")

    r = h_results[0.5]
    tones = r['tones']
    good = r['good']
    snrs_arr = r['snrs']
    freqs_arr = r['freqs']

    print(f"\n  72-symbol tone sequence (NSPS=360, h=0.5):")
    print(f"  Tone spacing: {0.5 * baud:.2f} Hz")
    print(f"  Base freq: {r['base']:.2f} Hz")
    print()

    for s in range(72):
        snr_str = f"{snrs_arr[s]:6.1f}" if snrs_arr[s] > 0 else "  ---"
        marker = "***" if snrs_arr[s] > 10 else "** " if snrs_arr[s] > 3 else "*  " if snrs_arr[s] > 2 else "   "
        print(f"  [{s:2d}] tone={tones[s]:2d}  f={freqs_arr[s]:7.1f}  snr={snr_str} {marker}")

    print(f"\n  Tone indices only (good symbols):")
    print(f"  {tones[good].tolist()}")

    # ================================================================
    # Same for h=0.75
    # ================================================================
    print(f"\n{'='*60}")
    print(f"  STEP 7: COMPLETE TONE SEQUENCE (h=0.75)")
    print(f"{'='*60}")

    r = h_results[0.75]
    tones = r['tones']
    good = r['good']
    snrs_arr = r['snrs']
    freqs_arr = r['freqs']

    print(f"\n  72-symbol tone sequence (NSPS=360, h=0.75):")
    print(f"  Tone spacing: {0.75 * baud:.2f} Hz")
    print(f"  Base freq: {r['base']:.2f} Hz")
    print()

    for s in range(72):
        snr_str = f"{snrs_arr[s]:6.1f}" if snrs_arr[s] > 0 else "  ---"
        marker = "***" if snrs_arr[s] > 10 else "** " if snrs_arr[s] > 3 else "*  " if snrs_arr[s] > 2 else "   "
        print(f"  [{s:2d}] tone={tones[s]:2d}  f={freqs_arr[s]:7.1f}  snr={snr_str} {marker}")

    print(f"\n  Tone indices only (good symbols):")
    print(f"  {tones[good].tolist()}")

    print(f"\n{'='*60}")
    print(f"  DONE — check output_360/")
    print(f"{'='*60}")


if __name__ == "__main__":
    main()
