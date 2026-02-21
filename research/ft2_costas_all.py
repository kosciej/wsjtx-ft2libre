#!/usr/bin/env python3
"""
FT2 Costas Array Search — enumerate ALL 7x7 Costas arrays and test each one.

There are exactly 200 Costas arrays of order 7. We test all of them
at all reasonable sync positions within the frame, using both quantized
tones and frequency-domain correlation.

Also tests with fuzzy matching (±1 tone tolerance) to handle quantization errors.

Usage:
    uv run python ft2_costas_all.py ft2_capture.wav
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


def is_costas(perm):
    """Check if a permutation is a Costas array.
    A Costas array has all distinct displacement vectors."""
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
    """Generate all Costas arrays of order n."""
    # For n=7, brute-force is feasible (7! = 5040 permutations)
    costas_arrays = []
    for perm in permutations(range(n)):
        if is_costas(perm):
            costas_arrays.append(np.array(perm))
    return costas_arrays


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


def extract_tones(data, sr, nsps, fmin, fmax, t_start, t_end, h, power_pct=30):
    """Extract quantized tone indices from a burst."""
    i0 = int(t_start * sr)
    i1 = int(t_end * sr)
    seg = data[i0:i1]
    nfft = nsps * 8
    df = sr / nfft
    baud = sr / nsps
    spacing = h * baud
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
        pk = np.argmax(sb)
        pk_abs = np.where(band)[0][pk]
        freqs_out[s] = parabolic_peak(spec, pk_abs, df)
        powers_out[s] = sb[pk]

    good = powers_out > np.percentile(powers_out, power_pct)
    if np.sum(good) < 5:
        return freqs_out, powers_out, np.zeros(n_sym, dtype=int), good

    gf = freqs_out[good]
    best_rms = 1e20
    best_base = gf.min()
    for base in np.arange(gf.min() - spacing, gf.min() + spacing, 0.5):
        q = np.round((gf - base) / spacing) * spacing + base
        rms = np.sqrt(np.mean((gf - q) ** 2))
        if rms < best_rms:
            best_rms = rms
            best_base = base

    tone_idx = np.round((freqs_out - best_base) / spacing).astype(int)
    return freqs_out, powers_out, tone_idx, good


def test_costas_array(tone_idx, costas, positions, n_frame, tolerance=0):
    """Test a specific Costas array at specific positions.
    Returns (best_score, best_offset, best_base_tone)."""
    n = len(tone_idx)
    costas_len = len(costas)
    best_score = 0
    best_off = 0
    best_base = 0

    for off in range(max(1, n - n_frame)):
        for tb in range(int(tone_idx[off:off + n_frame].min()) - 8,
                        int(tone_idx[off:off + n_frame].max()) + 2):
            score = 0
            for sp in positions:
                seg = tone_idx[off + sp:off + sp + costas_len]
                if len(seg) < costas_len:
                    continue
                for k in range(costas_len):
                    expected = tb + costas[k]
                    actual = seg[k]
                    if tolerance == 0:
                        if actual == expected:
                            score += 1
                    else:
                        if abs(actual - expected) <= tolerance:
                            score += 1

            if score > best_score:
                best_score = score
                best_off = off
                best_base = tb

    return best_score, best_off, best_base


def frequency_domain_sync(data, sr, nsps, fmin, fmax, t_start, t_end,
                          costas, h, sync_positions, n_frame):
    """
    Frequency-domain Costas detection: for each candidate time offset,
    compute the spectral power at expected Costas tone frequencies and sum.
    This avoids quantization errors.
    """
    i0 = int(t_start * sr)
    i1 = int(t_end * sr)
    seg = data[i0:i1]
    baud = sr / nsps
    spacing = h * baud
    nfft = nsps * 2  # standard oversampling
    df = sr / nfft
    n_sym = len(seg) // nsps

    if n_sym < n_frame:
        return 0, 0, 0

    # Compute all symbol spectra
    spectra = np.zeros((n_sym, nfft // 2))
    for s in range(n_sym):
        chunk = seg[s * nsps:(s + 1) * nsps]
        w = np.hanning(nsps)
        padded = np.zeros(nfft)
        padded[:nsps] = chunk * w
        spectra[s] = np.abs(fft(padded))[:nfft // 2] ** 2

    freqs = np.arange(nfft // 2) * df
    f_lo = int(fmin / df)
    f_hi = int(fmax / df) + 1
    f_lo = max(0, f_lo)
    f_hi = min(nfft // 2, f_hi)

    best_score = 0
    best_off = 0
    best_f0 = 0
    costas_len = len(costas)

    # Scan over time offsets and base frequency
    for off in range(n_sym - n_frame + 1):
        for f0_bin in range(f_lo, f_hi):
            score = 0
            for sp in sync_positions:
                for k in range(costas_len):
                    sym_idx = off + sp + k
                    if sym_idx >= n_sym:
                        continue
                    tone_bin = f0_bin + int(round(costas[k] * spacing / df))
                    if 0 <= tone_bin < nfft // 2:
                        score += spectra[sym_idx, tone_bin]

            if score > best_score:
                best_score = score
                best_off = off
                best_f0 = f0_bin * df

    # Normalize
    noise = np.median(spectra[:, f_lo:f_hi])
    n_sync_tones = costas_len * len(sync_positions)
    if noise > 0:
        snr = best_score / (n_sync_tones * noise)
    else:
        snr = 0

    return snr, best_off, best_f0


def main():
    wav_path = sys.argv[1] if len(sys.argv) > 1 else "ft2_capture.wav"
    Path("output_costas").mkdir(exist_ok=True)

    data, sr = load_wav(wav_path)
    print(f"Loaded: {sr} Hz, {len(data)/sr:.2f}s")

    # Generate all 7x7 Costas arrays
    print("\nGenerating all 7x7 Costas arrays...")
    all_costas_7 = generate_all_costas(7)
    print(f"Found {len(all_costas_7)} Costas arrays of order 7")

    # Also check the FT8 Costas for sanity
    ft8_costas = np.array([3, 1, 4, 0, 6, 5, 2])
    assert any(np.array_equal(c, ft8_costas) for c in all_costas_7), "FT8 Costas should be in the list!"

    fmin, fmax = 1440, 1720

    # Use burst1 which has cleanest signal (39 good symbols with NSPS=576 h=1.0)
    bursts = [
        (0.2, 5.5, "burst0"),
        (5.8, 8.5, "burst1"),
    ]

    # Primary parameter set
    nsps = 576
    h = 1.0
    baud = sr / nsps
    spacing = h * baud

    # === STEP 1: Quantized tone search across ALL Costas arrays ===
    print(f"\n{'='*60}")
    print(f"  STEP 1: ALL {len(all_costas_7)} COSTAS ARRAYS — QUANTIZED TONE SEARCH")
    print(f"  NSPS={nsps}, h={h}, baud={baud:.2f}")
    print(f"{'='*60}")

    # Standard FT8 positions in 79-symbol frame
    sync_positions_79 = [0, 36, 72]

    for bt0, bt1, blabel in bursts:
        print(f"\n  {blabel} ({bt0}-{bt1}s):")
        freqs, powers, tones, good = extract_tones(
            data, sr, nsps, fmin, fmax, bt0, bt1, h, power_pct=20)

        n_good = int(np.sum(good))
        print(f"    {len(tones)} symbols, {n_good} good")

        if n_good < 20:
            print("    Not enough good symbols")
            continue

        results = []
        for ci, costas in enumerate(all_costas_7):
            score, off, base = test_costas_array(
                tones, costas, sync_positions_79, 79, tolerance=0)
            results.append((score, ci, off, base, costas.tolist()))

        results.sort(reverse=True)
        print(f"\n    Top 10 (exact match, 79-sym frame, positions [0,36,72]):")
        for score, ci, off, base, pat in results[:10]:
            print(f"      {score:2d}/21  array #{ci:3d}: {pat}  off={off} base={base}")

        # Also try with ±1 tolerance
        results_fuzzy = []
        for ci, costas in enumerate(all_costas_7):
            score, off, base = test_costas_array(
                tones, costas, sync_positions_79, 79, tolerance=1)
            results_fuzzy.append((score, ci, off, base, costas.tolist()))

        results_fuzzy.sort(reverse=True)
        print(f"\n    Top 10 (±1 tone tolerance, 79-sym frame):")
        for score, ci, off, base, pat in results_fuzzy[:10]:
            print(f"      {score:2d}/21  array #{ci:3d}: {pat}  off={off} base={base}")

    # === STEP 2: Try different sync position patterns ===
    print(f"\n{'='*60}")
    print(f"  STEP 2: ALTERNATIVE SYNC POSITIONS")
    print(f"{'='*60}")

    # Various position patterns to try
    position_sets = {
        "FT8 [0,36,72] 79sym": ([0, 36, 72], 79),
        "FT8 [0,36,72] 77sym": ([0, 36, 72], 77),
        "evenly_3x [0,26,52] 79sym": ([0, 26, 52], 79),
        "evenly_3x [0,24,48] 79sym": ([0, 24, 48], 79),
        "start_mid_end [0,36,72] 85sym": ([0, 36, 72], 85),
        "start_end [0,72] 79sym": ([0, 72], 79),
        "4_costas [0,18,39,57] 77sym": ([0, 18, 39, 57], 77),  # FT4-like
        "2_costas [0,36] 79sym": ([0, 36], 79),
    }

    bt0, bt1, blabel = bursts[1]  # Use burst1
    freqs, powers, tones, good = extract_tones(
        data, sr, nsps, fmin, fmax, bt0, bt1, h, power_pct=20)

    for pos_label, (positions, n_frame) in position_sets.items():
        max_possible = 7 * len(positions)
        best_overall = (0, -1, 0, 0, [])

        for ci, costas in enumerate(all_costas_7):
            score, off, base = test_costas_array(
                tones, costas, positions, n_frame, tolerance=1)
            if score > best_overall[0]:
                best_overall = (score, ci, off, base, costas.tolist())

        score, ci, off, base, pat = best_overall
        print(f"  {pos_label}: best {score}/{max_possible}  "
              f"array #{ci}: {pat}  off={off} base={base}")

    # === STEP 3: Frequency-domain search with top candidates ===
    print(f"\n{'='*60}")
    print(f"  STEP 3: FREQUENCY-DOMAIN SYNC CORRELATION")
    print(f"  (Top 5 Costas arrays from quantized search)")
    print(f"{'='*60}")

    # Get top candidates from step 1
    bt0, bt1, blabel = bursts[0]  # Use longest burst for freq-domain
    freqs_b0, powers_b0, tones_b0, good_b0 = extract_tones(
        data, sr, nsps, fmin, fmax, bt0, bt1, h, power_pct=20)

    # Get top 5 unique Costas arrays from step 1 results
    seen = set()
    top_costas = []
    for score, ci, off, base, pat in results_fuzzy[:20]:
        key = tuple(pat)
        if key not in seen:
            seen.add(key)
            top_costas.append((ci, np.array(pat)))
            if len(top_costas) >= 10:
                break

    # Always include FT8 Costas
    ft8_key = tuple(ft8_costas)
    if ft8_key not in seen:
        top_costas.append((-1, ft8_costas))

    for ci, costas in top_costas:
        snr, off, f0 = frequency_domain_sync(
            data, sr, nsps, fmin, fmax, bt0, bt1,
            costas, h, sync_positions_79, 79)
        print(f"  Array #{ci:3d} {costas.tolist()}: SNR={snr:.2f}, off={off}, f0={f0:.1f}")

    # === STEP 4: Also try NSPS=576 h=0.5 and h=0.75 ===
    print(f"\n{'='*60}")
    print(f"  STEP 4: ALTERNATIVE h VALUES")
    print(f"{'='*60}")

    for h_test in [0.5, 0.75]:
        spacing_test = h_test * baud
        print(f"\n  h={h_test} (spacing={spacing_test:.2f} Hz)")

        for bt0, bt1, blabel in bursts[:1]:
            freqs_t, powers_t, tones_t, good_t = extract_tones(
                data, sr, nsps, fmin, fmax, bt0, bt1, h_test, power_pct=20)

            best = (0, -1, [])
            for ci, costas in enumerate(all_costas_7):
                score, off, base = test_costas_array(
                    tones_t, costas, sync_positions_79, 79, tolerance=1)
                if score > best[0]:
                    best = (score, ci, costas.tolist())

            print(f"    {blabel}: best {best[0]}/21 array #{best[1]}: {best[2]}")

    # === STEP 5: Try Costas arrays of different orders ===
    print(f"\n{'='*60}")
    print(f"  STEP 5: COSTAS ARRAYS OF OTHER ORDERS (4, 5, 6)")
    print(f"{'='*60}")

    for order in [4, 5, 6]:
        costas_arrays = generate_all_costas(order)
        print(f"\n  Order {order}: {len(costas_arrays)} arrays")

        bt0, bt1, blabel = bursts[1]
        freqs_t, powers_t, tones_t, good_t = extract_tones(
            data, sr, nsps, fmin, fmax, bt0, bt1, 1.0, power_pct=20)

        # For order 4, try FT4-like 4 positions
        if order == 4:
            positions_list = [
                ([0, 18, 39, 57], 77, "FT4-like"),
                ([0, 20, 40, 60], 79, "evenly_4x"),
                ([0, 18, 36, 72], 79, "mixed"),
            ]
        elif order == 5:
            positions_list = [
                ([0, 16, 32, 48, 64], 79, "evenly_5x"),
                ([0, 15, 36, 57, 72], 79, "mixed_5x"),
            ]
        else:  # order 6
            positions_list = [
                ([0, 12, 24, 36, 48, 60], 79, "evenly_6x"),
                ([0, 12, 36, 48, 60, 72], 79, "mixed_6x"),
            ]

        for positions, n_frame, pos_label in positions_list:
            max_possible = order * len(positions)
            best = (0, -1, [])
            for ci, costas in enumerate(costas_arrays):
                score, off, base = test_costas_array(
                    tones_t, costas, positions, n_frame, tolerance=1)
                if score > best[0]:
                    best = (score, ci, costas.tolist())

            frac = best[0] / max_possible if max_possible > 0 else 0
            print(f"    {pos_label} pos={positions}: "
                  f"best {best[0]}/{max_possible} ({frac:.0%}) "
                  f"array #{best[1]}: {best[2]}")

    print(f"\n{'='*60}")
    print(f"  DONE — check output_costas/")
    print(f"{'='*60}")


if __name__ == "__main__":
    main()
