#!/usr/bin/env python3
"""
FT2 Test NSPS=288 h=0.5 — the baud/h combination that gives:
  - baud=41.67 Bd, 90 symbols in 2.16s
  - tone_spacing=20.83 Hz (same as measured)
  - BW = 7×20.83 = 145.8 Hz ≈ 150 Hz (matches ft2.it)
  - 90 symbols = enough for LDPC(174,91) + sync

Also tests NSPS=320 h=0.5 (37.5 Bd, 81 symbols) and other candidates.

Usage:
    uv run python ft2_test_nsps288.py ft2_capture.wav
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
from scipy import signal as sig


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


def extract_symbols(data, sr, nsps, fmin, fmax, t_start, t_end, h):
    """Extract per-symbol frequencies and metadata."""
    i0 = int(t_start * sr)
    i1 = int(t_end * sr)
    seg_data = data[i0:i1]
    nfft = nsps * 8
    df = sr / nfft
    baud = sr / nsps
    spacing = h * baud
    n_sym = len(seg_data) // nsps

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
    """Quantize frequencies to tone indices."""
    good = snrs > snr_threshold
    if np.sum(good) < 5:
        return np.zeros(len(freqs), dtype=int), good, 1e20, 0

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
    return tone_idx, good, best_rms, best_base


def scan_alignment(data, sr, nsps, fmin, fmax, h, t_center, n_target_sym,
                   scan_range_ms=100):
    """Fine-scan start time to find best symbol alignment."""
    baud = sr / nsps
    spacing = h * baud
    sym_dur = nsps / sr

    offsets = np.arange(-scan_range_ms, scan_range_ms + 1, 0.5) / 1000.0
    scores = []

    for offset in offsets:
        t0 = t_center + offset
        t1 = t0 + n_target_sym * sym_dur
        if t0 < 0 or t1 > len(data) / sr:
            scores.append(0)
            continue

        freqs, powers, snrs, n_sym = extract_symbols(
            data, sr, nsps, fmin, fmax, t0, t1, h)
        # Score: average SNR of all symbols (higher = better alignment)
        scores.append(np.mean(snrs))

    scores = np.array(scores)
    best_idx = np.argmax(scores)
    return offsets[best_idx], offsets, scores


def costas_search_on_tones(tone_idx, good_mask, costas_arrays, positions, n_frame):
    """Search for best Costas array match."""
    n = len(tone_idx)
    best = (0, -1, 0, 0, [])

    for ci, costas in enumerate(costas_arrays):
        for off in range(max(1, n - n_frame)):
            for tb in range(int(tone_idx[off:off + n_frame].min()) - 8,
                            int(tone_idx[off:off + n_frame].max()) + 2):
                score = 0
                for sp in positions:
                    seg = tone_idx[off + sp:off + sp + len(costas)]
                    if len(seg) < len(costas):
                        continue
                    for k in range(len(costas)):
                        if abs(seg[k] - (tb + costas[k])) <= 1:
                            score += 1
                if score > best[0]:
                    best = (score, ci, off, tb, costas.tolist())

    return best


def main():
    wav_path = sys.argv[1] if len(sys.argv) > 1 else "ft2_capture.wav"
    Path("output_288").mkdir(exist_ok=True)

    data, sr = load_wav(wav_path)
    total_dur = len(data) / sr
    print(f"Loaded: {sr} Hz, {total_dur:.2f}s")

    fmin, fmax = 1440, 1720

    # Middle burst TX boundaries (from previous analysis)
    tx_on = 5.879
    tx_off = 8.041
    tx_dur = tx_off - tx_on

    print(f"\nMiddle burst: {tx_on:.3f}-{tx_off:.3f}s ({tx_dur:.3f}s)")

    # ================================================================
    # Test ALL NSPS values that give integer symbol counts in ~2.16s
    # and have enough capacity for LDPC(174,91)
    # ================================================================
    print(f"\n{'='*70}")
    print(f"  SYSTEMATIC NSPS SCAN (NSPS values giving enough symbols)")
    print(f"{'='*70}")

    # For LDPC(174,91) with 8-GFSK (3 bits/sym):
    # Min data symbols = ceil(174/3) = 58
    # With 21 sync symbols: min total = 79
    # Max from 2.16s: baud × 2.16 symbols

    candidates = []
    for nsps in range(120, 700):
        baud = sr / nsps
        n_sym = tx_dur * baud
        # We want integer-ish symbol count
        if abs(n_sym - round(n_sym)) < 0.15 and round(n_sym) >= 58:
            n_int = round(n_sym)
            # Test h values
            for h in [0.25, 0.5, 0.75, 1.0]:
                spacing = h * baud
                n_tones_in_bw = 150 / spacing  # expected number of tones in 150 Hz
                # 8-GFSK has 8 tones spanning 7×spacing
                bw_8gfsk = 7 * spacing
                if 100 < bw_8gfsk < 200:  # reasonable bandwidth
                    candidates.append((nsps, baud, h, n_int, spacing, bw_8gfsk))

    # Sort by how well n_sym fits an integer
    candidates.sort(key=lambda x: abs(tx_dur * x[1] - x[3]))

    print(f"\n  Found {len(candidates)} viable candidates")
    print(f"  {'NSPS':>5} {'Baud':>7} {'h':>5} {'N_sym':>6} {'Spacing':>8} "
          f"{'BW_8tone':>9} {'Data_sym':>9} {'Data_bits':>10}")
    print(f"  {'-'*70}")

    for nsps, baud, h, n_sym, spacing, bw in candidates[:30]:
        data_sym = n_sym - 21  # assuming 21 sync symbols
        data_bits = data_sym * 3
        enough = "✓" if data_bits >= 174 else "✗"
        print(f"  {nsps:5d} {baud:7.2f} {h:5.2f} {n_sym:6d} {spacing:8.2f} "
              f"{bw:9.2f} {data_sym:9d} {data_bits:10d} {enough}")

    # ================================================================
    # Detailed analysis of top candidates
    # ================================================================
    print(f"\n{'='*70}")
    print(f"  DETAILED ANALYSIS OF TOP CANDIDATES")
    print(f"{'='*70}")

    # Focus on the most promising: NSPS=288, 240, 320, 360, 480, 576
    top_candidates = [
        (288, 0.5, 90, "NSPS288_h0.5"),
        (240, 0.5, 108, "NSPS240_h0.5"),
        (320, 0.5, 81, "NSPS320_h0.5"),
        (360, 0.5, 72, "NSPS360_h0.5"),
        (192, 0.5, 135, "NSPS192_h0.5"),
        (288, 0.75, 90, "NSPS288_h0.75"),
        (576, 1.0, 45, "NSPS576_h1.0"),
    ]

    all_results = {}

    for nsps, h, expected_sym, label in top_candidates:
        baud = sr / nsps
        spacing = h * baud

        # Fine-tune alignment
        best_off, offsets, scores = scan_alignment(
            data, sr, nsps, fmin, fmax, h, tx_on,
            expected_sym, scan_range_ms=int(1000 * nsps / sr * 2))

        t_start = tx_on + best_off

        # Extract symbols
        freqs, powers, snrs, n_sym = extract_symbols(
            data, sr, nsps, fmin, fmax, t_start,
            t_start + expected_sym * nsps / sr + 0.01, h)

        tones, good, rms, base = quantize_tones(freqs, snrs, spacing)

        n_good = int(np.sum(good))
        good_tones = tones[good] if n_good > 0 else np.array([0])
        t_range = int(good_tones.max() - good_tones.min()) + 1

        # Compute average SNR of good symbols
        avg_snr = np.mean(snrs[good]) if n_good > 0 else 0

        print(f"\n  {label} (baud={baud:.2f}, spacing={spacing:.2f} Hz):")
        print(f"    Symbols: {n_sym} total, {n_good} good (SNR>2)")
        print(f"    RMS quant error: {rms:.2f} Hz")
        print(f"    Tone range: {t_range} tones ({good_tones.min()}-{good_tones.max()})")
        print(f"    Avg good-symbol SNR: {avg_snr:.1f}")
        print(f"    Data capacity: {(n_sym - 21) * 3} bits (need ≥174)")

        all_results[label] = {
            'nsps': nsps, 'h': h, 'baud': baud, 'spacing': spacing,
            'tones': tones, 'good': good, 'freqs': freqs, 'snrs': snrs,
            'rms': rms, 'n_sym': n_sym, 'n_good': n_good,
            'avg_snr': avg_snr, 't_start': t_start,
        }

        # Print tone sequence
        good_seq = tones[good]
        print(f"    Good-symbol tones: {good_seq[:60].tolist()}")

    # ================================================================
    # Compare: which NSPS gives the most consistent tone quantization?
    # ================================================================
    print(f"\n{'='*70}")
    print(f"  COMPARISON: QUANTIZATION QUALITY")
    print(f"{'='*70}")

    print(f"\n  {'Label':<20} {'RMS':>6} {'AvgSNR':>8} {'Tones':>6} "
          f"{'Good%':>6} {'DataBits':>9}")
    print(f"  {'-'*60}")

    for label in sorted(all_results.keys()):
        r = all_results[label]
        good_pct = 100 * r['n_good'] / r['n_sym'] if r['n_sym'] > 0 else 0
        data_bits = (r['n_sym'] - 21) * 3
        gt = r['tones'][r['good']]
        t_range = int(gt.max() - gt.min()) + 1 if len(gt) > 0 else 0
        print(f"  {label:<20} {r['rms']:6.2f} {r['avg_snr']:8.1f} "
              f"{t_range:6d} {good_pct:5.1f}% {data_bits:9d}")

    # ================================================================
    # Symbol-aligned waterfall for NSPS=288 h=0.5
    # ================================================================
    print(f"\n{'='*70}")
    print(f"  WATERFALLS FOR TOP 3 CANDIDATES")
    print(f"{'='*70}")

    for label in ["NSPS288_h0.5", "NSPS240_h0.5", "NSPS576_h1.0"]:
        r = all_results[label]
        nsps, h, baud = r['nsps'], r['h'], r['baud']
        t_start = r['t_start']
        sym_dur = nsps / sr
        t_end = t_start + r['n_sym'] * sym_dur

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

        fig, ax = plt.subplots(figsize=(16, max(8, n_sym * 0.12)))
        ax.imshow(wf_db, aspect='auto', origin='lower',
                  extent=[wf_freqs[mask][0], wf_freqs[mask][-1], 0, n_sym],
                  cmap='inferno', vmin=vmin, vmax=vmax)

        for i in range(n_sym):
            ax.axhline(i, color='white', alpha=0.06, linewidth=0.3)

        # Mark tone grid
        spacing = h * baud
        if r['n_good'] > 0:
            base_freq = r['freqs'][r['good']].min()
            for t in range(-1, 9):
                f = base_freq + t * spacing
                if fmin <= f <= fmax:
                    ax.axvline(f, color='cyan', alpha=0.15, linewidth=0.5)

        ax.set_xlabel("Frequency (Hz)")
        ax.set_ylabel("Symbol index")
        ax.set_title(f"{label}\nNSPS={nsps}, baud={baud:.2f}, h={h}, "
                     f"spacing={spacing:.2f} Hz, {n_sym} sym")
        plt.colorbar(ax.images[0], ax=ax, label='dB')
        plt.tight_layout()
        plt.savefig(f"output_288/waterfall_{label}.png", dpi=200)
        plt.close()
        print(f"  Saved: output_288/waterfall_{label}.png")

    # ================================================================
    # Costas search on NSPS=288 h=0.5
    # ================================================================
    print(f"\n{'='*70}")
    print(f"  COSTAS SEARCH ON NSPS=288 h=0.5 (90 symbols)")
    print(f"{'='*70}")

    r = all_results["NSPS288_h0.5"]
    tones = r['tones']
    good = r['good']

    # Generate all 7x7 Costas arrays
    all_costas_7 = generate_all_costas(7)
    print(f"  Testing {len(all_costas_7)} Costas arrays of order 7")

    # Try at standard FT8-like positions scaled to 90 symbols
    # For 90 symbols with 3 Costas blocks:
    # 90 = 3×7 + 69 data
    # Possible positions: evenly spaced
    position_sets = [
        ([0, 36, 72], 90, "FT8-like [0,36,72]"),
        ([0, 42, 83], 90, "Evenly spaced"),
        ([0, 30, 60], 90, "[0,30,60]"),
        ([0, 36, 83], 90, "[0,36,83]"),
        ([0, 45, 83], 90, "[0,45,83]"),
    ]

    for positions, n_frame, pos_label in position_sets:
        max_possible = 7 * len(positions)
        best = costas_search_on_tones(tones, good, all_costas_7, positions, n_frame)
        score, ci, off, tb, pat = best
        frac = score / max_possible
        print(f"  {pos_label}: best {score}/{max_possible} ({frac:.0%}) "
              f"array #{ci}: {pat} off={off}")

    # Also try order 4 Costas (FT4-style)
    all_costas_4 = generate_all_costas(4)
    print(f"\n  Testing {len(all_costas_4)} Costas arrays of order 4")

    position_sets_4 = [
        ([0, 22, 45, 67], 90, "Evenly spaced 4x"),
        ([0, 18, 39, 57], 77, "FT4-like [0,18,39,57]"),
        ([0, 23, 46, 69], 90, "[0,23,46,69]"),
    ]

    for positions, n_frame, pos_label in position_sets_4:
        max_possible = 4 * len(positions)
        best = costas_search_on_tones(tones, good, all_costas_4, positions, n_frame)
        score, ci, off, tb, pat = best
        frac = score / max_possible
        print(f"  {pos_label}: best {score}/{max_possible} ({frac:.0%}) "
              f"array #{ci}: {pat} off={off}")

    # ================================================================
    # Key metric: peak sharpness comparison
    # ================================================================
    print(f"\n{'='*70}")
    print(f"  PEAK SHARPNESS COMPARISON")
    print(f"  (Higher = better symbol alignment)")
    print(f"{'='*70}")

    for label in sorted(all_results.keys()):
        r = all_results[label]
        nsps = r['nsps']
        t_start = r['t_start']
        seg = data[int(t_start * sr):int((t_start + r['n_sym'] * nsps / sr) * sr)]

        # Compute peak-to-average ratio per symbol
        nfft = nsps * 4
        df = sr / nfft
        ratios = []

        for s in range(r['n_sym']):
            chunk = seg[s * nsps:(s + 1) * nsps]
            if len(chunk) < nsps:
                break
            w = np.hanning(nsps)
            spec = np.abs(fft(chunk * w, n=nfft))[:nfft // 2]
            f = np.arange(nfft // 2) * df
            band = (f >= fmin) & (f <= fmax)
            if not np.any(band):
                continue
            sb = spec[band]
            if sb.mean() > 0:
                ratios.append(sb.max() / sb.mean())

        ratios = np.array(ratios)
        print(f"  {label:<20}: median peak/avg = {np.median(ratios):.2f}, "
              f"mean = {np.mean(ratios):.2f}")

    print(f"\n{'='*70}")
    print(f"  DONE — check output_288/")
    print(f"{'='*70}")


if __name__ == "__main__":
    main()
