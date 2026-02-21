#!/usr/bin/env python3
"""
FT2 Sync Pattern Discovery — find repeating patterns across bursts.

Since Costas 7x7 sync detection failed, FT2 must use a different sync structure.
This script:
1. Extracts tone sequences from all three bursts
2. Cross-correlates between bursts to find common patterns
3. Looks for any repeating subsequences that could be sync markers
4. Tests FT4-style sync (different Costas positions or different pattern)

Usage:
    uv run python ft2_find_sync.py ft2_capture.wav
"""

import sys
from pathlib import Path

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


def extract_tones(data, sr, nsps, fmin, fmax, t_start, t_end, h):
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

    # Quantize to tone indices
    good = powers_out > np.percentile(powers_out, 30)
    if np.sum(good) < 5:
        return freqs_out, powers_out, np.zeros(n_sym, dtype=int), good

    # Find best base frequency
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


def cross_correlate_tones(tones_a, tones_b, max_lag=None):
    """Cross-correlate two tone sequences (integer-valued).
    Returns correlation as fraction of matching tones at each lag."""
    na, nb = len(tones_a), len(tones_b)
    if max_lag is None:
        max_lag = max(na, nb)

    results = []
    for lag in range(-max_lag, max_lag + 1):
        # Overlap region
        if lag >= 0:
            a_start, b_start = lag, 0
        else:
            a_start, b_start = 0, -lag

        overlap = min(na - a_start, nb - b_start)
        if overlap < 7:
            results.append((lag, 0, 0))
            continue

        a_seg = tones_a[a_start:a_start + overlap]
        b_seg = tones_b[b_start:b_start + overlap]

        # Count exact matches (allowing for base offset)
        best_matches = 0
        for offset in range(-3, 4):
            matches = np.sum(a_seg == (b_seg + offset))
            best_matches = max(best_matches, matches)

        results.append((lag, best_matches, overlap))

    return results


def find_repeating_subsequences(tone_seq, good_mask, min_len=4, max_len=12):
    """Find subsequences that appear multiple times in the tone sequence."""
    good_idx = np.where(good_mask)[0]
    if len(good_idx) < min_len * 2:
        return []

    # Use good-power symbols only
    gt = tone_seq[good_idx]
    # Normalize to start from 0
    gt = gt - gt.min()

    found = {}
    for length in range(min_len, max_len + 1):
        if length >= len(gt):
            break
        for i in range(len(gt) - length):
            pattern = tuple(gt[i:i + length])
            if pattern not in found:
                found[pattern] = []
            found[pattern].append(good_idx[i])

    # Filter to patterns that appear 2+ times
    repeating = [(pat, positions) for pat, positions in found.items()
                 if len(positions) >= 2]
    # Sort by length * count
    repeating.sort(key=lambda x: len(x[0]) * len(x[1]), reverse=True)
    return repeating[:20]


def test_ft4_sync(tone_idx, good_mask, baud, label):
    """Test FT4 sync structure: 4 Costas arrays at positions [0, 18, 39, 57] in 77-symbol frame."""
    FT4_COSTAS = np.array([0, 1, 3, 2])  # FT4 uses 4x4 Costas
    FT4_SYNC_POS = [0, 18, 39, 57]
    FT4_NSYM = 77  # FT4: 4+18+4+18+4+18+4+7 = 77

    print(f"\n  --- FT4-style Sync Test ({label}) ---")
    good = np.where(good_mask)[0]
    if len(good) < 20:
        print("    Not enough good symbols")
        return

    gt = tone_idx.copy()
    n = len(gt)
    best_score = 0
    best_off = 0

    for off in range(max(1, n - FT4_NSYM)):
        score = 0
        for sp in FT4_SYNC_POS:
            seg = gt[off + sp:off + sp + 4]
            if len(seg) < 4:
                continue
            for tb in range(int(seg.min()) - 4, int(seg.max()) + 2):
                m = int(np.sum((seg - tb) == FT4_COSTAS))
                score = max(score, m)
        if score > best_score:
            best_score = score
            best_off = off

    print(f"    Best: {best_score}/16 at offset {best_off}")


def test_generic_costas(tone_idx, good_mask, n_syms_frame, costas_len,
                        costas_positions, costas_pattern, baud, label):
    """Test a generic Costas sync pattern."""
    n = len(tone_idx)
    best_score = 0
    best_off = 0

    for off in range(max(1, n - n_syms_frame)):
        score = 0
        for sp in costas_positions:
            seg = tone_idx[off + sp:off + sp + costas_len]
            if len(seg) < costas_len:
                continue
            for tb in range(int(seg.min()) - 8, int(seg.max()) + 2):
                m = int(np.sum((seg - tb) == costas_pattern))
                score = max(score, m)
        if score > best_score:
            best_score = score
            best_off = off

    max_possible = costas_len * len(costas_positions)
    print(f"    {label}: Best {best_score}/{max_possible} at offset {best_off}")
    return best_score, best_off


def exhaustive_sync_search(tone_idx, good_mask, baud, label):
    """
    Brute-force search: try ALL possible 7-tone patterns at ALL possible
    position triples within a 79-symbol frame.
    """
    print(f"\n  --- Exhaustive Sync Search ({label}) ---")

    n = len(tone_idx)
    if n < 79:
        print(f"    Only {n} symbols, need 79+")
        return

    # First, try the standard Costas at different frame lengths
    print("  Testing standard Costas [3,1,4,0,6,5,2] at various frame sizes:")
    COSTAS = np.array([3, 1, 4, 0, 6, 5, 2])

    # Try different total frame lengths (not just 79)
    for n_frame in [77, 79, 85, 87, 90, 103, 105]:
        if n_frame > n:
            continue
        # Try Costas at start, middle, end
        positions_to_try = [
            [0, n_frame // 2 - 3, n_frame - 7],  # evenly spaced
            [0, (n_frame - 7) // 2, n_frame - 7],  # start, mid, end
        ]
        # Also try positions based on FT8 ratios
        if n_frame >= 72 + 7:
            positions_to_try.append([0, 36, 72])

        for positions in positions_to_try:
            score, off = test_generic_costas(
                tone_idx, good_mask, n_frame, 7, positions, COSTAS,
                baud, f"frame={n_frame}, pos={positions}")

    # Try without the standard Costas — look for ANY 7-tone pattern
    # that repeats at regular intervals
    print("\n  Searching for ANY repeating 7-symbol pattern:")
    best_pattern_score = 0
    best_pattern = None
    best_pattern_info = None

    for off in range(min(20, n - 79)):
        for frame_len in [77, 79, 85, 87, 90, 103, 105]:
            if off + frame_len > n:
                continue
            # Extract the first 7 symbols as candidate pattern
            candidate = tone_idx[off:off + 7]
            base = candidate.min()
            candidate_norm = candidate - base

            # Check if this pattern repeats at regular positions
            for mid_pos in range(20, frame_len - 14):
                seg_mid = tone_idx[off + mid_pos:off + mid_pos + 7]
                for end_pos in range(mid_pos + 14, frame_len - 6):
                    seg_end = tone_idx[off + end_pos:off + end_pos + 7]

                    # Count matches with base offset search
                    total_matches = 7  # first one matches by definition
                    for seg in [seg_mid, seg_end]:
                        best_m = 0
                        for tb in range(int(seg.min()) - 8, int(seg.max()) + 2):
                            m = int(np.sum((seg - tb) == candidate_norm))
                            best_m = max(best_m, m)
                        total_matches += best_m

                    if total_matches > best_pattern_score:
                        best_pattern_score = total_matches
                        best_pattern = candidate_norm.copy()
                        best_pattern_info = (off, mid_pos, end_pos, frame_len)

    if best_pattern is not None:
        print(f"    Best: {best_pattern_score}/21, pattern={best_pattern.tolist()}")
        print(f"    At offset={best_pattern_info[0]}, positions=[0, {best_pattern_info[1]}, {best_pattern_info[2]}], frame_len={best_pattern_info[3]}")


def analyze_frame_structure(tone_idx, powers, good_mask, baud, label):
    """
    Analyze the structure by looking at power envelope and tone transitions.
    If there are sync tones, they should have distinctive power or pattern characteristics.
    """
    print(f"\n  --- Frame Structure Analysis ({label}) ---")

    n = len(tone_idx)
    good = good_mask

    # Power profile — look for regular structure
    # Smooth power to reduce noise
    from scipy.ndimage import uniform_filter1d
    smooth_power = uniform_filter1d(powers.astype(float), size=5)

    # Autocorrelation of power envelope
    p_centered = smooth_power - smooth_power.mean()
    autocorr = np.correlate(p_centered, p_centered, mode="full")
    autocorr = autocorr[len(autocorr) // 2:]  # positive lags only
    autocorr /= autocorr[0] + 1e-20

    # Find peaks in autocorrelation
    from scipy.signal import find_peaks
    peaks, props = find_peaks(autocorr[5:], height=0.1, distance=10)
    peaks += 5  # offset

    if len(peaks) > 0:
        print(f"    Power autocorrelation peaks (symbols): {peaks[:10].tolist()}")
        print(f"    Corresponding times: {[f'{p/baud:.3f}s' for p in peaks[:10]]}")
        # The first strong peak might be the frame length
        for p in peaks[:5]:
            print(f"      Period={p} symbols = {p/baud:.3f}s")
    else:
        print("    No significant autocorrelation peaks in power")

    # Tone transition autocorrelation
    # Create binary signal: 1 where tone changes, 0 where it stays
    transitions = np.abs(np.diff(tone_idx)).astype(float)
    transitions[~good[1:]] = 0  # mask bad symbols

    if len(transitions) > 20:
        t_centered = transitions - transitions.mean()
        t_autocorr = np.correlate(t_centered, t_centered, mode="full")
        t_autocorr = t_autocorr[len(t_autocorr) // 2:]
        t_autocorr /= t_autocorr[0] + 1e-20

        t_peaks, _ = find_peaks(t_autocorr[5:], height=0.08, distance=10)
        t_peaks += 5
        if len(t_peaks) > 0:
            print(f"    Transition autocorrelation peaks: {t_peaks[:10].tolist()}")

    return autocorr


def main():
    wav_path = sys.argv[1] if len(sys.argv) > 1 else "ft2_capture.wav"
    Path("output_sync").mkdir(exist_ok=True)

    data, sr = load_wav(wav_path)
    print(f"Loaded: {sr} Hz, {len(data)/sr:.2f}s")

    fmin, fmax = 1440, 1720

    bursts = [
        (0.2, 5.5, "burst0"),
        (5.8, 8.5, "burst1"),
        (8.5, 13.0, "burst2"),
    ]

    # Test the most promising parameter combinations
    param_sets = [
        (576, 1.0, "NSPS576_h1.0"),
        (576, 0.75, "NSPS576_h0.75"),
        (480, 0.75, "NSPS480_h0.75"),
        (360, 0.5, "NSPS360_h0.5"),
    ]

    all_tones = {}  # (params, burst) -> (tone_idx, good_mask)

    # === STEP 1: Extract tone sequences ===
    print(f"\n{'='*60}")
    print(f"  STEP 1: EXTRACT TONE SEQUENCES")
    print(f"{'='*60}")

    for nsps, h, plabel in param_sets:
        baud = sr / nsps
        spacing = h * baud
        print(f"\n  {plabel} (baud={baud:.2f}, spacing={spacing:.2f} Hz)")

        for bt0, bt1, blabel in bursts:
            freqs, powers, tones, good = extract_tones(
                data, sr, nsps, fmin, fmax, bt0, bt1, h)
            all_tones[(plabel, blabel)] = (tones, good, freqs, powers)

            gt = tones[good]
            n_good = int(np.sum(good))
            n_tones = int(gt.max() - gt.min()) + 1 if n_good > 0 else 0
            print(f"    {blabel}: {len(tones)} syms, {n_good} good, "
                  f"{n_tones} tones, range [{gt.min()}-{gt.max()}]" if n_good > 0 else
                  f"    {blabel}: {len(tones)} syms, {n_good} good")

    # === STEP 2: Cross-correlate between bursts ===
    print(f"\n{'='*60}")
    print(f"  STEP 2: CROSS-CORRELATION BETWEEN BURSTS")
    print(f"{'='*60}")

    for nsps, h, plabel in param_sets:
        baud = sr / nsps
        print(f"\n  {plabel}")

        burst_keys = [(plabel, b[2]) for b in bursts]
        for i in range(len(burst_keys)):
            for j in range(i + 1, len(burst_keys)):
                ki, kj = burst_keys[i], burst_keys[j]
                ti, gi = all_tones[ki][:2]
                tj, gj = all_tones[kj][:2]

                # Use good symbols only
                ti_good = ti[gi]
                tj_good = tj[gj]

                if len(ti_good) < 10 or len(tj_good) < 10:
                    continue

                results = cross_correlate_tones(ti_good, tj_good,
                                                max_lag=max(len(ti_good), len(tj_good)))
                # Find best
                best = max(results, key=lambda x: x[1] / max(x[2], 1))
                lag, matches, overlap = best
                if overlap > 0:
                    frac = matches / overlap
                    print(f"    {ki[1]} vs {kj[1]}: best lag={lag}, "
                          f"matches={matches}/{overlap} ({frac:.1%})")

    # === STEP 3: Look for repeating subsequences within each burst ===
    print(f"\n{'='*60}")
    print(f"  STEP 3: REPEATING SUBSEQUENCES WITHIN BURSTS")
    print(f"{'='*60}")

    for nsps, h, plabel in param_sets[:2]:  # top 2 candidates
        baud = sr / nsps
        for bt0, bt1, blabel in bursts[:1]:  # longest burst only
            key = (plabel, blabel)
            tones, good = all_tones[key][:2]
            print(f"\n  {plabel}, {blabel}:")
            repeating = find_repeating_subsequences(tones, good, min_len=4, max_len=8)
            for pat, positions in repeating[:10]:
                print(f"    Pattern {list(pat)}: found at symbol positions "
                      f"{positions[:6]}{'...' if len(positions) > 6 else ''}")

    # === STEP 4: Frame structure from autocorrelation ===
    print(f"\n{'='*60}")
    print(f"  STEP 4: FRAME STRUCTURE (AUTOCORRELATION)")
    print(f"{'='*60}")

    fig, axes = plt.subplots(len(param_sets), 1,
                              figsize=(16, 4 * len(param_sets)))
    for ax, (nsps, h, plabel) in zip(axes, param_sets):
        baud = sr / nsps
        key = (plabel, bursts[0][2])  # longest burst
        tones, good, freqs, powers = all_tones[key]

        autocorr = analyze_frame_structure(tones, powers, good, baud, plabel)

        ax.plot(autocorr[:min(200, len(autocorr))])
        ax.set_xlabel("Lag (symbols)")
        ax.set_ylabel("Autocorrelation")
        ax.set_title(f"{plabel} — Power Autocorrelation")
        ax.grid(True, alpha=0.3)
        ax.axhline(0, color="gray", linewidth=0.5)

    plt.tight_layout()
    plt.savefig("output_sync/autocorrelation.png", dpi=150)
    plt.close()
    print(f"\n  Saved: output_sync/autocorrelation.png")

    # === STEP 5: Exhaustive sync pattern search ===
    print(f"\n{'='*60}")
    print(f"  STEP 5: EXHAUSTIVE SYNC SEARCH")
    print(f"{'='*60}")

    for nsps, h, plabel in param_sets[:2]:  # top 2
        baud = sr / nsps
        for bt0, bt1, blabel in bursts[:1]:  # longest burst
            key = (plabel, blabel)
            tones, good = all_tones[key][:2]
            exhaustive_sync_search(tones, good, baud, f"{plabel} {blabel}")

    # === STEP 6: Test FT4 sync structure ===
    print(f"\n{'='*60}")
    print(f"  STEP 6: FT4-STYLE SYNC TEST")
    print(f"{'='*60}")

    for nsps, h, plabel in param_sets:
        baud = sr / nsps
        for bt0, bt1, blabel in bursts:
            key = (plabel, blabel)
            tones, good = all_tones[key][:2]
            test_ft4_sync(tones, good, baud, f"{plabel} {blabel}")

    # === STEP 7: Visualize tone sequences ===
    print(f"\n{'='*60}")
    print(f"  STEP 7: TONE SEQUENCE PLOTS")
    print(f"{'='*60}")

    for nsps, h, plabel in param_sets[:2]:
        baud = sr / nsps
        fig, axes = plt.subplots(len(bursts), 1, figsize=(20, 4 * len(bursts)))
        for ax, (bt0, bt1, blabel) in zip(axes, bursts):
            key = (plabel, blabel)
            tones, good, freqs, powers = all_tones[key]

            x = np.arange(len(tones))
            ax.plot(x[good], tones[good], "b.-", markersize=4, linewidth=0.8)
            ax.plot(x[~good], tones[~good], "r.", markersize=2, alpha=0.3)
            ax.set_ylabel("Tone index")
            ax.set_title(f"{plabel} — {blabel} ({np.sum(good)} good symbols)")
            ax.grid(True, alpha=0.3)

            # Mark every 79 symbols (FT8 frame length)
            for i in range(0, len(tones), 79):
                ax.axvline(i, color="green", alpha=0.3, linewidth=0.5)
            # Mark every 77 symbols (FT4 frame length)
            for i in range(0, len(tones), 77):
                ax.axvline(i, color="orange", alpha=0.3, linewidth=0.5,
                          linestyle="--")

        axes[-1].set_xlabel("Symbol index")
        plt.suptitle(f"Tone Sequences — {plabel}", fontsize=13)
        plt.tight_layout()
        plt.savefig(f"output_sync/tone_sequences_{plabel}.png", dpi=150)
        plt.close()
        print(f"  Saved: output_sync/tone_sequences_{plabel}.png")

    print(f"\n{'='*60}")
    print(f"  DONE — check output_sync/")
    print(f"{'='*60}")


if __name__ == "__main__":
    main()
