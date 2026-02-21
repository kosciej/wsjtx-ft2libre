#!/usr/bin/env python3
"""
FT2 Verify Costas — focused test of FT8 Costas [3,1,4,0,6,5,2] with h=0.5
on the 72-symbol frame, assuming 2 Costas blocks.

With correct h=0.5 quantization, check if the Costas pattern emerges at
positions that divide 72 symbols into 2 sync + 2 data blocks.

Also test with the FT4 Costas pattern.

Usage:
    uv run python ft2_verify_costas.py ft2_capture.wav
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


def extract_and_quantize(data, sr, nsps, fmin, fmax, t_start, n_sym, h):
    """Extract and quantize symbols."""
    i0 = int(t_start * sr)
    seg = data[i0:]
    nfft = nsps * 8
    df = sr / nfft
    baud = sr / nsps
    spacing = h * baud
    n = min(len(seg) // nsps, n_sym)

    freqs = np.zeros(n)
    powers = np.zeros(n)
    snrs = np.zeros(n)

    for s in range(n):
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

    # Quantize
    good = snrs > 1.5
    if np.sum(good) < 5:
        return freqs, powers, snrs, np.zeros(n, dtype=int), good, 1e20

    gf = freqs[good]
    best_rms = 1e20
    best_base = gf.min()
    for base in np.arange(gf.min() - spacing, gf.min() + spacing, 0.2):
        q = np.round((gf - base) / spacing) * spacing + base
        rms = np.sqrt(np.mean((gf - q) ** 2))
        if rms < best_rms:
            best_rms = rms
            best_base = base

    tones = np.round((freqs - best_base) / spacing).astype(int)
    return freqs, powers, snrs, tones, good, best_rms


def check_costas_at_position(tones, position, costas, base_tone):
    """Check exact match of Costas array at a specific position with a specific base."""
    seg = tones[position:position + len(costas)]
    if len(seg) < len(costas):
        return 0, []
    matches = []
    score = 0
    for k in range(len(costas)):
        expected = base_tone + costas[k]
        actual = seg[k]
        match = actual == expected
        fuzzy = abs(actual - expected) <= 1
        if match:
            score += 1
        matches.append((k, actual, expected, match, fuzzy))
    return score, matches


def main():
    wav_path = sys.argv[1] if len(sys.argv) > 1 else "ft2_capture.wav"
    Path("output_verify").mkdir(exist_ok=True)

    data, sr = load_wav(wav_path)
    print(f"Loaded: {sr} Hz, {len(data)/sr:.2f}s")

    fmin, fmax = 1440, 1720
    nsps = 360
    baud = sr / nsps
    h = 0.5
    spacing = h * baud

    # TX boundary
    tx_on = 5.857  # From alignment scan

    # FT8 Costas array
    FT8_COSTAS = np.array([3, 1, 4, 0, 6, 5, 2])

    print(f"NSPS={nsps}, h={h}, baud={baud:.2f}, spacing={spacing:.2f} Hz")
    print(f"FT8 Costas: {FT8_COSTAS.tolist()}")

    # ================================================================
    # Try multiple start offsets (fine scan)
    # ================================================================
    print(f"\n{'='*60}")
    print(f"  FINE ALIGNMENT SCAN WITH COSTAS DETECTION")
    print(f"{'='*60}")

    # Scan start time in 1ms steps
    offsets = np.arange(-60, 60, 0.5) / 1000.0
    best_costas_score = 0
    best_config = None

    # For each offset, extract 72 symbols and check Costas at all possible positions
    for offset in offsets:
        t_start = tx_on + offset
        if t_start < 0:
            continue

        freqs, powers, snrs, tones, good, rms = extract_and_quantize(
            data, sr, nsps, fmin, fmax, t_start, 80, h)

        if len(tones) < 72:
            continue

        # Check Costas at all position pairs within 72 symbols
        for pos1 in range(0, 72 - 7):
            for pos2 in range(pos1 + 14, 72 - 6):
                # Try different base tones
                for base in range(-1, 4):
                    score1, _ = check_costas_at_position(tones, pos1, FT8_COSTAS, base)
                    score2, _ = check_costas_at_position(tones, pos2, FT8_COSTAS, base)
                    total = score1 + score2

                    if total > best_costas_score:
                        best_costas_score = total
                        best_config = {
                            'offset': offset, 'pos1': pos1, 'pos2': pos2,
                            'base': base, 'score': total,
                            'score1': score1, 'score2': score2,
                            'tones': tones[:80].copy(), 'rms': rms,
                        }

    if best_config:
        c = best_config
        print(f"\n  Best Costas match: {c['score']}/14")
        print(f"    Offset: {c['offset']*1000:.1f}ms")
        print(f"    Positions: [{c['pos1']}, {c['pos2']}]")
        print(f"    Base tone: {c['base']}")
        print(f"    Score1: {c['score1']}/7, Score2: {c['score2']}/7")
        print(f"    RMS: {c['rms']:.2f} Hz")

        # Show the actual vs expected at best positions
        for pos, label in [(c['pos1'], "Block1"), (c['pos2'], "Block2")]:
            seg = c['tones'][pos:pos + 7]
            expected = c['base'] + FT8_COSTAS
            print(f"\n    {label} at position {pos}:")
            print(f"      Expected: {expected.tolist()}")
            print(f"      Actual:   {seg.tolist()}")
            print(f"      Match:    {['✓' if a==e else '✗' for a,e in zip(seg, expected)]}")

        # Show full tone sequence with sync positions marked
        print(f"\n    Full 72-symbol tone sequence:")
        tones72 = c['tones'][:72]
        for s in range(72):
            in_sync1 = c['pos1'] <= s < c['pos1'] + 7
            in_sync2 = c['pos2'] <= s < c['pos2'] + 7
            marker = "S1" if in_sync1 else "S2" if in_sync2 else "  "
            print(f"      [{s:2d}] tone={tones72[s]:2d}  {marker}")

    # ================================================================
    # Also check: maybe it's 79 symbols and we're missing the edges
    # ================================================================
    print(f"\n{'='*60}")
    print(f"  CHECK: 79-SYMBOL FRAME (wider window)")
    print(f"{'='*60}")

    # Start earlier to capture potential first Costas
    for early_start in [tx_on - 0.21, tx_on - 0.15, tx_on - 0.10, tx_on - 0.05, tx_on]:
        freqs, powers, snrs, tones, good, rms = extract_and_quantize(
            data, sr, nsps, fmin, fmax, early_start, 85, h)

        n_good_first7 = int(np.sum(snrs[:7] > 2.0))
        n_good_total = int(np.sum(snrs[:79] > 2.0)) if len(snrs) >= 79 else 0

        # Check 3-Costas at standard FT8 positions
        if len(tones) >= 79:
            total_score = 0
            for pos in [0, 36, 72]:
                for base in range(-1, 4):
                    score, _ = check_costas_at_position(tones, pos, FT8_COSTAS, base)
                    total_score = max(total_score, score)

            # Check 3-Costas with best-base search
            best_3costas = 0
            best_3base = 0
            for base in range(-2, 5):
                s1, _ = check_costas_at_position(tones, 0, FT8_COSTAS, base)
                s2, _ = check_costas_at_position(tones, 36, FT8_COSTAS, base)
                s3, _ = check_costas_at_position(tones, 72, FT8_COSTAS, base)
                total = s1 + s2 + s3
                if total > best_3costas:
                    best_3costas = total
                    best_3base = base

            print(f"  Start={early_start:.3f}s: "
                  f"first7_good={n_good_first7}, total_good={n_good_total}, "
                  f"3-Costas best={best_3costas}/21 (base={best_3base})")

            if best_3costas >= 10:
                for pos in [0, 36, 72]:
                    seg = tones[pos:pos + 7]
                    expected = best_3base + FT8_COSTAS
                    matches = sum(1 for a, e in zip(seg, expected) if a == e)
                    print(f"    pos={pos}: expected={expected.tolist()} "
                          f"actual={seg.tolist()} matches={matches}/7")

    # ================================================================
    # Frequency-domain Costas correlation (more robust)
    # ================================================================
    print(f"\n{'='*60}")
    print(f"  FREQUENCY-DOMAIN COSTAS CORRELATION")
    print(f"{'='*60}")

    # Compute symbol spectra for the middle burst with wider window
    t_start = tx_on - 0.25
    i0 = int(t_start * sr)
    seg = data[i0:i0 + 85 * nsps]
    nfft = 2 * nsps  # Standard for sync detection
    df = sr / nfft
    nfos = nfft // nsps  # = 2
    tone_bins = nfos * h  # = 1.0 bin per tone for h=0.5

    n_sym = len(seg) // nsps
    spectra = np.zeros((n_sym, nfft // 2))

    for s in range(n_sym):
        chunk = seg[s * nsps:(s + 1) * nsps]
        w = np.hanning(nsps)
        padded = np.zeros(nfft)
        padded[:nsps] = chunk * w
        spectra[s] = np.abs(fft(padded))[:nfft // 2] ** 2

    freqs_axis = np.arange(nfft // 2) * df
    f_lo = int(fmin / df)
    f_hi = min(int(fmax / df) + 1, nfft // 2)

    print(f"  nfft={nfft}, df={df:.2f} Hz, tone_bins={tone_bins:.1f}")
    print(f"  {n_sym} symbol spectra computed")

    # 2-Costas search (positions within 72 symbols)
    best_snr = 0
    best_freq_config = None

    for pos1 in [0]:  # Costas at start
        for pos2 in range(25, 66):  # Second Costas somewhere in the middle
            for f0_bin in range(f_lo, f_hi - 8):
                score = 0
                for pos in [pos1, pos2]:
                    for k in range(7):
                        sym_idx = pos + k
                        if sym_idx >= n_sym:
                            continue
                        tone_bin = f0_bin + int(round(FT8_COSTAS[k] * tone_bins))
                        if 0 <= tone_bin < nfft // 2:
                            score += spectra[sym_idx, tone_bin]

                # Normalize by noise
                noise = np.median(spectra[pos1:pos2 + 7, f0_bin:f0_bin + 8])
                snr = score / (14 * noise + 1e-20)

                if snr > best_snr:
                    best_snr = snr
                    best_freq_config = (pos1, pos2, f0_bin, f0_bin * df)

    if best_freq_config:
        p1, p2, fb, f0 = best_freq_config
        print(f"\n  2-Costas freq-domain: SNR={best_snr:.2f}")
        print(f"    Positions: [{p1}, {p2}], base freq bin={fb} ({f0:.1f} Hz)")
        print(f"    Data symbols: {p2 - p1 - 7} between sync blocks")

    # 3-Costas search (79 symbols, extended window)
    best_snr_3 = 0
    best_freq_config_3 = None

    for f0_bin in range(f_lo, f_hi - 8):
        score = 0
        for pos in [0, 36, 72]:
            for k in range(7):
                sym_idx = pos + k
                if sym_idx >= n_sym:
                    continue
                tone_bin = f0_bin + int(round(FT8_COSTAS[k] * tone_bins))
                if 0 <= tone_bin < nfft // 2:
                    score += spectra[sym_idx, tone_bin]

        noise = np.median(spectra[:79, f0_bin:f0_bin + 8]) if n_sym >= 79 else 1
        snr = score / (21 * noise + 1e-20)

        if snr > best_snr_3:
            best_snr_3 = snr
            best_freq_config_3 = (f0_bin, f0_bin * df)

    if best_freq_config_3:
        fb, f0 = best_freq_config_3
        print(f"\n  3-Costas freq-domain: SNR={best_snr_3:.2f}")
        print(f"    Base freq bin={fb} ({f0:.1f} Hz)")

    # ================================================================
    # Print definitive tone sequence
    # ================================================================
    print(f"\n{'='*60}")
    print(f"  DEFINITIVE TONE SEQUENCE (NSPS=360, h=0.5)")
    print(f"{'='*60}")

    t_start = tx_on
    freqs, powers, snrs, tones, good, rms = extract_and_quantize(
        data, sr, nsps, fmin, fmax, t_start, 72, h)

    print(f"  Start: {t_start:.4f}s")
    print(f"  RMS quantization error: {rms:.2f} Hz")
    print(f"  Good symbols: {int(np.sum(good))}/{len(tones)}")
    print()

    # Find the base tone (minimum tone in the sequence)
    gt = tones[good]
    base = gt.min()
    print(f"  Base tone: {base}")
    print(f"  Tone range: {gt.min()}-{gt.max()} ({gt.max()-gt.min()+1} levels)")

    # Normalize so tones start at 0
    tones_norm = tones - base
    print(f"\n  Normalized tone sequence (0-indexed):")
    for s in range(len(tones)):
        snr_str = f"{snrs[s]:6.1f}" if snrs[s] > 0 else "  ---"
        marker = "***" if snrs[s] > 10 else "** " if snrs[s] > 3 else "*  " if snrs[s] > 2 else "   "
        print(f"  [{s:2d}] tone={tones_norm[s]:2d}  (raw={tones[s]:2d})  "
              f"f={freqs[s]:7.1f}  snr={snr_str} {marker}")

    print(f"\n  Compact: {tones_norm[:72].tolist()}")

    print(f"\n{'='*60}")
    print(f"  DONE")
    print(f"{'='*60}")


if __name__ == "__main__":
    main()
