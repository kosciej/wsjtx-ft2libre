#!/usr/bin/env python3
"""
Cross-reference v5: Fine-grained alignment with fractional NSPS support.
Fix capture 1 alignment at best RMS, then sweep capture 2 with 1-sample steps.
Try NSPS=359 and NSPS=360.
"""

import numpy as np
import soundfile as sf
import matplotlib.pyplot as plt
from pathlib import Path
from collections import Counter

OUTPUT_DIR = Path("output_xref")
OUTPUT_DIR.mkdir(exist_ok=True)

FS = 12000
NSYMBOLS = 72


def extract_and_quantize(signal, start, nsps, h=0.5):
    """Extract tones with given NSPS and h."""
    tone_sp = h * FS / nsps
    nfft = int(nsps) * 4
    freqs = np.zeros(NSYMBOLS)

    for i in range(NSYMBOLS):
        s = start + int(round(i * nsps))
        e = s + int(nsps)
        if e > len(signal):
            freqs[i] = np.nan
            continue
        seg_len = e - s
        seg = signal[s:e] * np.hanning(seg_len)
        spec = np.abs(np.fft.rfft(seg, n=nfft))
        lo = int(200 * nfft / FS)
        hi = int(2500 * nfft / FS)
        pk = lo + np.argmax(spec[lo:hi])
        if 0 < pk < len(spec) - 1:
            a = np.log(spec[pk-1]+1e-30)
            b = np.log(spec[pk]+1e-30)
            c = np.log(spec[pk+1]+1e-30)
            d = a - 2*b + c
            delta = 0.5*(a-c)/d if abs(d) > 1e-10 else 0
            freqs[i] = (pk+delta)*FS/nfft
        else:
            freqs[i] = pk*FS/nfft

    # Quantize
    valid = ~np.isnan(freqs)
    f = freqs[valid]
    best_rms = 999
    best_base = 0
    best_tones = None

    for base in np.arange(np.min(f) - 2*tone_sp, np.min(f) + tone_sp, 0.3):
        tf = (f - base) / tone_sp
        ti = np.round(tf).astype(int)
        in_range = (ti >= 0) & (ti <= 7)
        if np.sum(~in_range) > 4:
            continue
        if np.sum(in_range) > 0:
            rms = np.sqrt(np.mean((tf[in_range] - ti[in_range])**2))
            if rms < best_rms:
                best_rms = rms
                best_base = base
                tones = np.full(NSYMBOLS, -99, dtype=int)
                tones[valid] = ti
                best_tones = tones

    return best_tones, best_rms, best_base, freqs


def main():
    sig1, _ = sf.read("ft2_capture.wav")
    if len(sig1.shape) > 1: sig1 = sig1[:, 0]
    sig2, _ = sf.read("ft2_capture2.wav")
    if len(sig2.shape) > 1: sig2 = sig2[:, 0]

    # Burst detection
    def detect(signal, skip_first=False):
        wl = int(0.05 * FS)
        pw = np.convolve(signal**2, np.ones(wl)/wl, mode='same')
        th = np.max(pw) * 0.01
        ab = pw > th
        tr = np.diff(ab.astype(int))
        starts = np.where(tr == 1)[0]
        ends = np.where(tr == -1)[0]
        if ab[0]: starts = np.concatenate([[0], starts])
        if ab[-1]: ends = np.concatenate([ends, [len(signal)-1]])
        n = min(len(starts), len(ends))
        bursts = [(starts[i], ends[i]) for i in range(n)
                  if (ends[i]-starts[i])/FS > 1.5]
        if skip_first and len(bursts) > 1:
            return bursts[1][0]
        return max(bursts, key=lambda b: b[1]-b[0])[0] if bursts else None

    b1 = detect(sig1) + int(0.05*FS)
    b2 = detect(sig2, skip_first=True) + int(0.05*FS)
    print(f"Burst centers: C1={b1/FS:.4f}s, C2={b2/FS:.4f}s")

    # For each NSPS value, fix C1 at best alignment, sweep C2
    for nsps in [359, 360]:
        print(f"\n{'='*70}")
        print(f"NSPS = {nsps} (baud = {FS/nsps:.4f} Hz)")
        print(f"{'='*70}")

        # Find best C1 alignment
        best_c1_rms = 999
        best_c1_off = 0
        best_c1_tones = None
        best_c1_base = 0

        for off in range(-400, 400, 2):
            start = b1 + off
            if start < 0: continue
            tones, rms, base, _ = extract_and_quantize(sig1, start, nsps)
            if tones is not None and rms < best_c1_rms:
                best_c1_rms = rms
                best_c1_off = off
                best_c1_tones = tones
                best_c1_base = base

        print(f"\nCapture 1: best off={best_c1_off} RMS={best_c1_rms:.4f} base={best_c1_base:.1f}")
        print(f"  Tones: {best_c1_tones.tolist()}")

        # Sweep C2 with 1-sample steps, compute match count
        c1_valid = (best_c1_tones >= 0) & (best_c1_tones <= 7)
        match_curve = []

        for off in range(-400, 400, 1):
            start = b2 + off
            if start < 0: continue
            tones, rms, base, _ = extract_and_quantize(sig2, start, nsps)
            if tones is None:
                match_curve.append((off, 0, 999, None))
                continue

            c2_valid = (tones >= 0) & (tones <= 7)
            both_valid = c1_valid & c2_valid
            matches = np.sum((best_c1_tones == tones) & both_valid)
            match_curve.append((off, matches, rms, tones))

        # Find top match alignments
        match_curve.sort(key=lambda x: -x[1])
        print(f"\nTop 10 C2 alignments by match count:")
        for off, matches, rms, tones in match_curve[:10]:
            print(f"  off={off:+4d}: matches={matches:2d} rms={rms:.4f}")

        # Use the best match
        best_off2, best_matches, best_rms2, best_c2_tones = match_curve[0]
        print(f"\nBest: C2 off={best_off2:+d}, matches={best_matches}, rms={best_rms2:.4f}")
        print(f"  Tones: {best_c2_tones.tolist()}")

        # Cross-reference
        t1 = best_c1_tones
        t2 = best_c2_tones
        v1 = (t1 >= 0) & (t1 <= 7)
        v2 = (t2 >= 0) & (t2 <= 7)
        both_v = v1 & v2

        classification = []
        for i in range(NSYMBOLS):
            if not both_v[i]:
                classification.append('?')
            elif t1[i] == t2[i]:
                classification.append('S')
            else:
                classification.append('D')

        n_sync = classification.count('S')
        n_data = classification.count('D')
        n_unc = classification.count('?')

        print(f"\nSYNC={n_sync} DATA={n_data} UNC={n_unc}")
        print(f"Data bits: {n_data}×3 = {n_data*3}")

        # Frame map
        for row in range(0, NSYMBOLS, 36):
            end = min(row + 36, NSYMBOLS)
            print(f"\n  Pos: {' '.join(f'{i:2d}' for i in range(row, end))}")
            print(f"  Map: {' '.join(f' {classification[i]}' for i in range(row, end))}")
            print(f"  C1:  {' '.join(f'{t1[i]:2d}' for i in range(row, end))}")
            print(f"  C2:  {' '.join(f'{t2[i]:2d}' for i in range(row, end))}")

        sync_pos = [i for i in range(NSYMBOLS) if classification[i] == 'S']
        print(f"\nSync positions: {sync_pos}")
        print(f"Sync tones: {[int(t1[i]) for i in sync_pos]}")

        # Look for contiguous blocks
        print(f"\nContiguous sync blocks:")
        blocks = []
        in_b = False
        bs_start = 0
        for i in range(NSYMBOLS + 1):
            if i < NSYMBOLS and classification[i] == 'S':
                if not in_b:
                    bs_start = i
                    in_b = True
            else:
                if in_b:
                    bl = i - bs_start
                    bt = t1[bs_start:i]
                    blocks.append((bs_start, i, bl, bt))
                    norm = bt - np.min(bt)
                    in_b = False
                    if bl >= 3:
                        is_perm = len(set(norm)) == len(norm)
                        print(f"  [{bs_start}:{i}] len={bl} "
                              f"tones={bt.tolist()} norm={norm.tolist()}"
                              f"{' (permutation)' if is_perm else ''}")

        # Best 2×7 pattern search
        print(f"\nBest 2×7 sync block positions:")
        best_score = 0
        best_2x7 = None
        for s1 in range(NSYMBOLS - 13):
            for s2 in range(s1+7, NSYMBOLS - 6):
                positions = list(range(s1, s1+7)) + list(range(s2, s2+7))
                sc = sum(1 for p in positions if classification[p] == 'S')
                dw = sum(1 for p in positions if classification[p] == 'D')
                score = sc - 2*dw  # Penalize data in sync positions
                if score > best_score:
                    best_score = score
                    best_2x7 = (s1, s2, sc, dw)

        if best_2x7:
            s1, s2, sc, dw = best_2x7
            b1t = t1[s1:s1+7]
            b2t = t1[s2:s2+7]
            n1 = b1t - np.min(b1t)
            n2 = b2t - np.min(b2t)
            nd = NSYMBOLS - 14
            print(f"  [{s1}:{s1+7}] + [{s2}:{s2+7}]")
            print(f"  Sync correct: {sc}/14, Data in sync: {dw}")
            print(f"  Data capacity: {nd}×3 = {nd*3} bits")
            print(f"  Block1: tones={b1t.tolist()} norm={n1.tolist()}")
            print(f"  Block2: tones={b2t.tolist()} norm={n2.tolist()}")
            # Check if blocks are identical
            if np.array_equal(n1, n2):
                print(f"  -> Blocks have SAME normalized pattern!")
            # Check Costas
            for bi, bn in [(1, n1), (2, n2)]:
                if len(set(bn)) == 7 and set(bn) == set(range(7)):
                    diffs = set()
                    is_costas = True
                    for a in range(7):
                        for b in range(a+1, 7):
                            d = (int(bn[b])-int(bn[a]), b-a)
                            if d in diffs:
                                is_costas = False
                            diffs.add(d)
                    print(f"  Block{bi}: {'COSTAS ARRAY' if is_costas else 'not Costas'}")

        # Plot match count curve
        offs_plot = [x[0] for x in sorted(match_curve, key=lambda x: x[0])
                     if x[3] is not None]
        matches_plot = [x[1] for x in sorted(match_curve, key=lambda x: x[0])
                        if x[3] is not None]

        plt.figure(figsize=(14, 4))
        plt.plot(offs_plot, matches_plot, 'b-', linewidth=0.5)
        plt.axhline(NSYMBOLS/8, color='r', linestyle='--', alpha=0.5,
                     label=f'Random chance ({NSYMBOLS/8:.1f})')
        plt.xlabel("Capture 2 offset (samples)")
        plt.ylabel("Matching symbols")
        plt.title(f"Match count vs C2 alignment offset (NSPS={nsps})")
        plt.legend()
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.savefig(OUTPUT_DIR / f"match_curve_nsps{nsps}.png", dpi=150)
        plt.close()

    print("\nDone!")


if __name__ == "__main__":
    main()
