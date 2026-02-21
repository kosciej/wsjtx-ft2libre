#!/usr/bin/env python3
"""
Cross-reference using raw frequency comparison instead of quantized tones.

For sync symbols, the raw frequencies should be nearly identical between captures
(since they transmit the same sync tones). For data symbols, the frequencies
should generally differ (different messages → different data tones).

We compute |f1[i] - f2[i]| mod tone_spacing for each symbol position and look
for positions where the difference is small.
"""

import numpy as np
import soundfile as sf
import matplotlib.pyplot as plt
from pathlib import Path

OUTPUT_DIR = Path("output_xref")
OUTPUT_DIR.mkdir(exist_ok=True)

FS = 12000
NSPS = 360
H = 0.5
TONE_SPACING = H * FS / NSPS  # 16.667 Hz
NSYMBOLS = 72
NFFT = NSPS * 4


def extract_freqs(signal, start):
    freqs = np.zeros(NSYMBOLS)
    for i in range(NSYMBOLS):
        s = start + i * NSPS
        e = s + NSPS
        if e > len(signal):
            freqs[i] = np.nan
            continue
        seg = signal[s:e] * np.hanning(NSPS)
        spec = np.abs(np.fft.rfft(seg, n=NFFT))
        lo = int(200 * NFFT / FS)
        hi = int(2500 * NFFT / FS)
        pk = lo + np.argmax(spec[lo:hi])
        if 0 < pk < len(spec) - 1:
            a = np.log(spec[pk-1]+1e-30)
            b = np.log(spec[pk]+1e-30)
            c = np.log(spec[pk+1]+1e-30)
            d = a - 2*b + c
            delta = 0.5*(a-c)/d if abs(d) > 1e-10 else 0
            freqs[i] = (pk+delta)*FS/NFFT
        else:
            freqs[i] = pk*FS/NFFT
    return freqs


def detect_burst(signal, skip_first=False):
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


def main():
    sig1, _ = sf.read("ft2_capture.wav")
    if len(sig1.shape) > 1: sig1 = sig1[:, 0]
    sig2, _ = sf.read("ft2_capture2.wav")
    if len(sig2.shape) > 1: sig2 = sig2[:, 0]

    b1 = detect_burst(sig1) + int(0.05*FS)
    b2 = detect_burst(sig2, skip_first=True) + int(0.05*FS)

    # Scan: for each (off1, off2) pair, compute frequency difference metric
    print("Computing frequency-difference metric for all alignment pairs...")

    step = 2  # 2-sample steps for C1
    offsets_1 = list(range(-400, 400, step))
    offsets_2 = list(range(-400, 400, 1))  # 1-sample for C2

    # Pre-compute C1 frequencies at each offset
    c1_freqs = {}
    for off1 in offsets_1:
        start = b1 + off1
        if start < 0 or start + NSYMBOLS * NSPS > len(sig1):
            continue
        c1_freqs[off1] = extract_freqs(sig1, start)

    # For each C1 offset, sweep C2 and compute frequency difference
    best_result = None
    best_score = -999

    # Also collect per-offset2 best scores for plotting
    c2_match_curve = {}

    for off1 in c1_freqs:
        f1 = c1_freqs[off1]
        if np.any(np.isnan(f1)):
            continue

        for off2 in offsets_2:
            start2 = b2 + off2
            if start2 < 0 or start2 + NSYMBOLS * NSPS > len(sig2):
                continue

            # Only compute for a few C1 offsets near the best RMS region
            if off1 < -380 or off1 > -300:
                continue

            f2 = extract_freqs(sig2, start2)
            if np.any(np.isnan(f2)):
                continue

            # Frequency difference modulo tone_spacing
            # For sync symbols, diff should be near 0 (mod tone_spacing)
            # For data symbols, diff should be random
            diff = f1 - f2

            # Compute "residual" relative to nearest multiple of tone_spacing
            # A small residual means f1 and f2 are on the same tone
            residual = diff / TONE_SPACING
            residual_frac = residual - np.round(residual)
            abs_frac = np.abs(residual_frac)

            # Score: number of positions with small fractional residual
            # A position is "matching" if |residual_frac| < threshold
            threshold = 0.2
            n_match = np.sum(abs_frac < threshold)

            # Also check: matching positions should have residual near integer
            # (i.e., same tone number), not just same modular class
            same_tone = np.abs(diff) < TONE_SPACING * 0.3
            n_same_tone = np.sum(same_tone)

            if n_same_tone > best_score:
                best_score = n_same_tone
                best_result = {
                    'off1': off1, 'off2': off2,
                    'f1': f1.copy(), 'f2': f2.copy(),
                    'diff': diff.copy(),
                    'same_tone': same_tone.copy(),
                    'n_same': n_same_tone,
                    'n_mod_match': n_match,
                }

    if best_result is None:
        print("No valid alignment found!")
        return

    r = best_result
    print(f"\nBest alignment: off1={r['off1']:+d} off2={r['off2']:+d}")
    print(f"  Same-tone positions: {r['n_same']}/{NSYMBOLS}")
    print(f"  Mod-match positions: {r['n_mod_match']}/{NSYMBOLS}")

    # Now do a FINER search around this best alignment
    print(f"\nRefining around best alignment...")
    fine_best = best_result
    for off1 in range(r['off1']-10, r['off1']+11, 1):
        start1 = b1 + off1
        if start1 < 0 or start1 + NSYMBOLS * NSPS > len(sig1):
            continue
        f1 = extract_freqs(sig1, start1)
        if np.any(np.isnan(f1)):
            continue

        for off2 in range(r['off2']-20, r['off2']+21, 1):
            start2 = b2 + off2
            if start2 < 0 or start2 + NSYMBOLS * NSPS > len(sig2):
                continue
            f2 = extract_freqs(sig2, start2)
            if np.any(np.isnan(f2)):
                continue

            diff = f1 - f2
            same_tone = np.abs(diff) < TONE_SPACING * 0.3
            n_same = np.sum(same_tone)

            if n_same > fine_best['n_same']:
                fine_best = {
                    'off1': off1, 'off2': off2,
                    'f1': f1.copy(), 'f2': f2.copy(),
                    'diff': diff.copy(),
                    'same_tone': same_tone.copy(),
                    'n_same': n_same,
                }

    r = fine_best
    print(f"  Refined: off1={r['off1']:+d} off2={r['off2']:+d}, same={r['n_same']}")

    # Quantize using a COMMON base frequency
    all_freqs = np.concatenate([r['f1'], r['f2']])
    global_base = np.min(all_freqs)
    best_base = global_base
    best_rms = 999
    for base in np.arange(global_base - 2*TONE_SPACING, global_base + TONE_SPACING, 0.2):
        tf = (all_freqs - base) / TONE_SPACING
        ti = np.round(tf).astype(int)
        valid = (ti >= 0) & (ti <= 7)
        if np.sum(valid) > len(all_freqs) * 0.9:
            rms = np.sqrt(np.mean((tf[valid] - ti[valid])**2))
            if rms < best_rms:
                best_rms = rms
                best_base = base

    print(f"\n  Common base: {best_base:.1f} Hz (RMS: {best_rms:.4f})")

    t1 = np.round((r['f1'] - best_base) / TONE_SPACING).astype(int)
    t2 = np.round((r['f2'] - best_base) / TONE_SPACING).astype(int)

    # Classify
    v1 = (t1 >= 0) & (t1 <= 7)
    v2 = (t2 >= 0) & (t2 <= 7)

    classification = []
    for i in range(NSYMBOLS):
        if not (v1[i] and v2[i]):
            classification.append('?')
        elif t1[i] == t2[i]:
            classification.append('S')
        else:
            classification.append('D')

    n_sync = classification.count('S')
    n_data = classification.count('D')
    n_unc = classification.count('?')

    print(f"\n  SYNC={n_sync} DATA={n_data} UNC={n_unc}")
    print(f"  Data bits: {n_data}×3 = {n_data*3}")

    # Frame map
    print(f"\n  Frame map:")
    for row in range(0, NSYMBOLS, 36):
        end = min(row + 36, NSYMBOLS)
        print(f"    Pos: {' '.join(f'{i:2d}' for i in range(row, end))}")
        print(f"    Map: {' '.join(f' {classification[i]}' for i in range(row, end))}")
        print(f"    C1:  {' '.join(f'{t1[i]:2d}' for i in range(row, end))}")
        print(f"    C2:  {' '.join(f'{t2[i]:2d}' for i in range(row, end))}")
        print()

    sync_pos = [i for i in range(NSYMBOLS) if classification[i] == 'S']
    print(f"  Sync positions: {sync_pos}")
    print(f"  Sync tones: {[int(t1[i]) for i in sync_pos]}")

    # Contiguous blocks
    print(f"\n  Contiguous sync blocks:")
    in_b = False
    bs = 0
    for i in range(NSYMBOLS + 1):
        if i < NSYMBOLS and classification[i] == 'S':
            if not in_b:
                bs = i
                in_b = True
        else:
            if in_b:
                bl = i - bs
                bt = t1[bs:i]
                norm = bt - np.min(bt)
                tag = ""
                if bl >= 4 and len(set(norm)) == len(norm):
                    tag = " (permutation)"
                    if bl == 7 and set(norm) == set(range(7)):
                        # Check Costas
                        diffs = set()
                        ok = True
                        for a in range(7):
                            for b in range(a+1, 7):
                                d = (int(norm[b])-int(norm[a]), b-a)
                                if d in diffs: ok = False
                                diffs.add(d)
                        tag = " COSTAS!" if ok else " (perm, not Costas)"
                print(f"    [{bs}:{i}] len={bl} tones={bt.tolist()} norm={norm.tolist()}{tag}")
                in_b = False

    # ================================================================
    # VISUALIZATION
    # ================================================================
    fig, axes = plt.subplots(4, 1, figsize=(22, 16))
    x = np.arange(NSYMBOLS)

    # 1. Frequency difference
    diff = r['f1'] - r['f2']
    colors_diff = ['green' if classification[i]=='S' else
                   ('orange' if classification[i]=='?' else 'red')
                   for i in range(NSYMBOLS)]
    axes[0].bar(x, diff, color=colors_diff, edgecolor='black', linewidth=0.3)
    for t in range(-4, 5):
        axes[0].axhline(t * TONE_SPACING, color='gray', alpha=0.2, linewidth=0.5)
    axes[0].axhline(0, color='black', linewidth=1)
    axes[0].set_ylabel("Freq difference (Hz)")
    axes[0].set_title("Frequency Difference (C1 - C2)\nGreen = sync (near 0), Red = data")

    # 2. Absolute frequency difference
    abs_diff = np.abs(diff)
    axes[1].bar(x, abs_diff, color=colors_diff, edgecolor='black', linewidth=0.3)
    axes[1].axhline(TONE_SPACING * 0.3, color='blue', linestyle='--',
                     label=f'Threshold ({TONE_SPACING*0.3:.1f} Hz)')
    axes[1].set_ylabel("|Freq diff| (Hz)")
    axes[1].set_title("Absolute Frequency Difference")
    axes[1].legend()

    # 3. Tones comparison
    w = 0.35
    c1c = ['green' if classification[i]=='S' else ('orange' if classification[i]=='?' else '#cc4444')
           for i in range(NSYMBOLS)]
    c2c = ['green' if classification[i]=='S' else ('orange' if classification[i]=='?' else '#4444cc')
           for i in range(NSYMBOLS)]
    axes[2].bar(x-w/2, np.clip(t1, -1, 8), w, color=c1c, alpha=0.8, label='C1')
    axes[2].bar(x+w/2, np.clip(t2, -1, 8), w, color=c2c, alpha=0.8, label='C2')
    axes[2].set_ylabel("Tone (0-7)")
    axes[2].legend()
    axes[2].set_ylim(-1.5, 8.5)

    # 4. Frame structure
    fc = ['green' if classification[i]=='S' else ('orange' if classification[i]=='?' else 'red')
          for i in range(NSYMBOLS)]
    axes[3].bar(x, [1]*NSYMBOLS, color=fc, edgecolor='black', linewidth=0.5)
    for i in range(NSYMBOLS):
        if classification[i] == 'S':
            axes[3].text(i, 0.5, str(int(t1[i])), ha='center', va='center',
                        fontsize=6, fontweight='bold', color='white')
    axes[3].set_title("Frame: Green=SYNC, Red=DATA, Orange=UNCERTAIN")
    axes[3].set_yticks([])
    axes[3].set_xlabel("Symbol Position")

    for ax in axes:
        ax.set_xlim(-0.5, NSYMBOLS-0.5)
        ax.set_xticks(range(0, NSYMBOLS, 2))

    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / "xref_freq.png", dpi=150)
    plt.close()
    print(f"\n  Saved xref_freq.png")


if __name__ == "__main__":
    main()
