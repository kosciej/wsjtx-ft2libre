#!/usr/bin/env python3
"""
Cross-reference v3: Outlier-tolerant tone extraction.
Allows a few symbols to fall outside 0-7 range and marks them as uncertain.
"""

import numpy as np
import soundfile as sf
import matplotlib.pyplot as plt
from pathlib import Path

OUTPUT_DIR = Path("output_xref")
OUTPUT_DIR.mkdir(exist_ok=True)

FS = 12000
NSPS = 360
BAUD = FS / NSPS
H = 0.5
TONE_SPACING = H * BAUD  # 16.667 Hz
NSYMBOLS = 72
NFFT = NSPS * 4


def extract_freqs(signal, start, nsymbols=NSYMBOLS, nsps=NSPS, nfft=NFFT):
    freqs = np.zeros(nsymbols)
    for i in range(nsymbols):
        s = start + i * nsps
        e = s + nsps
        seg = signal[s:e] * np.hanning(nsps)
        spec = np.abs(np.fft.rfft(seg, n=nfft))
        lo = int(200 * nfft / FS)
        hi = int(2500 * nfft / FS)
        pk = lo + np.argmax(spec[lo:hi])
        if 0 < pk < len(spec) - 1:
            a, b, c = np.log(spec[pk-1]+1e-30), np.log(spec[pk]+1e-30), np.log(spec[pk+1]+1e-30)
            d = a - 2*b + c
            delta = 0.5*(a-c)/d if abs(d) > 1e-10 else 0
            freqs[i] = (pk+delta)*FS/nfft
        else:
            freqs[i] = pk*FS/nfft
    return freqs


def quantize_with_outliers(freqs, tone_spacing=TONE_SPACING, max_outliers=5):
    """Quantize to 0-7, allowing some outlier symbols."""
    best_rms = 999
    best_base = 0
    best_tones = None

    for base in np.arange(np.min(freqs) - 2*tone_spacing,
                          np.min(freqs) + tone_spacing, 0.2):
        tf = (freqs - base) / tone_spacing
        ti = np.round(tf).astype(int)
        valid = (ti >= 0) & (ti <= 7)
        n_outliers = np.sum(~valid)

        if n_outliers > max_outliers:
            continue
        if np.sum(valid) < NSYMBOLS - max_outliers:
            continue

        rms = np.sqrt(np.mean((tf[valid] - ti[valid])**2))
        if rms < best_rms:
            best_rms = rms
            best_base = base
            best_tones = ti.copy()

    return best_tones, best_rms, best_base


def find_best_alignment(signal, burst_start_sample, label=""):
    """Scan offsets to find best symbol alignment."""
    margin = int(0.05 * FS)
    center = burst_start_sample + margin

    best_rms = 999
    best_result = None

    for offset in range(-NSPS, NSPS + 1, NSPS // 36):
        start = center + offset
        if start < 0 or start + NSYMBOLS * NSPS > len(signal):
            continue

        freqs = extract_freqs(signal, start)
        tones, rms, base = quantize_with_outliers(freqs)

        if tones is not None and rms < best_rms:
            best_rms = rms
            valid = (tones >= 0) & (tones <= 7)
            outliers = np.where(~valid)[0]
            best_result = {
                'tones': tones,
                'freqs': freqs,
                'rms': rms,
                'base': base,
                'offset': offset,
                'start': start,
                'outlier_pos': outliers,
            }

    if best_result:
        r = best_result
        valid_mask = (r['tones'] >= 0) & (r['tones'] <= 7)
        unique = len(np.unique(r['tones'][valid_mask]))
        print(f"  {label} best alignment:")
        print(f"    Offset: {r['offset']} ({r['offset']/FS*1000:.1f}ms)")
        print(f"    RMS: {r['rms']:.4f}, Base: {r['base']:.1f} Hz, Unique: {unique}")
        print(f"    Outliers ({len(r['outlier_pos'])}): {r['outlier_pos'].tolist()}")
        print(f"    Tones: {r['tones'].tolist()}")
    else:
        print(f"  {label}: FAILED")

    return best_result


def main():
    # Load both captures
    sig1, _ = sf.read("ft2_capture.wav")
    if len(sig1.shape) > 1: sig1 = sig1[:, 0]
    sig2, _ = sf.read("ft2_capture2.wav")
    if len(sig2.shape) > 1: sig2 = sig2[:, 0]

    print(f"Capture 1: {len(sig1)/FS:.2f}s")
    print(f"Capture 2: {len(sig2)/FS:.2f}s")

    # Detect bursts for both captures
    def find_middle_burst(signal, skip_first=False):
        window_len = int(0.05 * FS)
        power = np.convolve(signal**2, np.ones(window_len)/window_len, mode='same')
        threshold = np.max(power) * 0.01
        above = power > threshold
        trans = np.diff(above.astype(int))
        starts = np.where(trans == 1)[0]
        ends = np.where(trans == -1)[0]
        if above[0]:
            starts = np.concatenate([[0], starts])
        if above[-1]:
            ends = np.concatenate([ends, [len(signal)-1]])
        n = min(len(starts), len(ends))
        bursts = [(starts[i], ends[i]) for i in range(n)
                  if (ends[i]-starts[i])/FS > 1.5]

        for i, (s, e) in enumerate(bursts):
            tag = " [SKIP]" if (skip_first and i == 0) else ""
            print(f"    Burst {i}: {s/FS:.3f}-{e/FS:.3f}s ({(e-s)/FS:.3f}s){tag}")

        if skip_first and len(bursts) > 1:
            return bursts[1][0]
        elif bursts:
            longest = max(bursts, key=lambda b: b[1]-b[0])
            return longest[0]
        return None

    print("\nCapture 1 bursts:")
    burst1_start = find_middle_burst(sig1, skip_first=False)
    print("\nCapture 2 bursts:")
    burst2_start = find_middle_burst(sig2, skip_first=True)

    if burst1_start is None or burst2_start is None:
        print("Cannot find bursts!")
        return

    # Extract tones
    print(f"\n{'='*70}")
    print("TONE EXTRACTION")
    print(f"{'='*70}")

    r1 = find_best_alignment(sig1, burst1_start, "Capture 1")
    r2 = find_best_alignment(sig2, burst2_start, "Capture 2")

    if r1 is None or r2 is None:
        print("\nExtraction failed!")
        return

    # ================================================================
    # CROSS-REFERENCE
    # ================================================================
    t1 = r1['tones']
    t2 = r2['tones']
    outlier1 = set(r1['outlier_pos'])
    outlier2 = set(r2['outlier_pos'])
    all_outliers = outlier1 | outlier2

    print(f"\n{'='*70}")
    print("CROSS-REFERENCE")
    print(f"{'='*70}")

    # Classification:
    # - SYNC: same tone in both captures (both valid)
    # - DATA: different tone in both captures (both valid)
    # - UNCERTAIN: one or both are outliers
    classification = []
    for i in range(NSYMBOLS):
        if i in all_outliers:
            classification.append('?')
        elif t1[i] == t2[i]:
            classification.append('S')  # sync
        else:
            classification.append('D')  # data

    n_sync = classification.count('S')
    n_data = classification.count('D')
    n_uncertain = classification.count('?')

    print(f"\nSYNC: {n_sync}, DATA: {n_data}, UNCERTAIN: {n_uncertain}")
    print(f"Data capacity: {n_data} × 3 = {n_data * 3} bits (need 174 for LDPC)")

    # Print detailed frame map
    print(f"\nDetailed frame map:")
    for row in range(0, NSYMBOLS, 36):
        end = min(row + 36, NSYMBOLS)
        print(f"  Pos:  {' '.join(f'{i:2d}' for i in range(row, end))}")
        print(f"  Map:  {' '.join(f' {classification[i]}' for i in range(row, end))}")
        print(f"  C1:   {' '.join(f'{t1[i]:2d}' for i in range(row, end))}")
        print(f"  C2:   {' '.join(f'{t2[i]:2d}' for i in range(row, end))}")
        print()

    # Sync positions and values
    sync_pos = [i for i in range(NSYMBOLS) if classification[i] == 'S']
    data_pos = [i for i in range(NSYMBOLS) if classification[i] == 'D']
    uncertain_pos = [i for i in range(NSYMBOLS) if classification[i] == '?']

    print(f"SYNC positions: {sync_pos}")
    print(f"Sync tone values: {[t1[i] for i in sync_pos]}")
    print(f"\nDATA positions: {data_pos}")
    print(f"UNCERTAIN positions: {uncertain_pos}")

    # Look for contiguous sync blocks
    print(f"\nContiguous sync blocks:")
    blocks = []
    in_block = False
    block_start = 0
    for i in range(NSYMBOLS + 1):
        if i < NSYMBOLS and classification[i] == 'S':
            if not in_block:
                block_start = i
                in_block = True
        else:
            if in_block:
                block_len = i - block_start
                block_tones = t1[block_start:i]
                blocks.append((block_start, i, block_len, block_tones))
                norm = block_tones - np.min(block_tones)
                is_perm = (len(set(norm)) == len(norm) and
                           set(norm) == set(range(len(norm))))
                print(f"  [{block_start}:{i}] len={block_len} "
                      f"tones={block_tones.tolist()} "
                      f"norm={norm.tolist()}"
                      f"{' PERMUTATION!' if is_perm else ''}")
                in_block = False

    # Check if uncertain positions could extend sync blocks
    print(f"\nChecking if uncertain positions extend sync blocks...")
    for pos in uncertain_pos:
        neighbors_sync = []
        if pos > 0 and classification[pos-1] == 'S':
            neighbors_sync.append(pos-1)
        if pos < NSYMBOLS-1 and classification[pos+1] == 'S':
            neighbors_sync.append(pos+1)
        if neighbors_sync:
            print(f"  Position {pos}: adjacent to sync at {neighbors_sync}")
            print(f"    C1 tone: {t1[pos]}, C2 tone: {t2[pos]}")
            # If we treat uncertain as sync, what values would match?
            # The outlier tone value might still be meaningful
            if t1[pos] in range(0, 8) and t2[pos] in range(0, 8):
                print(f"    Both in range but differ: t1={t1[pos]}, t2={t2[pos]}")
            else:
                print(f"    At least one out of range")

    # ================================================================
    # COSTAS CHECK
    # ================================================================
    print(f"\n{'='*70}")
    print("COSTAS ARRAY CHECK ON SYNC BLOCKS")
    print(f"{'='*70}")

    for start, end, length, tones in blocks:
        if length >= 4:
            norm = tones - np.min(tones)
            max_val = np.max(norm)
            is_perm = (len(set(norm)) == len(norm) and
                       set(norm) == set(range(len(norm))))
            if is_perm:
                # Check Costas property
                diffs = set()
                is_costas = True
                for a in range(length):
                    for b in range(a+1, length):
                        d = (int(norm[b]) - int(norm[a]), b - a)
                        if d in diffs:
                            is_costas = False
                            break
                        diffs.add(d)
                    if not is_costas:
                        break
                print(f"  Block [{start}:{end}] ({length} symbols): {norm.tolist()}")
                print(f"    Costas: {is_costas}")
            else:
                print(f"  Block [{start}:{end}] ({length} symbols): {norm.tolist()} - not a permutation")

    # ================================================================
    # VISUALIZATION
    # ================================================================
    fig, axes = plt.subplots(3, 1, figsize=(22, 12))

    x = np.arange(NSYMBOLS)

    # 1. Tone comparison
    w = 0.35
    c1_colors = []
    c2_colors = []
    for i in range(NSYMBOLS):
        if classification[i] == 'S':
            c1_colors.append('green')
            c2_colors.append('green')
        elif classification[i] == '?':
            c1_colors.append('orange')
            c2_colors.append('orange')
        else:
            c1_colors.append('#cc4444')
            c2_colors.append('#4444cc')

    axes[0].bar(x - w/2, np.clip(t1, 0, 7), w, color=c1_colors, alpha=0.8,
                label='Capture 1', edgecolor='black', linewidth=0.3)
    axes[0].bar(x + w/2, np.clip(t2, 0, 7), w, color=c2_colors, alpha=0.8,
                label='Capture 2', edgecolor='black', linewidth=0.3)
    axes[0].set_ylabel("Tone (0-7)")
    axes[0].set_title("Cross-Reference: Green=SYNC, Red/Blue=DATA, Orange=UNCERTAIN")
    axes[0].legend(loc='upper right')
    axes[0].set_ylim(-0.5, 8.5)

    # 2. Frame structure
    fc = []
    for i in range(NSYMBOLS):
        if classification[i] == 'S':
            fc.append('green')
        elif classification[i] == '?':
            fc.append('orange')
        else:
            fc.append('red')
    axes[1].bar(x, [1]*NSYMBOLS, color=fc, edgecolor='black', linewidth=0.5)
    for i in range(NSYMBOLS):
        if classification[i] == 'S':
            axes[1].text(i, 0.5, str(t1[i]), ha='center', va='center',
                        fontsize=6, fontweight='bold', color='white')
        elif classification[i] == '?':
            axes[1].text(i, 0.5, '?', ha='center', va='center',
                        fontsize=6, fontweight='bold', color='black')
    axes[1].set_title("Frame: Green=SYNC (tone shown), Red=DATA, Orange=UNCERTAIN")
    axes[1].set_yticks([])

    # 3. Raw frequencies
    axes[2].plot(x, r1['freqs'], 'r.-', markersize=4, linewidth=0.8, label='Capture 1')
    axes[2].plot(x, r2['freqs'], 'b.-', markersize=4, linewidth=0.8, label='Capture 2')
    for pos in uncertain_pos:
        axes[2].axvline(pos, color='orange', alpha=0.3, linewidth=2)
    axes[2].set_ylabel("Frequency (Hz)")
    axes[2].set_xlabel("Symbol Position")
    axes[2].set_title("Peak Frequencies (orange bands = uncertain positions)")
    axes[2].legend()
    axes[2].grid(True, alpha=0.2)

    for ax in axes:
        ax.set_xlim(-0.5, NSYMBOLS - 0.5)
        ax.set_xticks(range(0, NSYMBOLS, 2))

    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / "xref_v3.png", dpi=150)
    plt.close()
    print(f"\nSaved xref_v3.png")


if __name__ == "__main__":
    main()
