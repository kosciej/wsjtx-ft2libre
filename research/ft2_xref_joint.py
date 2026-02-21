#!/usr/bin/env python3
"""
Cross-reference v4: Joint alignment optimization.
Instead of optimizing each capture's alignment independently (which may give
different symbol boundaries), we scan all alignment combinations and find
the one that maximizes matching symbol count.
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
        if e > len(signal):
            freqs[i] = np.nan
            continue
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


def quantize(freqs, tone_spacing=TONE_SPACING, max_outliers=3):
    """Quantize frequencies to 0-7 tones with best base frequency."""
    valid_mask = ~np.isnan(freqs)
    f = freqs[valid_mask]
    if len(f) == 0:
        return None, 999, 0

    best_rms = 999
    best_base = 0
    best_tones_valid = None

    for base in np.arange(np.min(f) - 2*tone_spacing,
                          np.min(f) + tone_spacing, 0.3):
        tf = (f - base) / tone_spacing
        ti = np.round(tf).astype(int)
        in_range = (ti >= 0) & (ti <= 7)
        n_outliers = np.sum(~in_range)

        if n_outliers > max_outliers:
            continue

        # RMS only on in-range symbols
        if np.sum(in_range) > 0:
            rms = np.sqrt(np.mean((tf[in_range] - ti[in_range])**2))
            if rms < best_rms:
                best_rms = rms
                best_base = base
                best_tones_valid = ti.copy()

    if best_tones_valid is None:
        return None, 999, 0

    # Map back to full array
    tones = np.full(len(freqs), -99, dtype=int)
    tones[valid_mask] = best_tones_valid
    return tones, best_rms, best_base


def find_burst_start(signal, skip_first=False):
    """Find the start sample of the target burst."""
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
        print(f"    Burst {i}: {s/FS:.3f}-{e/FS:.3f}s ({(e-s)/FS:.3f}s)")

    if skip_first and len(bursts) > 1:
        return bursts[1][0]
    elif bursts:
        return max(bursts, key=lambda b: b[1]-b[0])[0]
    return None


def main():
    sig1, _ = sf.read("ft2_capture.wav")
    if len(sig1.shape) > 1: sig1 = sig1[:, 0]
    sig2, _ = sf.read("ft2_capture2.wav")
    if len(sig2.shape) > 1: sig2 = sig2[:, 0]

    print("Capture 1 bursts:")
    burst1 = find_burst_start(sig1, skip_first=False)
    print("Capture 2 bursts:")
    burst2 = find_burst_start(sig2, skip_first=True)

    margin = int(0.05 * FS)
    center1 = burst1 + margin
    center2 = burst2 + margin

    # ================================================================
    # PHASE 1: Pre-compute tone sequences for many alignments
    # ================================================================
    print(f"\nPre-computing tone sequences...")
    step = NSPS // 36  # 10-sample steps

    offsets = list(range(-NSPS, NSPS + 1, step))

    # Cache: offset -> (tones, rms, base, freqs)
    cache1 = {}
    cache2 = {}

    for offset in offsets:
        start1 = center1 + offset
        if 0 <= start1 and start1 + NSYMBOLS * NSPS <= len(sig1):
            freqs1 = extract_freqs(sig1, start1)
            tones1, rms1, base1 = quantize(freqs1)
            if tones1 is not None:
                cache1[offset] = (tones1, rms1, base1, freqs1)

        start2 = center2 + offset
        if 0 <= start2 and start2 + NSYMBOLS * NSPS <= len(sig2):
            freqs2 = extract_freqs(sig2, start2)
            tones2, rms2, base2 = quantize(freqs2)
            if tones2 is not None:
                cache2[offset] = (tones2, rms2, base2, freqs2)

    print(f"  Capture 1: {len(cache1)} valid alignments")
    print(f"  Capture 2: {len(cache2)} valid alignments")

    # ================================================================
    # PHASE 2: Find alignment pair that maximizes matches
    # ================================================================
    print(f"\nSearching for best alignment pair...")

    best_matches = 0
    best_pair = None
    best_details = None

    for off1 in cache1:
        t1, rms1, base1, f1 = cache1[off1]
        for off2 in cache2:
            t2, rms2, base2, f2 = cache2[off2]

            # Count matches (only where both are in 0-7 range)
            both_valid = (t1 >= 0) & (t1 <= 7) & (t2 >= 0) & (t2 <= 7)
            matches = np.sum((t1 == t2) & both_valid)

            if matches > best_matches:
                best_matches = matches
                best_pair = (off1, off2)
                best_details = {
                    't1': t1.copy(), 't2': t2.copy(),
                    'f1': f1.copy(), 'f2': f2.copy(),
                    'rms1': rms1, 'rms2': rms2,
                    'base1': base1, 'base2': base2,
                }

    # Also check top candidates
    results = []
    for off1 in cache1:
        t1, rms1, _, _ = cache1[off1]
        for off2 in cache2:
            t2, rms2, _, _ = cache2[off2]
            both_valid = (t1 >= 0) & (t1 <= 7) & (t2 >= 0) & (t2 <= 7)
            matches = np.sum((t1 == t2) & both_valid)
            results.append((matches, rms1, rms2, off1, off2))

    results.sort(key=lambda r: (-r[0], r[1]+r[2]))

    print(f"\nTop 10 alignment pairs by match count:")
    for matches, rms1, rms2, off1, off2 in results[:10]:
        print(f"  matches={matches:2d}  off1={off1:+4d} off2={off2:+4d}  "
              f"rms1={rms1:.4f} rms2={rms2:.4f}")

    if best_pair is None:
        print("No valid alignment pair found!")
        return

    off1, off2 = best_pair
    d = best_details
    t1, t2 = d['t1'], d['t2']

    print(f"\n{'='*70}")
    print(f"BEST ALIGNMENT: off1={off1:+d} off2={off2:+d}, matches={best_matches}")
    print(f"{'='*70}")
    print(f"  C1: RMS={d['rms1']:.4f}, Base={d['base1']:.1f} Hz")
    print(f"  C2: RMS={d['rms2']:.4f}, Base={d['base2']:.1f} Hz")
    print(f"  C1 tones: {t1.tolist()}")
    print(f"  C2 tones: {t2.tolist()}")

    # Classification
    both_valid = (t1 >= 0) & (t1 <= 7) & (t2 >= 0) & (t2 <= 7)
    match_mask = (t1 == t2) & both_valid

    classification = []
    for i in range(NSYMBOLS):
        if not both_valid[i]:
            classification.append('?')
        elif match_mask[i]:
            classification.append('S')
        else:
            classification.append('D')

    n_sync = classification.count('S')
    n_data = classification.count('D')
    n_unc = classification.count('?')

    print(f"\nSYNC: {n_sync}, DATA: {n_data}, UNCERTAIN: {n_unc}")
    print(f"Data bits: {n_data} × 3 = {n_data * 3}")

    # Frame map
    print(f"\nFrame map:")
    for row in range(0, NSYMBOLS, 36):
        end = min(row + 36, NSYMBOLS)
        print(f"  Pos:  {' '.join(f'{i:2d}' for i in range(row, end))}")
        print(f"  Map:  {' '.join(f' {classification[i]}' for i in range(row, end))}")
        print(f"  C1:   {' '.join(f'{t1[i]:2d}' for i in range(row, end))}")
        print(f"  C2:   {' '.join(f'{t2[i]:2d}' for i in range(row, end))}")
        print()

    sync_pos = [i for i in range(NSYMBOLS) if classification[i] == 'S']
    data_pos = [i for i in range(NSYMBOLS) if classification[i] == 'D']
    unc_pos = [i for i in range(NSYMBOLS) if classification[i] == '?']

    print(f"SYNC positions: {sync_pos}")
    print(f"Sync values: {[int(t1[i]) for i in sync_pos]}")
    print(f"DATA positions: {data_pos}")
    print(f"UNCERTAIN positions: {unc_pos}")

    # Contiguous sync blocks
    print(f"\nContiguous sync blocks:")
    blocks = []
    in_block = False
    bs = 0
    for i in range(NSYMBOLS + 1):
        if i < NSYMBOLS and classification[i] == 'S':
            if not in_block:
                bs = i
                in_block = True
        else:
            if in_block:
                bl = i - bs
                bt = t1[bs:i]
                blocks.append((bs, i, bl, bt))
                norm = bt - np.min(bt)
                is_perm = len(set(norm)) == len(norm) and set(norm) == set(range(len(norm)))
                print(f"  [{bs}:{i}] len={bl} tones={bt.tolist()} norm={norm.tolist()}"
                      f"{' PERM' if is_perm else ''}")
                in_block = False

    # Check: what if we also include uncertain positions as potential sync?
    print(f"\nExtended sync (including uncertain as potential sync):")
    for pos in unc_pos:
        # Check if adjacent to sync
        adj = []
        if pos > 0 and classification[pos-1] == 'S':
            adj.append(pos-1)
        if pos < NSYMBOLS-1 and classification[pos+1] == 'S':
            adj.append(pos+1)
        t1v = int(t1[pos]) if 0 <= t1[pos] <= 7 else '?'
        t2v = int(t2[pos]) if 0 <= t2[pos] <= 7 else '?'
        print(f"  Pos {pos}: C1={t1v} C2={t2v} adj_sync={adj}")

    # ================================================================
    # Check known sync structures
    # ================================================================
    print(f"\n{'='*70}")
    print("PATTERN ANALYSIS")
    print(f"{'='*70}")

    # With n_sync sync positions, what patterns could they form?
    # Expected: 14 sync symbols for 58 data (174 bits)
    # or 20 sync for 52 data (156 bits) - not enough
    # 14 sync = 2 × 7 Costas or some other arrangement

    # Check if the sync positions happen to be in two groups of 7
    if n_sync >= 12:
        print(f"\nLooking for two 7-symbol sync blocks among the {n_sync} sync positions:")
        # Try all possible positions for two 7-symbol blocks
        best_score = 0
        best_config = None
        for s1 in range(NSYMBOLS - 6):
            for s2 in range(s1 + 7, NSYMBOLS - 6):
                positions = list(range(s1, s1+7)) + list(range(s2, s2+7))
                sync_correct = sum(1 for p in positions if classification[p] == 'S')
                # Also check: data positions should NOT be in this set
                data_wrong = sum(1 for p in positions if classification[p] == 'D')
                score = sync_correct - data_wrong
                if score > best_score:
                    best_score = score
                    n_data_in = NSYMBOLS - 14
                    best_config = (s1, s2, sync_correct, data_wrong, n_data_in)

        if best_config:
            s1, s2, sc, dw, nd = best_config
            print(f"  Best 2×7: [{s1}:{s1+7}] + [{s2}:{s2+7}]")
            print(f"    Sync correct: {sc}/14, Data wrongly included: {dw}")
            print(f"    Data symbols: {nd} × 3 = {nd*3} bits")
            # Show the sync tones at these positions
            block1 = t1[s1:s1+7]
            block2 = t1[s2:s2+7]
            n1 = block1 - np.min(block1)
            n2 = block2 - np.min(block2)
            print(f"    Block 1 tones: {block1.tolist()} norm: {n1.tolist()}")
            print(f"    Block 2 tones: {block2.tolist()} norm: {n2.tolist()}")

        # Also try 3×7
        print(f"\n  Looking for three 7-symbol sync blocks...")
        best_score = 0
        best_config = None
        for s1 in range(NSYMBOLS - 20):
            for s2 in range(s1 + 7, NSYMBOLS - 13):
                for s3 in range(s2 + 7, NSYMBOLS - 6):
                    positions = (list(range(s1, s1+7)) +
                                list(range(s2, s2+7)) +
                                list(range(s3, s3+7)))
                    sc = sum(1 for p in positions if classification[p] == 'S')
                    dw = sum(1 for p in positions if classification[p] == 'D')
                    score = sc - dw
                    if score > best_score:
                        best_score = score
                        best_config = (s1, s2, s3, sc, dw)

        if best_config:
            s1, s2, s3, sc, dw = best_config
            nd = NSYMBOLS - 21
            print(f"  Best 3×7: [{s1}:{s1+7}] + [{s2}:{s2+7}] + [{s3}:{s3+7}]")
            print(f"    Sync correct: {sc}/21, Data wrongly included: {dw}")
            print(f"    Data symbols: {nd} × 3 = {nd*3} bits")
            for si, label in [(s1, "B1"), (s2, "B2"), (s3, "B3")]:
                block = t1[si:si+7]
                norm = block - np.min(block)
                match_status = [classification[p] for p in range(si, si+7)]
                print(f"    {label} [{si}:{si+7}]: tones={block.tolist()} "
                      f"norm={norm.tolist()} match={''.join(match_status)}")

    # ================================================================
    # VISUALIZATION
    # ================================================================
    fig, axes = plt.subplots(3, 1, figsize=(22, 12))
    x = np.arange(NSYMBOLS)

    w = 0.35
    c1c = ['green' if classification[i]=='S' else ('orange' if classification[i]=='?' else '#cc4444')
           for i in range(NSYMBOLS)]
    c2c = ['green' if classification[i]=='S' else ('orange' if classification[i]=='?' else '#4444cc')
           for i in range(NSYMBOLS)]

    axes[0].bar(x-w/2, np.clip(t1, -1, 8), w, color=c1c, alpha=0.8, label='Capture 1')
    axes[0].bar(x+w/2, np.clip(t2, -1, 8), w, color=c2c, alpha=0.8, label='Capture 2')
    axes[0].set_ylabel("Tone")
    axes[0].set_title(f"Joint-Aligned Cross-Reference (matches={best_matches})\n"
                      f"Green=SYNC({n_sync}), Red/Blue=DATA({n_data}), Orange=UNCERTAIN({n_unc})")
    axes[0].legend()
    axes[0].set_ylim(-1.5, 8.5)

    fc = ['green' if classification[i]=='S' else ('orange' if classification[i]=='?' else 'red')
          for i in range(NSYMBOLS)]
    axes[1].bar(x, [1]*NSYMBOLS, color=fc, edgecolor='black', linewidth=0.5)
    for i in range(NSYMBOLS):
        if classification[i] == 'S':
            axes[1].text(i, 0.5, str(int(t1[i])), ha='center', va='center',
                        fontsize=6, fontweight='bold', color='white')
    axes[1].set_title("Frame: Green=SYNC, Red=DATA, Orange=?")
    axes[1].set_yticks([])

    axes[2].plot(x, d['f1'], 'r.-', ms=4, lw=0.8, label='Capture 1')
    axes[2].plot(x, d['f2'], 'b.-', ms=4, lw=0.8, label='Capture 2')
    axes[2].set_ylabel("Frequency (Hz)")
    axes[2].set_xlabel("Symbol Position")
    axes[2].legend()
    axes[2].grid(True, alpha=0.2)

    for ax in axes:
        ax.set_xlim(-0.5, NSYMBOLS-0.5)
        ax.set_xticks(range(0, NSYMBOLS, 2))

    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / "xref_joint.png", dpi=150)
    plt.close()
    print(f"\nSaved xref_joint.png")


if __name__ == "__main__":
    main()
