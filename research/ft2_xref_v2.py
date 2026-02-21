#!/usr/bin/env python3
"""
Cross-reference v2: More robust tone extraction with diagnostic output.
Analyzes middle burst from both captures, cross-references to find sync.
"""

import numpy as np
import soundfile as sf
import matplotlib.pyplot as plt
from pathlib import Path

OUTPUT_DIR = Path("output_xref")
OUTPUT_DIR.mkdir(exist_ok=True)

FS = 12000
NSPS = 360
BAUD = FS / NSPS  # 33.333 Hz
H = 0.5
TONE_SPACING = H * BAUD  # 16.667 Hz
NSYMBOLS = 72
NFFT = NSPS * 4  # 1440


def get_power_envelope(signal, window_ms=50):
    window_len = int(window_ms / 1000 * FS)
    return np.convolve(signal**2, np.ones(window_len)/window_len, mode='same')


def extract_peak_freqs(signal, start, nsymbols=NSYMBOLS, nsps=NSPS, nfft=NFFT):
    """Extract peak frequency per symbol using parabolic interpolation."""
    freqs = np.zeros(nsymbols)
    powers = np.zeros(nsymbols)
    freq_bins = np.fft.rfftfreq(nfft, 1/FS)

    for i in range(nsymbols):
        s = start + i * nsps
        e = s + nsps
        if e > len(signal):
            freqs[i] = np.nan
            continue

        seg = signal[s:e] * np.hanning(nsps)
        spec = np.abs(np.fft.rfft(seg, n=nfft))

        # Search in signal band
        lo_bin = int(200 * nfft / FS)
        hi_bin = int(2500 * nfft / FS)
        band = spec[lo_bin:hi_bin]
        pk = np.argmax(band) + lo_bin

        # Parabolic interpolation
        if 0 < pk < len(spec) - 1:
            a = np.log(spec[pk-1] + 1e-30)
            b = np.log(spec[pk] + 1e-30)
            c = np.log(spec[pk+1] + 1e-30)
            d = a - 2*b + c
            if abs(d) > 1e-10:
                delta = 0.5 * (a - c) / d
            else:
                delta = 0
            freqs[i] = (pk + delta) * FS / nfft
        else:
            freqs[i] = pk * FS / nfft

        powers[i] = spec[pk]

    return freqs, powers


def quantize_best(freqs, tone_spacing=TONE_SPACING, max_tones=8):
    """Find best quantization to 0..max_tones-1."""
    valid = ~np.isnan(freqs)
    f = freqs[valid]

    best_score = -1
    best_tones = None
    best_base = 0
    best_rms = 999

    for base in np.arange(np.min(f) - 2*tone_spacing, np.min(f) + tone_spacing, 0.2):
        tf = (f - base) / tone_spacing
        ti = np.round(tf).astype(int)
        res = tf - ti

        if np.any(ti < 0) or np.any(ti >= max_tones):
            continue

        rms = np.sqrt(np.mean(res**2))
        unique = len(np.unique(ti))

        # Score: prefer low RMS and high unique count
        score = unique / (rms + 0.01)

        if score > best_score:
            best_score = score
            best_rms = rms
            best_base = base
            result = np.full(len(freqs), -1, dtype=int)
            result[valid] = ti
            best_tones = result.copy()

    return best_tones, best_rms, best_base


def analyze_capture(filename, label, skip_first=False):
    """Analyze a capture: detect bursts, extract tones from middle burst."""
    print(f"\n{'='*60}")
    print(f"Analyzing {label}: {filename}")
    print(f"{'='*60}")

    signal, sr = sf.read(filename)
    if len(signal.shape) > 1:
        signal = signal[:, 0]
    print(f"Duration: {len(signal)/FS:.2f}s, max={np.max(np.abs(signal)):.4f}")

    # Detect bursts
    power = get_power_envelope(signal)
    threshold = np.max(power) * 0.01
    above = power > threshold
    transitions = np.diff(above.astype(int))
    starts = np.where(transitions == 1)[0]
    ends = np.where(transitions == -1)[0]
    if above[0]:
        starts = np.concatenate([[0], starts])
    if above[-1]:
        ends = np.concatenate([ends, [len(signal)-1]])

    n = min(len(starts), len(ends))
    bursts = [(starts[i], ends[i]) for i in range(n)
              if (ends[i] - starts[i]) / FS > 1.5]

    print(f"Bursts found: {len(bursts)}")
    for i, (s, e) in enumerate(bursts):
        print(f"  Burst {i}: {s/FS:.3f} - {e/FS:.3f}s ({(e-s)/FS:.3f}s)")

    # Select burst: skip first if requested, use longest otherwise
    if skip_first and len(bursts) > 1:
        target = bursts[1]
        print(f"  -> Using burst 1 (skipping first)")
    elif len(bursts) > 0:
        # Use the longest burst
        longest = max(bursts, key=lambda b: b[1]-b[0])
        idx = bursts.index(longest)
        print(f"  -> Using burst {idx} (longest)")
        target = longest
    else:
        print("  ERROR: No bursts found!")
        return None

    burst_start, burst_end = target
    burst_dur = (burst_end - burst_start) / FS
    print(f"  Burst: {burst_start/FS:.3f} - {burst_end/FS:.3f}s ({burst_dur:.3f}s)")

    # Fine alignment scan: try many offsets around the burst
    print(f"\nFine alignment scan...")
    scan_results = []
    scan_center = burst_start + int(0.05 * FS)  # 50ms after burst start

    for offset in range(-NSPS, NSPS + 1, NSPS // 36):
        start = scan_center + offset
        if start < 0 or start + NSYMBOLS * NSPS > len(signal):
            continue

        freqs, powers = extract_peak_freqs(signal, start)
        tones, rms, base = quantize_best(freqs)

        if tones is not None:
            valid_mask = tones >= 0
            unique = len(np.unique(tones[valid_mask]))
            scan_results.append({
                'offset': offset,
                'start': start,
                'tones': tones,
                'freqs': freqs,
                'rms': rms,
                'base': base,
                'unique': unique,
            })

    if not scan_results:
        print("  FAILED: no valid alignment found")
        return None

    # Pick best by lowest RMS
    scan_results.sort(key=lambda r: r['rms'])
    best = scan_results[0]

    print(f"  Best alignment:")
    print(f"    Offset: {best['offset']} samples ({best['offset']/FS*1000:.1f}ms)")
    print(f"    Start sample: {best['start']} ({best['start']/FS:.4f}s)")
    print(f"    RMS: {best['rms']:.4f}")
    print(f"    Base freq: {best['base']:.1f} Hz")
    print(f"    Unique tones: {best['unique']}")
    print(f"    Freq range: {np.nanmin(best['freqs']):.1f} - {np.nanmax(best['freqs']):.1f} Hz")
    print(f"    Tones: {best['tones'].tolist()}")

    # Diagnostic: show top 5 alignments
    print(f"\n  Top 5 alignments by RMS:")
    for i, r in enumerate(scan_results[:5]):
        print(f"    #{i}: offset={r['offset']:+5d} RMS={r['rms']:.4f} "
              f"base={r['base']:.1f} unique={r['unique']}")

    return best


def main():
    # Analyze both captures
    r1 = analyze_capture("ft2_capture.wav", "Capture 1 (unknown msg)", skip_first=False)
    r2 = analyze_capture("ft2_capture2.wav", "Capture 2 (CQ CQ SP6TLS KO02)", skip_first=True)

    if r1 is None or r2 is None:
        print("\nCannot cross-reference: tone extraction failed for one or both captures.")
        return

    # ================================================================
    # CROSS-REFERENCE
    # ================================================================
    t1 = r1['tones']
    t2 = r2['tones']

    print(f"\n{'='*80}")
    print(f"CROSS-REFERENCE")
    print(f"{'='*80}")

    print(f"\nCapture 1: {t1.tolist()}")
    print(f"Capture 2: {t2.tolist()}")

    matches = (t1 == t2)
    n_match = np.sum(matches)
    print(f"\nMatching positions: {n_match}/{NSYMBOLS} ({n_match/NSYMBOLS*100:.1f}%)")

    sync_pos = np.where(matches)[0]
    data_pos = np.where(~matches)[0]

    print(f"\nSYNC positions ({len(sync_pos)}): {sync_pos.tolist()}")
    print(f"Sync tone values: {t1[sync_pos].tolist()}")
    print(f"\nDATA positions ({len(data_pos)}): {data_pos.tolist()}")
    print(f"  Cap1 data tones: {t1[data_pos].tolist()}")
    print(f"  Cap2 data tones: {t2[data_pos].tolist()}")

    # Frame map
    print(f"\nFrame structure (S=sync, D=data):")
    for row in range(0, NSYMBOLS, 36):
        end = min(row + 36, NSYMBOLS)
        pos_str = ' '.join(f'{i:2d}' for i in range(row, end))
        map_str = ' '.join(f' {"S" if matches[i] else "D"}' for i in range(row, end))
        t1_str = ' '.join(f'{t1[i]:2d}' for i in range(row, end))
        t2_str = ' '.join(f'{t2[i]:2d}' for i in range(row, end))
        print(f"  Pos:  {pos_str}")
        print(f"  Map:  {map_str}")
        print(f"  C1:   {t1_str}")
        print(f"  C2:   {t2_str}")
        print()

    # Look for contiguous sync blocks
    print("Contiguous sync blocks:")
    in_block = False
    block_start = 0
    blocks = []
    for i in range(NSYMBOLS + 1):
        if i < NSYMBOLS and matches[i]:
            if not in_block:
                block_start = i
                in_block = True
        else:
            if in_block:
                block_len = i - block_start
                block_tones = t1[block_start:i]
                blocks.append((block_start, i, block_tones))
                norm = block_tones - np.min(block_tones)
                is_perm = (len(set(norm)) == len(norm) and
                           set(norm) == set(range(len(norm))))
                print(f"  [{block_start}:{i}] len={block_len} "
                      f"tones={block_tones.tolist()} norm={norm.tolist()}"
                      f"{' <- PERMUTATION!' if is_perm else ''}")
                in_block = False

    # Check for 7-symbol Costas-like blocks
    print(f"\n7-symbol blocks that are permutations of 0-6:")
    for start, end, tones in blocks:
        if end - start == 7:
            norm = tones - np.min(tones)
            if set(norm) == set(range(7)):
                # Check Costas property
                diffs = set()
                is_costas = True
                for a in range(7):
                    for b in range(a+1, 7):
                        d = (norm[b] - norm[a], b - a)
                        if d in diffs:
                            is_costas = False
                            break
                        diffs.add(d)
                    if not is_costas:
                        break
                print(f"  [{start}:{end}] {norm.tolist()} "
                      f"{'COSTAS!' if is_costas else 'not Costas'}")

    # Data capacity
    n_sync = len(sync_pos)
    n_data = len(data_pos)
    n_bits = n_data * 3
    print(f"\nCapacity: {n_sync} sync + {n_data} data = {n_bits} coded bits")
    print(f"  LDPC(174,91) needs 174 bits = 58 data symbols")

    # ================================================================
    # VISUALIZATION
    # ================================================================
    fig, axes = plt.subplots(4, 1, figsize=(22, 16))

    x = np.arange(NSYMBOLS)

    # 1. Side-by-side tone bars
    w = 0.35
    c1 = ['green' if matches[i] else '#cc4444' for i in range(NSYMBOLS)]
    c2 = ['green' if matches[i] else '#4444cc' for i in range(NSYMBOLS)]
    axes[0].bar(x - w/2, t1, w, color=c1, alpha=0.8, label='Capture 1')
    axes[0].bar(x + w/2, t2, w, color=c2, alpha=0.8, label='Capture 2')
    axes[0].set_ylabel("Tone (0-7)")
    axes[0].set_title("Tone Cross-Reference (Green=SYNC, Red/Blue=DATA)")
    axes[0].legend(loc='upper right')
    axes[0].set_ylim(-0.5, 8.5)

    # 2. Frame structure
    fc = ['green' if matches[i] else 'red' for i in range(NSYMBOLS)]
    axes[1].bar(x, [1]*NSYMBOLS, color=fc, edgecolor='black', linewidth=0.5)
    for i in range(NSYMBOLS):
        if matches[i]:
            axes[1].text(i, 0.5, str(t1[i]), ha='center', va='center',
                        fontsize=6, fontweight='bold', color='white')
    axes[1].set_title("Frame Map: Green=SYNC (tone value shown), Red=DATA")
    axes[1].set_yticks([])

    # 3. Raw frequencies
    f1 = r1['freqs']
    f2 = r2['freqs']
    axes[2].plot(x, f1, 'r.-', markersize=4, linewidth=0.8, label='Capture 1')
    axes[2].plot(x, f2, 'b.-', markersize=4, linewidth=0.8, label='Capture 2')
    # Grid lines at tone frequencies
    for cap_label, result, color in [("C1", r1, 'red'), ("C2", r2, 'blue')]:
        for tone in range(8):
            f_tone = result['base'] + tone * TONE_SPACING
            axes[2].axhline(f_tone, color=color, alpha=0.15, linewidth=0.5)
    axes[2].set_ylabel("Frequency (Hz)")
    axes[2].set_title("Peak Frequencies per Symbol")
    axes[2].legend()
    axes[2].grid(True, alpha=0.2)

    # 4. Difference (t1 - t2) at each position
    diff = t1.astype(float) - t2.astype(float)
    dc = ['green' if matches[i] else 'red' for i in range(NSYMBOLS)]
    axes[3].bar(x, diff, color=dc)
    axes[3].axhline(0, color='black', linewidth=0.5)
    axes[3].set_ylabel("Tone difference (C1 - C2)")
    axes[3].set_xlabel("Symbol Position")
    axes[3].set_title("Tone Difference (0 = matching = likely sync)")

    for ax in axes:
        ax.set_xlim(-0.5, NSYMBOLS - 0.5)
        ax.set_xticks(range(0, NSYMBOLS, 2))

    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / "xref_analysis.png", dpi=150)
    plt.close()
    print(f"\nSaved xref_analysis.png")

    # ================================================================
    # KNOWN-PATTERN CHECK
    # ================================================================
    print(f"\n{'='*80}")
    print(f"KNOWN PATTERN CHECKS")
    print(f"{'='*80}")

    # What if sync is 2 blocks of 7 at start and end?
    for sync_len in [4, 7]:
        for n_blocks in [2, 3]:
            if n_blocks == 2:
                positions_list = []
                for gap in range(NSYMBOLS - 2*sync_len + 1):
                    p1 = list(range(0, sync_len))
                    p2 = list(range(NSYMBOLS - sync_len, NSYMBOLS))
                    positions_list.append((f"[0:{sync_len}]+[{NSYMBOLS-sync_len}:{NSYMBOLS}]",
                                          p1 + p2))
                    # Also try middle positions
                    mid = NSYMBOLS // 2
                    p_mid = list(range(mid - sync_len//2, mid - sync_len//2 + sync_len))
                    positions_list.append((f"[0:{sync_len}]+mid", p1 + p_mid))
            elif n_blocks == 3:
                positions_list = []
                for mid_start in range(sync_len, NSYMBOLS - 2*sync_len):
                    p1 = list(range(0, sync_len))
                    p2 = list(range(mid_start, mid_start + sync_len))
                    p3 = list(range(NSYMBOLS - sync_len, NSYMBOLS))
                    name = f"[0:{sync_len}]+[{mid_start}:{mid_start+sync_len}]+[{NSYMBOLS-sync_len}:{NSYMBOLS}]"
                    positions_list.append((name, p1 + p2 + p3))

            best_pattern = None
            best_score = 0
            for name, positions in positions_list:
                sync_correct = sum(1 for p in positions if matches[p])
                data_correct = sum(1 for p in range(NSYMBOLS)
                                   if p not in positions and not matches[p])
                total_sync = len(positions)
                total_data = NSYMBOLS - total_sync
                score = sync_correct + data_correct
                if score > best_score:
                    best_score = score
                    best_pattern = (name, sync_correct, total_sync, data_correct, total_data)

            if best_pattern:
                name, sc, ts, dc_val, td = best_pattern
                print(f"\n  Best {n_blocks}×{sync_len} pattern: {name}")
                print(f"    Sync: {sc}/{ts}, Data: {dc_val}/{td}, "
                      f"Total score: {best_score}/{NSYMBOLS}")


if __name__ == "__main__":
    main()
