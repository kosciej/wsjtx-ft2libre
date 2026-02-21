#!/usr/bin/env python3
"""
Cross-reference with 4-GFSK model.
4 tones at ~41.3 Hz spacing, 2 bits per symbol.
Scan NSPS range to find best cross-reference match.
"""

import numpy as np
import soundfile as sf
import matplotlib.pyplot as plt
from pathlib import Path

OUTPUT_DIR = Path("output_raw")
OUTPUT_DIR.mkdir(exist_ok=True)

FS = 12000
TONE_FREQS = np.array([1500.6, 1541.9, 1583.2, 1624.5])


def load_burst(fn, skip_first=False):
    sig, _ = sf.read(fn)
    if len(sig.shape) > 1: sig = sig[:, 0]
    wl = int(0.05 * FS)
    pw = np.convolve(sig**2, np.ones(wl)/wl, mode='same')
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
        s, e = bursts[1]
    else:
        s, e = max(bursts, key=lambda b: b[1]-b[0])
    return sig, s, e


def extract_tones_4(signal, start, nsymbols, nsps):
    """Extract 4-tone sequence."""
    nfft = max(nsps * 4, 1024)
    tones = np.zeros(nsymbols, dtype=int)
    freqs = np.zeros(nsymbols)

    for i in range(nsymbols):
        s = start + i * nsps
        e = s + nsps
        if e > len(signal):
            tones[i] = -1
            freqs[i] = np.nan
            continue
        seg = signal[s:e] * np.hanning(nsps)
        spec = np.abs(np.fft.rfft(seg, n=nfft))
        lo = int(1400 * nfft / FS)
        hi = int(1700 * nfft / FS)
        pk = lo + np.argmax(spec[lo:hi])

        if 0 < pk < len(spec) - 1:
            a = np.log(spec[pk-1]+1e-30)
            b = np.log(spec[pk]+1e-30)
            c = np.log(spec[pk+1]+1e-30)
            d = a - 2*b + c
            delta = 0.5*(a-c)/d if abs(d) > 1e-10 else 0
            freq = (pk+delta)*FS/nfft
        else:
            freq = pk*FS/nfft

        freqs[i] = freq
        tones[i] = np.argmin(np.abs(TONE_FREQS - freq))

    return tones, freqs


def main():
    sig1, _ = sf.read("ft2_capture.wav")
    if len(sig1.shape) > 1: sig1 = sig1[:, 0]
    sig2, _ = sf.read("ft2_capture2.wav")
    if len(sig2.shape) > 1: sig2 = sig2[:, 0]

    # Detect bursts
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
            return bursts[1]
        return max(bursts, key=lambda b: b[1]-b[0]) if bursts else None

    b1s, b1e = detect(sig1)
    b2s, b2e = detect(sig2, skip_first=True)
    print(f"Burst 1: {b1s/FS:.3f}-{b1e/FS:.3f}s ({(b1e-b1s)/FS:.3f}s)")
    print(f"Burst 2: {b2s/FS:.3f}-{b2e/FS:.3f}s ({(b2e-b2s)/FS:.3f}s)")

    # Scan NSPS values and find which gives best cross-reference
    print(f"\nScanning NSPS with 4-tone cross-reference:")
    print(f"{'NSPS':>5} {'Baud':>7} {'h':>5} {'Nsym':>5} "
          f"{'MaxMatch':>9} {'BestOff1':>9} {'BestOff2':>9} "
          f"{'Data':>5} {'Bits':>5}")

    best_overall = None
    best_overall_matches = 0

    for nsps in range(260, 400, 2):
        baud = FS / nsps
        h = np.mean(np.diff(TONE_FREQS)) / baud
        nsymbols = int((b1e - b1s - 0.1*FS) / nsps)  # conservative symbol count
        if nsymbols < 40:
            continue

        margin1 = int(0.03 * FS)
        margin2 = int(0.03 * FS)

        best_matches = 0
        best_off1 = 0
        best_off2 = 0
        best_t1 = None
        best_t2 = None

        # Coarse scan
        for off1 in range(0, nsps, nsps//6):
            t1, f1 = extract_tones_4(sig1, b1s + margin1 + off1, nsymbols, nsps)

            for off2 in range(0, nsps, nsps//6):
                t2, f2 = extract_tones_4(sig2, b2s + margin2 + off2, nsymbols, nsps)

                valid = (t1 >= 0) & (t2 >= 0)
                matches = np.sum((t1 == t2) & valid)

                if matches > best_matches:
                    best_matches = matches
                    best_off1 = off1
                    best_off2 = off2
                    best_t1 = t1.copy()
                    best_t2 = t2.copy()

        # Fine scan around best
        for off1 in range(max(0, best_off1-nsps//6), best_off1+nsps//6, 2):
            for off2 in range(max(0, best_off2-nsps//6), best_off2+nsps//6, 2):
                t1, f1 = extract_tones_4(sig1, b1s + margin1 + off1, nsymbols, nsps)
                t2, f2 = extract_tones_4(sig2, b2s + margin2 + off2, nsymbols, nsps)
                valid = (t1 >= 0) & (t2 >= 0)
                matches = np.sum((t1 == t2) & valid)
                if matches > best_matches:
                    best_matches = matches
                    best_off1 = off1
                    best_off2 = off2
                    best_t1 = t1.copy()
                    best_t2 = t2.copy()

        n_data = nsymbols - best_matches
        bits = n_data * 2

        print(f"{nsps:5d} {baud:7.2f} {h:5.2f} {nsymbols:5d} "
              f"{best_matches:9d} {best_off1:9d} {best_off2:9d} "
              f"{n_data:5d} {bits:5d}")

        if best_matches > best_overall_matches:
            best_overall_matches = best_matches
            best_overall = {
                'nsps': nsps, 'baud': baud, 'h': h,
                'nsymbols': nsymbols,
                'matches': best_matches,
                'off1': best_off1, 'off2': best_off2,
                't1': best_t1, 't2': best_t2,
            }

    # ================================================================
    # Detailed analysis of best NSPS
    # ================================================================
    if best_overall is None:
        print("No good result found!")
        return

    r = best_overall
    print(f"\n{'='*70}")
    print(f"BEST: NSPS={r['nsps']}, baud={r['baud']:.2f}, h={r['h']:.3f}")
    print(f"  {r['nsymbols']} symbols, {r['matches']} matches")
    print(f"{'='*70}")

    t1 = r['t1']
    t2 = r['t2']
    nsymbols = r['nsymbols']

    matches = (t1 == t2)
    sync_pos = np.where(matches)[0]
    data_pos = np.where(~matches)[0]

    print(f"\nCapture 1 tones: {t1.tolist()}")
    print(f"Capture 2 tones: {t2.tolist()}")
    print(f"\nSYNC positions ({len(sync_pos)}): {sync_pos.tolist()}")
    print(f"Sync tones: {t1[sync_pos].tolist()}")
    print(f"DATA positions ({len(data_pos)}): {data_pos.tolist()}")
    print(f"Data capacity: {len(data_pos)} × 2 = {len(data_pos)*2} bits")

    # Frame map
    print(f"\nFrame map (S=sync, D=data):")
    cls = ['S' if matches[i] else 'D' for i in range(nsymbols)]
    for row in range(0, nsymbols, 40):
        end = min(row + 40, nsymbols)
        print(f"  Pos: {' '.join(f'{i:2d}' for i in range(row, end))}")
        print(f"  Map: {' '.join(f' {cls[i]}' for i in range(row, end))}")
        print(f"  C1:  {' '.join(f'{t1[i]:2d}' for i in range(row, end))}")
        print(f"  C2:  {' '.join(f'{t2[i]:2d}' for i in range(row, end))}")
        print()

    # Contiguous sync blocks
    print("Contiguous sync blocks (len >= 3):")
    in_b = False
    bs = 0
    for i in range(nsymbols + 1):
        if i < nsymbols and matches[i]:
            if not in_b:
                bs = i
                in_b = True
        else:
            if in_b:
                bl = i - bs
                if bl >= 3:
                    bt = t1[bs:i]
                    print(f"  [{bs}:{i}] len={bl} tones={bt.tolist()}")
                in_b = False

    # Visualization
    fig, axes = plt.subplots(2, 1, figsize=(22, 8))
    x = np.arange(nsymbols)

    w = 0.35
    c1c = ['green' if matches[i] else '#cc4444' for i in range(nsymbols)]
    c2c = ['green' if matches[i] else '#4444cc' for i in range(nsymbols)]
    axes[0].bar(x-w/2, t1, w, color=c1c, alpha=0.8, label='Capture 1')
    axes[0].bar(x+w/2, t2, w, color=c2c, alpha=0.8, label='Capture 2')
    axes[0].set_ylabel("Tone (0-3)")
    axes[0].set_title(f"4-GFSK Cross-Reference (NSPS={r['nsps']}, baud={r['baud']:.1f})\n"
                      f"Green=SYNC({len(sync_pos)}), Red/Blue=DATA({len(data_pos)})")
    axes[0].legend()
    axes[0].set_ylim(-0.5, 3.5)
    axes[0].set_yticks([0, 1, 2, 3])

    fc = ['green' if matches[i] else 'red' for i in range(nsymbols)]
    axes[1].bar(x, [1]*nsymbols, color=fc, edgecolor='black', linewidth=0.5)
    for i in range(nsymbols):
        if matches[i]:
            axes[1].text(i, 0.5, str(t1[i]), ha='center', va='center',
                        fontsize=5, fontweight='bold', color='white')
    axes[1].set_title("Frame: Green=SYNC, Red=DATA")
    axes[1].set_yticks([])
    axes[1].set_xlabel("Symbol Position")

    for ax in axes:
        ax.set_xlim(-0.5, nsymbols-0.5)
        ax.set_xticks(range(0, nsymbols, 5))

    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / "xref_4tone.png", dpi=150)
    plt.close()
    print(f"\nSaved xref_4tone.png")


if __name__ == "__main__":
    main()
