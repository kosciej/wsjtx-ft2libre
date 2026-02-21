#!/usr/bin/env python3
"""Quick examination of ft2_capture2.wav to understand signal characteristics."""

import numpy as np
import soundfile as sf
import matplotlib.pyplot as plt
from pathlib import Path

OUTPUT_DIR = Path("output_xref")
OUTPUT_DIR.mkdir(exist_ok=True)

signal, fs = sf.read("ft2_capture2.wav")
if len(signal.shape) > 1:
    signal = signal[:, 0]

print(f"Duration: {len(signal)/fs:.2f}s")
print(f"Sample rate: {fs} Hz")
print(f"Samples: {len(signal)}")
print(f"Max amplitude: {np.max(np.abs(signal)):.6f}")
print(f"RMS: {np.sqrt(np.mean(signal**2)):.6f}")

# Plot waveform
fig, axes = plt.subplots(3, 1, figsize=(20, 12))

t = np.arange(len(signal)) / fs
axes[0].plot(t, signal, linewidth=0.2)
axes[0].set_title(f"Waveform (duration={len(signal)/fs:.2f}s, fs={fs}Hz)")
axes[0].set_ylabel("Amplitude")
axes[0].set_xlabel("Time (s)")

# Power envelope - sliding window
window_ms = 50
window_len = int(window_ms / 1000 * fs)
power = np.convolve(signal**2, np.ones(window_len)/window_len, mode='same')
axes[1].plot(t, 10*np.log10(power + 1e-30), linewidth=0.5)
axes[1].set_title("Power envelope (50ms window)")
axes[1].set_ylabel("Power (dB)")
axes[1].set_xlabel("Time (s)")

# Show power percentiles
for pct in [10, 50, 90, 99]:
    val = np.percentile(power, pct)
    axes[1].axhline(10*np.log10(val + 1e-30), color='r', alpha=0.3, linestyle='--')
    axes[1].text(0, 10*np.log10(val + 1e-30), f'P{pct}', fontsize=8)

# Spectrogram
nfft = 2048
axes[2].specgram(signal, NFFT=nfft, Fs=fs, noverlap=nfft//2, cmap='viridis')
axes[2].set_title("Spectrogram")
axes[2].set_ylabel("Frequency (Hz)")
axes[2].set_xlabel("Time (s)")
axes[2].set_ylim(0, 3000)

plt.tight_layout()
plt.savefig(OUTPUT_DIR / "capture2_overview.png", dpi=150)
plt.close()

# Also zoom into first few seconds
fig, ax = plt.subplots(figsize=(20, 6))
ax.specgram(signal[:int(10*fs)], NFFT=1024, Fs=fs, noverlap=512, cmap='viridis')
ax.set_title("First 10 seconds - Spectrogram")
ax.set_ylabel("Frequency (Hz)")
ax.set_xlabel("Time (s)")
ax.set_ylim(0, 3000)
plt.tight_layout()
plt.savefig(OUTPUT_DIR / "capture2_first10s.png", dpi=150)
plt.close()

print("\nSaved overview and spectrogram plots.")
