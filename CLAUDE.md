# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

WSJT-X 3.0.0-rc1 — a Qt5 application for weak-signal digital amateur radio communication. Supports modes: FT8, FT4, JT4, JT9, JT65, Q65, MSK144, FST4, WSPR, Echo, and SuperFox. Written in C++11 (Qt5 GUI/networking) and Fortran (DSP/decoding algorithms), with C utility code. Depends on a bundled fork of Hamlib for rig control.

## Repository Layout

- **Top level** — superbuild CMake wrapper (`CMakeLists.txt`) that builds Hamlib then WSJT-X
- **`wsjtx/`** — the full WSJT-X source tree (extracted from `src/wsjtx.tgz`)
- **`hamlib-4.7/`** — the Hamlib library source (extracted from `src/hamlib-4.7.tar.gz`)
- **`hamlib.patch` / `wsjtx.patch`** — patches applied before build (currently empty)

## Build Commands

Out-of-source builds are mandatory. Requires: CMake >= 3.7.2, gcc/g++/gfortran >= 4.8.2 (or clang 3.4+), Qt5 (Widgets, SerialPort, Multimedia, PrintSupport, Sql, WebSockets, LinguistTools), Boost >= 1.62 (log), FFTW3 (single-precision + threads), libusb, and optionally OpenMP.

### Superbuild (builds hamlib + wsjtx together)

```bash
mkdir build && cd build
cmake -DWSJT_SKIP_MANPAGES=ON -DWSJT_GENERATE_DOCS=OFF ..
cmake --build .
cmake --build . --target install    # default prefix: /usr/local
cmake --build . --target package    # CPack package
```

### Building WSJT-X directly (if hamlib is already installed)

```bash
mkdir build && cd build
cmake -DCMAKE_PREFIX_PATH=/path/to/hamlib \
      -DWSJT_SKIP_MANPAGES=ON -DWSJT_GENERATE_DOCS=OFF \
      ../wsjtx
cmake --build .
```

### Useful CMake options

| Option | Default | Description |
|--------|---------|-------------|
| `WSJT_BUILD_UTILS` | ON | Build simulator/codec utilities (ft8sim, jt65sim, etc.) |
| `WSJT_SKIP_MANPAGES` | OFF | Skip manpage generation |
| `WSJT_GENERATE_DOCS` | ON | Build documentation |
| `WSJT_ENABLE_EXPERIMENTAL_FEATURES` | ON (Debug) | Enable experimental features |
| `WSJT_FOX_OTP` | ON | Enable Fox OTP verification |
| `CMAKE_BUILD_TYPE` | RELEASE | `Debug` enables bounds-checking and backtraces in Fortran |

### Running tests

```bash
cd build
ctest                              # runs all tests
ctest -R test_qt_helpers           # run specific test
```

The test suite is minimal — only `tests/test_qt_helpers.cpp` (Qt5Test-based, tests `wsjtx/qt_helpers`).

## WSJT-X Architecture

### Build artifacts — static libraries and executables

The build produces several static libraries linked into the final executables:

- **`wsjt_cxx`** — C/C++ core: CRC, Reed-Solomon (Ka9Q), QRA codes, LDPC tables, math utilities
- **`wsjt_fort`** / **`wsjt_fort_omp`** — Fortran DSP/decoding library (with/without OpenMP)
- **`wsjt_qt`** — Qt integration: transceiver control, network messages, UI models, config, widgets
- **`wsjt_qtmm`** — Qt Multimedia: BWF audio file handling
- **`fort_qt`** — Shared memory bridge (`lib/shmem.cpp`) for Fortran↔Qt IPC
- **`qcp`** — QCustomPlot widget library (waterfall/spectrum plots)
- **`wsjtx_udp`** — UDP message protocol library for external integrations

Main executables: `wsjtx` (GUI app), `jt9` (decoder subprocess), `wsprd` (WSPR decoder), `message_aggregator` (UDP example app), `udp_daemon`.

### C++ / Fortran interop

The critical boundary between C++ and Fortran is the **shared memory segment** defined in `wsjtx/commons.h` (struct `dec_data_t`), which must stay in sync with `wsjtx/lib/jt9com.f90`. The main `wsjtx` process writes audio samples and parameters into shared memory; the `jt9` subprocess reads them and runs Fortran decoders.

Additional Fortran COMMON blocks: `spectra_`, `echocom_`, `foxcom_` (defined in `commons.h`).

Fortran functions are called from C++ via `extern "C"` declarations (e.g., `four2a_` in `main.cpp`). The `FortranCInterface` CMake module generates `FC.h` for name-mangling compatibility.

### Key source directories (under `wsjtx/`)

| Directory | Role |
|-----------|------|
| `widgets/` | Qt UI — `mainwindow.cpp` (17k lines, main app window), widegraph, echograph, fastgraph, logqso, astro, etc. |
| `lib/` | Fortran DSP core — 400+ files for signal processing, encoding/decoding |
| `lib/ft8/`, `lib/ft4/`, `lib/fst4/`, `lib/superfox/` | Mode-specific Fortran codecs |
| `lib/ft8var/` | FT8 variant/experimental decoder code |
| `lib/77bit/` | Message packing/unpacking for 77-bit protocol (packjt77.f90) |
| `lib/qra/` | Q65 and QRA code implementation (C + Fortran) |
| `lib/ftrsd/` | Reed-Solomon FT codec (Ka9Q-derived C code) |
| `lib/wsprd/` | WSPR decoder (standalone C program) |
| `Audio/` | Sound card I/O (soundin/soundout via Qt Multimedia), BWF file support |
| `Detector/` | Audio signal detector (feeds samples to DSP) |
| `Modulator/` | TX audio waveform generation |
| `Decoder/` | Decoded text parsing (`decodedtext.cpp`) |
| `Network/` | UDP protocol (`MessageClient`/`NetworkMessage`), PSKReporter, LotW, Cloudlog, FoxVerifier |
| `Transceiver/` | Rig control abstraction — `TransceiverFactory` creates backends: Hamlib, HRD, DXLab Commander, OmniRig (Win), TCI (WebSocket) |
| `models/` | Qt model classes — Bands, Modes, FrequencyList, StationList, FoxLog, CabrilloLog, DecodeHighlighting |
| `logbook/` | AD1C country data, WorkedBefore tracking, contest multiplier |
| `WSPR/` | WSPR band hopping and TX scheduling |
| `UDPExamples/` | Example UDP client apps (message_aggregator, udp_daemon) |
| `Configuration.cpp/.ui` | Settings dialog (5.6k lines) — all user preferences |
| `translations/` | Qt i18n files (.ts) for ca, da, en, es, it, ja, ru, zh, zh_TW |

### Data flow

1. **Audio in**: `soundin` → `Detector` → samples written to shared memory (`dec_data.d2[]`)
2. **Decoding**: `jt9` subprocess reads shared memory, runs Fortran decoders (ft8_decode, ft4_decode, jt65_decode, q65_decode, etc.), writes results back
3. **Display**: `MainWindow` reads decoded messages, displays in text widgets, updates waterfall (`plotter`)
4. **TX**: User composes message → Fortran encoder generates waveform → `Modulator` → `soundout`
5. **Rig control**: `TransceiverFactory` → backend (HamlibTransceiver, etc.) → serial/network to radio
6. **Reporting**: Decoded spots sent to PSKReporter and optionally Cloudlog/WSPRnet via `Network/`

### Hamlib (under `hamlib-4.7/`)

C library providing a unified API for controlling amateur radio transceivers, rotators, and amplifiers. WSJT-X uses it primarily through `Transceiver/HamlibTransceiver.cpp`. Key directories: `rigs/` (per-manufacturer backends like `icom/`, `kenwood/`, `yaesu/`), `rotators/`, `amplifiers/`, `include/` (public API in `hamlib/rig.h`), `src/` (core implementation). Built as a static library via autotools (`./configure && make`).

## Coding Conventions

- C++11 standard (`-std=gnu++11` on Linux, `-std=c++11 -stdlib=libc++` on macOS)
- C++ compiled with `-Werror -Wall -Wextra`
- Fortran uses `-Wall -Wno-conversion -fno-second-underscore`; Debug adds `-fbounds-check -fbacktrace -fcheck=all`
- Qt5 with `CMAKE_AUTOMOC ON` and `CMAKE_AUTOUIC ON`
- UI files (`.ui`) in `widgets/` and top-level `Configuration.ui`
- Pimpl pattern used in some classes (see `pimpl_h.hpp` / `pimpl_impl.hpp`)
- Fortran modules listed before non-module sources in CMakeLists.txt (build order matters)
- `wsjtx_config.h` is force-included in all C/C++ translation units
- Shared memory struct in `commons.h` must match `lib/jt9com.f90` exactly

## Superbuild Notes

- The top-level `CMakeLists.txt` uses `ExternalProject_Add` to build hamlib (autotools) then wsjtx (CMake)
- Dependency chain: `hamlib-install` → `wsjtx-configure` → `wsjtx-build`
- The superbuild does not work on Windows; macOS support is incomplete
- Extra CMake args pass through to the inner WSJT-X build via command line or `WSJTX_EXTRA_CMAKE_OPTIONS`
