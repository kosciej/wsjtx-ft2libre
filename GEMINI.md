# WSJT-X 3.0.0-rc1 & FT2 Research

## Project Overview

This workspace contains the source code for **WSJT-X 3.0.0-rc1**, a weak-signal amateur radio communication software, along with a dedicated **research environment** for reverse-engineering the **FT2** digital protocol.

The project is structured as a **CMake superbuild** that orchestrates the build of:
1.  **Hamlib** (Rig control library, version 4.7)
2.  **WSJT-X** (The main application)

Additionally, the `research/` directory contains a Python-based analysis suite used to characterize the FT2 protocol.

## Directory Structure

*   **`wsjtx/`**: Main C++/Fortran source code for WSJT-X.
    *   `CMakeLists.txt`: Main build configuration.
    *   `src/`: Source files.
*   **`hamlib-4.7/`**: Source code for the Hamlib library.
*   **`research/`**: Python environment for FT2 signal analysis.
    *   `README.md`: **CRITICAL**. Detailed log of the reverse-engineering process and protocol findings.
    *   `protocol.md`: Technical specification of the FT2 protocol.
    *   `*.py`: Analysis scripts (e.g., `ft2_raw_analysis.py`, `analyze_capture.py`).
    *   `*.wav`: Audio captures of FT2 signals.
*   **`build/`**: Directory for build artifacts (contains `wsjtx-build`, `hamlib-build`).
*   **`CMakeLists.txt`**: The root superbuild script.

## Build Instructions

### Prerequisites
*   CMake
*   Qt 5 (Components: Core, Gui, Widgets, Network, PrintSupport, etc.)
*   Fortran Compiler (gfortran)
*   C++ Compiler (Clang/GCC)
*   Python 3 (for research scripts)

### Building WSJT-X (Superbuild)

The root `CMakeLists.txt` is a superbuild wrapper.

```bash
mkdir build
cd build
cmake ..
cmake --build .
```

This process will:
1.  Build `hamlib` (static).
2.  Build `wsjtx` linked against the local `hamlib`.

### Research Environment (Python)

The `research/` directory uses a `pyproject.toml` for dependencies.

```bash
cd research
# Install dependencies (using pip or uv)
pip install .
# OR if using uv
uv sync
```

**Key Python Dependencies:**
*   `numpy`
*   `scipy`
*   `matplotlib`
*   `soundfile`

## Research Context: FT2 Protocol

**Objective:** Implement FT2 support in WSJT-X by reverse-engineering the protocol from the closed-source "Decodium" application and on-air signals.

**Key Findings (from `research/README.md`):**
*   **Protocol:** FT2 is **4-GFSK** (not 8-GFSK as originally claimed).
*   **Relation:** It is essentially **FT4 running at 2x speed**.
*   **Symbol Rate:** 41.667 baud (vs 20.833 for FT4).
*   **Tones:** 4 tones spaced ~41.6 Hz apart.
*   **NSPS:** 288 samples per symbol (at 12kHz).

**Important Scripts:**
*   `research/ft2_raw_analysis.py`: The script that revealed the 4-tone structure.
*   `research/analyze_capture.py`: General signal analysis.

## Development Conventions

*   **Languages:** C++ (Qt framework) for GUI/Control, Fortran for DSP/Modulation/Demodulation.
*   **Style:** Follow existing Qt-style C++ naming (CamelCase for classes, mixedCase for variables).
*   **Documentation:** Maintain `research/README.md` with new findings.
*   **Git:** Do not commit large `.wav` files or build artifacts.
