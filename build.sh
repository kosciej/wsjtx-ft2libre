#!/bin/bash
set -e

# 1. Clean previous build
echo "Cleaning old build..."
rm -rf build

# 2. Setup build directory
mkdir build
cd build

# 3. Configure Superbuild
echo "Configuring project..."
# We use the Qt5 path from Homebrew
QT5_PATH=$(brew --prefix qt@5)
cmake .. -DCMAKE_PREFIX_PATH="$QT5_PATH"

# 4. Build
echo "Starting build (this may take a few minutes)..."
cmake --build .

# 5. Fix macOS App Bundle
echo "Setting up wsjtx.app bundle links..."
BUILD_ROOT=".."
APP_DIR="wsjtx-prefix/src/wsjtx-build/wsjtx.app"
MACOS_DIR="$APP_DIR/Contents/MacOS"
RES_DIR="$APP_DIR/Contents/Resources"

# Link helper binaries
ln -sf "../../../jt9" "$MACOS_DIR/jt9"
ln -sf "../../../../../hamlib-prefix/bin/rigctl" "$MACOS_DIR/rigctl-wsjtx"
ln -sf "../../../../../hamlib-prefix/bin/rigctld" "$MACOS_DIR/rigctld-wsjtx"
ln -sf "../../../../../hamlib-prefix/bin/rigctlcom" "$MACOS_DIR/rigctlcom-wsjtx"

# Link resource files
ln -sf "../../../../../wsjtx/cty.dat" "$RES_DIR/cty.dat"
ln -sf "../../../../../wsjtx/CALL3.TXT" "$RES_DIR/CALL3.TXT"
ln -sf "../../../../../wsjtx/ALLCALL7.TXT" "$RES_DIR/ALLCALL7.TXT"

# Link translation files
for f in wsjtx-prefix/src/wsjtx-build/*.qm; do
  if [ -f "$f" ]; then
    ln -sf "../../../$(basename "$f")" "$RES_DIR/$(basename "$f")"
  fi
done

echo "------------------------------------------------"
echo "Build complete!"
echo "To run the application, use:"
echo "open build/$APP_DIR"
echo "------------------------------------------------"
