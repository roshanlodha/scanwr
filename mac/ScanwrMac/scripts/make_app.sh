#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DIST="$ROOT/dist"
APP_NAME="scanwr"
BUNDLE_ID="com.roshanlodha.scanwr"
VERSION="0.0.4"

PY_RUNTIME_DIR="${SCANWR_PY_RUNTIME_DIR:-$DIST/python-runtime/python}"

if [[ -z "$PY_RUNTIME_DIR" ]]; then
  echo "ERROR: Set SCANWR_PY_RUNTIME_DIR to a relocatable Python runtime directory to bundle." >&2
  echo "Expected layout: \$SCANWR_PY_RUNTIME_DIR/bin/python3 (and its adjacent libs)." >&2
  echo "Example: export SCANWR_PY_RUNTIME_DIR=/path/to/python-runtime" >&2
  exit 2
fi

if [[ ! -x "$PY_RUNTIME_DIR/bin/python3" ]]; then
  echo "ERROR: Not executable: $PY_RUNTIME_DIR/bin/python3" >&2
  exit 2
fi

mkdir -p "$DIST"

echo "Building release binary via SwiftPM…"
cd "$ROOT"
export TMPDIR="${TMPDIR:-/tmp}"
export CLANG_MODULE_CACHE_PATH="${CLANG_MODULE_CACHE_PATH:-/tmp/scanwr-clang-cache}"
rm -rf "$TMPDIR/scanwr-swift-scratch" "$TMPDIR/scanwr-swift-cache" || true
swift build -c release \
  --disable-sandbox \
  --scratch-path "$TMPDIR/scanwr-swift-scratch" \
  --cache-path "$TMPDIR/scanwr-swift-cache" \
  --manifest-cache local

# For Apple Silicon only (arm64), SwiftPM emits into this bin dir when using the scratch path above.
BIN="$TMPDIR/scanwr-swift-scratch/arm64-apple-macosx/release/ScanwrMacApp"

if [[ ! -f "$BIN" ]]; then
  echo "ERROR: Build output not found: $BIN" >&2
  exit 3
fi

APP="$DIST/$APP_NAME.app"
CONTENTS="$APP/Contents"
MACOS="$CONTENTS/MacOS"
RES="$CONTENTS/Resources"

rm -rf "$APP"
mkdir -p "$MACOS" "$RES"

echo "Creating app bundle at: $APP"
cp "$BIN" "$MACOS/$APP_NAME"

# Put the python server script in Resources (Swift looks in Bundle.main first).
cp "$ROOT/Sources/ScanwrMacApp/Resources/scanwr_rpc_server.py" "$RES/scanwr_rpc_server.py"

echo "Bundling Python runtime…"
rm -rf "$RES/python"
mkdir -p "$RES"
ditto --noqtn "$PY_RUNTIME_DIR" "$RES/python"

# Minimal Info.plist
cat > "$CONTENTS/Info.plist" <<EOF
<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE plist PUBLIC "-//Apple//DTD PLIST 1.0//EN" "http://www.apple.com/DTDs/PropertyList-1.0.dtd">
<plist version="1.0">
<dict>
  <key>CFBundleDevelopmentRegion</key>
  <string>en</string>
  <key>CFBundleExecutable</key>
  <string>$APP_NAME</string>
  <key>CFBundleIdentifier</key>
  <string>$BUNDLE_ID</string>
  <key>CFBundleInfoDictionaryVersion</key>
  <string>6.0</string>
  <key>CFBundleName</key>
  <string>$APP_NAME</string>
  <key>CFBundlePackageType</key>
  <string>APPL</string>
  <key>CFBundleShortVersionString</key>
  <string>$VERSION</string>
  <key>CFBundleVersion</key>
  <string>$VERSION</string>
  <key>LSMinimumSystemVersion</key>
  <string>14.0</string>
  <key>NSHighResolutionCapable</key>
  <true/>
</dict>
</plist>
EOF

echo "OK: $APP"
