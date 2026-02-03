#!/usr/bin/env bash
set -euo pipefail

# Build a self-contained Python runtime suitable for bundling into scanwr.app.
#
# This script intentionally runs at *build time* (developer machine), not on the end user's machine.
# End users should only download the .dmg, drag the app, and run it.
#
# Inputs:
# - PYTHON_BS_TARBALL: path to a python-build-standalone tarball for macOS arm64
#   OR
# - PYTHON_BS_URL: URL to download that tarball
#
# Output:
# - dist/python-runtime/python/bin/python3  (and its libs + site-packages)
#
# Then package with:
#   export SCANWR_PY_RUNTIME_DIR="$(pwd)/dist/python-runtime/python"
#   ./scripts/make_app.sh

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
OUT="${SCANWR_PY_RUNTIME_OUT:-$ROOT/dist/python-runtime}"

PYTHON_BS_TARBALL="${PYTHON_BS_TARBALL:-}"
PYTHON_BS_URL="${PYTHON_BS_URL:-}"

# Prefer a lockfile (best for reproducible packaging).
# Default points at repo-root requirements-lock.txt if present.
REPO_ROOT="$(cd "$ROOT/../.." && pwd)"
DEFAULT_REQ="$REPO_ROOT/requirements-lock.txt"
SCANWR_PY_REQUIREMENTS="${SCANWR_PY_REQUIREMENTS:-}"
if [[ -z "$SCANWR_PY_REQUIREMENTS" && -f "$DEFAULT_REQ" ]]; then
  SCANWR_PY_REQUIREMENTS="$DEFAULT_REQ"
fi

# Build mode:
# - venv: copy site-packages from a prebuilt venv (offline-friendly; best in this repo)
# - requirements: pip install -r requirements (needs network or a local wheelhouse)
DEFAULT_MODE="requirements"
DEFAULT_VENV_PY="$REPO_ROOT/venv/bin/python"
if [[ -x "$DEFAULT_VENV_PY" ]]; then
  DEFAULT_MODE="venv"
fi
SCANWR_PY_MODE="${SCANWR_PY_MODE:-$DEFAULT_MODE}"
SCANWR_VENV_PYTHON="${SCANWR_VENV_PYTHON:-$DEFAULT_VENV_PY}"

# Fallback: space-separated pip specs to install into the embedded runtime.
# Examples:
#   export SCANWR_PY_PACKAGES="scanpy==1.12.* squidpy==1.8.*"
#   export SCANWR_PY_PACKAGES="scanpy squidpy"
SCANWR_PY_PACKAGES="${SCANWR_PY_PACKAGES:-scanpy squidpy}"

mkdir -p "$ROOT/dist"
rm -rf "$OUT"
mkdir -p "$OUT"

_need() {
  command -v "$1" >/dev/null 2>&1 || { echo "ERROR: missing dependency: $1" >&2; exit 2; }
}

_need tar

TARBALL=""
if [[ -n "$PYTHON_BS_TARBALL" ]]; then
  TARBALL="$PYTHON_BS_TARBALL"
elif [[ -z "$PYTHON_BS_URL" ]]; then
  # Convenience: if the repo already contains a suitable standalone tarball, use it.
  # Example: cpython-3.12.12+20260127-aarch64-apple-darwin-install_only.tar
  CANDIDATE="$(ls -1 "$REPO_ROOT"/cpython-*-aarch64-apple-darwin-install_only.tar 2>/dev/null | head -n 1 || true)"
  if [[ -n "$CANDIDATE" ]]; then
    TARBALL="$CANDIDATE"
  fi
elif [[ -n "$PYTHON_BS_URL" ]]; then
  _need curl
  TARBALL="$ROOT/dist/python-build-standalone.tar.zst"
  echo "Downloading python-build-standalone…"
  curl -L "$PYTHON_BS_URL" -o "$TARBALL"
else
  echo "ERROR: set PYTHON_BS_TARBALL or PYTHON_BS_URL" >&2
  exit 2
fi

if [[ ! -f "$TARBALL" ]]; then
  echo "ERROR: tarball not found: $TARBALL" >&2
  exit 2
fi

echo "Extracting Python runtime…"
# Strip quarantine so embedded dylibs can load.
xattr -dr com.apple.quarantine "$TARBALL" >/dev/null 2>&1 || true
if [[ "$TARBALL" == *.tar.zst ]]; then
  _need zstd
  # Many python-build-standalone releases are .tar.zst
  zstd -dc "$TARBALL" | tar -x -C "$OUT"
else
  tar -x -C "$OUT" -f "$TARBALL"
fi
xattr -dr com.apple.quarantine "$OUT/python" >/dev/null 2>&1 || true

# We expect the extracted archive to contain a top-level 'python' directory.
if [[ ! -x "$OUT/python/bin/python3" ]]; then
  echo "ERROR: expected extracted python at: $OUT/python/bin/python3" >&2
  echo "Contents:" >&2
  ls -la "$OUT" >&2 || true
  exit 3
fi

PY="$OUT/python/bin/python3"

echo "Installing packages into embedded runtime…"
if [[ "$SCANWR_PY_MODE" == "venv" ]]; then
  if [[ ! -x "$SCANWR_VENV_PYTHON" ]]; then
    echo "ERROR: SCANWR_VENV_PYTHON not found/executable: $SCANWR_VENV_PYTHON" >&2
    exit 2
  fi
  echo "Mode=venv (offline-friendly)"
  echo "Copying site-packages from: $SCANWR_VENV_PYTHON"

  RT_PURELIB="$("$PY" -c "import sysconfig; print(sysconfig.get_paths()['purelib'])")"
  RT_PLATLIB="$("$PY" -c "import sysconfig; print(sysconfig.get_paths()['platlib'])")"
  VENV_PURELIB="$("$SCANWR_VENV_PYTHON" -c "import sysconfig; print(sysconfig.get_paths()['purelib'])")"
  VENV_PLATLIB="$("$SCANWR_VENV_PYTHON" -c "import sysconfig; print(sysconfig.get_paths()['platlib'])")"

  echo "Runtime purelib: $RT_PURELIB"
  echo "Venv purelib:    $VENV_PURELIB"

  _need rsync

  # Clear the runtime's site-packages and replace with the venv's.
  if [[ "$RT_PURELIB" == "$RT_PLATLIB" ]]; then
    rm -rf "$RT_PURELIB"
    mkdir -p "$RT_PURELIB"
    rsync -a "$VENV_PURELIB/" "$RT_PURELIB/"
  else
    rm -rf "$RT_PURELIB" "$RT_PLATLIB"
    mkdir -p "$RT_PURELIB" "$RT_PLATLIB"
    rsync -a "$VENV_PURELIB/" "$RT_PURELIB/"
    rsync -a "$VENV_PLATLIB/" "$RT_PLATLIB/"
  fi

  echo "Verifying imports…"
  "$PY" -c "import scanpy as sc, squidpy as sq; print('scanpy', getattr(sc,'__version__','?')); print('squidpy', getattr(sq,'__version__','?'))"
elif [[ "$SCANWR_PY_MODE" == "requirements" ]]; then
  echo "Bootstrapping pip…"
  "$PY" -m ensurepip --upgrade >/dev/null 2>&1 || true
  "$PY" -m pip install --upgrade pip setuptools wheel

  if [[ -n "$SCANWR_PY_REQUIREMENTS" ]]; then
    if [[ ! -f "$SCANWR_PY_REQUIREMENTS" ]]; then
      echo "ERROR: SCANWR_PY_REQUIREMENTS not found: $SCANWR_PY_REQUIREMENTS" >&2
      exit 2
    fi
    echo "Mode=requirements"
    echo "Using requirements file: $SCANWR_PY_REQUIREMENTS"
    "$PY" -m pip install -r "$SCANWR_PY_REQUIREMENTS"
  else
    echo "Mode=requirements"
    echo "Using SCANWR_PY_PACKAGES=$SCANWR_PY_PACKAGES"
    # shellcheck disable=SC2086
    "$PY" -m pip install $SCANWR_PY_PACKAGES
  fi
else
  echo "ERROR: Unknown SCANWR_PY_MODE: $SCANWR_PY_MODE (use 'venv' or 'requirements')" >&2
  exit 2
fi

echo "OK: Built embedded runtime at: $OUT/python"
echo "Next:"
echo "  export SCANWR_PY_RUNTIME_DIR=\"$OUT/python\""
echo "  cd \"$ROOT\" && ./scripts/make_app.sh"
