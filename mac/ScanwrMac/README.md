# scGUI (macOS app)

scGUI, formerly scAnWr, is a native SwiftUI macOS front-end that talks to a Python JSON-RPC backend.

## Notes

- Project state is persisted under `Project/.scanwr/` (including the current canvas workflow in `template.json`).
- Multi-sample runs are concatenated into one AnnData (`Project/.scanwr/checkpoints/cohort.h5ad`) with `adata.obs["sample"]` and `adata.obs["group"]` populated from metadata.
- The Pipeline panel includes a **Cohort (concat)** toggle to switch between concatenated cohort analysis and per-sample execution.
- The Welcome screen includes buttons for Settings (verbosity slider) and clearing app cache.
- Custom modules are exposed as `custom.<module>` (orange badge in the UI). `custom.select_group` filters cells by `adata.obs[group] == value` (e.g. `group=leiden`, `value=3` keeps cluster "3").

## Run (development)

1. Ensure your Python env has `scanpy` installed (your repo `venv/` already does).
2. Build/run with an env var pointing at that Python:

```bash
cd mac/ScanwrMac
SCANWR_PYTHON=/Users/roshanlodha/Documents/scanwr/venv/bin/python swift run ScanwrMacApp
```

## Packaging (.app + .dmg)

The app bundle id is `com.roshanlodha.scanwr`.

Current release version: `0.3.1`

Note: Leiden clustering (e.g. `scanpy.tl.leiden`) may require optional Python deps `leidenalg` + `igraph`.
The embedded runtime build script installs these automatically; if you maintain your own `venv/`,
ensure they are installed there too.

### 1) Provide a relocatable Python runtime

You must provide a self-contained Python runtime directory with:

- `bin/python3`
- its adjacent `lib/` (or framework) so it runs without Homebrew/Python.org installs
- `site-packages` containing `scanpy` and all dependencies

Set it as:

```bash
export SCANWR_PY_RUNTIME_DIR=/path/to/python-runtime
```

If you use `python-build-standalone`, you can build this automatically:

```bash
cd mac/ScanwrMac
# Choose either:
export PYTHON_BS_TARBALL=/path/to/python-build-standalone.tar.zst
# or:
export PYTHON_BS_URL="https://…/python-build-standalone-…-macos-arm64-….tar.zst"
# Prefer the repo lockfile:
export SCANWR_PY_REQUIREMENTS="/Users/roshanlodha/Documents/scanwr/requirements-lock.txt"
# Or install ad-hoc:
# export SCANWR_PY_PACKAGES="scanpy==1.12.* squidpy==1.8.*"
./scripts/build_python_runtime.sh
```

### Offline-friendly (recommended for this repo)

If you already have `scanpy` + `squidpy` installed in the repo `venv/`, you can build the embedded runtime
without downloading any wheels by copying the venv’s `site-packages`:

```bash
cd mac/ScanwrMac
export PYTHON_BS_TARBALL="/Users/roshanlodha/Documents/scanwr/cpython-3.12.12+20260127-aarch64-apple-darwin-install_only.tar"
export SCANWR_PY_MODE="venv"
./scripts/build_python_runtime.sh
```

### 2) Build the `.app`

```bash
cd mac/ScanwrMac
./scripts/make_app.sh
```

### 3) Create the `.dmg`

```bash
cd mac/ScanwrMac
./scripts/make_dmg.sh
```

Outputs land in `mac/ScanwrMac/dist/` (e.g. `scGUI-0.3.1.dmg`).
