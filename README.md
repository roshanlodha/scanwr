# scGUI (scanwr), formerly scAnWr

scGUI is a native macOS app (SwiftUI) for building and running single‑cell analysis pipelines powered by [Scanpy](https://scanpy.readthedocs.io/). It uses a bundled Python JSON‑RPC backend so end users can run analyses without setting up a local Python environment.

## What scGUI does

- **Metadata management**: add/import sample rows (`sample`, `group`, `path`, optional reader override).
- **Pipeline builder**: drag & drop modules onto a canvas and connect them into a workflow.
- **Multi-sample runs**: run the pipeline across a cohort of samples, with per-sample progress and logs.
- **Incremental re-runs**: caches step signatures and resumes from checkpoints when possible.
- **Visualization**: interactive plotting UI backed by seaborn/matplotlib; export plots as SVG.
- **Project = a folder on disk**: scGUI projects are portable and easy to inspect/version.

## Project format on disk

Each scGUI “project” is just a folder on disk. The app stores state in a single hidden directory:

- `Project/.scanwr/metadata.txt`: sample table (sample, group, path, reader)
- `Project/.scanwr/template.json`: current pipeline canvas workflow
- `Project/.scanwr/checkpoints/`: per-sample `.h5ad` outputs
- `Project/.scanwr/history/`: per-sample cached step signatures (for incremental re-runs)
- `Project/plots/`: plots (organized as `plots/{sample}/*.png`)
- `Project/.scanwr/templates/`: optional workflow templates

## Repo layout

- `mac/ScanwrMac/`: the macOS SwiftUI app + bundled Python RPC server.
- `mac/ScanwrMac/Sources/ScanwrMacApp/Resources/scanwr_rpc_server.py`: the backend’s module registry + execution.
- `requirements-lock.txt`: pinned Python deps for the backend/runtime.

## Release packaging

- Current version: `0.2.5`
- Build outputs: `mac/ScanwrMac/dist/scGUI.app` and `mac/ScanwrMac/dist/scGUI-0.2.5.dmg`

## Development & contributing

Development/build instructions, packaging details, and contributing notes live in `mac/ScanwrMac/README.md`.
