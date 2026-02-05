# scGUI (scanwr), formerly scAnWr

scGUI is a native macOS app (SwiftUI) for building and running single‑cell analysis pipelines powered by [Scanpy](https://scanpy.readthedocs.io/). It uses a bundled Python JSON‑RPC backend so end users can run analyses without setting up a local Python environment.

## What scGUI does

- **Metadata management**: add/import sample rows (`sample`, `group`, `path`, optional reader override).
- **Pipeline builder**: drag & drop modules onto a canvas and connect them into a workflow.
- **Multi-sample (cohort) analysis**: loads all samples, concatenates them, and runs the pipeline once with `adata.obs['sample']` and `adata.obs['group']` attached.
- **Incremental re-runs**: caches step signatures and resumes from checkpoints when possible.
- **Visualization**: interactive plotting UI backed by seaborn/matplotlib; export plots as SVG.
- **Project = a folder on disk**: scGUI projects are portable and easy to inspect/version.

## Project format on disk

Each scGUI “project” is just a folder on disk. The app stores state in a single hidden directory:

- `Project/.scanwr/metadata.txt`: sample table (sample, group, path, reader)
- `Project/.scanwr/template.json`: current pipeline canvas workflow
- `Project/.scanwr/checkpoints/`: analysis outputs (e.g. `cohort.h5ad` for concatenated multi-sample runs)
- `Project/.scanwr/history/`: cached step signatures (for incremental re-runs)
- `Project/plots/`: plots (organized as `plots/{sample}/*.png`; cohort runs use `plots/__cohort__/`)
- `Project/.scanwr/templates/`: optional workflow templates

## Repo layout

- `mac/ScanwrMac/`: the macOS SwiftUI app + bundled Python RPC server.
- `mac/ScanwrMac/Sources/ScanwrMacApp/Resources/scanwr_rpc_server.py`: the backend’s module registry + execution.
- `requirements-lock.txt`: pinned Python deps for the backend/runtime.

## Multi-sample marker genes (groups)

When you run a multi-sample project, scGUI concatenates all samples into one AnnData and adds:

- `adata.obs["sample"]`: the sample id from metadata
- `adata.obs["group"]`: the group label from metadata

To find marker genes between groups, add **Rank Genes Groups** and set `groupby = group` (this runs `scanpy.tl.rank_genes_groups` and saves results into `adata.uns["rank_genes_groups"]`).

You can switch between **Cohort (concat)** and **Per-sample** runs using the toggle at the top of the Pipeline panel.

## Release packaging

- Current version: `0.3.1`
- Build outputs: `mac/ScanwrMac/dist/scGUI.app` and `mac/ScanwrMac/dist/scGUI-0.3.1.dmg`

## Development & contributing

Development/build instructions, packaging details, and contributing notes live in `mac/ScanwrMac/README.md`.
