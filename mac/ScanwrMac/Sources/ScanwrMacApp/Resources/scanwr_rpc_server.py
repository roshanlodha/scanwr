#!/usr/bin/env python3
from __future__ import annotations

import json
import os
import sys
import tempfile
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple


def _set_safe_env() -> None:
    # Scanpy pulls matplotlib/fontconfig; route caches to writable places.
    os.environ.setdefault("MPLCONFIGDIR", tempfile.mkdtemp(prefix="scanwr-mpl-"))
    os.environ.setdefault("XDG_CACHE_HOME", tempfile.gettempdir())
    os.environ.setdefault("PYTHONUNBUFFERED", "1")


def _notify_log(message: str) -> None:
    sys.stdout.write(json.dumps({"method": "log", "params": {"message": message}}) + "\n")
    sys.stdout.flush()


def _notify_progress(
    percent: float,
    message: str,
    sample: str | None = None,
    step_index: int | None = None,
    step_count: int | None = None,
) -> None:
    params: Dict[str, Any] = {"percent": max(0.0, min(1.0, float(percent))), "message": message}
    if sample is not None:
        params["sample"] = sample
    if step_index is not None:
        params["stepIndex"] = step_index
    if step_count is not None:
        params["stepCount"] = step_count
    sys.stdout.write(json.dumps({"method": "progress", "params": params}) + "\n")
    sys.stdout.flush()


def _is_10x_mtx_dir(path: Path) -> bool:
    required_any = {"matrix.mtx", "matrix.mtx.gz"}
    barcodes_any = {"barcodes.tsv", "barcodes.tsv.gz"}
    features_any = {"features.tsv", "features.tsv.gz", "genes.tsv", "genes.tsv.gz"}
    names = {p.name for p in path.iterdir() if p.is_file()}
    return bool(names & required_any) and bool(names & barcodes_any) and bool(names & features_any)


def _find_first(path: Path, candidates: List[str]) -> Optional[Path]:
    for name in candidates:
        p = path / name
        if p.exists():
            return p
    return None


def _detect_reader(path: str) -> Tuple[str, str]:
    p = Path(path).expanduser().resolve()
    if not p.exists():
        raise FileNotFoundError(str(p))

    if p.is_file():
        suffix = "".join(p.suffixes).lower()
        if suffix.endswith(".h5ad"):
            return ("scanpy.read_h5ad", str(p))
        if suffix.endswith(".loom"):
            return ("scanpy.read_loom", str(p))
        if suffix.endswith(".mtx") or suffix.endswith(".mtx.gz"):
            return ("scanpy.read_mtx", str(p))
        return ("scanpy.read", str(p))

    tenx_h5 = _find_first(
        p,
        [
            "filtered_feature_bc_matrix.h5",
            "raw_feature_bc_matrix.h5",
            "filtered_gene_bc_matrices_h5.h5",
        ],
    )
    if tenx_h5 is not None:
        return ("scanpy.read_10x_h5", str(tenx_h5))

    if _is_10x_mtx_dir(p):
        return ("scanpy.read_10x_mtx", str(p))

    h5ads = [x for x in p.iterdir() if x.is_file() and x.name.lower().endswith(".h5ad")]
    if len(h5ads) == 1:
        return ("scanpy.read_h5ad", str(h5ads[0]))

    raise ValueError(f"Unsupported input path; cannot detect reader for: {p}")


def detect_reader(path: str) -> Dict[str, Any]:
    reader, arg = _detect_reader(path)
    reason = "auto-detected from file extension" if Path(path).is_file() else "auto-detected from directory contents"
    return {"suggested": reader, "reason": reason}


def _read_with_reader(reader: str, path: str):
    import scanpy as sc  # local import after env set

    p = Path(path).expanduser().resolve()
    if not p.exists():
        raise FileNotFoundError(str(p))

    # If user didn't override, detect.
    if not reader:
        reader, arg = _detect_reader(str(p))
    else:
        arg = str(p)

    if reader == "scanpy.read_h5ad":
        return reader, sc.read_h5ad(arg)
    if reader == "scanpy.read_loom":
        return reader, sc.read_loom(arg)
    if reader == "scanpy.read_mtx":
        return reader, sc.read_mtx(arg)
    if reader == "scanpy.read_10x_h5":
        return reader, sc.read_10x_h5(arg)
    if reader == "scanpy.read_10x_mtx":
        return reader, sc.read_10x_mtx(arg, var_names="gene_symbols")
    if reader == "scanpy.read":
        return reader, sc.read(arg)

    raise ValueError(f"Unsupported reader override: {reader}")


def _read_any(path: str):
    import scanpy as sc  # local import after env set

    reader, arg = _detect_reader(path)
    if reader == "scanpy.read_h5ad":
        return reader, sc.read_h5ad(arg)
    if reader == "scanpy.read_loom":
        return reader, sc.read_loom(arg)
    if reader == "scanpy.read_mtx":
        return reader, sc.read_mtx(arg)
    if reader == "scanpy.read_10x_h5":
        return reader, sc.read_10x_h5(arg)
    if reader == "scanpy.read_10x_mtx":
        return reader, sc.read_10x_mtx(arg, var_names="gene_symbols")
    return reader, sc.read(arg)


def _parse_csv_list(value: str) -> List[str]:
    return [x.strip() for x in (value or "").split(",") if x.strip()]


def _parse_int_list(value: str) -> Optional[List[int]]:
    items = _parse_csv_list(value)
    if not items:
        return None
    return [int(x) for x in items]


def _run_module(adata, spec_id: str, params: Dict[str, Any]) -> None:
    import scanpy as sc  # local import after env set

    # Backwards compatibility (older ids)
    if spec_id == "pp.calculate_qc_metrics":
        spec_id = "scanpy.pp.calculate_qc_metrics"

    def _opt_int(key: str) -> Optional[int]:
        v = (params or {}).get(key)
        if v is None:
            return None
        if isinstance(v, (int, float)):
            return int(v)
        if isinstance(v, str):
            if not v.strip():
                return None
            return int(v)
        return int(v)

    if spec_id == "scanpy.pp.filter_cells":
        # scanpy.pp.filter_cells allows only ONE of these thresholds per call.
        for key in ["min_counts", "min_genes", "max_counts", "max_genes"]:
            val = _opt_int(key)
            if val is None:
                continue
            sc.pp.filter_cells(adata, **{key: val}, inplace=True)
        return

    if spec_id == "scanpy.pp.filter_genes":
        # scanpy.pp.filter_genes allows only ONE of these thresholds per call.
        for key in ["min_counts", "min_cells", "max_counts", "max_cells"]:
            val = _opt_int(key)
            if val is None:
                continue
            sc.pp.filter_genes(adata, **{key: val}, inplace=True)
        return

    if spec_id == "scanpy.pp.calculate_qc_metrics":
        expr_type = str(params.get("expr_type", "counts") or "counts")
        var_type = str(params.get("var_type", "genes") or "genes")
        qc_vars = _parse_csv_list(str(params.get("qc_vars", "") or "")) or ()
        percent_top = tuple(_parse_int_list(str(params.get("percent_top", "") or "")) or [50, 100, 200, 500])
        layer_raw = str(params.get("layer", "") or "")
        layer = layer_raw if layer_raw else None
        use_raw = bool(params.get("use_raw", False))
        inplace = bool(params.get("inplace", False))
        log1p = bool(params.get("log1p", True))
        parallel_raw = str(params.get("parallel", "") or "")
        parallel = int(parallel_raw) if parallel_raw.strip() else None

        sc.pp.calculate_qc_metrics(
            adata,
            expr_type=expr_type,
            var_type=var_type,
            qc_vars=qc_vars,
            percent_top=percent_top,
            layer=layer,
            use_raw=use_raw,
            inplace=inplace,
            log1p=log1p,
            parallel=parallel,
        )
        return

    raise NotImplementedError(spec_id)


def list_modules() -> List[Dict[str, Any]]:
    # For the next version, expose only these two preprocessing filters for focused testing.
    return [
        {
            "id": "scanpy.pp.filter_cells",
            "group": "pp",
            "namespace": "core",
            "title": "Filter Cells",
            "scanpyQualname": "scanpy.pp.filter_cells",
        },
        {
            "id": "scanpy.pp.filter_genes",
            "group": "pp",
            "namespace": "core",
            "title": "Filter Genes",
            "scanpyQualname": "scanpy.pp.filter_genes",
        },
    ]


def _sanitize_filename(s: str) -> str:
    out = []
    for ch in s:
        if ch.isalnum() or ch in ("-", "_", "."):
            out.append(ch)
        else:
            out.append("_")
    return "".join(out).strip("_") or "step"


def run_pipeline(input_path: str, output_dir: str, steps: List[Dict[str, Any]]) -> Dict[str, Any]:
    import scanpy as sc  # local import after env set

    in_path = Path(input_path).expanduser().resolve()
    out_dir = Path(output_dir).expanduser().resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    if not in_path.exists():
        raise FileNotFoundError(str(in_path))
    if in_path.suffix.lower() != ".h5ad":
        raise ValueError(f"Expected .h5ad input, got: {in_path.name}")

    _notify_log("Loading input .h5ad…")
    adata = sc.read_h5ad(str(in_path))

    checkpoints: List[str] = []
    for i, step in enumerate(steps, start=1):
        spec_id = str(step.get("specId"))
        params = step.get("params") or {}
        _notify_log(f"Step {i}/{len(steps)}: {spec_id}")
        _run_module(adata, spec_id, params)

        ck_name = f"{i:02d}_{_sanitize_filename(spec_id)}.h5ad"
        ck_path = out_dir / ck_name
        _notify_log(f"Saving checkpoint: {ck_path.name}")
        adata.write_h5ad(str(ck_path))
        checkpoints.append(str(ck_path))

        # Reload from checkpoint to ensure the next step operates on an on-disk checkpointed state.
        adata = sc.read_h5ad(str(ck_path))

    final_path = out_dir / "final.h5ad"
    _notify_log(f"Saving final: {final_path.name}")
    adata.write_h5ad(str(final_path))

    shape = list(getattr(adata, "shape", (0, 0)))
    return {
        "inputPath": str(in_path),
        "outputDir": str(out_dir),
        "checkpoints": checkpoints,
        "finalPath": str(final_path),
        "shape": shape,
    }


def _params_sig(params: Dict[str, Any]) -> str:
    import hashlib

    blob = json.dumps(params, sort_keys=True, separators=(",", ":"), ensure_ascii=True)
    return hashlib.sha256(blob.encode("utf-8")).hexdigest()


def _step_sig(step: Dict[str, Any]) -> Dict[str, Any]:
    spec_id = str(step.get("specId"))
    params = step.get("params") or {}
    return {"specId": spec_id, "paramsSha256": _params_sig(params)}


def run_pipeline_multi(output_dir: str, project_name: str, samples: List[Dict[str, Any]], steps: List[Dict[str, Any]]) -> Dict[str, Any]:
    import scanpy as sc  # local import after env set

    out_dir = Path(output_dir).expanduser().resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    project = out_dir / _sanitize_filename(project_name)
    project.mkdir(parents=True, exist_ok=True)
    scanwr_dir = project / ".scanwr"
    scanwr_dir.mkdir(parents=True, exist_ok=True)
    checkpoints_dir = scanwr_dir / "checkpoints"
    history_dir = scanwr_dir / "history"
    checkpoints_dir.mkdir(parents=True, exist_ok=True)
    history_dir.mkdir(parents=True, exist_ok=True)

    # Persist metadata for reload/debug.
    meta_path = scanwr_dir / "metadata.txt"
    lines = ["sample\tgroup\tpath\treader"]
    for s in samples:
        lines.append(
            "\t".join(
                [
                    str(s.get("sample", "")).strip(),
                    str(s.get("group", "")).strip(),
                    str(s.get("path", "")).strip(),
                    str(s.get("reader", "")).strip(),
                ]
            )
        )
    meta_path.write_text("\n".join(lines) + "\n")

    requested_sigs = [_step_sig(s) for s in steps]

    results: List[Dict[str, Any]] = []
    total_samples = max(1, len(samples))
    for s_idx, s in enumerate(samples, start=1):
        sample = str(s.get("sample", "")).strip()
        group = str(s.get("group", "")).strip()
        path = str(s.get("path", "")).strip()
        reader_override = str(s.get("reader", "")).strip()

        if not sample or not group or not path:
            raise ValueError("Each sample must include sample, group, path")

        checkpoint_h5ad = checkpoints_dir / f"{_sanitize_filename(sample)}.h5ad"
        checkpoint_steps = history_dir / f"{_sanitize_filename(sample)}.json"

        def _load_cached_prefix() -> List[Dict[str, Any]]:
            if not checkpoint_steps.exists():
                return []
            try:
                return json.loads(checkpoint_steps.read_text())
            except Exception:
                return []

        def _write_cached_prefix(prefix: List[Dict[str, Any]]) -> None:
            checkpoint_steps.write_text(json.dumps(prefix, indent=2, sort_keys=True))

        cached_sigs = _load_cached_prefix()
        prefix_len = 0
        for a, b in zip(cached_sigs, requested_sigs):
            if a == b:
                prefix_len += 1
            else:
                break

        if prefix_len > 0 and checkpoint_h5ad.exists():
            _notify_log(f"[{sample}] Cache hit: skipping {prefix_len}/{len(requested_sigs)} step(s)")
            _notify_progress(
                percent=(s_idx - 1) / total_samples,
                message=f"{sample}: resuming from checkpoint",
                sample=sample,
                step_index=prefix_len,
                step_count=len(requested_sigs),
            )
            adata = sc.read_h5ad(str(checkpoint_h5ad))
            reader = reader_override or "auto"
        else:
            _notify_log(f"Reading {sample}…")
            _notify_progress(
                percent=(s_idx - 1) / total_samples,
                message=f"{sample}: reading input",
                sample=sample,
                step_index=0,
                step_count=len(requested_sigs),
            )
            reader, adata = _read_with_reader(reader_override, path)
            # Reset cache if it doesn't match
            _write_cached_prefix([])

        checkpoints: List[str] = []
        for i, step in enumerate(steps, start=1):
            if i <= prefix_len:
                continue
            spec_id = str(step.get("specId"))
            params = step.get("params") or {}

            overall = ((s_idx - 1) + (i / max(1, len(steps)))) / total_samples
            _notify_progress(
                percent=overall,
                message=f"{sample}: {spec_id}",
                sample=sample,
                step_index=i,
                step_count=len(steps),
            )
            _notify_log(f"[{sample}] Step {i}/{len(steps)}: {spec_id}")
            _run_module(adata, spec_id, params)

            _notify_log(f"[{sample}] Updating checkpoint.h5ad")
            adata.write_h5ad(str(checkpoint_h5ad))
            # Update cached signature list in order.
            done = requested_sigs[:i]
            _write_cached_prefix(done)
            checkpoints.append(str(checkpoint_h5ad))
            adata = sc.read_h5ad(str(checkpoint_h5ad))

        final_path = checkpoint_h5ad
        _notify_progress(
            percent=min(1.0, s_idx / total_samples),
            message=f"{sample}: done",
            sample=sample,
            step_index=len(steps),
            step_count=len(steps),
        )

        shape = list(getattr(adata, "shape", (0, 0)))
        results.append(
            {
                "sample": sample,
                "group": group,
                "path": path,
                "reader": reader,
                "outputDir": str(project),
                "checkpoints": checkpoints,
                "finalPath": str(final_path),
                "shape": shape,
            }
        )

    return {"outputDir": str(project), "results": results}


def _handle(method: str, params: Any) -> Any:
    if method == "ping":
        return {"ok": True}
    if method == "list_modules":
        return list_modules()
    if method == "detect_reader":
        return detect_reader(path=params["path"])
    if method == "run_pipeline":
        # New API: per-sample raw input paths
        if "samples" in params:
            return run_pipeline_multi(
                output_dir=params["outputDir"],
                project_name=params.get("projectName") or "scanwr-project",
                samples=params["samples"],
                steps=params["steps"],
            )
        # Legacy (kept for compatibility during iteration)
        return run_pipeline(input_path=params["inputPath"], output_dir=params["outputDir"], steps=params["steps"])
    raise ValueError(f"Unknown method: {method}")


def main() -> int:
    _set_safe_env()
    _notify_log("scanwr backend ready")
    try:
        import logging

        logging.basicConfig(level=logging.DEBUG)
    except Exception:
        pass

    # Make Scanpy as verbose as possible.
    try:
        import scanpy as sc

        sc.settings.verbosity = 4
        # Force logging to stdout so it shows up in the app console.
        sc.settings.logfile = sys.stdout
        _notify_log(f"scanpy verbosity={sc.settings.verbosity}")
    except Exception as e:
        _notify_log(f"scanpy verbosity setup failed: {e}")

    for line in sys.stdin:
        line = line.strip()
        if not line:
            continue
        try:
            req = json.loads(line)
            req_id = req.get("id")
            method = req.get("method")
            params = req.get("params", {})

            result = _handle(method, params)
            resp = {"id": req_id, "result": result}
        except Exception as e:
            resp = {
                "id": req.get("id") if "req" in locals() else None,
                "error": {"code": 1, "message": str(e)},
            }
        sys.stdout.write(json.dumps(resp) + "\n")
        sys.stdout.flush()

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
