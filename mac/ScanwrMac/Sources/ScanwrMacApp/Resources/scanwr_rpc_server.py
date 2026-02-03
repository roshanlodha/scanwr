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

    if spec_id == "pp.calculate_qc_metrics":
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
    return [
        {
            "id": "pp.calculate_qc_metrics",
            "group": "pp",
            "title": "Calculate QC Metrics",
            "scanpyQualname": "scanpy.pp.calculate_qc_metrics",
        }
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


def run_pipeline_multi(output_dir: str, samples: List[Dict[str, Any]], steps: List[Dict[str, Any]]) -> Dict[str, Any]:
    import scanpy as sc  # local import after env set

    out_dir = Path(output_dir).expanduser().resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    results: List[Dict[str, Any]] = []
    for s in samples:
        sample = str(s.get("sample", "")).strip()
        group = str(s.get("group", "")).strip()
        path = str(s.get("path", "")).strip()
        reader_override = str(s.get("reader", "")).strip()

        if not sample or not group or not path:
            raise ValueError("Each sample must include sample, group, path")

        sample_out = out_dir / _sanitize_filename(sample)
        sample_out.mkdir(parents=True, exist_ok=True)

        _notify_log(f"Reading {sample}…")
        reader, adata = _read_with_reader(reader_override, path)

        checkpoints: List[str] = []
        for i, step in enumerate(steps, start=1):
            spec_id = str(step.get("specId"))
            params = step.get("params") or {}
            _notify_log(f"[{sample}] Step {i}/{len(steps)}: {spec_id}")
            _run_module(adata, spec_id, params)

            ck_name = f"{i:02d}_{_sanitize_filename(spec_id)}.h5ad"
            ck_path = sample_out / ck_name
            _notify_log(f"[{sample}] Saving checkpoint: {ck_path.name}")
            adata.write_h5ad(str(ck_path))
            checkpoints.append(str(ck_path))
            adata = sc.read_h5ad(str(ck_path))

        final_path = sample_out / "final.h5ad"
        _notify_log(f"[{sample}] Saving final: {final_path.name}")
        adata.write_h5ad(str(final_path))

        shape = list(getattr(adata, "shape", (0, 0)))
        results.append(
            {
                "sample": sample,
                "group": group,
                "path": path,
                "reader": reader,
                "outputDir": str(sample_out),
                "checkpoints": checkpoints,
                "finalPath": str(final_path),
                "shape": shape,
            }
        )

    return {"outputDir": str(out_dir), "results": results}


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
                samples=params["samples"],
                steps=params["steps"],
            )
        # Legacy (kept for compatibility during iteration)
        return run_pipeline(input_path=params["inputPath"], output_dir=params["outputDir"], steps=params["steps"])
    raise ValueError(f"Unknown method: {method}")


def main() -> int:
    _set_safe_env()
    _notify_log("scanwr backend ready")

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
