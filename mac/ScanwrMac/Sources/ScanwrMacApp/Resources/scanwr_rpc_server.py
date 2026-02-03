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
        qc_vars = _parse_csv_list(str(params.get("qc_vars", "") or "")) or None
        percent_top = _parse_int_list(str(params.get("percent_top", "") or ""))
        log1p = bool(params.get("log1p", True))
        sc.pp.calculate_qc_metrics(
            adata,
            qc_vars=qc_vars,
            percent_top=percent_top,
            log1p=log1p,
            inplace=True,
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


def run_pipeline(samples: List[Dict[str, Any]], pipeline: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    out: List[Dict[str, Any]] = []
    for s in samples:
        sample = str(s.get("sample", ""))
        group = str(s.get("group", ""))
        path = str(s.get("path", ""))
        _notify_log(f"Reading {sample}…")
        reader, adata = _read_any(path)
        for step in pipeline:
            spec_id = str(step.get("specId"))
            params = step.get("params") or {}
            _notify_log(f"Running {spec_id} on {sample}…")
            _run_module(adata, spec_id, params)
        shape = list(getattr(adata, "shape", (0, 0)))
        out.append(
            {
                "sample": sample,
                "group": group,
                "path": path,
                "reader": reader,
                "shape": shape,
            }
        )
    return out


def _handle(method: str, params: Any) -> Any:
    if method == "ping":
        return {"ok": True}
    if method == "list_modules":
        return list_modules()
    if method == "run_pipeline":
        return run_pipeline(samples=params["samples"], pipeline=params["pipeline"])
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

