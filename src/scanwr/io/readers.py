from __future__ import annotations

import os
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Optional


@dataclass(frozen=True)
class ReadResult:
    reader_name: str
    adata: object


def _is_10x_mtx_dir(path: Path) -> bool:
    required_any = {"matrix.mtx", "matrix.mtx.gz"}
    barcodes_any = {"barcodes.tsv", "barcodes.tsv.gz"}
    features_any = {"features.tsv", "features.tsv.gz", "genes.tsv", "genes.tsv.gz"}
    names = {p.name for p in path.iterdir() if p.is_file()}
    return bool(names & required_any) and bool(names & barcodes_any) and bool(names & features_any)


def _find_first(path: Path, candidates: list[str]) -> Optional[Path]:
    for name in candidates:
        p = path / name
        if p.exists():
            return p
    return None


def detect_scanpy_reader(path: str | os.PathLike[str]) -> tuple[str, Callable[[], object]]:
    p = Path(path).expanduser().resolve()
    if not p.exists():
        raise FileNotFoundError(str(p))

    def _sc():  # type: ignore[no-redef]
        import scanpy as sc  # local import so env vars can be set first

        return sc

    if p.is_file():
        suffix = "".join(p.suffixes).lower()
        if suffix.endswith(".h5ad"):
            return "scanpy.read_h5ad", lambda: _sc().read_h5ad(str(p))
        if suffix.endswith(".loom"):
            return "scanpy.read_loom", lambda: _sc().read_loom(str(p))
        if suffix.endswith(".mtx") or suffix.endswith(".mtx.gz"):
            return "scanpy.read_mtx", lambda: _sc().read_mtx(str(p))
        # Generic fallback (scanpy.read can handle many formats)
        return "scanpy.read", lambda: _sc().read(str(p))

    # Directory cases
    tenx_h5 = _find_first(
        p,
        [
            "filtered_feature_bc_matrix.h5",
            "raw_feature_bc_matrix.h5",
            "filtered_gene_bc_matrices_h5.h5",
        ],
    )
    if tenx_h5 is not None:
        return "scanpy.read_10x_h5", lambda: _sc().read_10x_h5(str(tenx_h5))

    if _is_10x_mtx_dir(p):
        return "scanpy.read_10x_mtx", lambda: _sc().read_10x_mtx(str(p), var_names="gene_symbols")

    # If it's a folder with a single .h5ad inside, use it.
    h5ads = [x for x in p.iterdir() if x.is_file() and x.name.lower().endswith(".h5ad")]
    if len(h5ads) == 1:
        return "scanpy.read_h5ad", lambda: _sc().read_h5ad(str(h5ads[0]))

    raise ValueError(f"Unsupported input path; cannot detect reader for: {p}")


def read_any(path: str | os.PathLike[str]) -> ReadResult:
    reader_name, fn = detect_scanpy_reader(path)
    adata = fn()
    return ReadResult(reader_name=reader_name, adata=adata)
