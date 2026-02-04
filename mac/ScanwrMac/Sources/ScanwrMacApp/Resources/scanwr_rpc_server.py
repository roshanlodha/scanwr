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
    os.environ.setdefault("MPLBACKEND", "Agg")
    os.environ.setdefault("MPLCONFIGDIR", tempfile.mkdtemp(prefix="scanwr-mpl-"))
    os.environ.setdefault("XDG_CACHE_HOME", tempfile.gettempdir())
    os.environ.setdefault("NUMBA_CACHE_DIR", tempfile.mkdtemp(prefix="scanwr-numba-"))
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


def _configure_python_logging(verbosity: int) -> None:
    try:
        import logging

        # Keep non-Scanpy libraries quiet by default. Scanpy verbosity controls Scanpy's chatter.
        base = logging.WARNING if verbosity <= 3 else logging.INFO
        logging.getLogger().setLevel(base)
        for noisy in [
            "numba",
            "numba.core",
            "numba.core.byteflow",
            "h5py",
            "anndata",
            "matplotlib",
            "fontTools",
        ]:
            logging.getLogger(noisy).setLevel(logging.WARNING)
    except Exception:
        pass


def _json_to_py(v: Any) -> Any:
    # Match the app's convention: empty strings mean "None" for backend args.
    if v is None:
        return None
    if isinstance(v, str):
        s = v.strip()
        return None if s == "" else v
    return v


def inspect_h5ad(path: str, var_names_limit: int = 5000) -> Dict[str, Any]:
    import anndata as ad  # local import after env set

    p = Path(path).expanduser().resolve()
    if not p.exists():
        raise FileNotFoundError(str(p))
    if p.suffix.lower() != ".h5ad":
        raise ValueError(f"Expected .h5ad input, got: {p.name}")

    adata = ad.read_h5ad(str(p))

    try:
        obs_cols = list(getattr(adata, "obs", {}).columns)
    except Exception:
        obs_cols = []

    groupby_candidates: List[str] = []
    numeric_candidates: List[str] = []
    try:
        import pandas as pd  # type: ignore

        obs = adata.obs
        for k in obs_cols:
            try:
                s = obs[k]
                try:
                    if pd.api.types.is_numeric_dtype(s) and not pd.api.types.is_bool_dtype(s):
                        numeric_candidates.append(k)
                except Exception:
                    pass
                if str(getattr(s.dtype, "name", "")) == "category":
                    groupby_candidates.append(k)
                    continue
                nunique = int(getattr(s, "nunique")()) if hasattr(s, "nunique") else None
                if nunique is not None and 1 <= nunique <= 50:
                    groupby_candidates.append(k)
            except Exception:
                continue
    except Exception:
        pass

    try:
        layers = list(getattr(adata, "layers", {}).keys())
    except Exception:
        layers = []

    obsm_dims: Dict[str, int] = {}
    try:
        import numpy as np  # type: ignore

        obsm = getattr(adata, "obsm", None)
        if obsm is not None:
            keys = list(getattr(obsm, "keys", lambda: [])())
            for k in keys:
                try:
                    arr = np.asarray(obsm[k])
                    if arr.ndim == 1:
                        obsm_dims[str(k)] = 1
                    elif arr.ndim == 2:
                        obsm_dims[str(k)] = int(arr.shape[1])
                except Exception:
                    continue
    except Exception:
        obsm_dims = {}

    try:
        var_names_all = list(getattr(adata, "var_names", []))
    except Exception:
        var_names_all = []
    total = len(var_names_all)
    limit = max(0, int(var_names_limit or 0))
    truncated = total > limit if limit else False
    var_names = var_names_all[:limit] if limit else var_names_all

    has_raw = False
    try:
        has_raw = getattr(adata, "raw", None) is not None
    except Exception:
        has_raw = False

    try:
        n_obs, n_vars = map(int, getattr(adata, "shape", (0, 0)))
    except Exception:
        n_obs, n_vars = (0, 0)

    return {
        "path": str(p),
        "nObs": n_obs,
        "nVars": n_vars,
        "obsColumns": sorted(map(str, obs_cols)),
        "groupbyCandidates": sorted(set(map(str, groupby_candidates))),
        "numericObsColumns": sorted(set(map(str, numeric_candidates))),
        "varNames": list(map(str, var_names)),
        "varNamesTotal": int(total),
        "varNamesTruncated": bool(truncated),
        "layers": sorted(map(str, layers)),
        "hasRaw": bool(has_raw),
        "obsmDims": {k: int(v) for k, v in sorted(obsm_dims.items())},
    }


def plot_violin(params: Dict[str, Any]) -> Dict[str, Any]:
    import anndata as ad  # local import after env set
    import matplotlib as mpl

    h5ad_path = str(params.get("h5adPath") or params.get("path") or "").strip()
    if not h5ad_path:
        raise ValueError("Missing h5adPath")
    p = Path(h5ad_path).expanduser().resolve()
    if not p.exists():
        raise FileNotFoundError(str(p))

    keys = params.get("keys")
    if not keys:
        raise ValueError("Missing keys")
    if isinstance(keys, str):
        keys = [keys]
    keys = [str(k) for k in keys]

    output_path = str(params.get("outputPath") or "").strip()
    if not output_path:
        raise ValueError("Missing outputPath")
    out = Path(output_path).expanduser().resolve()
    if out.suffix.lower() != ".svg":
        out = out.with_suffix(".svg")
    out.parent.mkdir(parents=True, exist_ok=True)

    # Ensure SVG text stays editable.
    # Source - https://stackoverflow.com/a/50563208
    new_rc_params = {'text.usetex': False,
    "svg.fonttype": 'none'
    }
    mpl.rcParams.update(new_rc_params)

    adata = ad.read_h5ad(str(p))

    kwargs: Dict[str, Any] = {}
    groupby = _json_to_py(params.get("groupby"))
    if groupby is not None:
        kwargs["groupby"] = str(groupby)

    # Map json keys to Scanpy arg names (snake_case).
    if "log" in params:
        kwargs["log"] = bool(params.get("log"))
    if "useRaw" in params or "use_raw" in params:
        use_raw = params.get("useRaw", params.get("use_raw"))
        if use_raw is not None:
            kwargs["use_raw"] = bool(use_raw)
    if "stripplot" in params:
        kwargs["stripplot"] = bool(params.get("stripplot"))
    if "jitter" in params:
        jitter = params.get("jitter")
        if isinstance(jitter, bool):
            kwargs["jitter"] = jitter
        elif jitter is None:
            pass
        else:
            try:
                kwargs["jitter"] = float(jitter)
            except Exception:
                kwargs["jitter"] = True
    if "size" in params:
        try:
            kwargs["size"] = int(params.get("size"))
        except Exception:
            kwargs["size"] = 1
    layer = _json_to_py(params.get("layer"))
    if layer is not None:
        kwargs["layer"] = str(layer)
    density_norm = _json_to_py(params.get("densityNorm", params.get("density_norm")))
    if density_norm is not None:
        kwargs["density_norm"] = density_norm
    order = params.get("order")
    if isinstance(order, list) and order:
        kwargs["order"] = [str(x) for x in order]
    if "multiPanel" in params or "multi_panel" in params:
        kwargs["multi_panel"] = bool(params.get("multiPanel", params.get("multi_panel")))
    if "xlabel" in params:
        kwargs["xlabel"] = str(params.get("xlabel") or "")
    ylabel = _json_to_py(params.get("ylabel"))
    if ylabel is not None:
        if isinstance(ylabel, str) and "," in ylabel:
            parts = [x.strip() for x in ylabel.split(",") if x.strip()]
            kwargs["ylabel"] = parts if parts else ylabel
        else:
            kwargs["ylabel"] = ylabel
    if "rotation" in params and params.get("rotation") is not None:
        kwargs["rotation"] = float(params.get("rotation"))
    if "show" in params and params.get("show") is not None:
        kwargs["show"] = bool(params.get("show"))
    scale = _json_to_py(params.get("scale"))
    if scale is not None:
        kwargs["scale"] = scale

    kwds = params.get("kwds") or {}
    if isinstance(kwds, dict):
        for k, v in kwds.items():
            kwargs[str(k)] = _json_to_py(v)

    # Always run headlessly for the GUI.
    kwargs["show"] = False
    kwargs["save"] = None

    ret = sc.pl.violin(adata, keys=keys, **kwargs)

    import matplotlib.pyplot as plt

    fig = None
    try:
        fig = getattr(ret, "figure", None)
        if fig is None:
            fig = getattr(ret, "fig", None)
    except Exception:
        fig = None
    if fig is None:
        fig = plt.gcf()

    fig.savefig(str(out), format="svg", bbox_inches="tight")
    plt.close(fig)

    return {"svgPath": str(out)}


def _parse_keyref(keyref: str) -> Tuple[str, str]:
    s = str(keyref or "").strip()
    if ":" not in s:
        raise ValueError(f"Invalid key reference (expected obs:<col>, gene:<name>, or obsm:<key>:<dim>): {s}")
    kind, name = s.split(":", 1)
    kind = kind.strip().lower()
    name = name.strip()
    if kind not in ("obs", "gene", "obsm") or not name:
        raise ValueError(f"Invalid key reference: {s}")
    return kind, name


def _vector_from_keyref(adata: Any, keyref: str, layer: str | None, use_raw: bool | None) -> Any:
    import numpy as np

    kind, name = _parse_keyref(keyref)
    if kind == "obs":
        try:
            return adata.obs[name].to_numpy()
        except Exception:
            return adata.obs[name].values

    if kind == "obsm":
        # Expected: "obsm:<key>:<dim>" (dim is 0-based). If dim is omitted, default to 0.
        if ":" in name:
            obsm_key, dim_raw = name.split(":", 1)
            obsm_key = obsm_key.strip()
            dim_raw = dim_raw.strip()
            if not dim_raw:
                dim = 0
            else:
                dim = int(dim_raw)
        else:
            obsm_key = name.strip()
            dim = 0
        if not obsm_key:
            raise ValueError(f"Invalid obsm key reference: {keyref}")
        obsm = getattr(adata, "obsm", None)
        if obsm is None or obsm_key not in obsm:
            raise ValueError(f"obsm key not found: {obsm_key}")
        arr = np.asarray(obsm[obsm_key])
        if arr.ndim == 1:
            if dim != 0:
                raise ValueError(f"obsm key is 1D; dim must be 0: {obsm_key}")
            return arr
        if arr.ndim == 2:
            if dim < 0 or dim >= arr.shape[1]:
                raise ValueError(f"obsm dim out of range: {obsm_key}:{dim} (n={arr.shape[1]})")
            return arr[:, dim]
        raise ValueError(f"Unsupported obsm value for key: {obsm_key}")

    # gene
    view = adata
    if use_raw:
        raw = getattr(adata, "raw", None)
        if raw is not None:
            view = raw

    if layer:
        layers = getattr(adata, "layers", {})
        if layer not in layers:
            raise ValueError(f"Layer not found: {layer}")
        X = adata[:, name].layers[layer]
    else:
        X = view[:, name].X

    try:
        # sparse
        if hasattr(X, "toarray"):
            X = X.toarray()
    except Exception:
        pass

    arr = np.asarray(X)
    if arr.ndim == 2 and arr.shape[1] == 1:
        arr = arr[:, 0]
    return arr


def plot_custom(params: Dict[str, Any]) -> Dict[str, Any]:
    import anndata as ad  # local import after env set
    import matplotlib as mpl
    import matplotlib.pyplot as plt

    # Optional seaborn (preferred)
    sns = None
    try:
        import seaborn as seaborn  # type: ignore

        sns = seaborn
    except Exception:
        sns = None

    h5ad_path = str(params.get("h5adPath") or "").strip()
    if not h5ad_path:
        raise ValueError("Missing h5adPath")
    p = Path(h5ad_path).expanduser().resolve()
    if not p.exists():
        raise FileNotFoundError(str(p))

    plot_type = str(params.get("plotType") or "").strip().lower()
    if plot_type not in ("scatter", "violin", "box", "density"):
        raise ValueError(f"Unsupported plotType: {plot_type}")

    output_path = str(params.get("outputPath") or "").strip()
    if not output_path:
        raise ValueError("Missing outputPath")
    out = Path(output_path).expanduser().resolve()
    if out.suffix.lower() != ".svg":
        out = out.with_suffix(".svg")
    out.parent.mkdir(parents=True, exist_ok=True)

    # Keep SVG text editable
    mpl.rcParams.update({"text.usetex": False, "svg.fonttype": "none"})

    _notify_log(f"Plot: loading {p.name}…")
    adata = ad.read_h5ad(str(p))

    layer = str(params.get("layer") or "").strip() or None
    use_raw = params.get("useRaw")
    if use_raw is not None:
        use_raw = bool(use_raw)

    xref = params.get("x")
    yref = params.get("y")
    cref = params.get("color")

    title = str(params.get("title") or "").strip()
    subtitle = str(params.get("subtitle") or "").strip()
    legend_title = str(params.get("legendTitle") or "").strip()
    xlabel = str(params.get("xLabel") or "").strip()
    ylabel = str(params.get("yLabel") or "").strip()
    xrot = params.get("xTickRotation")
    try:
        xrot_f = float(xrot) if xrot is not None else None
    except Exception:
        xrot_f = None

    point_size = params.get("pointSize")
    alpha = params.get("alpha")
    try:
        point_size_f = float(point_size) if point_size is not None else None
    except Exception:
        point_size_f = None
    try:
        alpha_f = float(alpha) if alpha is not None else None
    except Exception:
        alpha_f = None

    def _default_label(keyref: Any) -> str:
        if not keyref:
            return ""
        kind, name = _parse_keyref(str(keyref))
        if kind == "obs":
            return name
        if kind == "gene":
            return f"{name} (expr)"
        # obsm: "X_umap:0"
        if ":" in name:
            k, d = name.split(":", 1)
            try:
                return f"{k}[{int(d) + 1}]"
            except Exception:
                return k
        return name

    if not xlabel:
        xlabel = _default_label(xref)
    if not ylabel:
        ylabel = _default_label(yref)

    import pandas as pd
    import numpy as np

    df = pd.DataFrame(index=getattr(adata, "obs_names", None))
    if xref:
        df["x"] = _vector_from_keyref(adata, str(xref), layer=layer, use_raw=use_raw)
    if yref:
        df["y"] = _vector_from_keyref(adata, str(yref), layer=layer, use_raw=use_raw)
    if cref:
        df["c"] = _vector_from_keyref(adata, str(cref), layer=layer, use_raw=use_raw)

    # Drop rows with missing values in required columns.
    required_cols: List[str] = []
    if plot_type in ("scatter",):
        required_cols = ["x", "y"]
    if plot_type in ("violin", "box", "density"):
        required_cols = ["y"]
    df = df.dropna(subset=required_cols)

    if plot_type == "scatter":
        if "x" not in df.columns or "y" not in df.columns:
            raise ValueError("Scatter plot requires x and y.")

        fig, ax = plt.subplots(figsize=(7.4, 5.2), dpi=120)
        if "c" in df.columns:
            cvals = df["c"].to_numpy()
            # heuristic: treat as continuous if numeric with many unique values
            is_numeric = np.issubdtype(cvals.dtype, np.number)
            nunique = int(pd.Series(cvals).nunique(dropna=True))
            if is_numeric and nunique > 30:
                sca = ax.scatter(
                    df["x"].to_numpy(),
                    df["y"].to_numpy(),
                    c=cvals,
                    s=point_size_f or 10,
                    alpha=alpha_f if alpha_f is not None else 0.8,
                    cmap="viridis",
                    linewidths=0,
                )
                cb = fig.colorbar(sca, ax=ax, fraction=0.05, pad=0.04)
                cb.set_label(legend_title or _default_label(cref))
            else:
                if sns is None:
                    # matplotlib fallback (categorical)
                    cats = pd.Series(cvals).astype(str)
                    uniq = list(dict.fromkeys(cats.tolist()))
                    palette = plt.get_cmap("tab20")
                    for i, u in enumerate(uniq):
                        mask = cats == u
                        ax.scatter(
                            df.loc[mask, "x"].to_numpy(),
                            df.loc[mask, "y"].to_numpy(),
                            s=point_size_f or 10,
                            alpha=alpha_f if alpha_f is not None else 0.8,
                            label=u,
                            color=palette(i % 20),
                            linewidths=0,
                        )
                    ax.legend(title=legend_title or _default_label(cref), loc="best", frameon=False)
                else:
                    sns.scatterplot(
                        data=df,
                        x="x",
                        y="y",
                        hue="c",
                        ax=ax,
                        s=point_size_f or 30,
                        alpha=alpha_f if alpha_f is not None else 0.8,
                        linewidth=0,
                    )
                    leg = ax.get_legend()
                    if leg is not None:
                        leg.set_title(legend_title or _default_label(cref))
        else:
            ax.scatter(
                df["x"].to_numpy(),
                df["y"].to_numpy(),
                s=point_size_f or 10,
                alpha=alpha_f if alpha_f is not None else 0.8,
                linewidths=0,
            )

        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.grid(True, alpha=0.18)

    if plot_type in ("violin", "box"):
        if "y" not in df.columns:
            raise ValueError(f"{plot_type} plot requires y.")

        fig, ax = plt.subplots(figsize=(7.4, 5.2), dpi=120)
        if sns is None:
            raise RuntimeError("seaborn is required for violin/box plots (not installed in runtime).")

        has_x = "x" in df.columns
        hue = "c" if "c" in df.columns else None
        if not has_x and hue is not None:
            # Seaborn's behavior with only y+hue is a bit inconsistent; keep it simple for now.
            _notify_log("Note: ignoring color because x/category is not set.")
            hue = None

        if plot_type == "violin":
            if has_x:
                sns.violinplot(data=df, x="x", y="y", hue=hue, ax=ax, cut=0, inner="quart")
            else:
                sns.violinplot(data=df, y="y", ax=ax, cut=0, inner="quart")
        else:
            if has_x:
                sns.boxplot(data=df, x="x", y="y", hue=hue, ax=ax)
            else:
                sns.boxplot(data=df, y="y", ax=ax)

        if hue is not None:
            leg = ax.get_legend()
            if leg is not None:
                leg.set_title(legend_title or _default_label(cref))
        else:
            try:
                ax.get_legend().remove()
            except Exception:
                pass

        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)

    if plot_type == "density":
        if "y" not in df.columns:
            raise ValueError("Density plot requires a value (y).")
        fig, ax = plt.subplots(figsize=(7.4, 5.2), dpi=120)
        if sns is None:
            raise RuntimeError("seaborn is required for density plots (not installed in runtime).")

        fill = bool(params.get("densityFill")) if params.get("densityFill") is not None else False

        if "x" in df.columns:
            sns.kdeplot(data=df, x="y", hue="x", ax=ax, common_norm=False, fill=fill)
            leg = ax.get_legend()
            if leg is not None:
                leg.set_title(legend_title or _default_label(xref))
        else:
            sns.kdeplot(data=df, x="y", ax=ax, fill=fill)

        ax.set_xlabel(xlabel or _default_label(yref))
        ax.set_ylabel(ylabel or "Density")
        ax.grid(True, alpha=0.18)

    if title:
        fig.suptitle(title, y=0.985, fontsize=14, fontweight="bold")
    if subtitle:
        fig.text(0.5, 0.945 if title else 0.98, subtitle, ha="center", va="top", fontsize=10, color="#666666")

    if xrot_f is not None:
        plt.setp(ax.get_xticklabels(), rotation=xrot_f, ha="right")

    fig.tight_layout()
    _notify_log(f"Saving SVG: {out.name}")
    fig.savefig(str(out), format="svg", bbox_inches="tight")
    plt.close(fig)
    return {"svgPath": str(out)}

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


def _run_module(
    adata,
    spec_id: str,
    params: Dict[str, Any],
    *,
    plots_dir: Path | None = None,
    sample: str | None = None,
) -> None:
    import scanpy as sc  # local import after env set

    # Backwards compatibility: map RAPIDS ids (and older shorthand) back to Scanpy ids.
    legacy_to_scanpy = {
        "pp.calculate_qc_metrics": "scanpy.pp.calculate_qc_metrics",
        "rapids_singlecell.pp.filter_cells": "scanpy.pp.filter_cells",
        "rapids_singlecell.pp.filter_genes": "scanpy.pp.filter_genes",
        "rapids_singlecell.pp.calculate_qc_metrics": "scanpy.pp.calculate_qc_metrics",
        "rapids_singlecell.pp.scrublet": "scanpy.pp.scrublet",
        "rapids_singlecell.pp.highly_variable_genes": "scanpy.pp.highly_variable_genes",
        "rapids_singlecell.pp.normalize_total": "scanpy.pp.normalize_total",
        "rapids_singlecell.pp.log1p": "scanpy.pp.log1p",
        "rapids_singlecell.pp.pca": "scanpy.tl.pca",
        "rapids_singlecell.pp.neighbors": "scanpy.pp.neighbors",
        "rapids_singlecell.tl.umap": "scanpy.tl.umap",
        "rapids_singlecell.tl.leiden": "scanpy.tl.leiden",
        "rapids_singlecell.tl.rank_genes_groups": "scanpy.tl.rank_genes_groups",
    }
    spec_id = legacy_to_scanpy.get(spec_id, spec_id)

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

    def _opt_float(key: str) -> Optional[float]:
        v = (params or {}).get(key)
        if v is None:
            return None
        if isinstance(v, (int, float)):
            return float(v)
        if isinstance(v, str):
            if not v.strip():
                return None
            return float(v)
        return float(v)

    def _opt_str(key: str) -> Optional[str]:
        v = (params or {}).get(key)
        if v is None:
            return None
        if isinstance(v, str):
            s = v.strip()
            return s if s else None
        return str(v).strip() or None

    def _opt_bool(key: str, default: bool) -> bool:
        v = (params or {}).get(key)
        if v is None:
            return default
        if isinstance(v, bool):
            return v
        if isinstance(v, (int, float)):
            return bool(v)
        if isinstance(v, str):
            s = v.strip().lower()
            if s in ("true", "1", "yes", "y", "t", "on"):
                return True
            if s in ("false", "0", "no", "n", "f", "off"):
                return False
        return bool(v)

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

    if spec_id == "scanpy.pp.scrublet":
        batch_key = _opt_str("batch_key")
        sc.pp.scrublet(adata, batch_key=batch_key)
        return

    if spec_id == "scanpy.pp.highly_variable_genes":
        n_top_genes = _opt_int("n_top_genes")
        batch_key = _opt_str("batch_key")
        kwargs: Dict[str, Any] = {"adata": adata}
        if n_top_genes is not None:
            kwargs["n_top_genes"] = n_top_genes
        if batch_key is not None:
            kwargs["batch_key"] = batch_key
        sc.pp.highly_variable_genes(**kwargs)
        return

    if spec_id == "scanpy.pp.calculate_qc_metrics":
        # Prefer the simplified, checkbox-driven interface (mt/ribo/hb). Fall back to
        # advanced parameters if those aren't provided.
        uses_gene_flags = any(k in (params or {}) for k in ("use_mt", "use_ribo", "use_hb"))
        log1p = _opt_bool("log1p", True)

        percent_top_list = _parse_int_list(str((params or {}).get("percent_top", "") or ""))
        percent_top = tuple(percent_top_list) if percent_top_list else None

        kwargs: Dict[str, Any] = {
            "adata": adata,
            "inplace": True,
            "log1p": log1p,
        }

        if uses_gene_flags:
            qc_vars: List[str] = []
            if _opt_bool("use_mt", True):
                adata.var["mt"] = adata.var_names.str.startswith("MT-")
                qc_vars.append("mt")
            if _opt_bool("use_ribo", True):
                adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
                qc_vars.append("ribo")
            if _opt_bool("use_hb", True):
                adata.var["hb"] = adata.var_names.str.contains(r"^HB[^(P)]", regex=True, na=False)
                qc_vars.append("hb")
            kwargs["qc_vars"] = qc_vars
        else:
            # Advanced mode (legacy): pass through Scanpy arguments.
            kwargs["expr_type"] = str((params or {}).get("expr_type", "counts") or "counts")
            kwargs["var_type"] = str((params or {}).get("var_type", "genes") or "genes")
            kwargs["qc_vars"] = _parse_csv_list(str((params or {}).get("qc_vars", "") or "")) or ()
            layer = _opt_str("layer")
            if layer is not None:
                kwargs["layer"] = layer
            kwargs["use_raw"] = _opt_bool("use_raw", False)
            parallel = _opt_int("parallel")
            if parallel is not None:
                kwargs["parallel"] = parallel

        if percent_top is not None:
            kwargs["percent_top"] = percent_top

        sc.pp.calculate_qc_metrics(**kwargs)
        return

    if spec_id == "scanpy.pp.normalize_total":
        target_sum = _opt_float("target_sum")
        if target_sum is None:
            sc.pp.normalize_total(adata)
        else:
            sc.pp.normalize_total(adata, target_sum=target_sum)
        return

    if spec_id == "scanpy.pp.log1p":
        sc.pp.log1p(adata)
        return

    if spec_id == "scanpy.tl.pca":
        sc.tl.pca(adata)
        create_scatterplot = _opt_bool("create_scatterplot", True)
        if create_scatterplot and plots_dir is not None:
            import matplotlib.pyplot as plt

            sample_dir = _sample_plots_dir(plots_dir, sample)
            if sample_dir is not None:
                out = sample_dir / "pca_scatter.png"
                sc.pl.pca(adata, show=False)
                plt.gcf().savefig(str(out), format="png", dpi=160, bbox_inches="tight")
                plt.close("all")
        return

    if spec_id == "scanpy.pp.neighbors":
        sc.pp.neighbors(adata)
        return

    if spec_id == "scanpy.tl.umap":
        sc.tl.umap(adata)
        create_scatterplot = _opt_bool("create_scatterplot", True)
        if create_scatterplot and plots_dir is not None:
            import matplotlib.pyplot as plt

            sample_dir = _sample_plots_dir(plots_dir, sample)
            if sample_dir is not None:
                out = sample_dir / "umap_scatter.png"
                sc.pl.umap(adata, show=False)
                plt.gcf().savefig(str(out), format="png", dpi=160, bbox_inches="tight")
                plt.close("all")
        return

    if spec_id == "scanpy.tl.leiden":
        try:
            import leidenalg  # noqa: F401
            import igraph  # noqa: F401
        except Exception as e:
            raise RuntimeError(
                "Leiden clustering requires the optional deps `leidenalg` + `igraph`, "
                "but they are missing from the app's embedded Python runtime. "
                "Rebuild the embedded runtime (mac/ScanwrMac/scripts/build_python_runtime.sh) "
                "and repackage the .app."
            ) from e
        # UI uses "res"; allow "resolution" too for advanced usage.
        res = _opt_float("res")
        if res is None:
            res = _opt_float("resolution")
        if res is None:
            sc.tl.leiden(adata)
        else:
            sc.tl.leiden(adata, resolution=res)
        return

    if spec_id == "scanpy.tl.rank_genes_groups":
        import matplotlib.pyplot as plt

        groupby_raw = _opt_str("groupby") or "leiden"
        groupby_norm = groupby_raw.strip().lower()
        groupby = {
            "leiden": "leiden",
            "louvain": "louvain",
            "kmeans": "kmeans",
            "k-means": "kmeans",
            "k_means": "kmeans",
        }.get(groupby_norm, groupby_raw.strip())

        try:
            obs = getattr(adata, "obs", None)
            if obs is None or groupby not in obs.columns:
                raise KeyError(groupby)
        except Exception:
            raise RuntimeError(
                f"rank_genes_groups: groupby='{groupby}' not found in adata.obs. "
                "Run clustering first (e.g. Leiden) or choose a different groupby."
            )

        sc.tl.rank_genes_groups(adata, groupby=groupby)

        create_dotplot = _opt_bool("create_dotplot", False)
        create_heatmap = _opt_bool("create_heatmap", True)

        if (create_dotplot or create_heatmap) and plots_dir is not None:
            sample_dir = _sample_plots_dir(plots_dir, sample)
            if sample_dir is None:
                return
            group_slug = _sanitize_filename(groupby)

            if create_dotplot:
                out = sample_dir / f"rank_genes_groups__{group_slug}__dotplot.png"
                sc.pl.rank_genes_groups_dotplot(adata, groupby=groupby, show=False)
                plt.gcf().savefig(str(out), format="png", dpi=160, bbox_inches="tight")
                plt.close("all")

            if create_heatmap:
                out = sample_dir / f"rank_genes_groups__{group_slug}__heatmap.png"
                sc.pl.rank_genes_groups_heatmap(adata, groupby=groupby, show=False)
                plt.gcf().savefig(str(out), format="png", dpi=160, bbox_inches="tight")
                plt.close("all")
        return

    raise NotImplementedError(spec_id)


def list_modules() -> List[Dict[str, Any]]:
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
        {
            "id": "scanpy.pp.calculate_qc_metrics",
            "group": "pp",
            "namespace": "core",
            "title": "Calculate QC Metrics",
            "scanpyQualname": "scanpy.pp.calculate_qc_metrics",
        },
        {
            "id": "scanpy.pp.scrublet",
            "group": "pp",
            "namespace": "core",
            "title": "Scrublet (Doublet Detection)",
            "scanpyQualname": "scanpy.pp.scrublet",
        },
        {
            "id": "scanpy.pp.highly_variable_genes",
            "group": "pp",
            "namespace": "core",
            "title": "Highly Variable Genes",
            "scanpyQualname": "scanpy.pp.highly_variable_genes",
        },
        {
            "id": "scanpy.pp.normalize_total",
            "group": "pp",
            "namespace": "core",
            "title": "Normalize Total Counts",
            "scanpyQualname": "scanpy.pp.normalize_total",
        },
        {
            "id": "scanpy.pp.log1p",
            "group": "pp",
            "namespace": "core",
            "title": "Log1p",
            "scanpyQualname": "scanpy.pp.log1p",
        },
        {
            "id": "scanpy.tl.pca",
            "group": "tl",
            "namespace": "core",
            "title": "PCA",
            "scanpyQualname": "scanpy.tl.pca",
        },
        {
            "id": "scanpy.pp.neighbors",
            "group": "pp",
            "namespace": "core",
            "title": "Neighbors",
            "scanpyQualname": "scanpy.pp.neighbors",
        },
        {
            "id": "scanpy.tl.umap",
            "group": "tl",
            "namespace": "core",
            "title": "UMAP",
            "scanpyQualname": "scanpy.tl.umap",
        },
        {
            "id": "scanpy.tl.leiden",
            "group": "tl",
            "namespace": "core",
            "title": "Leiden",
            "scanpyQualname": "scanpy.tl.leiden",
        },
        {
            "id": "scanpy.tl.rank_genes_groups",
            "group": "tl",
            "namespace": "core",
            "title": "Rank Genes Groups",
            "scanpyQualname": "scanpy.tl.rank_genes_groups",
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


def _sample_plots_dir(plots_root: Path | None, sample: str | None) -> Path | None:
    if plots_root is None:
        return None
    sample_slug = _sanitize_filename((sample or "").strip() or "sample")
    out = plots_root / sample_slug
    out.mkdir(parents=True, exist_ok=True)
    return out


def run_pipeline(input_path: str, output_dir: str, steps: List[Dict[str, Any]]) -> Dict[str, Any]:
    import anndata as ad  # local import after env set

    in_path = Path(input_path).expanduser().resolve()
    out_dir = Path(output_dir).expanduser().resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    plots_root = out_dir / "plots"
    plots_root.mkdir(parents=True, exist_ok=True)

    if not in_path.exists():
        raise FileNotFoundError(str(in_path))
    if in_path.suffix.lower() != ".h5ad":
        raise ValueError(f"Expected .h5ad input, got: {in_path.name}")

    _notify_log("Loading input .h5ad…")
    adata = ad.read_h5ad(str(in_path))

    checkpoints: List[str] = []
    for i, step in enumerate(steps, start=1):
        spec_id = str(step.get("specId"))
        params = step.get("params") or {}
        _notify_log(f"Step {i}/{len(steps)}: {spec_id}")
        _run_module(adata, spec_id, params, plots_dir=plots_root, sample="single")

        ck_name = f"{i:02d}_{_sanitize_filename(spec_id)}.h5ad"
        ck_path = out_dir / ck_name
        _notify_log(f"Saving checkpoint: {ck_path.name}")
        adata.write_h5ad(str(ck_path))
        checkpoints.append(str(ck_path))

        # Reload from checkpoint to ensure the next step operates on an on-disk checkpointed state.
        adata = ad.read_h5ad(str(ck_path))

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
    import anndata as ad  # local import after env set

    out_dir = Path(output_dir).expanduser().resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    project_name = str(project_name or "").strip()
    if project_name:
        project = out_dir / _sanitize_filename(project_name)
        project.mkdir(parents=True, exist_ok=True)
    else:
        project = out_dir
    scanwr_dir = project / ".scanwr"
    scanwr_dir.mkdir(parents=True, exist_ok=True)
    checkpoints_dir = scanwr_dir / "checkpoints"
    history_dir = scanwr_dir / "history"
    plots_root = project / "plots"
    checkpoints_dir.mkdir(parents=True, exist_ok=True)
    history_dir.mkdir(parents=True, exist_ok=True)
    plots_root.mkdir(parents=True, exist_ok=True)

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
            adata = ad.read_h5ad(str(checkpoint_h5ad))
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
            _run_module(adata, spec_id, params, plots_dir=plots_root, sample=sample)

            _notify_log(f"[{sample}] Updating checkpoint.h5ad")
            adata.write_h5ad(str(checkpoint_h5ad))
            # Update cached signature list in order.
            done = requested_sigs[:i]
            _write_cached_prefix(done)
            checkpoints.append(str(checkpoint_h5ad))
            adata = ad.read_h5ad(str(checkpoint_h5ad))

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
    if method == "set_verbosity":
        try:
            import logging

            level = int(params.get("level", 3))
            level = max(0, min(4, level))
            _configure_python_logging(level)
            # Keep root logger conservative; backend functions should emit explicit logs to the app.
            base = logging.WARNING if level <= 3 else logging.INFO
            logging.getLogger().setLevel(base)
            _notify_log(f"verbosity={level}")
            return {"ok": True}
        except Exception as e:
            return {"ok": False, "error": str(e)}
    if method == "detect_reader":
        return detect_reader(path=params["path"])
    if method == "inspect_h5ad":
        limit = int(params.get("varNamesLimit", 5000))
        return inspect_h5ad(path=params["path"], var_names_limit=limit)
    if method == "plot_violin":
        return plot_violin(params=params)
    if method == "plot_custom":
        return plot_custom(params=params)
    if method == "run_pipeline":
        # New API: per-sample raw input paths
        if "samples" in params:
            return run_pipeline_multi(
                output_dir=params["outputDir"],
                project_name=params.get("projectName", ""),
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

        # Default to a conservative log level; the app console is the primary signal channel.
        logging.basicConfig(level=logging.WARNING)
    except Exception:
        pass

    try:
        v_raw = os.environ.get("SCANWR_VERBOSITY", "3")
        try:
            v = int(v_raw)
        except Exception:
            v = 3
        v = max(0, min(4, v))
        _configure_python_logging(v)
    except Exception:
        pass

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
