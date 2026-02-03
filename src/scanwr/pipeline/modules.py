from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Dict, List, Literal, Optional

ModuleGroup = Literal["pp", "tl", "pl"]


@dataclass(frozen=True)
class ModuleSpec:
    id: str
    group: ModuleGroup
    title: str
    scanpy_qualname: str
    color_hex: str


@dataclass
class ModuleInstance:
    spec: ModuleSpec
    params: Dict[str, Any] = field(default_factory=dict)

    def run(self, adata: Any) -> None:
        if self.spec.id == "pp.calculate_qc_metrics":
            _run_calculate_qc_metrics(adata, self.params)
            return
        raise NotImplementedError(self.spec.id)


def _parse_csv_list(value: str) -> List[str]:
    return [x.strip() for x in value.split(",") if x.strip()]


def _parse_int_list(value: str) -> Optional[List[int]]:
    items = _parse_csv_list(value)
    if not items:
        return None
    return [int(x) for x in items]


def _run_calculate_qc_metrics(adata: Any, params: Dict[str, Any]) -> None:
    import scanpy as sc  # local import so env vars can be set first

    qc_vars_raw = str(params.get("qc_vars", "") or "")
    percent_top_raw = str(params.get("percent_top", "") or "")
    log1p = bool(params.get("log1p", True))

    qc_vars = _parse_csv_list(qc_vars_raw) or None
    percent_top = _parse_int_list(percent_top_raw)

    sc.pp.calculate_qc_metrics(
        adata,
        qc_vars=qc_vars,
        percent_top=percent_top,
        log1p=log1p,
        inplace=True,
    )


def available_modules() -> List[ModuleSpec]:
    return [
        ModuleSpec(
            id="pp.calculate_qc_metrics",
            group="pp",
            title="Calculate QC Metrics",
            scanpy_qualname="scanpy.pp.calculate_qc_metrics",
            color_hex="#2D6CDF",
        ),
    ]
