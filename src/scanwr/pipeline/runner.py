from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Iterable, List

from scanwr.io.readers import ReadResult, read_any
from scanwr.pipeline.modules import ModuleInstance


@dataclass(frozen=True)
class SampleRow:
    sample: str
    group: str
    path: str


@dataclass(frozen=True)
class RunResult:
    sample: SampleRow
    read: ReadResult
    adata: Any


def run_pipeline(samples: Iterable[SampleRow], pipeline: List[ModuleInstance]) -> List[RunResult]:
    results: List[RunResult] = []
    for s in samples:
        read = read_any(s.path)
        adata = read.adata
        for step in pipeline:
            step.run(adata)
        results.append(RunResult(sample=s, read=read, adata=adata))
    return results

