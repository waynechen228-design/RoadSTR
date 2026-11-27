"""End-to-end analysis of one run."""
from __future__ import annotations

import platform
from dataclasses import dataclass, field
from pathlib import Path

import pandas as pd

from . import __version__
from .calling import Bins, SampleResult, call_sample, control_concordance, negative_signals
from .compare import cross_matrix, pair_verdicts
from .io import read_fsa
from .kit import Kit
from .ladder import LadderResult, build_bins, load_bins


@dataclass
class RunResult:
    kit: Kit
    bins: Bins
    ladder: LadderResult | None
    samples: dict                                  # id -> SampleResult
    positive: SampleResult | None = None
    negative: SampleResult | None = None
    pairs: pd.DataFrame = field(default_factory=pd.DataFrame)
    matrix: pd.DataFrame = field(default_factory=pd.DataFrame)
    positive_table: pd.DataFrame = field(default_factory=pd.DataFrame)
    negative_table: pd.DataFrame = field(default_factory=pd.DataFrame)
    run_notes: list = field(default_factory=list)
    settings: dict = field(default_factory=dict)


def run(kit: Kit, samples: list, ladder: str | None = None, bins_csv: str | None = None,
        positive: str | None = None, negative: str | None = None, pairs: list | None = None) -> RunResult:
    if not ladder and not bins_csv:
        raise ValueError("an allelic ladder (.fsa) or a bins file (.csv) is required")
    notes = []
    lres = None
    if ladder:
        lt = read_fsa(ladder); lres = build_bins(lt, kit)
        if not lres.ok:
            raise RuntimeError("allelic ladder does not fit the kit file:\n" + "\n".join(lres.report))
        if lres.warnings:
            notes.append(f"ladder check: {lres.warnings} warning(s) - see the Ladder sheet")
        table = lres.bins; ladder_run = lt.meta["run_name"]
    else:
        table = load_bins(bins_csv); ladder_run = None
        notes.append(f"bins loaded from {Path(bins_csv).name} (not measured in this run)")
    bins = Bins(table, kit)
    neg = call_sample(read_fsa(negative), kit, bins) if negative else None
    pos = call_sample(read_fsa(positive), kit, bins, negative=neg) if positive else None
    res = {}
    for f in samples:
        t = read_fsa(f)
        if ladder_run and t.meta["run_name"] != ladder_run:
            notes.append(f"{t.id}: run '{t.meta['run_name']}' differs from the ladder run '{ladder_run}'")
        res[t.id] = call_sample(t, kit, bins, negative=neg)
    markers = list(kit.markers)
    profiles = {k: v.genotype() for k, v in res.items()}
    if pos is not None:
        profiles[f"{pos.id} (positive control)"] = pos.genotype()
    out = RunResult(kit, bins, lres, res, pos, neg, run_notes=notes)
    if pairs:
        out.pairs = pair_verdicts(pairs, profiles, markers, kit.analysis["match_threshold_percent"])
    out.matrix = cross_matrix(profiles)
    if pos is not None and kit.control_genotype:
        out.positive_table = control_concordance(pos, kit.control_genotype)
    if neg is not None:
        out.negative_table = negative_signals(neg, kit.analysis["analytical_threshold_rfu"])
    out.settings = {"fsa-str-typer version": __version__, "python": platform.python_version(),
                    "kit file": str(kit.path), "kit": kit.name, "ladder": ladder or "", "bins file": bins_csv or "",
                    **{f"analysis.{k}": v for k, v in kit.analysis.items()}}
    return out
