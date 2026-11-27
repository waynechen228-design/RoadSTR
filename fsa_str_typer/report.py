"""Excel report."""
from __future__ import annotations

import pandas as pd
from openpyxl import Workbook
from openpyxl.styles import Alignment, Font, PatternFill

from .pipeline import RunResult

FILL = {"Mismatch": "F4CCCC", "No match": "F4CCCC", "Fail": "F4CCCC", "FAIL": "F4CCCC", False: "F4CCCC",
        "Check": "FFF2CC", "REVIEW": "FFF2CC"}


def _sheet(wb, name, df):
    ws = wb.create_sheet(name[:31]); ws.append([str(c) for c in df.columns])
    for r in df.itertuples(index=False):
        ws.append([None if (isinstance(v, float) and pd.isna(v)) else (str(v) if isinstance(v, (list, tuple)) else v) for v in r])
    for c in ws[1]:
        c.font = Font(name="Arial", bold=True, color="FFFFFF"); c.fill = PatternFill("solid", fgColor="305496")
        c.alignment = Alignment(wrap_text=True, horizontal="center")
    for row in ws.iter_rows(min_row=2):
        for c in row:
            c.font = Font(name="Arial", size=10)
            if c.value in FILL: c.fill = PatternFill("solid", fgColor=FILL[c.value])
    for col in ws.columns:
        ws.column_dimensions[col[0].column_letter].width = min(70, max(9, max(len(str(c.value or "")) for c in col) + 2))
    ws.freeze_panes = "A2"


def genotype_table(rr: RunResult) -> pd.DataFrame:
    cols = {sid: [", ".join(s.loci[m].genotype) or "—" for m in rr.kit.markers] for sid, s in rr.samples.items()}
    if rr.positive is not None:
        cols[f"{rr.positive.id} (positive)"] = [", ".join(rr.positive.loci[m].genotype) or "—" for m in rr.kit.markers]
    return pd.DataFrame({"Marker": list(rr.kit.markers), **cols})


def allele_table(rr: RunResult) -> pd.DataFrame:
    rows = []
    allres = dict(rr.samples)
    if rr.positive is not None: allres[rr.positive.id] = rr.positive
    for sid, s in allres.items():
        for m, c in s.loci.items():
            if not c.alleles:
                rows.append(dict(Sample=sid, Marker=m, Allele="", Size_bp=None, Deviation_bp=None, Height_RFU=None,
                                 Flags="", Status=c.status, Notes="; ".join(c.notes)))
            for a, bp, dv, h, fl in c.alleles:
                rows.append(dict(Sample=sid, Marker=m, Allele=a, Size_bp=bp, Deviation_bp=dv, Height_RFU=h, Flags=fl,
                                 Status=c.status, Notes="; ".join(c.notes)))
    return pd.DataFrame(rows)


def qc_table(rr: RunResult) -> pd.DataFrame:
    rows = []
    allres = dict(rr.samples)
    if rr.positive is not None: allres[rr.positive.id] = rr.positive
    if rr.negative is not None: allres[rr.negative.id] = rr.negative
    for sid, s in allres.items():
        t = s.trace; a = [x for c in s.loci.values() for x in c.alleles]
        hb = [min(x[3] for x in c.alleles) / max(x[3] for x in c.alleles) for c in s.loci.values() if len(c.alleles) == 2]
        rows.append(dict(Sample=sid, File=t.path.name, Run=t.meta["run_name"], Well=t.meta["well"], Status=s.status,
                         Size_std_found=None if s.sizing is None else len(s.sizing.found),
                         Size_std_residual_bp=None if s.sizing is None else s.sizing.max_residual_bp,
                         Loci_called=sum(1 for c in s.loci.values() if c.alleles),
                         Loci_check=sum(1 for c in s.loci.values() if c.status == "Check"),
                         Max_abs_dev_bp=round(max(abs(x[2]) for x in a), 2) if a else None,
                         Min_height_RFU=min(x[3] for x in a) if a else None,
                         Min_het_balance=round(min(hb), 2) if hb else None,
                         Saturated_peaks=int((s.peaks.artifact == "saturated").sum()) if len(s.peaks) else 0,
                         Notes="; ".join(s.notes)))
    return pd.DataFrame(rows)


def write_excel(rr: RunResult, path: str):
    wb = Workbook(); wb.remove(wb.active)
    if len(rr.pairs): _sheet(wb, "Pairs", rr.pairs)
    _sheet(wb, "Genotypes", genotype_table(rr))
    _sheet(wb, "Alleles", allele_table(rr))
    _sheet(wb, "QC", qc_table(rr))
    if len(rr.positive_table): _sheet(wb, "Positive_control", rr.positive_table)
    if rr.negative is not None:
        nt = rr.negative_table if len(rr.negative_table) else pd.DataFrame({"result": ["no in-marker signal above the analytical threshold"]})
        _sheet(wb, "Negative_control", nt)
    m = rr.matrix.reset_index().rename(columns={"index": "Tanabe %"}); _sheet(wb, "Cross_match", m)
    peaks = []
    allres = dict(rr.samples)
    for x in (rr.positive, rr.negative):
        if x is not None: allres[x.id] = x
    for sid, s in allres.items():
        if len(s.peaks):
            p = s.peaks.copy(); p.insert(0, "Sample", sid); peaks.append(p)
    if peaks: _sheet(wb, "Peaks", pd.concat(peaks).round(3))
    if rr.ladder is not None:
        _sheet(wb, "Ladder", pd.DataFrame({f"Ladder check: {rr.ladder.status}": rr.ladder.report}))
    _sheet(wb, "Bins", pd.DataFrame([(mk, a, round(b, 3)) for mk, d in rr.bins.by.items() for a, b in d.items()],
                                     columns=["Marker", "Allele", "bp"]))
    _sheet(wb, "Settings", pd.DataFrame({"Setting": list(rr.settings), "Value": [str(v) for v in rr.settings.values()]}
                                         ).pipe(lambda d: pd.concat([d, pd.DataFrame({"Setting": ["run note"] * len(rr.run_notes),
                                                                                       "Value": rr.run_notes})])))
    wb.save(path)
