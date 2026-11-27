"""Allele calling for one sample."""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import pandas as pd

from .io import Trace
from .kit import Kit, allele_sort_key
from .peaks import NON_ALLELIC, classify, detect_peaks
from .sizing import SizingError, match_size_standard


@dataclass
class LocusCall:
    marker: str
    alleles: list                       # [(allele, bp, deviation_bp, height_rfu, flags)]
    status: str = "OK"                  # OK | Check | Fail
    notes: list = field(default_factory=list)

    @property
    def genotype(self) -> list[str]:
        return [a[0] for a in self.alleles]


@dataclass
class SampleResult:
    trace: Trace
    loci: dict                          # marker -> LocusCall
    peaks: pd.DataFrame
    sizing: object = None
    status: str = "OK"
    notes: list = field(default_factory=list)

    @property
    def id(self):
        return self.trace.id

    def genotype(self) -> dict:
        return {m: c.genotype for m, c in self.loci.items()}


class Bins:
    """Ladder bins plus virtual positions for off-ladder alleles (microvariants, alleles beyond the ladder)."""

    def __init__(self, bins: pd.DataFrame, kit: Kit):
        self.kit = kit
        self.by = {m: dict(zip(g.Allele, g.bp)) for m, g in bins.groupby("Marker")}
        self.windows = {}
        for dye in kit.sample_dyes:
            ms = [m for m in kit.markers_in_dye(dye) if m in self.by]
            spans = [(min(self.by[m].values()), max(self.by[m].values())) for m in ms]
            for i, (m, (lo, hi)) in enumerate(zip(ms, spans)):
                ru = kit.markers[m].get("repeat_unit") or 4
                wl, wh = lo - ru, hi + ru
                if i > 0: wl = max(wl, (spans[i - 1][1] + lo) / 2)
                if i < len(ms) - 1: wh = min(wh, (hi + spans[i + 1][0]) / 2)
                self.windows[m] = (wl, wh)

    def integer_position(self, m: str, n: int) -> float:
        """Expected bp of integer allele n, interpolated/extrapolated from the integer ladder alleles."""
        pts = sorted((int(a), b) for a, b in self.by[m].items() if a.isdigit())
        x = np.array([p[0] for p in pts], float); y = np.array([p[1] for p in pts], float)
        if x[0] <= n <= x[-1]:
            return float(np.interp(n, x, y))
        sl = np.polyfit(x[-3:], y[-3:], 1)[0] if n > x[-1] else np.polyfit(x[:3], y[:3], 1)[0]
        return float(y[-1] + (n - x[-1]) * sl) if n > x[-1] else float(y[0] - (x[0] - n) * sl)

    def assign(self, m: str, bp: float, tol: float, microvariants: bool = True):
        """Return (allele, deviation, on_ladder) or (None, nearest_dev, False)."""
        b = self.by[m]
        a = min(b, key=lambda k: abs(b[k] - bp)); dev = bp - b[a]
        if abs(dev) <= tol:
            return a, dev, True
        ru = self.kit.markers[m].get("repeat_unit")
        if not (microvariants and ru):
            return None, dev, False
        ints = [int(k) for k in b if k.isdigit()]
        best = None
        for n in range(min(ints) - 3, max(ints) + 4):
            base = self.integer_position(m, n)
            for x in range(0, ru):
                pos = base + x * (self.integer_position(m, n + 1) - base) / ru
                d = bp - pos
                if abs(d) <= tol and (best is None or abs(d) < abs(best[1])):
                    best = (f"{n}" if x == 0 else f"{n}.{x}", d)
        if best and best[0] not in b:
            return best[0], best[1], False
        return None, dev, False


def call_sample(trace: Trace, kit: Kit, bins: Bins, negative: SampleResult | None = None) -> SampleResult:
    prm = kit.analysis; tol = prm["bin_tolerance_bp"]; thr = prm["analytical_threshold_rfu"]
    ss = kit.size_standard
    try:
        sizing = match_size_standard(trace.signal[-1], ss["fragments"], ss.get("exclude", []))
    except SizingError as e:
        loci = {m: LocusCall(m, [], "Fail", ["sizing failed"]) for m in kit.markers}
        return SampleResult(trace, loci, pd.DataFrame(), None, "Fail", [str(e)])
    chans = [kit.channel_of(d, trace.dyes) for d in kit.sample_dyes]
    P = detect_peaks(trace, sizing.sizer, chans, kit.sample_dyes)
    P = classify(P, trace, trace.n_channels - 1, prm)
    P["marker"] = ""; P["allele"] = ""; P["deviation"] = np.nan; P["called"] = False
    neg_peaks = negative.peaks if negative is not None and len(negative.peaks) else None
    P["in_negative"] = False
    if neg_peaks is not None:
        # Signals of the negative control recur at the same scan in other capillaries of the run
        # (fluorescent debris, array artifacts); they are scan-locked, so they are matched by scan.
        # Such peaks are never used to remove a main allele - only to keep them out of the
        # "additional peak" / mixture warnings.
        npk = neg_peaks[(neg_peaks.height >= 100) & (neg_peaks.artifact != "spike")]
        for i, r in P[P.artifact == ""].iterrows():
            n = npk[(npk.dye == r.dye) & ((npk.scan - r.scan).abs() <= 3)]
            if len(n) and r.height < 5 * n.height.max():
                P.loc[i, "in_negative"] = True
    loci = {}
    for m, v in kit.markers.items():
        if m not in bins.by:
            loci[m] = LocusCall(m, [], "Fail", ["marker missing from ladder bins"]); continue
        lo, hi = bins.windows[m]; ru = v.get("repeat_unit")
        sel = P[(P.dye == v["dye"]) & (P.bp >= lo) & (P.bp <= hi)]
        notes = []; cands = []
        for i, r in sel.iterrows():
            a, dev, on = bins.assign(m, r.bp, tol, prm["microvariant_detection"])
            P.loc[i, ["marker", "deviation"]] = [m, dev]
            if a is not None:
                P.loc[i, "allele"] = a
            if (r.artifact == "under saturated peak (undetermined)" and a is not None
                    and r.height >= max(thr, 0.2 * sel.height.max()) and not any(f"allele {a} (" in n for n in notes)):
                notes.append(f"possible allele {a} ({int(r.height)} RFU) hidden under an off-scale peak of another dye; "
                             "re-inject a diluted sample to decide")
            if r.artifact in NON_ALLELIC or r.height < thr:
                continue
            if a is None:
                if r.height >= 0.2 * max(sel.height.max(), thr):
                    notes.append(f"unassigned off-ladder peak {r.bp:.2f} bp ({int(r.height)} RFU)")
                continue
            shadow = bool(len(P[(P.dye != r.dye) & ((P.scan - r.scan).abs() <= 4) & (P.height >= 1.5 * r.height)]))
            if not on and r.in_negative:
                P.loc[i, "artifact"] = "run artifact (off-ladder, in negative control)"; continue
            if not on and shadow:
                # off-ladder alleles are rare and cross-dye bleed is common: an off-ladder peak that sits
                # under a taller peak of another dye is treated as bleed
                P.loc[i, "artifact"] = "pull-up (off-ladder)"; continue
            if not on:
                # an off-ladder name is only given to a stand-alone peak: within 1.6 bp of a clearly taller
                # peak it is a shoulder (incomplete +A addition, split or over-loaded peak), not a microvariant
                near = sel[(sel.index != i) & ((sel.bp - r.bp).abs() <= 1.6) & (sel.height >= 1.5 * r.height)]
                if len(near):
                    P.loc[i, "artifact"] = "shoulder"; continue
            cands.append(dict(i=i, allele=a, bp=r.bp, dev=dev, h=r.height, on=on, sat=r.artifact == "saturated",
                              comig=r.artifact == "co-migrates with saturated peak", neg=bool(r.in_negative), shadow=shadow))
        # one peak per allele (the tallest)
        best = {}
        for c in cands:
            if c["allele"] not in best or c["h"] > best[c["allele"]]["h"]:
                best[c["allele"]] = c
        cands = sorted(best.values(), key=lambda c: -c["h"])
        # stutter: one repeat shorter than a taller candidate and below the stutter ratio
        if ru:
            keep = []
            for c in cands:
                parent = [q for q in cands if q["h"] > c["h"] and abs((q["bp"] - c["bp"]) - ru) < 0.8]
                if parent and c["h"] < prm["stutter_max_ratio"] * parent[0]["h"]:
                    P.loc[c["i"], "artifact"] = "stutter"; continue
                keep.append(c)
            cands = keep
        if not cands:
            loci[m] = LocusCall(m, [], "Fail", notes + ["no allele above the analytical threshold"]); continue
        hidden = any(n.startswith("possible allele") for n in notes)
        top = cands[0]["h"]
        major = [c for c in cands if c["h"] >= prm["min_allele_fraction"] * top]
        minor = [c for c in cands if c["h"] < prm["min_allele_fraction"] * top]
        calls = major[:2]
        status = "Check" if hidden else "OK"
        runart = [c for c in major[2:] + minor if c["neg"]]
        if runart:
            notes.append("ignored run artifact(s) also present in the negative control: "
                         + ", ".join(f"{c['allele']} ({int(c['h'])} RFU)" for c in runart))
        bleed = [c for c in major[2:] + minor if c["shadow"] and not c["neg"]]
        if bleed:
            notes.append("ignored peak(s) under taller peaks of another dye (likely bleed): "
                         + ", ".join(f"{c['allele']} ({int(c['h'])} RFU)" for c in bleed))
        major = major[:2] + [c for c in major[2:] if not (c["neg"] or c["shadow"])]
        minor = [c for c in minor if not (c["neg"] or c["shadow"])]
        if len(major) > 2:
            status = "Check"
            notes.append("additional peak(s) >= {:.0%} of the main allele: ".format(prm["min_allele_fraction"])
                         + ", ".join(f"{c['allele']} ({int(c['h'])} RFU)" for c in major[2:])
                         + " - possible mixture, trisomy or artifact; inspect the trace")
        minor = [c for c in minor if c["h"] >= 0.1 * top]
        if minor:
            notes.append("minor peak(s) not called: " + ", ".join(f"{c['allele']} ({int(c['h'])} RFU)" for c in minor[:5])
                         + (f" and {len(minor) - 5} more" if len(minor) > 5 else ""))
        if len(calls) == 2:
            hb = min(c["h"] for c in calls) / max(c["h"] for c in calls)
            if hb < prm["heterozygote_balance_min"] and not any(c["sat"] for c in calls):
                status = "Check"; notes.append(f"heterozygote imbalance {hb:.2f}")
        st = prm.get("stochastic_threshold_rfu")
        if len(calls) == 1 and st and calls[0]["h"] < st:
            status = "Check"; notes.append(f"homozygote below the stochastic threshold ({st} RFU)")
        for c in calls:
            if not c["on"]:
                status = "Check"; notes.append(f"off-ladder allele {c['allele']} (virtual bin, {c['dev']:+.2f} bp)")
            if c["sat"]:
                notes.append(f"allele {c['allele']} saturated (height/balance not measurable)")
            if c.get("comig"):
                status = "Check"; notes.append(f"allele {c['allele']} co-migrates with an off-scale peak of another dye "
                                               "(kept: far taller than the bleed in the other channels); confirm on the trace")
            if neg_peaks is not None:
                n = neg_peaks[(neg_peaks.dye == v["dye"]) & ((neg_peaks.bp - c["bp"]).abs() <= tol) & (neg_peaks.height >= thr)
                              & (~neg_peaks.artifact.isin(NON_ALLELIC))]
                if len(n) and c["h"] < 3 * n.height.max():
                    status = "Check"; notes.append(f"allele {c['allele']} < 3x a signal at the same position in the negative control")
            P.loc[c["i"], "called"] = True
        calls = sorted(calls, key=lambda c: allele_sort_key(c["allele"]))
        flags = lambda c: ",".join(f for f, on in (("off-ladder", not c["on"]), ("saturated", c["sat"])) if on)
        loci[m] = LocusCall(m, [(c["allele"], round(float(c["bp"]), 2), round(float(c["dev"]), 2), int(c["h"]), flags(c)) for c in calls],
                            status, notes)
    res = SampleResult(trace, loci, P, sizing)
    res.notes = list(sizing.notes)
    nsat = int((P.artifact == "saturated").sum())
    if nsat:
        res.notes.append(f"{nsat} saturated peak(s): sample over-loaded; genotype usable, peak heights and balance are not - "
                         "consider re-injection of a diluted sample")
    if any(c.status == "Fail" for c in loci.values()):
        res.status = "Check"
    if sizing.status != "OK":
        res.status = "Check"
    return res


def control_concordance(res: SampleResult, expected: dict) -> pd.DataFrame:
    rows = []
    for m, exp in expected.items():
        got = res.loci[m].genotype if m in res.loci else []
        exp = sorted([str(a) for a in exp], key=allele_sort_key)
        rows.append(dict(Marker=m, Expected=", ".join(exp), Observed=", ".join(got), Concordant=got == exp))
    return pd.DataFrame(rows)


def negative_signals(res: SampleResult, thr: float) -> pd.DataFrame:
    P = res.peaks
    if not len(P):
        return P
    x = P[(P.height >= thr) & (P.marker != "") & (~P.artifact.isin(NON_ALLELIC))]
    return x[["dye", "marker", "allele", "bp", "height", "deviation", "artifact"]].round(2)
