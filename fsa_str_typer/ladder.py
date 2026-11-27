"""Allele bins measured from the allelic ladder of the same run, and kit-file validation."""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import pandas as pd
from scipy.signal import find_peaks, peak_widths

from .io import Trace
from .kit import Kit, nominal_offset
from .sizing import match_size_standard


@dataclass
class LadderResult:
    bins: pd.DataFrame                    # Marker, Allele, bp
    ok: bool
    report: list = field(default_factory=list)   # one line per dye / marker
    warnings: int = 0

    @property
    def status(self) -> str:
        if not self.ok:
            return "FAIL"
        return "REVIEW" if self.warnings else "PASS"


def _all_peaks(trace: Trace):
    out = {}
    for c in range(trace.n_channels):
        y = trace.signal[c]
        q, qr = find_peaks(y, height=100, distance=4, prominence=60)
        fw = peak_widths(y, q, rel_height=0.5)[0] if len(q) else np.array([])
        keep = fw <= 30                                     # drop broad humps (dye blobs, bleed)
        out[c] = (q[keep], qr["peak_heights"][keep])
    return out


def build_bins(trace: Trace, kit: Kit, tol: float = 0.8) -> LadderResult:
    """Detect ladder peaks per dye and fit each marker's allele pattern (from the allele names
    and repeat unit) to them; markers of one dye are placed in size order without overlap.
    Every listed allele must be found. Strong unexplained peaks next to a marker are reported,
    because they usually mean an allele is missing from the kit file."""
    ss = kit.size_standard
    sizing = match_size_standard(trace.signal[-1], ss["fragments"], ss.get("exclude", []))
    sz = sizing.sizer
    P = _all_peaks(trace)
    ok = True; nwarn = 0; lines = []; BINS = []
    lines += kit.check_file_dyes(trace.dyes)
    for dye in kit.sample_dyes:
        ci = kit.channel_of(dye, trace.dyes); q, h = P[ci]
        flag = np.array([any(((abs(P[o][0] - s) <= 1) & (P[o][1] > 2 * hh)).any() for o in P if o != ci)
                         for s, hh in zip(q, h)], bool)     # >2x taller peak in another dye = likely pull-up
        bp = sz(q.astype(float)); hh = h.astype(float)
        m_ = np.isfinite(bp); bp, hh, flag, qq = bp[m_], hh[m_], flag[m_], q[m_]
        mk = kit.markers_in_dye(dye)
        n = sum(len(kit.markers[m]["ladder_alleles"]) for m in mk)
        if not len(hh):
            lines.append(f"{dye}: no peaks"); ok = False; continue
        med = np.median(np.sort(hh[~flag])[-n:]) if (~flag).any() else np.median(hh)
        cand = hh >= 0.25 * med
        cb, ch, cf, cs = bp[cand], hh[cand], flag[cand], qq[cand]

        def pick(x, t):
            w = np.where(np.abs(cb - x) <= t)[0]
            if not len(w):
                return -1
            good = [i for i in w if not cf[i]]
            return max(good or list(w), key=lambda i: ch[i])

        opts = []
        for m in mk:
            v = kit.markers[m]; al = v["ladder_alleles"]; ru = v.get("repeat_unit")
            nom = np.array([nominal_offset(a, ru) for a in al], float)
            o = []
            for j in range(len(cb)):
                idx = [j]; pos = [cb[j]]; slope = 1.0
                for kk in range(1, len(nom)):
                    mk_ = [(nom[t], cb[i]) for t, i in enumerate(idx) if i >= 0]
                    if len(mk_) >= 2:
                        xs_, ys_ = zip(*mk_[-6:]); slope = np.polyfit(xs_, ys_, 1)[0]
                    ln, lb = mk_[-1]; gap = nom[kk] - ln; x = lb + gap * slope; pos.append(x)
                    idx.append(pick(x, min(tol + 0.03 * gap, 2.0)))
                hit = [(nom[t], cb[i]) for t, i in enumerate(idx) if i >= 0]
                if len(hit) >= 3:
                    xs_, ys_ = map(np.array, zip(*hit))
                    resid = float(np.abs(ys_ - np.polyval(np.polyfit(xs_, ys_, 2), xs_)).max())
                else:
                    resid = 0.0
                nflag = sum(bool(cf[i]) for i in idx if i >= 0)
                found = [cb[i] for i in idx if i >= 0]
                o.append(dict(idx=idx, pos=np.array(pos), hits=len(hit), resid=resid, lo=min(found), hi=max(found),
                              score=len(hit) * 1000 - nflag * 10 - resid))
            opts.append(o)
        dp = [{i: (o["score"], [i]) for i, o in enumerate(opts[0])}]
        for t in range(1, len(opts)):
            cur = {}
            for i, o in enumerate(opts[t]):
                prev = [(sc, pth) for pi, (sc, pth) in dp[-1].items() if opts[t - 1][pi]["hi"] < o["lo"]]
                if prev:
                    sc, pth = max(prev, key=lambda x: x[0]); cur[i] = (sc + o["score"], pth + [i])
            dp.append(cur)
        if not dp[-1]:
            lines.append(f"{dye}: markers could not be placed in size order"); ok = False; continue
        path = max(dp[-1].values(), key=lambda x: x[0])[1]
        used = np.zeros(len(cb), bool)
        for o_list, choice in zip(opts, path):
            for i in o_list[choice]["idx"]:
                if i >= 0: used[i] = True
        shadow = np.array([any(((abs(P[o][0] - s_) <= 1) & (P[o][1] > 1.2 * x)).any() for o in P if o != ci)
                           for s_, x in zip(cs, ch)], bool)
        lines.append(dye)
        for m, o_list, choice in zip(mk, opts, path):
            v = kit.markers[m]; al = v["ladder_alleles"]; ru = v.get("repeat_unit") or 5.4; o = o_list[choice]
            idx, pos = o["idx"], o["pos"]
            miss = [a for a, i in zip(al, idx) if i < 0]
            for a, i in zip(al, idx):
                if i >= 0:
                    BINS.append((m, a, float(cb[i])))
            span = (pos.min() - ru, pos.max() + 2.5 * ru)
            extra = [(round(float(b), 1), int(x)) for b, x, u, f in zip(cb, ch, used, shadow)
                     if not u and not f and span[0] < b < span[1] and x >= 0.6 * med]
            good = not miss and o["resid"] <= tol
            ok &= good; nwarn += bool(extra)
            lines.append(f"   {m:9s}{o['hits']:3d}/{len(al):<3d} alleles {pos.min():6.1f}-{pos.max():6.1f} bp  "
                         f"mobility residual {o['resid']:.2f} bp  {'ok' if good else 'FAIL'}"
                         + (f"  missing: {miss}" if miss else "")
                         + (f"  WARNING strong unassigned peaks (allele missing from kit file?): {extra}" if extra else ""))
    bins = pd.DataFrame(BINS, columns=["Marker", "Allele", "bp"])
    return LadderResult(bins=bins, ok=bool(ok), report=lines, warnings=nwarn)


def load_bins(path) -> pd.DataFrame:
    b = pd.read_csv(path, dtype={"Allele": str})
    return b[["Marker", "Allele", "bp"]]
