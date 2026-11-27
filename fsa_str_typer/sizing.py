"""Size-standard detection and Local Southern sizing (scan number -> base pairs)."""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy.signal import find_peaks


class SizingError(RuntimeError):
    """Raised when the size standard cannot be matched reliably. The sample must not be called."""


def _southern3(xs, ys):
    """Fit  L = L0 + K / (m - M0)  exactly through three points.
    Rearranged as  m*L = M0*L + L0*m - C  (C = L0*M0 - K), which is linear in (M0, L0, C)."""
    x = np.asarray(xs, float); y = np.asarray(ys, float)
    A = np.column_stack([y, x, -np.ones(3)])
    try:
        M0, L0, C = np.linalg.solve(A, x * y)
    except np.linalg.LinAlgError:
        p = np.polyfit(x, y, 2); return lambda m: np.polyval(p, m)
    K = L0 * M0 - C
    return lambda m: L0 + K / (np.asarray(m, float) - M0)


class LocalSouthern:
    """Local Southern sizing: for a scan between standards i and i+1, average the
    3-point Southern fits through (i-1, i, i+1) and (i, i+1, i+2). Undefined outside the standards."""

    def __init__(self, scans, sizes):
        o = np.argsort(scans)
        self.s = np.asarray(scans, float)[o]; self.b = np.asarray(sizes, float)[o]
        if len(self.s) < 4:
            raise SizingError("fewer than 4 size-standard peaks")
        n = len(self.s)
        self.fits = [_southern3(self.s[i:i + 3], self.b[i:i + 3]) for i in range(n - 2)]

    def __call__(self, m):
        m = np.atleast_1d(np.asarray(m, float)); out = np.full(m.shape, np.nan)
        n = len(self.s)
        for i in range(n - 1):
            sel = (m >= self.s[i]) & ((m < self.s[i + 1]) if i < n - 2 else (m <= self.s[i + 1]))
            if not sel.any():
                continue
            fs = []
            if i >= 1: fs.append(self.fits[i - 1])
            if i + 2 < n: fs.append(self.fits[i])
            out[sel] = np.mean([f(m[sel]) for f in fs], axis=0)
        return out


@dataclass
class SizingResult:
    sizer: LocalSouthern
    pairs: list                          # (scan, bp) used for sizing
    found: list                          # fragment sizes found (incl. excluded ones)
    missing: list
    extra_peaks: int                     # size-standard-channel peaks not assigned (bleed, primer front)
    max_residual_bp: float               # standards vs. smooth global fit (mobility consistency)
    status: str = "OK"
    notes: list = field(default_factory=list)


def match_size_standard(liz: np.ndarray, fragments, exclude=(), min_fraction: float = 0.8) -> SizingResult:
    """Assign size-standard peaks to fragment sizes.

    Walks back from the largest fragment and predicts every next position from the local
    scans-per-bp trend, so extra peaks are skipped instead of shifting the assignment.
    Extra peaks come from bleed-through of over-loaded samples into the size-standard
    channel and from the primer front; counting peaks from either end fails on them.
    """
    p, pr = find_peaks(liz, height=100, distance=20, prominence=80)
    h = pr["peak_heights"]
    if len(p) < 4:
        raise SizingError("size-standard channel: fewer than 4 peaks")
    frags = sorted(fragments, reverse=True)
    # typical size-standard peak height, measured in the second half of the run, away from the
    # primer front (which dominates the start of the size-standard channel in sample injections)
    late = h[p > 0.5 * len(liz)]
    ref = late if len(late) >= 3 else h
    med = np.median(np.sort(ref)[-max(3, len(frags) // 2):])
    ok = (h > 0.3 * med) & (h < 3 * med)            # size-standard peaks have similar heights
    # fragments hidden in the primer front (off-scale peaks at the start of the run) are not assigned:
    # small peaks there are mostly front debris and would distort sizing of the first interval
    huge = np.sort(p[(h >= 3 * med) & (p < 0.5 * len(liz))])
    front_end = 0
    if len(huge):                                   # contiguous cluster of off-scale peaks = primer front
        front_end = int(huge[0])
        for x in huge[1:]:
            if x - front_end > 150:
                break                               # a later isolated off-scale peak is bleed, not the front
            front_end = int(x)
        front_end += 20
    ok &= p > front_end
    cand = p[ok]
    best = None
    for i0 in range(len(cand) - 1, max(len(cand) - 5, 0), -1):        # last few peaks as the largest fragment
        for i1 in range(i0 - 1, max(i0 - 5, -1), -1):
            pairs = [(cand[i0], frags[0]), (cand[i1], frags[1])]
            for f in frags[2:]:
                (s1, b1), (s2, b2) = pairs[-1], pairs[-2]
                rate = (s2 - s1) / (b2 - b1)                            # scans per bp, locally
                if len(pairs) >= 3:
                    s3, b3 = pairs[-3]
                    r2 = (s1 - s3) / (b1 - b3)
                    if abs(rate - r2) < 0.3 * rate:
                        rate = 2 * rate - r2                            # follow the trend
                x = s1 - (b1 - f) * rate
                j = int(np.argmin(np.abs(cand - x)))
                if abs(cand[j] - x) <= max(8.0, 0.25 * rate * min(b1 - f, 50)) and cand[j] < s1:
                    pairs.append((int(cand[j]), f))
            rate_ok = all(np.diff([s for s, _ in sorted(pairs)]) > 0)
            score = (len(pairs) if rate_ok else 0, -abs(pairs[0][0] - cand[-1]))
            if best is None or score > best[0]:
                best = (score, pairs)
    pairs = sorted(best[1])
    found = [b for _, b in pairs]
    missing = [f for f in sorted(fragments) if f not in found]
    used = [(s, b) for s, b in pairs if b not in exclude]
    # mobility consistency: residual of the standards against a smooth (cubic) fit
    s_, b_ = np.array(used, float).T
    resid = float(np.abs(b_ - np.polyval(np.polyfit(s_, b_, 3), s_)).max()) if len(used) >= 6 else np.nan
    res = SizingResult(LocalSouthern(s_, b_), used, found, missing, int(len(cand) - len(pairs)), round(resid, 2))
    n_needed = int(np.ceil(min_fraction * len(fragments)))
    if len(found) < n_needed:
        raise SizingError(f"only {len(found)}/{len(fragments)} size-standard fragments matched "
                          f"(missing {missing}); sizing is not reliable")
    if missing:
        res.notes.append(f"size-standard fragments not found: {missing}")
    if res.extra_peaks:
        res.notes.append(f"{res.extra_peaks} extra peak(s) in the size-standard channel skipped "
                         f"(bleed-through from an over-loaded sample or primer front)")
    if np.isfinite(resid) and resid > 2.0:
        res.status = "Check"; res.notes.append(f"irregular size-standard mobility (max residual {resid:.2f} bp)")
    return res
