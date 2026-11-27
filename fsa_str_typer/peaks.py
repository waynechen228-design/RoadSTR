"""Peak detection in the sample dye channels and classification of non-allelic signals."""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.signal import find_peaks, peak_widths

from .io import Trace

DETECT_FLOOR_RFU = 50      # peaks are detected down to this height so that artifacts can be recognised


def detect_peaks(trace: Trace, sizer, channels: list[int], dye_names: list[str], floor: float = DETECT_FLOOR_RFU) -> pd.DataFrame:
    rows = []
    for c, dye in zip(channels, dye_names):
        y = trace.signal[c]
        p, pr = find_peaks(y, height=floor, distance=4, prominence=floor * 0.6)
        if not len(p):
            continue
        fw = peak_widths(y, p, rel_height=0.5)[0]
        for s, h, w in zip(p, pr["peak_heights"], fw):
            if 0 < s < len(y) - 1:                               # parabolic refinement of the apex
                y0, y1, y2 = y[s - 1], y[s], y[s + 1]; den = y0 - 2 * y1 + y2
                off = 0.5 * (y0 - y2) / den if den != 0 else 0.0
            else:
                off = 0.0
            bp = float(sizer(s + off)[0])
            if np.isfinite(bp):
                rows.append(dict(dye=dye, channel=c, scan=int(s), bp=bp, height=float(h), fwhm=float(w)))
    return pd.DataFrame(rows, columns=["dye", "channel", "scan", "bp", "height", "fwhm"])


def bleed_ratios(P: pd.DataFrame, min_source: float = 2000) -> dict:
    """Cross-dye bleed (pull-up) ratio measured in this sample, per (source dye, target dye):
    90th percentile of target height at the same scan / source height, over strong source peaks."""
    out = {}
    for a in P.dye.unique():
        src = P[(P.dye == a) & (P.height >= min_source)]
        for b in P.dye.unique():
            if a == b or not len(src):
                continue
            v = []
            for r in src.itertuples():
                o = P[(P.dye == b) & ((P.scan - r.scan).abs() <= 1)]
                v.append(o.height.max() / r.height if len(o) else 0.0)
            out[(a, b)] = float(np.percentile(v, 90))
    return out


def classify(P: pd.DataFrame, trace: Trace, size_std_channel: int, params: dict) -> pd.DataFrame:
    """Label every peak: '' (candidate allele) or the artifact type.

    saturated                  peak reaches the detector ceiling (height and balance unreliable; still a candidate)
    co-migrates with saturated peak   under an off-scale peak of another dye but far taller than the bleed
                               into this dye under that dye's other off-scale peaks (candidate, flagged)
    pull-up (saturated source) under an off-scale peak of another dye and consistent with its bleed
    under saturated peak (undetermined)  not enough other off-scale peaks to judge (not called, flagged)
    spike                      narrower than any DNA fragment (electrical spike / bubble)
    size-std bleed             under a much taller size-standard peak
    pull-up                    under a taller peak of another dye and no taller than this sample's
                               measured bleed ratio allows; a clearly taller peak is a co-migrating allele
    """
    P = P.copy()
    if not len(P):
        P["artifact"] = []; return P
    sat = trace.saturated_scans(params["saturation_rfu"])
    ratio = bleed_ratios(P)
    ls, lpr = find_peaks(trace.signal[size_std_channel], height=300, distance=20, prominence=200)
    lh = lpr["peak_heights"]
    k = params["pullup_call_ratio"]
    sat_apex = {}
    for r in P.itertuples():
        if len(sat[r.channel]) and np.abs(sat[r.channel] - r.scan).min() <= 3 and r.height > 0.5 * params["saturation_rfu"]:
            sat_apex.setdefault(r.channel, []).append(r.scan)
    labels = []
    for r in P.itertuples():
        own_sat = len(sat[r.channel]) and np.abs(sat[r.channel] - r.scan).min() <= 3
        sat_src = any(len(sat[c]) and np.abs(sat[c] - r.scan).min() <= 3
                      for c in set(P.channel) if c != r.channel)
        taller = P[(P.channel != r.channel) & ((P.scan - r.scan).abs() <= 4) & (P.height > 2.5 * r.height)]
        pull = any(r.height <= k * ratio.get((q.dye, r.dye), 0.0) * q.height + 0.05 * q.height for q in taller.itertuples())
        lz = (np.abs(ls - r.scan) <= 3) & (lh > 2.5 * r.height)
        if own_sat:
            lab = "saturated"
        elif sat_src:
            # Bleed under an off-scale peak is disproportionate, so the measured ratio does not apply.
            # Compare with the bleed into this dye under the OTHER off-scale peaks of the same source dye:
            # a peak far taller than all of them is a real allele co-migrating with the saturated peak.
            src = [c for c in set(P.channel) if c != r.channel and len(sat[c]) and np.abs(sat[c] - r.scan).min() <= 3]
            refs = []
            for c in src:
                for ap in sat_apex.get(c, []):
                    if abs(ap - r.scan) <= 6:
                        continue
                    t = P[(P.channel == r.channel) & ((P.scan - ap).abs() <= 4)]
                    refs.append(t.height.max() if len(t) else 0.0)
            if len(refs) >= 2:
                ref = float(np.percentile(refs, 75))
                lab = "co-migrates with saturated peak" if r.height > 2 * max(ref, 250.0) else "pull-up (saturated source)"
            else:
                lab = "under saturated peak (undetermined)"
        elif r.fwhm < params["spike_max_fwhm_scans"]:
            lab = "spike"
        elif lz.any():
            lab = "size-std bleed"
        elif pull:
            lab = "pull-up"
        else:
            lab = ""
        labels.append(lab)
    P["artifact"] = labels
    return P


NON_ALLELIC = {"pull-up (saturated source)", "spike", "size-std bleed", "pull-up", "shoulder", "stutter",
               "under saturated peak (undetermined)", "pull-up (off-ladder)"}
