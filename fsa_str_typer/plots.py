"""Figures at final print size (6.5 in wide, 7 pt text, 600 dpi), following Elsevier artwork guidance."""
from __future__ import annotations

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import font_manager
from matplotlib.ticker import FuncFormatter

from .calling import Bins, SampleResult

_fonts = {f.name for f in font_manager.fontManager.ttflist}
FONT = next((f for f in ("Arial", "Helvetica", "Liberation Sans") if f in _fonts), "DejaVu Sans")
plt.rcParams.update({"font.family": FONT, "font.size": 7, "axes.linewidth": 0.5, "xtick.major.width": 0.5,
                     "ytick.major.width": 0.5, "xtick.major.size": 2, "ytick.major.size": 2})
WIDTH_IN = 6.5
DPI = 600
COLORS = ["#1f4fd1", "#1a9a2a", "#000000", "#d62020", "#8a2be2"]     # blue, green, yellow (drawn black), red, 5th dye
OVERVIEW_COLORS = ["#0000ff", "#00a000", "#000000", "#ff0000", "#8a2be2", "#ff9900"]
_k = FuncFormatter(lambda v, p: f"{v / 1000:g}k" if v else "0")


def _save(fig, out, formats):
    for f in formats:
        kw = {"pil_kwargs": {"compression": "tiff_lzw"}} if f == "tif" else {}
        fig.savefig(f"{out}.{f}", dpi=DPI, **kw)
    plt.close(fig)


def overview(res: SampleResult, out: str, title: str | None = None, project: str = "", formats=("png", "tif")):
    """Full-run overlay of all channels against scan number (as in the instrument software)."""
    t = res.trace; n = t.signal.shape[1]; x = np.arange(n)
    fig = plt.figure(figsize=(WIDTH_IN, 2.05))
    ax = fig.add_axes([0.085, 0.06, 0.9, 0.62])
    for c in range(t.n_channels):
        col = OVERVIEW_COLORS[-1] if c == t.n_channels - 1 else OVERVIEW_COLORS[c]
        ax.plot(x, t.signal[c], color=col, lw=0.35)
    ax.set_xlim(0, n); ax.set_ylim(0, t.signal.max() * 1.05)
    ax.xaxis.tick_top(); ax.set_xticks(range(0, n, 1000)); ax.yaxis.set_major_formatter(_k); ax.set_ylabel("RFU", labelpad=2)
    for s in ax.spines.values(): s.set_color("#999999"); s.set_linewidth(0.4)
    ax.grid(True, ls=":", color="#cccccc", lw=0.4)
    L = [f"Project Name : {project}", f"Sample Name : {title or t.sample_name}", f"Sample Type : {t.meta.get('sample_type', '')}"]
    R = [f"Sample File Name : {t.path.name}", f"Size Standard : {t.meta.get('size_standard_name', '')}", "Plot Type : raw (baseline-subtracted)"]
    for i, (l, r) in enumerate(zip(L, R)):
        fig.text(0.085, 0.985 - i * 0.065, l, fontsize=7, fontweight="bold", va="top")
        fig.text(0.58, 0.985 - i * 0.065, r, fontsize=7, fontweight="bold", va="top")
    _save(fig, out, formats)


def rfu_trace(res: SampleResult, bins: Bins, out: str, xlim=(100, 365), formats=("png", "tif")):
    """Allele-labelled RFU traces against fragment size, one panel per sample dye."""
    kit = bins.kit; t = res.trace; dyes = kit.sample_dyes
    bp = res.sizing.sizer(np.arange(t.signal.shape[1]))
    fig, axs = plt.subplots(len(dyes), 1, figsize=(WIDTH_IN, 1.03 * len(dyes)), sharex=True)
    for k, (ax, dye) in enumerate(zip(np.atleast_1d(axs), dyes)):
        c = kit.channel_of(dye, t.dyes); col = COLORS[k % len(COLORS)]
        x = bp; y = t.signal[c]; w = np.isfinite(x) & (x >= xlim[0]) & (x <= xlim[1])
        ymax = max(y[w & (x > xlim[0] + 3)].max(), 100) * 1.45
        ax.plot(x[w], y[w], color=col, lw=0.45)
        calls = [(m, a) for m in kit.markers_in_dye(dye) if m in res.loci for a in res.loci[m].alleles]
        calls.sort(key=lambda z: z[1][1])
        for j, (m, a) in enumerate(calls):
            dx = 0
            if j + 1 < len(calls) and calls[j + 1][1][1] - a[1] < 5.5: dx = -(6 if len(a[0]) > 2 else 4)
            if j > 0 and a[1] - calls[j - 1][1][1] < 5.5: dx = 6 if len(a[0]) > 2 else 4
            ww = w & (np.abs(x - a[1]) < 0.6); jj = np.argmax(np.where(ww, y, -1e9))
            ax.annotate(a[0], (x[jj], y[jj]), xytext=(dx, 1.5), textcoords="offset points", ha="center", va="bottom",
                        fontsize=7, fontweight="bold", bbox=dict(boxstyle="square,pad=0.12", fc="white", ec=col, lw=0.4))
        for m in kit.markers_in_dye(dye):
            if m not in bins.by: continue
            lo, hi = min(bins.by[m].values()), max(bins.by[m].values())
            ax.axvspan(lo - 2, hi + 2, ymin=0.90, ymax=1, color=col, alpha=0.15, lw=0)
            ax.text((lo + hi) / 2, ymax * 0.95, m, ha="center", va="center", fontsize=7, color=col if k != 2 else "#333")
        ax.set_ylim(-0.04 * ymax, ymax); ax.set_xlim(*xlim)
        ax.set_ylabel(f"{dye}\nRFU", color=col if k != 2 else "#333", labelpad=2)
        ax.spines[["top", "right"]].set_visible(False); ax.yaxis.set_major_formatter(_k)
    np.atleast_1d(axs)[-1].set_xlabel("Fragment size (bp)")
    fig.tight_layout(h_pad=0.25, pad=0.3)
    _save(fig, out, formats)
