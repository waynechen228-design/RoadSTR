"""Reading ABIF (.fsa / .hid) capillary-electrophoresis files."""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
from Bio import SeqIO
from scipy.ndimage import percentile_filter, uniform_filter1d

# Data tags of the dye channels in ABIF files: 4 sample dyes + 1 or 2 more (5- or 6-dye chemistry).
DATA_TAGS = ("DATA1", "DATA2", "DATA3", "DATA4", "DATA105", "DATA106")


def _txt(v):
    return v.decode(errors="replace") if isinstance(v, (bytes, bytearray)) else v


@dataclass
class Trace:
    """One capillary injection."""
    path: Path
    sample_name: str
    dyes: list[str]                 # dye names in channel order (last = size standard)
    raw: np.ndarray                 # channels x scans, as stored (instrument-corrected, no baseline)
    signal: np.ndarray              # channels x scans, baseline-subtracted
    meta: dict = field(default_factory=dict)

    @property
    def id(self) -> str:
        return self.path.stem

    @property
    def n_channels(self) -> int:
        return self.raw.shape[0]

    def saturated_scans(self, level: float) -> list[np.ndarray]:
        """Scans at or above the detector ceiling, per channel."""
        return [np.where(ch >= level)[0] for ch in self.raw]


def baseline(y: np.ndarray, window: int = 501, smooth: int = 301, pct: float = 10) -> np.ndarray:
    """Rolling low-percentile baseline (robust to peaks), lightly smoothed."""
    return uniform_filter1d(percentile_filter(y, pct, size=window), smooth)


def read_fsa(path: str | Path) -> Trace:
    path = Path(path)
    raw_tags = SeqIO.read(str(path), "abi").annotations["abif_raw"]
    keys = [k for k in DATA_TAGS if k in raw_tags]
    if len(keys) < 5:
        raise ValueError(f"{path.name}: expected at least 5 dye channels, found {keys}")
    raw = np.array([np.asarray(raw_tags[k], dtype=float) for k in keys])
    sig = np.array([ch - baseline(ch) for ch in raw])
    dyes = [_txt(raw_tags.get(f"DyeN{i}", f"dye{i}".encode())) for i in range(1, len(keys) + 1)]
    meta = {
        "sample_name": _txt(raw_tags.get("SpNm1", b"")),
        "sample_type": _txt(raw_tags.get("STYP1", b"")),
        "well": _txt(raw_tags.get("TUBE1", b"")),
        "capillary": raw_tags.get("LANE1"),
        "run_name": _txt(raw_tags.get("RunN1", b"")),
        "instrument": _txt(raw_tags.get("MODL1", b"")),
        "dye_set": _txt(raw_tags.get("DySN1", b"")),
        "size_standard_name": _txt(raw_tags.get("StdF1", b"")),
        "n_scans": raw.shape[1],
    }
    return Trace(path=path, sample_name=meta["sample_name"] or path.stem, dyes=dyes, raw=raw, signal=sig, meta=meta)
