"""Kit definition files (YAML)."""
from __future__ import annotations

from dataclasses import dataclass
from importlib import resources
from pathlib import Path

import yaml

DEFAULT_ANALYSIS = {
    "analytical_threshold_rfu": 150,
    "stochastic_threshold_rfu": None,
    "bin_tolerance_bp": 0.5,
    "stutter_max_ratio": 0.15,
    "min_allele_fraction": 0.2,
    "heterozygote_balance_min": 0.5,
    "pullup_call_ratio": 2.0,
    "saturation_rfu": 32000,
    "spike_max_fwhm_scans": 2.0,
    "microvariant_detection": True,
    "match_threshold_percent": 80,
}


def nominal_offset(allele: str, repeat_unit: int | None) -> float:
    """Length of an allele relative to allele '0' in bp: '9.3' -> 9*4+3."""
    if allele in ("X", "Y"):
        return 0.0 if allele == "X" else 6.0
    base, _, frac = allele.partition(".")
    return int(base) * (repeat_unit or 4) + int(frac or 0)


def allele_sort_key(a: str):
    if a in ("X", "Y"):
        return (-1, a)
    base, _, frac = a.partition(".")
    return (int(base), int(frac or 0))


@dataclass
class Kit:
    name: str
    data: dict
    path: Path

    @property
    def markers(self) -> dict:
        return self.data["markers"]

    @property
    def analysis(self) -> dict:
        a = dict(DEFAULT_ANALYSIS); a.update(self.data.get("analysis") or {}); return a

    @property
    def size_standard(self) -> dict:
        return self.data["size_standards"][self.data["size_standard"]]

    @property
    def sample_dyes(self) -> list[str]:
        return [d["name"] for d in self.data["dyes"] if d.get("role") != "size_standard"]

    @property
    def control_genotype(self) -> dict | None:
        c = self.data.get("control_dna") or {}
        return c.get("genotype")

    def markers_in_dye(self, dye: str) -> list[str]:
        return [m for m, v in self.markers.items() if v["dye"] == dye]

    def channel_of(self, dye: str, trace_dyes: list[str]) -> int:
        """Channel index of a kit dye in a trace: by dye name, else by the channel number in the kit file."""
        if dye in trace_dyes:
            return trace_dyes.index(dye)
        for d in self.data["dyes"]:
            if d["name"] == dye:
                return int(d["channel"]) - 1
        raise KeyError(dye)

    def check_file_dyes(self, trace_dyes: list[str]) -> list[str]:
        missing = [d["name"] for d in self.data["dyes"] if d["name"] not in trace_dyes]
        return [f"kit dye(s) {missing} not named in the file (file dyes: {trace_dyes}); using channel numbers"] if missing else []


def load_kit(name_or_path: str) -> Kit:
    """Load a kit by file path or by the name of a bundled kit (e.g. 'cla_identifiler_plus')."""
    p = Path(name_or_path)
    if not p.exists():
        cand = resources.files("fsa_str_typer") / "kits" / f"{name_or_path}.yaml"
        if not cand.is_file():
            raise FileNotFoundError(f"kit '{name_or_path}' not found (bundled kits: {bundled_kits()})")
        p = Path(str(cand))
    data = yaml.safe_load(p.read_text())
    for m, v in data["markers"].items():
        v["ladder_alleles"] = [str(a) for a in v["ladder_alleles"]]
    return Kit(name=data["kit"]["name"], data=data, path=p)


def bundled_kits() -> list[str]:
    d = resources.files("fsa_str_typer") / "kits"
    return sorted(Path(str(f)).stem for f in d.iterdir() if str(f).endswith(".yaml") and "template" not in str(f))
