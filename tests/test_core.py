from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from fsa_str_typer.calling import Bins, call_sample
from fsa_str_typer.cli import main
from fsa_str_typer.compare import pair_verdicts, tanabe
from fsa_str_typer.io import read_fsa
from fsa_str_typer.kit import load_kit
from fsa_str_typer.ladder import build_bins
from fsa_str_typer.sizing import LocalSouthern, match_size_standard

EX = Path(__file__).resolve().parents[1] / "examples"


@pytest.fixture(scope="module")
def kit():
    return load_kit("cla_identifiler_plus")


@pytest.fixture(scope="module")
def ladder(kit):
    return build_bins(read_fsa(EX / "allelic_ladder.fsa"), kit)


def test_local_southern_reproduces_standards():
    s = np.array([1600, 1850, 2250, 2350, 2450, 2900, 4050, 4500, 4600, 5200])
    b = np.array([75, 100, 139, 150, 160, 200, 300, 340, 350, 400])
    f = LocalSouthern(s, b)
    assert np.allclose(f(s), b, atol=1e-6)


def test_size_standard_matched(kit):
    t = read_fsa(EX / "positive_control_9947A.fsa")
    ss = kit.size_standard
    r = match_size_standard(t.signal[-1], ss["fragments"], ss["exclude"])
    assert {75, 100, 139, 150, 160, 200, 300, 340, 350, 400, 450, 490, 500} <= set(r.found)


def test_ladder_fits_kit(kit, ladder):
    assert ladder.ok
    n = sum(len(v["ladder_alleles"]) for v in kit.markers.values())
    assert len(ladder.bins) == n


def test_positive_control_concordant(kit, ladder):
    bins = Bins(ladder.bins, kit)
    res = call_sample(read_fsa(EX / "positive_control_9947A.fsa"), kit, bins)
    exp = kit.control_genotype
    for m, a in exp.items():
        assert sorted(res.loci[m].genotype) == sorted(str(x) for x in a), m


def test_negative_control_has_no_profile(kit, ladder):
    bins = Bins(ladder.bins, kit)
    res = call_sample(read_fsa(EX / "negative_control.fsa"), kit, bins)
    # sporadic low-level signals only: nothing resembling a DNA profile. Off-ladder artifacts may be
    # strong in a negative control, but they are flagged as off-ladder, never reported as normal alleles.
    on_ladder = [a for c in res.loci.values() for a in c.alleles if "off-ladder" not in a[4]]
    assert max((a[3] for a in on_ladder), default=0) < 500
    assert sum(1 for c in res.loci.values() if any("off-ladder" not in a[4] for a in c.alleles)) <= 5


def test_microvariant_virtual_bin(kit, ladder):
    bins = Bins(ladder.bins, kit)
    b22 = bins.by["FGA"]["22"]
    a, dev, on = bins.assign("FGA", b22 + 2.0, 0.5)
    assert a == "22.2" and not on


def test_tanabe_and_swap_detection():
    p = {"A": ["10", "12"], "B": ["8"], "AMEL": ["X"]}
    q = {"A": ["10", "12"], "B": ["8"], "AMEL": ["X", "Y"]}
    r = {"A": ["11", "13"], "B": ["9"], "AMEL": ["X"]}
    assert tanabe(p, q) == 100.0 and tanabe(p, r) == 0.0
    prof = {"pb1": p, "ips1": r, "pb2": r, "ips2": p}
    v = pair_verdicts([("pb1", "ips1"), ("pb2", "ips2")], prof, ["A", "B", "AMEL"])
    assert list(v.Verdict) == ["No match", "No match"] and "label swap" in v.Note[0]


def test_cli_run(tmp_path):
    rc = main(["run", "--kit", "cla_identifiler_plus", "--ladder", str(EX / "allelic_ladder.fsa"),
               "--positive", str(EX / "positive_control_9947A.fsa"), "--negative", str(EX / "negative_control.fsa"),
               "--samples", str(EX / "positive_control_9947A.fsa"), "--out", str(tmp_path)])
    assert rc == 0
    x = pd.read_excel(tmp_path / "str_results.xlsx", sheet_name="Positive_control")
    assert x.Concordant.all()
