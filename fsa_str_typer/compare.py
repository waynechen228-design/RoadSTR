"""Profile comparison: locus-by-locus match and percent match (Tanabe), cross-sample matrix."""
from __future__ import annotations

import pandas as pd

SEX_MARKERS = {"AMEL", "Amelogenin", "DYS391", "Yindel"}


def tanabe(p: dict, q: dict, markers=None) -> float | None:
    """Tanabe percent match: 2 x shared alleles / (alleles in p + alleles in q) x 100, STR loci only,
    over loci typed in both profiles (ANSI/ATCC ASN-0002)."""
    ms = [m for m in (markers or p) if m not in SEX_MARKERS and p.get(m) and q.get(m)]
    tot = sum(len(set(p[m])) + len(set(q[m])) for m in ms)
    if not tot:
        return None
    sh = sum(len(set(p[m]) & set(q[m])) for m in ms)
    return 200.0 * sh / tot


def compare_pair(ref: dict, test: dict, markers) -> pd.DataFrame:
    rows = []
    for m in markers:
        a, b = ref.get(m, []), test.get(m, [])
        st = "Not typed" if not a or not b else ("Match" if sorted(a) == sorted(b) else "Mismatch")
        rows.append(dict(Marker=m, Reference=", ".join(a) or "—", Test=", ".join(b) or "—", Result=st))
    return pd.DataFrame(rows)


def cross_matrix(profiles: dict) -> pd.DataFrame:
    ids = list(profiles)
    return pd.DataFrame([[None if tanabe(profiles[a], profiles[b]) is None else round(tanabe(profiles[a], profiles[b]), 1)
                          for b in ids] for a in ids], index=ids, columns=ids)


def pair_verdicts(pairs, profiles: dict, markers, threshold: float = 80) -> pd.DataFrame:
    """For each (reference, test) pair: loci matched, Tanabe %, verdict; if a pair fails, look for a
    better partner among all profiles (detects sample label swaps)."""
    rows = []
    for ref, test in pairs:
        cmp = compare_pair(profiles[ref], profiles[test], markers)
        n_typed = int((cmp.Result != "Not typed").sum()); n_match = int((cmp.Result == "Match").sum())
        t = tanabe(profiles[ref], profiles[test], markers)
        verdict = "Match" if t is not None and t >= threshold else "No match"
        note = ""
        if verdict != "Match":
            alt = []
            for who, other in ((test, ref), (ref, test)):
                for k, prof in profiles.items():
                    if k in (ref, test):
                        continue
                    tt = tanabe(profiles[who], prof, markers)
                    if tt is not None and tt >= threshold:
                        alt.append(f"{who} matches {k} ({tt:.0f}%)")
            if alt:
                note = "; ".join(dict.fromkeys(alt)) + " - possible sample label swap"
        rows.append(dict(Reference=ref, Test=test, Loci_typed=n_typed, Loci_matched=n_match,
                         Tanabe_percent=None if t is None else round(t, 1), Verdict=verdict, Note=note))
    return pd.DataFrame(rows)
