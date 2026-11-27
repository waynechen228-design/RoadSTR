"""Command-line interface:  str-typer {inspect, validate-kit, run}"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

from . import __version__


def _set(kit, assignments):
    """Override analysis settings:  --set analytical_threshold_rfu=100"""
    a = kit.data.setdefault("analysis", {})
    for s in assignments or []:
        k, _, v = s.partition("=")
        try:
            v = float(v) if "." in v or v.isdigit() else ({"true": True, "false": False, "none": None}.get(v.lower(), v))
        except ValueError:
            pass
        a[k] = v


def cmd_inspect(args):
    from .io import read_fsa
    from .sizing import SizingError, match_size_standard
    for f in args.files:
        t = read_fsa(f)
        print(f"{t.path.name}\n  sample: {t.sample_name}   type: {t.meta['sample_type']}   well: {t.meta['well']}   capillary: {t.meta['capillary']}")
        print(f"  run: {t.meta['run_name']}   instrument: {t.meta['instrument']}   dye set: {t.meta['dye_set']}")
        print(f"  channels: {t.dyes}   size standard: {t.meta['size_standard_name']}   scans: {t.meta['n_scans']}")
        sat = [len(s) for s in t.saturated_scans(32000)]
        print(f"  saturated scans per channel: {sat}")
        if args.kit:
            from .kit import load_kit
            k = load_kit(args.kit); ss = k.size_standard
            try:
                r = match_size_standard(t.signal[-1], ss["fragments"], ss.get("exclude", []))
                print(f"  size standard: {len(r.found)}/{len(ss['fragments'])} fragments matched, residual {r.max_residual_bp} bp"
                      + (f"; {'; '.join(r.notes)}" if r.notes else ""))
            except SizingError as e:
                print(f"  size standard: FAILED - {e}")


def cmd_validate(args):
    from .io import read_fsa
    from .kit import load_kit
    from .ladder import build_bins
    kit = load_kit(args.kit); t = read_fsa(args.ladder)
    print(f"Kit: {kit.name}   Ladder: {t.path.name} ({t.meta['instrument']}, {t.meta['run_name']})\n")
    r = build_bins(t, kit)
    print("\n".join(r.report)); print(f"\nRESULT: {r.status}")
    if args.bins_out:
        r.bins.to_csv(args.bins_out, index=False); print(f"bins written to {args.bins_out}")
    return 0 if r.ok else 1


def _pairs(specs):
    out = []
    for s in specs or []:
        a, sep, b = s.partition(":")
        if not sep:
            sys.exit(f"--pair must look like REFERENCE:TEST (file names without .fsa), got '{s}'")
        out.append((Path(a).stem, Path(b).stem))
    return out


def cmd_run(args):
    from .kit import load_kit
    from .pipeline import run
    from .report import write_excel
    kit = load_kit(args.kit); _set(kit, args.set)
    out = Path(args.out); out.mkdir(parents=True, exist_ok=True)
    rr = run(kit, args.samples, ladder=args.ladder, bins_csv=args.bins, positive=args.positive,
             negative=args.negative, pairs=_pairs(args.pair))
    write_excel(rr, out / "str_results.xlsx")
    print(f"Kit: {kit.name}   ladder: {rr.ladder.status if rr.ladder else 'bins file'}")
    for n in rr.run_notes: print(f"  note: {n}")
    if len(rr.positive_table):
        c = int(rr.positive_table.Concordant.sum()); print(f"Positive control: {c}/{len(rr.positive_table)} loci concordant")
    if rr.negative is not None:
        nt = rr.negative_table
        tol = kit.analysis["bin_tolerance_bp"]
        onl = nt[(nt.allele.astype(str) != "") & (nt.deviation.abs() <= tol)] if len(nt) else nt
        at_allele = len(onl); loci = sorted(set(onl.marker)) if len(onl) else []
        print(f"Negative control: {len(nt)} signal(s) above the analytical threshold inside marker ranges; "
              f"{at_allele} at allele or off-ladder (virtual) allele positions{(' (' + ', '.join(loci) + ')') if loci else ''} - see the Negative_control sheet")
    for sid, s in rr.samples.items():
        chk = [m for m, c in s.loci.items() if c.status != "OK"]
        print(f"  {sid:20s} {s.status:5s} loci needing review: {chk if chk else 'none'}")
    if len(rr.pairs):
        print("\nPairs (reference : test)")
        for r in rr.pairs.itertuples():
            print(f"  {r.Reference} : {r.Test}   {r.Loci_matched}/{r.Loci_typed} loci, Tanabe {r.Tanabe_percent}%  -> {r.Verdict}"
                  + (f"   [{r.Note}]" if r.Note else ""))
    if args.plots:
        from .plots import overview, rfu_trace
        figs = out / "figures"; figs.mkdir(exist_ok=True)
        allres = dict(rr.samples)
        for x in (rr.positive, rr.negative):
            if x is not None: allres[x.id] = x
        for sid, s in allres.items():
            if s.sizing is None: continue
            overview(s, str(figs / f"{sid}_overview"), project=args.project or "")
            rfu_trace(s, rr.bins, str(figs / f"{sid}_RFU_trace"))
    print(f"\nResults: {out / 'str_results.xlsx'}")
    return 0


def main(argv=None):
    p = argparse.ArgumentParser(prog="str-typer", description="STR genotyping and cell-line authentication from .fsa files")
    p.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    sub = p.add_subparsers(dest="cmd", required=True)
    a = sub.add_parser("inspect", help="show what is stored in .fsa files (dyes, size standard, run)")
    a.add_argument("files", nargs="+"); a.add_argument("--kit", help="also test size-standard matching with this kit")
    a.set_defaults(fn=cmd_inspect)
    a = sub.add_parser("validate-kit", help="check a kit file against an allelic ladder")
    a.add_argument("--kit", required=True); a.add_argument("--ladder", required=True)
    a.add_argument("--bins-out", help="write the measured bins to this CSV"); a.set_defaults(fn=cmd_validate)
    a = sub.add_parser("run", help="genotype samples and compare pairs")
    a.add_argument("--kit", required=True, help="bundled kit name (e.g. cla_identifiler_plus) or path to a kit YAML")
    g = a.add_mutually_exclusive_group(required=True)
    g.add_argument("--ladder", help="allelic ladder .fsa from the same run (recommended)")
    g.add_argument("--bins", help="bins CSV (Marker,Allele,bp) e.g. from validate-kit --bins-out")
    a.add_argument("--samples", nargs="+", required=True)
    a.add_argument("--positive", help="positive control .fsa (checked against control_dna in the kit file)")
    a.add_argument("--negative", help="negative control .fsa")
    a.add_argument("--pair", action="append", help="REFERENCE:TEST, e.g. 3371p:3371 (repeatable)")
    a.add_argument("--set", action="append", help="override an analysis setting, e.g. analytical_threshold_rfu=100")
    a.add_argument("--plots", action="store_true", help="write overview and RFU-trace figures (PNG + 600 dpi TIFF)")
    a.add_argument("--project", help="project name printed on overview figures")
    a.add_argument("--out", default="str_results"); a.set_defaults(fn=cmd_run)
    args = p.parse_args(argv)
    return args.fn(args) or 0


if __name__ == "__main__":
    sys.exit(main())
