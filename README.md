# fsa-str-typer

STR genotyping and cell-line authentication directly from raw capillary-electrophoresis files
(`.fsa` / `.hid`, ABIF format), without the instrument vendor's analysis software.

Typical use: an iPSC line and the donor's PBMC are typed by a service provider, you receive the
`.fsa` files, and you want an independent, reproducible answer to *"is this line derived from this
donor?"* together with publication-ready figures and tables.

```text
str-typer run --kit cla_identifiler_plus --ladder Ladder.fsa \
    --positive P.fsa --negative N.fsa \
    --samples PBMC_01.fsa iPSC_01.fsa PBMC_02.fsa iPSC_02.fsa \
    --pair PBMC_01:iPSC_01 --pair PBMC_02:iPSC_02 --plots --out results/

Positive control: 16/16 loci concordant
Pairs (reference : test)
  PBMC_01 : iPSC_01   16/16 loci, Tanabe 100.0%  -> Match
  PBMC_02 : iPSC_02    1/16 loci, Tanabe 46.4%   -> No match  [iPSC_02 matches PBMC_03 (100%) - possible sample label swap]
```

Everything runs locally; no data leave your computer.

---

## Installation

```bash
git clone https://github.com/waynechen228-design/RoadSTR.git
cd RoadSTR
conda env create -f environment.yml      # or:  pip install -e .
conda activate fsa-str-typer
pytest                                    # runs the tests on the bundled example files
```

Requires Python ≥ 3.9 and numpy, scipy, pandas, biopython, pyyaml, openpyxl, matplotlib.

## What you need

| Input | Why |
|---|---|
| sample `.fsa` files | the injections to genotype |
| **allelic ladder `.fsa` of the same run** | allele bins are measured from it, so run- and instrument-specific mobility is corrected automatically. Without it (`--bins file.csv` from another run) calls still work but are less precise |
| positive control `.fsa` (recommended) | checked against the control genotype in the kit file (e.g. 9947A) |
| negative control `.fsa` (recommended) | contamination check; also used to recognise run artifacts |
| a kit file | which loci, dyes and ladder alleles the kit has (bundled: `cla_identifiler_plus`) |

## Commands

```bash
str-typer inspect file.fsa [--kit KIT]          # what is in a file: dyes, size standard, run, saturation
str-typer validate-kit --kit KIT --ladder Ladder.fsa [--bins-out bins.csv]
str-typer run --kit KIT --ladder Ladder.fsa --samples *.fsa [--positive P.fsa] [--negative N.fsa]
              [--pair REF:TEST ...] [--plots] [--set name=value ...] [--out DIR]
```

`--pair` uses file names without `.fsa`. `--set` overrides an analysis setting for this run
(e.g. `--set analytical_threshold_rfu=100`); every setting used is written to the report.

## Output

`str_results.xlsx`

| Sheet | Content |
|---|---|
| Pairs | per pair: loci matched, Tanabe %, verdict, label-swap hint |
| Genotypes | one column per sample |
| Alleles | every called allele: size, deviation from the bin, height, flags, status, notes |
| QC | per injection: sizing, loci called / to check, max deviation, min height, heterozygote balance, saturation |
| Positive_control / Negative_control | concordance with the kit's control DNA; signals in the negative control |
| Cross_match | Tanabe % between all profiles of the run |
| Peaks | every detected peak with its classification (allele, pull-up, spike, …) |
| Ladder, Bins | ladder check and the bins measured from it |
| Settings | program version, kit file, all analysis settings, run notes |

With `--plots`, for every injection (PNG and 600-dpi LZW TIFF, 6.5 in wide, 7-pt text):
`*_overview` (all channels vs. scan number, as in the vendor software) and
`*_RFU_trace` (one panel per dye vs. bp, alleles labelled).

Every locus gets a status: **OK**, **Check** (inspect the trace; the note says why) or **Fail**.

## How it works

1. **Reading** – raw dye channels are read from the ABIF file and a rolling low-percentile baseline is subtracted.
2. **Sizing** – size-standard peaks are matched to the fragment list by their spacing, walking back from the
   largest fragment and predicting each next position from the local scans-per-bp trend, so extra peaks
   (bleed-through from over-loaded samples, primer front) are skipped. Fragments inside the primer front are
   not used. Sizes are computed by the Local Southern method; the 250-bp fragment of GS500 is excluded, as in
   the vendor software. If too few fragments can be matched the sample is not called.
3. **Bins** – in the allelic ladder, each marker's allele pattern (from the allele names and the repeat unit)
   is fitted to the detected peaks, markers of one dye in size order. Every listed allele must be found;
   strong unassigned peaks are reported (usually an allele missing from the kit file).
4. **Peak classification** – every peak is labelled before calling:
   - *pull-up*: under a taller peak of another dye and no taller than the bleed ratio measured in this sample allows
     (a much taller peak is kept as a co-migrating allele);
   - *saturated*, *pull-up (saturated source)*, *co-migrates with saturated peak*, *under saturated peak (undetermined)*:
     off-scale peaks are detected from the raw signal; bleed under them is judged against the bleed under the other
     off-scale peaks of the same dye;
   - *spike* (narrower than 2 scans), *size-std bleed*, *stutter* (one repeat shorter, < 15 % of the parent),
     *shoulder* (within 1.6 bp of a much taller peak).
5. **Calling** – peaks above the analytical threshold are assigned to the nearest ladder allele within ±0.5 bp.
   Off-ladder alleles (e.g. FGA 22.2) are named from virtual bins (repeat-unit interpolation) and flagged.
   Up to two alleles are called per locus; further peaks ≥ 20 % of the main allele are reported (possible mixture,
   trisomy or artifact), heterozygote balance < 0.5 is flagged.
6. **Comparison** – Tanabe percent match over the STR loci (Amelogenin excluded), ≥ 80 % = match
   (ANSI/ATCC ASN-0002); every failing pair is compared with all other profiles to detect label swaps.

A suggested methods sentence:

> STR genotypes were called from the raw .fsa data with fsa-str-typer v0.1.0 (Local Southern sizing against
> GS500 LIZ; allele bins measured from the allelic ladder injected in the same run, ±0.5 bp). Profiles were
> compared by Tanabe percent match (ANSI/ATCC ASN-0002).

## Validation so far

Two 3730xl runs, CLA Identifiler Plus, GS500 LIZ, POP-7 (19 injections):

- allelic ladder: all 205 ladder alleles of the kit file found (`validate-kit`);
- 9947A positive control: 16/16 loci concordant in both runs;
- peak sizes agree with the vendor software within 0.2 bp;
- 7 PBMC/iPSC pairs (272 locus calls) identical to an independent manual review of the traces; one iPSC
  sample whose sizing failed in the vendor software was called completely; one sample label swap between
  two lines was detected automatically;
- over-loaded injections (signal off-scale in all four dyes) gave correct genotypes, with the affected loci flagged.

Other kits are supported through kit files but are **unverified** until checked with `validate-kit` on a real
ladder and with a positive control. Please report results (issues / pull requests welcome).

## Adding a kit

Copy `fsa_str_typer/kits/kit_template.yaml`, fill in markers, dyes, repeat units and ladder alleles from the
kit user guide (comments in the template say where to find each item), then:

```bash
str-typer validate-kit --kit my_kit.yaml --ladder Ladder.fsa
```

`PASS` = every ladder allele was found in the stated dye and order; `REVIEW` = all found, but strong
unexplained peaks exist (often an allele missing from the list). No base-pair positions are needed.

## Limitations

- For research use (cell-line authentication). Not validated for forensic, clinical or diagnostic use.
- Mixtures are detected (more than two alleles), not deconvolved.
- Over-loaded samples: genotypes are usually correct, but heights and balance are not measurable; re-inject diluted.
- Tested on 5-dye chemistry (4 dyes + size standard). 6-dye files are read, but no 6-dye kit has been validated yet.
- Always look at the traces of every locus marked **Check**.

## Data protection

STR profiles identify individuals. This repository contains only an allelic ladder and commercial control DNA.
Do not commit donor or patient files (`.gitignore` excludes `*.fsa` outside `examples/`).

## License

MIT
