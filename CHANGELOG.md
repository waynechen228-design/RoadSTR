# Changelog

## 0.1.0 (2026-11-25)
First public version. Developed and tested on two 3730xl runs (CLA Identifiler Plus, GS500 LIZ):
19 injections incl. allelic ladders, 9947A positive controls, negative controls and 7 PBMC/iPSC pairs.

Robustness rules added during testing:
- **Size-standard matching by spacing**, walking back from the largest fragment; extra peaks in the
  size-standard channel (bleed from over-loaded samples, primer front) no longer shift the assignment.
  Sizing that cannot be matched stops with an error instead of producing wrong sizes.
- **Bins measured from the allelic ladder of the same run** by fitting each marker's allele pattern
  (allele names + repeat unit); kit files are validated against a real ladder (`validate-kit`).
- **Cross-dye pull-up judged from the bleed ratio measured in each sample**, so real alleles that
  co-migrate with alleles of another dye (e.g. D3S1358 / D19S433 in 9947A) are kept.
- **Saturation handling**: off-scale peaks are detected from the raw signal; peaks under them are
  judged against the bleed under the other off-scale peaks of the same dye; undecidable cases are
  flagged (re-inject diluted) instead of being called or silently removed.
- **Electrical spikes** (narrower than any DNA fragment), **size-standard bleed-through**, **stutter**
  and **shoulders** (incomplete +A addition / split peaks) are classified as non-allelic.
- **Off-ladder alleles** (e.g. FGA 22.2) are named from virtual bins built from the repeat unit and
  flagged for review; off-ladder peaks under taller peaks of another dye or present in the negative
  control are not named.
- **Run artifacts**: signals of the negative control recur at the same scan in other capillaries; they
  are matched by scan and kept out of mixture warnings (never used to remove a main allele).
- **Pair comparison** with Tanabe % (ANSI/ATCC ASN-0002) and automatic detection of sample label swaps
  across all profiles of the run.
