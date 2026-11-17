# Matched (in-sample) LD vs 1000 Genomes reference LD: h2_gwas = 0.003, h2_eqtl = 0.015

Loci in both runs: 9600 (100 replicates x 12 windows x 8 scenarios). Reference: `outputs/bulk_h2ref_v1/scores_main.tsv.gz`; in-sample: `outputs/bulk_h2ref_v1_insample/scores_main.tsv.gz`.

**Check.** Whole-region coloc.abf PP4, reference vs in-sample: max |diff| = 0; NA pattern identical: True.

## A. Matched 5% FPR (cut-offs from the calibration windows of each LD), test windows

Detection rate / false-positive rate (single_distinct + double_distinct) / gwas_only, %; cut-off in brackets.

| method | in-sample LD | reference LD | detection: reference − in-sample |
|---|---|---|---|
| coloc.abf | 73.2 / 3.5 / 5.0 [0.141] | 73.2 / 3.5 / 5.0 [0.141] | +0.0 |
| block-coloc (conditional leads, 1e-3) | 77.0 / 3.8 / 4.8 [0.143] | 76.8 / 3.8 / 4.8 [0.145] | -0.3 |
| block-coloc (marginal leads, 1e-3) | 76.0 / 5.0 / 4.3 [0.156] | 74.9 / 4.9 / 3.7 [0.168] | -1.1 |
| block-coloc (marginal leads, 1e-5) | 74.0 / 3.7 / 5.0 [0.141] | 74.0 / 3.6 / 5.0 [0.142] | +0.0 |
| block, no gate (marginal leads, 1e-3) | 76.0 / 5.0 / 4.3 [0.156] | 74.9 / 4.9 / 3.7 [0.168] | -1.1 |
| block-coloc (original handoff, 1e-4) | 77.7 / 4.6 / 5.0 [0.154] | 77.3 / 4.6 / 4.3 [0.159] | -0.4 |
| coloc.susie | 17.9 / 3.7 / 0.2 [0.011] | 18.3 / 3.6 / 0.0 [0.011] | +0.4 |
| coloc.susie + fallback | 72.5 / 3.5 / 5.3 [0.139] | 72.6 / 3.4 / 5.3 [0.139] | +0.0 |
| coloc.signals (mask) | 78.5 / 4.8 / 3.5 [0.170] | 78.2 / 4.6 / 3.5 [0.171] | -0.3 |
| coloc.signals (cond) | 79.5 / 4.4 / 3.5 [0.176] | 76.2 / 3.8 / 2.5 [0.210] | -3.3 |

By scenario (detection or false-positive rate, in-sample → reference):

| method | single_shared | single_double_shared | double_shared2 | double_shared1 | single_distinct | double_distinct |
|---|---|---|---|---|---|---|
| coloc.abf | 96.3 → 96.3 | 74.3 → 74.3 | 76.8 → 76.8 | 45.3 → 45.3 | 3.5 → 3.5 | 3.5 → 3.5 |
| block-coloc (conditional leads, 1e-3) | 96.3 → 96.2 | 74.5 → 74.2 | 82.2 → 81.8 | 55.0 → 54.8 | 3.5 → 3.7 | 4.2 → 3.8 |
| block-coloc (marginal leads, 1e-3) | 95.7 → 95.7 | 72.0 → 70.5 | 81.2 → 78.8 | 55.2 → 54.5 | 5.5 → 5.5 | 4.5 → 4.3 |
| block-coloc (marginal leads, 1e-5) | 96.3 → 96.3 | 74.8 → 74.8 | 78.2 → 78.0 | 46.8 → 47.0 | 3.7 → 3.7 | 3.7 → 3.5 |
| block, no gate (marginal leads, 1e-3) | 95.7 → 95.7 | 72.0 → 70.5 | 81.2 → 78.8 | 55.2 → 54.5 | 5.5 → 5.5 | 4.5 → 4.3 |
| block-coloc (original handoff, 1e-4) | 95.8 → 95.7 | 73.3 → 72.3 | 84.5 → 84.2 | 57.2 → 57.0 | 4.2 → 4.3 | 5.0 → 4.8 |
| coloc.susie | 36.7 → 37.5 | 10.8 → 10.8 | 13.2 → 14.0 | 11.0 → 11.0 | 4.0 → 4.2 | 3.3 → 3.0 |
| coloc.susie + fallback | 96.2 → 96.0 | 70.8 → 71.8 | 78.0 → 77.3 | 45.2 → 45.2 | 3.5 → 3.3 | 3.5 → 3.5 |
| coloc.signals (mask) | 95.5 → 95.5 | 71.3 → 71.3 | 84.7 → 84.3 | 62.5 → 61.8 | 3.7 → 3.7 | 5.8 → 5.5 |
| coloc.signals (cond) | 95.5 → 94.8 | 71.3 → 67.7 | 86.5 → 82.0 | 64.5 → 60.2 | 3.3 → 3.2 | 5.5 → 4.3 |

## B. Fixed rule PP4 > 0.7, test windows

| method | colocalised loci called, in-sample → reference | not-colocalised loci called, in-sample → reference |
|---|---|---|
| coloc.abf | 36.5 → 36.5 | 0.1 → 0.1 |
| block-coloc (conditional leads, 1e-3) | 44.1 → 44.1 | 0.2 → 0.2 |
| block-coloc (marginal leads, 1e-3) | 45.4 → 45.4 | 0.5 → 0.6 |
| block-coloc (marginal leads, 1e-5) | 38.1 → 38.0 | 0.1 → 0.1 |
| block, no gate (marginal leads, 1e-3) | 45.4 → 45.5 | 0.5 → 0.6 |
| block-coloc (original handoff, 1e-4) | 44.4 → 44.5 | 0.2 → 0.3 |
| coloc.susie | 14.9 → 15.5 | 0.0 → 0.0 |
| coloc.susie + fallback | 37.6 → 37.5 | 0.1 → 0.1 |
| coloc.signals (mask) | 42.8 → 42.6 | 0.1 → 0.1 |
| coloc.signals (cond) | 44.5 → 44.9 | 0.1 → 0.2 |

## C. Paired change in detection, reference − in-sample (same positive test loci, matched 5% FPR)

| method | locus-level diff, McNemar p | window mean (t 95% CI) | windows < 0 | per-window |
|---|---|---|---|---|
| coloc.abf | +0.0, p = 1 | +0.0 (+0.0, +0.0) | 0/6 | +0.0, +0.0, +0.0, +0.0, +0.0, +0.0 |
| block-coloc (conditional leads, 1e-3) | -0.2, p = 0.18 | -0.3 (-0.8, +0.3) | 3/6 | +0.0, -0.5, -1.0, +0.5, -0.5, +0.0 |
| block-coloc (marginal leads, 1e-3) | -1.1, p = 4.2e-05 | -1.1 (-1.7, -0.5) | 6/6 | -0.8, -1.8, -0.7, -0.8, -2.0, -0.7 |
| block-coloc (marginal leads, 1e-5) | +0.0, p = 1 | +0.0 (-0.2, +0.2) | 1/6 | +0.0, +0.0, -0.2, +0.2, +0.0, +0.0 |
| block, no gate (marginal leads, 1e-3) | -1.1, p = 4.2e-05 | -1.1 (-1.7, -0.5) | 6/6 | -0.8, -1.8, -0.7, -0.8, -2.0, -0.7 |
| block-coloc (original handoff, 1e-4) | -0.4, p = 0.087 | -0.4 (-1.2, +0.4) | 5/6 | -0.5, -1.2, -0.2, +1.0, -0.5, -1.0 |
| coloc.susie | +0.4, p = 0.27 | +0.4 (-0.2, +1.0) | 1/6 | +0.0, +1.3, +0.5, +0.3, -0.2, +0.8 |
| coloc.susie + fallback | +0.0, p = 1 | +0.0 (-0.3, +0.3) | 1/6 | -0.5, +0.2, +0.3, +0.0, +0.2, +0.0 |
| coloc.signals (mask) | -0.2, p = 0.34 | -0.3 (-0.8, +0.3) | 4/6 | -0.8, -0.5, -0.5, +0.5, +0.2, -0.5 |
| coloc.signals (cond) | -3.3, p = 3e-20 | -3.3 (-4.4, -2.2) | 6/6 | -3.7, -2.5, -4.0, -4.7, -2.2, -2.5 |

## D. Difference in differences: (block-coloc (conditional leads, 1e-3) − method) under reference LD minus the same under in-sample LD

Positive: the method loses more than block-coloc when the matched LD is replaced by the reference LD.

| method | block-coloc − method, in-sample | block-coloc − method, reference | difference (window t 95% CI) | windows > 0 |
|---|---|---|---|---|
| block-coloc (marginal leads, 1e-3) | +1.0 | +1.9 | +0.9 (+0.2, +1.5) | 5/6 |
| block-coloc (marginal leads, 1e-5) | +3.0 | +2.7 | -0.3 (-0.7, +0.2) | 1/6 |
| block, no gate (marginal leads, 1e-3) | +1.0 | +1.9 | +0.9 (+0.2, +1.5) | 5/6 |
| block-coloc (original handoff, 1e-4) | -0.7 | -0.5 | +0.2 (-0.6, +0.9) | 3/6 |
| coloc.susie | +59.1 | +58.4 | -0.7 (-1.5, +0.2) | 1/6 |
| coloc.susie + fallback | +4.5 | +4.2 | -0.3 (-1.1, +0.5) | 2/6 |
| coloc.signals (mask) | -1.5 | -1.5 | +0.0 (-0.6, +0.6) | 2/6 |
| coloc.signals (cond) | -2.5 | +0.6 | +3.0 (+1.7, +4.4) | 6/6 |

## E. Signal counts and no-call rates (all windows)

Conditional lead selection (lead threshold 1e-3): share of loci where the number of leads equals the number of causal variants, in-sample → reference; and share taking the block branch.

| scenario | trait 1: leads = causal | trait 2: leads = causal | block branch | no valid block |
|---|---|---|---|---|
| single_shared | 90.2 → 90.5 | 76.8 → 76.9 | 14.6 → 14.2 | 1.0 → 1.1 |
| single_double_shared | 91.2 → 90.0 | 13.9 → 15.9 | 19.8 → 22.7 | 6.4 → 7.4 |
| double_shared2 | 87.5 → 84.0 | 15.4 → 15.4 | 67.0 → 66.8 | 12.3 → 12.8 |
| double_shared1 | 86.4 → 82.8 | 13.4 → 13.4 | 66.1 → 65.7 | 34.1 → 33.8 |
| single_distinct | 90.7 → 88.1 | 76.6 → 75.2 | 12.8 → 16.8 | 12.2 → 15.9 |
| double_distinct | 86.4 → 85.8 | 15.4 → 15.8 | 65.8 → 65.8 | 63.7 → 64.0 |
| gwas_only | 92.2 → 91.2 | 86.8 → 86.8 | 1.8 → 1.9 | 1.8 → 1.9 |
| null_both | 86.1 → 86.1 | 85.5 → 85.5 | 0.4 → 0.2 | 0.4 → 0.2 |

No result (PP4 = NA), positives on test windows, in-sample → reference: coloc.susie 79.7 → 79.2; coloc.susie + fallback 0.0 → 0.0; coloc.signals (mask) 0.0 → 0.0; coloc.signals (cond) 0.0 → 0.0
