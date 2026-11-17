# Matched (in-sample) LD vs 1000 Genomes reference LD: h2_gwas = 0.02, h2_eqtl = 0.015

Loci in both runs: 4800 (50 replicates x 12 windows x 8 scenarios). Reference: `outputs/bulk_h2ref_v1/scores_grid.tsv.gz`; in-sample: `outputs/bulk_h2ref_v1_insample_g0.02_e0.015/scores_main.tsv.gz`.

**Check.** Whole-region coloc.abf PP4, reference vs in-sample: max |diff| = 0; NA pattern identical: True.

## A. Matched 5% FPR (cut-offs from the calibration windows of each LD), test windows

Detection rate / false-positive rate (single_distinct + double_distinct) / gwas_only, %; cut-off in brackets.

| method | in-sample LD | reference LD | detection: reference − in-sample |
|---|---|---|---|
| coloc.abf | 70.1 / 3.0 / 5.7 [0.139] | 70.1 / 3.0 / 5.7 [0.139] | +0.0 |
| block-coloc (conditional leads, 1e-3) | 76.8 / 3.7 / 5.3 [0.143] | 76.3 / 4.3 / 5.0 [0.147] | -0.4 |
| block-coloc (marginal leads, 1e-5) | 70.6 / 5.0 / 5.0 [0.153] | 70.6 / 4.7 / 5.0 [0.153] | +0.0 |
| block-coloc (original handoff, 1e-4) | 66.2 / 6.5 / 4.7 [0.391] | 66.4 / 6.0 / 4.7 [0.393] | +0.3 |
| coloc.susie | 18.0 / 4.0 / 0.0 [0.012] | 19.1 / 2.3 / 0.0 [0.015] | +1.1 |
| coloc.susie + fallback | 72.2 / 4.0 / 6.7 [0.130] | 70.9 / 3.0 / 5.7 [0.139] | -1.3 |
| coloc.signals (mask) | 78.7 / 4.5 / 3.3 [0.196] | 75.2 / 4.5 / 3.3 [0.238] | -3.5 |
| coloc.signals (cond) | 82.9 / 5.2 / 4.3 [0.163] | 76.2 / 6.5 / 2.0 [0.270] | -6.7 |

By scenario (detection or false-positive rate, in-sample → reference):

| method | single_shared | single_double_shared | double_shared2 | double_shared1 | single_distinct | double_distinct |
|---|---|---|---|---|---|---|
| coloc.abf | 96.3 → 96.3 | 73.3 → 73.3 | 72.0 → 72.0 | 38.7 → 38.7 | 2.3 → 2.3 | 3.7 → 3.7 |
| block-coloc (conditional leads, 1e-3) | 96.3 → 96.3 | 73.3 → 72.7 | 85.7 → 84.3 | 51.7 → 52.0 | 2.3 → 3.0 | 5.0 → 5.7 |
| block-coloc (marginal leads, 1e-5) | 96.3 → 96.3 | 71.7 → 72.3 | 74.7 → 74.3 | 39.7 → 39.3 | 4.3 → 4.3 | 5.7 → 5.0 |
| block-coloc (original handoff, 1e-4) | 92.3 → 92.0 | 53.0 → 53.0 | 74.3 → 74.3 | 45.0 → 46.3 | 6.0 → 5.3 | 7.0 → 6.7 |
| coloc.susie | 34.7 → 35.7 | 11.0 → 11.7 | 15.0 → 16.7 | 11.3 → 12.3 | 5.0 → 2.3 | 3.0 → 2.3 |
| coloc.susie + fallback | 96.0 → 95.7 | 71.7 → 71.7 | 79.3 → 76.0 | 41.7 → 40.3 | 3.7 → 2.0 | 4.3 → 4.0 |
| coloc.signals (mask) | 95.0 → 94.7 | 69.7 → 64.7 | 87.3 → 84.3 | 62.7 → 57.0 | 3.0 → 3.0 | 6.0 → 6.0 |
| coloc.signals (cond) | 95.7 → 94.3 | 73.3 → 61.7 | 94.0 → 86.7 | 68.7 → 62.3 | 2.3 → 6.7 | 8.0 → 6.3 |

## B. Fixed rule PP4 > 0.7, test windows

| method | colocalised loci called, in-sample → reference | not-colocalised loci called, in-sample → reference |
|---|---|---|
| coloc.abf | 36.6 → 36.6 | 0.0 → 0.0 |
| block-coloc (conditional leads, 1e-3) | 47.1 → 47.2 | 0.0 → 0.2 |
| block-coloc (marginal leads, 1e-5) | 38.9 → 38.9 | 0.7 → 0.7 |
| block-coloc (original handoff, 1e-4) | 49.2 → 49.2 | 2.3 → 2.7 |
| coloc.susie | 15.2 → 15.9 | 0.0 → 0.0 |
| coloc.susie + fallback | 38.8 → 38.7 | 0.0 → 0.0 |
| coloc.signals (mask) | 46.2 → 45.5 | 0.0 → 0.2 |
| coloc.signals (cond) | 48.2 → 48.2 | 0.0 → 1.5 |

## C. Paired change in detection, reference − in-sample (same positive test loci, matched 5% FPR)

| method | locus-level diff, McNemar p | window mean (t 95% CI) | windows < 0 | per-window |
|---|---|---|---|---|
| coloc.abf | +0.0, p = 1 | +0.0 (+0.0, +0.0) | 0/6 | +0.0, +0.0, +0.0, +0.0, +0.0, +0.0 |
| block-coloc (conditional leads, 1e-3) | -0.4, p = 0.27 | -0.4 (-1.3, +0.4) | 4/6 | +0.5, -1.5, -1.0, -0.5, -0.5, +0.5 |
| block-coloc (marginal leads, 1e-5) | +0.0, p = 1 | +0.0 (-0.5, +0.5) | 2/6 | -0.5, -0.5, +0.5, +0.0, +0.5, +0.0 |
| block-coloc (original handoff, 1e-4) | +0.2, p = 0.68 | +0.2 (-1.8, +2.3) | 2/6 | +1.0, -0.5, -3.5, +1.5, +1.5, +1.5 |
| coloc.susie | +1.1, p = 0.053 | +1.1 (-0.4, +2.5) | 2/6 | +1.0, +2.0, -0.5, -0.5, +1.5, +3.0 |
| coloc.susie + fallback | -1.2, p = 0.014 | -1.3 (-2.2, -0.3) | 5/6 | -1.0, -2.5, -2.0, -1.0, -1.0, +0.0 |
| coloc.signals (mask) | -3.5, p = 1.6e-08 | -3.5 (-4.9, -2.1) | 6/6 | -5.0, -4.5, -2.5, -2.0, -2.5, -4.5 |
| coloc.signals (cond) | -6.7, p = 2.7e-21 | -6.7 (-7.8, -5.5) | 6/6 | -7.5, -6.0, -5.5, -8.0, -5.5, -7.5 |

## D. Difference in differences: (block-coloc (conditional leads, 1e-3) − method) under reference LD minus the same under in-sample LD

Positive: the method loses more than block-coloc when the matched LD is replaced by the reference LD.

| method | block-coloc − method, in-sample | block-coloc − method, reference | difference (window t 95% CI) | windows > 0 |
|---|---|---|---|---|
| block-coloc (marginal leads, 1e-5) | +6.2 | +5.8 | -0.4 (-1.4, +0.6) | 2/6 |
| block-coloc (original handoff, 1e-4) | +10.6 | +9.9 | -0.7 (-2.4, +1.1) | 1/6 |
| coloc.susie | +58.7 | +57.2 | -1.5 (-2.9, -0.1) | 0/6 |
| coloc.susie + fallback | +4.6 | +5.4 | +0.8 (+0.4, +1.3) | 6/6 |
| coloc.signals (mask) | -1.9 | +1.2 | +3.1 (+1.2, +4.9) | 6/6 |
| coloc.signals (cond) | -6.2 | +0.1 | +6.2 (+4.4, +8.1) | 6/6 |

## E. Signal counts and no-call rates (all windows)

Conditional lead selection (lead threshold 1e-3): share of loci where the number of leads equals the number of causal variants, in-sample → reference; and share taking the block branch.

| scenario | trait 1: leads = causal | trait 2: leads = causal | block branch | no valid block |
|---|---|---|---|---|
| single_shared | 92.5 → 76.3 | 76.5 → 76.7 | 13.0 → 24.2 | 0.7 → 1.7 |
| single_double_shared | 93.3 → 71.2 | 14.0 → 15.3 | 18.3 → 31.0 | 6.2 → 12.2 |
| double_shared2 | 90.0 → 67.8 | 15.0 → 15.5 | 70.7 → 71.5 | 10.8 → 11.2 |
| double_shared1 | 91.7 → 66.7 | 12.5 → 12.7 | 70.2 → 69.7 | 37.7 → 36.5 |
| single_distinct | 94.2 → 73.2 | 77.0 → 76.2 | 11.0 → 27.3 | 10.5 → 26.0 |
| double_distinct | 90.5 → 69.3 | 15.7 → 15.7 | 69.0 → 69.2 | 66.2 → 66.2 |
| gwas_only | 93.7 → 75.2 | 86.2 → 86.2 | 1.2 → 3.5 | 1.2 → 3.5 |
| null_both | 85.3 → 85.3 | 86.5 → 86.5 | 0.3 → 0.3 | 0.3 → 0.3 |

No result (PP4 = NA), positives on test windows, in-sample → reference: coloc.susie 79.3 → 78.4; coloc.susie + fallback 0.0 → 0.0; coloc.signals (mask) 0.0 → 0.1; coloc.signals (cond) 0.0 → 0.0
