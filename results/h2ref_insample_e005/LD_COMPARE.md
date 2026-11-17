# Matched (in-sample) LD vs 1000 Genomes reference LD: h2_gwas = 0.003, h2_eqtl = 0.05

Loci in both runs: 4800 (50 replicates x 12 windows x 8 scenarios). Reference: `outputs/bulk_h2ref_v1/scores_grid.tsv.gz`; in-sample: `outputs/bulk_h2ref_v1_insample_g0.003_e0.05/scores_main.tsv.gz`.

**Check.** Whole-region coloc.abf PP4, reference vs in-sample: max |diff| = 0; NA pattern identical: True.

## A. Matched 5% FPR (cut-offs from the calibration windows of each LD), test windows

Detection rate / false-positive rate (single_distinct + double_distinct) / gwas_only, %; cut-off in brackets.

| method | in-sample LD | reference LD | detection: reference − in-sample |
|---|---|---|---|
| coloc.abf | 89.8 / 5.0 / 100.0 [0.012] | 89.8 / 5.0 / 100.0 [0.012] | +0.0 |
| block-coloc (conditional leads, 1e-3) | 97.5 / 6.2 / 100.0 [0.014] | 97.7 / 6.8 / 100.0 [0.014] | +0.2 |
| block-coloc (marginal leads, 1e-5) | 95.3 / 5.5 / 100.0 [0.014] | 95.2 / 5.7 / 100.0 [0.014] | -0.1 |
| block-coloc (original handoff, 1e-4) | 96.6 / 6.0 / 99.0 [0.022] | 96.4 / 4.5 / 91.3 [0.028] | -0.2 |
| coloc.susie | 81.2 / 4.2 / 0.0 [0.015] | 81.2 / 3.2 / 0.0 [0.015] | -0.1 |
| coloc.susie + fallback | 87.5 / 4.8 / 100.0 [0.020] | 87.2 / 5.2 / 99.7 [0.020] | -0.2 |
| coloc.signals (mask) | 95.8 / 7.7 / 21.0 [0.072] | 94.8 / 6.2 / 14.0 [0.092] | -1.1 |
| coloc.signals (cond) | 97.9 / 5.5 / 11.0 [0.101] | 98.2 / 7.5 / 4.3 [0.152] | +0.3 |

By scenario (detection or false-positive rate, in-sample → reference):

| method | single_shared | single_double_shared | double_shared2 | double_shared1 | single_distinct | double_distinct |
|---|---|---|---|---|---|---|
| coloc.abf | 100.0 → 100.0 | 90.3 → 90.3 | 94.7 → 94.7 | 74.3 → 74.3 | 0.0 → 0.0 | 10.0 → 10.0 |
| block-coloc (conditional leads, 1e-3) | 100.0 → 100.0 | 97.3 → 97.3 | 98.3 → 98.0 | 94.3 → 95.3 | 1.0 → 1.7 | 11.3 → 12.0 |
| block-coloc (marginal leads, 1e-5) | 100.0 → 100.0 | 96.0 → 96.0 | 97.3 → 97.0 | 88.0 → 88.0 | 0.3 → 0.7 | 10.7 → 10.7 |
| block-coloc (original handoff, 1e-4) | 100.0 → 100.0 | 97.0 → 96.7 | 98.3 → 98.0 | 91.0 → 91.0 | 3.0 → 3.0 | 9.0 → 6.0 |
| coloc.susie | 99.7 → 99.7 | 76.3 → 75.3 | 84.7 → 86.7 | 64.3 → 63.0 | 0.0 → 0.0 | 8.3 → 6.3 |
| coloc.susie + fallback | 100.0 → 100.0 | 81.0 → 82.0 | 96.0 → 96.3 | 73.0 → 70.7 | 0.0 → 0.0 | 9.7 → 10.3 |
| coloc.signals (mask) | 100.0 → 99.7 | 95.3 → 92.3 | 97.0 → 97.7 | 91.0 → 89.3 | 3.3 → 4.3 | 12.0 → 8.0 |
| coloc.signals (cond) | 99.7 → 99.7 | 96.7 → 97.3 | 99.3 → 99.3 | 96.0 → 96.3 | 1.7 → 5.7 | 9.3 → 9.3 |

## B. Fixed rule PP4 > 0.7, test windows

| method | colocalised loci called, in-sample → reference | not-colocalised loci called, in-sample → reference |
|---|---|---|
| coloc.abf | 68.2 → 68.2 | 0.0 → 0.0 |
| block-coloc (conditional leads, 1e-3) | 91.9 → 92.3 | 0.0 → 0.2 |
| block-coloc (marginal leads, 1e-5) | 85.6 → 85.8 | 0.3 → 0.3 |
| block-coloc (original handoff, 1e-4) | 90.9 → 91.3 | 0.8 → 1.0 |
| coloc.susie | 77.4 → 77.3 | 0.0 → 0.0 |
| coloc.susie + fallback | 82.2 → 81.8 | 0.0 → 0.0 |
| coloc.signals (mask) | 88.2 → 86.8 | 0.0 → 0.2 |
| coloc.signals (cond) | 93.2 → 92.8 | 0.0 → 1.2 |

## C. Paired change in detection, reference − in-sample (same positive test loci, matched 5% FPR)

| method | locus-level diff, McNemar p | window mean (t 95% CI) | windows < 0 | per-window |
|---|---|---|---|---|
| coloc.abf | +0.0, p = 1 | +0.0 (+0.0, +0.0) | 0/6 | +0.0, +0.0, +0.0, +0.0, +0.0, +0.0 |
| block-coloc (conditional leads, 1e-3) | +0.2, p = 0.73 | +0.2 (-0.5, +0.8) | 2/6 | +0.0, +1.0, +0.5, -0.5, -0.5, +0.5 |
| block-coloc (marginal leads, 1e-5) | -0.1, p = 1 | -0.1 (-0.3, +0.1) | 1/6 | +0.0, +0.0, +0.0, -0.5, +0.0, +0.0 |
| block-coloc (original handoff, 1e-4) | -0.2, p = 0.62 | -0.2 (-0.4, +0.1) | 2/6 | +0.0, -0.5, +0.0, -0.5, +0.0, +0.0 |
| coloc.susie | -0.1, p = 1 | -0.1 (-2.0, +1.8) | 3/6 | -2.0, -1.5, +3.0, -1.0, +0.5, +0.5 |
| coloc.susie + fallback | -0.2, p = 0.81 | -0.3 (-2.0, +1.5) | 2/6 | -3.5, -0.5, +1.0, +0.0, +0.5, +1.0 |
| coloc.signals (mask) | -1.1, p = 0.011 | -1.1 (-2.3, +0.2) | 4/6 | -3.0, -0.5, +0.0, -2.0, +0.0, -1.0 |
| coloc.signals (cond) | +0.2, p = 0.58 | +0.3 (-0.6, +1.1) | 2/6 | -0.5, +0.0, +1.5, +0.0, -0.5, +1.0 |

## D. Difference in differences: (block-coloc (conditional leads, 1e-3) − method) under reference LD minus the same under in-sample LD

Positive: the method loses more than block-coloc when the matched LD is replaced by the reference LD.

| method | block-coloc − method, in-sample | block-coloc − method, reference | difference (window t 95% CI) | windows > 0 |
|---|---|---|---|---|
| block-coloc (marginal leads, 1e-5) | +2.2 | +2.4 | +0.2 (-0.3, +0.8) | 3/6 |
| block-coloc (original handoff, 1e-4) | +0.9 | +1.3 | +0.3 (-0.4, +1.1) | 3/6 |
| coloc.susie | +16.2 | +16.5 | +0.2 (-1.7, +2.2) | 3/6 |
| coloc.susie + fallback | +10.0 | +10.4 | +0.4 (-1.4, +2.2) | 2/6 |
| coloc.signals (mask) | +1.7 | +2.9 | +1.2 (+0.0, +2.5) | 5/6 |
| coloc.signals (cond) | -0.4 | -0.5 | -0.1 (-0.9, +0.7) | 2/6 |

## E. Signal counts and no-call rates (all windows)

Conditional lead selection (lead threshold 1e-3): share of loci where the number of leads equals the number of causal variants, in-sample → reference; and share taking the block branch.

| scenario | trait 1: leads = causal | trait 2: leads = causal | block branch | no valid block |
|---|---|---|---|---|
| single_shared | 89.3 → 90.2 | 93.5 → 90.8 | 16.2 → 18.0 | 0.2 → 0.2 |
| single_double_shared | 91.2 → 90.0 | 84.8 → 83.3 | 90.2 → 91.0 | 3.5 → 4.5 |
| double_shared2 | 88.5 → 85.5 | 82.0 → 78.5 | 98.5 → 98.0 | 2.5 → 1.8 |
| double_shared1 | 85.7 → 82.2 | 85.7 → 81.0 | 99.2 → 99.3 | 11.5 → 11.5 |
| single_distinct | 91.3 → 90.0 | 94.5 → 91.0 | 14.0 → 18.5 | 13.0 → 17.5 |
| double_distinct | 85.7 → 84.2 | 85.5 → 81.7 | 99.0 → 99.3 | 97.5 → 97.3 |
| gwas_only | 91.7 → 91.3 | 86.2 → 86.2 | 1.7 → 1.3 | 1.7 → 1.3 |
| null_both | 85.3 → 85.3 | 86.5 → 86.5 | 0.3 → 0.3 | 0.3 → 0.3 |

No result (PP4 = NA), positives on test windows, in-sample → reference: coloc.susie 7.2 → 7.2; coloc.susie + fallback 0.0 → 0.0; coloc.signals (mask) 0.1 → 0.0; coloc.signals (cond) 0.0 → 0.0
