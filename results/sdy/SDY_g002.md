# sdY sensitivity: h2_gwas = 0.02, h2_eqtl = 0.015

4800 loci, arms: `est` = sdY estimated on each subset (current); `whole` = sdY estimated once on the whole region; `one` = sdY = 1 (true value).
Leads and blocks are the same in every arm; each method and arm gets its own cut-off at 5% FPR on the calibration windows.

**Reproduction check** (arm est vs `outputs/bulk_h2ref_v1/scores_grid.tsv.gz`): 4800 loci matched; block PP4 max |diff| = 0; NA pattern identical: True; whole-region PP4 max |diff| = 0.

**Internal check**: whole-region PP4, arm whole vs arm est, max |diff| = 0 (should be ~0: same sdY); own sdY estimate vs coloc:::sdY.est, max |diff| = 0; loci with NA sdy_check: 0.

## Matched 5% FPR, test windows

| arm | method | cut-off | detection | FPR (distinct) | gwas_only | block-coloc − coloc.abf (pp), McNemar p | window mean (t 95% CI), windows > 0 |
|---|---|---|---|---|---|---|---|
| est | block-coloc | 0.1471 | 76.3 | 4.3 | 5.0 |  |  |
| est | coloc.abf | 0.1391 | 70.1 | 3.0 | 5.7 | +6.2, p = 7.6e-17 | +6.3 (+3.0, +9.5), 6/6 |
| whole | block-coloc | 0.1471 | 76.3 | 4.3 | 5.0 |  |  |
| whole | coloc.abf | 0.1391 | 70.1 | 3.0 | 5.7 | +6.2, p = 7.6e-17 | +6.3 (+3.0, +9.5), 6/6 |
| one | block-coloc | 0.1469 | 76.3 | 4.3 | 5.0 |  |  |
| one | coloc.abf | 0.1378 | 70.4 | 3.0 | 6.0 | +5.9, p = 1.4e-14 | +5.9 (+2.5, +9.3), 6/6 |

Block-coloc, arm vs arm est (same positive test loci):

| arm | diff vs est (pp), McNemar p | window mean (t 95% CI) |
|---|---|---|
| whole | +0.0, p = 1 | +0.0 (+0.0, +0.0), 0/6 |
| one | +0.0, p = 1 | +0.0 (+0.0, +0.0), 0/6 |

## By scenario (block-coloc / coloc.abf, %)

| arm | single_shared | single_double_shared | double_shared2 | double_shared1 | single_distinct | double_distinct |
|---|---|---|---|---|---|---|
| est | 96.3 / 96.3 | 72.7 / 73.3 | 84.3 / 72.0 | 52.0 / 38.7 | 3.0 / 2.3 | 5.7 / 3.7 |
| whole | 96.3 / 96.3 | 72.7 / 73.3 | 84.3 / 72.0 | 52.0 / 38.7 | 3.0 / 2.3 | 5.7 / 3.7 |
| one | 96.3 / 96.3 | 72.7 / 73.3 | 84.3 / 72.7 | 52.0 / 39.3 | 3.0 / 2.3 | 5.7 / 3.7 |

## Fixed rule (PP4 > 0.7 colocalised; PP3 > 0.7 independent), test windows, block-coloc

| arm | colocalised: called coloc / indep | single_distinct: indep (coloc) | double_distinct: indep (coloc) | gwas_only: coloc |
|---|---|---|---|---|
| est | 47.2 / 3.2 | 34.7 (0.0) | 9.0 (0.3) | 0.0 |
| whole | 47.2 / 3.2 | 34.7 (0.0) | 9.0 (0.3) | 0.0 |
| one | 47.2 / 3.2 | 35.0 (0.0) | 9.0 (0.3) | 0.0 |

## Diagnostics (all windows)

| quantity | median / p90 / p99 |
|---|---|
| sdY trait 1, whole region | 1.002 / 1.008 / 1.014 |
| sdY trait 2, whole region | 1.002 / 1.031 / 1.054 |
| sdY trait 1, best block / whole region (block branch, arm est) | 0.9962 / 1.001 / 1.005 |
| sdY trait 2, best block / whole region (block branch, arm est) | 0.9957 / 1.013 / 1.025 |
| abs(log sdY2 ratio), block branch | 0.009093 / 0.02221 / 0.03246 |
| block-coloc abs(PP4 whole − est), all loci | 0 / 9.078e-05 / 0.001802 |
| coloc.abf abs(PP4 whole − est), all loci | 0 / 0 / 0 |
| best block differs from arm est (whole), block branch | 0.0% |
| block-coloc abs(PP4 one − est), all loci | 0.0002944 / 0.001549 / 0.005253 |
| coloc.abf abs(PP4 one − est), all loci | 0.0003031 / 0.001598 / 0.005282 |
| best block differs from arm est (one), block branch | 0.0% |
| block subsets where the sdY estimate is NA (negative coefficient), block branch | 0.0% |
