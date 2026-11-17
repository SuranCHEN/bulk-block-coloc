# sdY sensitivity: h2_gwas = 0.003, h2_eqtl = 0.015

9600 loci, arms: `est` = sdY estimated on each subset (current); `whole` = sdY estimated once on the whole region; `one` = sdY = 1 (true value).
Leads and blocks are the same in every arm; each method and arm gets its own cut-off at 5% FPR on the calibration windows.

**Reproduction check** (arm est vs `outputs/bulk_h2ref_v1/scores_main.tsv.gz`): 9600 loci matched; block PP4 max |diff| = 0; NA pattern identical: True; whole-region PP4 max |diff| = 0.

**Internal check**: whole-region PP4, arm whole vs arm est, max |diff| = 0 (should be ~0: same sdY); own sdY estimate vs coloc:::sdY.est, max |diff| = 0; loci with NA sdy_check: 0.

## Matched 5% FPR, test windows

| arm | method | cut-off | detection | FPR (distinct) | gwas_only | block-coloc − coloc.abf (pp), McNemar p | window mean (t 95% CI), windows > 0 |
|---|---|---|---|---|---|---|---|
| est | block-coloc | 0.1448 | 76.8 | 3.8 | 4.8 |  |  |
| est | coloc.abf | 0.1412 | 73.2 | 3.5 | 5.0 | +3.5, p = 2.7e-17 | +3.5 (+2.2, +4.9), 6/6 |
| whole | block-coloc | 0.1448 | 76.8 | 3.8 | 4.8 |  |  |
| whole | coloc.abf | 0.1412 | 73.2 | 3.5 | 5.0 | +3.5, p = 2.7e-17 | +3.5 (+2.2, +4.9), 6/6 |
| one | block-coloc | 0.1462 | 76.7 | 3.8 | 4.8 |  |  |
| one | coloc.abf | 0.1404 | 73.3 | 3.5 | 5.2 | +3.4, p = 2.5e-15 | +3.4 (+1.9, +4.8), 6/6 |

Block-coloc, arm vs arm est (same positive test loci):

| arm | diff vs est (pp), McNemar p | window mean (t 95% CI) |
|---|---|---|
| whole | +0.0, p = 1 | +0.0 (+0.0, +0.0), 0/6 |
| one | -0.1, p = 0.5 | -0.1 (-0.2, +0.1), 0/6 |

## By scenario (block-coloc / coloc.abf, %)

| arm | single_shared | single_double_shared | double_shared2 | double_shared1 | single_distinct | double_distinct |
|---|---|---|---|---|---|---|
| est | 96.2 / 96.3 | 74.2 / 74.3 | 81.8 / 76.8 | 54.8 / 45.3 | 3.7 / 3.5 | 3.8 / 3.5 |
| whole | 96.2 / 96.3 | 74.2 / 74.3 | 81.8 / 76.8 | 54.8 / 45.3 | 3.7 / 3.5 | 3.8 / 3.5 |
| one | 96.2 / 96.3 | 74.0 / 74.3 | 81.8 / 77.0 | 54.7 / 45.5 | 3.7 / 3.5 | 3.8 / 3.5 |

## Fixed rule (PP4 > 0.7 colocalised; PP3 > 0.7 independent), test windows, block-coloc

| arm | colocalised: called coloc / indep | single_distinct: indep (coloc) | double_distinct: indep (coloc) | gwas_only: coloc |
|---|---|---|---|---|
| est | 44.1 / 2.7 | 35.5 (0.2) | 9.8 (0.2) | 0.2 |
| whole | 44.1 / 2.7 | 35.5 (0.2) | 9.8 (0.2) | 0.2 |
| one | 44.1 / 2.7 | 36.0 (0.2) | 9.8 (0.2) | 0.2 |

## Diagnostics (all windows)

| quantity | median / p90 / p99 |
|---|---|
| sdY trait 1, whole region | 1.002 / 1.009 / 1.014 |
| sdY trait 2, whole region | 1.002 / 1.031 / 1.054 |
| sdY trait 1, best block / whole region (block branch, arm est) | 0.9992 / 1.003 / 1.007 |
| sdY trait 2, best block / whole region (block branch, arm est) | 0.9956 / 1.013 / 1.027 |
| abs(log sdY2 ratio), block branch | 0.009306 / 0.02346 / 0.03348 |
| block-coloc abs(PP4 whole − est), all loci | 0 / 6.377e-05 / 0.001723 |
| coloc.abf abs(PP4 whole − est), all loci | 0 / 0 / 0 |
| best block differs from arm est (whole), block branch | 0.0% |
| block-coloc abs(PP4 one − est), all loci | 0.000293 / 0.001559 / 0.005449 |
| coloc.abf abs(PP4 one − est), all loci | 0.0002962 / 0.001554 / 0.005472 |
| best block differs from arm est (one), block branch | 0.0% |
| block subsets where the sdY estimate is NA (negative coefficient), block branch | 0.0% |
