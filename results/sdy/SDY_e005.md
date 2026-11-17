# sdY sensitivity: h2_gwas = 0.003, h2_eqtl = 0.05

4800 loci, arms: `est` = sdY estimated on each subset (current); `whole` = sdY estimated once on the whole region; `one` = sdY = 1 (true value).
Leads and blocks are the same in every arm; each method and arm gets its own cut-off at 5% FPR on the calibration windows.

**Reproduction check** (arm est vs `outputs/bulk_h2ref_v1/scores_grid.tsv.gz`): 4800 loci matched; block PP4 max |diff| = 0; NA pattern identical: True; whole-region PP4 max |diff| = 0.

**Internal check**: whole-region PP4, arm whole vs arm est, max |diff| = 0 (should be ~0: same sdY); own sdY estimate vs coloc:::sdY.est, max |diff| = 0; loci with NA sdy_check: 0.

## Matched 5% FPR, test windows

| arm | method | cut-off | detection | FPR (distinct) | gwas_only | block-coloc − coloc.abf (pp), McNemar p | window mean (t 95% CI), windows > 0 |
|---|---|---|---|---|---|---|---|
| est | block-coloc | 0.0140 | 97.7 | 6.8 | 100.0 |  |  |
| est | coloc.abf | 0.0116 | 89.8 | 5.0 | 100.0 | +7.8, p = 3.1e-26 | +7.8 (+5.7, +9.9), 6/6 |
| whole | block-coloc | 0.0140 | 97.7 | 6.8 | 100.0 |  |  |
| whole | coloc.abf | 0.0116 | 89.8 | 5.0 | 100.0 | +7.8, p = 3.1e-26 | +7.8 (+5.7, +9.9), 6/6 |
| one | block-coloc | 0.0135 | 97.7 | 6.8 | 100.0 |  |  |
| one | coloc.abf | 0.0118 | 89.8 | 5.0 | 100.0 | +7.8, p = 3.1e-26 | +7.8 (+5.7, +9.9), 6/6 |

Block-coloc, arm vs arm est (same positive test loci):

| arm | diff vs est (pp), McNemar p | window mean (t 95% CI) |
|---|---|---|
| whole | +0.0, p = 1 | +0.0 (+0.0, +0.0), 0/6 |
| one | +0.0, p = 1 | +0.0 (+0.0, +0.0), 0/6 |

## By scenario (block-coloc / coloc.abf, %)

| arm | single_shared | single_double_shared | double_shared2 | double_shared1 | single_distinct | double_distinct |
|---|---|---|---|---|---|---|
| est | 100.0 / 100.0 | 97.3 / 90.3 | 98.0 / 94.7 | 95.3 / 74.3 | 1.7 / 0.0 | 12.0 / 10.0 |
| whole | 100.0 / 100.0 | 97.3 / 90.3 | 98.0 / 94.7 | 95.3 / 74.3 | 1.7 / 0.0 | 12.0 / 10.0 |
| one | 100.0 / 100.0 | 97.3 / 90.3 | 98.0 / 94.7 | 95.3 / 74.3 | 1.7 / 0.0 | 12.0 / 10.0 |

## Fixed rule (PP4 > 0.7 colocalised; PP3 > 0.7 independent), test windows, block-coloc

| arm | colocalised: called coloc / indep | single_distinct: indep (coloc) | double_distinct: indep (coloc) | gwas_only: coloc |
|---|---|---|---|---|
| est | 92.3 / 4.8 | 98.3 (0.3) | 88.0 (0.0) | 0.0 |
| whole | 92.3 / 4.8 | 98.3 (0.3) | 88.0 (0.0) | 0.0 |
| one | 92.3 / 4.8 | 98.3 (0.3) | 88.0 (0.0) | 0.0 |

## Diagnostics (all windows)

| quantity | median / p90 / p99 |
|---|---|
| sdY trait 1, whole region | 1.002 / 1.009 / 1.014 |
| sdY trait 2, whole region | 1.001 / 1.03 / 1.054 |
| sdY trait 1, best block / whole region (block branch, arm est) | 0.9991 / 1.003 / 1.008 |
| sdY trait 2, best block / whole region (block branch, arm est) | 0.9922 / 1.01 / 1.025 |
| abs(log sdY2 ratio), block branch | 0.01071 / 0.02522 / 0.03957 |
| block-coloc abs(PP4 whole − est), all loci | 0 / 7.311e-05 / 0.001079 |
| coloc.abf abs(PP4 whole − est), all loci | 0 / 0 / 0 |
| best block differs from arm est (whole), block branch | 0.1% |
| block-coloc abs(PP4 one − est), all loci | 1.702e-05 / 0.0006161 / 0.003342 |
| coloc.abf abs(PP4 one − est), all loci | 2.25e-05 / 0.001016 / 0.009817 |
| best block differs from arm est (one), block branch | 0.1% |
| block subsets where the sdY estimate is NA (negative coefficient), block branch | 0.0% |
