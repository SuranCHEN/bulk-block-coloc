# block-coloc construction ablation: h2_gwas = 0.02, h2_eqtl = 0.015

4800 loci, 9 configurations. Same loci for every configuration; each configuration gets its own cut-off at 5% FPR on the calibration windows.

**Reproduction check** (config base vs `outputs/bulk_h2ref_v1/scores_grid.tsv.gz`): 4800 loci matched; block PP4 max |diff| = 0; NA pattern identical: True; whole-region PP4 max |diff| = 0.

## Matched 5% FPR, test windows

| configuration | cut-off | detection | FPR (distinct) | gwas_only | block branch (positives) | best block has shared variant | diff vs base (pp), McNemar p | window mean (t 95% CI), windows > 0 |
|---|---|---|---|---|---|---|---|---|
| base: stage-1 rule (clump 0.1, neighbourhood 0.4, sign split, gate) | 0.147 | 76.3 | 4.3 | 5.0 | 34.4% | 98.3% | — | — |
| nosplit: no sign split | 0.147 | 76.3 | 4.3 | 5.0 | 34.4% | 98.3% | +0.0, p = 1 | +0.0 (+0.0, +0.0), 0/6 |
| rb0.2: neighbourhood r2 >= 0.2 | 0.149 | 76.8 | 4.3 | 5.0 | 39.0% | 95.9% | +0.5, p = 0.31 | +0.5 (-0.4, +1.4), 3/6 |
| rb0.3: neighbourhood r2 >= 0.3 | 0.149 | 76.6 | 4.2 | 5.0 | 36.2% | 97.9% | +0.2, p = 0.55 | +0.3 (-0.2, +0.7), 2/6 |
| rb0.5: neighbourhood r2 >= 0.5 | 0.144 | 76.1 | 4.0 | 5.3 | 32.5% | 98.7% | -0.2, p = 0.45 | -0.3 (-0.5, +0.0), 0/6 |
| rb0.6: neighbourhood r2 >= 0.6 | 0.144 | 75.2 | 3.5 | 5.3 | 30.1% | 98.3% | -1.2, p = 0.0013 | -1.2 (-1.9, -0.4), 0/6 |
| rc0.05: clump r2 < 0.05 | 0.148 | 76.2 | 3.8 | 5.0 | 33.5% | 99.3% | -0.2, p = 0.73 | -0.2 (-0.6, +0.3), 0/6 |
| rc0.2: clump r2 < 0.2 | 0.147 | 76.4 | 4.0 | 5.0 | 35.9% | 97.9% | +0.1, p = 1 | +0.1 (-0.1, +0.3), 1/6 |
| nogate: no gate (block whenever both traits have a lead) | 0.147 | 76.3 | 4.3 | 5.0 | 52.3% | 98.9% | +0.0, p = 1 | +0.0 (+0.0, +0.0), 0/6 |

## By scenario (detection or false-positive rate, %)

| configuration | single_shared | single_double_shared | double_shared2 | double_shared1 | single_distinct | double_distinct |
|---|---|---|---|---|---|---|
| base | 96.3 | 72.7 | 84.3 | 52.0 | 3.0 | 5.7 |
| nosplit | 96.3 | 72.7 | 84.3 | 52.0 | 3.0 | 5.7 |
| rb0.2 | 96.3 | 72.3 | 85.7 | 53.0 | 3.3 | 5.3 |
| rb0.3 | 96.3 | 72.3 | 85.3 | 52.3 | 3.0 | 5.3 |
| rb0.5 | 96.3 | 73.3 | 83.7 | 51.0 | 3.0 | 5.0 |
| rb0.6 | 96.3 | 73.0 | 83.0 | 48.3 | 2.7 | 4.3 |
| rc0.05 | 96.3 | 72.7 | 84.7 | 51.0 | 3.0 | 4.7 |
| rc0.2 | 96.3 | 72.7 | 85.0 | 51.7 | 3.0 | 5.0 |
| nogate | 96.3 | 72.7 | 84.3 | 52.0 | 3.0 | 5.7 |

## Fixed rule (PP4 > 0.7 colocalised; PP3 > 0.7 independent), test windows

| configuration | colocalised: called coloc / indep | single_distinct: indep (coloc) | double_distinct: indep (coloc) |
|---|---|---|---|
| base | 47.2 / 3.2 | 34.7 (0.0) | 9.0 (0.3) |
| nosplit | 47.2 / 3.2 | 34.7 (0.0) | 9.0 (0.3) |
| rb0.2 | 47.3 / 2.8 | 34.3 (0.3) | 7.0 (0.3) |
| rb0.3 | 47.2 / 3.0 | 34.7 (0.0) | 8.3 (0.3) |
| rb0.5 | 46.8 / 3.3 | 34.7 (0.0) | 9.0 (0.3) |
| rb0.6 | 46.0 / 3.6 | 34.7 (0.0) | 9.0 (0.3) |
| rc0.05 | 46.8 / 3.3 | 34.7 (0.3) | 9.0 (0.0) |
| rc0.2 | 47.5 / 3.2 | 34.7 (0.0) | 9.0 (0.3) |
| nogate | 47.4 / 3.2 | 34.7 (0.0) | 9.0 (0.3) |

## Leads and blocks (all windows)

| configuration | single-signal loci with >= 2 leads in a trait | two-causal trait with >= 2 leads (trait 1 / trait 2) | median blocks (block branch) | median best-block size |
|---|---|---|---|---|
| base | 30.7% | 99.1% / 15.9% | 1 | 46 |
| nosplit | 30.7% | 99.1% / 15.9% | 1 | 46 |
| rb0.2 | 30.7% | 99.1% / 15.9% | 1 | 94 |
| rb0.3 | 30.7% | 99.1% / 15.9% | 1 | 65 |
| rb0.5 | 30.7% | 99.1% / 15.9% | 1 | 36 |
| rb0.6 | 30.7% | 99.1% / 15.9% | 1 | 29 |
| rc0.05 | 29.5% | 98.9% / 15.5% | 1 | 46 |
| rc0.2 | 34.6% | 99.7% / 15.9% | 1 | 46 |
| nogate | 30.7% | 99.1% / 15.9% | 1 | 46 |
