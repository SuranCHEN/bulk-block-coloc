# block-coloc construction ablation: h2_gwas = 0.003, h2_eqtl = 0.015

9600 loci, 9 configurations. Same loci for every configuration; each configuration gets its own cut-off at 5% FPR on the calibration windows.

**Reproduction check** (config base vs `outputs/bulk_h2ref_v1/scores_main.tsv.gz`): 9600 loci matched; block PP4 max |diff| = 0; NA pattern identical: True; whole-region PP4 max |diff| = 0.

## Matched 5% FPR, test windows

| configuration | cut-off | detection | FPR (distinct) | gwas_only | block branch (positives) | best block has shared variant | diff vs base (pp), McNemar p | window mean (t 95% CI), windows > 0 |
|---|---|---|---|---|---|---|---|---|
| base: stage-1 rule (clump 0.1, neighbourhood 0.4, sign split, gate) | 0.145 | 76.8 | 3.8 | 4.8 | 29.2% | 98.9% | — | — |
| nosplit: no sign split | 0.145 | 76.8 | 3.8 | 4.8 | 29.2% | 98.9% | +0.0, p = 1 | +0.0 (+0.0, +0.0), 0/6 |
| rb0.2: neighbourhood r2 >= 0.2 | 0.144 | 77.3 | 4.2 | 4.8 | 33.3% | 97.2% | +0.6, p = 0.016 | +0.6 (-0.3, +1.5), 4/6 |
| rb0.3: neighbourhood r2 >= 0.3 | 0.142 | 77.3 | 4.0 | 5.0 | 30.6% | 98.4% | +0.6, p = 0.00052 | +0.6 (+0.3, +0.9), 6/6 |
| rb0.5: neighbourhood r2 >= 0.5 | 0.147 | 76.3 | 3.8 | 4.8 | 27.4% | 98.9% | -0.5, p = 0.0034 | -0.5 (-0.8, -0.2), 0/6 |
| rb0.6: neighbourhood r2 >= 0.6 | 0.147 | 75.9 | 3.7 | 4.8 | 25.0% | 98.7% | -0.8, p = 0.00018 | -0.8 (-1.6, -0.1), 0/6 |
| rc0.05: clump r2 < 0.05 | 0.145 | 76.8 | 3.8 | 4.8 | 28.6% | 99.0% | +0.0, p = 1 | +0.0 (-0.3, +0.3), 2/6 |
| rc0.2: clump r2 < 0.2 | 0.144 | 76.9 | 3.8 | 4.8 | 29.5% | 99.0% | +0.2, p = 0.12 | +0.2 (+0.0, +0.3), 4/6 |
| nogate: no gate (block whenever both traits have a lead) | 0.145 | 76.8 | 3.8 | 4.8 | 50.5% | 99.3% | +0.0, p = 1 | +0.0 (+0.0, +0.0), 0/6 |

## By scenario (detection or false-positive rate, %)

| configuration | single_shared | single_double_shared | double_shared2 | double_shared1 | single_distinct | double_distinct |
|---|---|---|---|---|---|---|
| base | 96.2 | 74.2 | 81.8 | 54.8 | 3.7 | 3.8 |
| nosplit | 96.2 | 74.2 | 81.8 | 54.8 | 3.7 | 3.8 |
| rb0.2 | 96.2 | 74.7 | 82.0 | 56.5 | 3.7 | 4.8 |
| rb0.3 | 96.3 | 74.7 | 82.2 | 56.2 | 3.7 | 4.3 |
| rb0.5 | 96.2 | 73.8 | 81.7 | 53.5 | 3.7 | 3.8 |
| rb0.6 | 96.2 | 73.7 | 81.3 | 52.5 | 3.7 | 3.7 |
| rc0.05 | 96.2 | 74.3 | 81.8 | 54.7 | 3.7 | 3.8 |
| rc0.2 | 96.2 | 74.3 | 82.0 | 55.2 | 3.7 | 4.0 |
| nogate | 96.2 | 74.2 | 81.8 | 54.8 | 3.7 | 3.8 |

## Fixed rule (PP4 > 0.7 colocalised; PP3 > 0.7 independent), test windows

| configuration | colocalised: called coloc / indep | single_distinct: indep (coloc) | double_distinct: indep (coloc) |
|---|---|---|---|
| base | 44.1 / 2.7 | 35.5 (0.2) | 9.8 (0.2) |
| nosplit | 44.1 / 2.7 | 35.5 (0.2) | 9.8 (0.2) |
| rb0.2 | 44.2 / 2.5 | 35.0 (0.2) | 8.0 (0.2) |
| rb0.3 | 44.1 / 2.7 | 35.5 (0.2) | 9.3 (0.2) |
| rb0.5 | 43.7 / 2.8 | 35.5 (0.2) | 9.8 (0.2) |
| rb0.6 | 43.2 / 2.9 | 35.5 (0.2) | 10.0 (0.2) |
| rc0.05 | 43.9 / 2.7 | 35.5 (0.2) | 9.8 (0.2) |
| rc0.2 | 44.2 / 2.6 | 35.5 (0.2) | 9.8 (0.2) |
| nogate | 44.3 / 2.7 | 35.5 (0.2) | 9.8 (0.2) |

## Leads and blocks (all windows)

| configuration | single-signal loci with >= 2 leads in a trait | two-causal trait with >= 2 leads (trait 1 / trait 2) | median blocks (block branch) | median best-block size |
|---|---|---|---|---|
| base | 17.4% | 93.7% / 16.4% | 1 | 42 |
| nosplit | 17.4% | 93.7% / 16.4% | 1 | 42 |
| rb0.2 | 17.4% | 93.7% / 16.4% | 1 | 93 |
| rb0.3 | 17.4% | 93.7% / 16.4% | 1 | 59 |
| rb0.5 | 17.4% | 93.7% / 16.4% | 1 | 34 |
| rb0.6 | 17.4% | 93.7% / 16.4% | 1 | 27 |
| rc0.05 | 17.0% | 93.0% / 16.2% | 1 | 42 |
| rc0.2 | 17.8% | 94.2% / 16.5% | 1 | 43 |
| nogate | 17.4% | 93.7% / 16.4% | 1 | 46 |
