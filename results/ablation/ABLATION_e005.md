# block-coloc construction ablation: h2_gwas = 0.003, h2_eqtl = 0.05

4800 loci, 9 configurations. Same loci for every configuration; each configuration gets its own cut-off at 5% FPR on the calibration windows.

**Reproduction check** (config base vs `outputs/bulk_h2ref_v1/scores_grid.tsv.gz`): 4800 loci matched; block PP4 max |diff| = 0; NA pattern identical: True; whole-region PP4 max |diff| = 0.

## Matched 5% FPR, test windows

| configuration | cut-off | detection | FPR (distinct) | gwas_only | block branch (positives) | best block has shared variant | diff vs base (pp), McNemar p | window mean (t 95% CI), windows > 0 |
|---|---|---|---|---|---|---|---|---|
| base: stage-1 rule (clump 0.1, neighbourhood 0.4, sign split, gate) | 0.014 | 97.7 | 6.8 | 100.0 | 73.0% | 99.8% | — | — |
| nosplit: no sign split | 0.014 | 97.7 | 6.8 | 100.0 | 73.0% | 99.8% | +0.0, p = 1 | +0.0 (+0.0, +0.0), 0/6 |
| rb0.2: neighbourhood r2 >= 0.2 | 0.027 | 97.4 | 6.0 | 92.7 | 74.6% | 99.9% | -0.2, p = 0.58 | -0.3 (-1.1, +0.6), 2/6 |
| rb0.3: neighbourhood r2 >= 0.3 | 0.017 | 97.7 | 6.5 | 100.0 | 73.8% | 99.9% | +0.0, p = 1 | +0.0 (-0.5, +0.5), 2/6 |
| rb0.5: neighbourhood r2 >= 0.5 | 0.013 | 97.2 | 6.3 | 100.0 | 71.2% | 99.8% | -0.4, p = 0.12 | -0.4 (-0.8, -0.0), 0/6 |
| rb0.6: neighbourhood r2 >= 0.6 | 0.013 | 97.1 | 6.2 | 100.0 | 68.8% | 99.4% | -0.6, p = 0.065 | -0.6 (-1.4, +0.2), 0/6 |
| rc0.05: clump r2 < 0.05 | 0.016 | 97.3 | 6.3 | 100.0 | 71.6% | 99.8% | -0.3, p = 0.12 | -0.3 (-0.9, +0.2), 0/6 |
| rc0.2: clump r2 < 0.2 | 0.016 | 97.8 | 6.2 | 100.0 | 73.3% | 99.7% | +0.1, p = 1 | +0.1 (-0.1, +0.3), 1/6 |
| nogate: no gate (block whenever both traits have a lead) | 0.014 | 97.7 | 6.8 | 100.0 | 94.8% | 99.7% | +0.0, p = 1 | +0.0 (+0.0, +0.0), 0/6 |

## By scenario (detection or false-positive rate, %)

| configuration | single_shared | single_double_shared | double_shared2 | double_shared1 | single_distinct | double_distinct |
|---|---|---|---|---|---|---|
| base | 100.0 | 97.3 | 98.0 | 95.3 | 1.7 | 12.0 |
| nosplit | 100.0 | 97.3 | 98.0 | 95.3 | 1.7 | 12.0 |
| rb0.2 | 100.0 | 96.3 | 98.7 | 94.7 | 1.7 | 10.3 |
| rb0.3 | 100.0 | 97.3 | 98.7 | 94.7 | 1.7 | 11.3 |
| rb0.5 | 100.0 | 97.0 | 98.3 | 93.7 | 1.3 | 11.3 |
| rb0.6 | 100.0 | 97.0 | 98.7 | 92.7 | 1.3 | 11.0 |
| rc0.05 | 100.0 | 97.3 | 98.0 | 94.0 | 1.7 | 11.0 |
| rc0.2 | 100.0 | 97.3 | 98.3 | 95.3 | 1.7 | 10.7 |
| nogate | 100.0 | 97.3 | 98.0 | 95.3 | 1.7 | 12.0 |

## Fixed rule (PP4 > 0.7 colocalised; PP3 > 0.7 independent), test windows

| configuration | colocalised: called coloc / indep | single_distinct: indep (coloc) | double_distinct: indep (coloc) |
|---|---|---|---|
| base | 92.3 / 4.8 | 98.3 (0.3) | 88.0 (0.0) |
| nosplit | 92.2 / 4.9 | 98.3 (0.3) | 88.0 (0.0) |
| rb0.2 | 92.2 / 4.4 | 98.0 (0.3) | 81.0 (0.0) |
| rb0.3 | 92.4 / 4.5 | 98.3 (0.3) | 86.7 (0.0) |
| rb0.5 | 91.7 / 5.7 | 98.7 (0.3) | 88.7 (0.0) |
| rb0.6 | 91.0 / 5.8 | 98.7 (0.3) | 89.0 (0.0) |
| rc0.05 | 91.6 / 5.2 | 98.3 (0.3) | 88.0 (0.0) |
| rc0.2 | 92.6 / 4.4 | 98.3 (0.3) | 88.3 (0.0) |
| nogate | 92.3 / 4.8 | 98.3 (0.3) | 88.0 (0.0) |

## Leads and blocks (all windows)

| configuration | single-signal loci with >= 2 leads in a trait | two-causal trait with >= 2 leads (trait 1 / trait 2) | median blocks (block branch) | median best-block size |
|---|---|---|---|---|
| base | 18.2% | 93.6% / 89.3% | 1 | 32 |
| nosplit | 18.2% | 93.6% / 89.3% | 1 | 32 |
| rb0.2 | 18.2% | 93.6% / 89.3% | 1 | 71 |
| rb0.3 | 18.2% | 93.6% / 89.3% | 1 | 47 |
| rb0.5 | 18.2% | 93.6% / 89.3% | 1 | 26 |
| rb0.6 | 18.2% | 93.6% / 89.3% | 1 | 21 |
| rc0.05 | 17.7% | 93.3% / 88.0% | 1 | 32 |
| rc0.2 | 18.8% | 94.1% / 89.6% | 1 | 32 |
| nogate | 18.2% | 93.6% / 89.3% | 1 | 33 |
