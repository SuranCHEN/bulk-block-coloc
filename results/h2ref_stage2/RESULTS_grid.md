# Bulk block-coloc: stage 2, one factor at a time (test windows, frozen cut-offs)

Generated 2026-10-29 03:30. Replicates per cell: 50. Cut-offs and lead thresholds from stage 1 (thresholds.tsv), not re-tuned.
Strength = local h2 of each trait; each trait's h2 split equally over its causal SNPs (per-SNP h2 = h2 / K). LD given to the methods: 1000 Genomes EUR reference LD (489 people) for both traits.

## Three-line conclusion

1. The gain of block-coloc over whole-region coloc.abf grows with eQTL strength, as predicted by the lead-count calculation (Supplementary Note, Proposition 5): block-coloc (whole-region fallback when no block contains leads of both traits, the adopted rule) detects +4.0 points more than coloc.abf at eQTL h2 1.5%, +7.1 at 3% and +13.7 at 5% (95.9% vs 82.2%; all McNemar p < 1e-9). The share of two-causal eQTL loci with two conditional leads was 4.0 / 15.9 / 55.9 / 89.3% (predicted 4.0 / 14.9 / 56.2 / 90.2%). coloc.abf itself loses power at the strongest eQTL (83.1% to 82.2%) because two eQTL signals violate its single-causal-variant assumption.
2. coloc.signals (cond) detects 1.9-3.8 points more than block-coloc once either trait is strong, but with the frozen stage-1 cut-offs its false-positive rate rises to 6.0% (eQTL h2 5%), 6.2% and 8.7% (GWAS h2 1% and 2%), while block-coloc stays at or below 4.3%, so these cells are not compared at matched error. With weak signals (GWAS h2 0.1% or eQTL h2 0.75%) block-coloc detects 2.9-4.7 points more than coloc.signals (cond).
3. When the two causal variants of a trait are in LD (|r| >= 0.32), all methods lose power and every coloc.abf-based method calls 19-29% (mid tier) and 70-79% (high tier) of single_distinct loci colocalised: distinct causal variants in strong LD are not separable at these sample sizes, and block-coloc rarely splits (double_shared2: 70% -> 38% -> 17% of loci), as expected from Proposition 3. gwas_only false positives of block-coloc are 3.3-5.3%, the same range as coloc.abf; the version without fallback reaches 8.7-9.3% because of its lower cut-off.

## Detection rate vs GWAS h2 (eQTL h2 1.5%, low LD)

| h2_gwas | block-coloc, conditional (proposed) | block-coloc (no fallback) | block-coloc, marginal | block-coloc, handoff | coloc.abf | coloc.susie | coloc.susie + fallback | signals (mask) | signals (cond) |
|---|---|---|---|---|---|---|---|---|---|
| 0.001 | 53.0% | 52.9% | 52.1% | 52.9% | 52.2% | 10.7% | 51.6% | 50.7% | 48.3% |
| 0.003 | 76.4% | 72.5% | 72.8% | 76.5% | 72.4% | 18.4% | 71.9% | 77.6% | 75.7% |
| 0.010 | 76.9% | 72.7% | 71.0% | 76.5% | 70.8% | 19.5% | 71.7% | 79.8% | 79.3% |
| 0.020 | 76.4% | 72.2% | 69.8% | 78.2% | 69.8% | 19.3% | 70.9% | 79.7% | 80.2% |

## False-positive rate vs GWAS h2

| h2_gwas | block-coloc, conditional (proposed) | block-coloc (no fallback) | block-coloc, marginal | block-coloc, handoff | coloc.abf | coloc.susie | coloc.susie + fallback | signals (mask) | signals (cond) |
|---|---|---|---|---|---|---|---|---|---|
| 0.001 | 2.7% | 2.7% | 2.5% | 2.5% | 2.5% | 3.2% | 2.5% | 2.5% | 1.8% |
| 0.003 | 3.2% | 3.8% | 3.0% | 3.3% | 3.0% | 3.3% | 2.8% | 5.0% | 4.2% |
| 0.010 | 3.7% | 4.2% | 3.8% | 8.0% | 3.0% | 4.5% | 3.0% | 5.7% | 6.2% |
| 0.020 | 4.3% | 4.5% | 4.7% | 13.0% | 3.0% | 4.7% | 3.0% | 6.5% | 8.7% |

## Detection rate vs eQTL h2 (GWAS h2 0.3%, low LD)

| h2_eqtl | block-coloc, conditional (proposed) | block-coloc (no fallback) | block-coloc, marginal | block-coloc, handoff | coloc.abf | coloc.susie | coloc.susie + fallback | signals (mask) | signals (cond) |
|---|---|---|---|---|---|---|---|---|---|
| 0.0075 | 52.9% | 54.7% | 51.7% | 54.1% | 51.4% | 4.6% | 51.9% | 54.1% | 50.0% |
| 0.0150 | 76.4% | 72.5% | 72.8% | 76.5% | 72.4% | 18.4% | 71.9% | 77.6% | 75.7% |
| 0.0300 | 90.2% | 86.2% | 83.2% | 89.2% | 83.1% | 55.6% | 80.8% | 90.7% | 93.1% |
| 0.0500 | 95.9% | 94.3% | 89.1% | 94.8% | 82.2% | 81.9% | 85.6% | 94.2% | 97.8% |

## False-positive rate vs eQTL h2

| h2_eqtl | block-coloc, conditional (proposed) | block-coloc (no fallback) | block-coloc, marginal | block-coloc, handoff | coloc.abf | coloc.susie | coloc.susie + fallback | signals (mask) | signals (cond) |
|---|---|---|---|---|---|---|---|---|---|
| 0.0075 | 3.7% | 6.3% | 3.8% | 5.2% | 4.0% | 1.2% | 4.0% | 5.3% | 4.2% |
| 0.0150 | 3.2% | 3.8% | 3.0% | 3.3% | 3.0% | 3.3% | 2.8% | 5.0% | 4.2% |
| 0.0300 | 2.0% | 1.5% | 1.0% | 1.7% | 1.2% | 7.0% | 1.2% | 3.7% | 4.3% |
| 0.0500 | 1.3% | 1.3% | 0.8% | 1.8% | 0.3% | 5.5% | 0.7% | 3.5% | 6.0% |

## False-positive rate when at most one trait has a signal (extra negatives, not used for calibration)

gwas_only: GWAS has one causal variant, the eQTL has none. null_both: neither trait has a signal.

### vs GWAS h2 (eQTL h2 1.5%)

| scenario | h2_gwas | n | block-coloc, conditional (proposed) | block-coloc (no fallback) | block-coloc, marginal | block-coloc, handoff | coloc.abf | coloc.susie | coloc.susie + fallback | signals (mask) | signals (cond) |
|---|---|---|---|---|---|---|---|---|---|---|---|
| gwas_only | 0.001 | 300 | 3.3% | 4.7% | 3.3% | 3.3% | 3.3% | 0.0% | 3.3% | 2.3% | 1.3% |
| gwas_only | 0.003 | 300 | 4.7% | 8.7% | 5.0% | 5.3% | 5.0% | 0.0% | 5.3% | 3.7% | 2.7% |
| gwas_only | 0.010 | 300 | 5.3% | 9.3% | 5.3% | 9.3% | 5.3% | 0.0% | 5.3% | 5.0% | 2.7% |
| gwas_only | 0.020 | 300 | 5.3% | 9.0% | 5.7% | 12.7% | 5.3% | 0.0% | 5.7% | 5.3% | 2.7% |
| null_both | 0.001 | 300 | 0.0% | 0.0% | 0.0% | 0.0% | 0.0% | 0.0% | 0.0% | 0.0% | 0.0% |
| null_both | 0.003 | 300 | 0.0% | 0.0% | 0.0% | 0.0% | 0.0% | 0.0% | 0.0% | 0.0% | 0.0% |
| null_both | 0.010 | 300 | 0.0% | 0.0% | 0.0% | 0.0% | 0.0% | 0.0% | 0.0% | 0.0% | 0.0% |
| null_both | 0.020 | 300 | 0.0% | 0.0% | 0.0% | 0.0% | 0.0% | 0.0% | 0.0% | 0.0% | 0.0% |

### vs eQTL h2 (GWAS h2 0.3%)

| scenario | h2_eqtl | n | block-coloc, conditional (proposed) | block-coloc (no fallback) | block-coloc, marginal | block-coloc, handoff | coloc.abf | coloc.susie | coloc.susie + fallback | signals (mask) | signals (cond) |
|---|---|---|---|---|---|---|---|---|---|---|---|
| gwas_only | 0.0075 | 300 | 4.7% | 8.7% | 5.0% | 5.3% | 5.0% | 0.0% | 5.3% | 3.7% | 2.7% |
| gwas_only | 0.0150 | 300 | 4.7% | 8.7% | 5.0% | 5.3% | 5.0% | 0.0% | 5.3% | 3.7% | 2.7% |
| gwas_only | 0.0300 | 300 | 4.7% | 8.7% | 5.0% | 5.3% | 5.0% | 0.0% | 5.3% | 3.7% | 2.7% |
| gwas_only | 0.0500 | 300 | 4.7% | 8.7% | 5.0% | 5.3% | 5.0% | 0.0% | 5.3% | 3.7% | 2.7% |
| null_both | 0.0075 | 300 | 0.0% | 0.0% | 0.0% | 0.0% | 0.0% | 0.0% | 0.0% | 0.0% | 0.0% |
| null_both | 0.0150 | 300 | 0.0% | 0.0% | 0.0% | 0.0% | 0.0% | 0.0% | 0.0% | 0.0% | 0.0% |
| null_both | 0.0300 | 300 | 0.0% | 0.0% | 0.0% | 0.0% | 0.0% | 0.0% | 0.0% | 0.0% | 0.0% |
| null_both | 0.0500 | 300 | 0.0% | 0.0% | 0.0% | 0.0% | 0.0% | 0.0% | 0.0% | 0.0% | 0.0% |

## LD between causal SNPs (main strengths; single_double_shared, double_shared2, double_shared1 pooled for TPR; single_distinct, double_distinct for FPR)

low |r| < 0.2; mid 0.32-0.63; high 0.63-0.9.

| tier | block-coloc, conditional (proposed) | block-coloc (no fallback) | block-coloc, marginal | block-coloc, handoff | coloc.abf | coloc.susie | coloc.susie + fallback | signals (mask) | signals (cond) |
|---|---|---|---|---|---|---|---|---|---|
| low | 69.9% | 64.8% | 65.0% | 70.1% | 64.6% | 12.7% | 64.1% | 71.6% | 69.2% |
| mid | 60.8% | 57.9% | 59.9% | 61.8% | 60.6% | 12.3% | 60.4% | 56.2% | 58.6% |
| high | 46.0% | 49.9% | 46.8% | 47.7% | 46.6% | 17.0% | 46.4% | 44.0% | 44.1% |

| tier | block-coloc, conditional (proposed) | block-coloc (no fallback) | block-coloc, marginal | block-coloc, handoff | coloc.abf | coloc.susie | coloc.susie + fallback | signals (mask) | signals (cond) |
|---|---|---|---|---|---|---|---|---|---|
| low | 3.2% | 3.8% | 3.0% | 3.3% | 3.0% | 3.3% | 2.8% | 5.0% | 4.2% |
| mid | 14.8% | 16.7% | 14.5% | 16.5% | 14.2% | 9.8% | 13.7% | 12.0% | 11.0% |
| high | 39.2% | 41.3% | 39.2% | 39.3% | 39.0% | 16.5% | 38.5% | 36.8% | 35.5% |

### By scenario

| tier | scenario | n | block-coloc, conditional (proposed) | block-coloc (no fallback) | block-coloc, marginal | block-coloc, handoff | coloc.abf | coloc.susie | coloc.susie + fallback | signals (mask) | signals (cond) |
|---|---|---|---|---|---|---|---|---|---|---|---|
| low | single_double_shared | 300 | 72.7% | 72.0% | 74.0% | 71.3% | 73.3% | 12.0% | 71.7% | 71.3% | 66.0% |
| low | double_shared2 | 300 | 84.7% | 79.3% | 79.0% | 84.3% | 78.0% | 15.3% | 79.3% | 85.7% | 85.3% |
| low | double_shared1 | 300 | 52.3% | 43.0% | 42.0% | 54.7% | 42.3% | 10.7% | 41.3% | 57.7% | 56.3% |
| low | single_distinct | 300 | 2.7% | 4.7% | 2.3% | 2.7% | 2.3% | 4.3% | 2.0% | 3.0% | 3.3% |
| low | double_distinct | 300 | 3.7% | 3.0% | 3.7% | 4.0% | 3.7% | 2.3% | 3.7% | 7.0% | 5.0% |
| mid | single_double_shared | 300 | 64.3% | 65.3% | 65.0% | 64.3% | 65.0% | 10.7% | 64.7% | 61.0% | 56.7% |
| mid | double_shared2 | 300 | 68.0% | 63.3% | 66.0% | 68.3% | 67.3% | 16.7% | 67.3% | 63.3% | 70.3% |
| mid | double_shared1 | 300 | 50.0% | 45.0% | 48.7% | 52.7% | 49.3% | 9.7% | 49.3% | 44.3% | 48.7% |
| mid | single_distinct | 300 | 26.0% | 28.7% | 25.7% | 29.0% | 25.0% | 17.0% | 23.7% | 21.3% | 19.3% |
| mid | double_distinct | 300 | 3.7% | 4.7% | 3.3% | 4.0% | 3.3% | 2.7% | 3.7% | 2.7% | 2.7% |
| high | single_double_shared | 300 | 59.0% | 64.3% | 60.0% | 59.3% | 59.7% | 23.7% | 59.7% | 56.0% | 55.0% |
| high | double_shared2 | 300 | 46.0% | 49.0% | 46.3% | 47.7% | 46.3% | 15.0% | 45.7% | 44.7% | 46.0% |
| high | double_shared1 | 300 | 33.0% | 36.3% | 34.0% | 36.0% | 33.7% | 12.3% | 34.0% | 31.3% | 31.3% |
| high | single_distinct | 300 | 75.7% | 78.7% | 75.7% | 75.7% | 75.3% | 31.3% | 74.3% | 72.3% | 70.3% |
| high | double_distinct | 300 | 2.7% | 4.0% | 2.7% | 3.0% | 2.7% | 1.7% | 2.7% | 1.3% | 0.7% |

### Branch choice by LD tier

| tier | method | scenario | true_signals | should_go | n | branch_block | correct_branch | leads_match_truth | median_leads |
|---|---|---|---|---|---|---|---|---|---|
| low | block-coloc (conditional lead selection, proposed) | single_double_shared | 1+2 | block | 300 | 23.7% | 23.7% | 15.0% | 1+1 |
| low | block-coloc (conditional lead selection, proposed) | double_shared2 | 2+2 | block | 300 | 70.3% | 70.3% | 13.3% | 2+1 |
| low | block-coloc (conditional lead selection, proposed) | double_shared1 | 2+2 | block | 300 | 68.7% | 68.7% | 10.7% | 2+1 |
| low | block-coloc (conditional lead selection, proposed) | single_distinct | 1+1 | whole | 300 | 19.0% | 81.0% | 65.7% | 1+1 |
| low | block-coloc (conditional lead selection, proposed) | double_distinct | 2+2 | block | 300 | 64.0% | 64.0% | 13.7% | 2+1 |
| low | block-coloc (no fallback) | single_double_shared | 1+2 | block | 300 | 23.7% | 23.7% | 15.0% | 1+1 |
| low | block-coloc (no fallback) | double_shared2 | 2+2 | block | 300 | 70.3% | 70.3% | 13.3% | 2+1 |
| low | block-coloc (no fallback) | double_shared1 | 2+2 | block | 300 | 68.7% | 68.7% | 10.7% | 2+1 |
| low | block-coloc (no fallback) | single_distinct | 1+1 | whole | 300 | 19.0% | 81.0% | 65.7% | 1+1 |
| low | block-coloc (no fallback) | double_distinct | 2+2 | block | 300 | 64.0% | 64.0% | 13.7% | 2+1 |
| low | block-coloc (marginal lead selection) | single_double_shared | 1+2 | block | 300 | 1.7% | 1.7% | 0.7% | 1+0 |
| low | block-coloc (marginal lead selection) | double_shared2 | 2+2 | block | 300 | 14.7% | 14.7% | 0.3% | 2+0 |
| low | block-coloc (marginal lead selection) | double_shared1 | 2+2 | block | 300 | 16.3% | 16.3% | 0.3% | 2+0 |
| low | block-coloc (marginal lead selection) | single_distinct | 1+1 | whole | 300 | 5.7% | 94.3% | 30.7% | 1+0 |
| low | block-coloc (marginal lead selection) | double_distinct | 2+2 | block | 300 | 7.3% | 7.3% | 0.7% | 2+0 |
| low | block-coloc (original handoff) | single_double_shared | 1+2 | block | 300 | 48.3% | 48.3% | 18.3% | 1+1 |
| low | block-coloc (original handoff) | double_shared2 | 2+2 | block | 300 | 82.0% | 82.0% | 15.7% | 2+1 |
| low | block-coloc (original handoff) | double_shared1 | 2+2 | block | 300 | 55.7% | 55.7% | 17.7% | 2+1 |
| low | block-coloc (original handoff) | single_distinct | 1+1 | whole | 300 | 7.0% | 93.0% | 57.7% | 1+1 |
| low | block-coloc (original handoff) | double_distinct | 2+2 | block | 300 | 33.0% | 33.0% | 18.3% | 2+1 |
| mid | block-coloc (conditional lead selection, proposed) | single_double_shared | 1+2 | block | 300 | 14.7% | 14.7% | 7.0% | 1+1 |
| mid | block-coloc (conditional lead selection, proposed) | double_shared2 | 2+2 | block | 300 | 38.0% | 38.0% | 2.7% | 1+1 |
| mid | block-coloc (conditional lead selection, proposed) | double_shared1 | 2+2 | block | 300 | 33.3% | 33.3% | 2.7% | 1+1 |
| mid | block-coloc (conditional lead selection, proposed) | single_distinct | 1+1 | whole | 300 | 12.0% | 88.0% | 72.7% | 1+1 |
| mid | block-coloc (conditional lead selection, proposed) | double_distinct | 2+2 | block | 300 | 36.7% | 36.7% | 3.7% | 1+1 |
| mid | block-coloc (no fallback) | single_double_shared | 1+2 | block | 300 | 14.7% | 14.7% | 7.0% | 1+1 |
| mid | block-coloc (no fallback) | double_shared2 | 2+2 | block | 300 | 38.0% | 38.0% | 2.7% | 1+1 |
| mid | block-coloc (no fallback) | double_shared1 | 2+2 | block | 300 | 33.3% | 33.3% | 2.7% | 1+1 |
| mid | block-coloc (no fallback) | single_distinct | 1+1 | whole | 300 | 12.0% | 88.0% | 72.7% | 1+1 |
| mid | block-coloc (no fallback) | double_distinct | 2+2 | block | 300 | 36.7% | 36.7% | 3.7% | 1+1 |
| mid | block-coloc (marginal lead selection) | single_double_shared | 1+2 | block | 300 | 2.3% | 2.3% | 0.0% | 1+0 |
| mid | block-coloc (marginal lead selection) | double_shared2 | 2+2 | block | 300 | 12.3% | 12.3% | 0.0% | 1+0 |
| mid | block-coloc (marginal lead selection) | double_shared1 | 2+2 | block | 300 | 9.3% | 9.3% | 0.3% | 1+0 |
| mid | block-coloc (marginal lead selection) | single_distinct | 1+1 | whole | 300 | 6.7% | 93.3% | 31.3% | 1+0 |
| mid | block-coloc (marginal lead selection) | double_distinct | 2+2 | block | 300 | 11.3% | 11.3% | 0.3% | 1+0 |
| mid | block-coloc (original handoff) | single_double_shared | 1+2 | block | 300 | 51.7% | 51.7% | 12.7% | 1+1 |
| mid | block-coloc (original handoff) | double_shared2 | 2+2 | block | 300 | 69.3% | 69.3% | 11.0% | 2+1 |
| mid | block-coloc (original handoff) | double_shared1 | 2+2 | block | 300 | 48.0% | 48.0% | 7.3% | 2+1 |
| mid | block-coloc (original handoff) | single_distinct | 1+1 | whole | 300 | 34.7% | 65.3% | 59.3% | 1+1 |
| mid | block-coloc (original handoff) | double_distinct | 2+2 | block | 300 | 17.0% | 17.0% | 10.0% | 2+1 |
| high | block-coloc (conditional lead selection, proposed) | single_double_shared | 1+2 | block | 300 | 15.7% | 15.7% | 7.3% | 1+1 |
| high | block-coloc (conditional lead selection, proposed) | double_shared2 | 2+2 | block | 300 | 17.0% | 17.0% | 0.0% | 1+1 |
| high | block-coloc (conditional lead selection, proposed) | double_shared1 | 2+2 | block | 300 | 12.7% | 12.7% | 0.7% | 1+1 |
| high | block-coloc (conditional lead selection, proposed) | single_distinct | 1+1 | whole | 300 | 16.3% | 83.7% | 63.7% | 1+1 |
| high | block-coloc (conditional lead selection, proposed) | double_distinct | 2+2 | block | 300 | 17.0% | 17.0% | 1.0% | 1+1 |
| high | block-coloc (no fallback) | single_double_shared | 1+2 | block | 300 | 15.7% | 15.7% | 7.3% | 1+1 |
| high | block-coloc (no fallback) | double_shared2 | 2+2 | block | 300 | 17.0% | 17.0% | 0.0% | 1+1 |
| high | block-coloc (no fallback) | double_shared1 | 2+2 | block | 300 | 12.7% | 12.7% | 0.7% | 1+1 |
| high | block-coloc (no fallback) | single_distinct | 1+1 | whole | 300 | 16.3% | 83.7% | 63.7% | 1+1 |
| high | block-coloc (no fallback) | double_distinct | 2+2 | block | 300 | 17.0% | 17.0% | 1.0% | 1+1 |
| high | block-coloc (marginal lead selection) | single_double_shared | 1+2 | block | 300 | 4.0% | 4.0% | 0.7% | 1+0 |
| high | block-coloc (marginal lead selection) | double_shared2 | 2+2 | block | 300 | 2.0% | 2.0% | 0.0% | 1+0 |
| high | block-coloc (marginal lead selection) | double_shared1 | 2+2 | block | 300 | 1.7% | 1.7% | 0.0% | 1+0 |
| high | block-coloc (marginal lead selection) | single_distinct | 1+1 | whole | 300 | 4.3% | 95.7% | 32.7% | 1+0 |
| high | block-coloc (marginal lead selection) | double_distinct | 2+2 | block | 300 | 1.7% | 1.7% | 0.0% | 1+0 |
| high | block-coloc (original handoff) | single_double_shared | 1+2 | block | 300 | 60.3% | 60.3% | 11.7% | 1+1 |
| high | block-coloc (original handoff) | double_shared2 | 2+2 | block | 300 | 52.7% | 52.7% | 4.0% | 1+1 |
| high | block-coloc (original handoff) | double_shared1 | 2+2 | block | 300 | 45.0% | 45.0% | 7.3% | 1+1 |
| high | block-coloc (original handoff) | single_distinct | 1+1 | whole | 300 | 89.0% | 11.0% | 56.0% | 1+1 |
| high | block-coloc (original handoff) | double_distinct | 2+2 | block | 300 | 12.0% | 12.0% | 3.7% | 1+1 |

Figure: `grid_lines.png`
