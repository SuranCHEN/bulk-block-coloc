# Ablation: where block-coloc's gain over coloc.abf comes from

Existing per-locus scores, reference LD. Exact decomposition of the detection-rate difference on the same positive test loci (6 test windows) by the branch block-coloc takes. In G0 and G2 block-coloc returns the whole-region PP4, so any difference there comes only from the two methods' different cut-offs.

Check: block-coloc PP4 equals whole-region PP4 outside G1, max |diff| = 0.

## Stage 1 (default setting)

### Default: h2_gwas 0.3%, h2_eqtl 1.5%

Cut-offs (5% FPR on calibration windows): coloc.abf 0.1412, block-coloc 0.1448, block-coloc without fallback 0.1067. Positive test loci: 2400. Detection: coloc.abf 73.2%, block-coloc 76.8% (difference +3.5), without fallback 73.8%.

| branch | share of positives | coloc.abf | block-coloc | contribution to the difference (pp) | of which: PP4 change at common cut-off | of which: cut-off change |
|---|---|---|---|---|---|---|
| G0 gate not passed (PP4 identical) | 56.8% | 74.4 | 74.1 | -0.17 | +0.00 | -0.17 |
| G1 valid block (PP4 = max block) | 29.2% | 84.7 | 97.7 | +3.79 | +3.79 | +0.00 |
| G2 no valid block (PP4 identical, fallback) | 14.1% | 44.7 | 44.1 | -0.08 | +0.00 | -0.08 |
| total | 100% | 73.2 | 76.8 | +3.54 | | |

Negatives (calibration windows) above coloc.abf's cut-off, by branch — this sets block-coloc's own cut-off:

- G0: share 59.3%, above t_abf coloc.abf 4.4% / block-coloc 4.4%; G1: share 1.3%, above t_abf coloc.abf 37.5% / block-coloc 62.5%; G2: share 39.3%, above t_abf coloc.abf 4.9% / block-coloc 4.9%

G1 by trigger: trait 2 has >= 2 leads: 41.4% of G1, gain +9.7 pp within; only trait 1 has >= 2 leads: 58.6% of G1, gain +15.4 pp within.

G1 flips: 98 loci missed by coloc.abf and called by block-coloc (selected block contains the shared causal variant in 100.0%); 7 called by coloc.abf and missed by block-coloc. Median block size 46 SNPs.

| scenario | G0 / G1 / G2 share | contribution G1 (pp, within scenario) | contribution G0+G2 (pp) | total |
|---|---|---|---|---|
| single_shared | 84 / 14 / 1 | +0.0 | -0.2 | -0.2 |
| single_double_shared | 78 / 15 / 8 | +0.3 | -0.5 | -0.2 |
| double_shared2 | 33 / 55 / 12 | +5.2 | -0.2 | +5.0 |
| double_shared1 | 32 / 33 / 36 | +9.7 | -0.2 | +9.5 |

Without fallback: G2 positives get no result; 6.2% of positives are called by block-coloc in G2 and lost. Its lower cut-off (0.107) comes from the 0.0 / NA negatives in G2 (no-valid-block negatives are removed from the top of the ranking).

## Stage 2 (stage-1 cut-offs)

### h2_gwas 0.001, h2_eqtl 0.015

Cut-offs (5% FPR on calibration windows): coloc.abf 0.1412, block-coloc 0.1448, block-coloc without fallback 0.1067. Positive test loci: 1200. Detection: coloc.abf 52.2%, block-coloc 53.0% (difference +0.8), without fallback 52.9%.

| branch | share of positives | coloc.abf | block-coloc | contribution to the difference (pp) | of which: PP4 change at common cut-off | of which: cut-off change |
|---|---|---|---|---|---|---|
| G0 gate not passed (PP4 identical) | 74.4% | 49.0 | 48.7 | -0.25 | +0.00 | -0.25 |
| G1 valid block (PP4 = max block) | 14.6% | 81.1 | 88.6 | +1.08 | +1.08 | +0.00 |
| G2 no valid block (PP4 identical, fallback) | 11.0% | 34.8 | 34.8 | +0.00 | +0.00 | +0.00 |
| total | 100% | 52.2 | 53.0 | +0.83 | | |

Negatives (calibration windows) above coloc.abf's cut-off, by branch — this sets block-coloc's own cut-off:

- G0: share 80.5%, above t_abf coloc.abf 1.9% / block-coloc 1.9%; G1: share 0.7%, above t_abf coloc.abf 25.0% / block-coloc 75.0%; G2: share 18.8%, above t_abf coloc.abf 4.4% / block-coloc 4.4%

G1 by trigger: trait 2 has >= 2 leads: 62.9% of G1, gain +6.4 pp within; only trait 1 has >= 2 leads: 37.1% of G1, gain +9.2 pp within.

G1 flips: 13 loci missed by coloc.abf and called by block-coloc (selected block contains the shared causal variant in 92.3%); 0 called by coloc.abf and missed by block-coloc. Median block size 46 SNPs.

| scenario | G0 / G1 / G2 share | contribution G1 (pp, within scenario) | contribution G0+G2 (pp) | total |
|---|---|---|---|---|
| single_shared | 86 / 13 / 1 | +0.0 | +0.0 | +0.0 |
| single_double_shared | 77 / 14 / 9 | +0.7 | -0.7 | +0.0 |
| double_shared2 | 67 / 19 / 14 | +1.3 | -0.3 | +1.0 |
| double_shared1 | 68 / 12 / 20 | +2.3 | +0.0 | +2.3 |

Without fallback: G2 positives get no result; 3.8% of positives are called by block-coloc in G2 and lost. Its lower cut-off (0.107) comes from the 0.0 / NA negatives in G2 (no-valid-block negatives are removed from the top of the ranking).

### h2_gwas 0.003, h2_eqtl 0.0075

Cut-offs (5% FPR on calibration windows): coloc.abf 0.1412, block-coloc 0.1448, block-coloc without fallback 0.1067. Positive test loci: 1200. Detection: coloc.abf 51.4%, block-coloc 52.9% (difference +1.5), without fallback 54.7%.

| branch | share of positives | coloc.abf | block-coloc | contribution to the difference (pp) | of which: PP4 change at common cut-off | of which: cut-off change |
|---|---|---|---|---|---|---|
| G0 gate not passed (PP4 identical) | 73.6% | 50.5 | 49.8 | -0.50 | +0.00 | -0.50 |
| G1 valid block (PP4 = max block) | 12.9% | 76.1 | 91.6 | +2.00 | +2.08 | -0.08 |
| G2 no valid block (PP4 identical, fallback) | 13.5% | 32.7 | 32.7 | +0.00 | +0.00 | +0.00 |
| total | 100% | 51.4 | 52.9 | +1.50 | | |

Negatives (calibration windows) above coloc.abf's cut-off, by branch — this sets block-coloc's own cut-off:

- G0: share 80.3%, above t_abf coloc.abf 5.4% / block-coloc 5.4%; G1: share 0.5%, above t_abf coloc.abf 66.7% / block-coloc 33.3%; G2: share 19.2%, above t_abf coloc.abf 7.8% / block-coloc 7.8%

G1 by trigger: trait 2 has >= 2 leads: 23.2% of G1, gain +5.6 pp within; only trait 1 has >= 2 leads: 76.8% of G1, gain +18.5 pp within.

G1 flips: 27 loci missed by coloc.abf and called by block-coloc (selected block contains the shared causal variant in 92.6%); 3 called by coloc.abf and missed by block-coloc. Median block size 56 SNPs.

| scenario | G0 / G1 / G2 share | contribution G1 (pp, within scenario) | contribution G0+G2 (pp) | total |
|---|---|---|---|---|
| single_shared | 89 / 7 / 3 | +0.3 | -0.7 | -0.3 |
| single_double_shared | 92 / 3 / 5 | +0.0 | -1.0 | -1.0 |
| double_shared2 | 57 / 28 / 15 | +5.3 | -0.3 | +5.0 |
| double_shared1 | 56 / 13 / 30 | +2.3 | +0.0 | +2.3 |

Without fallback: G2 positives get no result; 4.4% of positives are called by block-coloc in G2 and lost. Its lower cut-off (0.107) comes from the 0.0 / NA negatives in G2 (no-valid-block negatives are removed from the top of the ranking).

### h2_gwas 0.003, h2_eqtl 0.03

Cut-offs (5% FPR on calibration windows): coloc.abf 0.1412, block-coloc 0.1448, block-coloc without fallback 0.1067. Positive test loci: 1200. Detection: coloc.abf 83.1%, block-coloc 90.2% (difference +7.1), without fallback 86.2%.

| branch | share of positives | coloc.abf | block-coloc | contribution to the difference (pp) | of which: PP4 change at common cut-off | of which: cut-off change |
|---|---|---|---|---|---|---|
| G0 gate not passed (PP4 identical) | 32.9% | 89.9 | 89.6 | -0.08 | +0.00 | -0.08 |
| G1 valid block (PP4 = max block) | 56.3% | 85.7 | 98.5 | +7.25 | +7.25 | +0.00 |
| G2 no valid block (PP4 identical, fallback) | 10.8% | 48.8 | 48.1 | -0.08 | +0.00 | -0.08 |
| total | 100% | 83.1 | 90.2 | +7.08 | | |

Negatives (calibration windows) above coloc.abf's cut-off, by branch — this sets block-coloc's own cut-off:

- G0: share 45.3%, above t_abf coloc.abf 0.4% / block-coloc 0.4%; G1: share 1.2%, above t_abf coloc.abf 14.3% / block-coloc 42.9%; G2: share 53.5%, above t_abf coloc.abf 1.6% / block-coloc 1.6%

G1 by trigger: trait 2 has >= 2 leads: 75.4% of G1, gain +13.1 pp within; only trait 1 has >= 2 leads: 24.6% of G1, gain +12.0 pp within.

G1 flips: 89 loci missed by coloc.abf and called by block-coloc (selected block contains the shared causal variant in 100.0%); 2 called by coloc.abf and missed by block-coloc. Median block size 45 SNPs.

| scenario | G0 / G1 / G2 share | contribution G1 (pp, within scenario) | contribution G0+G2 (pp) | total |
|---|---|---|---|---|
| single_shared | 83 / 17 / 0 | +0.0 | +0.0 | +0.0 |
| single_double_shared | 38 / 54 / 8 | +4.7 | -0.3 | +4.3 |
| double_shared2 | 5 / 90 / 5 | +6.3 | +0.0 | +6.3 |
| double_shared1 | 6 / 65 / 29 | +18.0 | -0.3 | +17.7 |

Without fallback: G2 positives get no result; 5.2% of positives are called by block-coloc in G2 and lost. Its lower cut-off (0.107) comes from the 0.0 / NA negatives in G2 (no-valid-block negatives are removed from the top of the ranking).

### h2_gwas 0.003, h2_eqtl 0.05

Cut-offs (5% FPR on calibration windows): coloc.abf 0.1412, block-coloc 0.1448, block-coloc without fallback 0.1067. Positive test loci: 1200. Detection: coloc.abf 82.2%, block-coloc 95.9% (difference +13.7), without fallback 94.3%.

| branch | share of positives | coloc.abf | block-coloc | contribution to the difference (pp) | of which: PP4 change at common cut-off | of which: cut-off change |
|---|---|---|---|---|---|---|
| G0 gate not passed (PP4 identical) | 23.2% | 96.0 | 96.0 | +0.00 | +0.00 | +0.00 |
| G1 valid block (PP4 = max block) | 73.0% | 80.0 | 98.7 | +13.67 | +13.67 | +0.00 |
| G2 no valid block (PP4 identical, fallback) | 3.8% | 41.3 | 41.3 | +0.00 | +0.00 | +0.00 |
| total | 100% | 82.2 | 95.9 | +13.67 | | |

Negatives (calibration windows) above coloc.abf's cut-off, by branch — this sets block-coloc's own cut-off:

- G0: share 42.7%, above t_abf coloc.abf 0.0% / block-coloc 0.0%; G1: share 1.0%, above t_abf coloc.abf 0.0% / block-coloc 50.0%; G2: share 56.3%, above t_abf coloc.abf 0.0% / block-coloc 0.0%

G1 by trigger: trait 2 has >= 2 leads: 92.9% of G1, gain +19.7 pp within; only trait 1 has >= 2 leads: 7.1% of G1, gain +6.5 pp within.

G1 flips: 165 loci missed by coloc.abf and called by block-coloc (selected block contains the shared causal variant in 100.0%); 1 called by coloc.abf and missed by block-coloc. Median block size 33 SNPs.

| scenario | G0 / G1 / G2 share | contribution G1 (pp, within scenario) | contribution G0+G2 (pp) | total |
|---|---|---|---|---|
| single_shared | 82 / 18 / 0 | +0.0 | +0.0 | +0.0 |
| single_double_shared | 9 / 87 / 4 | +13.3 | +0.0 | +13.3 |
| double_shared2 | 2 / 97 / 2 | +9.7 | +0.0 | +9.7 |
| double_shared1 | 1 / 90 / 9 | +31.7 | +0.0 | +31.7 |

Without fallback: G2 positives get no result; 1.6% of positives are called by block-coloc in G2 and lost. Its lower cut-off (0.107) comes from the 0.0 / NA negatives in G2 (no-valid-block negatives are removed from the top of the ranking).

### h2_gwas 0.01, h2_eqtl 0.015

Cut-offs (5% FPR on calibration windows): coloc.abf 0.1412, block-coloc 0.1448, block-coloc without fallback 0.1067. Positive test loci: 1200. Detection: coloc.abf 70.8%, block-coloc 76.9% (difference +6.2), without fallback 72.7%.

| branch | share of positives | coloc.abf | block-coloc | contribution to the difference (pp) | of which: PP4 change at common cut-off | of which: cut-off change |
|---|---|---|---|---|---|---|
| G0 gate not passed (PP4 identical) | 52.8% | 74.1 | 73.8 | -0.17 | +0.00 | -0.17 |
| G1 valid block (PP4 = max block) | 31.6% | 77.8 | 97.9 | +6.33 | +6.33 | +0.00 |
| G2 no valid block (PP4 identical, fallback) | 15.6% | 44.9 | 44.9 | +0.00 | +0.00 | +0.00 |
| total | 100% | 70.8 | 76.9 | +6.17 | | |

Negatives (calibration windows) above coloc.abf's cut-off, by branch — this sets block-coloc's own cut-off:

- G0: share 56.0%, above t_abf coloc.abf 4.8% / block-coloc 4.8%; G1: share 1.2%, above t_abf coloc.abf 14.3% / block-coloc 28.6%; G2: share 42.8%, above t_abf coloc.abf 3.9% / block-coloc 3.9%

G1 by trigger: trait 2 has >= 2 leads: 38.0% of G1, gain +14.6 pp within; only trait 1 has >= 2 leads: 62.0% of G1, gain +23.4 pp within.

G1 flips: 76 loci missed by coloc.abf and called by block-coloc (selected block contains the shared causal variant in 100.0%); 0 called by coloc.abf and missed by block-coloc. Median block size 49 SNPs.

| scenario | G0 / G1 / G2 share | contribution G1 (pp, within scenario) | contribution G0+G2 (pp) | total |
|---|---|---|---|---|
| single_shared | 82 / 16 / 2 | +0.0 | +0.0 | +0.0 |
| single_double_shared | 74 / 16 / 10 | +0.3 | -0.3 | +0.0 |
| double_shared2 | 27 / 61 / 11 | +12.0 | -0.3 | +11.7 |
| double_shared1 | 28 / 33 / 39 | +13.0 | +0.0 | +13.0 |

Without fallback: G2 positives get no result; 7.0% of positives are called by block-coloc in G2 and lost. Its lower cut-off (0.107) comes from the 0.0 / NA negatives in G2 (no-valid-block negatives are removed from the top of the ranking).

### h2_gwas 0.02, h2_eqtl 0.015

Cut-offs (5% FPR on calibration windows): coloc.abf 0.1412, block-coloc 0.1448, block-coloc without fallback 0.1067. Positive test loci: 1200. Detection: coloc.abf 69.8%, block-coloc 76.4% (difference +6.7), without fallback 72.2%.

| branch | share of positives | coloc.abf | block-coloc | contribution to the difference (pp) | of which: PP4 change at common cut-off | of which: cut-off change |
|---|---|---|---|---|---|---|
| G0 gate not passed (PP4 identical) | 49.7% | 72.1 | 72.0 | -0.08 | +0.00 | -0.08 |
| G1 valid block (PP4 = max block) | 34.4% | 78.2 | 97.8 | +6.75 | +6.75 | +0.00 |
| G2 no valid block (PP4 identical, fallback) | 15.9% | 44.0 | 44.0 | +0.00 | +0.00 | +0.00 |
| total | 100% | 69.8 | 76.4 | +6.67 | | |

Negatives (calibration windows) above coloc.abf's cut-off, by branch — this sets block-coloc's own cut-off:

- G0: share 50.7%, above t_abf coloc.abf 4.6% / block-coloc 4.6%; G1: share 1.8%, above t_abf coloc.abf 18.2% / block-coloc 54.5%; G2: share 47.5%, above t_abf coloc.abf 4.6% / block-coloc 4.6%

G1 by trigger: trait 2 has >= 2 leads: 35.1% of G1, gain +13.8 pp within; only trait 1 has >= 2 leads: 64.9% of G1, gain +22.8 pp within.

G1 flips: 83 loci missed by coloc.abf and called by block-coloc (selected block contains the shared causal variant in 97.6%); 2 called by coloc.abf and missed by block-coloc. Median block size 49 SNPs.

| scenario | G0 / G1 / G2 share | contribution G1 (pp, within scenario) | contribution G0+G2 (pp) | total |
|---|---|---|---|---|
| single_shared | 74 / 24 / 2 | +0.0 | +0.0 | +0.0 |
| single_double_shared | 69 / 18 / 12 | +0.0 | +0.0 | +0.0 |
| double_shared2 | 27 / 62 / 11 | +13.0 | +0.0 | +13.0 |
| double_shared1 | 28 / 34 / 38 | +14.0 | -0.3 | +13.7 |

Without fallback: G2 positives get no result; 7.0% of positives are called by block-coloc in G2 and lost. Its lower cut-off (0.107) comes from the 0.0 / NA negatives in G2 (no-valid-block negatives are removed from the top of the ranking).

## Normal effects, default (own cut-offs)

### Normal effects, default

Cut-offs (5% FPR on calibration windows): coloc.abf 0.0893, block-coloc 0.0900, block-coloc without fallback 0.0867. Positive test loci: 2400. Detection: coloc.abf 37.8%, block-coloc 39.7% (difference +1.9), without fallback 38.2%.

| branch | share of positives | coloc.abf | block-coloc | contribution to the difference (pp) | of which: PP4 change at common cut-off | of which: cut-off change |
|---|---|---|---|---|---|---|
| G0 gate not passed (PP4 identical) | 81.2% | 32.6 | 32.6 | -0.04 | +0.00 | -0.04 |
| G1 valid block (PP4 = max block) | 11.4% | 81.0 | 98.5 | +2.00 | +2.00 | +0.00 |
| G2 no valid block (PP4 identical, fallback) | 7.3% | 28.4 | 27.3 | -0.08 | +0.00 | -0.08 |
| total | 100% | 37.8 | 39.7 | +1.88 | | |

Negatives (calibration windows) above coloc.abf's cut-off, by branch — this sets block-coloc's own cut-off:

- G0: share 84.7%, above t_abf coloc.abf 5.2% / block-coloc 5.2%; G1: share 0.5%, above t_abf coloc.abf 33.3% / block-coloc 66.7%; G2: share 14.8%, above t_abf coloc.abf 2.8% / block-coloc 2.8%

G1 by trigger: trait 2 has >= 2 leads: 48.9% of G1, gain +14.2 pp within; only trait 1 has >= 2 leads: 51.1% of G1, gain +20.7 pp within.

G1 flips: 48 loci missed by coloc.abf and called by block-coloc (selected block contains the shared causal variant in 100.0%); 0 called by coloc.abf and missed by block-coloc. Median block size 44 SNPs.

| scenario | G0 / G1 / G2 share | contribution G1 (pp, within scenario) | contribution G0+G2 (pp) | total |
|---|---|---|---|---|
| single_shared | 92 / 7 / 1 | +0.0 | +0.0 | +0.0 |
| single_double_shared | 88 / 7 / 4 | +1.0 | -0.2 | +0.8 |
| double_shared2 | 74 / 18 / 8 | +3.0 | -0.2 | +2.8 |
| double_shared1 | 71 / 14 / 15 | +4.0 | -0.2 | +3.8 |

Without fallback: G2 positives get no result; 2.0% of positives are called by block-coloc in G2 and lost. Its lower cut-off (0.087) comes from the 0.0 / NA negatives in G2 (no-valid-block negatives are removed from the top of the ranking).

## Normal effects, h2_eqtl 5% (own cut-offs)

### Normal effects, h2_eqtl 5%

Cut-offs (5% FPR on calibration windows): coloc.abf 0.0697, block-coloc 0.0702, block-coloc without fallback 0.0676. Positive test loci: 1200. Detection: coloc.abf 46.2%, block-coloc 51.7% (difference +5.5), without fallback 50.0%.

| branch | share of positives | coloc.abf | block-coloc | contribution to the difference (pp) | of which: PP4 change at common cut-off | of which: cut-off change |
|---|---|---|---|---|---|---|
| G0 gate not passed (PP4 identical) | 67.0% | 41.3 | 41.3 | +0.00 | +0.00 | +0.00 |
| G1 valid block (PP4 = max block) | 22.3% | 73.1 | 97.8 | +5.50 | +5.50 | +0.00 |
| G2 no valid block (PP4 identical, fallback) | 10.7% | 20.3 | 20.3 | +0.00 | +0.00 | +0.00 |
| total | 100% | 46.2 | 51.7 | +5.50 | | |

Negatives (calibration windows) above coloc.abf's cut-off, by branch — this sets block-coloc's own cut-off:

- G0: share 73.7%, above t_abf coloc.abf 5.7% / block-coloc 5.7%; G1: share 0.8%, above t_abf coloc.abf 20.0% / block-coloc 60.0%; G2: share 25.5%, above t_abf coloc.abf 2.6% / block-coloc 2.6%

G1 by trigger: trait 2 has >= 2 leads: 66.8% of G1, gain +27.9 pp within; only trait 1 has >= 2 leads: 33.2% of G1, gain +18.0 pp within.

G1 flips: 67 loci missed by coloc.abf and called by block-coloc (selected block contains the shared causal variant in 100.0%); 1 called by coloc.abf and missed by block-coloc. Median block size 47 SNPs.

| scenario | G0 / G1 / G2 share | contribution G1 (pp, within scenario) | contribution G0+G2 (pp) | total |
|---|---|---|---|---|
| single_shared | 89 / 10 / 1 | +0.0 | +0.0 | +0.0 |
| single_double_shared | 73 / 16 / 10 | +4.7 | +0.0 | +4.7 |
| double_shared2 | 53 / 37 / 10 | +10.0 | +0.0 | +10.0 |
| double_shared1 | 53 / 26 / 22 | +7.3 | +0.0 | +7.3 |

Without fallback: G2 positives get no result; 2.2% of positives are called by block-coloc in G2 and lost. Its lower cut-off (0.068) comes from the 0.0 / NA negatives in G2 (no-valid-block negatives are removed from the top of the ranking).

