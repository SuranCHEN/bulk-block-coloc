# Bulk block-coloc: main comparison (stage 1)

Generated 2026-10-26 18:10. Setting: GWAS N 20,000, q 0.3% per causal SNP; eQTL N 1,000, q 1.5% per causal SNP; causal |r| < 0.2.
Calibration windows: 31263287, 37199593, 39694225, 42313324, 44324730, 46518082. Test windows: 21059630, 24009444, 26904965, 29033189, 33343871, 34894283.
Loci: 3000 calibration, 3000 test (per method). Target FPR 5% on calibration windows (negatives: single_distinct and double_distinct).
Extra negatives gwas_only and null_both are reported for the false-positive rate only (not calibrated, not pooled).

## Three-line conclusion

1. **检出率**：误报率在校准窗口上匹配到 5%（测试窗口实际为 4.1%–4.8%，coloc.susie 为 2.8%）时，block-coloc（条件 lead 选择，本文方法）检出 93.3% 的共定位位点，高于整区域 coloc.abf 的 91.2%（同一批位点配对差 +2.2 个百分点，95% CI +1.1 至 +3.2）和 coloc.susie + coloc.abf fallback 的 90.1%（+3.2，+2.1 至 +4.4）。提升主要来自 double_shared2（94.2% 对 89.5%），single_shared 上没有损失（98.2% 对 97.8%）。但它没有超过 coloc.signals (cond)（95.1%），与 block-coloc（交接包原版，94.2%）也相近。
2. **是否切块**：条件 lead 选择基本消除了单信号位点被误切块：在 single_shared、single_distinct 上走 block 的比例为 1.5% 和 0.3%；边际 lead 选择为 3.7% 和 4.3%；交接包原版在 single_shared 上为 84.8%。但多信号常被漏判：single_double_shared 只有 26.0% 走 block，double_shared2 为 78.8%，因为 eQTL 的第二个信号（N = 1,000，因果 |z| ≈ 3.9）经常过不了 lead 阈值。single_double_shared 中当 eQTL 特有的信号比共享信号更显著时，所有方法检出率都明显下降（整区域 coloc.abf 73.0%，block-coloc（条件 lead 选择）76.2%，coloc.susie + fallback 54.7%）。
3. **只有一个性状有信号时的误报**（额外阴性，未用于定判定线）：gwas_only 上，基于 coloc.abf 的方法误报 8.8%–11.8%，超过目标 5%；block-coloc（条件 lead 选择）为 11.8%，略高于整区域 coloc.abf 的 8.8%，原因是它的判定线更低（0.094 对 0.109），这些位点上两者给出的 PP4 在 600 个中有 599 个完全相同。coloc.signals 为 3.5%，coloc.susie 为 0.2%；null_both 上所有方法均为 0%。

## 1. Main table (test windows, frozen cut-offs)

n_positive = loci called colocalised (PP4 above the method's cut-off). For shared scenarios rate is the detection rate;
for not-shared scenarios it is the false-positive rate. Denominators include every locus; no result counts as not called.
'not shared (extra)' rows (gwas_only, null_both) use the same frozen cut-offs but were not used to set them.

| scenario | truth | method | lead_p | n_total | n_positive | rate | 95% CI | no_call_rate | branch_block |
|---|---|---|---|---|---|---|---|---|---|
| single_shared | shared | block-coloc (conditional lead selection, proposed) | 1e-04 | 600 | 589 | 98.2% | [96.7%, 99.0%] | 0.2% | 1.5% |
| single_shared | shared | block-coloc (marginal lead selection) | 1e-05 | 600 | 587 | 97.8% | [96.3%, 98.7%] | 0.0% | 3.7% |
| single_shared | shared | block-coloc (original handoff) | 1e-05 | 600 | 585 | 97.5% | [95.9%, 98.5%] | 0.3% | 84.8% |
| single_shared | shared | coloc.abf, whole region | - | 600 | 587 | 97.8% | [96.3%, 98.7%] | 0.0% | 0.0% |
| single_shared | shared | coloc.susie | - | 600 | 215 | 35.8% | [32.1%, 39.7%] | 64.0% | 0.0% |
| single_shared | shared | coloc.susie + coloc.abf fallback | - | 600 | 586 | 97.7% | [96.1%, 98.6%] | 0.0% | 0.0% |
| single_shared | shared | coloc.signals (mask) | - | 600 | 573 | 95.5% | [93.5%, 96.9%] | 0.0% | 0.0% |
| single_shared | shared | coloc.signals (cond) | - | 600 | 573 | 95.5% | [93.5%, 96.9%] | 0.0% | 0.0% |
| single_double_shared | shared | block-coloc (conditional lead selection, proposed) | 1e-04 | 600 | 526 | 87.7% | [84.8%, 90.1%] | 1.3% | 26.0% |
| single_double_shared | shared | block-coloc (marginal lead selection) | 1e-05 | 600 | 508 | 84.7% | [81.6%, 87.3%] | 3.3% | 21.7% |
| single_double_shared | shared | block-coloc (original handoff) | 1e-05 | 600 | 526 | 87.7% | [84.8%, 90.1%] | 0.2% | 56.8% |
| single_double_shared | shared | coloc.abf, whole region | - | 600 | 517 | 86.2% | [83.2%, 88.7%] | 0.0% | 0.0% |
| single_double_shared | shared | coloc.susie | - | 600 | 228 | 38.0% | [34.2%, 41.9%] | 45.2% | 0.0% |
| single_double_shared | shared | coloc.susie + coloc.abf fallback | - | 600 | 455 | 75.8% | [72.3%, 79.1%] | 0.0% | 0.0% |
| single_double_shared | shared | coloc.signals (mask) | - | 600 | 529 | 88.2% | [85.3%, 90.5%] | 0.0% | 0.0% |
| single_double_shared | shared | coloc.signals (cond) | - | 600 | 544 | 90.7% | [88.1%, 92.7%] | 0.0% | 0.0% |
| double_shared2 | shared | block-coloc (conditional lead selection, proposed) | 1e-04 | 600 | 565 | 94.2% | [92.0%, 95.8%] | 3.3% | 78.8% |
| double_shared2 | shared | block-coloc (marginal lead selection) | 1e-05 | 600 | 569 | 94.8% | [92.8%, 96.3%] | 1.8% | 56.2% |
| double_shared2 | shared | block-coloc (original handoff) | 1e-05 | 600 | 585 | 97.5% | [95.9%, 98.5%] | 0.3% | 94.2% |
| double_shared2 | shared | coloc.abf, whole region | - | 600 | 537 | 89.5% | [86.8%, 91.7%] | 0.0% | 0.0% |
| double_shared2 | shared | coloc.susie | - | 600 | 301 | 50.2% | [46.2%, 54.2%] | 49.3% | 0.0% |
| double_shared2 | shared | coloc.susie + coloc.abf fallback | - | 600 | 581 | 96.8% | [95.1%, 98.0%] | 0.0% | 0.0% |
| double_shared2 | shared | coloc.signals (mask) | - | 600 | 579 | 96.5% | [94.7%, 97.7%] | 0.0% | 0.0% |
| double_shared2 | shared | coloc.signals (cond) | - | 600 | 594 | 99.0% | [97.8%, 99.5%] | 0.0% | 0.0% |
| single_distinct | not shared | block-coloc (conditional lead selection, proposed) | 1e-04 | 600 |  41 | 6.8% | [5.1%, 9.1%] | 0.3% | 0.3% |
| single_distinct | not shared | block-coloc (marginal lead selection) | 1e-05 | 600 |  35 | 5.8% | [4.2%, 8.0%] | 4.0% | 4.3% |
| single_distinct | not shared | block-coloc (original handoff) | 1e-05 | 600 |  39 | 6.5% | [4.8%, 8.8%] | 0.8% | 3.2% |
| single_distinct | not shared | coloc.abf, whole region | - | 600 |  37 | 6.2% | [4.5%, 8.4%] | 0.0% | 0.0% |
| single_distinct | not shared | coloc.susie | - | 600 |   6 | 1.0% | [0.5%, 2.2%] | 67.2% | 0.0% |
| single_distinct | not shared | coloc.susie + coloc.abf fallback | - | 600 |  36 | 6.0% | [4.4%, 8.2%] | 0.0% | 0.0% |
| single_distinct | not shared | coloc.signals (mask) | - | 600 |  21 | 3.5% | [2.3%, 5.3%] | 0.0% | 0.0% |
| single_distinct | not shared | coloc.signals (cond) | - | 600 |  24 | 4.0% | [2.7%, 5.9%] | 0.0% | 0.0% |
| double_distinct | not shared | block-coloc (conditional lead selection, proposed) | 1e-04 | 600 |  14 | 2.3% | [1.4%, 3.9%] | 80.2% | 81.0% |
| double_distinct | not shared | block-coloc (marginal lead selection) | 1e-05 | 600 |  20 | 3.3% | [2.2%, 5.1%] | 57.8% | 59.8% |
| double_distinct | not shared | block-coloc (original handoff) | 1e-05 | 600 |  18 | 3.0% | [1.9%, 4.7%] | 24.5% | 28.0% |
| double_distinct | not shared | coloc.abf, whole region | - | 600 |  17 | 2.8% | [1.8%, 4.5%] | 0.0% | 0.0% |
| double_distinct | not shared | coloc.susie | - | 600 |  28 | 4.7% | [3.2%, 6.7%] | 47.3% | 0.0% |
| double_distinct | not shared | coloc.susie + coloc.abf fallback | - | 600 |  19 | 3.2% | [2.0%, 4.9%] | 0.0% | 0.0% |
| double_distinct | not shared | coloc.signals (mask) | - | 600 |  28 | 4.7% | [3.2%, 6.7%] | 0.2% | 0.0% |
| double_distinct | not shared | coloc.signals (cond) | - | 600 |  34 | 5.7% | [4.1%, 7.8%] | 0.0% | 0.0% |
| gwas_only | not shared (extra) | block-coloc (conditional lead selection, proposed) | 1e-04 | 600 |  71 | 11.8% | [9.5%, 14.7%] | 0.2% | 0.2% |
| gwas_only | not shared (extra) | block-coloc (marginal lead selection) | 1e-05 | 600 |  53 | 8.8% | [6.8%, 11.4%] | 0.0% | 0.0% |
| gwas_only | not shared (extra) | block-coloc (original handoff) | 1e-05 | 600 |  53 | 8.8% | [6.8%, 11.4%] | 6.0% | 13.5% |
| gwas_only | not shared (extra) | coloc.abf, whole region | - | 600 |  53 | 8.8% | [6.8%, 11.4%] | 0.0% | 0.0% |
| gwas_only | not shared (extra) | coloc.susie | - | 600 |   1 | 0.2% | [0.0%, 0.9%] | 99.8% | 0.0% |
| gwas_only | not shared (extra) | coloc.susie + coloc.abf fallback | - | 600 |  56 | 9.3% | [7.3%, 11.9%] | 0.0% | 0.0% |
| gwas_only | not shared (extra) | coloc.signals (mask) | - | 600 |  21 | 3.5% | [2.3%, 5.3%] | 0.0% | 0.0% |
| gwas_only | not shared (extra) | coloc.signals (cond) | - | 600 |  21 | 3.5% | [2.3%, 5.3%] | 0.0% | 0.0% |
| null_both | not shared (extra) | block-coloc (conditional lead selection, proposed) | 1e-04 | 600 |   0 | 0.0% | [0.0%, 0.6%] | 0.0% | 0.0% |
| null_both | not shared (extra) | block-coloc (marginal lead selection) | 1e-05 | 600 |   0 | 0.0% | [0.0%, 0.6%] | 0.0% | 0.0% |
| null_both | not shared (extra) | block-coloc (original handoff) | 1e-05 | 600 |   0 | 0.0% | [0.0%, 0.6%] | 26.3% | 32.8% |
| null_both | not shared (extra) | coloc.abf, whole region | - | 600 |   0 | 0.0% | [0.0%, 0.6%] | 0.0% | 0.0% |
| null_both | not shared (extra) | coloc.susie | - | 600 |   0 | 0.0% | [0.0%, 0.6%] | 100.0% | 0.0% |
| null_both | not shared (extra) | coloc.susie + coloc.abf fallback | - | 600 |   0 | 0.0% | [0.0%, 0.6%] | 0.0% | 0.0% |
| null_both | not shared (extra) | coloc.signals (mask) | - | 600 |   0 | 0.0% | [0.0%, 0.6%] | 0.0% | 0.0% |
| null_both | not shared (extra) | coloc.signals (cond) | - | 600 |   0 | 0.0% | [0.0%, 0.6%] | 0.0% | 0.0% |

## 1b. False-positive rate when at most one trait has a signal (test windows)

gwas_only: GWAS has one causal variant, the eQTL has none. null_both: neither trait has a signal.
Same frozen cut-offs as above; these loci were not used to set them.

| method | gwas_only | null_both |
|---|---|---|
| block-coloc (conditional lead selection, proposed) | 11.8% [9.5%, 14.7%] | 0.0% [0.0%, 0.6%] |
| block-coloc (marginal lead selection) | 8.8% [6.8%, 11.4%] | 0.0% [0.0%, 0.6%] |
| block-coloc (original handoff) | 8.8% [6.8%, 11.4%] | 0.0% [0.0%, 0.6%] |
| coloc.abf, whole region | 8.8% [6.8%, 11.4%] | 0.0% [0.0%, 0.6%] |
| coloc.susie | 0.2% [0.0%, 0.9%] | 0.0% [0.0%, 0.6%] |
| coloc.susie + coloc.abf fallback | 9.3% [7.3%, 11.9%] | 0.0% [0.0%, 0.6%] |
| coloc.signals (mask) | 3.5% [2.3%, 5.3%] | 0.0% [0.0%, 0.6%] |
| coloc.signals (cond) | 3.5% [2.3%, 5.3%] | 0.0% [0.0%, 0.6%] |

## 2. Pooled over scenarios (test windows)

TPR pools single_shared, single_double_shared, double_shared2; FPR pools single_distinct and double_distinct. AUC does not depend on any cut-off.

| method | lead_p | n_pos | TPR | n_neg | FPR | AUC | no_call | branch_block |
|---|---|---|---|---|---|---|---|---|
| block-coloc (conditional lead selection, proposed) | 1e-04 | 1800 | 93.3% [92.1%, 94.4%] | 1200 | 4.6% [3.5%, 5.9%] | 0.971 [0.965, 0.977] | 17.1% | 37.5% |
| block-coloc (marginal lead selection) | 1e-05 | 1800 | 92.4% [91.1%, 93.6%] | 1200 | 4.6% [3.5%, 5.9%] | 0.966 [0.959, 0.973] | 13.4% | 29.1% |
| block-coloc (original handoff) | 1e-05 | 1800 | 94.2% [93.0%, 95.2%] | 1200 | 4.8% [3.7%, 6.1%] | 0.977 [0.972, 0.982] | 5.2% | 53.4% |
| coloc.abf, whole region | - | 1800 | 91.2% [89.8%, 92.4%] | 1200 | 4.5% [3.5%, 5.8%] | 0.964 [0.957, 0.972] | 0.0% | 0.0% |
| coloc.susie | - | 1800 | 41.3% [39.1%, 43.6%] | 1200 | 2.8% [2.0%, 3.9%] | 0.610 [0.593, 0.627] | 54.6% | 0.0% |
| coloc.susie + coloc.abf fallback | - | 1800 | 90.1% [88.6%, 91.4%] | 1200 | 4.6% [3.5%, 5.9%] | 0.946 [0.936, 0.954] | 0.0% | 0.0% |
| coloc.signals (mask) | - | 1800 | 93.4% [92.1%, 94.4%] | 1200 | 4.1% [3.1%, 5.4%] | 0.974 [0.968, 0.980] | 0.0% | 0.0% |
| coloc.signals (cond) | - | 1800 | 95.1% [94.0%, 96.0%] | 1200 | 4.8% [3.8%, 6.2%] | 0.983 [0.978, 0.987] | 0.0% | 0.0% |

## 3. Cut-offs chosen on calibration windows

lead_p = lead threshold (GWAS / eQTL when they differ). '*' marks the lead threshold kept for each block-coloc version.

| method | lead_p | cutoff | calib_FPR | calib_TPR |
|---|---|---|---|---|
| block-coloc (conditional lead selection, proposed) | 1e-03 | 0.0778 | 5.0% | 93.3% |
| block-coloc (conditional lead selection, proposed) | 1e-04 * | 0.0936 | 5.0% | 93.7% |
| block-coloc (conditional lead selection, proposed) | 1e-05 | 0.1055 | 5.0% | 93.7% |
| block-coloc (conditional lead selection, proposed) [earlier bulk thresholds] | 5e-08 / 1e-05 | 0.1055 | 5.0% | 93.6% |
| block-coloc (conditional lead selection, proposed) + whole-region fallback | 1e-03 | 0.1171 | 5.0% | 94.8% |
| block-coloc (conditional lead selection, proposed) + whole-region fallback | 1e-04 * | 0.1112 | 5.0% | 94.2% |
| block-coloc (conditional lead selection, proposed) + whole-region fallback | 1e-05 | 0.1094 | 5.0% | 93.8% |
| block-coloc (conditional lead selection, proposed) + whole-region fallback [earlier bulk thresholds] | 5e-08 / 1e-05 | 0.1094 | 5.0% | 93.6% |
| block-coloc (marginal lead selection) | 1e-03 | 0.2037 | 5.0% | 89.4% |
| block-coloc (marginal lead selection) | 1e-04 | 0.1094 | 5.0% | 92.8% |
| block-coloc (marginal lead selection) | 1e-05 * | 0.1116 | 5.0% | 93.6% |
| block-coloc (marginal lead selection) [earlier bulk thresholds] | 5e-08 / 1e-05 | 0.1067 | 5.0% | 93.9% |
| block-coloc (marginal lead selection) + whole-region fallback | 1e-03 | 0.2451 | 5.0% | 92.7% |
| block-coloc (marginal lead selection) + whole-region fallback | 1e-04 | 0.1396 | 5.0% | 94.1% |
| block-coloc (marginal lead selection) + whole-region fallback | 1e-05 * | 0.1178 | 5.0% | 93.9% |
| block-coloc (marginal lead selection) + whole-region fallback [earlier bulk thresholds] | 5e-08 / 1e-05 | 0.1116 | 5.0% | 94.2% |
| block-coloc (original handoff) | 1e-03 | 0.2125 | 5.0% | 92.8% |
| block-coloc (original handoff) | 1e-04 | 0.1319 | 5.0% | 94.0% |
| block-coloc (original handoff) | 1e-05 * | 0.1074 | 5.0% | 94.6% |
| block-coloc (original handoff) [earlier bulk thresholds] | 5e-08 / 1e-05 | 0.1028 | 5.0% | 94.7% |
| block-coloc (original handoff) + whole-region fallback | 1e-03 | 0.2347 | 5.0% | 93.3% |
| block-coloc (original handoff) + whole-region fallback | 1e-04 | 0.1469 | 5.0% | 93.9% |
| block-coloc (original handoff) + whole-region fallback | 1e-05 * | 0.1244 | 5.0% | 94.3% |
| block-coloc (original handoff) + whole-region fallback [earlier bulk thresholds] | 5e-08 / 1e-05 | 0.1171 | 5.0% | 94.2% |
| coloc.abf, whole region | - | 0.1094 | 5.0% | 91.1% |
| coloc.susie | - | 0.0185 | 5.0% | 42.9% |
| coloc.susie + coloc.abf fallback | - | 0.1028 | 5.0% | 90.8% |
| coloc.signals (mask) | - | 0.1706 | 5.0% | 92.8% |
| coloc.signals (cond) | - | 0.1714 | 5.0% | 95.2% |

## 4. Paired differences on the same test loci (McNemar)

Difference = method minus reference, in percentage points, with 95% CI and exact McNemar p-value.

| method | vs | TPR_difference | FPR_difference |
|---|---|---|---|
| block-coloc (conditional lead selection, proposed) | coloc.susie + coloc.abf fallback | +3.2 pts [+2.1, +4.4], p=6.5e-08 | +0.0 pts [-0.7, +0.7], p=1 |
| block-coloc (marginal lead selection) | coloc.susie + coloc.abf fallback | +2.3 pts [+1.3, +3.4], p=1.4e-05 | +0.0 pts [-0.6, +0.6], p=1 |
| block-coloc (original handoff) | coloc.susie + coloc.abf fallback | +4.1 pts [+3.1, +5.1], p=3.9e-16 | +0.2 pts [-0.6, +0.9], p=0.83 |
| coloc.abf, whole region | coloc.susie + coloc.abf fallback | +1.1 pts [-0.1, +2.3], p=0.1 | -0.1 pts [-0.6, +0.4], p=1 |
| coloc.susie | coloc.susie + coloc.abf fallback | -48.8 pts [-51.2, -46.3], p=<2e-16 | -1.8 pts [-3.3, -0.2], p=0.029 |
| coloc.signals (mask) | coloc.susie + coloc.abf fallback | +3.3 pts [+2.0, +4.6], p=6e-07 | -0.5 pts [-1.5, +0.5], p=0.43 |
| coloc.signals (cond) | coloc.susie + coloc.abf fallback | +4.9 pts [+3.7, +6.2], p=2.4e-14 | +0.2 pts [-0.8, +1.3], p=0.77 |
| block-coloc (conditional lead selection, proposed) | coloc.abf, whole region | +2.2 pts [+1.1, +3.2], p=5.3e-05 | +0.1 pts [-0.6, +0.8], p=1 |
| block-coloc (marginal lead selection) | coloc.abf, whole region | +1.3 pts [+0.3, +2.2], p=0.011 | +0.1 pts [-0.4, +0.6], p=1 |
| block-coloc (original handoff) | coloc.abf, whole region | +3.1 pts [+2.2, +3.9], p=1.4e-13 | +0.2 pts [-0.4, +0.9], p=0.61 |
| coloc.susie | coloc.abf, whole region | -49.8 pts [-52.4, -47.2], p=<2e-16 | -1.7 pts [-3.1, -0.2], p=0.035 |
| coloc.susie + coloc.abf fallback | coloc.abf, whole region | -1.1 pts [-2.3, +0.1], p=0.1 | +0.1 pts [-0.4, +0.6], p=1 |
| coloc.signals (mask) | coloc.abf, whole region | +2.2 pts [+1.1, +3.3], p=0.00011 | -0.4 pts [-1.5, +0.6], p=0.53 |
| coloc.signals (cond) | coloc.abf, whole region | +3.9 pts [+2.7, +5.1], p=3.8e-10 | +0.3 pts [-0.8, +1.4], p=0.66 |

## 5. Branch choice: split into blocks or not (test windows)

should_go: block when either trait has 2 causal SNPs, whole otherwise. leads_match_truth: the method's lead counts
(GWAS, eQTL) equal the true numbers of causal SNPs.

| method | scenario | true_signals | should_go | n | branch_block | correct_branch | leads_match_truth | median_leads |
|---|---|---|---|---|---|---|---|---|
| block-coloc (conditional lead selection, proposed) | single_shared | 1+1 | whole | 600 | 1.5% | 98.5% | 59.3% | 1+1 |
| block-coloc (conditional lead selection, proposed) | single_double_shared | 1+2 | block | 600 | 26.0% | 26.0% | 24.5% | 1+1 |
| block-coloc (conditional lead selection, proposed) | double_shared2 | 2+2 | block | 600 | 78.8% | 78.8% | 22.5% | 2+1 |
| block-coloc (conditional lead selection, proposed) | single_distinct | 1+1 | whole | 600 | 0.3% | 99.7% | 58.2% | 1+1 |
| block-coloc (conditional lead selection, proposed) | double_distinct | 2+2 | block | 600 | 81.0% | 81.0% | 23.7% | 2+1 |
| block-coloc (conditional lead selection, proposed) | gwas_only | 1+0 | whole | 600 | 0.2% | 99.8% | 97.7% | 1+0 |
| block-coloc (conditional lead selection, proposed) | null_both | 0+0 | whole | 600 | 0.0% | 100.0% | 96.3% | 0+0 |
| block-coloc (marginal lead selection) | single_shared | 1+1 | whole | 600 | 3.7% | 96.3% | 35.8% | 1+0 |
| block-coloc (marginal lead selection) | single_double_shared | 1+2 | block | 600 | 21.7% | 21.7% | 14.5% | 1+1 |
| block-coloc (marginal lead selection) | double_shared2 | 2+2 | block | 600 | 56.2% | 56.2% | 9.2% | 2+1 |
| block-coloc (marginal lead selection) | single_distinct | 1+1 | whole | 600 | 4.3% | 95.7% | 31.3% | 1+0 |
| block-coloc (marginal lead selection) | double_distinct | 2+2 | block | 600 | 59.8% | 59.8% | 9.3% | 2+1 |
| block-coloc (marginal lead selection) | gwas_only | 1+0 | whole | 600 | 0.0% | 100.0% | 89.5% | 1+0 |
| block-coloc (marginal lead selection) | null_both | 0+0 | whole | 600 | 0.0% | 100.0% | 99.5% | 0+0 |
| block-coloc (original handoff) | single_shared | 1+1 | whole | 600 | 84.8% | 15.2% | 78.3% | 1+1 |
| block-coloc (original handoff) | single_double_shared | 1+2 | block | 600 | 56.8% | 56.8% | 25.7% | 1+1 |
| block-coloc (original handoff) | double_shared2 | 2+2 | block | 600 | 94.2% | 94.2% | 17.8% | 2+1 |
| block-coloc (original handoff) | single_distinct | 1+1 | whole | 600 | 3.2% | 96.8% | 76.0% | 1+1 |
| block-coloc (original handoff) | double_distinct | 2+2 | block | 600 | 28.0% | 28.0% | 17.7% | 2+1 |
| block-coloc (original handoff) | gwas_only | 1+0 | whole | 600 | 13.5% | 86.5% | 0.0% | 1+2 |
| block-coloc (original handoff) | null_both | 0+0 | whole | 600 | 32.8% | 67.2% | 0.0% | 2+2 |

## 6. single_double_shared split by which eQTL signal is stronger (test windows)

The eQTL has two causal variants: one shared with the GWAS and one eQTL-only. 'stronger' compares their observed eQTL |z|.

| method | n_shared_stronger | TPR_shared_stronger | n_other_stronger | TPR_other_stronger |
|---|---|---|---|---|
| block-coloc (conditional lead selection, proposed) | 293 | 99.7% [98.1%, 99.9%] | 307 | 76.2% [71.2%, 80.6%] |
| block-coloc (marginal lead selection) | 293 | 99.3% [97.5%, 99.8%] | 307 | 70.7% [65.4%, 75.5%] |
| block-coloc (original handoff) | 293 | 100.0% [98.7%, 100.0%] | 307 | 75.9% [70.8%, 80.3%] |
| coloc.abf, whole region | 293 | 100.0% [98.7%, 100.0%] | 307 | 73.0% [67.7%, 77.6%] |
| coloc.susie | 293 | 50.9% [45.2%, 56.5%] | 307 | 25.7% [21.2%, 30.9%] |
| coloc.susie + coloc.abf fallback | 293 | 98.0% [95.6%, 99.1%] | 307 | 54.7% [49.1%, 60.2%] |
| coloc.signals (mask) | 293 | 99.7% [98.1%, 99.9%] | 307 | 77.2% [72.2%, 81.5%] |
| coloc.signals (cond) | 293 | 99.7% [98.1%, 99.9%] | 307 | 82.1% [77.4%, 86.0%] |

## 7. Localisation among detected test positives

Sub-region = best block (block branch) or the two credible sets of the best SuSiE signal pair.
Whole-locus calls and coloc.signals calls have no sub-region and are not counted here.

| method | detected | with_subregion | contains_shared | median_size_snps |
|---|---|---|---|---|
| block-coloc (conditional lead selection, proposed) | 1680 |  605 | 100.0% [99.4%, 100.0%] | 36 |
| block-coloc (marginal lead selection) | 1664 |  454 | 99.3% [98.1%, 99.8%] | 48 |
| block-coloc (original handoff) | 1696 | 1398 | 99.8% [99.4%, 99.9%] | 36 |
| coloc.abf, whole region | 1641 |    0 | NA | NA |
| coloc.susie |  744 |  744 | 99.6% [98.8%, 99.9%] | 29 |
| coloc.susie + coloc.abf fallback | 1622 |  716 | 99.6% [98.8%, 99.9%] | 28 |
| coloc.signals (mask) | 1681 |    0 | NA | NA |
| coloc.signals (cond) | 1711 |    0 | NA | NA |

## 8. Sensitivity: fallback policy and the earlier bulk lead thresholds (test windows)

'+ whole-region fallback' = when the block branch builds no valid block, use the whole-region PP4 instead of no result.
'[earlier bulk thresholds]' = lead thresholds GWAS 5e-8, eQTL 1e-5 (August bulk setting), each with its own calibrated cut-off.

| method | lead_p | cutoff | n_pos | TPR | n_neg | FPR | AUC | no_call | branch_block |
|---|---|---|---|---|---|---|---|---|---|
| block-coloc (conditional lead selection, proposed) | 1e-04 | 0.0936 | 1800 | 93.3% [92.1%, 94.4%] | 1200 | 4.6% [3.5%, 5.9%] | 0.971 | 17.1% | 37.5% |
| block-coloc (marginal lead selection) | 1e-05 | 0.1116 | 1800 | 92.4% [91.1%, 93.6%] | 1200 | 4.6% [3.5%, 5.9%] | 0.966 | 13.4% | 29.1% |
| block-coloc (original handoff) | 1e-05 | 0.1074 | 1800 | 94.2% [93.0%, 95.2%] | 1200 | 4.8% [3.7%, 6.1%] | 0.977 | 5.2% | 53.4% |
| block-coloc (conditional lead selection, proposed) + whole-region fallback | 1e-04 | 0.1112 | 1800 | 93.8% [92.6%, 94.9%] | 1200 | 4.5% [3.5%, 5.8%] | 0.977 | 0.0% | 37.5% |
| block-coloc (marginal lead selection) + whole-region fallback | 1e-05 | 0.1178 | 1800 | 93.5% [92.3%, 94.5%] | 1200 | 4.5% [3.5%, 5.8%] | 0.976 | 0.0% | 29.1% |
| block-coloc (original handoff) + whole-region fallback | 1e-05 | 0.1244 | 1800 | 93.8% [92.6%, 94.8%] | 1200 | 4.8% [3.8%, 6.2%] | 0.977 | 0.0% | 53.4% |
| block-coloc (original handoff) [earlier bulk thresholds] | 5e-08 / 1e-05 | 0.1028 | 1800 | 94.1% [92.9%, 95.1%] | 1200 | 4.3% [3.3%, 5.6%] | 0.978 | 5.0% | 52.5% |
| block-coloc (marginal lead selection) [earlier bulk thresholds] | 5e-08 / 1e-05 | 0.1067 | 1800 | 93.2% [91.9%, 94.2%] | 1200 | 4.2% [3.2%, 5.5%] | 0.973 | 12.2% | 26.3% |
| block-coloc (conditional lead selection, proposed) [earlier bulk thresholds] | 5e-08 / 1e-05 | 0.1055 | 1800 | 92.8% [91.5%, 93.9%] | 1200 | 4.3% [3.3%, 5.6%] | 0.973 | 11.8% | 23.7% |

## 9. Inputs (test windows)

Expected causal |z| about sqrt(N q): GWAS 7.7, eQTL 3.9.

| scenario | n | causal_z_gwas | causal_z_eqtl | pair_abs_r_median | gwas_no_cs | eqtl_no_cs |
|---|---|---|---|---|---|---|
| single_shared | 600 | 7.8 | 4.0 | NA | 0.0% | 64.0% |
| single_double_shared | 600 | 7.8 | 3.9 | 0.03 | 0.3% | 45.2% |
| double_shared2 | 600 | 7.8 | 3.9 | 0.03 | 0.0% | 49.3% |
| single_distinct | 600 | 7.7 | 3.9 | 0.03 | 0.0% | 67.2% |
| double_distinct | 600 | 7.7 | 3.8 | 0.03 | 0.0% | 47.3% |
| gwas_only | 600 | 7.8 | NA | NA | 0.0% | 99.8% |
| null_both | 600 | NA | NA | NA | 99.7% | 99.8% |

## Figures

- Main figure: `matched_fpr_by_scenario.png`
- PP4 distribution by scenario and method: `pp4_distribution.png`
- `roc_test.png`

## Reproduce

```
Rscript scripts/bulk_gated_v1/run.R --stage main --cores 4 && Rscript scripts/bulk_gated_v1/evaluate.R --stage main
```

## Files changed (git diff --stat)

(Codex pastes the output of `git diff --stat HEAD~1` here after committing)
