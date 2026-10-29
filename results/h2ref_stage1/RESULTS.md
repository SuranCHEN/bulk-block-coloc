# Bulk block-coloc: main comparison (stage 1)

Generated 2026-10-29 03:28. Setting: GWAS N 20,000, local h2 0.3%; eQTL N 1,000, local h2 1.5%; each trait's h2 split equally over its causal SNPs (per-SNP h2 = h2 / K); causal |r| < 0.2.
LD given to the methods: 1000 Genomes EUR reference LD (489 people) for both traits.
Calibration windows: 31263287, 37199593, 39694225, 42313324, 44324730, 46518082. Test windows: 21059630, 24009444, 26904965, 29033189, 33343871, 34894283.
Loci: 3600 calibration, 3600 test (per method). Target FPR 5% on calibration windows (negatives: single_distinct and double_distinct).
Extra negatives gwas_only and null_both are reported for the false-positive rate only (not calibrated, not pooled).

## Three-line conclusion

1. **检出率**：误报率在校准窗口上匹配到 5%（测试窗口实际为 3.4%–4.6%）时，block-coloc（条件 lead 选择，本文方法；没有合格的块时退回整区域 PP4，2026-10-28 定为正式规则）检出 76.8% 的共定位位点（误报率 3.8%，AUC 0.934），比整区域 coloc.abf（73.2%）高 3.5 个百分点（同一批位点配对，95% CI +2.7 至 +4.4，p < 1e-15），比 coloc.susie + coloc.abf fallback（72.6%）高 4.2 个百分点，与 coloc.signals (cond)（76.2%，+0.6，p = 0.38）相当，低于 coloc.signals (mask)（78.2%，−1.5，p = 0.007）。不含退回规则的版本（无合格块记为无结果）检出 73.8%，与整区域 coloc.abf 无差别（+0.5，p = 0.5），作为敏感性分析报告；配对比较见 4b 节。
2. **为什么提升有限**：h² 按因果 SNP 个数平分后，eQTL 的每个信号只有 |z| ≈ 2.75，在选定的 lead 阈值 1e-3 下，多信号场景里 eQTL 只有 15%–17% 能找到 2 个 lead（与理论预测 14.9% 一致）；切块主要由 GWAS 的两个 lead 触发，而 eQTL 唯一的 lead 常和 GWAS 的 lead 不在同一个块里，于是没有合格的块、退回整区域 PP4（总体 22.9%，double_shared1 35.5%，double_distinct 63.8%）。阈值 1e-3 下单信号位点也有 15.5%（single_shared）走切块，来自噪声造成的额外 lead（无信号区域出现 lead 的理论概率 13.4%）；同样的规则在阈值 1e-4 下只有 1.5%，两个阈值在校准窗口上的检出率相近，按事先定的规则取 1e-3。
3. **最难的场景与额外阴性**：double_shared1（两个性状各有一个不共享的信号）上所有方法都变差（block-coloc 54.8%，整区域 coloc.abf 45.3%，coloc.signals (mask) 61.8%）；共享变异在两个性状里都不是较强的那个时，整区域 coloc.abf 只检出 11.8%，block-coloc 20.7%，coloc.signals 36.7%。只有 GWAS 有信号时（gwas_only），block-coloc 误报 4.8%，整区域 coloc.abf 5.0%，coloc.signals 2.5%–3.5%，coloc.susie 0%（不含退回规则的版本为 7.7%，来自其较低的判定线）；两性状都没有信号时（null_both）所有方法均为 0%。coloc.susie 在 eQTL 上 62%–89% 的位点找不到 credible set，单独使用时检出率只有 18.3%。

## 1. Main table (test windows, frozen cut-offs)

n_positive = loci called colocalised (PP4 above the method's cut-off). For shared scenarios rate is the detection rate;
for not-shared scenarios it is the false-positive rate. Denominators include every locus; no result counts as not called.
'not shared (extra)' rows (gwas_only, null_both) use the same frozen cut-offs but were not used to set them.

| scenario | truth | method | lead_p | n_total | n_positive | rate | 95% CI | no_call_rate | branch_block |
|---|---|---|---|---|---|---|---|---|---|
| single_shared | shared | block-coloc (conditional lead selection, proposed) | 1e-03 | 600 | 577 | 96.2% | [94.3%, 97.4%] | 0.0% | 15.5% |
| single_shared | shared | block-coloc (no fallback) | 1e-03 | 600 | 577 | 96.2% | [94.3%, 97.4%] | 1.2% | 15.5% |
| single_shared | shared | block-coloc (marginal lead selection) | 1e-05 | 600 | 578 | 96.3% | [94.5%, 97.6%] | 0.0% | 4.0% |
| single_shared | shared | block-coloc (original handoff) | 1e-04 | 600 | 576 | 96.0% | [94.1%, 97.3%] | 0.3% | 86.7% |
| single_shared | shared | coloc.abf, whole region | - | 600 | 578 | 96.3% | [94.5%, 97.6%] | 0.0% | 0.0% |
| single_shared | shared | coloc.susie | - | 600 | 225 | 37.5% | [33.7%, 41.4%] | 62.3% | 0.0% |
| single_shared | shared | coloc.susie + coloc.abf fallback | - | 600 | 576 | 96.0% | [94.1%, 97.3%] | 0.0% | 0.0% |
| single_shared | shared | coloc.signals (mask) | - | 600 | 573 | 95.5% | [93.5%, 96.9%] | 0.0% | 0.0% |
| single_shared | shared | coloc.signals (cond) | - | 600 | 569 | 94.8% | [92.8%, 96.3%] | 0.0% | 0.0% |
| single_double_shared | shared | block-coloc (conditional lead selection, proposed) | 1e-03 | 600 | 445 | 74.2% | [70.5%, 77.5%] | 0.0% | 22.2% |
| single_double_shared | shared | block-coloc (no fallback) | 1e-03 | 600 | 446 | 74.3% | [70.7%, 77.7%] | 7.5% | 22.2% |
| single_double_shared | shared | block-coloc (marginal lead selection) | 1e-05 | 600 | 448 | 74.7% | [71.0%, 78.0%] | 0.5% | 2.2% |
| single_double_shared | shared | block-coloc (original handoff) | 1e-04 | 600 | 435 | 72.5% | [68.8%, 75.9%] | 3.8% | 50.0% |
| single_double_shared | shared | coloc.abf, whole region | - | 600 | 446 | 74.3% | [70.7%, 77.7%] | 0.0% | 0.0% |
| single_double_shared | shared | coloc.susie | - | 600 |  65 | 10.8% | [8.6%, 13.6%] | 85.7% | 0.0% |
| single_double_shared | shared | coloc.susie + coloc.abf fallback | - | 600 | 431 | 71.8% | [68.1%, 75.3%] | 0.0% | 0.0% |
| single_double_shared | shared | coloc.signals (mask) | - | 600 | 428 | 71.3% | [67.6%, 74.8%] | 0.0% | 0.0% |
| single_double_shared | shared | coloc.signals (cond) | - | 600 | 406 | 67.7% | [63.8%, 71.3%] | 0.0% | 0.0% |
| double_shared2 | shared | block-coloc (conditional lead selection, proposed) | 1e-03 | 600 | 491 | 81.8% | [78.5%, 84.7%] | 0.0% | 66.8% |
| double_shared2 | shared | block-coloc (no fallback) | 1e-03 | 600 | 471 | 78.5% | [75.0%, 81.6%] | 12.2% | 66.8% |
| double_shared2 | shared | block-coloc (marginal lead selection) | 1e-05 | 600 | 466 | 77.7% | [74.2%, 80.8%] | 0.8% | 12.0% |
| double_shared2 | shared | block-coloc (original handoff) | 1e-04 | 600 | 494 | 82.3% | [79.1%, 85.2%] | 7.2% | 80.3% |
| double_shared2 | shared | coloc.abf, whole region | - | 600 | 461 | 76.8% | [73.3%, 80.0%] | 0.0% | 0.0% |
| double_shared2 | shared | coloc.susie | - | 600 |  84 | 14.0% | [11.5%, 17.0%] | 85.3% | 0.0% |
| double_shared2 | shared | coloc.susie + coloc.abf fallback | - | 600 | 464 | 77.3% | [73.8%, 80.5%] | 0.0% | 0.0% |
| double_shared2 | shared | coloc.signals (mask) | - | 600 | 506 | 84.3% | [81.2%, 87.0%] | 0.2% | 0.0% |
| double_shared2 | shared | coloc.signals (cond) | - | 600 | 492 | 82.0% | [78.7%, 84.9%] | 0.0% | 0.0% |
| double_shared1 | shared | block-coloc (conditional lead selection, proposed) | 1e-03 | 600 | 329 | 54.8% | [50.8%, 58.8%] | 0.0% | 68.5% |
| double_shared1 | shared | block-coloc (no fallback) | 1e-03 | 600 | 276 | 46.0% | [42.1%, 50.0%] | 35.5% | 68.5% |
| double_shared1 | shared | block-coloc (marginal lead selection) | 1e-05 | 600 | 267 | 44.5% | [40.6%, 48.5%] | 8.5% | 14.5% |
| double_shared1 | shared | block-coloc (original handoff) | 1e-04 | 600 | 335 | 55.8% | [51.8%, 59.8%] | 11.5% | 56.7% |
| double_shared1 | shared | coloc.abf, whole region | - | 600 | 272 | 45.3% | [41.4%, 49.3%] | 0.0% | 0.0% |
| double_shared1 | shared | coloc.susie | - | 600 |  66 | 11.0% | [8.7%, 13.8%] | 83.5% | 0.0% |
| double_shared1 | shared | coloc.susie + coloc.abf fallback | - | 600 | 271 | 45.2% | [41.2%, 49.2%] | 0.0% | 0.0% |
| double_shared1 | shared | coloc.signals (mask) | - | 600 | 371 | 61.8% | [57.9%, 65.6%] | 0.0% | 0.0% |
| double_shared1 | shared | coloc.signals (cond) | - | 600 | 361 | 60.2% | [56.2%, 64.0%] | 0.0% | 0.0% |
| single_distinct | not shared | block-coloc (conditional lead selection, proposed) | 1e-03 | 600 |  22 | 3.7% | [2.4%, 5.5%] | 0.0% | 18.0% |
| single_distinct | not shared | block-coloc (no fallback) | 1e-03 | 600 |  32 | 5.3% | [3.8%, 7.4%] | 17.3% | 18.0% |
| single_distinct | not shared | block-coloc (marginal lead selection) | 1e-05 | 600 |  22 | 3.7% | [2.4%, 5.5%] | 4.8% | 5.0% |
| single_distinct | not shared | block-coloc (original handoff) | 1e-04 | 600 |  25 | 4.2% | [2.8%, 6.1%] | 3.0% | 6.2% |
| single_distinct | not shared | coloc.abf, whole region | - | 600 |  21 | 3.5% | [2.3%, 5.3%] | 0.0% | 0.0% |
| single_distinct | not shared | coloc.susie | - | 600 |  25 | 4.2% | [2.8%, 6.1%] | 64.2% | 0.0% |
| single_distinct | not shared | coloc.susie + coloc.abf fallback | - | 600 |  20 | 3.3% | [2.2%, 5.1%] | 0.0% | 0.0% |
| single_distinct | not shared | coloc.signals (mask) | - | 600 |  22 | 3.7% | [2.4%, 5.5%] | 0.0% | 0.0% |
| single_distinct | not shared | coloc.signals (cond) | - | 600 |  19 | 3.2% | [2.0%, 4.9%] | 0.0% | 0.0% |
| double_distinct | not shared | block-coloc (conditional lead selection, proposed) | 1e-03 | 600 |  23 | 3.8% | [2.6%, 5.7%] | 0.0% | 65.8% |
| double_distinct | not shared | block-coloc (no fallback) | 1e-03 | 600 |  20 | 3.3% | [2.2%, 5.1%] | 63.8% | 65.8% |
| double_distinct | not shared | block-coloc (marginal lead selection) | 1e-05 | 600 |  20 | 3.3% | [2.2%, 5.1%] | 9.8% | 10.0% |
| double_distinct | not shared | block-coloc (original handoff) | 1e-04 | 600 |  27 | 4.5% | [3.1%, 6.5%] | 24.5% | 30.8% |
| double_distinct | not shared | coloc.abf, whole region | - | 600 |  21 | 3.5% | [2.3%, 5.3%] | 0.0% | 0.0% |
| double_distinct | not shared | coloc.susie | - | 600 |  18 | 3.0% | [1.9%, 4.7%] | 89.0% | 0.0% |
| double_distinct | not shared | coloc.susie + coloc.abf fallback | - | 600 |  21 | 3.5% | [2.3%, 5.3%] | 0.0% | 0.0% |
| double_distinct | not shared | coloc.signals (mask) | - | 600 |  33 | 5.5% | [3.9%, 7.6%] | 0.0% | 0.0% |
| double_distinct | not shared | coloc.signals (cond) | - | 600 |  26 | 4.3% | [3.0%, 6.3%] | 0.0% | 0.0% |
| gwas_only | not shared (extra) | block-coloc (conditional lead selection, proposed) | 1e-03 | 600 |  29 | 4.8% | [3.4%, 6.9%] | 0.0% | 2.2% |
| gwas_only | not shared (extra) | block-coloc (no fallback) | 1e-03 | 600 |  46 | 7.7% | [5.8%, 10.1%] | 2.2% | 2.2% |
| gwas_only | not shared (extra) | block-coloc (marginal lead selection) | 1e-05 | 600 |  31 | 5.2% | [3.7%, 7.2%] | 0.2% | 0.2% |
| gwas_only | not shared (extra) | block-coloc (original handoff) | 1e-04 | 600 |  29 | 4.8% | [3.4%, 6.9%] | 11.7% | 20.8% |
| gwas_only | not shared (extra) | coloc.abf, whole region | - | 600 |  30 | 5.0% | [3.5%, 7.0%] | 0.0% | 0.0% |
| gwas_only | not shared (extra) | coloc.susie | - | 600 |   0 | 0.0% | [0.0%, 0.6%] | 100.0% | 0.0% |
| gwas_only | not shared (extra) | coloc.susie + coloc.abf fallback | - | 600 |  32 | 5.3% | [3.8%, 7.4%] | 0.0% | 0.0% |
| gwas_only | not shared (extra) | coloc.signals (mask) | - | 600 |  21 | 3.5% | [2.3%, 5.3%] | 0.0% | 0.0% |
| gwas_only | not shared (extra) | coloc.signals (cond) | - | 600 |  15 | 2.5% | [1.5%, 4.1%] | 0.0% | 0.0% |
| null_both | not shared (extra) | block-coloc (conditional lead selection, proposed) | 1e-03 | 600 |   0 | 0.0% | [0.0%, 0.6%] | 0.0% | 0.3% |
| null_both | not shared (extra) | block-coloc (no fallback) | 1e-03 | 600 |   0 | 0.0% | [0.0%, 0.6%] | 0.3% | 0.3% |
| null_both | not shared (extra) | block-coloc (marginal lead selection) | 1e-05 | 600 |   0 | 0.0% | [0.0%, 0.6%] | 0.0% | 0.0% |
| null_both | not shared (extra) | block-coloc (original handoff) | 1e-04 | 600 |   0 | 0.0% | [0.0%, 0.6%] | 26.7% | 33.0% |
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
| block-coloc (conditional lead selection, proposed) | 4.8% [3.4%, 6.9%] | 0.0% [0.0%, 0.6%] |
| block-coloc (no fallback) | 7.7% [5.8%, 10.1%] | 0.0% [0.0%, 0.6%] |
| block-coloc (marginal lead selection) | 5.2% [3.7%, 7.2%] | 0.0% [0.0%, 0.6%] |
| block-coloc (original handoff) | 4.8% [3.4%, 6.9%] | 0.0% [0.0%, 0.6%] |
| coloc.abf, whole region | 5.0% [3.5%, 7.0%] | 0.0% [0.0%, 0.6%] |
| coloc.susie | 0.0% [0.0%, 0.6%] | 0.0% [0.0%, 0.6%] |
| coloc.susie + coloc.abf fallback | 5.3% [3.8%, 7.4%] | 0.0% [0.0%, 0.6%] |
| coloc.signals (mask) | 3.5% [2.3%, 5.3%] | 0.0% [0.0%, 0.6%] |
| coloc.signals (cond) | 2.5% [1.5%, 4.1%] | 0.0% [0.0%, 0.6%] |

## 2. Pooled over scenarios (test windows)

TPR pools the shared scenarios (single_shared, single_double_shared, double_shared2, double_shared1); FPR pools single_distinct and double_distinct. AUC does not depend on any cut-off.

| method | lead_p | n_pos | TPR | n_neg | FPR | AUC | no_call | branch_block |
|---|---|---|---|---|---|---|---|---|
| block-coloc (conditional lead selection, proposed) | 1e-03 | 2400 | 76.8% [75.0%, 78.4%] | 1200 | 3.8% [2.8%, 5.0%] | 0.934 [0.926, 0.942] | 0.0% | 42.8% |
| block-coloc (no fallback) | 1e-03 | 2400 | 73.8% [72.0%, 75.5%] | 1200 | 4.3% [3.3%, 5.6%] | 0.864 [0.851, 0.875] | 22.9% | 42.8% |
| block-coloc (marginal lead selection) | 1e-05 | 2400 | 73.3% [71.5%, 75.0%] | 1200 | 3.5% [2.6%, 4.7%] | 0.914 [0.904, 0.922] | 4.1% | 7.9% |
| block-coloc (original handoff) | 1e-04 | 2400 | 76.7% [74.9%, 78.3%] | 1200 | 4.3% [3.3%, 5.6%] | 0.900 [0.889, 0.910] | 8.4% | 51.8% |
| coloc.abf, whole region | - | 2400 | 73.2% [71.4%, 74.9%] | 1200 | 3.5% [2.6%, 4.7%] | 0.922 [0.913, 0.930] | 0.0% | 0.0% |
| coloc.susie | - | 2400 | 18.3% [16.8%, 19.9%] | 1200 | 3.6% [2.7%, 4.8%] | 0.509 [0.494, 0.522] | 78.3% | 0.0% |
| coloc.susie + coloc.abf fallback | - | 2400 | 72.6% [70.8%, 74.3%] | 1200 | 3.4% [2.5%, 4.6%] | 0.914 [0.904, 0.923] | 0.0% | 0.0% |
| coloc.signals (mask) | - | 2400 | 78.2% [76.6%, 79.9%] | 1200 | 4.6% [3.5%, 5.9%] | 0.940 [0.933, 0.947] | 0.0% | 0.0% |
| coloc.signals (cond) | - | 2400 | 76.2% [74.4%, 77.8%] | 1200 | 3.8% [2.8%, 5.0%] | 0.943 [0.935, 0.950] | 0.0% | 0.0% |

## 3. Cut-offs chosen on calibration windows

lead_p = lead threshold (GWAS / eQTL when they differ). '*' marks the lead threshold kept for each block-coloc version.

| method | lead_p | cutoff | calib_FPR | calib_TPR |
|---|---|---|---|---|
| block-coloc (conditional lead selection, proposed) | 1e-03 * | 0.1448 | 5.0% | 76.5% |
| block-coloc (conditional lead selection, proposed) | 1e-04 | 0.1412 | 5.0% | 74.4% |
| block-coloc (conditional lead selection, proposed) | 1e-05 | 0.1412 | 5.0% | 73.3% |
| block-coloc (conditional lead selection, proposed) [earlier bulk thresholds] | 5e-08 / 1e-05 | 0.1412 | 5.0% | 72.9% |
| block-coloc (no fallback) | 1e-03 * | 0.1067 | 5.0% | 73.9% |
| block-coloc (no fallback) | 1e-04 | 0.1300 | 5.0% | 73.8% |
| block-coloc (no fallback) | 1e-05 | 0.1396 | 5.0% | 73.0% |
| block-coloc (no fallback) [earlier bulk thresholds] | 5e-08 / 1e-05 | 0.1412 | 5.0% | 72.8% |
| block-coloc (marginal lead selection) | 1e-03 | 0.1210 | 5.0% | 70.7% |
| block-coloc (marginal lead selection) | 1e-04 | 0.1344 | 5.0% | 72.6% |
| block-coloc (marginal lead selection) | 1e-05 * | 0.1396 | 5.0% | 73.0% |
| block-coloc (marginal lead selection) [earlier bulk thresholds] | 5e-08 / 1e-05 | 0.1396 | 5.0% | 73.1% |
| block-coloc (original handoff) | 1e-03 | 0.1599 | 5.0% | 74.6% |
| block-coloc (original handoff) | 1e-04 * | 0.1416 | 5.0% | 76.1% |
| block-coloc (original handoff) | 1e-05 | 0.1431 | 5.0% | 75.8% |
| block-coloc (original handoff) [earlier bulk thresholds] | 5e-08 / 1e-05 | 0.1431 | 5.0% | 73.8% |
| block-coloc (marginal lead selection) + whole-region fallback | 1e-03 | 0.1683 | 5.0% | 74.6% |
| block-coloc (marginal lead selection) + whole-region fallback | 1e-04 | 0.1471 | 5.0% | 74.3% |
| block-coloc (marginal lead selection) + whole-region fallback | 1e-05 * | 0.1416 | 5.0% | 73.3% |
| block-coloc (marginal lead selection) + whole-region fallback [earlier bulk thresholds] | 5e-08 / 1e-05 | 0.1412 | 5.0% | 73.1% |
| block-coloc (original handoff) + whole-region fallback | 1e-03 | 0.1797 | 5.0% | 75.9% |
| block-coloc (original handoff) + whole-region fallback | 1e-04 * | 0.1591 | 5.0% | 76.5% |
| block-coloc (original handoff) + whole-region fallback | 1e-05 | 0.1504 | 5.0% | 76.7% |
| block-coloc (original handoff) + whole-region fallback [earlier bulk thresholds] | 5e-08 / 1e-05 | 0.1435 | 5.0% | 74.6% |
| coloc.abf, whole region | - | 0.1412 | 5.0% | 72.5% |
| coloc.susie | - | 0.0114 | 5.0% | 18.5% |
| coloc.susie + coloc.abf fallback | - | 0.1389 | 5.0% | 72.9% |
| coloc.signals (mask) | - | 0.1715 | 5.0% | 77.6% |
| coloc.signals (cond) | - | 0.2103 | 5.0% | 75.5% |

## 4. Paired differences on the same test loci (McNemar)

Difference = method minus reference, in percentage points, with 95% CI and exact McNemar p-value.

| method | vs | TPR_difference | FPR_difference |
|---|---|---|---|
| block-coloc (conditional lead selection, proposed) | coloc.susie + coloc.abf fallback | +4.2 pts [+3.2, +5.1], p=<2e-16 | +0.3 pts [-0.1, +0.7], p=0.22 |
| block-coloc (no fallback) | coloc.susie + coloc.abf fallback | +1.2 pts [-0.3, +2.6], p=0.14 | +0.9 pts [-0.1, +1.9], p=0.099 |
| block-coloc (marginal lead selection) | coloc.susie + coloc.abf fallback | +0.7 pts [+0.0, +1.4], p=0.053 | +0.1 pts [-0.2, +0.4], p=1 |
| block-coloc (original handoff) | coloc.susie + coloc.abf fallback | +4.1 pts [+2.9, +5.3], p=3.6e-11 | +0.9 pts [+0.1, +1.7], p=0.043 |
| coloc.abf, whole region | coloc.susie + coloc.abf fallback | +0.6 pts [-0.0, +1.3], p=0.082 | +0.1 pts [-0.1, +0.2], p=1 |
| coloc.susie | coloc.susie + coloc.abf fallback | -54.2 pts [-56.4, -52.1], p=<2e-16 | +0.2 pts [-1.3, +1.7], p=0.91 |
| coloc.signals (mask) | coloc.susie + coloc.abf fallback | +5.7 pts [+4.4, +6.9], p=<2e-16 | +1.2 pts [+0.4, +2.0], p=0.0066 |
| coloc.signals (cond) | coloc.susie + coloc.abf fallback | +3.6 pts [+2.2, +5.0], p=1.1e-06 | +0.3 pts [-0.7, +1.4], p=0.64 |
| block-coloc (conditional lead selection, proposed) | coloc.abf, whole region | +3.5 pts [+2.7, +4.4], p=<2e-16 | +0.2 pts [-0.1, +0.6], p=0.37 |
| block-coloc (no fallback) | coloc.abf, whole region | +0.5 pts [-0.9, +2.0], p=0.5 | +0.8 pts [-0.1, +1.8], p=0.13 |
| block-coloc (marginal lead selection) | coloc.abf, whole region | +0.1 pts [-0.4, +0.6], p=0.87 | +0.0 pts [-0.2, +0.2], p=1 |
| block-coloc (original handoff) | coloc.abf, whole region | +3.5 pts [+2.3, +4.6], p=2e-09 | +0.8 pts [+0.0, +1.6], p=0.064 |
| coloc.susie | coloc.abf, whole region | -54.9 pts [-57.0, -52.7], p=<2e-16 | +0.1 pts [-1.4, +1.6], p=1 |
| coloc.susie + coloc.abf fallback | coloc.abf, whole region | -0.6 pts [-1.3, +0.0], p=0.082 | -0.1 pts [-0.2, +0.1], p=1 |
| coloc.signals (mask) | coloc.abf, whole region | +5.0 pts [+3.8, +6.3], p=1.1e-15 | +1.1 pts [+0.3, +1.9], p=0.011 |
| coloc.signals (cond) | coloc.abf, whole region | +3.0 pts [+1.6, +4.4], p=4.5e-05 | +0.2 pts [-0.8, +1.3], p=0.75 |

### 4b. Adopted block-coloc rule and the no-fallback sensitivity analysis

The proposed method uses whole-region PP4 when no valid block is available; 'block-coloc (no fallback)' keeps those loci as no result (section 8).

| method | vs | TPR_difference | FPR_difference |
|---|---|---|---|
| block-coloc (conditional lead selection, proposed) | coloc.abf, whole region | +3.5 pts [+2.7, +4.4], p=<2e-16 | +0.2 pts [-0.1, +0.6], p=0.37 |
| block-coloc (conditional lead selection, proposed) | coloc.susie + coloc.abf fallback | +4.2 pts [+3.2, +5.1], p=<2e-16 | +0.3 pts [-0.1, +0.7], p=0.22 |
| block-coloc (conditional lead selection, proposed) | coloc.signals (cond) | +0.6 pts [-0.6, +1.8], p=0.38 | +0.0 pts [-1.0, +1.0], p=1 |
| block-coloc (conditional lead selection, proposed) | coloc.signals (mask) | -1.5 pts [-2.6, -0.4], p=0.0074 | -0.8 pts [-1.6, -0.0], p=0.064 |
| block-coloc (no fallback) | coloc.signals (cond) | -2.4 pts [-4.0, -0.8], p=0.0037 | +0.6 pts [-0.6, +1.8], p=0.43 |
| block-coloc (no fallback) | coloc.signals (mask) | -4.5 pts [-6.0, -3.0], p=1.5e-08 | -0.2 pts [-1.4, +0.9], p=0.78 |

## 5. Branch choice: split into blocks or not (test windows)

should_go: block when either trait has 2 causal SNPs, whole otherwise. leads_match_truth: the method's lead counts
(GWAS, eQTL) equal the true numbers of causal SNPs.

| method | scenario | true_signals | should_go | n | branch_block | correct_branch | leads_match_truth | median_leads |
|---|---|---|---|---|---|---|---|---|
| block-coloc (conditional lead selection, proposed) | single_shared | 1+1 | whole | 600 | 15.5% | 84.5% | 69.2% | 1+1 |
| block-coloc (conditional lead selection, proposed) | single_double_shared | 1+2 | block | 600 | 22.2% | 22.2% | 14.2% | 1+1 |
| block-coloc (conditional lead selection, proposed) | double_shared2 | 2+2 | block | 600 | 66.8% | 66.8% | 13.5% | 2+1 |
| block-coloc (conditional lead selection, proposed) | double_shared1 | 2+2 | block | 600 | 68.5% | 68.5% | 11.3% | 2+1 |
| block-coloc (conditional lead selection, proposed) | single_distinct | 1+1 | whole | 600 | 18.0% | 82.0% | 66.2% | 1+1 |
| block-coloc (conditional lead selection, proposed) | double_distinct | 2+2 | block | 600 | 65.8% | 65.8% | 13.0% | 2+1 |
| block-coloc (conditional lead selection, proposed) | gwas_only | 1+0 | whole | 600 | 2.2% | 97.8% | 76.7% | 1+0 |
| block-coloc (conditional lead selection, proposed) | null_both | 0+0 | whole | 600 | 0.3% | 99.7% | 74.5% | 0+0 |
| block-coloc (no fallback) | single_shared | 1+1 | whole | 600 | 15.5% | 84.5% | 69.2% | 1+1 |
| block-coloc (no fallback) | single_double_shared | 1+2 | block | 600 | 22.2% | 22.2% | 14.2% | 1+1 |
| block-coloc (no fallback) | double_shared2 | 2+2 | block | 600 | 66.8% | 66.8% | 13.5% | 2+1 |
| block-coloc (no fallback) | double_shared1 | 2+2 | block | 600 | 68.5% | 68.5% | 11.3% | 2+1 |
| block-coloc (no fallback) | single_distinct | 1+1 | whole | 600 | 18.0% | 82.0% | 66.2% | 1+1 |
| block-coloc (no fallback) | double_distinct | 2+2 | block | 600 | 65.8% | 65.8% | 13.0% | 2+1 |
| block-coloc (no fallback) | gwas_only | 1+0 | whole | 600 | 2.2% | 97.8% | 76.7% | 1+0 |
| block-coloc (no fallback) | null_both | 0+0 | whole | 600 | 0.3% | 99.7% | 74.5% | 0+0 |
| block-coloc (marginal lead selection) | single_shared | 1+1 | whole | 600 | 4.0% | 96.0% | 36.5% | 1+0 |
| block-coloc (marginal lead selection) | single_double_shared | 1+2 | block | 600 | 2.2% | 2.2% | 0.8% | 1+0 |
| block-coloc (marginal lead selection) | double_shared2 | 2+2 | block | 600 | 12.0% | 12.0% | 0.7% | 2+0 |
| block-coloc (marginal lead selection) | double_shared1 | 2+2 | block | 600 | 14.5% | 14.5% | 0.5% | 2+0 |
| block-coloc (marginal lead selection) | single_distinct | 1+1 | whole | 600 | 5.0% | 95.0% | 33.5% | 1+0 |
| block-coloc (marginal lead selection) | double_distinct | 2+2 | block | 600 | 10.0% | 10.0% | 0.8% | 2+0 |
| block-coloc (marginal lead selection) | gwas_only | 1+0 | whole | 600 | 0.2% | 99.8% | 89.3% | 1+0 |
| block-coloc (marginal lead selection) | null_both | 0+0 | whole | 600 | 0.0% | 100.0% | 99.5% | 0+0 |
| block-coloc (original handoff) | single_shared | 1+1 | whole | 600 | 86.7% | 13.3% | 59.8% | 1+1 |
| block-coloc (original handoff) | single_double_shared | 1+2 | block | 600 | 50.0% | 50.0% | 17.8% | 1+1 |
| block-coloc (original handoff) | double_shared2 | 2+2 | block | 600 | 80.3% | 80.3% | 17.5% | 2+1 |
| block-coloc (original handoff) | double_shared1 | 2+2 | block | 600 | 56.7% | 56.7% | 19.8% | 2+1 |
| block-coloc (original handoff) | single_distinct | 1+1 | whole | 600 | 6.2% | 93.8% | 59.0% | 1+1 |
| block-coloc (original handoff) | double_distinct | 2+2 | block | 600 | 30.8% | 30.8% | 15.7% | 2+1 |
| block-coloc (original handoff) | gwas_only | 1+0 | whole | 600 | 20.8% | 79.2% | 0.0% | 1+2 |
| block-coloc (original handoff) | null_both | 0+0 | whole | 600 | 33.0% | 67.0% | 0.0% | 2+2 |

## 6. Detection when the shared signal is not the stronger one (test windows)

single_double_shared: the eQTL has two causal variants, one shared with the GWAS and one eQTL-only. 'stronger' compares their observed eQTL |z|.

| method | n_shared_stronger | TPR_shared_stronger | n_other_stronger | TPR_other_stronger |
|---|---|---|---|---|
| block-coloc (conditional lead selection, proposed) | 300 | 95.0% [91.9%, 96.9%] | 300 | 53.3% [47.7%, 58.9%] |
| block-coloc (no fallback) | 300 | 95.7% [92.7%, 97.5%] | 300 | 53.0% [47.4%, 58.6%] |
| block-coloc (marginal lead selection) | 300 | 95.3% [92.3%, 97.2%] | 300 | 54.0% [48.3%, 59.6%] |
| block-coloc (original handoff) | 300 | 94.0% [90.7%, 96.2%] | 300 | 51.0% [45.4%, 56.6%] |
| coloc.abf, whole region | 300 | 95.3% [92.3%, 97.2%] | 300 | 53.3% [47.7%, 58.9%] |
| coloc.susie | 300 | 14.7% [11.1%, 19.1%] | 300 | 7.0% [4.6%, 10.5%] |
| coloc.susie + coloc.abf fallback | 300 | 94.0% [90.7%, 96.2%] | 300 | 49.7% [44.0%, 55.3%] |
| coloc.signals (mask) | 300 | 93.0% [89.5%, 95.4%] | 300 | 49.7% [44.0%, 55.3%] |
| coloc.signals (cond) | 300 | 91.0% [87.2%, 93.7%] | 300 | 44.3% [38.8%, 50.0%] |

double_shared1: each trait has one shared and one trait-specific causal variant. Columns: in how many traits the shared variant has the larger observed |z|.

| method | shared_stronger_in_both | in_one_trait | in_neither |
|---|---|---|---|
| block-coloc (conditional lead selection, proposed) | 91.4% [85.5%, 95.0%] (n=139) | 57.2% [51.5%, 62.7%] (n=292) | 20.7% [15.3%, 27.4%] (n=169) |
| block-coloc (no fallback) | 79.9% [72.4%, 85.7%] (n=139) | 46.2% [40.6%, 52.0%] (n=292) | 17.8% [12.7%, 24.2%] (n=169) |
| block-coloc (marginal lead selection) | 88.5% [82.1%, 92.8%] (n=139) | 42.1% [36.6%, 47.9%] (n=292) | 12.4% [8.3%, 18.2%] (n=169) |
| block-coloc (original handoff) | 92.1% [86.4%, 95.5%] (n=139) | 59.2% [53.5%, 64.7%] (n=292) | 20.1% [14.8%, 26.8%] (n=169) |
| coloc.abf, whole region | 91.4% [85.5%, 95.0%] (n=139) | 42.8% [37.3%, 48.5%] (n=292) | 11.8% [7.8%, 17.6%] (n=169) |
| coloc.susie | 15.8% [10.7%, 22.8%] (n=139) | 10.3% [7.3%, 14.3%] (n=292) | 8.3% [5.0%, 13.4%] (n=169) |
| coloc.susie + coloc.abf fallback | 92.8% [87.3%, 96.0%] (n=139) | 42.1% [36.6%, 47.9%] (n=292) | 11.2% [7.3%, 16.9%] (n=169) |
| coloc.signals (mask) | 90.6% [84.7%, 94.5%] (n=139) | 62.7% [57.0%, 68.0%] (n=292) | 36.7% [29.8%, 44.2%] (n=169) |
| coloc.signals (cond) | 87.8% [81.3%, 92.2%] (n=139) | 60.6% [54.9%, 66.0%] (n=292) | 36.7% [29.8%, 44.2%] (n=169) |

## 7. Localisation among detected test positives

Sub-region = best block (block branch) or the two credible sets of the best SuSiE signal pair.
Whole-locus calls and coloc.signals calls have no sub-region and are not counted here.

| method | detected | with_subregion | contains_shared | median_size_snps |
|---|---|---|---|---|
| block-coloc (conditional lead selection, proposed) | 1842 |  684 | 99.4% [98.5%, 99.8%] | 46 |
| block-coloc (no fallback) | 1770 |  690 | 99.3% [98.3%, 99.7%] | 46 |
| block-coloc (marginal lead selection) | 1759 |  136 | 98.5% [94.8%, 99.6%] | 56 |
| block-coloc (original handoff) | 1840 | 1456 | 98.8% [98.1%, 99.3%] | 48 |
| coloc.abf, whole region | 1757 |    0 | NA | NA |
| coloc.susie |  440 |  440 | 98.9% [97.4%, 99.5%] | 38 |
| coloc.susie + coloc.abf fallback | 1742 |  393 | 99.5% [98.2%, 99.9%] | 37 |
| coloc.signals (mask) | 1878 |    0 | NA | NA |
| coloc.signals (cond) | 1828 |    0 | NA | NA |

## 8. Sensitivity: fallback policy and the earlier bulk lead thresholds (test windows)

The proposed conditional method includes whole-region fallback; 'block-coloc (no fallback)' omits it.
'+ whole-region fallback' = when the block branch builds no valid block, use the whole-region PP4 instead of no result.
'[earlier bulk thresholds]' = lead thresholds GWAS 5e-8, eQTL 1e-5 (August bulk setting), each with its own calibrated cut-off.

| method | lead_p | cutoff | n_pos | TPR | n_neg | FPR | AUC | no_call | branch_block |
|---|---|---|---|---|---|---|---|---|---|
| block-coloc (conditional lead selection, proposed) | 1e-03 | 0.1448 | 2400 | 76.8% [75.0%, 78.4%] | 1200 | 3.8% [2.8%, 5.0%] | 0.934 | 0.0% | 42.8% |
| block-coloc (no fallback) | 1e-03 | 0.1067 | 2400 | 73.8% [72.0%, 75.5%] | 1200 | 4.3% [3.3%, 5.6%] | 0.864 | 22.9% | 42.8% |
| block-coloc (marginal lead selection) | 1e-05 | 0.1396 | 2400 | 73.3% [71.5%, 75.0%] | 1200 | 3.5% [2.6%, 4.7%] | 0.914 | 4.1% | 7.9% |
| block-coloc (original handoff) | 1e-04 | 0.1416 | 2400 | 76.7% [74.9%, 78.3%] | 1200 | 4.3% [3.3%, 5.6%] | 0.900 | 8.4% | 51.8% |
| block-coloc (marginal lead selection) + whole-region fallback | 1e-05 | 0.1416 | 2400 | 74.0% [72.3%, 75.8%] | 1200 | 3.6% [2.7%, 4.8%] | 0.926 | 0.0% | 7.9% |
| block-coloc (original handoff) + whole-region fallback | 1e-04 | 0.1591 | 2400 | 77.3% [75.6%, 78.9%] | 1200 | 4.6% [3.5%, 5.9%] | 0.934 | 0.0% | 51.8% |
| block-coloc (original handoff) [earlier bulk thresholds] | 5e-08 / 1e-05 | 0.1431 | 2400 | 74.7% [72.9%, 76.4%] | 1200 | 3.7% [2.7%, 4.9%] | 0.916 | 2.6% | 39.5% |
| block-coloc (marginal lead selection) [earlier bulk thresholds] | 5e-08 / 1e-05 | 0.1396 | 2400 | 73.4% [71.6%, 75.1%] | 1200 | 3.4% [2.5%, 4.6%] | 0.919 | 1.7% | 3.2% |
| block-coloc (no fallback) [earlier bulk thresholds] | 5e-08 / 1e-05 | 0.1412 | 2400 | 73.2% [71.4%, 74.9%] | 1200 | 3.4% [2.5%, 4.6%] | 0.919 | 1.0% | 1.8% |
| block-coloc (conditional lead selection, proposed) [earlier bulk thresholds] | 5e-08 / 1e-05 | 0.1412 | 2400 | 73.3% [71.5%, 75.1%] | 1200 | 3.5% [2.6%, 4.7%] | 0.922 | 0.0% | 1.8% |

## 9. Inputs (test windows)

Expected causal |z| (sqrt(N q / (1 - q)), q = h2 / K): GWAS 7.8 (K = 1) / 5.5 (K = 2); eQTL 3.9 (K = 1) / 2.7 (K = 2). Table: medians over loci and causal SNPs.

| scenario | n | causal_z_gwas | causal_z_eqtl | pair_abs_r_median | gwas_no_cs | eqtl_no_cs |
|---|---|---|---|---|---|---|
| single_shared | 600 | 7.8 | 4.0 | NA | 0.2% | 62.3% |
| single_double_shared | 600 | 7.8 | 2.8 | 0.03 | 0.0% | 85.7% |
| double_shared2 | 600 | 5.5 | 2.7 | 0.03 | 3.8% | 84.7% |
| double_shared1 | 600 | 5.5 | 2.8 | 0.03 | 3.7% | 83.2% |
| single_distinct | 600 | 7.7 | 4.0 | 0.03 | 0.0% | 64.2% |
| double_distinct | 600 | 5.4 | 2.7 | 0.03 | 3.3% | 88.8% |
| gwas_only | 600 | 7.8 | NA | NA | 0.3% | 99.8% |
| null_both | 600 | NA | NA | NA | 99.7% | 99.8% |

## Figures

- Main figure: `matched_fpr_by_scenario.png`
- PP4 distribution by scenario and method: `pp4_distribution.png`
- `roc_test.png`

## Reproduce

```
Rscript scripts/bulk_h2ref_v1/run.R --stage main --cores 4 && Rscript scripts/bulk_h2ref_v1/evaluate.R --stage main
```

## Files changed (git diff --stat)

(Codex pastes the output of `git diff --stat HEAD~1` here after committing)
