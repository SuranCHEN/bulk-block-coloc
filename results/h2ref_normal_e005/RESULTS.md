# Bulk block-coloc: main comparison (stage 1)

Generated 2026-10-29 03:29. Setting: GWAS N 20,000, local h2 0.3%; eQTL N 1,000, local h2 5%; the effect of each causal SNP is drawn from N(0, h2 / K), independently for the two traits, without rescaling (Var(g) = h2 in expectation only; sensitivity analysis); causal |r| < 0.2.
LD given to the methods: 1000 Genomes EUR reference LD (489 people) for both traits.
Calibration windows: 31263287, 37199593, 39694225, 42313324, 44324730, 46518082. Test windows: 21059630, 24009444, 26904965, 29033189, 33343871, 34894283.
Loci: 1800 calibration, 1800 test (per method). Target FPR 5% on calibration windows (negatives: single_distinct and double_distinct).
Extra negatives gwas_only and null_both are reported for the false-positive rate only (not calibrated, not pooled).

## Three-line conclusion

1. **检出率**：正态效应、eQTL h² 5%（50 次重复）、误报率在本设计的校准窗口上匹配到 5%（测试窗口实际 7.0%–8.7%，所有方法相近）时，block-coloc 检出 51.7%，比整区域 coloc.abf（46.2%）高 5.5 个百分点（配对 95% CI +4.2 至 +6.8，p < 1e-15，6 个测试窗口全部为正，+4.0 至 +7.0），比 coloc.susie + fallback 高 2.8 个百分点，与 coloc.signals (cond)（+0.8，p = 0.42）和 (mask)（+0.3，p = 0.74）无差别；AUC 0.786，为所有方法中最高（coloc.abf 0.755，signals 0.750–0.763）。
2. **与固定效应的比较**：固定效应下同一强度的增益为 +13.7 个百分点（第二段，沿用第一段判定线），正态效应下为 +5.5；方向不变，幅度减小。原因是正态效应下两个 eQTL 因果变异常常一强一弱，两因果位点中只有 31.1% 找到两个条件独立的 lead（固定效应 89%），切块比例 31.4%。整区域 coloc.abf 在共享信号较弱的位点上最差（double_shared1 中共享变异在两性状都不是较强者时 2.8%，block-coloc 11.3%）。
3. **额外阴性**：gwas_only 误报率 block-coloc 11.7%、整区域 coloc.abf 12.0%，高于 coloc.signals（2.7%–4.3%）：本设计的判定线（0.070）低于只有 GWAS 信号时 PP4 的典型值（约 0.09，Model 式 5.3），这对所有基于整区域 coloc.abf 的方法同样成立；null_both 所有方法为 0%。

## 1. Main table (test windows, frozen cut-offs)

n_positive = loci called colocalised (PP4 above the method's cut-off). For shared scenarios rate is the detection rate;
for not-shared scenarios it is the false-positive rate. Denominators include every locus; no result counts as not called.
'not shared (extra)' rows (gwas_only, null_both) use the same frozen cut-offs but were not used to set them.

| scenario | truth | method | lead_p | n_total | n_positive | rate | 95% CI | no_call_rate | branch_block |
|---|---|---|---|---|---|---|---|---|---|
| single_shared | shared | block-coloc (conditional lead selection, proposed) | 1e-03 | 300 | 193 | 64.3% | [58.8%, 69.5%] | 0.0% | 11.0% |
| single_shared | shared | block-coloc (no fallback) | 1e-03 | 300 | 193 | 64.3% | [58.8%, 69.5%] | 1.0% | 11.0% |
| single_shared | shared | block-coloc (marginal lead selection) | 1e-04 | 300 | 189 | 63.0% | [57.4%, 68.3%] | 1.0% | 26.7% |
| single_shared | shared | block-coloc (original handoff) | 1e-04 | 300 | 186 | 62.0% | [56.4%, 67.3%] | 6.7% | 56.7% |
| single_shared | shared | coloc.abf, whole region | - | 300 | 193 | 64.3% | [58.8%, 69.5%] | 0.0% | 0.0% |
| single_shared | shared | coloc.susie | - | 300 | 111 | 37.0% | [31.7%, 42.6%] | 63.0% | 0.0% |
| single_shared | shared | coloc.susie + coloc.abf fallback | - | 300 | 192 | 64.0% | [58.4%, 69.2%] | 0.0% | 0.0% |
| single_shared | shared | coloc.signals (mask) | - | 300 | 183 | 61.0% | [55.4%, 66.3%] | 0.0% | 0.0% |
| single_shared | shared | coloc.signals (cond) | - | 300 | 173 | 57.7% | [52.0%, 63.1%] | 0.0% | 0.0% |
| single_double_shared | shared | block-coloc (conditional lead selection, proposed) | 1e-03 | 300 | 127 | 42.3% | [36.9%, 48.0%] | 0.0% | 26.7% |
| single_double_shared | shared | block-coloc (no fallback) | 1e-03 | 300 | 123 | 41.0% | [35.6%, 46.6%] | 10.3% | 26.7% |
| single_double_shared | shared | block-coloc (marginal lead selection) | 1e-04 | 300 | 122 | 40.7% | [35.3%, 46.3%] | 11.3% | 30.3% |
| single_double_shared | shared | block-coloc (original handoff) | 1e-04 | 300 | 130 | 43.3% | [37.8%, 49.0%] | 8.0% | 43.7% |
| single_double_shared | shared | coloc.abf, whole region | - | 300 | 113 | 37.7% | [32.4%, 43.3%] | 0.0% | 0.0% |
| single_double_shared | shared | coloc.susie | - | 300 |  65 | 21.7% | [17.4%, 26.7%] | 65.0% | 0.0% |
| single_double_shared | shared | coloc.susie + coloc.abf fallback | - | 300 | 119 | 39.7% | [34.3%, 45.3%] | 0.0% | 0.0% |
| single_double_shared | shared | coloc.signals (mask) | - | 300 | 128 | 42.7% | [37.2%, 48.3%] | 0.0% | 0.0% |
| single_double_shared | shared | coloc.signals (cond) | - | 300 | 131 | 43.7% | [38.2%, 49.3%] | 0.0% | 0.0% |
| double_shared2 | shared | block-coloc (conditional lead selection, proposed) | 1e-03 | 300 | 183 | 61.0% | [55.4%, 66.3%] | 0.0% | 47.0% |
| double_shared2 | shared | block-coloc (no fallback) | 1e-03 | 300 | 172 | 57.3% | [51.7%, 62.8%] | 9.7% | 47.0% |
| double_shared2 | shared | block-coloc (marginal lead selection) | 1e-04 | 300 | 171 | 57.0% | [51.3%, 62.5%] | 7.3% | 40.7% |
| double_shared2 | shared | block-coloc (original handoff) | 1e-04 | 300 | 178 | 59.3% | [53.7%, 64.7%] | 6.0% | 55.0% |
| double_shared2 | shared | coloc.abf, whole region | - | 300 | 153 | 51.0% | [45.4%, 56.6%] | 0.0% | 0.0% |
| double_shared2 | shared | coloc.susie | - | 300 | 101 | 33.7% | [28.6%, 39.2%] | 57.0% | 0.0% |
| double_shared2 | shared | coloc.susie + coloc.abf fallback | - | 300 | 167 | 55.7% | [50.0%, 61.2%] | 0.0% | 0.0% |
| double_shared2 | shared | coloc.signals (mask) | - | 300 | 175 | 58.3% | [52.7%, 63.8%] | 0.0% | 0.0% |
| double_shared2 | shared | coloc.signals (cond) | - | 300 | 179 | 59.7% | [54.0%, 65.1%] | 0.0% | 0.0% |
| double_shared1 | shared | block-coloc (conditional lead selection, proposed) | 1e-03 | 300 | 117 | 39.0% | [33.7%, 44.6%] | 0.0% | 47.3% |
| double_shared1 | shared | block-coloc (no fallback) | 1e-03 | 300 | 112 | 37.3% | [32.1%, 42.9%] | 21.7% | 47.3% |
| double_shared1 | shared | block-coloc (marginal lead selection) | 1e-04 | 300 | 109 | 36.3% | [31.1%, 41.9%] | 19.0% | 44.3% |
| double_shared1 | shared | block-coloc (original handoff) | 1e-04 | 300 | 117 | 39.0% | [33.7%, 44.6%] | 8.3% | 43.0% |
| double_shared1 | shared | coloc.abf, whole region | - | 300 |  95 | 31.7% | [26.7%, 37.1%] | 0.0% | 0.0% |
| double_shared1 | shared | coloc.susie | - | 300 |  85 | 28.3% | [23.5%, 33.7%] | 52.7% | 0.0% |
| double_shared1 | shared | coloc.susie + coloc.abf fallback | - | 300 | 108 | 36.0% | [30.8%, 41.6%] | 0.0% | 0.0% |
| double_shared1 | shared | coloc.signals (mask) | - | 300 | 130 | 43.3% | [37.8%, 49.0%] | 0.0% | 0.0% |
| double_shared1 | shared | coloc.signals (cond) | - | 300 | 128 | 42.7% | [37.2%, 48.3%] | 0.0% | 0.0% |
| single_distinct | not shared | block-coloc (conditional lead selection, proposed) | 1e-03 | 300 |  25 | 8.3% | [5.7%, 12.0%] | 0.0% | 11.7% |
| single_distinct | not shared | block-coloc (no fallback) | 1e-03 | 300 |  25 | 8.3% | [5.7%, 12.0%] | 10.7% | 11.7% |
| single_distinct | not shared | block-coloc (marginal lead selection) | 1e-04 | 300 |  27 | 9.0% | [6.3%, 12.8%] | 18.3% | 22.0% |
| single_distinct | not shared | block-coloc (original handoff) | 1e-04 | 300 |  29 | 9.7% | [6.8%, 13.5%] | 14.3% | 23.0% |
| single_distinct | not shared | coloc.abf, whole region | - | 300 |  24 | 8.0% | [5.4%, 11.6%] | 0.0% | 0.0% |
| single_distinct | not shared | coloc.susie | - | 300 |   4 | 1.3% | [0.5%, 3.4%] | 70.3% | 0.0% |
| single_distinct | not shared | coloc.susie + coloc.abf fallback | - | 300 |  25 | 8.3% | [5.7%, 12.0%] | 0.0% | 0.0% |
| single_distinct | not shared | coloc.signals (mask) | - | 300 |  24 | 8.0% | [5.4%, 11.6%] | 0.0% | 0.0% |
| single_distinct | not shared | coloc.signals (cond) | - | 300 |  18 | 6.0% | [3.8%, 9.3%] | 0.0% | 0.0% |
| double_distinct | not shared | block-coloc (conditional lead selection, proposed) | 1e-03 | 300 |  19 | 6.3% | [4.1%, 9.7%] | 0.0% | 45.0% |
| double_distinct | not shared | block-coloc (no fallback) | 1e-03 | 300 |  17 | 5.7% | [3.6%, 8.9%] | 43.0% | 45.0% |
| double_distinct | not shared | block-coloc (marginal lead selection) | 1e-04 | 300 |  22 | 7.3% | [4.9%, 10.9%] | 36.7% | 40.7% |
| double_distinct | not shared | block-coloc (original handoff) | 1e-04 | 300 |  23 | 7.7% | [5.2%, 11.2%] | 18.0% | 25.7% |
| double_distinct | not shared | coloc.abf, whole region | - | 300 |  18 | 6.0% | [3.8%, 9.3%] | 0.0% | 0.0% |
| double_distinct | not shared | coloc.susie | - | 300 |  21 | 7.0% | [4.6%, 10.5%] | 56.3% | 0.0% |
| double_distinct | not shared | coloc.susie + coloc.abf fallback | - | 300 |  20 | 6.7% | [4.4%, 10.1%] | 0.0% | 0.0% |
| double_distinct | not shared | coloc.signals (mask) | - | 300 |  21 | 7.0% | [4.6%, 10.5%] | 0.0% | 0.0% |
| double_distinct | not shared | coloc.signals (cond) | - | 300 |  25 | 8.3% | [5.7%, 12.0%] | 0.0% | 0.0% |
| gwas_only | not shared (extra) | block-coloc (conditional lead selection, proposed) | 1e-03 | 300 |  35 | 11.7% | [8.5%, 15.8%] | 0.0% | 1.0% |
| gwas_only | not shared (extra) | block-coloc (no fallback) | 1e-03 | 300 |  38 | 12.7% | [9.4%, 16.9%] | 1.0% | 1.0% |
| gwas_only | not shared (extra) | block-coloc (marginal lead selection) | 1e-04 | 300 |  34 | 11.3% | [8.2%, 15.4%] | 0.0% | 0.0% |
| gwas_only | not shared (extra) | block-coloc (original handoff) | 1e-04 | 300 |  34 | 11.3% | [8.2%, 15.4%] | 20.0% | 31.3% |
| gwas_only | not shared (extra) | coloc.abf, whole region | - | 300 |  36 | 12.0% | [8.8%, 16.2%] | 0.0% | 0.0% |
| gwas_only | not shared (extra) | coloc.susie | - | 300 |   0 | 0.0% | [0.0%, 1.3%] | 100.0% | 0.0% |
| gwas_only | not shared (extra) | coloc.susie + coloc.abf fallback | - | 300 |  36 | 12.0% | [8.8%, 16.2%] | 0.0% | 0.0% |
| gwas_only | not shared (extra) | coloc.signals (mask) | - | 300 |  13 | 4.3% | [2.5%, 7.3%] | 0.0% | 0.0% |
| gwas_only | not shared (extra) | coloc.signals (cond) | - | 300 |   8 | 2.7% | [1.4%, 5.2%] | 0.0% | 0.0% |
| null_both | not shared (extra) | block-coloc (conditional lead selection, proposed) | 1e-03 | 300 |   0 | 0.0% | [0.0%, 1.3%] | 0.0% | 0.3% |
| null_both | not shared (extra) | block-coloc (no fallback) | 1e-03 | 300 |   0 | 0.0% | [0.0%, 1.3%] | 0.3% | 0.3% |
| null_both | not shared (extra) | block-coloc (marginal lead selection) | 1e-04 | 300 |   0 | 0.0% | [0.0%, 1.3%] | 0.0% | 0.0% |
| null_both | not shared (extra) | block-coloc (original handoff) | 1e-04 | 300 |   0 | 0.0% | [0.0%, 1.3%] | 24.0% | 30.7% |
| null_both | not shared (extra) | coloc.abf, whole region | - | 300 |   0 | 0.0% | [0.0%, 1.3%] | 0.0% | 0.0% |
| null_both | not shared (extra) | coloc.susie | - | 300 |   0 | 0.0% | [0.0%, 1.3%] | 100.0% | 0.0% |
| null_both | not shared (extra) | coloc.susie + coloc.abf fallback | - | 300 |   0 | 0.0% | [0.0%, 1.3%] | 0.0% | 0.0% |
| null_both | not shared (extra) | coloc.signals (mask) | - | 300 |   0 | 0.0% | [0.0%, 1.3%] | 0.0% | 0.0% |
| null_both | not shared (extra) | coloc.signals (cond) | - | 300 |   0 | 0.0% | [0.0%, 1.3%] | 0.0% | 0.0% |

## 1b. False-positive rate when at most one trait has a signal (test windows)

gwas_only: GWAS has one causal variant, the eQTL has none. null_both: neither trait has a signal.
Same frozen cut-offs as above; these loci were not used to set them.

| method | gwas_only | null_both |
|---|---|---|
| block-coloc (conditional lead selection, proposed) | 11.7% [8.5%, 15.8%] | 0.0% [0.0%, 1.3%] |
| block-coloc (no fallback) | 12.7% [9.4%, 16.9%] | 0.0% [0.0%, 1.3%] |
| block-coloc (marginal lead selection) | 11.3% [8.2%, 15.4%] | 0.0% [0.0%, 1.3%] |
| block-coloc (original handoff) | 11.3% [8.2%, 15.4%] | 0.0% [0.0%, 1.3%] |
| coloc.abf, whole region | 12.0% [8.8%, 16.2%] | 0.0% [0.0%, 1.3%] |
| coloc.susie | 0.0% [0.0%, 1.3%] | 0.0% [0.0%, 1.3%] |
| coloc.susie + coloc.abf fallback | 12.0% [8.8%, 16.2%] | 0.0% [0.0%, 1.3%] |
| coloc.signals (mask) | 4.3% [2.5%, 7.3%] | 0.0% [0.0%, 1.3%] |
| coloc.signals (cond) | 2.7% [1.4%, 5.2%] | 0.0% [0.0%, 1.3%] |

## 2. Pooled over scenarios (test windows)

TPR pools the shared scenarios (single_shared, single_double_shared, double_shared2, double_shared1); FPR pools single_distinct and double_distinct. AUC does not depend on any cut-off.

| method | lead_p | n_pos | TPR | n_neg | FPR | AUC | no_call | branch_block |
|---|---|---|---|---|---|---|---|---|
| block-coloc (conditional lead selection, proposed) | 1e-03 | 1200 | 51.7% [48.8%, 54.5%] | 600 | 7.3% [5.5%, 9.7%] | 0.786 [0.766, 0.807] | 0.0% | 31.4% |
| block-coloc (no fallback) | 1e-03 | 1200 | 50.0% [47.2%, 52.8%] | 600 | 7.0% [5.2%, 9.3%] | 0.770 [0.750, 0.792] | 16.1% | 31.4% |
| block-coloc (marginal lead selection) | 1e-04 | 1200 | 49.2% [46.4%, 52.1%] | 600 | 8.2% [6.2%, 10.6%] | 0.765 [0.744, 0.785] | 15.6% | 34.1% |
| block-coloc (original handoff) | 1e-04 | 1200 | 50.9% [48.1%, 53.7%] | 600 | 8.7% [6.7%, 11.2%] | 0.762 [0.741, 0.784] | 10.2% | 41.2% |
| coloc.abf, whole region | - | 1200 | 46.2% [43.4%, 49.0%] | 600 | 7.0% [5.2%, 9.3%] | 0.755 [0.732, 0.775] | 0.0% | 0.0% |
| coloc.susie | - | 1200 | 30.2% [27.6%, 32.8%] | 600 | 4.2% [2.8%, 6.1%] | 0.573 [0.550, 0.595] | 60.7% | 0.0% |
| coloc.susie + coloc.abf fallback | - | 1200 | 48.8% [46.0%, 51.7%] | 600 | 7.5% [5.7%, 9.9%] | 0.766 [0.745, 0.787] | 0.0% | 0.0% |
| coloc.signals (mask) | - | 1200 | 51.3% [48.5%, 54.2%] | 600 | 7.5% [5.7%, 9.9%] | 0.763 [0.742, 0.785] | 0.0% | 0.0% |
| coloc.signals (cond) | - | 1200 | 50.9% [48.1%, 53.7%] | 600 | 7.2% [5.4%, 9.5%] | 0.750 [0.727, 0.772] | 0.0% | 0.0% |

## 3. Cut-offs chosen on calibration windows

lead_p = lead threshold (GWAS / eQTL when they differ). '*' marks the lead threshold kept for each block-coloc version.

| method | lead_p | cutoff | calib_FPR | calib_TPR |
|---|---|---|---|---|
| block-coloc (conditional lead selection, proposed) | 1e-03 * | 0.0702 | 5.0% | 55.1% |
| block-coloc (conditional lead selection, proposed) | 1e-04 | 0.0697 | 5.0% | 53.8% |
| block-coloc (conditional lead selection, proposed) | 1e-05 | 0.0697 | 5.0% | 52.0% |
| block-coloc (conditional lead selection, proposed) [earlier bulk thresholds] | 5e-08 / 1e-05 | 0.0697 | 5.0% | 51.7% |
| block-coloc (no fallback) | 1e-03 * | 0.0676 | 5.0% | 53.6% |
| block-coloc (no fallback) | 1e-04 | 0.0697 | 5.0% | 53.2% |
| block-coloc (no fallback) | 1e-05 | 0.0697 | 5.0% | 51.9% |
| block-coloc (no fallback) [earlier bulk thresholds] | 5e-08 / 1e-05 | 0.0697 | 5.0% | 51.7% |
| block-coloc (marginal lead selection) | 1e-03 | 0.1354 | 5.0% | 50.0% |
| block-coloc (marginal lead selection) | 1e-04 * | 0.0721 | 5.0% | 53.6% |
| block-coloc (marginal lead selection) | 1e-05 | 0.0697 | 5.0% | 52.5% |
| block-coloc (marginal lead selection) [earlier bulk thresholds] | 5e-08 / 1e-05 | 0.0702 | 5.0% | 52.0% |
| block-coloc (original handoff) | 1e-03 | 0.1442 | 5.0% | 51.1% |
| block-coloc (original handoff) | 1e-04 * | 0.0749 | 5.0% | 53.8% |
| block-coloc (original handoff) | 1e-05 | 0.0702 | 5.0% | 52.8% |
| block-coloc (original handoff) [earlier bulk thresholds] | 5e-08 / 1e-05 | 0.0704 | 5.0% | 52.8% |
| block-coloc (marginal lead selection) + whole-region fallback | 1e-03 | 0.1373 | 5.0% | 51.8% |
| block-coloc (marginal lead selection) + whole-region fallback | 1e-04 * | 0.0749 | 5.0% | 54.2% |
| block-coloc (marginal lead selection) + whole-region fallback | 1e-05 | 0.0702 | 5.0% | 52.8% |
| block-coloc (marginal lead selection) + whole-region fallback [earlier bulk thresholds] | 5e-08 / 1e-05 | 0.0702 | 5.0% | 52.2% |
| block-coloc (original handoff) + whole-region fallback | 1e-03 | 0.1604 | 5.0% | 51.3% |
| block-coloc (original handoff) + whole-region fallback | 1e-04 * | 0.0755 | 5.0% | 54.8% |
| block-coloc (original handoff) + whole-region fallback | 1e-05 | 0.0704 | 5.0% | 53.9% |
| block-coloc (original handoff) + whole-region fallback [earlier bulk thresholds] | 5e-08 / 1e-05 | 0.0706 | 5.0% | 53.3% |
| coloc.abf, whole region | - | 0.0697 | 5.0% | 49.6% |
| coloc.susie | - | 0.0076 | 5.0% | 29.9% |
| coloc.susie + coloc.abf fallback | - | 0.0697 | 5.0% | 51.8% |
| coloc.signals (mask) | - | 0.1021 | 5.0% | 53.2% |
| coloc.signals (cond) | - | 0.1361 | 5.0% | 53.8% |

## 4. Paired differences on the same test loci (McNemar)

Difference = method minus reference, in percentage points, with 95% CI and exact McNemar p-value.

| method | vs | TPR_difference | FPR_difference |
|---|---|---|---|
| block-coloc (conditional lead selection, proposed) | coloc.susie + coloc.abf fallback | +2.8 pts [+1.7, +4.0], p=1.2e-06 | -0.2 pts [-1.0, +0.7], p=1 |
| block-coloc (no fallback) | coloc.susie + coloc.abf fallback | +1.2 pts [-0.3, +2.6], p=0.14 | -0.5 pts [-1.5, +0.5], p=0.51 |
| block-coloc (marginal lead selection) | coloc.susie + coloc.abf fallback | +0.4 pts [-0.7, +1.6], p=0.58 | +0.7 pts [-1.0, +2.3], p=0.56 |
| block-coloc (original handoff) | coloc.susie + coloc.abf fallback | +2.1 pts [+0.6, +3.5], p=0.0073 | +1.2 pts [-0.8, +3.1], p=0.31 |
| coloc.abf, whole region | coloc.susie + coloc.abf fallback | -2.7 pts [-3.8, -1.5], p=3.3e-06 | -0.5 pts [-1.2, +0.2], p=0.37 |
| coloc.susie | coloc.susie + coloc.abf fallback | -18.7 pts [-21.1, -16.2], p=<2e-16 | -3.3 pts [-5.9, -0.8], p=0.015 |
| coloc.signals (mask) | coloc.susie + coloc.abf fallback | +2.5 pts [+0.8, +4.2], p=0.0042 | +0.0 pts [-1.8, +1.8], p=1 |
| coloc.signals (cond) | coloc.susie + coloc.abf fallback | +2.1 pts [+0.2, +3.9], p=0.034 | -0.3 pts [-2.4, +1.7], p=0.87 |
| block-coloc (conditional lead selection, proposed) | coloc.abf, whole region | +5.5 pts [+4.2, +6.8], p=<2e-16 | +0.3 pts [-0.1, +0.8], p=0.5 |
| block-coloc (no fallback) | coloc.abf, whole region | +3.8 pts [+2.2, +5.5], p=4.7e-06 | +0.0 pts [-0.7, +0.7], p=1 |
| block-coloc (marginal lead selection) | coloc.abf, whole region | +3.1 pts [+1.7, +4.5], p=2.2e-05 | +1.2 pts [-0.4, +2.7], p=0.21 |
| block-coloc (original handoff) | coloc.abf, whole region | +4.8 pts [+3.1, +6.4], p=6.9e-09 | +1.7 pts [-0.1, +3.5], p=0.099 |
| coloc.susie | coloc.abf, whole region | -16.0 pts [-18.7, -13.3], p=<2e-16 | -2.8 pts [-5.5, -0.2], p=0.046 |
| coloc.susie + coloc.abf fallback | coloc.abf, whole region | +2.7 pts [+1.5, +3.8], p=3.3e-06 | +0.5 pts [-0.2, +1.2], p=0.37 |
| coloc.signals (mask) | coloc.abf, whole region | +5.2 pts [+3.4, +6.9], p=6.4e-09 | +0.5 pts [-1.4, +2.4], p=0.73 |
| coloc.signals (cond) | coloc.abf, whole region | +4.8 pts [+2.7, +6.8], p=8.3e-06 | +0.2 pts [-2.0, +2.4], p=1 |

### 4b. Adopted block-coloc rule and the no-fallback sensitivity analysis

The proposed method uses whole-region PP4 when no valid block is available; 'block-coloc (no fallback)' keeps those loci as no result (section 8).

| method | vs | TPR_difference | FPR_difference |
|---|---|---|---|
| block-coloc (conditional lead selection, proposed) | coloc.abf, whole region | +5.5 pts [+4.2, +6.8], p=<2e-16 | +0.3 pts [-0.1, +0.8], p=0.5 |
| block-coloc (conditional lead selection, proposed) | coloc.susie + coloc.abf fallback | +2.8 pts [+1.7, +4.0], p=1.2e-06 | -0.2 pts [-1.0, +0.7], p=1 |
| block-coloc (conditional lead selection, proposed) | coloc.signals (cond) | +0.8 pts [-0.9, +2.4], p=0.42 | +0.2 pts [-2.0, +2.3], p=1 |
| block-coloc (conditional lead selection, proposed) | coloc.signals (mask) | +0.3 pts [-1.1, +1.8], p=0.74 | -0.2 pts [-2.0, +1.7], p=1 |
| block-coloc (no fallback) | coloc.signals (cond) | -0.9 pts [-2.7, +0.9], p=0.37 | -0.2 pts [-2.3, +2.0], p=1 |
| block-coloc (no fallback) | coloc.signals (mask) | -1.3 pts [-3.0, +0.4], p=0.15 | -0.5 pts [-2.4, +1.4], p=0.73 |

## 5. Branch choice: split into blocks or not (test windows)

should_go: block when either trait has 2 causal SNPs, whole otherwise. leads_match_truth: the method's lead counts
(GWAS, eQTL) equal the true numbers of causal SNPs.

| method | scenario | true_signals | should_go | n | branch_block | correct_branch | leads_match_truth | median_leads |
|---|---|---|---|---|---|---|---|---|
| block-coloc (conditional lead selection, proposed) | single_shared | 1+1 | whole | 300 | 11.0% | 89.0% | 43.3% | 1+1 |
| block-coloc (conditional lead selection, proposed) | single_double_shared | 1+2 | block | 300 | 26.7% | 26.7% | 19.7% | 1+1 |
| block-coloc (conditional lead selection, proposed) | double_shared2 | 2+2 | block | 300 | 47.0% | 47.0% | 9.3% | 1+1 |
| block-coloc (conditional lead selection, proposed) | double_shared1 | 2+2 | block | 300 | 47.3% | 47.3% | 10.0% | 1+1 |
| block-coloc (conditional lead selection, proposed) | single_distinct | 1+1 | whole | 300 | 11.7% | 88.3% | 39.3% | 1+1 |
| block-coloc (conditional lead selection, proposed) | double_distinct | 2+2 | block | 300 | 45.0% | 45.0% | 9.3% | 1+1 |
| block-coloc (conditional lead selection, proposed) | gwas_only | 1+0 | whole | 300 | 1.0% | 99.0% | 58.0% | 1+0 |
| block-coloc (conditional lead selection, proposed) | null_both | 0+0 | whole | 300 | 0.3% | 99.7% | 74.0% | 0+0 |
| block-coloc (no fallback) | single_shared | 1+1 | whole | 300 | 11.0% | 89.0% | 43.3% | 1+1 |
| block-coloc (no fallback) | single_double_shared | 1+2 | block | 300 | 26.7% | 26.7% | 19.7% | 1+1 |
| block-coloc (no fallback) | double_shared2 | 2+2 | block | 300 | 47.0% | 47.0% | 9.3% | 1+1 |
| block-coloc (no fallback) | double_shared1 | 2+2 | block | 300 | 47.3% | 47.3% | 10.0% | 1+1 |
| block-coloc (no fallback) | single_distinct | 1+1 | whole | 300 | 11.7% | 88.3% | 39.3% | 1+1 |
| block-coloc (no fallback) | double_distinct | 2+2 | block | 300 | 45.0% | 45.0% | 9.3% | 1+1 |
| block-coloc (no fallback) | gwas_only | 1+0 | whole | 300 | 1.0% | 99.0% | 58.0% | 1+0 |
| block-coloc (no fallback) | null_both | 0+0 | whole | 300 | 0.3% | 99.7% | 74.0% | 0+0 |
| block-coloc (marginal lead selection) | single_shared | 1+1 | whole | 300 | 26.7% | 73.3% | 16.0% | 1+1 |
| block-coloc (marginal lead selection) | single_double_shared | 1+2 | block | 300 | 30.3% | 30.3% | 8.7% | 1+1 |
| block-coloc (marginal lead selection) | double_shared2 | 2+2 | block | 300 | 40.7% | 40.7% | 4.0% | 1+1 |
| block-coloc (marginal lead selection) | double_shared1 | 2+2 | block | 300 | 44.3% | 44.3% | 3.7% | 1+1 |
| block-coloc (marginal lead selection) | single_distinct | 1+1 | whole | 300 | 22.0% | 78.0% | 14.0% | 1+1 |
| block-coloc (marginal lead selection) | double_distinct | 2+2 | block | 300 | 40.7% | 40.7% | 4.3% | 1+1 |
| block-coloc (marginal lead selection) | gwas_only | 1+0 | whole | 300 | 0.0% | 100.0% | 37.7% | 1+0 |
| block-coloc (marginal lead selection) | null_both | 0+0 | whole | 300 | 0.0% | 100.0% | 96.0% | 0+0 |
| block-coloc (original handoff) | single_shared | 1+1 | whole | 300 | 56.7% | 43.3% | 36.7% | 1+1 |
| block-coloc (original handoff) | single_double_shared | 1+2 | block | 300 | 43.7% | 43.7% | 17.0% | 1+1 |
| block-coloc (original handoff) | double_shared2 | 2+2 | block | 300 | 55.0% | 55.0% | 7.7% | 2+1 |
| block-coloc (original handoff) | double_shared1 | 2+2 | block | 300 | 43.0% | 43.0% | 7.7% | 1+1 |
| block-coloc (original handoff) | single_distinct | 1+1 | whole | 300 | 23.0% | 77.0% | 35.3% | 1+1 |
| block-coloc (original handoff) | double_distinct | 2+2 | block | 300 | 25.7% | 25.7% | 7.3% | 1+1 |
| block-coloc (original handoff) | gwas_only | 1+0 | whole | 300 | 31.3% | 68.7% | 0.0% | 1+2 |
| block-coloc (original handoff) | null_both | 0+0 | whole | 300 | 30.7% | 69.3% | 0.0% | 2+2 |

## 6. Detection when the shared signal is not the stronger one (test windows)

single_double_shared: the eQTL has two causal variants, one shared with the GWAS and one eQTL-only. 'stronger' compares their observed eQTL |z|.

| method | n_shared_stronger | TPR_shared_stronger | n_other_stronger | TPR_other_stronger |
|---|---|---|---|---|
| block-coloc (conditional lead selection, proposed) | 138 | 67.4% [59.2%, 74.6%] | 162 | 21.0% [15.4%, 27.9%] |
| block-coloc (no fallback) | 138 | 65.9% [57.7%, 73.3%] | 162 | 19.8% [14.4%, 26.6%] |
| block-coloc (marginal lead selection) | 138 | 65.9% [57.7%, 73.3%] | 162 | 19.1% [13.8%, 25.9%] |
| block-coloc (original handoff) | 138 | 68.1% [59.9%, 75.3%] | 162 | 22.2% [16.5%, 29.2%] |
| coloc.abf, whole region | 138 | 67.4% [59.2%, 74.6%] | 162 | 12.3% [8.1%, 18.3%] |
| coloc.susie | 138 | 35.5% [28.0%, 43.8%] | 162 | 9.9% [6.2%, 15.4%] |
| coloc.susie + coloc.abf fallback | 138 | 67.4% [59.2%, 74.6%] | 162 | 16.0% [11.2%, 22.5%] |
| coloc.signals (mask) | 138 | 66.7% [58.4%, 74.0%] | 162 | 22.2% [16.5%, 29.2%] |
| coloc.signals (cond) | 138 | 62.3% [54.0%, 70.0%] | 162 | 27.8% [21.5%, 35.1%] |

double_shared1: each trait has one shared and one trait-specific causal variant. Columns: in how many traits the shared variant has the larger observed |z|.

| method | shared_stronger_in_both | in_one_trait | in_neither |
|---|---|---|---|
| block-coloc (conditional lead selection, proposed) | 86.9% [78.1%, 92.5%] (n=84) | 24.8% [18.5%, 32.4%] (n=145) | 11.3% [5.8%, 20.7%] (n=71) |
| block-coloc (no fallback) | 84.5% [75.3%, 90.7%] (n=84) | 22.8% [16.7%, 30.2%] (n=145) | 11.3% [5.8%, 20.7%] (n=71) |
| block-coloc (marginal lead selection) | 86.9% [78.1%, 92.5%] (n=84) | 19.3% [13.7%, 26.5%] (n=145) | 11.3% [5.8%, 20.7%] (n=71) |
| block-coloc (original handoff) | 88.1% [79.5%, 93.4%] (n=84) | 24.8% [18.5%, 32.4%] (n=145) | 9.9% [4.9%, 19.0%] (n=71) |
| coloc.abf, whole region | 88.1% [79.5%, 93.4%] (n=84) | 13.1% [8.6%, 19.6%] (n=145) | 2.8% [0.8%, 9.7%] (n=71) |
| coloc.susie | 57.1% [46.5%, 67.2%] (n=84) | 20.0% [14.3%, 27.2%] (n=145) | 11.3% [5.8%, 20.7%] (n=71) |
| coloc.susie + coloc.abf fallback | 86.9% [78.1%, 92.5%] (n=84) | 20.7% [14.9%, 28.0%] (n=145) | 7.0% [3.0%, 15.4%] (n=71) |
| coloc.signals (mask) | 86.9% [78.1%, 92.5%] (n=84) | 33.8% [26.6%, 41.8%] (n=145) | 11.3% [5.8%, 20.7%] (n=71) |
| coloc.signals (cond) | 86.9% [78.1%, 92.5%] (n=84) | 31.7% [24.7%, 39.7%] (n=145) | 12.7% [6.8%, 22.4%] (n=71) |

## 7. Localisation among detected test positives

Sub-region = best block (block branch) or the two credible sets of the best SuSiE signal pair.
Whole-locus calls and coloc.signals calls have no sub-region and are not counted here.

| method | detected | with_subregion | contains_shared | median_size_snps |
|---|---|---|---|---|
| block-coloc (conditional lead selection, proposed) | 620 | 262 | 99.6% [97.9%, 99.9%] | 46 |
| block-coloc (no fallback) | 600 | 263 | 99.6% [97.9%, 99.9%] | 46 |
| block-coloc (marginal lead selection) | 591 | 305 | 95.4% [92.4%, 97.2%] | 46 |
| block-coloc (original handoff) | 611 | 465 | 94.8% [92.4%, 96.5%] | 46 |
| coloc.abf, whole region | 554 |   0 | NA | NA |
| coloc.susie | 362 | 362 | 98.1% [96.1%, 99.1%] | 20 |
| coloc.susie + coloc.abf fallback | 586 | 338 | 100.0% [98.9%, 100.0%] | 18 |
| coloc.signals (mask) | 616 |   0 | NA | NA |
| coloc.signals (cond) | 611 |   0 | NA | NA |

## 8. Sensitivity: fallback policy and the earlier bulk lead thresholds (test windows)

The proposed conditional method includes whole-region fallback; 'block-coloc (no fallback)' omits it.
'+ whole-region fallback' = when the block branch builds no valid block, use the whole-region PP4 instead of no result.
'[earlier bulk thresholds]' = lead thresholds GWAS 5e-8, eQTL 1e-5 (August bulk setting), each with its own calibrated cut-off.

| method | lead_p | cutoff | n_pos | TPR | n_neg | FPR | AUC | no_call | branch_block |
|---|---|---|---|---|---|---|---|---|---|
| block-coloc (conditional lead selection, proposed) | 1e-03 | 0.0702 | 1200 | 51.7% [48.8%, 54.5%] | 600 | 7.3% [5.5%, 9.7%] | 0.786 | 0.0% | 31.4% |
| block-coloc (no fallback) | 1e-03 | 0.0676 | 1200 | 50.0% [47.2%, 52.8%] | 600 | 7.0% [5.2%, 9.3%] | 0.770 | 16.1% | 31.4% |
| block-coloc (marginal lead selection) | 1e-04 | 0.0721 | 1200 | 49.2% [46.4%, 52.1%] | 600 | 8.2% [6.2%, 10.6%] | 0.765 | 15.6% | 34.1% |
| block-coloc (original handoff) | 1e-04 | 0.0749 | 1200 | 50.9% [48.1%, 53.7%] | 600 | 8.7% [6.7%, 11.2%] | 0.762 | 10.2% | 41.2% |
| block-coloc (marginal lead selection) + whole-region fallback | 1e-04 | 0.0749 | 1200 | 50.2% [47.3%, 53.0%] | 600 | 8.7% [6.7%, 11.2%] | 0.773 | 0.0% | 34.1% |
| block-coloc (original handoff) + whole-region fallback | 1e-04 | 0.0755 | 1200 | 52.2% [49.3%, 55.0%] | 600 | 9.5% [7.4%, 12.1%] | 0.777 | 0.0% | 41.2% |
| block-coloc (original handoff) [earlier bulk thresholds] | 5e-08 / 1e-05 | 0.0704 | 1200 | 50.1% [47.3%, 52.9%] | 600 | 7.7% [5.8%, 10.1%] | 0.764 | 5.3% | 31.2% |
| block-coloc (marginal lead selection) [earlier bulk thresholds] | 5e-08 / 1e-05 | 0.0702 | 1200 | 48.5% [45.7%, 51.3%] | 600 | 7.5% [5.7%, 9.9%] | 0.772 | 7.1% | 15.8% |
| block-coloc (no fallback) [earlier bulk thresholds] | 5e-08 / 1e-05 | 0.0697 | 1200 | 48.0% [45.2%, 50.8%] | 600 | 7.0% [5.2%, 9.3%] | 0.771 | 3.2% | 7.5% |
| block-coloc (conditional lead selection, proposed) [earlier bulk thresholds] | 5e-08 / 1e-05 | 0.0697 | 1200 | 48.2% [45.3%, 51.0%] | 600 | 7.0% [5.2%, 9.3%] | 0.770 | 0.0% | 7.5% |

## 9. Inputs (test windows)

Expected causal |z| (normal effects: E|z| = sqrt(2 / pi * (1 + N q / (1 - q))), q = h2 / K): GWAS 6.2 (K = 1) / 4.4 (K = 2); eQTL 5.8 (K = 1) / 4.1 (K = 2). Table: medians over loci and causal SNPs.

| scenario | n | causal_z_gwas | causal_z_eqtl | pair_abs_r_median | gwas_no_cs | eqtl_no_cs |
|---|---|---|---|---|---|---|
| single_shared | 300 | 5.6 | 5.8 | NA | 40.0% | 38.0% |
| single_double_shared | 300 | 5.3 | 3.4 | 0.03 | 41.0% | 38.0% |
| double_shared2 | 300 | 4.0 | 3.2 | 0.03 | 29.0% | 37.3% |
| double_shared1 | 300 | 4.0 | 3.6 | 0.03 | 29.3% | 35.0% |
| single_distinct | 300 | 5.4 | 4.6 | 0.03 | 43.7% | 47.3% |
| double_distinct | 300 | 3.9 | 3.4 | 0.03 | 30.7% | 36.3% |
| gwas_only | 300 | 5.3 | NA | NA | 39.0% | 100.0% |
| null_both | 300 | NA | NA | NA | 99.3% | 99.7% |

## Figures

- Main figure: `matched_fpr_by_scenario.png`
- PP4 distribution by scenario and method: `pp4_distribution.png`
- `roc_test.png`

## Reproduce

```
Rscript scripts/bulk_h2ref_v1/run.R --stage main --effect normal --h2_eqtl 0.05 --reps 50 --cores 4 && Rscript scripts/bulk_h2ref_v1/evaluate.R --stage main --out outputs/bulk_h2ref_v1_normal_g0.003_e0.05
```

## Files changed (git diff --stat)

(Codex pastes the output of `git diff --stat HEAD~1` here after committing)
