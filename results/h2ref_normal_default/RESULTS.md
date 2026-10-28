# Bulk block-coloc: main comparison (stage 1)

Generated 2026-10-29 03:29. Setting: GWAS N 20,000, local h2 0.3%; eQTL N 1,000, local h2 1.5%; the effect of each causal SNP is drawn from N(0, h2 / K), independently for the two traits, without rescaling (Var(g) = h2 in expectation only; sensitivity analysis); causal |r| < 0.2.
LD given to the methods: 1000 Genomes EUR reference LD (489 people) for both traits.
Calibration windows: 31263287, 37199593, 39694225, 42313324, 44324730, 46518082. Test windows: 21059630, 24009444, 26904965, 29033189, 33343871, 34894283.
Loci: 3600 calibration, 3600 test (per method). Target FPR 5% on calibration windows (negatives: single_distinct and double_distinct).
Extra negatives gwas_only and null_both are reported for the false-positive rate only (not calibrated, not pooled).

## Three-line conclusion

1. **检出率**：效应改为从 N(0, h²/K) 抽取（敏感性分析）后，在默认强度（GWAS h² 0.3%、eQTL h² 1.5%）、误报率在本设计的校准窗口上匹配到 5%（测试窗口实际 5.8%–6.3%）时，所有基于 coloc.abf 的方法检出率都降到 38%–40%（固定效应下为 73%–77%）。block-coloc 检出 39.7%，比整区域 coloc.abf（37.8%）高 1.9 个百分点（配对 95% CI +1.3 至 +2.5，p = 2e-11，6 个测试窗口全部为正），比 coloc.susie + fallback 高 1.5 个百分点，与 coloc.signals (cond)（+0.5，p = 0.36）和 (mask)（+0.2，p = 0.59）无差别；AUC 0.708，整区域 coloc.abf 0.698。
2. **为什么整体变低、增益变小**：正态效应下因果变异强弱不一，eQTL 因果变异的 |z| 有 36%（单因果）到 50%（双因果）低于 2；两因果的 eQTL 只有 11.4% 找到两个条件独立的 lead（固定效应下 15%–17%），block-coloc 只在 18.3% 的阳性位点上切块，其余位点与整区域 coloc.abf 相同。lead 阈值按本文方法的校准检出率选为 1e-3（41.8%，1e-4 为 41.0%）。
3. **额外阴性**：gwas_only 误报率 block-coloc 6.8%，整区域 coloc.abf 7.0%，coloc.signals 3.7%–4.0%；null_both 所有方法为 0%。

## 1. Main table (test windows, frozen cut-offs)

n_positive = loci called colocalised (PP4 above the method's cut-off). For shared scenarios rate is the detection rate;
for not-shared scenarios it is the false-positive rate. Denominators include every locus; no result counts as not called.
'not shared (extra)' rows (gwas_only, null_both) use the same frozen cut-offs but were not used to set them.

| scenario | truth | method | lead_p | n_total | n_positive | rate | 95% CI | no_call_rate | branch_block |
|---|---|---|---|---|---|---|---|---|---|
| single_shared | shared | block-coloc (conditional lead selection, proposed) | 1e-03 | 600 | 290 | 48.3% | [44.4%, 52.3%] | 0.0% | 8.2% |
| single_shared | shared | block-coloc (no fallback) | 1e-03 | 600 | 291 | 48.5% | [44.5%, 52.5%] | 1.3% | 8.2% |
| single_shared | shared | block-coloc (marginal lead selection) | 1e-04 | 600 | 289 | 48.2% | [44.2%, 52.2%] | 0.5% | 10.8% |
| single_shared | shared | block-coloc (original handoff) | 1e-04 | 600 | 277 | 46.2% | [42.2%, 50.2%] | 9.3% | 50.5% |
| single_shared | shared | coloc.abf, whole region | - | 600 | 290 | 48.3% | [44.4%, 52.3%] | 0.0% | 0.0% |
| single_shared | shared | coloc.susie | - | 600 | 106 | 17.7% | [14.8%, 20.9%] | 82.3% | 0.0% |
| single_shared | shared | coloc.susie + coloc.abf fallback | - | 600 | 289 | 48.2% | [44.2%, 52.2%] | 0.0% | 0.0% |
| single_shared | shared | coloc.signals (mask) | - | 600 | 279 | 46.5% | [42.5%, 50.5%] | 0.0% | 0.0% |
| single_shared | shared | coloc.signals (cond) | - | 600 | 272 | 45.3% | [41.4%, 49.3%] | 0.0% | 0.0% |
| single_double_shared | shared | block-coloc (conditional lead selection, proposed) | 1e-03 | 600 | 209 | 34.8% | [31.1%, 38.7%] | 0.0% | 11.7% |
| single_double_shared | shared | block-coloc (no fallback) | 1e-03 | 600 | 201 | 33.5% | [29.8%, 37.4%] | 4.3% | 11.7% |
| single_double_shared | shared | block-coloc (marginal lead selection) | 1e-04 | 600 | 203 | 33.8% | [30.2%, 37.7%] | 4.7% | 10.8% |
| single_double_shared | shared | block-coloc (original handoff) | 1e-04 | 600 | 203 | 33.8% | [30.2%, 37.7%] | 10.8% | 37.7% |
| single_double_shared | shared | coloc.abf, whole region | - | 600 | 204 | 34.0% | [30.3%, 37.9%] | 0.0% | 0.0% |
| single_double_shared | shared | coloc.susie | - | 600 |  56 | 9.3% | [7.3%, 11.9%] | 86.3% | 0.0% |
| single_double_shared | shared | coloc.susie + coloc.abf fallback | - | 600 | 204 | 34.0% | [30.3%, 37.9%] | 0.0% | 0.0% |
| single_double_shared | shared | coloc.signals (mask) | - | 600 | 202 | 33.7% | [30.0%, 37.5%] | 0.0% | 0.0% |
| single_double_shared | shared | coloc.signals (cond) | - | 600 | 193 | 32.2% | [28.6%, 36.0%] | 0.0% | 0.0% |
| double_shared2 | shared | block-coloc (conditional lead selection, proposed) | 1e-03 | 600 | 265 | 44.2% | [40.2%, 48.2%] | 0.0% | 26.0% |
| double_shared2 | shared | block-coloc (no fallback) | 1e-03 | 600 | 251 | 41.8% | [38.0%, 45.8%] | 8.5% | 26.0% |
| double_shared2 | shared | block-coloc (marginal lead selection) | 1e-04 | 600 | 252 | 42.0% | [38.1%, 46.0%] | 3.0% | 16.3% |
| double_shared2 | shared | block-coloc (original handoff) | 1e-04 | 600 | 266 | 44.3% | [40.4%, 48.3%] | 8.5% | 46.7% |
| double_shared2 | shared | coloc.abf, whole region | - | 600 | 248 | 41.3% | [37.5%, 45.3%] | 0.0% | 0.0% |
| double_shared2 | shared | coloc.susie | - | 600 |  82 | 13.7% | [11.1%, 16.6%] | 83.0% | 0.0% |
| double_shared2 | shared | coloc.susie + coloc.abf fallback | - | 600 | 250 | 41.7% | [37.8%, 45.7%] | 0.0% | 0.0% |
| double_shared2 | shared | coloc.signals (mask) | - | 600 | 269 | 44.8% | [40.9%, 48.8%] | 0.0% | 0.0% |
| double_shared2 | shared | coloc.signals (cond) | - | 600 | 278 | 46.3% | [42.4%, 50.3%] | 0.0% | 0.0% |
| double_shared1 | shared | block-coloc (conditional lead selection, proposed) | 1e-03 | 600 | 189 | 31.5% | [27.9%, 35.3%] | 0.0% | 29.2% |
| double_shared1 | shared | block-coloc (no fallback) | 1e-03 | 600 | 173 | 28.8% | [25.4%, 32.6%] | 15.2% | 29.2% |
| double_shared1 | shared | block-coloc (marginal lead selection) | 1e-04 | 600 | 176 | 29.3% | [25.8%, 33.1%] | 9.0% | 19.3% |
| double_shared1 | shared | block-coloc (original handoff) | 1e-04 | 600 | 188 | 31.3% | [27.8%, 35.2%] | 11.0% | 39.0% |
| double_shared1 | shared | coloc.abf, whole region | - | 600 | 166 | 27.7% | [24.2%, 31.4%] | 0.0% | 0.0% |
| double_shared1 | shared | coloc.susie | - | 600 |  75 | 12.5% | [10.1%, 15.4%] | 81.0% | 0.0% |
| double_shared1 | shared | coloc.susie + coloc.abf fallback | - | 600 | 175 | 29.2% | [25.7%, 32.9%] | 0.0% | 0.0% |
| double_shared1 | shared | coloc.signals (mask) | - | 600 | 197 | 32.8% | [29.2%, 36.7%] | 0.0% | 0.0% |
| double_shared1 | shared | coloc.signals (cond) | - | 600 | 199 | 33.2% | [29.5%, 37.0%] | 0.0% | 0.0% |
| single_distinct | not shared | block-coloc (conditional lead selection, proposed) | 1e-03 | 600 |  39 | 6.5% | [4.8%, 8.8%] | 0.0% | 8.8% |
| single_distinct | not shared | block-coloc (no fallback) | 1e-03 | 600 |  37 | 6.2% | [4.5%, 8.4%] | 8.5% | 8.8% |
| single_distinct | not shared | block-coloc (marginal lead selection) | 1e-04 | 600 |  38 | 6.3% | [4.6%, 8.6%] | 11.7% | 12.2% |
| single_distinct | not shared | block-coloc (original handoff) | 1e-04 | 600 |  44 | 7.3% | [5.5%, 9.7%] | 10.8% | 19.3% |
| single_distinct | not shared | coloc.abf, whole region | - | 600 |  38 | 6.3% | [4.6%, 8.6%] | 0.0% | 0.0% |
| single_distinct | not shared | coloc.susie | - | 600 |  19 | 3.2% | [2.0%, 4.9%] | 83.0% | 0.0% |
| single_distinct | not shared | coloc.susie + coloc.abf fallback | - | 600 |  37 | 6.2% | [4.5%, 8.4%] | 0.0% | 0.0% |
| single_distinct | not shared | coloc.signals (mask) | - | 600 |  34 | 5.7% | [4.1%, 7.8%] | 0.0% | 0.0% |
| single_distinct | not shared | coloc.signals (cond) | - | 600 |  30 | 5.0% | [3.5%, 7.0%] | 0.0% | 0.0% |
| double_distinct | not shared | block-coloc (conditional lead selection, proposed) | 1e-03 | 600 |  36 | 6.0% | [4.4%, 8.2%] | 0.0% | 26.0% |
| double_distinct | not shared | block-coloc (no fallback) | 1e-03 | 600 |  29 | 4.8% | [3.4%, 6.9%] | 24.7% | 26.0% |
| double_distinct | not shared | block-coloc (marginal lead selection) | 1e-04 | 600 |  32 | 5.3% | [3.8%, 7.4%] | 15.5% | 16.5% |
| double_distinct | not shared | block-coloc (original handoff) | 1e-04 | 600 |  37 | 6.2% | [4.5%, 8.4%] | 15.5% | 21.5% |
| double_distinct | not shared | coloc.abf, whole region | - | 600 |  35 | 5.8% | [4.2%, 8.0%] | 0.0% | 0.0% |
| double_distinct | not shared | coloc.susie | - | 600 |  34 | 5.7% | [4.1%, 7.8%] | 84.0% | 0.0% |
| double_distinct | not shared | coloc.susie + coloc.abf fallback | - | 600 |  33 | 5.5% | [3.9%, 7.6%] | 0.0% | 0.0% |
| double_distinct | not shared | coloc.signals (mask) | - | 600 |  38 | 6.3% | [4.6%, 8.6%] | 0.0% | 0.0% |
| double_distinct | not shared | coloc.signals (cond) | - | 600 |  42 | 7.0% | [5.2%, 9.3%] | 0.0% | 0.0% |
| gwas_only | not shared (extra) | block-coloc (conditional lead selection, proposed) | 1e-03 | 600 |  41 | 6.8% | [5.1%, 9.1%] | 0.0% | 2.0% |
| gwas_only | not shared (extra) | block-coloc (no fallback) | 1e-03 | 600 |  44 | 7.3% | [5.5%, 9.7%] | 2.0% | 2.0% |
| gwas_only | not shared (extra) | block-coloc (marginal lead selection) | 1e-04 | 600 |  41 | 6.8% | [5.1%, 9.1%] | 0.5% | 0.7% |
| gwas_only | not shared (extra) | block-coloc (original handoff) | 1e-04 | 600 |  44 | 7.3% | [5.5%, 9.7%] | 19.3% | 28.5% |
| gwas_only | not shared (extra) | coloc.abf, whole region | - | 600 |  42 | 7.0% | [5.2%, 9.3%] | 0.0% | 0.0% |
| gwas_only | not shared (extra) | coloc.susie | - | 600 |   1 | 0.2% | [0.0%, 0.9%] | 99.8% | 0.0% |
| gwas_only | not shared (extra) | coloc.susie + coloc.abf fallback | - | 600 |  43 | 7.2% | [5.4%, 9.5%] | 0.0% | 0.0% |
| gwas_only | not shared (extra) | coloc.signals (mask) | - | 600 |  24 | 4.0% | [2.7%, 5.9%] | 0.0% | 0.0% |
| gwas_only | not shared (extra) | coloc.signals (cond) | - | 600 |  22 | 3.7% | [2.4%, 5.5%] | 0.0% | 0.0% |
| null_both | not shared (extra) | block-coloc (conditional lead selection, proposed) | 1e-03 | 600 |   0 | 0.0% | [0.0%, 0.6%] | 0.0% | 0.3% |
| null_both | not shared (extra) | block-coloc (no fallback) | 1e-03 | 600 |   0 | 0.0% | [0.0%, 0.6%] | 0.3% | 0.3% |
| null_both | not shared (extra) | block-coloc (marginal lead selection) | 1e-04 | 600 |   0 | 0.0% | [0.0%, 0.6%] | 0.0% | 0.0% |
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
| block-coloc (conditional lead selection, proposed) | 6.8% [5.1%, 9.1%] | 0.0% [0.0%, 0.6%] |
| block-coloc (no fallback) | 7.3% [5.5%, 9.7%] | 0.0% [0.0%, 0.6%] |
| block-coloc (marginal lead selection) | 6.8% [5.1%, 9.1%] | 0.0% [0.0%, 0.6%] |
| block-coloc (original handoff) | 7.3% [5.5%, 9.7%] | 0.0% [0.0%, 0.6%] |
| coloc.abf, whole region | 7.0% [5.2%, 9.3%] | 0.0% [0.0%, 0.6%] |
| coloc.susie | 0.2% [0.0%, 0.9%] | 0.0% [0.0%, 0.6%] |
| coloc.susie + coloc.abf fallback | 7.2% [5.4%, 9.5%] | 0.0% [0.0%, 0.6%] |
| coloc.signals (mask) | 4.0% [2.7%, 5.9%] | 0.0% [0.0%, 0.6%] |
| coloc.signals (cond) | 3.7% [2.4%, 5.5%] | 0.0% [0.0%, 0.6%] |

## 2. Pooled over scenarios (test windows)

TPR pools the shared scenarios (single_shared, single_double_shared, double_shared2, double_shared1); FPR pools single_distinct and double_distinct. AUC does not depend on any cut-off.

| method | lead_p | n_pos | TPR | n_neg | FPR | AUC | no_call | branch_block |
|---|---|---|---|---|---|---|---|---|
| block-coloc (conditional lead selection, proposed) | 1e-03 | 2400 | 39.7% [37.8%, 41.7%] | 1200 | 6.2% [5.0%, 7.8%] | 0.708 [0.690, 0.723] | 0.0% | 18.3% |
| block-coloc (no fallback) | 1e-03 | 2400 | 38.2% [36.2%, 40.1%] | 1200 | 5.5% [4.3%, 6.9%] | 0.705 [0.688, 0.722] | 10.4% | 18.3% |
| block-coloc (marginal lead selection) | 1e-04 | 2400 | 38.3% [36.4%, 40.3%] | 1200 | 5.8% [4.6%, 7.3%] | 0.706 [0.690, 0.722] | 7.4% | 14.3% |
| block-coloc (original handoff) | 1e-04 | 2400 | 38.9% [37.0%, 40.9%] | 1200 | 6.8% [5.5%, 8.3%] | 0.688 [0.671, 0.704] | 11.0% | 35.8% |
| coloc.abf, whole region | - | 2400 | 37.8% [35.9%, 39.8%] | 1200 | 6.1% [4.9%, 7.6%] | 0.698 [0.680, 0.714] | 0.0% | 0.0% |
| coloc.susie | - | 2400 | 13.3% [12.0%, 14.7%] | 1200 | 4.4% [3.4%, 5.7%] | 0.512 [0.499, 0.524] | 83.3% | 0.0% |
| coloc.susie + coloc.abf fallback | - | 2400 | 38.2% [36.3%, 40.2%] | 1200 | 5.8% [4.6%, 7.3%] | 0.698 [0.680, 0.714] | 0.0% | 0.0% |
| coloc.signals (mask) | - | 2400 | 39.5% [37.5%, 41.4%] | 1200 | 6.0% [4.8%, 7.5%] | 0.696 [0.680, 0.713] | 0.0% | 0.0% |
| coloc.signals (cond) | - | 2400 | 39.2% [37.3%, 41.2%] | 1200 | 6.0% [4.8%, 7.5%] | 0.690 [0.672, 0.707] | 0.0% | 0.0% |

## 3. Cut-offs chosen on calibration windows

lead_p = lead threshold (GWAS / eQTL when they differ). '*' marks the lead threshold kept for each block-coloc version.

| method | lead_p | cutoff | calib_FPR | calib_TPR |
|---|---|---|---|---|
| block-coloc (conditional lead selection, proposed) | 1e-03 * | 0.0900 | 5.0% | 41.8% |
| block-coloc (conditional lead selection, proposed) | 1e-04 | 0.0893 | 5.0% | 41.0% |
| block-coloc (conditional lead selection, proposed) | 1e-05 | 0.0893 | 5.0% | 40.2% |
| block-coloc (conditional lead selection, proposed) [earlier bulk thresholds] | 5e-08 / 1e-05 | 0.0893 | 5.0% | 40.0% |
| block-coloc (no fallback) | 1e-03 * | 0.0867 | 5.0% | 40.2% |
| block-coloc (no fallback) | 1e-04 | 0.0884 | 5.0% | 40.7% |
| block-coloc (no fallback) | 1e-05 | 0.0893 | 5.0% | 40.2% |
| block-coloc (no fallback) [earlier bulk thresholds] | 5e-08 / 1e-05 | 0.0893 | 5.0% | 40.0% |
| block-coloc (marginal lead selection) | 1e-03 | 0.0974 | 5.0% | 38.4% |
| block-coloc (marginal lead selection) | 1e-04 * | 0.0898 | 5.0% | 40.6% |
| block-coloc (marginal lead selection) | 1e-05 | 0.0898 | 5.0% | 40.2% |
| block-coloc (marginal lead selection) [earlier bulk thresholds] | 5e-08 / 1e-05 | 0.0898 | 5.0% | 40.0% |
| block-coloc (original handoff) | 1e-03 | 0.1107 | 5.0% | 39.7% |
| block-coloc (original handoff) | 1e-04 * | 0.0939 | 5.0% | 40.2% |
| block-coloc (original handoff) | 1e-05 | 0.0938 | 5.0% | 39.8% |
| block-coloc (original handoff) [earlier bulk thresholds] | 5e-08 / 1e-05 | 0.0893 | 5.0% | 40.8% |
| block-coloc (marginal lead selection) + whole-region fallback | 1e-03 | 0.1095 | 5.0% | 39.9% |
| block-coloc (marginal lead selection) + whole-region fallback | 1e-04 * | 0.0916 | 5.0% | 40.9% |
| block-coloc (marginal lead selection) + whole-region fallback | 1e-05 | 0.0900 | 5.0% | 40.3% |
| block-coloc (marginal lead selection) + whole-region fallback [earlier bulk thresholds] | 5e-08 / 1e-05 | 0.0898 | 5.0% | 40.1% |
| block-coloc (original handoff) + whole-region fallback | 1e-03 | 0.1224 | 5.0% | 39.8% |
| block-coloc (original handoff) + whole-region fallback | 1e-04 * | 0.0992 | 5.0% | 40.9% |
| block-coloc (original handoff) + whole-region fallback | 1e-05 | 0.0974 | 5.0% | 40.2% |
| block-coloc (original handoff) + whole-region fallback [earlier bulk thresholds] | 5e-08 / 1e-05 | 0.0938 | 5.0% | 40.1% |
| coloc.abf, whole region | - | 0.0893 | 5.0% | 39.7% |
| coloc.susie | - | 0.0033 | 5.0% | 14.1% |
| coloc.susie + coloc.abf fallback | - | 0.0877 | 5.0% | 39.9% |
| coloc.signals (mask) | - | 0.1079 | 5.0% | 40.4% |
| coloc.signals (cond) | - | 0.1220 | 5.0% | 40.1% |

## 4. Paired differences on the same test loci (McNemar)

Difference = method minus reference, in percentage points, with 95% CI and exact McNemar p-value.

| method | vs | TPR_difference | FPR_difference |
|---|---|---|---|
| block-coloc (conditional lead selection, proposed) | coloc.susie + coloc.abf fallback | +1.5 pts [+0.9, +2.0], p=1.8e-07 | +0.4 pts [+0.1, +0.8], p=0.063 |
| block-coloc (no fallback) | coloc.susie + coloc.abf fallback | -0.1 pts [-0.9, +0.7], p=0.92 | -0.3 pts [-0.8, +0.2], p=0.34 |
| block-coloc (marginal lead selection) | coloc.susie + coloc.abf fallback | +0.1 pts [-0.5, +0.7], p=0.89 | +0.0 pts [-0.5, +0.5], p=1 |
| block-coloc (original handoff) | coloc.susie + coloc.abf fallback | +0.7 pts [-0.2, +1.5], p=0.16 | +0.9 pts [+0.0, +1.8], p=0.071 |
| coloc.abf, whole region | coloc.susie + coloc.abf fallback | -0.4 pts [-0.9, +0.0], p=0.099 | +0.2 pts [-0.0, +0.5], p=0.25 |
| coloc.susie | coloc.susie + coloc.abf fallback | -25.0 pts [-26.9, -23.1], p=<2e-16 | -1.4 pts [-3.2, +0.4], p=0.15 |
| coloc.signals (mask) | coloc.susie + coloc.abf fallback | +1.2 pts [+0.3, +2.1], p=0.01 | +0.2 pts [-0.8, +1.2], p=0.87 |
| coloc.signals (cond) | coloc.susie + coloc.abf fallback | +1.0 pts [-0.0, +2.0], p=0.067 | +0.2 pts [-1.0, +1.3], p=0.89 |
| block-coloc (conditional lead selection, proposed) | coloc.abf, whole region | +1.9 pts [+1.3, +2.5], p=2e-11 | +0.2 pts [-0.1, +0.4], p=0.5 |
| block-coloc (no fallback) | coloc.abf, whole region | +0.3 pts [-0.5, +1.2], p=0.5 | -0.6 pts [-1.2, +0.0], p=0.092 |
| block-coloc (marginal lead selection) | coloc.abf, whole region | +0.5 pts [-0.1, +1.1], p=0.16 | -0.2 pts [-0.8, +0.3], p=0.55 |
| block-coloc (original handoff) | coloc.abf, whole region | +1.1 pts [+0.2, +2.0], p=0.02 | +0.7 pts [-0.2, +1.5], p=0.18 |
| coloc.susie | coloc.abf, whole region | -24.5 pts [-26.5, -22.6], p=<2e-16 | -1.7 pts [-3.5, +0.1], p=0.082 |
| coloc.susie + coloc.abf fallback | coloc.abf, whole region | +0.4 pts [-0.0, +0.9], p=0.099 | -0.2 pts [-0.5, +0.0], p=0.25 |
| coloc.signals (mask) | coloc.abf, whole region | +1.6 pts [+0.8, +2.5], p=0.00035 | -0.1 pts [-1.0, +0.9], p=1 |
| coloc.signals (cond) | coloc.abf, whole region | +1.4 pts [+0.4, +2.5], p=0.0098 | -0.1 pts [-1.2, +1.1], p=1 |

### 4b. Adopted block-coloc rule and the no-fallback sensitivity analysis

The proposed method uses whole-region PP4 when no valid block is available; 'block-coloc (no fallback)' keeps those loci as no result (section 8).

| method | vs | TPR_difference | FPR_difference |
|---|---|---|---|
| block-coloc (conditional lead selection, proposed) | coloc.abf, whole region | +1.9 pts [+1.3, +2.5], p=2e-11 | +0.2 pts [-0.1, +0.4], p=0.5 |
| block-coloc (conditional lead selection, proposed) | coloc.susie + coloc.abf fallback | +1.5 pts [+0.9, +2.0], p=1.8e-07 | +0.4 pts [+0.1, +0.8], p=0.063 |
| block-coloc (conditional lead selection, proposed) | coloc.signals (cond) | +0.5 pts [-0.4, +1.3], p=0.36 | +0.2 pts [-0.9, +1.4], p=0.77 |
| block-coloc (conditional lead selection, proposed) | coloc.signals (mask) | +0.2 pts [-0.5, +1.0], p=0.59 | +0.2 pts [-0.7, +1.2], p=0.74 |
| block-coloc (no fallback) | coloc.signals (cond) | -1.1 pts [-2.1, -0.1], p=0.046 | -0.5 pts [-1.7, +0.7], p=0.5 |
| block-coloc (no fallback) | coloc.signals (mask) | -1.3 pts [-2.2, -0.3], p=0.0096 | -0.5 pts [-1.6, +0.6], p=0.45 |

## 5. Branch choice: split into blocks or not (test windows)

should_go: block when either trait has 2 causal SNPs, whole otherwise. leads_match_truth: the method's lead counts
(GWAS, eQTL) equal the true numbers of causal SNPs.

| method | scenario | true_signals | should_go | n | branch_block | correct_branch | leads_match_truth | median_leads |
|---|---|---|---|---|---|---|---|---|
| block-coloc (conditional lead selection, proposed) | single_shared | 1+1 | whole | 600 | 8.2% | 91.8% | 33.3% | 1+1 |
| block-coloc (conditional lead selection, proposed) | single_double_shared | 1+2 | block | 600 | 11.7% | 11.7% | 8.0% | 1+1 |
| block-coloc (conditional lead selection, proposed) | double_shared2 | 2+2 | block | 600 | 26.0% | 26.0% | 3.2% | 1+1 |
| block-coloc (conditional lead selection, proposed) | double_shared1 | 2+2 | block | 600 | 29.2% | 29.2% | 4.5% | 1+1 |
| block-coloc (conditional lead selection, proposed) | single_distinct | 1+1 | whole | 600 | 8.8% | 91.2% | 31.2% | 1+1 |
| block-coloc (conditional lead selection, proposed) | double_distinct | 2+2 | block | 600 | 26.0% | 26.0% | 3.5% | 1+1 |
| block-coloc (conditional lead selection, proposed) | gwas_only | 1+0 | whole | 600 | 2.0% | 98.0% | 57.2% | 1+0 |
| block-coloc (conditional lead selection, proposed) | null_both | 0+0 | whole | 600 | 0.3% | 99.7% | 74.5% | 0+0 |
| block-coloc (no fallback) | single_shared | 1+1 | whole | 600 | 8.2% | 91.8% | 33.3% | 1+1 |
| block-coloc (no fallback) | single_double_shared | 1+2 | block | 600 | 11.7% | 11.7% | 8.0% | 1+1 |
| block-coloc (no fallback) | double_shared2 | 2+2 | block | 600 | 26.0% | 26.0% | 3.2% | 1+1 |
| block-coloc (no fallback) | double_shared1 | 2+2 | block | 600 | 29.2% | 29.2% | 4.5% | 1+1 |
| block-coloc (no fallback) | single_distinct | 1+1 | whole | 600 | 8.8% | 91.2% | 31.2% | 1+1 |
| block-coloc (no fallback) | double_distinct | 2+2 | block | 600 | 26.0% | 26.0% | 3.5% | 1+1 |
| block-coloc (no fallback) | gwas_only | 1+0 | whole | 600 | 2.0% | 98.0% | 57.2% | 1+0 |
| block-coloc (no fallback) | null_both | 0+0 | whole | 600 | 0.3% | 99.7% | 74.5% | 0+0 |
| block-coloc (marginal lead selection) | single_shared | 1+1 | whole | 600 | 10.8% | 89.2% | 14.5% | 1+0 |
| block-coloc (marginal lead selection) | single_double_shared | 1+2 | block | 600 | 10.8% | 10.8% | 2.3% | 1+0 |
| block-coloc (marginal lead selection) | double_shared2 | 2+2 | block | 600 | 16.3% | 16.3% | 1.2% | 1+0 |
| block-coloc (marginal lead selection) | double_shared1 | 2+2 | block | 600 | 19.3% | 19.3% | 1.8% | 1+0 |
| block-coloc (marginal lead selection) | single_distinct | 1+1 | whole | 600 | 12.2% | 87.8% | 12.5% | 1+0 |
| block-coloc (marginal lead selection) | double_distinct | 2+2 | block | 600 | 16.5% | 16.5% | 1.3% | 1+0 |
| block-coloc (marginal lead selection) | gwas_only | 1+0 | whole | 600 | 0.7% | 99.3% | 36.8% | 1+0 |
| block-coloc (marginal lead selection) | null_both | 0+0 | whole | 600 | 0.0% | 100.0% | 96.3% | 0+0 |
| block-coloc (original handoff) | single_shared | 1+1 | whole | 600 | 50.5% | 49.5% | 40.8% | 1+1 |
| block-coloc (original handoff) | single_double_shared | 1+2 | block | 600 | 37.7% | 37.7% | 15.0% | 1+1 |
| block-coloc (original handoff) | double_shared2 | 2+2 | block | 600 | 46.7% | 46.7% | 7.5% | 1+1 |
| block-coloc (original handoff) | double_shared1 | 2+2 | block | 600 | 39.0% | 39.0% | 8.8% | 2+1 |
| block-coloc (original handoff) | single_distinct | 1+1 | whole | 600 | 19.3% | 80.7% | 34.7% | 1+1 |
| block-coloc (original handoff) | double_distinct | 2+2 | block | 600 | 21.5% | 21.5% | 8.5% | 1+1 |
| block-coloc (original handoff) | gwas_only | 1+0 | whole | 600 | 28.5% | 71.5% | 0.0% | 1+2 |
| block-coloc (original handoff) | null_both | 0+0 | whole | 600 | 33.0% | 67.0% | 0.0% | 2+2 |

## 6. Detection when the shared signal is not the stronger one (test windows)

single_double_shared: the eQTL has two causal variants, one shared with the GWAS and one eQTL-only. 'stronger' compares their observed eQTL |z|.

| method | n_shared_stronger | TPR_shared_stronger | n_other_stronger | TPR_other_stronger |
|---|---|---|---|---|
| block-coloc (conditional lead selection, proposed) | 282 | 56.7% [50.9%, 62.4%] | 318 | 15.4% [11.9%, 19.8%] |
| block-coloc (no fallback) | 282 | 55.7% [49.8%, 61.4%] | 318 | 13.8% [10.5%, 18.1%] |
| block-coloc (marginal lead selection) | 282 | 56.4% [50.5%, 62.0%] | 318 | 13.8% [10.5%, 18.1%] |
| block-coloc (original handoff) | 282 | 55.7% [49.8%, 61.4%] | 318 | 14.5% [11.0%, 18.8%] |
| coloc.abf, whole region | 282 | 56.7% [50.9%, 62.4%] | 318 | 13.8% [10.5%, 18.1%] |
| coloc.susie | 282 | 13.8% [10.3%, 18.3%] | 318 | 5.3% [3.4%, 8.4%] |
| coloc.susie + coloc.abf fallback | 282 | 56.7% [50.9%, 62.4%] | 318 | 13.8% [10.5%, 18.1%] |
| coloc.signals (mask) | 282 | 55.7% [49.8%, 61.4%] | 318 | 14.2% [10.7%, 18.4%] |
| coloc.signals (cond) | 282 | 53.2% [47.4%, 58.9%] | 318 | 13.5% [10.2%, 17.7%] |

double_shared1: each trait has one shared and one trait-specific causal variant. Columns: in how many traits the shared variant has the larger observed |z|.

| method | shared_stronger_in_both | in_one_trait | in_neither |
|---|---|---|---|
| block-coloc (conditional lead selection, proposed) | 70.6% [63.2%, 77.1%] (n=160) | 22.8% [18.4%, 28.0%] (n=289) | 6.6% [3.6%, 11.8%] (n=151) |
| block-coloc (no fallback) | 67.5% [59.9%, 74.3%] (n=160) | 19.7% [15.5%, 24.7%] (n=289) | 5.3% [2.7%, 10.1%] (n=151) |
| block-coloc (marginal lead selection) | 68.1% [60.6%, 74.8%] (n=160) | 20.4% [16.2%, 25.4%] (n=289) | 5.3% [2.7%, 10.1%] (n=151) |
| block-coloc (original handoff) | 68.1% [60.6%, 74.8%] (n=160) | 23.9% [19.3%, 29.1%] (n=289) | 6.6% [3.6%, 11.8%] (n=151) |
| coloc.abf, whole region | 70.6% [63.2%, 77.1%] (n=160) | 15.6% [11.8%, 20.2%] (n=289) | 5.3% [2.7%, 10.1%] (n=151) |
| coloc.susie | 21.9% [16.2%, 28.9%] (n=160) | 11.4% [8.2%, 15.6%] (n=289) | 4.6% [2.3%, 9.3%] (n=151) |
| coloc.susie + coloc.abf fallback | 70.6% [63.2%, 77.1%] (n=160) | 18.3% [14.3%, 23.2%] (n=289) | 6.0% [3.2%, 10.9%] (n=151) |
| coloc.signals (mask) | 68.8% [61.2%, 75.4%] (n=160) | 24.9% [20.3%, 30.2%] (n=289) | 9.9% [6.1%, 15.7%] (n=151) |
| coloc.signals (cond) | 69.4% [61.8%, 76.0%] (n=160) | 26.3% [21.6%, 31.7%] (n=289) | 7.9% [4.6%, 13.4%] (n=151) |

## 7. Localisation among detected test positives

Sub-region = best block (block branch) or the two credible sets of the best SuSiE signal pair.
Whole-locus calls and coloc.signals calls have no sub-region and are not counted here.

| method | detected | with_subregion | contains_shared | median_size_snps |
|---|---|---|---|---|
| block-coloc (conditional lead selection, proposed) | 953 | 270 | 100.0% [98.6%, 100.0%] | 43 |
| block-coloc (no fallback) | 916 | 271 | 100.0% [98.6%, 100.0%] | 43 |
| block-coloc (marginal lead selection) | 920 | 238 | 96.6% [93.5%, 98.3%] | 49 |
| block-coloc (original handoff) | 934 | 699 | 96.4% [94.8%, 97.6%] | 54 |
| coloc.abf, whole region | 908 |   0 | NA | NA |
| coloc.susie | 319 | 319 | 97.2% [94.7%, 98.5%] | 28 |
| coloc.susie + coloc.abf fallback | 918 | 272 | 99.6% [97.9%, 99.9%] | 24 |
| coloc.signals (mask) | 947 |   0 | NA | NA |
| coloc.signals (cond) | 942 |   0 | NA | NA |

## 8. Sensitivity: fallback policy and the earlier bulk lead thresholds (test windows)

The proposed conditional method includes whole-region fallback; 'block-coloc (no fallback)' omits it.
'+ whole-region fallback' = when the block branch builds no valid block, use the whole-region PP4 instead of no result.
'[earlier bulk thresholds]' = lead thresholds GWAS 5e-8, eQTL 1e-5 (August bulk setting), each with its own calibrated cut-off.

| method | lead_p | cutoff | n_pos | TPR | n_neg | FPR | AUC | no_call | branch_block |
|---|---|---|---|---|---|---|---|---|---|
| block-coloc (conditional lead selection, proposed) | 1e-03 | 0.0900 | 2400 | 39.7% [37.8%, 41.7%] | 1200 | 6.2% [5.0%, 7.8%] | 0.708 | 0.0% | 18.3% |
| block-coloc (no fallback) | 1e-03 | 0.0867 | 2400 | 38.2% [36.2%, 40.1%] | 1200 | 5.5% [4.3%, 6.9%] | 0.705 | 10.4% | 18.3% |
| block-coloc (marginal lead selection) | 1e-04 | 0.0898 | 2400 | 38.3% [36.4%, 40.3%] | 1200 | 5.8% [4.6%, 7.3%] | 0.706 | 7.4% | 14.3% |
| block-coloc (original handoff) | 1e-04 | 0.0939 | 2400 | 38.9% [37.0%, 40.9%] | 1200 | 6.8% [5.5%, 8.3%] | 0.688 | 11.0% | 35.8% |
| block-coloc (marginal lead selection) + whole-region fallback | 1e-04 | 0.0916 | 2400 | 39.0% [37.1%, 41.0%] | 1200 | 6.3% [5.1%, 7.9%] | 0.704 | 0.0% | 14.3% |
| block-coloc (original handoff) + whole-region fallback | 1e-04 | 0.0992 | 2400 | 39.8% [37.9%, 41.8%] | 1200 | 6.8% [5.5%, 8.3%] | 0.705 | 0.0% | 35.8% |
| block-coloc (original handoff) [earlier bulk thresholds] | 5e-08 / 1e-05 | 0.0893 | 2400 | 39.0% [37.0%, 40.9%] | 1200 | 6.7% [5.4%, 8.2%] | 0.697 | 5.8% | 28.2% |
| block-coloc (marginal lead selection) [earlier bulk thresholds] | 5e-08 / 1e-05 | 0.0898 | 2400 | 38.3% [36.4%, 40.3%] | 1200 | 6.1% [4.9%, 7.6%] | 0.702 | 1.8% | 4.0% |
| block-coloc (no fallback) [earlier bulk thresholds] | 5e-08 / 1e-05 | 0.0893 | 2400 | 38.1% [36.2%, 40.0%] | 1200 | 6.1% [4.9%, 7.6%] | 0.700 | 0.8% | 1.8% |
| block-coloc (conditional lead selection, proposed) [earlier bulk thresholds] | 5e-08 / 1e-05 | 0.0893 | 2400 | 38.2% [36.2%, 40.1%] | 1200 | 6.1% [4.9%, 7.6%] | 0.700 | 0.0% | 1.8% |

## 9. Inputs (test windows)

Expected causal |z| (normal effects: E|z| = sqrt(2 / pi * (1 + N q / (1 - q))), q = h2 / K): GWAS 6.2 (K = 1) / 4.4 (K = 2); eQTL 3.2 (K = 1) / 2.3 (K = 2). Table: medians over loci and causal SNPs.

| scenario | n | causal_z_gwas | causal_z_eqtl | pair_abs_r_median | gwas_no_cs | eqtl_no_cs |
|---|---|---|---|---|---|---|
| single_shared | 600 | 5.1 | 3.0 | NA | 41.8% | 68.3% |
| single_double_shared | 600 | 5.3 | 1.9 | 0.03 | 39.2% | 75.0% |
| double_shared2 | 600 | 3.7 | 2.0 | 0.03 | 32.2% | 74.7% |
| double_shared1 | 600 | 3.9 | 2.1 | 0.03 | 31.5% | 72.8% |
| single_distinct | 600 | 5.7 | 2.6 | 0.03 | 39.8% | 73.0% |
| double_distinct | 600 | 3.6 | 2.0 | 0.03 | 33.2% | 73.8% |
| gwas_only | 600 | 5.4 | NA | NA | 40.0% | 99.8% |
| null_both | 600 | NA | NA | NA | 99.5% | 99.8% |

## Figures

- Main figure: `matched_fpr_by_scenario.png`
- PP4 distribution by scenario and method: `pp4_distribution.png`
- `roc_test.png`

## Reproduce

```
Rscript scripts/bulk_h2ref_v1/run.R --stage main --effect normal --cores 4 && Rscript scripts/bulk_h2ref_v1/evaluate.R --stage main --out outputs/bulk_h2ref_v1_normal
```

## Files changed (git diff --stat)

(Codex pastes the output of `git diff --stat HEAD~1` here after committing)
