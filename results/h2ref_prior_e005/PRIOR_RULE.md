# Block priors and decision rules (GWAS h2 0.003, eQTL h2 0.05, fixed effects, low LD)

Loci: 4800 (2400 test). Same simulated loci as stage 1 (same seeds). block-coloc uses the whole-region fallback.
block-coloc, region-equivalent block priors: p1, p2, p12 multiplied by M / M_B inside a block of M_B SNPs (region of M SNPs).

## Three-line conclusion

1. **块内先验不是增益的来源**：eQTL h² 5%、误报率在本格校准窗口上匹配到 5% 时，区域等价先验下 block-coloc 比整区域 coloc.abf 高 7.7 个百分点（配对 95% CI +6.1 至 +9.2，p = 1e-24；原先验 +7.8），两种先验无差别（−0.2，p = 0.5）。本格的判定线很低（coloc.abf 0.012，block-coloc 0.014），因为 eQTL 很强时不共定位位点的 PP4 几乎为 0；因此 gwas_only 位点（PP4 约 0.09）在基于整区域 coloc.abf 的方法下全部越线（100%），coloc.signals 为 4.3%–14.0%。
2. **traceCB 规则下 block-coloc 的结论最明确**：原先验的 block-coloc 把 92.3% 的共定位位点判为共定位（coloc.abf 68.2%，并有 20.4% 被错判为独立；coloc.signals (cond) 92.8%），把 98.3% 的 single_distinct 和 88.0% 的 double_distinct 判为独立（coloc.signals (cond) 47.7% 与 19.0%，其余多为不确定），不共定位判为共定位 ≤ 0.3%，gwas_only 为 0%。区域等价先验下判为共定位的只有 49.8%（块内 PP4 中位数 0.66，整区域 0.88）。
3. **结论**：与默认格一致，块内保留每 SNP 先验；固定的 traceCB 规则在这一格比“按误报率匹配定线”更稳健（gwas_only 0% 对 100%）。

## 1. Matched 5% FPR (cut-off on calibration windows; test windows)

| method | cutoff | TPR | FPR | gwas_only | null_both | AUC | single_shared | single_double_shared | double_shared2 | double_shared1 | single_distinct | double_distinct |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| coloc.abf (whole region) | 0.0116 | 89.8% | 5.0% | 100.0% | 1.3% | 0.971 | 100.0% | 90.3% | 94.7% | 74.3% | 0.0% | 10.0% |
| block-coloc, per-SNP priors (stages 1-2) | 0.014 | 97.7% | 6.8% | 100.0% | 1.0% | 0.992 | 100.0% | 97.3% | 98.0% | 95.3% | 1.7% | 12.0% |
| block-coloc, region-equivalent block priors | 0.014 | 97.5% | 6.7% | 100.0% | 1.0% | 0.990 | 100.0% | 97.3% | 97.7% | 95.0% | 1.7% | 11.7% |
| coloc.signals (mask) | 0.0916 | 94.8% | 6.2% | 14.0% | 0.0% | 0.980 | 99.7% | 92.3% | 97.7% | 89.3% | 4.3% | 8.0% |
| coloc.signals (cond) | 0.1521 | 98.2% | 7.5% | 4.3% | 0.0% | 0.990 | 99.7% | 97.3% | 99.3% | 96.3% | 5.7% | 9.3% |

## 2. traceCB rule (PP4 > 0.7 colocalised, PP3 > 0.7 independent; test windows)

| method | group | n | colocalised | independent | undetermined |
|---|---|---|---|---|---|
| coloc.abf (whole region) | shared (4 scenarios) | 1200 | 68.2% | 20.4% | 11.3% |
| coloc.abf (whole region) | single_distinct | 300 | 0.0% | 100.0% | 0.0% |
| coloc.abf (whole region) | double_distinct | 300 | 0.0% | 90.0% | 10.0% |
| coloc.abf (whole region) | gwas_only | 300 | 0.0% | 0.0% | 100.0% |
| coloc.abf (whole region) | null_both | 300 | 0.0% | 0.0% | 100.0% |
| block-coloc, per-SNP priors (stages 1-2) | shared (4 scenarios) | 1200 | 92.3% | 4.8% | 2.8% |
| block-coloc, per-SNP priors (stages 1-2) | single_distinct | 300 | 0.3% | 98.3% | 1.3% |
| block-coloc, per-SNP priors (stages 1-2) | double_distinct | 300 | 0.0% | 88.0% | 12.0% |
| block-coloc, per-SNP priors (stages 1-2) | gwas_only | 300 | 0.0% | 0.0% | 100.0% |
| block-coloc, per-SNP priors (stages 1-2) | null_both | 300 | 0.0% | 0.0% | 100.0% |
| block-coloc, region-equivalent block priors | shared (4 scenarios) | 1200 | 49.8% | 6.8% | 43.4% |
| block-coloc, region-equivalent block priors | single_distinct | 300 | 0.0% | 98.3% | 1.7% |
| block-coloc, region-equivalent block priors | double_distinct | 300 | 0.0% | 88.7% | 11.3% |
| block-coloc, region-equivalent block priors | gwas_only | 300 | 0.0% | 0.0% | 100.0% |
| block-coloc, region-equivalent block priors | null_both | 300 | 0.0% | 0.0% | 100.0% |
| coloc.signals (mask) | shared (4 scenarios) | 1200 | 86.8% | 5.4% | 7.8% |
| coloc.signals (mask) | single_distinct | 300 | 0.3% | 70.7% | 29.0% |
| coloc.signals (mask) | double_distinct | 300 | 0.0% | 34.3% | 65.7% |
| coloc.signals (mask) | gwas_only | 300 | 0.0% | 0.0% | 100.0% |
| coloc.signals (mask) | null_both | 300 | 0.0% | 0.0% | 100.0% |
| coloc.signals (cond) | shared (4 scenarios) | 1200 | 92.8% | 1.3% | 5.9% |
| coloc.signals (cond) | single_distinct | 300 | 1.3% | 47.7% | 51.0% |
| coloc.signals (cond) | double_distinct | 300 | 1.0% | 19.0% | 80.0% |
| coloc.signals (cond) | gwas_only | 300 | 0.0% | 0.0% | 100.0% |
| coloc.signals (cond) | null_both | 300 | 0.0% | 0.0% | 100.0% |

## 3. Paired differences on the same test loci (percentage points, McNemar)

| rule | method | vs | part | diff | lo | hi | p |
|---|---|---|---|---|---|---|---|
| matched 5% FPR | block-coloc, region-equivalent block priors | coloc.abf (whole region) | TPR | +7.7 | +6.1 | +9.2 | 9.9e-25 |
| matched 5% FPR | block-coloc, region-equivalent block priors | coloc.abf (whole region) | FPR | +1.7 | +0.6 | +2.7 | 0.002 |
| matched 5% FPR | block-coloc, per-SNP priors (stages 1-2) | coloc.abf (whole region) | TPR | +7.8 | +6.3 | +9.4 | 3.1e-26 |
| matched 5% FPR | block-coloc, per-SNP priors (stages 1-2) | coloc.abf (whole region) | FPR | +1.8 | +0.8 | +2.9 | 0.00098 |
| matched 5% FPR | block-coloc, region-equivalent block priors | block-coloc, per-SNP priors (stages 1-2) | TPR | -0.2 | -0.4 | +0.1 | 0.5 |
| matched 5% FPR | block-coloc, region-equivalent block priors | block-coloc, per-SNP priors (stages 1-2) | FPR | -0.2 | -0.5 | +0.2 | 1 |
| matched 5% FPR | block-coloc, region-equivalent block priors | coloc.signals (cond) | TPR | -0.7 | -1.5 | +0.2 | 0.17 |
| matched 5% FPR | block-coloc, region-equivalent block priors | coloc.signals (cond) | FPR | -0.8 | -3.4 | +1.7 | 0.61 |
| matched 5% FPR | block-coloc, region-equivalent block priors | coloc.signals (mask) | TPR | +2.8 | +1.6 | +3.9 | 2e-06 |
| matched 5% FPR | block-coloc, region-equivalent block priors | coloc.signals (mask) | FPR | +0.5 | -1.8 | +2.8 | 0.78 |
| traceCB PP4 > 0.7 | block-coloc, region-equivalent block priors | coloc.abf (whole region) | TPR | -18.5 | -21.6 | -15.4 | 5.2e-29 |
| traceCB PP4 > 0.7 | block-coloc, region-equivalent block priors | coloc.abf (whole region) | FPR | +0.0 | +0.0 | +0.0 | 1 |
| traceCB PP4 > 0.7 | block-coloc, per-SNP priors (stages 1-2) | coloc.abf (whole region) | TPR | +24.1 | +21.7 | +26.5 | 1.5e-85 |
| traceCB PP4 > 0.7 | block-coloc, per-SNP priors (stages 1-2) | coloc.abf (whole region) | FPR | +0.2 | -0.2 | +0.5 | 1 |
| traceCB PP4 > 0.7 | block-coloc, region-equivalent block priors | block-coloc, per-SNP priors (stages 1-2) | TPR | -42.6 | -45.4 | -39.8 | 3e-154 |
| traceCB PP4 > 0.7 | block-coloc, region-equivalent block priors | block-coloc, per-SNP priors (stages 1-2) | FPR | -0.2 | -0.5 | +0.2 | 1 |
| traceCB PP4 > 0.7 | block-coloc, region-equivalent block priors | coloc.signals (cond) | TPR | -43.0 | -45.8 | -40.2 | 1.2e-153 |
| traceCB PP4 > 0.7 | block-coloc, region-equivalent block priors | coloc.signals (cond) | FPR | -1.2 | -2.0 | -0.3 | 0.016 |
| traceCB PP4 > 0.7 | block-coloc, region-equivalent block priors | coloc.signals (mask) | TPR | -37.0 | -39.9 | -34.1 | 1.1e-108 |
| traceCB PP4 > 0.7 | block-coloc, region-equivalent block priors | coloc.signals (mask) | FPR | -0.2 | -0.5 | +0.2 | 1 |

## 4. Prior diagnostics (loci that used a block)

block_loci 1748; median_block_snps 32; median_region_snps 938; median_pp4_orig 0.976; median_pp4_adj 0.664; median_pp4_whole 0.877

## 5. Consistency with stage 1 (PP4 should be identical)

branch (whole / block) mismatches for V3_gate_conditional: 0 of 4800; whole_abf: 4800 loci, max |diff| 0.00e+00, NA mismatches 0; V3_gate_conditional: 4800 loci, max |diff| 0.00e+00, NA mismatches 0; signals_mask: 4800 loci, max |diff| 0.00e+00, NA mismatches 0; signals_cond: 4800 loci, max |diff| 0.00e+00, NA mismatches 0
