# Block priors and decision rules (GWAS h2 0.003, eQTL h2 0.015, fixed effects, low LD)

Loci: 9600 (4800 test). Same simulated loci as stage 1 (same seeds). block-coloc uses the whole-region fallback.
block-coloc, region-equivalent block priors: p1, p2, p12 multiplied by M / M_B inside a block of M_B SNPs (region of M SNPs).

## Three-line conclusion

1. **块内先验不是增益的来源**：误报率在校准窗口上匹配到 5% 时，把块内先验调成与整区域等价（p1、p2、p12 乘以 M / M_B）后，block-coloc 仍比整区域 coloc.abf 高 3.0 个百分点（配对 95% CI +2.2 至 +3.9，p = 2e-12；原先验 +3.5），两种先验之间只差 0.5 个百分点（p = 0.02）。
2. **traceCB 规则（PP4 > 0.7 共定位、PP3 > 0.7 独立）下**：原先验的 block-coloc 把 44.1% 的共定位位点判为共定位（coloc.abf 36.5%，配对 +7.5，p = 1e-48；coloc.signals 42.6%–44.9%），不共定位位点判为共定位的比例与 coloc.abf 相同（0.0%–0.2%），gwas_only 为 0.2%。区域等价先验矫枉过正：它把整区域 H3 的先验集中到块内彼此连锁的 SNP 对上，块内 PP4 中位数从 0.87 降到 0.55，低于整区域的 0.72，判为共定位的只有 26.9%。
3. **结论**：块内保留每 SNP 先验；区域等价先验只作敏感性分析；真实数据可用 traceCB 的三分类规则。本格的阴性位点中因果变异间 |r| < 0.2，固定阈值在高 LD 阴性上的表现尚未检验。

## 1. Matched 5% FPR (cut-off on calibration windows; test windows)

| method | cutoff | TPR | FPR | gwas_only | null_both | AUC | single_shared | single_double_shared | double_shared2 | double_shared1 | single_distinct | double_distinct |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| coloc.abf (whole region) | 0.1412 | 73.2% | 3.5% | 5.0% | 0.0% | 0.922 | 96.3% | 74.3% | 76.8% | 45.3% | 3.5% | 3.5% |
| block-coloc, per-SNP priors (stages 1-2) | 0.1448 | 76.8% | 3.8% | 4.8% | 0.0% | 0.934 | 96.2% | 74.2% | 81.8% | 54.8% | 3.7% | 3.8% |
| block-coloc, region-equivalent block priors | 0.1448 | 76.2% | 3.6% | 4.8% | 0.0% | 0.932 | 96.0% | 74.3% | 81.2% | 53.5% | 3.5% | 3.7% |
| coloc.signals (mask) | 0.1715 | 78.2% | 4.6% | 3.5% | 0.0% | 0.940 | 95.5% | 71.3% | 84.3% | 61.8% | 3.7% | 5.5% |
| coloc.signals (cond) | 0.2103 | 76.2% | 3.8% | 2.5% | 0.0% | 0.943 | 94.8% | 67.7% | 82.0% | 60.2% | 3.2% | 4.3% |

## 2. traceCB rule (PP4 > 0.7 colocalised, PP3 > 0.7 independent; test windows)

| method | group | n | colocalised | independent | undetermined |
|---|---|---|---|---|---|
| coloc.abf (whole region) | shared (4 scenarios) | 2400 | 36.5% | 4.1% | 59.3% |
| coloc.abf (whole region) | single_distinct | 600 | 0.2% | 35.7% | 64.2% |
| coloc.abf (whole region) | double_distinct | 600 | 0.0% | 10.2% | 89.8% |
| coloc.abf (whole region) | gwas_only | 600 | 0.2% | 0.2% | 99.7% |
| coloc.abf (whole region) | null_both | 600 | 0.0% | 0.0% | 100.0% |
| block-coloc, per-SNP priors (stages 1-2) | shared (4 scenarios) | 2400 | 44.1% | 2.7% | 53.2% |
| block-coloc, per-SNP priors (stages 1-2) | single_distinct | 600 | 0.2% | 35.5% | 64.3% |
| block-coloc, per-SNP priors (stages 1-2) | double_distinct | 600 | 0.2% | 9.8% | 90.0% |
| block-coloc, per-SNP priors (stages 1-2) | gwas_only | 600 | 0.2% | 0.2% | 99.7% |
| block-coloc, per-SNP priors (stages 1-2) | null_both | 600 | 0.0% | 0.0% | 100.0% |
| block-coloc, region-equivalent block priors | shared (4 scenarios) | 2400 | 26.9% | 3.9% | 69.2% |
| block-coloc, region-equivalent block priors | single_distinct | 600 | 0.2% | 35.5% | 64.3% |
| block-coloc, region-equivalent block priors | double_distinct | 600 | 0.0% | 10.3% | 89.7% |
| block-coloc, region-equivalent block priors | gwas_only | 600 | 0.2% | 0.2% | 99.7% |
| block-coloc, region-equivalent block priors | null_both | 600 | 0.0% | 0.0% | 100.0% |
| coloc.signals (mask) | shared (4 scenarios) | 2400 | 42.6% | 2.1% | 55.3% |
| coloc.signals (mask) | single_distinct | 600 | 0.2% | 27.0% | 72.8% |
| coloc.signals (mask) | double_distinct | 600 | 0.0% | 4.7% | 95.3% |
| coloc.signals (mask) | gwas_only | 600 | 0.2% | 0.2% | 99.7% |
| coloc.signals (mask) | null_both | 600 | 0.0% | 0.0% | 100.0% |
| coloc.signals (cond) | shared (4 scenarios) | 2400 | 44.9% | 1.5% | 53.7% |
| coloc.signals (cond) | single_distinct | 600 | 0.3% | 21.2% | 78.5% |
| coloc.signals (cond) | double_distinct | 600 | 0.0% | 2.5% | 97.5% |
| coloc.signals (cond) | gwas_only | 600 | 0.2% | 0.0% | 99.8% |
| coloc.signals (cond) | null_both | 600 | 0.0% | 0.0% | 100.0% |

## 3. Paired differences on the same test loci (percentage points, McNemar)

| rule | method | vs | part | diff | lo | hi | p |
|---|---|---|---|---|---|---|---|
| matched 5% FPR | block-coloc, region-equivalent block priors | coloc.abf (whole region) | TPR | +3.0 | +2.2 | +3.9 | 1.9e-12 |
| matched 5% FPR | block-coloc, region-equivalent block priors | coloc.abf (whole region) | FPR | +0.1 | -0.3 | +0.4 | 1 |
| matched 5% FPR | block-coloc, per-SNP priors (stages 1-2) | coloc.abf (whole region) | TPR | +3.5 | +2.7 | +4.4 | 2.7e-17 |
| matched 5% FPR | block-coloc, per-SNP priors (stages 1-2) | coloc.abf (whole region) | FPR | +0.2 | -0.1 | +0.6 | 0.38 |
| matched 5% FPR | block-coloc, region-equivalent block priors | block-coloc, per-SNP priors (stages 1-2) | TPR | -0.5 | -0.9 | -0.1 | 0.023 |
| matched 5% FPR | block-coloc, region-equivalent block priors | block-coloc, per-SNP priors (stages 1-2) | FPR | -0.2 | -0.4 | +0.1 | 0.5 |
| matched 5% FPR | block-coloc, region-equivalent block priors | coloc.signals (cond) | TPR | +0.1 | -1.1 | +1.3 | 0.95 |
| matched 5% FPR | block-coloc, region-equivalent block priors | coloc.signals (cond) | FPR | -0.2 | -1.2 | +0.8 | 0.87 |
| matched 5% FPR | block-coloc, region-equivalent block priors | coloc.signals (mask) | TPR | -2.0 | -3.1 | -0.9 | 0.00046 |
| matched 5% FPR | block-coloc, region-equivalent block priors | coloc.signals (mask) | FPR | -1.0 | -1.8 | -0.2 | 0.029 |
| traceCB PP4 > 0.7 | block-coloc, region-equivalent block priors | coloc.abf (whole region) | TPR | -9.7 | -11.0 | -8.3 | 5.5e-47 |
| traceCB PP4 > 0.7 | block-coloc, region-equivalent block priors | coloc.abf (whole region) | FPR | +0.0 | +0.0 | +0.0 | 1 |
| traceCB PP4 > 0.7 | block-coloc, per-SNP priors (stages 1-2) | coloc.abf (whole region) | TPR | +7.5 | +6.5 | +8.6 | 1.3e-48 |
| traceCB PP4 > 0.7 | block-coloc, per-SNP priors (stages 1-2) | coloc.abf (whole region) | FPR | +0.1 | -0.1 | +0.2 | 1 |
| traceCB PP4 > 0.7 | block-coloc, region-equivalent block priors | block-coloc, per-SNP priors (stages 1-2) | TPR | -17.2 | -18.7 | -15.7 | 9.5e-125 |
| traceCB PP4 > 0.7 | block-coloc, region-equivalent block priors | block-coloc, per-SNP priors (stages 1-2) | FPR | -0.1 | -0.2 | +0.1 | 1 |
| traceCB PP4 > 0.7 | block-coloc, region-equivalent block priors | coloc.signals (cond) | TPR | -18.0 | -19.6 | -16.4 | 2.4e-122 |
| traceCB PP4 > 0.7 | block-coloc, region-equivalent block priors | coloc.signals (cond) | FPR | -0.1 | -0.2 | +0.1 | 1 |
| traceCB PP4 > 0.7 | block-coloc, region-equivalent block priors | coloc.signals (mask) | TPR | -15.8 | -17.2 | -14.3 | 7.6e-102 |
| traceCB PP4 > 0.7 | block-coloc, region-equivalent block priors | coloc.signals (mask) | FPR | +0.0 | +0.0 | +0.0 | 1 |

## 4. Prior diagnostics (loci that used a block)

block_loci 1402; median_block_snps 42; median_region_snps 938; median_pp4_orig 0.867; median_pp4_adj 0.552; median_pp4_whole 0.716

## 5. Consistency with stage 1 (PP4 should be identical)

branch (whole / block) mismatches for V3_gate_conditional: 0 of 9600; whole_abf: 9600 loci, max |diff| 0.00e+00, NA mismatches 0; V3_gate_conditional: 9600 loci, max |diff| 0.00e+00, NA mismatches 0; signals_mask: 9600 loci, max |diff| 0.00e+00, NA mismatches 0; signals_cond: 9600 loci, max |diff| 0.00e+00, NA mismatches 0
