# Bulk block-coloc: stage 2, one factor at a time (test windows, frozen cut-offs)

Generated 2026-10-26 16:25. Replicates per cell: 50. Cut-offs and lead thresholds from stage 1 (thresholds.tsv), not re-tuned.

## Three-line conclusion

1. **GWAS 变强时，条件 lead 选择能把误报率稳住**：GWAS q 从 0.1% 增大到 2% 时，block-coloc（条件 lead 选择）的误报率保持在 3.5%–4.3%，而 block-coloc（边际 lead 选择）升到 9.5%、block-coloc（交接包原版）升到 15.7%；整区域 coloc.abf 保持在 4.3%。q = 2% 时它的检出率为 92.7%，高于整区域 coloc.abf（88.0%）和 coloc.susie + fallback（89.7%），低于 coloc.signals (cond)（95.0%）。
2. **eQTL 弱时它最好，eQTL 强时被反超**：eQTL q = 0.75% 时，block-coloc（条件 lead 选择）检出率最高（82.4%；coloc.susie + fallback 80.3%，整区域 coloc.abf 79.2%，coloc.signals (cond) 78.4%）；eQTL q = 5% 时，coloc.signals (cond)（99.3%）和 coloc.susie + fallback（98.9%）超过它（97.3%），整区域 coloc.abf 反而降到 83.9%。
3. **局限**：不同因果变异之间的 LD 升高时，所有方法都明显变差：检出率从约 90%（|r| < 0.2）降到约 50%（|r| 0.63–0.9）；高 LD 时 single_distinct 被 block-coloc（条件 lead 选择）判为共定位的比例达 84.3%（coloc.susie 为 31.0%）。只有 GWAS 有信号时（gwas_only），基于 coloc.abf 的方法误报 4%–17%，block-coloc（条件 lead 选择）为 6.3%–13.3%，coloc.signals 为 1.3%–4.7%，coloc.susie 为 0%；两性状都没有信号时，所有方法误报均为 0%。

## Detection rate vs GWAS q (eQTL q 1.5%, low LD)

| q_gwas | block-coloc, conditional (proposed) | block-coloc, marginal | block-coloc, handoff | coloc.abf | coloc.susie | coloc.susie + fallback | signals (mask) | signals (cond) |
|---|---|---|---|---|---|---|---|---|
| 0.001 | 85.9% | 85.0% | 85.7% | 84.9% | 27.3% | 81.9% | 84.3% | 85.8% |
| 0.003 | 92.6% | 91.3% | 93.7% | 89.9% | 41.9% | 89.6% | 92.8% | 94.8% |
| 0.010 | 92.8% | 89.2% | 92.8% | 88.2% | 42.2% | 89.6% | 92.7% | 95.1% |
| 0.020 | 92.7% | 89.1% | 92.9% | 88.0% | 42.6% | 89.7% | 92.3% | 95.0% |

## False-positive rate vs GWAS q

| q_gwas | block-coloc, conditional (proposed) | block-coloc, marginal | block-coloc, handoff | coloc.abf | coloc.susie | coloc.susie + fallback | signals (mask) | signals (cond) |
|---|---|---|---|---|---|---|---|---|
| 0.001 | 3.5% | 2.5% | 2.5% | 2.8% | 2.7% | 3.3% | 2.8% | 2.7% |
| 0.003 | 4.3% | 4.5% | 4.7% | 4.7% | 3.8% | 5.2% | 4.3% | 4.8% |
| 0.010 | 4.2% | 6.0% | 7.8% | 4.3% | 3.8% | 4.7% | 5.3% | 5.2% |
| 0.020 | 4.2% | 9.5% | 15.7% | 4.3% | 4.0% | 4.8% | 6.8% | 5.5% |

## Detection rate vs eQTL q (GWAS q 0.3%, low LD)

| q_eqtl | block-coloc, conditional (proposed) | block-coloc, marginal | block-coloc, handoff | coloc.abf | coloc.susie | coloc.susie + fallback | signals (mask) | signals (cond) |
|---|---|---|---|---|---|---|---|---|
| 0.0075 | 82.4% | 79.9% | 81.7% | 79.2% | 12.0% | 80.3% | 77.4% | 78.4% |
| 0.0150 | 92.6% | 91.3% | 93.7% | 89.9% | 41.9% | 89.6% | 92.8% | 94.8% |
| 0.0300 | 96.3% | 95.4% | 96.3% | 87.3% | 90.6% | 95.1% | 96.9% | 99.3% |
| 0.0500 | 97.3% | 97.4% | 97.9% | 83.9% | 99.0% | 98.9% | 95.9% | 99.3% |

## False-positive rate vs eQTL q

| q_eqtl | block-coloc, conditional (proposed) | block-coloc, marginal | block-coloc, handoff | coloc.abf | coloc.susie | coloc.susie + fallback | signals (mask) | signals (cond) |
|---|---|---|---|---|---|---|---|---|
| 0.0075 | 8.7% | 7.7% | 7.7% | 7.7% | 0.3% | 8.2% | 5.2% | 5.0% |
| 0.0150 | 4.3% | 4.5% | 4.7% | 4.7% | 3.8% | 5.2% | 4.3% | 4.8% |
| 0.0300 | 0.5% | 1.0% | 1.0% | 0.7% | 2.5% | 0.8% | 2.7% | 2.3% |
| 0.0500 | 0.0% | 0.8% | 0.8% | 0.0% | 0.0% | 0.0% | 2.2% | 2.2% |

## False-positive rate when at most one trait has a signal (extra negatives, not used for calibration)

gwas_only: GWAS has one causal variant, the eQTL has none. null_both: neither trait has a signal.

### vs GWAS q (eQTL q 1.5%)

| scenario | q_gwas | n | block-coloc, conditional (proposed) | block-coloc, marginal | block-coloc, handoff | coloc.abf | coloc.susie | coloc.susie + fallback | signals (mask) | signals (cond) |
|---|---|---|---|---|---|---|---|---|---|---|
| gwas_only | 0.001 | 300 | 6.3% | 4.0% | 4.3% | 4.3% | 0.0% | 4.3% | 1.3% | 1.3% |
| gwas_only | 0.003 | 300 | 13.3% | 10.0% | 10.7% | 10.0% | 0.0% | 11.0% | 4.0% | 4.0% |
| gwas_only | 0.010 | 300 | 13.0% | 10.3% | 12.3% | 10.3% | 0.0% | 11.3% | 4.3% | 4.3% |
| gwas_only | 0.020 | 300 | 13.0% | 10.3% | 17.0% | 10.7% | 0.0% | 11.3% | 4.7% | 4.7% |
| null_both | 0.001 | 300 | 0.0% | 0.0% | 0.0% | 0.0% | 0.0% | 0.0% | 0.0% | 0.0% |
| null_both | 0.003 | 300 | 0.0% | 0.0% | 0.0% | 0.0% | 0.0% | 0.0% | 0.0% | 0.0% |
| null_both | 0.010 | 300 | 0.0% | 0.0% | 0.0% | 0.0% | 0.0% | 0.0% | 0.0% | 0.0% |
| null_both | 0.020 | 300 | 0.0% | 0.0% | 0.0% | 0.0% | 0.0% | 0.0% | 0.0% | 0.0% |

### vs eQTL q (GWAS q 0.3%)

| scenario | q_eqtl | n | block-coloc, conditional (proposed) | block-coloc, marginal | block-coloc, handoff | coloc.abf | coloc.susie | coloc.susie + fallback | signals (mask) | signals (cond) |
|---|---|---|---|---|---|---|---|---|---|---|
| gwas_only | 0.0075 | 300 | 13.3% | 10.0% | 10.7% | 10.0% | 0.0% | 11.0% | 4.0% | 4.0% |
| gwas_only | 0.0150 | 300 | 13.3% | 10.0% | 10.7% | 10.0% | 0.0% | 11.0% | 4.0% | 4.0% |
| gwas_only | 0.0300 | 300 | 13.3% | 10.0% | 10.7% | 10.0% | 0.0% | 11.0% | 4.0% | 4.0% |
| gwas_only | 0.0500 | 300 | 13.3% | 10.0% | 10.7% | 10.0% | 0.0% | 11.0% | 4.0% | 4.0% |
| null_both | 0.0075 | 300 | 0.0% | 0.0% | 0.0% | 0.0% | 0.0% | 0.0% | 0.0% | 0.0% |
| null_both | 0.0150 | 300 | 0.0% | 0.0% | 0.0% | 0.0% | 0.0% | 0.0% | 0.0% | 0.0% |
| null_both | 0.0300 | 300 | 0.0% | 0.0% | 0.0% | 0.0% | 0.0% | 0.0% | 0.0% | 0.0% |
| null_both | 0.0500 | 300 | 0.0% | 0.0% | 0.0% | 0.0% | 0.0% | 0.0% | 0.0% | 0.0% |

## LD between causal SNPs (main strengths; single_double_shared, double_shared2 pooled for TPR; single_distinct, double_distinct for FPR)

low |r| < 0.2; mid 0.32-0.63; high 0.63-0.9.

| tier | block-coloc, conditional (proposed) | block-coloc, marginal | block-coloc, handoff | coloc.abf | coloc.susie | coloc.susie + fallback | signals (mask) | signals (cond) |
|---|---|---|---|---|---|---|---|---|
| low | 89.7% | 88.0% | 91.7% | 85.8% | 46.7% | 85.5% | 91.3% | 94.3% |
| mid | 73.2% | 71.5% | 74.3% | 70.3% | 36.5% | 72.8% | 65.3% | 74.5% |
| high | 50.5% | 48.8% | 50.5% | 49.5% | 39.3% | 50.2% | 45.7% | 48.5% |

| tier | block-coloc, conditional (proposed) | block-coloc, marginal | block-coloc, handoff | coloc.abf | coloc.susie | coloc.susie + fallback | signals (mask) | signals (cond) |
|---|---|---|---|---|---|---|---|---|
| low | 4.3% | 4.5% | 4.7% | 4.7% | 3.8% | 5.2% | 4.3% | 4.8% |
| mid | 18.7% | 16.5% | 18.3% | 16.2% | 7.5% | 16.7% | 12.3% | 14.0% |
| high | 44.0% | 42.7% | 43.7% | 42.3% | 15.8% | 42.7% | 37.2% | 36.8% |

### By scenario

| tier | scenario | n | block-coloc, conditional (proposed) | block-coloc, marginal | block-coloc, handoff | coloc.abf | coloc.susie | coloc.susie + fallback | signals (mask) | signals (cond) |
|---|---|---|---|---|---|---|---|---|---|---|
| low | single_double_shared | 300 | 86.7% | 83.0% | 86.7% | 85.3% | 40.3% | 74.7% | 88.0% | 90.0% |
| low | double_shared2 | 300 | 92.7% | 93.0% | 96.7% | 86.3% | 53.0% | 96.3% | 94.7% | 98.7% |
| low | single_distinct | 300 | 5.7% | 5.0% | 6.0% | 5.3% | 1.3% | 5.3% | 3.0% | 3.3% |
| low | double_distinct | 300 | 3.0% | 4.0% | 3.3% | 4.0% | 6.3% | 5.0% | 5.7% | 6.3% |
| mid | single_double_shared | 300 | 73.0% | 70.3% | 71.7% | 69.7% | 34.0% | 68.7% | 61.7% | 66.0% |
| mid | double_shared2 | 300 | 73.3% | 72.7% | 77.0% | 71.0% | 39.0% | 77.0% | 69.0% | 83.0% |
| mid | single_distinct | 300 | 34.0% | 29.0% | 31.0% | 29.3% | 14.3% | 30.0% | 21.7% | 22.3% |
| mid | double_distinct | 300 | 3.3% | 4.0% | 5.7% | 3.0% | 0.7% | 3.3% | 3.0% | 5.7% |
| high | single_double_shared | 300 | 60.0% | 59.0% | 60.0% | 59.7% | 52.0% | 59.7% | 54.3% | 55.0% |
| high | double_shared2 | 300 | 41.0% | 38.7% | 41.0% | 39.3% | 26.7% | 40.7% | 37.0% | 42.0% |
| high | single_distinct | 300 | 84.3% | 82.3% | 83.3% | 82.0% | 31.0% | 82.0% | 72.3% | 72.0% |
| high | double_distinct | 300 | 3.7% | 3.0% | 4.0% | 2.7% | 0.7% | 3.3% | 2.0% | 1.7% |

### Branch choice by LD tier

| tier | method | scenario | true_signals | should_go | n | branch_block | correct_branch | leads_match_truth | median_leads |
|---|---|---|---|---|---|---|---|---|---|
| low | block-coloc (conditional lead selection, proposed) | single_double_shared | 1+2 | block | 300 | 25.7% | 25.7% | 24.0% | 1+1 |
| low | block-coloc (conditional lead selection, proposed) | double_shared2 | 2+2 | block | 300 | 81.0% | 81.0% | 22.0% | 2+1 |
| low | block-coloc (conditional lead selection, proposed) | single_distinct | 1+1 | whole | 300 | 0.3% | 99.7% | 58.0% | 1+1 |
| low | block-coloc (conditional lead selection, proposed) | double_distinct | 2+2 | block | 300 | 78.0% | 78.0% | 25.0% | 2+1 |
| low | block-coloc (marginal lead selection) | single_double_shared | 1+2 | block | 300 | 22.3% | 22.3% | 15.7% | 1+1 |
| low | block-coloc (marginal lead selection) | double_shared2 | 2+2 | block | 300 | 59.7% | 59.7% | 10.7% | 2+1 |
| low | block-coloc (marginal lead selection) | single_distinct | 1+1 | whole | 300 | 3.7% | 96.3% | 30.3% | 1+0 |
| low | block-coloc (marginal lead selection) | double_distinct | 2+2 | block | 300 | 55.7% | 55.7% | 8.3% | 2+1 |
| low | block-coloc (original handoff) | single_double_shared | 1+2 | block | 300 | 54.3% | 54.3% | 26.3% | 1+1 |
| low | block-coloc (original handoff) | double_shared2 | 2+2 | block | 300 | 91.3% | 91.3% | 17.7% | 2+1 |
| low | block-coloc (original handoff) | single_distinct | 1+1 | whole | 300 | 3.0% | 97.0% | 76.0% | 1+1 |
| low | block-coloc (original handoff) | double_distinct | 2+2 | block | 300 | 28.3% | 28.3% | 18.3% | 2+1 |
| mid | block-coloc (conditional lead selection, proposed) | single_double_shared | 1+2 | block | 300 | 1.7% | 1.7% | 1.3% | 1+1 |
| mid | block-coloc (conditional lead selection, proposed) | double_shared2 | 2+2 | block | 300 | 31.3% | 31.3% | 2.0% | 1+1 |
| mid | block-coloc (conditional lead selection, proposed) | single_distinct | 1+1 | whole | 300 | 1.0% | 99.0% | 61.3% | 1+1 |
| mid | block-coloc (conditional lead selection, proposed) | double_distinct | 2+2 | block | 300 | 30.7% | 30.7% | 1.7% | 1+1 |
| mid | block-coloc (marginal lead selection) | single_double_shared | 1+2 | block | 300 | 13.7% | 13.7% | 6.0% | 1+1 |
| mid | block-coloc (marginal lead selection) | double_shared2 | 2+2 | block | 300 | 41.7% | 41.7% | 3.3% | 2+1 |
| mid | block-coloc (marginal lead selection) | single_distinct | 1+1 | whole | 300 | 6.0% | 94.0% | 33.7% | 1+0 |
| mid | block-coloc (marginal lead selection) | double_distinct | 2+2 | block | 300 | 44.7% | 44.7% | 2.3% | 2+1 |
| mid | block-coloc (original handoff) | single_double_shared | 1+2 | block | 300 | 51.0% | 51.0% | 19.0% | 1+1 |
| mid | block-coloc (original handoff) | double_shared2 | 2+2 | block | 300 | 68.3% | 68.3% | 5.0% | 2+1 |
| mid | block-coloc (original handoff) | single_distinct | 1+1 | whole | 300 | 30.0% | 70.0% | 78.3% | 1+1 |
| mid | block-coloc (original handoff) | double_distinct | 2+2 | block | 300 | 18.7% | 18.7% | 7.7% | 2+1 |
| high | block-coloc (conditional lead selection, proposed) | single_double_shared | 1+2 | block | 300 | 1.3% | 1.3% | 0.7% | 1+1 |
| high | block-coloc (conditional lead selection, proposed) | double_shared2 | 2+2 | block | 300 | 2.7% | 2.7% | 0.0% | 1+1 |
| high | block-coloc (conditional lead selection, proposed) | single_distinct | 1+1 | whole | 300 | 1.7% | 98.3% | 57.3% | 1+1 |
| high | block-coloc (conditional lead selection, proposed) | double_distinct | 2+2 | block | 300 | 2.3% | 2.3% | 0.0% | 1+1 |
| high | block-coloc (marginal lead selection) | single_double_shared | 1+2 | block | 300 | 10.7% | 10.7% | 3.3% | 1+1 |
| high | block-coloc (marginal lead selection) | double_shared2 | 2+2 | block | 300 | 23.3% | 23.3% | 1.3% | 1+0 |
| high | block-coloc (marginal lead selection) | single_distinct | 1+1 | whole | 300 | 2.7% | 97.3% | 33.0% | 1+0 |
| high | block-coloc (marginal lead selection) | double_distinct | 2+2 | block | 300 | 22.0% | 22.0% | 0.7% | 1+0 |
| high | block-coloc (original handoff) | single_double_shared | 1+2 | block | 300 | 42.0% | 42.0% | 20.0% | 1+1 |
| high | block-coloc (original handoff) | double_shared2 | 2+2 | block | 300 | 46.7% | 46.7% | 6.0% | 2+1 |
| high | block-coloc (original handoff) | single_distinct | 1+1 | whole | 300 | 87.3% | 12.7% | 77.0% | 1+1 |
| high | block-coloc (original handoff) | double_distinct | 2+2 | block | 300 | 25.0% | 25.0% | 8.7% | 2+1 |

Figure: `grid_lines.png`
