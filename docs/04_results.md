# 4. 结果

> 第一段（主分析）已完成（2026-10-25）；第二段待用户看过第一段后再跑。
> 本页只摘要；完整表格见本仓库的 `results/stage1/RESULTS.md`。

## 4.1 评价规则（事先定好，不看结果调整）

- 每个方法的分数是 PP4；无结果计为“未判为共定位”，并单独报告无结果率。分母始终是全部位点。
- **固定误报率比较**：每个方法在 6 个校准窗口上定出一条 PP4 判定线，使不共定位场景（single_distinct、double_distinct）的误报率为 5%。
- block-coloc 各版本的 lead 阈值（1e-3 / 1e-4 / 1e-5）也只在校准窗口上选。
- 判定线冻结后，在另外 6 个测试窗口上报告：检出率、误报率、AUC、无结果率、走 block 分支的比例、是否正确判断要不要切块，以及与 `susie_fallback`、`whole_abf` 在同一批位点上的配对比较（McNemar 检验）。
- 区间：比例用 Wilson 95% 置信区间；AUC 用 1,000 次 bootstrap。

## 4.2 第一段：主设置

设定：GWAS q 0.3%（N 20,000）、eQTL q 1.5%（N 1,000）、因果间 |r| < 0.2；12 窗口 × 5 场景 × 100 次重复，校准、测试各 3,000 个位点。
lead 阈值在校准窗口上选出：block-coloc（条件 lead 选择）为 1e-4，block-coloc（边际 lead 选择）和 block-coloc（交接包原版）为 1e-5。

### 三行结论

1. **检出率**：误报率在校准窗口上匹配到 5%（测试窗口实际为 4.1%–4.8%，coloc.susie 为 2.8%）时，block-coloc（条件 lead 选择，本文方法）检出 93.3% 的共定位位点，高于整区域 coloc.abf 的 91.2%（同一批位点配对差 +2.2 个百分点，95% CI +1.1 至 +3.2）和 coloc.susie + coloc.abf fallback 的 90.1%（+3.2，+2.1 至 +4.4）。提升主要来自 double_shared2（94.2% 对 89.5%），single_shared 上没有损失（98.2% 对 97.8%）。但它没有超过 coloc.signals (cond)（95.1%），与 block-coloc（交接包原版，94.2%）也相近。
2. **是否切块**：条件 lead 选择基本消除了单信号位点被误切块：在 single_shared、single_distinct 上走 block 的比例为 1.5% 和 0.3%；边际 lead 选择为 3.7% 和 4.3%；交接包原版在 single_shared 上为 84.8%。但多信号常被漏判：single_double_shared 只有 26.0% 走 block，double_shared2 为 78.8%。原因是 eQTL 的第二个信号（N = 1,000，因果 |z| ≈ 3.9）经常过不了 lead 阈值。
3. **不共享信号更强时**：single_double_shared 中，当 eQTL 特有的因果变异比共享的因果变异更显著时，所有方法的检出率都明显下降（整区域 coloc.abf 73.0%、block-coloc（条件 lead 选择）76.2%、coloc.susie + fallback 54.7%、coloc.signals (cond) 82.1%；共享变异更显著时，除 coloc.susie 外各方法为 98.0%–100%）。eQTL 信号这么弱时，SuSiE 在 45%–67% 的位点上找不到 credible set，所以不加 fallback 的 coloc.susie 只检出 41.3%。

### 汇总（测试窗口）

| 方法 | 检出率 | 误报率 | AUC |
|---|---|---|---|
| **block-coloc（条件 lead 选择，本文方法）** | **93.3%** | 4.6% | 0.971 |
| block-coloc（边际 lead 选择） | 92.4% | 4.6% | 0.966 |
| block-coloc（交接包原版） | 94.2% | 4.8% | 0.977 |
| 整区域 coloc.abf | 91.2% | 4.5% | 0.964 |
| coloc.susie (raw) | 41.3% | 2.8% | 0.610 |
| coloc.susie + fallback | 90.1% | 4.6% | 0.946 |
| coloc.signals (mask) | 93.4% | 4.1% | 0.974 |
| coloc.signals (cond) | 95.1% | 4.8% | 0.983 |

| 场景 | block-coloc（条件 lead 选择） | coloc.abf | coloc.susie + fallback | coloc.signals (cond) |
|---|---|---|---|---|
| single_shared | 98.2% | 97.8% | 97.7% | 95.5% |
| single_double_shared | 87.7% | 86.2% | 75.8% | 90.7% |
| double_shared2 | 94.2% | 89.5% | 96.8% | 99.0% |
| single_distinct（误报） | 6.8% | 6.2% | 6.0% | 4.0% |
| double_distinct（误报） | 2.3% | 2.8% | 3.2% | 5.7% |

### 主图

![测试窗口上各方法在每个场景中被判为共定位的比例（误报率匹配到 5%）](matched_fpr_by_scenario.png)

左边 3 组是共定位场景（柱子越高越好），右边 2 组是不共定位场景（柱子越低越好），虚线是 5% 误报率。

### ROC 曲线

![测试窗口上的 ROC 曲线](roc_test.png)

### 读结果时要注意

- 这是**一个设定**下的结果；强度和 LD 改变后结论可能变化，见第二段。
- 测试窗口只有 6 个，各方法之间 2–3 个百分点的差距虽然配对检验显著，但幅度不大。
- 只有基因型和 LD 是真实的，结果不对应任何真实基因或疾病。

## 4.3 第二段：一次只变一个因素

（待补：GWAS 强度、eQTL 强度、因果间 LD 档位的影响，图 `grid_lines.png`）
