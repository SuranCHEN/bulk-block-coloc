# bulk 版 block-coloc

在 GWAS × bulk eQTL 共定位中，一个区域里可能有多个独立信号，这时经典 coloc.abf 的“每个性状至多一个因果变异”假设不成立。
本项目评估 **block-coloc**：只在确有多个独立信号时，用 LD 把区域切块、每块跑 coloc.abf；其余情况直接跑整区域 coloc.abf。
我们在真实 UKBB chr22 基因型上做模拟，与整区域 coloc.abf、coloc.susie、coloc.signals 比较，看在相同误报率下，block-coloc 在哪些条件下检出率更高或更低、定位更准或更差。

## 主要结论（第一段：主设置）

![测试窗口上各方法在每个场景中被判为共定位的比例（误报率匹配到 5%）](matched_fpr_by_scenario.png)

左边 3 组是共定位场景（越高越好），右边 2 组是不共定位场景（越低越好），虚线是 5% 误报率。

1. **检出率**：误报率都控制在 5% 时，block-coloc（条件 lead 选择）检出 93.3% 的共定位位点，高于整区域 coloc.abf（91.2%，配对差 +2.2 个百分点，95% CI +1.1 至 +3.2）和 coloc.susie + fallback（90.1%）。提升主要来自 double_shared2（94.2% 对 89.5%）。但它低于 coloc.signals (cond)（95.1%），与 block-coloc（交接包原版，94.2%）相近。
2. **是否切块**：条件 lead 选择几乎不会把单信号位点错误地切块（single_shared 1.5%，交接包原版为 84.8%）；但 eQTL 信号较弱（N = 1,000，|z| ≈ 3.9）时常漏掉第二个信号，single_double_shared 只有 26.0% 切块。
3. **最难的情况**：single_double_shared 中，当 eQTL 特有的信号比共享信号更显著时，所有方法的检出率都明显下降（coloc.abf 73.0%，block-coloc 条件 lead 选择 76.2%，coloc.susie + fallback 54.7%）。

这是一个设定下的结果，测试窗口只有 6 个；详细数字、置信区间和注意事项见 [结果](04_results.md)。第二段（变强度、变 LD）待定。

## 目录

1. [bulk 模型](01_bulk_model.md)：数据、线性模型、统计量、`z ~ N(Rλ, R)`、场景与 LD 分档
2. [block-coloc 流程与伪代码](02_block_coloc_pseudocode.md)：block-coloc 三个版本（条件 lead 选择、边际 lead 选择、交接包原版）各做了什么、改了什么、为什么改
3. [代码说明与复现](03_code.md)：文件结构、数据流、复现命令、自检
4. [结果](04_results.md)：第一段（主设置）已完成；第二段（变强度、变 LD）待定
5. [早期 bulk 结果](05_early_bulk_results.md)：10 月中按 8 月旧设定重建的结果（设定不同，不能直接和第 4 页比较）

## 一句话设定

| 项目 | 设定 |
|---|---|
| 基因型 | 真实 UKBB chr22，12 个约 1 Mb 窗口（6 个校准 + 6 个测试） |
| 样本 | GWAS 20,000 人，eQTL 1,000 人，互不重叠 |
| 模型 | 两个性状都是线性模型；每个因果 SNP 解释方差 q |
| 强度 | GWAS q 0.1%–2%（主 0.3%）；eQTL q 0.75%–5%（主 1.5%） |
| 场景 | single_shared、single_double_shared、double_shared2（共定位）；single_distinct、double_distinct（不共定位） |
| 因果间 LD | 主设置 \|r\| < 0.2；敏感性 0.32–0.63、0.63–0.9 |
| 方法 | block-coloc（条件 lead 选择，本文方法）、block-coloc（边际 lead 选择）、block-coloc（交接包原版）；对照：整区域 coloc.abf、coloc.susie（raw / fallback）、coloc.signals（mask / cond） |
| 评价 | 校准窗口上固定误报率 5%，测试窗口上报告检出率等 |

只有基因型和 LD 是真实的；因果 SNP、效应、表型都是模拟的，结果不对应任何真实基因或疾病。

## 主要参考文献

- Giambartolomei C, et al. (2014). *PLoS Genet* 10(5): e1004383.（coloc.abf）
- Wallace C (2020). *PLoS Genet* 16(4): e1008720.（coloc.signals）
- Wallace C (2021). *PLoS Genet* 17(9): e1009440.（coloc.susie；本项目 single_distinct、single_shared、single_double_shared、double_shared2 四个场景分别对应其 Fig 1 的第 1–4 行）
- Wang G, et al. (2020). *JRSS B* 82(5): 1273–1300.（SuSiE）
