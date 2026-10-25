# bulk 版 block-coloc

在 GWAS × bulk eQTL 共定位中，一个区域里可能有多个独立信号，这时经典 coloc.abf 的“每个性状至多一个因果变异”假设不成立。
本项目评估 **block-coloc**：只在确有多个独立信号时，用 LD 把区域切块、每块跑 coloc.abf；其余情况直接跑整区域 coloc.abf。
我们在真实 UKBB chr22 基因型上做模拟，与整区域 coloc.abf、coloc.susie、coloc.signals 比较，看在相同误报率下，block-coloc 在哪些条件下检出率更高或更低、定位更准或更差。

## 目录

1. [bulk 模型](01_bulk_model.md)：数据、线性模型、统计量、`z ~ N(Rλ, R)`、场景与 LD 分档
2. [block-coloc 流程与伪代码](02_block_coloc_pseudocode.md)：block-coloc 三个版本（条件 lead 选择、边际 lead 选择、交接包原版）各做了什么、改了什么、为什么改
3. [代码说明与复现](03_code.md)：文件结构、数据流、复现命令、自检
4. [结果](04_results.md)：第一段（主设置）已完成；第二段（变强度、变 LD）待定

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
