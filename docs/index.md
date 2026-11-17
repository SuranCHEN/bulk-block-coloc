# block-coloc: Block-Partitioned Genetic Colocalisation under Reference LD

多因果变异顺式区域中的汇总统计量共定位方法：仅当逐步条件分析在至少一个性状中检出多个独立信号时，利用参考 LD 邻域构造联合候选块，并在块内对未调整的原始边际统计量运行 `coloc.abf`，兼顾多因果信号局部分离与外部参考 LD 误差下的统计稳健性。

- [Documentation](#documentation)
- [Simulation](#simulation)
- [Main Figures](#main-figures)
- [Supplementary](#supplementary)
- [Useful Links](#useful-links)

## Documentation

- [Model, Algorithm, and Analytical Properties](method.md)

## Simulation

- [Simulation Design, Scenarios, and Evaluation Protocol](simulation.md)

## Main Figures

- [Figure 1: Workflow and Exact Branch Decomposition](fig1.md)
- [Figure 2: Detection Power, Type I Error, and Signal-Strength Grids](fig2.md)
- [Figure 3: Multi-Signal Behaviour of Whole-Region coloc.abf](fig3.md)
- [Figure 4: Matched In-Sample LD vs 1000 Genomes Reference LD](fig4.md)
- [Figure 5: Fixed Three-Way Calling Rule and Block-Prior Sensitivity](fig5.md)

## Supplementary

- [Formal Proofs, Sensitivity Analyses, and Extended Ablations](supplement.md)

## Useful Links

- [GitHub Repository](https://github.com/SuranCHEN/bulk-block-coloc)
- [Saved Summary Tables (`results/`)](https://github.com/SuranCHEN/bulk-block-coloc/tree/main/results)

## Change log

- **2026-11-18**：依据统计学同侪评审统一本文档符号体系与假设–定理–证明–注记结构，补全窗口聚类推断与匹配/参考 LD 平衡对照表，并升级全套 `plotmath` 论文主图。
- **2026-11-17**：精简重组文档结构为方法、模拟与五主图加补充材料，完成参考 LD 失配分布及先验赔率公式校核。
- **2026-11-16**：补充全概率三分支分解、分块规则消融、表型尺度估计敏感性及四类阴性联合校准分析。
- **2026-11-15**：完成样本内匹配 LD 与 1000 Genomes EUR 参考 LD 的逐位点配对比较。
