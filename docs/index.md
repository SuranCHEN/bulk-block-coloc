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
