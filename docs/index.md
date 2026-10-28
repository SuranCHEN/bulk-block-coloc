# block-coloc (bulk)

**Goal.** 评估 block-coloc 在“一个区域内有多个因果变异”时的共定位表现。block-coloc 只在条件分析发现多个独立信号时，按 LD 把区域切块并在块内做 coloc.abf；否则做整区域 coloc.abf。我们在真实 UK Biobank chr22 基因型上模拟两个连续性状（大样本 GWAS，$$N = 20{,}000$$；小样本 bulk eQTL，$$N = 1{,}000$$），所有方法使用 1000 Genomes EUR 参考 LD，在相同误报率下与整区域 coloc.abf、coloc.susie、coloc.signals 比较。

## Documentation

- [Model](01_bulk_model.md)：个体水平线性模型，summary statistics 的分布 $$\mathbf{z} \sim \mathcal{N}(\sqrt{N}\mathbf{R}\boldsymbol\beta, \mathbf{R})$$，coloc.abf 的假设与后验，多因果变异下的现有方法，条件 $$z$$，参考 LD
- [Method](02_block_coloc_pseudocode.md)：block-coloc 算法，比较的版本，与 coloc.signals / PWCoCo 的关系，对照方法
- [Simulation](03_simulation.md)：与已有模拟研究的对照，数据生成过程，参数，场景，评价方式
- [Supplementary Note](09_parameters.md)：模拟参数的推导。假设 A1–A5；Proposition 1–9：因果变异的信号强度、效应分布、因果间 LD、有无信号时的 lead、ABF 收缩、coloc 先验质量、参考 LD 误差、评价精度

## Results

- [Stage 1: main setting](04_results.md)：功效、误报率、机制分析、敏感性分析、理论与观察的对照
- [Stage 2: one factor at a time](04_results.md#6-stage-2-one-factor-at-a-time)：$$h^2_1$$、$$h^2_2$$、因果间 LD 的单因素变化
- [Previous results](07_results_v1.md)：上一版设定（每个因果变异固定解释比例、样本内 LD）
- [Early bulk results](05_early_bulk_results.md)：10 月中按 8 月旧设定重建的结果

## Code

- [Code](08_code.md)：文件结构、数据流和复现命令

## Main result

**Stage 1（默认设定）.** 误报率匹配到 5% 时，block-coloc 的检出率为 76.8%，比整区域 coloc.abf（73.2%）高 3.5 个百分点（$$p < 10^{-15}$$），与 coloc.signals (cond)（76.2%）相当。增益有限的原因是 eQTL 的两个因果变异平分 1.5% 的遗传度时，每个信号的期望 $$\lvert z\rvert$$ 只有 2.75，eQTL 只有约 15% 的位点找到两个条件独立的 lead。

**Stage 2（单因素变化）.** block-coloc 的优势随 eQTL 强度增大：比整区域 coloc.abf 高 4.0（$$h^2_2$$ = 1.5%）、7.1（3%）、13.7（5%）个百分点；eQTL 找到两个 lead 的比例（4%、16%、56%、89%）与运行前的理论预测一致。coloc.signals (cond) 在信号强时检出率更高，但沿用第一段判定线时误报率升到 6%–9%。因果变异间 $$\lvert r\rvert \ge 0.32$$ 时所有方法都无法区分不同的因果变异。详见 [Results](04_results.md)。

只有基因型和 LD 是真实数据；因果变异、效应和表型均为模拟，结果不对应任何真实基因或疾病。
