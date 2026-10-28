# Simulation

本页给出模拟的数据生成过程、参数、场景和评价方式。设计参照已有共定位方法学文章的模拟，尤其是 SharePro（Zhang et al. 2024）和 coloc.susie（Wallace 2021），差异和理由在 §1 列出。模型记号见 [Model](01_bulk_model.md)。代码：`scripts/`。

## 1. Design and relation to previous simulation studies

| 研究 | 基因型 / LD | 区域 | 样本量 | 因果变异与效应 | 重复 | 评价 |
|---|---|---|---|---|---|---|
| Giambartolomei et al. (2014), coloc.abf | 模拟 | cis 窗口 | eQTL 约 1,000；biomarker 2,000–10,000 | 1 个共享或不同的因果变异；按方差解释设定（eQTL 10%，biomarker 约 2%） | 每场景 1,000 | PP3 / PP4 > 0.9 的比例 |
| Wallace (2020), coloc.signals | 1000 Genomes EUR / YRI 单倍型 | 250–750 个 SNP | 100–10,000 | 1–2 个因果（MAF > 5%）；连续性状，残差 sd 1 | — | 各假设的平均后验；条件 vs 屏蔽 |
| Wallace (2021), coloc.susie | 1000 Genomes EUR 单倍型 | 1,000 个 SNP | 病例 / 对照各 10,000（simGWAS） | 每个性状 1–2 个因果（MAF > 5%），log OR 取自均值为 0 的正态分布；共享 0、1 或 2 个 | 10,000 | 按真实因果变异标注的平均后验 |
| Hukku et al. (2021), fastENLOC | GTEx 400 人，样本内 LD | 1,500 个 cis-SNP | 400 | 每区域约 3 个 eQTL，效应 $$\sim\mathcal{N}(0,1)$$ | — | 5% FDR 下的功效；LD 失配另做模拟 |
| Jiang et al. (2026), traceCB | UK Biobank 与东亚队列；1000 Genomes LD | 2,000 个 SNP（chr22） | 50–5,000 | 因果比例 $$p_{\mathrm{causal}} \in \{0.005, 0.01, 0.02\}$$，效应取自与 $$h^2$$ 挂钩的双变量正态 | 每设定 100 | 名义 0.05 下的 I 类错误与功效 |
| Zhang et al. (2024), SharePro | UK Biobank，不重叠的两组欧洲人 | 5 个 1 Mb 区域 | 25,000 / 1,000 | $$K_C$$ 个共享 + 各 $$K_S$$ 个特有因果；每个因果变异的遗传度固定为 0.01（性状 1）、0.05（性状 2） | 每设定 50 | PP > 0.8 时的 TPR / FPR |
| Zhang et al. (2025), jaxQTL（eQTL 定位，非共定位） | 单个 cis 变异，MAF 0.2；OneK1K 的文库大小 | 1 个 SNP | 1,000（默认） | $$h^2_{cis} = 0.05$$，$$\beta_g \sim \mathcal{N}(0, h^2_{cis})$$ | 每设定 500 | 名义 0.05 下的 I 类错误与功效 |
| **本项目** | UK Biobank chr22，不重叠的两组欧洲人；方法用 1000 Genomes EUR 参考 LD | 12 个 1 Mb 窗口（6 校准 + 6 测试） | 20,000 / 1,000 | 每性状 0–2 个因果，8 种共享结构，3 档因果间 LD；区域遗传度按因果个数平分 | 100（第一段）/ 50（第二段） | 误报率匹配到 5% 时的检出率；误报率；AUC |

设计上的几点选择：

- **样本与区域**：与 SharePro 相同，用真实 UK Biobank 基因型、两组不重叠的人、1 Mb 区域，“大样本性状 + 小样本 eQTL”。窗口数从 5 增加到 12，并分成校准和测试两组，用于在不重叠的窗口上定判定线和报告结果（§6）。
- **效应参数化**：本项目固定区域遗传度 $$h^2_k$$ 和因果个数 $$K_k$$，每个因果变异分到 $$h^2_k/K_k$$（沿用组内 8 月 bulk）。这与 traceCB 的 $$(h^2, p_{\mathrm{causal}})$$ 等价：$$K_k = M p_k$$，区域平均的 per-SNP 遗传度 $$\omega_k = h^2_k/M$$（[Model](01_bulk_model.md) 式 2.3）。SharePro 则固定每个因果变异的遗传度，因此多信号时总遗传度随因果个数增加；本项目多信号时总遗传度不变、每个信号更弱。
- **LD**：上述研究的方法大多使用样本内 LD（Hukku et al. 2021 另做了跨人群的失配模拟）。本项目让所有方法使用 1000 Genomes EUR 参考 LD，以接近只有 summary statistics 的实际情况（[Model](01_bulk_model.md) §7）。
- **评价**：上述研究多用固定后验阈值（0.8 或 0.9）。各方法 PP4 的刻度不同，固定阈值下误报率不同，比较检出率不公平；本项目先在校准窗口上把每个方法的误报率匹配到 5%，再在测试窗口上比较检出率（§6）。

## 2. Genotypes and reference LD

- **基因型**：本地 UK Biobank chr22（GRCh37），QC 为 MAF ≥ 0.05、缺失率 ≤ 2%、HWE $$p \ge 10^{-6}$$，并去掉 A/T、C/G 回文 SNP。12 个以已知位置为中心的 1 Mb 窗口，每窗口 546–1,537 个 SNP。
- **样本**：性状 1 用 20,000 人，性状 2 用另外 1,000 人，两组不重叠（已按样本编号核对）。
- **参考 LD**：1000 Genomes Phase 3 EUR，489 人。每个窗口只保留能与参考人群按 rsID、位置和等位基因对齐的 SNP（99.5%），模拟和分析都在这些 SNP 上进行。$$\widehat{\mathbf{R}}$$ 按计数等位基因对齐符号。

## 3. Data-generating process

对窗口 $$w$$、场景 $$s$$、LD 档位 $$\tau$$、重复 $$t$$：

1. **因果集合.** 按场景 $$s$$ 的共享结构（§4）从两组样本中 MAF 均 ≥ 0.05 的 SNP 里随机抽取 $$\mathcal{C}_1, \mathcal{C}_2$$，满足 §5 的 LD 约束（约束用 20,000 人的样本内 LD，即真实 LD）。
2. **效应.** 对 $$j \in \mathcal{C}_k$$，

   $$
   \beta_{kj} = c_k\, s_{kj} \sqrt{h^2_k / K_k}, \qquad s_{kj} \overset{\text{iid}}{\sim} \mathrm{Unif}\{-1, +1\},
   $$

   其中 $$c_k$$ 使 $$\operatorname{Var}(\mathbf{X}_k\boldsymbol\beta_k) = h^2_k$$ 在样本内严格成立（推导见 [Model](01_bulk_model.md) 式 2.2）。当 $$K_k = 2$$、两个因果变异 LD 为 $$r$$ 时，$$c_k = (1 + s_{k1}s_{k2} r)^{-1/2}$$。
3. **表型.** $$\mathbf{y}_k = \mathbf{X}_k\boldsymbol\beta_k + \boldsymbol\varepsilon_k$$，$$\boldsymbol\varepsilon_k \sim \mathcal{N}(\mathbf{0}, (1 - h^2_k)\mathbf{I})$$，$$\mathbf{X}_k$$ 为样本内标准化的真实基因型。
4. **Summary statistics.** 逐 SNP 的 OLS（[Model](01_bulk_model.md) 式 3.1、4.2），得到 $$(\hat b_{kj}, \hat s_{kj}, p_{kj})$$ 和样本内 MAF。
5. **分析.** 各方法输入 $$(\hat b, \hat s, \mathrm{MAF}, N)$$ 和 $$\widehat{\mathbf{R}}$$，输出 PP4（无结果记为 NA）。

**Common random numbers.** 噪声向量 $$\boldsymbol\varepsilon_k$$ 由 $$(w, s, t)$$ 决定，因果集合和效应方向由 $$(w, s, t, \tau)$$ 决定，与 $$h^2$$ 无关。因此改变 $$h^2$$ 时只有效应大小变化，不同强度之间的比较是配对的。

## 4. Parameters

| 记号 | 含义 | 默认值 | 第二段取值 | 依据 |
|---|---|---|---|---|
| $$N_1$$ / $$N_2$$ | 样本量 | 20,000 / 1,000 | — | 与 SharePro（25,000 / 1,000）同一量级；$$N_2$$ 与 OneK1K（982 人）相当 |
| $$h^2_1$$ | 性状 1 的区域遗传度 | 0.003 | 0.001, 0.003, 0.01, 0.02 | 组内 8 月 bulk；单因果时全基因组显著的概率 0.99（S2） |
| $$h^2_2$$ | 性状 2 的区域遗传度 | 0.015 | 0.0075, 0.015, 0.03, 0.05 | 组内 8 月 bulk；上端与 jaxQTL 的 $$h^2_{cis} = 0.05$$ 相同 |
| $$K_1, K_2$$ | 因果变异个数 | 0–2，由场景决定 | — | Wallace (2021)；SharePro 的 $$(K_C, K_S)$$ |
| 效应 | 每个因果 $$\pm\sqrt{h^2_k/K_k}$$，符号独立 | — | — | 信号强度只由 $$(N, h^2, K)$$ 决定（S3）；traceCB 用正态效应 |
| $$\sigma^2_k$$ | 残差方差 | $$1 - h^2_k$$ | — | $$\operatorname{Var}(y) = 1$$，同 traceCB |
| $$\tau$$ | 同一性状两个因果之间的 $$\lvert r\rvert$$ | low：[0, 0.2) | mid：[0.32, 0.63)，high：[0.63, 0.9) | 切点为 $$\sqrt{\rho_c}$$、$$\sqrt{\rho_b}$$（S4） |
| $$\mathrm{MAF}_{\min}$$ | 因果变异最小 MAF | 0.05 | — | Wallace (2020, 2021) |
| $$n_{\mathrm{ref}}$$ | 参考 LD 人数 | 489 | — | 1000 Genomes EUR（S8） |
| $$\sqrt{W}$$ | ABF 先验标准差 | 0.15 | — | coloc 默认（Giambartolomei et al. 2014；S6） |
| $$(p_1, p_2, p_{12})$$ | coloc 先验 | $$(10^{-4}, 10^{-4}, 10^{-5})$$ | — | coloc 默认（S7） |
| $$\alpha_{\mathrm{lead}}$$ | lead 阈值 | 在校准窗口上从 $$\{10^{-3}, 10^{-4}, 10^{-5}\}$$ 中选 | 沿用第一段 | [Method](02_block_coloc_pseudocode.md) §5；无信号时的假 lead 率见 S5 |
| 重复 | 每（窗口、场景） | 100 | 50 | 单场景检出率标准误 1.8 / 2.5 个百分点（S9） |

S2–S9 指 [Supplementary Note](09_parameters.md) 的各节。

**Signal strength.** 由 Proposition 1，因果变异 $$z$$ 的均值为 $$\lambda = \sqrt{N q/(1-q)}$$，$$q = h^2/K$$：

| | $$K = 1$$ | $$K = 2$$ | 第二段范围（$$K$$ = 1–2） |
|---|---|---|---|
| 性状 1，$$h^2_1 = 0.003$$ | 7.8 | 5.5 | 3.2–20.2 |
| 性状 2，$$h^2_2 = 0.015$$ | 3.9 | 2.75 | 1.9–7.3 |

按同一公式换算的已有研究：Giambartolomei et al. (2014) eQTL 10.5、biomarker 6.4–14.3；SharePro 性状 1 为 15.9、性状 2 为 7.3；jaxQTL 7.3；traceCB 每个因果 1.0–2.9（$$h^2$$ = 0.1–0.2 分给 10–40 个因果，$$N_2 = 400$$）。本项目默认的性状 2 介于 traceCB 与 SharePro / jaxQTL 之间，处于检测边缘（单个因果以 $$p < 10^{-3}$$ 检出的概率 0.73，双因果时每个 0.29）；第二段的上端与 SharePro / jaxQTL 相同。默认的性状 1 弱于 SharePro。

## 5. Scenarios

用 $$(K_C, K_{S1}, K_{S2})$$ 表示共享因果个数和两性状各自特有的因果个数（记法同 SharePro）。场景名沿用组内 8 月 bulk 的命名。

| 场景 | $$(K_C, K_{S1}, K_{S2})$$ | 性状 1 因果 | 性状 2 因果 | 真值 | 用途 |
|---|---|---|---|---|---|
| single_shared | (1, 0, 0) | A | A | 共定位 | 功效 |
| single_double_shared | (1, 0, 1) | A | A, B | 共定位 | 功效 |
| double_shared2 | (2, 0, 0) | A, B | A, B | 共定位 | 功效 |
| double_shared1 | (1, 1, 1) | A, B | A, C | 共定位 | 功效 |
| single_distinct | (0, 1, 1) | A | B | 不共定位 | 定判定线；误报率 |
| double_distinct | (0, 2, 2) | A, B | C, D | 不共定位 | 定判定线；误报率 |
| gwas_only | (0, 1, 0) | A | — | 不共定位 | 误报率（不参与定线） |
| null_both | (0, 0, 0) | — | — | 不共定位 | 误报率（不参与定线） |

single_distinct、single_shared、single_double_shared、double_shared2 分别对应 Wallace (2021) Fig 1 的四行；double_shared1、double_distinct 来自组内 8 月 bulk。

**因果间 LD.** 档位 $$\tau$$ 作用于同一性状的两个因果变异（single_double_shared、double_shared2、single_distinct 的 A–B；double_distinct 的 A–B 和 C–D；double_shared1 的 A–B 和 A–C），其余不同因果之间要求 $$\lvert r\rvert < 0.2$$。中、高档的 double_shared1 中，B 与 C 都和 A 强连锁，$$\lvert r_{BC}\rvert < 0.2$$ 在真实窗口中几乎无法满足（高档时每窗口只有 2–58 个 SNP 可作 A），故只要求 $$\lvert r_{BC}\rvert$$ 低于该档上限。档位切点 0.32 和 0.63 是方法阈值 $$\rho_c = 0.1$$、$$\rho_b = 0.4$$ 换算成 $$\lvert r\rvert$$ 的值，因此三档分别对应“两个 lead 各自成块”“第二个 lead 被 clump 去掉但不在第一个 lead 的邻域内”“两个因果落进同一块”。

## 6. Evaluation

**Calibration and test.** 12 个窗口固定分成 6 个校准窗口和 6 个测试窗口。所有依赖数据的选择（判定线、$$\alpha_{\mathrm{lead}}$$）只在校准窗口上做，结果只在测试窗口上报告。

**Matched type I error.** 对方法 $$m$$，取校准窗口上 single_distinct 和 double_distinct 的 PP4（NA 记为 −1），从大到小排序，判定线 $$t_m$$ 为第 $$\lfloor 0.05\, n \rfloor + 1$$ 个值，PP4 $$> t_m$$ 判为共定位。

**Metrics（测试窗口）.**

- 功效（检出率）：$$\Pr(\mathrm{PP4} > t_m \mid \text{共定位})$$，按场景和合并报告；
- 误报率：$$\Pr(\mathrm{PP4} > t_m \mid \text{不共定位})$$，分别报告 single_distinct、double_distinct 和两个额外阴性；
- AUC：不依赖判定线，1,000 次按真值分层的 bootstrap 给区间；
- 无结果率、走切块的比例、lead 个数是否等于真实因果个数、定位（选中的块是否包含共享因果）。

**Inference.** 比例用 Wilson 95% 区间。两个方法在同一批位点上的差异用 McNemar 精确检验：$$b_{10}, b_{01}$$ 为不一致对，在原假设下 $$b_{10} \sim \mathrm{Bin}(b_{10} + b_{01}, 1/2)$$；差值区间用配对比例的 Wald 区间。每个场景测试窗口 600 个位点，检出率 90% 时标准误约 1.2 个百分点。

## 7. Stages and replicates

| 阶段 | 设定 | 重复 | 位点数 |
|---|---|---|---|
| 第一段 | 默认参数，$$\tau$$ = low；8 个场景；$$\alpha_{\mathrm{lead}}$$ 三个候选及 8 月阈值都计算 | 100 | 12 × 8 × 100 = 9,600 |
| 第二段 | 从默认参数出发一次只改一个：$$h^2_1$$、$$h^2_2$$、$$\tau$$；判定线和 $$\alpha_{\mathrm{lead}}$$ 沿用第一段，不重新调 | 50 | 约 35,000 |

第二段沿用第一段的判定线，因此条件改变时误报率会偏离 5%，用来看判定线在条件变化下是否稳健。

## 8. Reproduction

```bash
Rscript scripts/test_synthetic.R          # 合成数据自检（22 项）
Rscript scripts/prepare_ref_ld.R          # 1000G EUR 参考 LD
Rscript scripts/run.R --stage main --cores 4
Rscript scripts/evaluate.R --stage main
Rscript scripts/run.R --stage grid --cores 3
Rscript scripts/evaluate.R --stage grid
```

R 4.5.1（macOS，arm64），R 包 coloc、susieR；第一段在 MacBook Air M4 上约 2.5 小时（3–4 个进程）。

## References

见 [Model](01_bulk_model.md#references)。
