# Simulation Design, Scenarios, and Evaluation Protocol

本页说明模拟研究的真实基因型与参考面板构建、数据生成过程（DGP）、参数网格、八种因果配置、独立校准—测试窗口划分以及复现命令。理论模型与记号沿用 [Model and Method](method.md)。

## 1. Genotypes, cohorts, and reference LD

- **个体基因型与独立双队列**：取自 UK Biobank 22 号染色体欧洲血统个体（GRCh37），质控标准为 $$\mathrm{MAF} \ge 0.05$$、缺失率 $$\le 2\%$$、Hardy–Weinberg 平衡检验 $$p \ge 10^{-6}$$，并剔除 A/T 与 C/G 回文 SNP。选取 12 个长度约 1 Mb 的顺式窗口（每窗口 $$M \in [546, 1{,}537]$$，中位数 $$938$$）。性状 1 随机抽取 $$N_1 = 20{,}000$$ 人，性状 2 抽取完全不重叠的另外 $$N_2 = 1{,}000$$ 人。
- **外部参考 LD 面板**：采用 1000 Genomes Project Phase 3 EUR 队列（$$n_{\mathrm{ref}} = 489$$），按 rsID、物理坐标及效应等位基因精确对齐（12 个窗口合计 11,850 个 SNP 中对齐 11,794 个，占 99.5%；等位基因链向互换时同步翻转相关矩阵符号）。主模拟中所有方法统一输入该外部参考相关矩阵 $$\widehat{\mathbf{R}}$$；样本内匹配 LD 对照实验（[Figure 4](fig4.md)）则分别输入各性状自身队列的真实样本内矩阵 $$\mathbf{R}_1, \mathbf{R}_2$$。

## 2. Data-generating process and parameter grid

对每个基因组窗口 $$w$$、因果场景 $$\mathcal{M}$$、因果间 LD 档位 $$\tau$$ 与蒙特卡洛重复 $$\iota \in \{1, \dots, R\}$$，按如下三步生成汇总统计量：

1. **抽取因果集合**：在两队列内均满足 $$\mathrm{MAF} \ge 0.05$$ 的多态性位点中，依场景 $$\mathcal{M}$$ 抽取 $$(\mathcal{C}_1, \mathcal{C}_2)$$，使指定因果变异对在性状 1 真实样本内 LD 上满足档位 $$\tau$$ 的相关性约束。
2. **生成真实联合效应**：当 $$K_k \ge 1$$ 时，主设计按 [Method](method.md) 式 (1.3)–(1.4) 生成 $$\beta_{kj} = \xi_k s_{kj}\sqrt{h^2_k / K_k}$$（$$s_{kj} \overset{\mathrm{iid}}{\sim} \mathrm{Unif}\{-1, +1\}$$），使经验遗传方差 $$\boldsymbol\beta_k^\top \mathbf{R}_k \boldsymbol\beta_k = h^2_k$$ 精确成立；正态效应灵敏度分支则生成不作重标定的 $$\beta_{kj} \sim \mathcal{N}(0, h^2_k / K_k)$$。当 $$K_k = 0$$（`gwas_only` 性状 2 或 `null_both`）时令 $$\boldsymbol\beta_k = \mathbf{0}$$、$$\boldsymbol\eta_k = \mathbf{X}_k \boldsymbol\beta_k = \mathbf{0}$$。
3. **生成表型与单变量 OLS 汇总统计量**：生成 $$\mathbf{y}_k = \boldsymbol\eta_k + \boldsymbol\varepsilon_k$$，其中 $$\boldsymbol\varepsilon_k = \sqrt{v_{ek}}\,\mathbf{u}_k$$、$$\mathbf{u}_k \overset{\mathrm{iid}}{\sim} \mathcal{N}(\mathbf{0}, \mathbf{I}_{N_k})$$。在固定效应主设计中，残差方差取 $$v_{ek} = 1 - \operatorname{Var}_{\mathrm{emp}}(\boldsymbol\eta_k)$$（故 $$K_k \ge 1$$ 时 $$v_{ek} = 1 - h^2_k$$，而当 $$K_k = 0$$ 时无论参数名义 $$h^2_k$$ 为何均精确有 $$v_{ek} = 1$$）；在正态效应灵敏度分支中，$$K_k \ge 1$$ 时取 $$v_{ek} = 1 - h^2_k$$（$$K_k = 0$$ 时亦为 $$1$$）。随后逐 SNP 运行单变量 OLS 得 $$(\hat b_{kj}, \hat s_{kj}, p_{kj})$$ 及样本内 MAF，连同 $$N_k$$ 与指定 LD 矩阵输入各方法。

**公共随机数（Common Random Numbers, CRN）配对设计.** 标准正态噪声抽样 $$\mathbf{u}_k$$ 由 $$(w, \mathcal{M}, \iota)$$ 唯一固定，因果集合与效应符号由 $$(w, \mathcal{M}, \tau, \iota)$$ 唯一固定，均与局部遗传度 $$h^2_k$$ 及分析所用 LD 来源无关（残差仅随 $$\sqrt{v_{ek}}$$ 缩放）。该设计使跨参数设定与跨 LD 面板的比较在逐位点配对样本上进行，正相关性有效降低了差值估计量的蒙特卡洛方差。

| 参数 | 统计含义 | 第一阶段（默认设定） | 第二阶段（单因素网格扫描） | 设定依据 |
|---|---|---|---|---|
| $$N_1, N_2$$ | 性状 1 / 性状 2 样本量 | 20,000 / 1,000 | 同左 | 对标大型 GWAS 与千人级 bulk eQTL 队列 |
| $$h^2_1$$ | 性状 1 局部遗传度 | 0.003（$$\lambda_1 = 7.76 / 5.48$$） | 0.001, 0.003, 0.01, 0.02 | 默认单因果全基因组显著率 99%；上端 $$\lambda_1 = 20.20$$ |
| $$h^2_2$$ | 性状 2 局部遗传度 | 0.015（$$\lambda_2 = 3.90 / 2.75$$） | 0.0075, 0.015, 0.03, 0.05 | 默认双因果处于检出边缘；上端 5% 与 SharePro (Zhang et al. 2024) 一致 |
| $$\tau$$ | 因果变异间 $$\lvert r\rvert$$ 档位 | `low`：$$[0, 0.2)$$ | `low`, `mid` $$[0.32, 0.63)$$, `high` $$[0.63, 0.9)$$ | 切点对齐 $$\sqrt{\rho_c} \approx 0.32$$ 与 $$\sqrt{\rho_b} \approx 0.63$$ |
| $$\sqrt{W}$$ | ABF 先验标准差 | $$0.15\,\sigma_y$$ | 同左 | `coloc` 连续性状默认先验尺度 |
| $$(p_1, p_2, p_{12})$$ | 每 SNP 先验概率 | $$(10^{-4}, 10^{-4}, 10^{-5})$$ | 同左 | `coloc` 默认先验；另设 $$(M/M_B)$$ 区域等价块先验对照 |
| 重复次数 $$R$$ | 每（窗口 $$\times$$ 场景）重复数 | 100（共 9,600 位点） | 50（每组参数设定 4,800 位点） | 第一阶段合并阳性率标准误约 0.9 个百分点 |

*Note.* 第一阶段默认与各扫描设定的非中心参数基准 $$\lambda_k = \sqrt{N_k q_k/(1 - q_k)}$$ 分别对应单因果（$$K_k = 1$$）与双因果（$$K_k = 2$$）。

## 3. Causal configurations and window-split evaluation protocol

记 $$(K_C, K_{S1}, K_{S2})$$ 为共享因果个数及两性状各自特有因果个数：

| 场景名称 | $$(K_C, K_{S1}, K_{S2})$$ | 性状 1 因果 | 性状 2 因果 | 真值类别 | 评价角色 |
|---|---|---|---|---|---|
| `single_shared` | $$(1, 0, 0)$$ | $$A$$ | $$A$$ | 共定位阳性（$$H_4^{\mathrm{R}}$$） | 检出率 |
| `single_double_shared` | $$(1, 0, 1)$$ | $$A$$ | $$A, B$$ | 共定位阳性（$$H_4^{\mathrm{R}}$$） | 检出率 |
| `double_shared2` | $$(2, 0, 0)$$ | $$A, B$$ | $$A, B$$ | 共定位阳性（$$H_4^{\mathrm{R}}$$） | 检出率 |
| `double_shared1` | $$(1, 1, 1)$$ | $$A, B$$ | $$A, C$$ | 共定位阳性（$$H_4^{\mathrm{R}}$$） | 检出率 |
| `single_distinct` | $$(0, 1, 1)$$ | $$A$$ | $$B$$ | 双性状独立信号（$$H_3^{\mathrm{R}}$$） | 校准临界值 + 测试误报率 |
| `double_distinct` | $$(0, 2, 2)$$ | $$A, B$$ | $$C, D$$ | 双性状独立信号（$$H_3^{\mathrm{R}}$$） | 校准临界值 + 测试误报率 |
| `gwas_only` | $$(0, 1, 0)$$ | $$A$$ | $$\varnothing$$ | 单性状信号（$$H_1^{\mathrm{R}}$$） | 额外阴性假阳性率 |
| `null_both` | $$(0, 0, 0)$$ | $$\varnothing$$ | $$\varnothing$$ | 全零假设（$$H_0^{\mathrm{R}}$$） | 额外阴性假阳性率 |

档位 $$\tau$$ 作用于同性状内的双因果对及 `single_distinct` 中的跨性状对 $$A$$–$$B$$；其余跨性状不同因果对在 `low` 档与多数场景中限制 $$\lvert r\rvert < 0.2$$，唯 `double_shared1` 的跨性状非共享对 $$B$$–$$C$$ 在 `mid` 与 `high` 档中允许低于该档位上界（跨性状中等 LD 情景另于 [Supplementary](supplement.md) §1.4 单独构造）。

**独立校准—测试窗口分割与统计推断.** 12 个顺式窗口固定划分为 6 个校准窗口与 6 个测试窗口：

- **按 5% 假阳性率校准临界值**：对每个方法 $$m$$，在 6 个校准窗口合并的 `single_distinct` 与 `double_distinct` 阴性位点（第一阶段 $$n = 1{,}200$$；无结果记为 $$-\infty$$，排序严格位于一切有限得分之下）上降序排列，取第 $$\lfloor 0.05\,n\rfloor + 1$$ 位值作为固定临界值 $$t_m$$。分类判决采用严格大于（$$> t_m$$），ROC 曲线则按 $$\ge$$ 扫描并完整保留并列点；第二阶段单因素扫描沿用第一阶段固定临界值 $$t_m$$ 与 $$\alpha_{\mathrm{lead}}$$ 以检验临界值迁移稳健性。
- **测试窗口聚类推断**：所有主基准检出率、误报率与 AUC（附 Wilson 95% 置信区间及自举 AUC 95% 区间）均在 6 个独立测试窗口上报告（唯一例外为 [Figure 3](fig3.md) 剖析全区域 `coloc.abf` 失效机理的分层补充分析，合并 12 个窗口以增大各层样本量）。针对同一窗口内位点共享 LD 结构的聚类相关性，配对方法比较以 **6 个独立测试窗口差值均值的自由度 5 的 $$t_5$$ 95% 置信区间**及同时对校准窗口与测试窗口进行整群重抽样的**校准—测试窗口自举（cluster bootstrap）95% 区间**作为主推断，并将合并位点层面的精确 McNemar 检验及 Wald 区间作为描述性补充（[Supplementary](supplement.md) §2.1）。

## 4. Reproducibility and core scripts

| 脚本文件 | 功能定位 |
|---|---|
| `R/methods.R`, `R/simulate_bulk.R` | 全部比较方法实现及基于真实基因型的表型与汇总统计量数据生成函数 |
| `scripts/test_synthetic.R` | 在合成基因型上执行 22 项自动化单元测试（LD 对齐、重标定、零分布与 CRN 验证） |
| `scripts/prepare_ref_ld.R` | 对齐 1000 Genomes EUR 参考面板并生成各窗口相关矩阵 |
| `scripts/run.R`, `scripts/evaluate.R` | 第一阶段（`main`）与第二阶段（`grid`）模拟主程序与窗口级评估脚本 |
| `scripts/prior_rule.R`, `scripts/prior_rule_eval.py` | 区块每 SNP 先验与区域等价先验对比（匹配误报率及固定三分类规则） |
| `scripts/gain_decomposition.py`, `scripts/ld_compare.py` | 全概率三分支精确分解（[Figure 1](fig1.md)）与样本内/参考 LD 逐位点配对比较（[Figure 4](fig4.md)） |
| `scripts/block_ablation.R`, `scripts/block_sdy.R` | 分块规则消融与表型尺度 $$\sigma_{yk}$$ 估计方式敏感性分析（[Supplementary](supplement.md)） |
| `scripts/extract_plot_tables.R`, `scripts/plot_main_figures.R` | 从原始输出提取校验汇总表并绘制全套统一格式主图与补充图 |

在仓库根目录（R $$\ge 4.4$$）运行复现流程：

```bash
# 1. 合成基因型 22 项自检（不依赖本地 UK Biobank 数据）
Rscript scripts/test_synthetic.R

# 2. 对齐参考 LD 并运行两阶段主模拟（需本地基因型目录）
Rscript scripts/prepare_ref_ld.R
Rscript scripts/run.R --stage main --cores 4 && Rscript scripts/evaluate.R --stage main
Rscript scripts/run.R --stage grid --cores 3 && Rscript scripts/evaluate.R --stage grid

# 3. 样本内匹配 LD、块先验、分块规则消融与表型尺度敏感性
Rscript scripts/run.R --stage main --ld insample --cores 4
Rscript scripts/prior_rule.R --cores 4
Rscript scripts/block_ablation.R --cores 4
Rscript scripts/block_sdy.R --cores 4

# 4. 由已保存汇总表重新绘制全部主图与补充图
Rscript scripts/plot_main_figures.R
```
