# Partitioned posterior（草稿）

本页把 block-coloc 从“对数据选出的块取最大 PP4”改写为一个分块模型下的区域后验。记号沿用 [Model](01_bulk_model.md) 与 [Method](02_block_coloc_pseudocode.md)；假设 (A1)–(A5) 见 [Supplementary Note](09_parameters.md)。本页是草稿。2026-11-13 已按 Codex 的独立审查（`collab/partition/REVIEW_P2_codex.md`）修订：命题的适用范围、先验的近似程度和文献先例。

## 1. Motivation

现行 block-coloc（[Method](02_block_coloc_pseudocode.md) Algorithm 1）有三处规则不来自模型：

1. **取最大.** 输出 $$\max_B \mathrm{PP4}(B)$$。它不是后验概率（Method §4 已说明），只能靠误报率匹配校准。
2. **no valid block 时退回整区域.** 两性状的 lead 邻域互不相交时，返回整区域 coloc.abf 的 PP4。但此时至少一个性状有两个信号，整区域的单因果假设恰好不成立；Results §6.1 中整区域 coloc.abf 在 $$h^2_2$$ 从 3% 增到 5% 时检出率反而下降，就是这一失效的表现。
3. **PP3 取自 PP4 最大的块.** 真实数据的 traceCB 规则需要 PP3，现在取的是“PP4 最大的那个块”的 PP3，没有区域层面的含义。

此外，块内用每 SNP 先验时，一个 $$M_B$$ 个 SNP 的块的先验比 $$H_4 : H_3 = p_{12}/\{(M_B-1)p_1p_2\}$$ 偏向共享（Proposition 7）；调成区域等价先验又矫枉过正（Results §8.2）。

下面的分块模型给出区域层面的 $$\mathrm{PP0}$$–$$\mathrm{PP4}$$；只有一块时退化为 coloc.abf；区域层面的先验与整区域 coloc.abf 在首项上相同（精确差别见 Remark 3）。数值比较（§6）显示它并不优于现行规则。

## 2. Region-level hypotheses

记 $$\mathcal{C}_1, \mathcal{C}_2$$ 为两性状在区域内的因果集合（个数不限）。区域层面的五个假设：

| | 定义 |
|---|---|
| $$H_0^{\mathrm{R}}$$ | $$\mathcal{C}_1 = \mathcal{C}_2 = \varnothing$$ |
| $$H_1^{\mathrm{R}}$$ | $$\mathcal{C}_1 \neq \varnothing,\ \mathcal{C}_2 = \varnothing$$ |
| $$H_2^{\mathrm{R}}$$ | $$\mathcal{C}_1 = \varnothing,\ \mathcal{C}_2 \neq \varnothing$$ |
| $$H_3^{\mathrm{R}}$$ | $$\mathcal{C}_1 \neq \varnothing,\ \mathcal{C}_2 \neq \varnothing,\ \mathcal{C}_1 \cap \mathcal{C}_2 = \varnothing$$ |
| $$H_4^{\mathrm{R}}$$ | $$\mathcal{C}_1 \cap \mathcal{C}_2 \neq \varnothing$$ |

当 $$\lvert\mathcal{C}_k\rvert \le 1$$ 时它们就是 coloc 的 $$H_0$$–$$H_4$$。模拟中的“共定位”真值即 $$H_4^{\mathrm{R}}$$，single_distinct、double_distinct 为 $$H_3^{\mathrm{R}}$$，gwas_only 为 $$H_1^{\mathrm{R}}$$，null_both 为 $$H_0^{\mathrm{R}}$$。block-coloc 检验的是 $$H_4^{\mathrm{R}}$$ 对其余四个假设。

## 3. Partitioned model

把区域的 SNP 分成互不相交的部分 $$\mathcal{P} = \{B_1, \dots, B_m\}$$，$$\bigcup_b B_b = \{1, \dots, M\}$$。假设：

- **(P1)** 划分 $$\mathcal{P}$$ 给定（不依赖数据；§5 讨论实际中的偏离）。
- **(P2)** 块间无 LD：$$\mathbf{R}_{B_a B_b} = \mathbf{0}$$，$$a \neq b$$。
- **(P3)** 每个性状在每个部分内至多一个因果变异。
- **(P4)** 先验在各部分之间独立；部分 $$b$$ 内采用 coloc 的每 SNP 先验 $$(p_1, p_2, p_{12})$$，即 $$H_0$$–$$H_4$$ 的先验质量与 $$1 : M_bp_1 : M_bp_2 : M_b(M_b-1)p_1p_2 : M_bp_{12}$$ 成比例。

**适用范围.** 以下命题针对一个分析用的工作模型：划分固定；$$\mathbf{z}_k \mid \boldsymbol\beta_k$$ 服从协方差已知为 $$\mathbf{R}$$ 的正态分布（式 4.3 本身是近似）；效应在标准化尺度上有独立的正态先验；$$\mathbf{R}_{BB}$$ 正定（参考 LD 只有 489 人，实际是秩亏的，严格的证明需改用广义逆并排除完全 LD）。模拟的生成模型（固定大小效应，A2）与此不同。实际的 coloc.abf 用每等位基因尺度的 $$\hat b$$、$$\hat s$$，并由 MAF、$$N$$ 和 $$\hat s$$ 估计 sdY（本项目未传入 sdY），只取块内 SNP 时估计值会变；因此对实际输出，下面的等式是把这些量当作已知的 plug-in 近似，不是精确后验。

**Proposition A (factorisation).** 在上述工作模型与 (P1)–(P4) 下，各部分的假设 $$(H^{(1)}, \dots, H^{(m)})$$ 的后验分解为

$$
\Pr\!\left(H^{(1)} = h_1, \dots, H^{(m)} = h_m \mid \mathbf{z}_1, \mathbf{z}_2\right) = \prod_{b=1}^m \mathrm{PP}_{h_b}(B_b),
$$

其中 $$\mathrm{PP}_h(B_b)$$ 是只在 $$B_b$$ 上运行 coloc.abf 得到的后验。

*Proof.* 由 [Model](01_bulk_model.md) 式 (4.3)，$$\mathbf{z}_k \mid \boldsymbol\beta_k \sim \mathcal{N}(\sqrt{N_k}\mathbf{R}\boldsymbol\beta_k, \mathbf{R})$$；两性状样本不重叠 (A3)，故给定效应时 $$\mathbf{z}_1 \perp \mathbf{z}_2$$。由 (P2)，$$\mathbf{R}$$ 分块对角，$$\mathbf{z}_{k, B_b}$$ 在各部分之间独立，且均值 $$\sqrt{N_k}\mathbf{R}_{B_bB_b}\boldsymbol\beta_{k,B_b}$$ 只依赖本部分的效应。于是似然分解为 $$\prod_b L_b$$。由 (P3)，部分 $$b$$ 内性状 $$k$$ 至多一个因果 $$j$$，均值为 $$\boldsymbol\mu = \sqrt{N_k}\beta_{kj}\mathbf{R}_{\cdot j}$$；因为 $$\mathbf{R}^{-1}\mathbf{R}_{\cdot j} = \mathbf{e}_j$$，

$$
\log\frac{L_b(\beta_{kj})}{L_b(0)} = \boldsymbol\mu^{\top}\mathbf{R}^{-1}\mathbf{z} - \tfrac12\boldsymbol\mu^{\top}\mathbf{R}^{-1}\boldsymbol\mu = \sqrt{N_k}\beta_{kj}z_{kj} - \tfrac12 N_k\beta_{kj}^2,
$$

只依赖 $$z_{kj}$$；对 $$\beta_{kj} \sim \mathcal{N}(0, W)$$ 积分即得 Wakefield 的 ABF（式 5.1）。因此部分 $$b$$ 内的 Bayes factor 正是 coloc.abf 在 $$B_b$$ 上所用的 Bayes factor。由 (P4) 先验也分解，后验随之分解，每个因子在 $$B_b$$ 内归一化后即 coloc.abf 的后验。$$\square$$

**Proposition B (region posterior).** 在 Proposition A 的条件下，记

$$
A = \prod_b \{1 - \mathrm{PP4}(B_b)\},\quad
C_1 = \prod_b \{\mathrm{PP0}(B_b) + \mathrm{PP2}(B_b)\},\quad
C_2 = \prod_b \{\mathrm{PP0}(B_b) + \mathrm{PP1}(B_b)\},\quad
C_0 = \prod_b \mathrm{PP0}(B_b),
$$

则

$$
\begin{aligned}
\Pr(H_0^{\mathrm{R}} \mid \mathbf{z}) &= C_0, &
\Pr(H_1^{\mathrm{R}} \mid \mathbf{z}) &= C_2 - C_0, &
\Pr(H_2^{\mathrm{R}} \mid \mathbf{z}) &= C_1 - C_0, \\
\Pr(H_3^{\mathrm{R}} \mid \mathbf{z}) &= A - C_1 - C_2 + C_0, &
\Pr(H_4^{\mathrm{R}} \mid \mathbf{z}) &= 1 - A.
\end{aligned}
$$

*Proof.* $$H_4^{\mathrm{R}}$$ 不成立当且仅当每个部分都不是 $$H_4$$，概率为 $$A$$。“性状 1 在区域内无因果”当且仅当每个部分为 $$H_0$$ 或 $$H_2$$，概率 $$C_1$$，且它蕴含“无共享”；同理 $$C_2$$；两者同时成立的概率为 $$C_0$$。$$H_3^{\mathrm{R}}$$ = “无共享”减去“至少一个性状无因果”，由容斥原理得 $$A - (C_1 + C_2 - C_0)$$。其余直接由定义得出。$$\square$$

**Remarks.**

1. **$$m = 1$$ 时退化为 coloc.abf.** 不切块的位点（gate 未通过）$$\mathcal{P} = \{\text{整区域}\}$$，区域后验就是整区域 coloc.abf 的后验。所以这个统计量与现行方法在 gate 上完全一致，只在切块的位点上不同。
2. **与取最大的关系.** 对同一组互不相交的部分，$$\max_b \mathrm{PP4}(B_b) \le \Pr(H_4^{\mathrm{R}} \mid \mathbf{z}) \le \min\{1, \sum_b \mathrm{PP4}(B_b)\}$$。只有一个部分有共享证据时，区域后验约等于取最大；有多个时（如 double_shared2 的两个共享信号各成一块）区域后验更大。现行 block-coloc 的块可以重叠，重叠的 SNP 重新分配后各部分的 PP4 会变，所以这个不等式不直接适用于现行方法的输出。
3. **区域先验与整区域 coloc.abf 在首项上相同.** 由 (P4) 与 Proposition B 计算区域层面的先验。在 $$Mp_1, Mp_2 \ll 1$$ 的一阶近似下，$$\Pr(H_4^{\mathrm{R}}) \approx \sum_b M_bp_{12} = Mp_{12}$$；$$\Pr(H_3^{\mathrm{R}})$$ 来自同一部分内的不同 SNP 对（$$\sum_b M_b(M_b-1)p_1p_2$$）与不同部分的 SNP 对（$$\sum_{a\ne b}M_aM_bp_1p_2$$），合计 $$M(M-1)p_1p_2$$。二者的首项都与整区域 coloc.abf 的先验相同（Proposition 7）；归一化和高阶项使它们不完全相等。精确计算（Codex，`p2_exact_priors.tsv`）：$$M = 940$$ 时区域先验 $$\Pr(H_4^{\mathrm{R}})$$ 比整区域高 1.7%（$$\{45, 895\}$$）、3.4%（$$\{45, 45, 850\}$$）、5.4%（$$\{45, 45, 30, 30, 790\}$$），均分四块时高 14.5%，每个 SNP 自成一块时高 20%。数值上，$$M = 940$$ 时整区域的先验比 $$H_4 : H_3 = 1.06$$；划分为 $$\{45, 895\}$$、$$\{45, 45, 850\}$$、$$\{45, 45, 30, 30, 790\}$$ 时区域层面为 1.07、1.08、1.10；单独一个 45 SNP 的块为 22.7。

   这为 Results §8 的两难提供了一个解释：每 SNP 先验在单个块内偏向共享，但把整个区域都纳入划分后，块内缺少的 $$H_3$$ 先验质量出现在“性状 1 的因果在一块、性状 2 的因果在另一块”的组合里，区域层面对共享的偏向只有上述几个百分点；区域等价先验则把整区域的 $$H_3$$ 质量全部压进一个块内。

4. **no valid block 不再需要特殊规则.** 两性状的 lead 邻域互不相交时，划分为“只含性状 1 lead 的邻域”“只含性状 2 lead 的邻域”和其余 SNP。若前者的后验集中在 $$H_1$$、后者集中在 $$H_2$$，由 Proposition B，$$H_3^{\mathrm{R}}$$ 的后验大。但 no valid block 本身不等于两性状各有真实信号且不在一起：lead 数不是因果数（默认格的 single_shared 中也有 7 个这样的测试位点）。实际表现见 §6.3。
5. **区域 PP3 有明确含义.** traceCB 规则中的 PP3 改用 $$\Pr(H_3^{\mathrm{R}} \mid \mathbf{z})$$。

## 4. Construction of the partition

划分取自现行 block-coloc 的中间结果，不引入新的调节参数：

1. lead 与 gate 同 Algorithm 1（Step 1–3）；gate 未通过时 $$\mathcal{P} = \{\text{整区域}\}$$。
2. 候选部分：Step 4 的联合块；不属于任何联合块的性状 1 lead 的邻域 $$\mathcal{W}(l)$$；同样的性状 2 lead 的邻域。
3. 候选部分可能重叠（同一 lead 出现在两个按符号拆分的块中，或同一性状两个 lead 的邻域相交）。重叠的 SNP 归入其 $$\hat r^2$$ 最大的 lead 所在的部分，使各部分互不相交。
4. 其余 SNP 组成一个部分（rest）。

对每个部分运行 coloc.abf（每 SNP 先验），按 Proposition B 合成。

## 5. When the assumptions fail

- **(P2) 块间 LD.** 邻域按与 lead 的 $$\hat r^2 \ge \rho_b = 0.4$$ 定义。这只限制 SNP 与 lead 的相关，不限制不同部分的两个 SNP 之间的相关：可以构造两个互不相关的 lead，各自邻域里的两个 SNP 彼此相关 0.999（REVIEW_P2 §4）。所以部分之间的 LD 没有 0.63 的上界。块外因果 $$c$$ 会使块内 SNP $$j$$ 的边际均值偏移 $$r_{jc}\lambda_c$$（式 3.1）。只有一个块外因果时偏移不超过 $$\max_{j \in B}\lvert r_{jc}\rvert\,\lambda_c$$，多个时为各项之和；默认设定下 $$\lambda_1 = 5.5$$–7.8，$$\lvert r\rvert = 0.63$$ 时即达 3.5–4.9。不含 rest 的版本作为敏感性分析报告，但去掉 rest 并不能恢复块间独立。给出这一偏移对 PP4 影响的界，是下一步的理论工作（对应 Model 页列出的“切块后 SCV 在多大程度上成立”）。
- **(P3) 同一部分内两个因果.** 同一性状两个因果 $$\lvert r\rvert \ge 0.63$$ 时落进同一部分（Results §6.4 的 high 档），部分内单因果假设失效，与 coloc.abf 的失效方式相同。
- **(P1) 划分依赖数据.** 划分由 lead 的位置决定，而 lead 来自同一组 $$\mathbf{z}$$。因此这里的“区域后验”是选定划分下工作模型的 plug-in 后验，不是某个包含选择机制的统一模型的后验；没有推导选择带来的偏差。误报率匹配检查的是判定的误报，不是概率的校准。

## 6. Comparison on the stage-1 loci

与第一段完全相同的位点（逐位点核对：现行 block-coloc 的 PP4 与 Results §8 的差为 0，默认格 9,600 个、$$h^2_2 = 5\%$$ 格 4,800 个位点；Codex 独立穷举 200 个多部分位点，区域 PP4 / PP3 与本页公式的最大差 $$< 10^{-13}$$）。脚本：`scripts/block_partition.R`、`partition_eval.py`；完整表格：`outputs/bulk_h2ref_v1_partition*/PARTITION.md`。

| 记号 | 定义 |
|---|---|
| whole | 整区域 coloc.abf |
| bc_stage1 | 现行 block-coloc：联合块取最大，gate 未通过或 no valid block 时退回整区域 |
| joint_prod | 只在联合块上按 $$1 - \prod(1 - \mathrm{PP4})$$ 合成，其余同 bc_stage1 |
| region | Proposition B，全部部分 |
| region_nr | 同上但去掉 rest |

### 6.1 Matched 5% false-positive rate

测试窗口检出率（%），判定线在校准窗口上各自重新定：

| 设定 | whole | bc_stage1 | joint_prod | region | region_nr |
|---|---|---|---|---|---|
| 默认，合并 | 73.2 | 76.8 | 76.8 | 72.5 | 74.1 |
| 默认，double_shared1 | 45.3 | 54.8 | 54.8 | 57.7 | 59.2 |
| 默认，single_double_shared | 74.3 | 74.2 | 74.2 | 62.5 | 65.2 |
| $$h^2_2 = 5\%$$，合并 | 89.8 | 97.7 | 97.7 | 97.6 | 97.8 |

配对比较（百分点；位点层面 95% CI；6 个测试窗口中差 > 0 的个数）：region − bc_stage1：默认 −4.3 (−5.7, −2.9)，0/6；$$h^2_2 = 5\%$$ −0.1 (−1.0, +0.9)，2/6。region − whole：默认 −0.8 (−2.3, +0.8)，2/6；$$h^2_2 = 5\%$$ +7.8 (+6.0, +9.5)，6/6。joint_prod 与 bc_stage1 在各自判定线下的判定在两格中完全相同（连续分数不一定相等）。

- **在相同误报率下，区域后验不优于取最大**：默认设定低 4.3 个百分点（6 个窗口全部为负），eQTL 强时相同。
- “取最大 vs 只在联合块上合成”没有改变任何判定。
- 额外阴性：$$h^2_2 = 5\%$$ 时 gwas_only 被判为共定位的比例，whole 与 bc_stage1 为 100%，region 为 0.3%（Results §8.3 指出的按误报率定线的盲区在区域后验下消失，因为它的判定线是 0.49 而不是 0.014）。

### 6.2 traceCB rule

测试窗口，%（区域方法的 PP3 为 $$\Pr(H_3^{\mathrm{R}}\mid\mathbf{z})$$）：

| 设定 | 统计量 | 共定位位点：判共定位 / 判独立 / 不确定 | single_distinct：判独立（判共定位） | double_distinct：判独立（判共定位） |
|---|---|---|---|---|
| 默认 | whole | 36.5 / 4.1 / 59.3 | 35.7 (0.2) | 10.2 (0.0) |
| | bc_stage1 | 44.1 / 2.7 / 53.2 | 35.5 (0.2) | 9.8 (0.2) |
| | region | 48.8 / 1.2 / 50.0 | 34.3 (0.3) | 7.3 (0.8) |
| $$h^2_2 = 5\%$$ | whole | 68.2 / 20.4 / 11.3 | 100 (0.0) | 90.0 (0.0) |
| | bc_stage1 | 92.3 / 4.8 / 2.8 | 98.3 (0.3) | 88.0 (0.0) |
| | region | 95.8 / 1.4 / 2.8 | 96.3 (1.7) | 72.7 (2.0) |

固定规则下区域后验多判出 3.5–4.7 个百分点的共定位，但不共定位位点被判为共定位的比例也升高（最高 2.0%，bc_stage1 最高 0.3%），double_distinct 被明确判为独立的比例下降（$$h^2_2 = 5\%$$：72.7% 对 88.0%）。

### 6.3 Why: residual shared mass in every single-trait part

§3 Remark 4 的预期只对了一半。no valid block 的位点上，区域 PP3 确实升高（$$h^2_2 = 5\%$$ 的 double_distinct：中位数 0.82，整区域 PP4 中位数 0.000），但区域 PP4 也升到 0.09–0.15。一个解释是 [Model](01_bulk_model.md) 式 (5.3) 的残留共享质量在每个部分里各出现一次：只有性状 1 有强信号的单个 SNP 上，$$\mathrm{PP4} = p_{12}\mathrm{ABF}_2/(p_1 + p_{12}\mathrm{ABF}_2)$$，性状 2 的 ABF 为 1 时是 0.09，为 0.1 时是 0.01。它不是下限：性状 2 的样本量足够大、数据支持零效应时它趋于 0；$$N_2 = 1{,}000$$ 时 ABF 降不了多少。测试窗口阴性位点上各类部分的 PP4：

| 部分 | 默认：中位数（95% 分位） | $$h^2_2 = 5\%$$：中位数（95% 分位） |
|---|---|---|
| 只含性状 1 lead 的邻域 | 0.04 (0.13–0.14) | 0.04 (0.18) |
| 只含性状 2 lead 的邻域 | 0.002–0.004 (0.03–0.05) | 0.01 (0.05–0.12) |
| rest | 0.02–0.03 (0.11) | 0.02 (0.15–0.23) |

double_distinct 通常有两个性状 1 邻域加 rest，合成后 $$1 - (1 - 0.04)^2(1 - 0.02) \approx 0.10$$，与观察到的中位数 0.12–0.15 一致。整区域 coloc.abf 在两性状都有强信号时没有这个下限，因为 $$H_3$$ 的质量把它挤掉了；分块模型的先验在部分之间独立（P4），每个部分各自保留它。它使阴性位点的 PP4 抬高、判定线升到 0.28–0.49，抵消了阳性位点上的增益。数据选块、块间 LD、部分内多因果和 sdY 的估计同时存在，现有结果不能把原因唯一归于先验独立。

### 6.4 Conclusion

1. **不建议用区域后验替换现行的取最大规则.** 相同误报率下它不更好（默认设定更差），固定规则下多出的检出伴随误报上升。
2. Proposition A、B 仍可作为 block-coloc 的模型解释（块独立近似下，取最大是区域后验的下界），以及说明块内每 SNP 先验在区域层面并不偏向共享（Remark 3）。
3. 一个可能的改进方向是不独立的先验（让一个性状的信号位置进入其他部分的先验）或用经验 Bayes 估计 $$p_{12}$$；二者都需要新的推导，尚未做。
4. **合成公式不是新的.** 见 §7。

## 7. Related work

**已有直接先例（Codex 核对，REVIEW_P2 §5）.** Hukku et al. (2022) 的式 (1) 用 $$\mathrm{GRCP} = 1 - \prod_r (1 - \mathrm{RCP}_r)$$ 把各独立 eQTL 信号的共定位概率合成为“基因至少有一个共定位因果变异”的概率，式 (12) 对 LCP 作同样的合成；fastENLOC 的输出也区分信号层面的 RCP 与基因层面的 GRCP。因此 $$\Pr(H_4^{\mathrm{R}}) = 1 - \prod_b(1 - \mathrm{PP4}_b)$$ 不能作为本文的贡献。coloc.signals（Wallace 2020）与 coloc.susie（Wallace 2021）是逐对信号比较，未发现把各块合成完整区域 $$H_0$$–$$H_4$$ 的同一公式，但这次检索不是系统综述；PWCoCo 的预印本全文未能取得，未核实。把区域按近似独立的 LD 块分解本身也是常见做法（Berisa & Pickrell 2016）。

## References

- Hukku A, Sampson MG, Luca F, Pique-Regi R, Wen X (2022). Analyzing and reconciling colocalization and transcriptome-wide association studies from the perspective of inferential reproducibility. *Am J Hum Genet* 109: 825–837.
- Berisa T, Pickrell JK (2016). Approximately independent linkage disequilibrium blocks in human populations. *Bioinformatics* 32(2): 283–285.
- 其余见 [Model](01_bulk_model.md#references)。
