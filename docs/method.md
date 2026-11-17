# Model, Algorithm, and Analytical Properties

本页从个体水平线性生成模型出发，推导单变量最小二乘统计量的分布、单因果贝叶斯因子的充分性与近似检验、面向多因果区域的 `block-coloc` 算法以及外部参考连锁不平衡（LD）矩阵失配下条件统计量的精确分布。所有形式化极限命题的证明及灵敏度分析见 [Supplementary](supplement.md)。

## 1. Individual-level generative model and summary-statistic distributions

考虑包含 $$M$$ 个单核苷酸多态性（SNP）的顺式（cis）基因组区域，设连续性状 $$k \in \{1, 2\}$$ 分别在互不重叠的独立队列中测量，样本量为 $$N_k$$（主模拟取 $$N_1 = 20{,}000$$ 对应 GWAS，$$N_2 = 1{,}000$$ 对应 bulk eQTL）。区域内基因型剂量 $$G_{k,ij} \in \{0, 1, 2\}$$（个体 $$i$$、位点 $$j$$）经列中心化得 $$\widetilde{G}_{k,ij} = G_{k,ij} - 2\hat f_{kj}$$，其中 $$\hat f_{kj} = (2 N_k)^{-1}\sum_{i=1}^{N_k} G_{k,ij}$$ 为计数等位基因频率，$$\mathrm{MAF}_{kj} = \min(\hat f_{kj}, 1 - \hat f_{kj})$$。记第 $$j$$ 列中心化剂量的平方和与样本标准差分别为 $$Q_{k,j} = \sum_{i=1}^{N_k} \widetilde{G}_{k,ij}^2$$ 与 $$d_{kj} = \sqrt{Q_{k,j}/(N_k - 1)}$$（在 Hardy–Weinberg 平衡下 $$d_{kj} \approx \sqrt{2\hat f_{kj}(1 - \hat f_{kj})}$$）。令样本内标准化基因型列为 $$\mathbf{x}_{kj} = \widetilde{\mathbf{G}}_{k,\cdot j} / d_{kj}$$，满足 $$\mathbf{1}^\top \mathbf{x}_{kj} = 0$$、$$\mathbf{x}_{kj}^\top \mathbf{x}_{kj} = N_k - 1$$，样本内 LD 相关矩阵精确为 $$\mathbf{R}_k = (N_k - 1)^{-1}\mathbf{X}_k^\top \mathbf{X}_k = (r_{k,jl})$$。

个体水平线性模型设为

$$
\mathbf{y}_k = \mathbf{X}_k \boldsymbol\beta_k + \boldsymbol\varepsilon_k, \qquad \boldsymbol\varepsilon_k \sim \mathcal{N}\!\left(\mathbf{0},\ \sigma^2_k \mathbf{I}_{N_k}\right), \qquad \boldsymbol\varepsilon_1 \perp \boldsymbol\varepsilon_2 \,. \tag{1.1}
$$

记因果指示变量 $$\gamma_{kj} \in \{0, 1\}$$、因果集合 $$\mathcal{C}_k = \{j : \gamma_{kj} = 1\}$$ 及因果个数 $$K_k = \lvert\mathcal{C}_k\rvert$$。区域水平共定位检验聚焦五类互斥穷尽假设：

$$
H_0^{\mathrm{R}} : \mathcal{C}_1 = \mathcal{C}_2 = \varnothing, \quad
H_1^{\mathrm{R}} : \mathcal{C}_1 \neq \varnothing,\, \mathcal{C}_2 = \varnothing, \quad
H_2^{\mathrm{R}} : \mathcal{C}_1 = \varnothing,\, \mathcal{C}_2 \neq \varnothing, \quad
H_3^{\mathrm{R}} : \mathcal{C}_1 \neq \varnothing,\, \mathcal{C}_2 \neq \varnothing,\, \mathcal{C}_1 \cap \mathcal{C}_2 = \varnothing, \quad
H_4^{\mathrm{R}} : \mathcal{C}_1 \cap \mathcal{C}_2 \neq \varnothing \,. \tag{1.2}
$$

在主模拟中，当 $$K_k \ge 1$$ 时，标准化因果效应定义为

$$
\beta_{kj} = \xi_k\, \gamma_{kj}\, s_{kj} \sqrt{q_k}, \qquad q_k = \frac{h^2_k}{K_k}, \qquad s_{kj} \overset{\mathrm{iid}}{\sim} \mathrm{Unif}\{-1, +1\} \,, \tag{1.3}
$$

两性状符号相互独立；线性遗传预测向量 $$\boldsymbol\eta_k = \mathbf{X}_k \boldsymbol\beta_k$$ 满足 $$\mathbf{1}^\top \boldsymbol\eta_k = 0$$，其样本内经验方差为

$$
\operatorname{Var}_{\mathrm{emp}}(\boldsymbol\eta_k) = \frac{1}{N_k - 1}\boldsymbol\beta_k^\top \mathbf{X}_k^\top \mathbf{X}_k \boldsymbol\beta_k = \boldsymbol\beta_k^\top \mathbf{R}_k \boldsymbol\beta_k = \xi_k^2 q_k\left(K_k + \sum_{j \neq l \in \mathcal{C}_k} s_{kj} s_{kl}\, r_{k,jl}\right) \,. \tag{1.4}
$$

令式 (1.4) 等于名义局部遗传度 $$h^2_k$$，即得精确重标定因子 $$\xi_k = \big(1 + K_k^{-1}\sum_{j \neq l \in \mathcal{C}_k} s_{kj} s_{kl}\, r_{k,jl}\big)^{-1/2}$$（$$K_k = 1$$ 时 $$\xi_k = 1$$；$$K_k = 2$$ 且 $$\lvert r_{k,12}\rvert < 0.2$$ 时 $$\xi_k \in [0.91, 1.12]$$）。配合残差方差 $$\sigma^2_k = 1 - h^2_k$$（当 $$K_k = 0$$ 时 $$\boldsymbol\eta_k = \mathbf{0}$$ 且 $$\sigma^2_k = 1$$），固定基因型下表型的条件期望样本方差为单位尺度：$$\mathbb{E}[\hat\sigma_{yk}^2 \mid \mathbf{X}_k] = \boldsymbol\beta_k^\top \mathbf{R}_k \boldsymbol\beta_k + \sigma^2_k = 1$$。

*Remark 1.* 在式 (1.3)–(1.4) 下，$$\xi_k$$ 依赖于因果符号与因果间 LD；对 $$K_k = 2$$（相关 $$r = r_{k,ab}$$）的四种符号配置求和可得 $$\mathbb{E}[\beta_{ka}^2 \mid \mathcal{C}_k] = q_k / (1 - r^2)$$ 及 $$\mathbb{E}[\beta_{ka}\beta_{kb} \mid \mathcal{C}_k] = -q_k r / (1 - r^2)$$。因此，区域平均每 SNP 遗传度关系 $$\omega_k := \mathbb{E}[\beta_{kj}^2] = (K_k/M)q_k = h^2_k / M$$ 与不相关 LD score 期望 $$\mathbb{E}[z_{kj}^2] \approx 1 + (N_k h^2_k / M)\sum_{l=1}^M r_{k,jl}^2$$（Bulik-Sullivan et al. 2015；Jiang et al. 2026）仅在**无重标定且互不相关的零均值随机效应模型**（例如 `effect = "normal"` 灵敏度分析）下精确成立；在固定重标定设计中，它仅作为因果近正交低 LD 情形下忽略 $$O(r^2)$$ 的启发式近似。一般地，无条件边际方差由 $$(\mathbf{R}_k \operatorname{Cov}(\boldsymbol\beta_k)\mathbf{R}_k)_{jj}$$ 给出。

**单变量最小二乘（OLS）估计与精确 $$t$$ 统计量.** 令去均值表型为 $$\tilde{\mathbf{y}}_k = (\mathbf{I}_{N_k} - N_k^{-1}\mathbf{1}\mathbf{1}^\top)\mathbf{y}_k$$，样本表型方差为 $$\hat\sigma_{yk}^2 = (N_k - 1)^{-1}\lVert\tilde{\mathbf{y}}_k\rVert_2^2$$。因 $$\mathbf{1}^\top \mathbf{x}_{kj} = 0$$，SNP $$j$$ 在每等位基因与标准化基因型尺度上的单变量 OLS 斜率分别为 $$\hat\beta^{\mathrm{allele}}_{kj} = \widetilde{\mathbf{G}}_{k,\cdot j}^\top \tilde{\mathbf{y}}_k / Q_{k,j}$$ 与 $$\hat b_{kj} = d_{kj}\hat\beta^{\mathrm{allele}}_{kj} = (N_k - 1)^{-1}\mathbf{x}_{kj}^\top \tilde{\mathbf{y}}_k$$，满足线性分解

$$
\hat b_{kj} = b_{kj} + e_{kj}, \qquad b_{kj} = (\mathbf{R}_k \boldsymbol\beta_k)_j = \sum_{c \in \mathcal{C}_k} r_{k,jc}\beta_{kc}, \qquad \mathbf{e}_k \mid \mathbf{X}_k \sim \mathcal{N}\!\left(\mathbf{0},\ \frac{\sigma^2_k}{N_k - 1}\mathbf{R}_k\right) \,. \tag{1.5}
$$

残差平方和为 $$\mathrm{RSS}_{kj} = \lVert\tilde{\mathbf{y}}_k\rVert_2^2 - (N_k - 1)\hat b_{kj}^2$$（自由度 $$\nu_k = N_k - 2$$）。两尺度上的 OLS 标准误平方估计量为

$$
\hat s_{kj}^2 = \frac{\mathrm{RSS}_{kj}}{(N_k - 2)(N_k - 1)} = \frac{\hat\sigma_{yk}^2 - \hat b_{kj}^2}{N_k - 2}, \qquad V_{kj} := \left(\hat s^{\mathrm{allele}}_{kj}\right)^2 = \frac{\hat s_{kj}^2}{d_{kj}^2} \approx \frac{\sigma_{yk}^2}{2 f_j (1 - f_j) N_k} \,, \tag{1.6}
$$

其中 $$\sigma_{yk}^2 := \mathbb{E}[\hat\sigma_{yk}^2 \mid \mathbf{X}_k] = 1$$。令样本 Pearson 相关系数为 $$\hat\rho_{kj} = \hat b_{kj} / \hat\sigma_{yk}$$，则检验统计量 $$z_{kj} := \hat\beta^{\mathrm{allele}}_{kj} / \hat s^{\mathrm{allele}}_{kj} = \hat b_{kj} / \hat s_{kj}$$ 与 $$\hat\rho_{kj}$$ 满足精确恒等式

$$
z_{kj} = \sqrt{N_k - 2}\,\frac{\hat\rho_{kj}}{\sqrt{1 - \hat\rho_{kj}^2}} =: \psi_{N_k}(\hat\rho_{kj}) \,, \tag{1.7}
$$

双侧 $$p$$ 值为 $$p_{kj} = 2 F_{t_{\nu_k}}(-\lvert z_{kj}\rvert)$$，在完全零假设 $$\boldsymbol\beta_k = \mathbf{0}$$ 下 $$z_{kj} \mid \mathbf{X}_k \sim t_{\nu_k}$$ 精确成立。

*Remark 2.* 仅当区域不存在与 $$j$$ 不共线的其他因果变异时，$$\mathrm{RSS}_{kj} / (N_k - 2)$$ 对生成噪声方差 $$\sigma^2_k$$ 无偏；若存在遗漏因果，残差吸收了确定性未解释遗传成分，$$\hat s_{kj}^2$$ 估计的是边际回归残差方差而非 $$\sigma^2_k$$。

**汇总统计量多元正态近似模型与因果非中心参数.** 在固定 $$\boldsymbol\beta_k$$ 下，由式 (1.6) 与 $$\hat\sigma_{yk}^2 = 1 + O_p(N_k^{-1/2})$$ 得

$$
\hat s_{kj}^2 = \frac{\sigma_{yk}^2 - b_{kj}^2}{N_k} + O_p\!\left(N_k^{-3/2}\right) \,. \tag{1.8}
$$

固定信号下 $$-b_{kj}^2 / N_k$$ 仍为 $$O(N_k^{-1})$$ 主项（仅在局部微弱信号 $$b_{kj}^2 = O(N_k^{-1/2})$$ 渐近区划下并入余项）。在本研究的小方差贡献区间（$$q_1 \le 0.02, q_2 \le 0.05$$）内，忽略 $$b_{kj}^2 \ll 1$$ 与 $$\sigma^2_k \approx 1$$ 作为额外小效应近似，结合式 (1.5)（将 $$N_k - 1$$ 替换为 $$N_k$$）得到汇总统计量多元正态近似模型（Zou et al. 2022）：

$$
\mathbf{z}_k \mid \boldsymbol\beta_k \;\dot\sim\; \mathcal{N}\!\left(\sqrt{N_k}\,\mathbf{R}_k \boldsymbol\beta_k,\ \sigma^2_k \mathbf{R}_k\right) \approx \mathcal{N}\!\left(\sqrt{N_k}\,\mathbf{R}_k \boldsymbol\beta_k,\ \mathbf{R}_k\right) \,. \tag{1.9}
$$

在本研究固定的样本内标准化基因型设计下，若性状仅含单个因果变异 $$c$$（或所有因果相互正交），式 (1.7) 的精确有限样本分布为非中心 $$t$$ 分布

$$
z_{kc} = \frac{Z + s_{kc}\delta_k}{\sqrt{U / \nu_k}}, \qquad Z \sim \mathcal{N}(0, 1) \perp U \sim \chi^2_{\nu_k}, \qquad \delta_k = \sqrt{\frac{(N_k - 1)\, q_k}{1 - q_k}} \approx \lambda_k := \sqrt{\frac{N_k\, q_k}{1 - q_k}} \,. \tag{1.10}
$$

该固定设计非中心 $$t$$ 统计量的渐近方差为 $$1 + q_k / \{2(1 - q_k)\}$$，不同于 i.i.d. 双变量正态随机设计下由 Pearson Delta 方法（$$\psi_{N_k}'(\rho)$$）导出的 $$(1 - q_k)^{-1}$$（完整推导与 100 万次 Monte Carlo 验算见 [Supplementary](supplement.md) Remark S1）。在 $$q_k \ll 1$$ 下，单点基准检验效能约为 $$\Phi(\lambda_k - z_{\alpha/2}) + \Phi(-\lambda_k - z_{\alpha/2})$$，其中 $$z_{\alpha/2} = \Phi^{-1}(1 - \alpha/2)$$。

---

## 2. Single-causal-variant colocalisation (`coloc.abf`)

Giambartolomei et al. (2014) 在每性状至多一个因果变异（SCV，$$K_k \le 1$$）假设下枚举五类互斥穷尽配置：$$H_0$$（$$\mathcal{C}_1 = \mathcal{C}_2 = \varnothing$$）、$$H_1$$（$$\lvert\mathcal{C}_1\rvert = 1, \mathcal{C}_2 = \varnothing$$）、$$H_2$$（$$\mathcal{C}_1 = \varnothing, \lvert\mathcal{C}_2\rvert = 1$$）、$$H_3$$（$$\mathcal{C}_1 = \{i\}, \mathcal{C}_2 = \{j\}, i \neq j$$）与 $$H_4$$（$$\mathcal{C}_1 = \mathcal{C}_2 = \{j\}$$）。

**Assumption 1.** *在汇总统计量多元正态近似模型 (1.9) 下，零假设 $$\boldsymbol\beta_k = \mathbf{0}$$ 与单因果备择假设 $$\boldsymbol\beta_k = \beta_{kj}\mathbf{e}_j$$ 具有相同的残差方差参数 $$\sigma^2_k$$（单变异微效近似）；当子矩阵 $$\mathbf{R}_{k,BB}$$ 秩亏（如参考样本量 $$n_{\mathrm{ref}} < \lvert B\rvert$$）时，边际统计量向量满足线性模型相容条件 $$\mathbf{z}_{k,B} \in \mathrm{col}(\mathbf{R}_{k,BB})$$。*

**Lemma 1（单因果假设下边际统计量的充分性）.** *在 Assumption 1 下，给定单一因果位点 $$j$$ 的边际统计量 $$z_{kj}$$，其余变异 $$\mathbf{z}_{k,-j}$$ 的条件分布不含效应参数 $$\beta_{kj}$$ 且与零假设下完全一致；当 $$\mathbf{z}_{k,B} \in \mathrm{col}(\mathbf{R}_{k,BB})$$ 时，基于 Moore–Penrose 广义逆 $$\mathbf{R}_{k,BB}^+$$ 的似然比精确退化为单点 $$z_{kj}$$ 的边际似然比。*

*Proof.* 在共享协方差下，$$\mathbb{E}[\mathbf{z}_k] = \sqrt{N_k}\beta_{kj}\mathbf{r}_{k,\cdot j}$$，故 $$\mathbf{z}_{k,-j} \mid z_{kj} \sim \mathcal{N}\big(\mathbf{r}_{k,-j,j} z_{kj},\ \sigma^2_k(\mathbf{R}_{k,-j,-j} - \mathbf{r}_{k,-j,j}\mathbf{r}_{k,-j,j}^\top)\big)$$ 与 $$\beta_{kj}$$ 无关。对秩亏情形，因 $$\mathbf{e}_j^\top \mathbf{R}_{k,BB}\mathbf{R}_{k,BB}^+ \mathbf{z}_{k,B} = z_{kj}$$ 对一切 $$\mathbf{z}_{k,B} \in \mathrm{col}(\mathbf{R}_{k,BB})$$ 恒成立，似然比二次型中的广义逆仅通过 $$z_{kj}$$ 起作用。$$\square$$

*Remark 3.* 当性状存在至少两个因果变异（$$K_k \ge 2$$）时，条件均值残留项 $$\sqrt{N_k}\sum_{l \neq j} \beta_{kl}(\mathbf{r}_{k,-j,l} - r_{k,jl}\mathbf{r}_{k,-j,j}) \neq \mathbf{0}$$，单点充分性失效。此外，由真实 LD 生成的统计量 $$\mathbf{z}_{k,B}$$ 并不自动落入外部参考面板矩阵 $$\widehat{\mathbf{R}}_{BB}$$ 的列空间中。

**Wakefield 近似贝叶斯因子（ABF）与先验尺度校准.** 为与 LD 相关系数 $$r_{k,jl}$$ 区分，记方差收缩因子为 $$a_{kj}$$。在每等位基因尺度下，令渐近抽样分布为 $$\hat\beta^{\mathrm{allele}}_{kj} \mid \beta^{\mathrm{allele}}_{kj} \sim \mathcal{N}(\beta^{\mathrm{allele}}_{kj}, V_{kj})$$，赋予共轭因果先验 $$\beta^{\mathrm{allele}}_{kj} \sim \mathcal{N}(0, W_k)$$。按线性尺度变换 $$\beta^{\mathrm{std}}_{kj} = (d_{kj}/\sigma_{yk})\beta^{\mathrm{allele}}_{kj}$$ 换算时，对应标准化先验方差乘以平方尺度比 $$(d_{kj}/\sigma_{yk})^2 \approx 2 f_j (1 - f_j)/\sigma_{yk}^2$$（密度变换雅可比因子为 $$\sigma_{yk}/d_{kj}$$）。对先验积分得备择边缘分布 $$\mathcal{N}(0, V_{kj} + W_k)$$，两正态密度求比并配方化简给出（Wakefield 2009）：

$$
\mathrm{ABF}_{kj} = \frac{\phi\!\left(\hat\beta^{\mathrm{allele}}_{kj};\, 0, V_{kj} + W_k\right)}{\phi\!\left(\hat\beta^{\mathrm{allele}}_{kj};\, 0, V_{kj}\right)} = \sqrt{1 - a_{kj}}\,\exp\!\left(\tfrac{1}{2} a_{kj} z_{kj}^2\right), \qquad a_{kj} = \frac{W_k}{V_{kj} + W_k} \in (0, 1) \,. \tag{2.1}
$$

*Remark 4.* 在固定 $$(V_{kj}, W_k)$$ 的标准正态近似 $$z_{kj} \sim \mathcal{N}(0, 1)$$ 下，由 $$\chi^2_1$$ 矩母函数立得 $$\mathbb{E}[\mathrm{ABF}_{kj}] = 1$$；在备择 $$z_{kj} \sim \mathcal{N}(\mu_{kj}, 1)$$ 下，$$\mathbb{E}[\log\mathrm{ABF}_{kj}] = \tfrac{1}{2}\log(1 - a_{kj}) + \tfrac{1}{2} a_{kj}(\mu_{kj}^2 + 1)$$。精确有限样本 OLS 统计量服从重尾 $$t_{\nu_k}$$ 分布（对任意固定 $$a_{kj} > 0$$ 有 $$\mathbb{E}[\exp(\tfrac{1}{2} a_{kj} z_{kj}^2)] = +\infty$$），且实际代入估计量（plug-in）$$\hat a_{kj}$$ 依赖估计方差与表型尺度 $$\widehat{\mathrm{sdY}}_k$$，故上述有限矩恒等式仅在固定参数正态近似下成立。数值量级上，由式 (1.3) 在 $$\mathrm{MAF} \in [0.05, 0.50]$$ 范围换算，默认 eQTL（$$q_2 = 0.015$$）每等位基因效应幅值为 $$0.17\,\sigma_{y2}$$–$$0.40\,\sigma_{y2}$$（$$a_{2j} \in [0.68, 0.92]$$），默认 GWAS（$$q_1 = 0.003$$）为 $$0.08\,\sigma_{y1}$$–$$0.18\,\sigma_{y1}$$（$$a_{1j} \in [0.977, 0.996]$$），均与 `coloc` 默认先验标准差 $$\sqrt{W_k} = 0.15\,\sigma_{yk}$$ 处于同一尺度。

赋予每 SNP 先验概率 $$(p_1, p_2, p_{12}) = (10^{-4}, 10^{-4}, 10^{-5})$$，记单点 ABF 求和 $$S_k = \sum_{j=1}^M \mathrm{ABF}_{kj}$$ 与跨性状乘积和 $$S_{12} = \sum_{j=1}^M \mathrm{ABF}_{1j}\mathrm{ABF}_{2j}$$，五假设后验与 $$H_4 : H_3$$ 先验赔率分别为

$$
\Pr(H_0 \mid D) \propto 1, \quad \Pr(H_1 \mid D) \propto p_1 S_1, \quad \Pr(H_2 \mid D) \propto p_2 S_2, \quad \Pr(H_3 \mid D) \propto p_1 p_2 (S_1 S_2 - S_{12}), \quad \Pr(H_4 \mid D) \propto p_{12} S_{12}, \tag{2.2}
$$

$$
\frac{\Pr(H_4)}{\Pr(H_3)} = \frac{p_{12}}{(M - 1)\, p_1 p_2} \,. \tag{2.3}
$$

在默认先验 $$p_{12}/(p_1 p_2) = 10^3$$ 下，全区域 $$M = 938$$（中位数）时先验赔率约为 $$1.07$$，而在中位大小 $$M_B = 42$$ 的局部块内先验赔率升至约 $$24.4$$（[Figure 5](fig5.md)）。

---

## 3. Block-partitioned colocalisation (`block-coloc`)

当至少一个性状含多个因果变异时，全区域 `coloc.abf` 渐近检验的是两性状最强边际信号是否重合，而非区域内是否存在共享因果变异（[Supplementary](supplement.md) Proposition S1；[Figure 3](fig3.md)）。`block-coloc` 将外部参考 LD $$\widehat{\mathbf{R}}$$ 仅用于离散的多信号计数与候选块边界划分，而块内共定位仍直接输入未经过条件调整的原始边际汇总统计量。

**分块多元正态条件分布与逐步条件统计量.** 设已选条件独立变异索引集为 $$\mathcal{S} \subset \{1, \dots, M\}$$，待检变异为 $$j \notin \mathcal{S}$$。在正态近似模型 (1.9) 下，记 $$\boldsymbol\mu_k = \sqrt{N_k}\,\mathbf{R}_k \boldsymbol\beta_k$$，分块联合分布及给定 $$\mathbf{z}_{k\mathcal{S}}$$ 的 Schur 补条件分布为

$$
\begin{pmatrix} \mathbf{z}_{k\mathcal{S}} \\ z_{kj} \end{pmatrix} \,\middle|\, \boldsymbol\beta_k \;\sim\; \mathcal{N}\!\left(\begin{pmatrix} \boldsymbol\mu_{k\mathcal{S}} \\ \mu_{kj} \end{pmatrix},\ \sigma^2_k \begin{pmatrix} \mathbf{R}_{k,\mathcal{S}\mathcal{S}} & \mathbf{R}_{k,\mathcal{S}j} \\ \mathbf{R}_{k,j\mathcal{S}} & 1 \end{pmatrix}\right), \tag{3.1}
$$

$$
z_{kj} \mid \mathbf{z}_{k\mathcal{S}} \;\sim\; \mathcal{N}\!\left(\mu_{kj} + \mathbf{R}_{k,j\mathcal{S}}\mathbf{R}_{k,\mathcal{S}\mathcal{S}}^{-1}(\mathbf{z}_{k\mathcal{S}} - \boldsymbol\mu_{k\mathcal{S}}),\ \sigma^2_k\left(1 - \mathbf{R}_{k,j\mathcal{S}}\mathbf{R}_{k,\mathcal{S}\mathcal{S}}^{-1}\mathbf{R}_{k,\mathcal{S}j}\right)\right) \,. \tag{3.2}
$$

当 $$\mathcal{C}_k \subseteq \mathcal{S}$$ 时，$$\mathbf{R}_{k,j\mathcal{S}}\mathbf{R}_{k,\mathcal{S}\mathcal{S}}^{-1}\boldsymbol\mu_{k\mathcal{S}} = \mu_{kj}$$，条件均值中来自已入选因果的偏移相消，从而消除中度连锁伴随信号的干扰；若 $$j \notin \mathcal{S}$$ 为新增独立因果，条件均值保留 $$\sqrt{N_k}\,\beta_{kj}\sqrt{1 - \mathbf{R}_{k,j\mathcal{S}}\mathbf{R}_{k,\mathcal{S}\mathcal{S}}^{-1}\mathbf{R}_{k,\mathcal{S}j}}$$。实际计算以参考面板 $$\widehat{\mathbf{R}}$$ 配合对角岭参数 $$\eta = 10^{-6}$$ 替代，记 $$\widehat{\mathbf{R}}_{\mathcal{S}\mathcal{S},\eta} = \widehat{\mathbf{R}}_{\mathcal{S}\mathcal{S}} + \eta \mathbf{I}_{\lvert\mathcal{S}\rvert}$$ 与权重行向量 $$\hat{\mathbf{w}}_{j \mid \mathcal{S}} = \widehat{\mathbf{R}}_{j\mathcal{S}}\widehat{\mathbf{R}}_{\mathcal{S}\mathcal{S},\eta}^{-1}$$，定义标准化条件统计量（Yang et al. 2012）：

$$
z_{kj \mid \mathcal{S}} = \frac{z_{kj} - \hat{\mathbf{w}}_{j \mid \mathcal{S}}\mathbf{z}_{k\mathcal{S}}}{\sqrt{1 - \hat{\mathbf{w}}_{j \mid \mathcal{S}}\widehat{\mathbf{R}}_{\mathcal{S}j}}} \,. \tag{3.3}
$$

*Remark 5.* 当条件方差因子 $$1 - \hat{\mathbf{w}}_{j \mid \mathcal{S}}\widehat{\mathbf{R}}_{\mathcal{S}j} \le 10^{-8}$$（近共线）时将 $$z_{kj \mid \mathcal{S}}$$ 置 $$0$$。由于集合 $$\mathcal{S}$$ 由同一样本通过边际显著性检验数据驱动选出，选择后残余统计量不再严格服从未经截断的标准正态分布；仅在给定固定超集 $$\mathcal{C}_k \subseteq \mathcal{S}$$、$$\widehat{\mathbf{R}} = \mathbf{R}_k$$、$$\eta \to 0$$ 及 $$\sigma^2_k \approx 1$$ 的理想设定下，$$z_{kj \mid \mathcal{S}} \sim \mathcal{N}(0, 1)$$ 精确成立。

> **Algorithm 1: `block-coloc`（逐步条件信号识别与局部 LD 邻域共定位）**
>
> - **输入**：两性状边际汇总统计量 $$(\hat\beta^{\mathrm{allele}}_{kj}, \hat s^{\mathrm{allele}}_{kj}, \mathrm{MAF}_{kj}, N_k)_{j=1}^M$$，参考 LD 矩阵 $$\widehat{\mathbf{R}}$$，超参数 $$\alpha_{\mathrm{lead}} = 10^{-3}$$，clump 阈值 $$\rho_c = 0.1$$，邻域阈值 $$\rho_b = 0.4$$，最小块规模 $$m_{\min} = 2$$。
> - **条件独立信号识别（对 $$k \in \{1, 2\}$$ 分别执行）**：
>   1. **边际初筛与 LD 去重（clumping）**：在 $$\{j : p_{kj} < \alpha_{\mathrm{lead}}\}$$ 上按 $$p_{kj}$$ 升序贪心聚类（仅保留与已选集合均满足 $$\hat r^2 < \rho_c$$ 的变异），得有序列表 $$\mathcal{L}^{\mathrm{marg}}_k = (j_1, \dots, j_{L_0})$$（$$L_0 = \lvert\mathcal{L}^{\mathrm{marg}}_k\rvert$$）。
>   2. **逐步条件回归筛选**：若 $$L_0 = 0$$ 令 $$\mathcal{L}_k = \varnothing$$；否则初始化 $$\mathcal{L}_k = \{j_1\}$$，依次对 $$m = 2, \dots, L_0$$ 计算式 (3.3) 的条件统计量 $$z_{k, j_m \mid \mathcal{L}_k}$$，若 $$2\Phi(-\lvert z_{k, j_m \mid \mathcal{L}_k}\rvert) < \alpha_{\mathrm{lead}}$$ 则更新 $$\mathcal{L}_k \leftarrow \mathcal{L}_k \cup \{j_m\}$$。
> - **多信号判定与块内共定位**：
>   3. **多信号筛选条件（事件 $$\mathcal{G}_0$$）**：若不满足 $$\lvert\mathcal{L}_1\rvert \ge 1,\ \lvert\mathcal{L}_2\rvert \ge 1$$ 且 $$\max(\lvert\mathcal{L}_1\rvert, \lvert\mathcal{L}_2\rvert) \ge 2$$，直接返回全区域 `coloc.abf` 统计量 $$\mathrm{PP4}(\{1, \dots, M\})$$。
>   4. **跨性状 LD 邻域合并分块**：对每个 $$l \in \mathcal{L}_k$$ 定义参考 LD 邻域 $$\mathcal{W}(l) = \{j : \hat r_{jl}^2 \ge \rho_b\}$$。当 $$a \in \mathcal{L}_1$$ 与 $$b \in \mathcal{L}_2$$ 满足 $$\mathcal{W}(a) \cap \mathcal{W}(b) \neq \varnothing$$ 时相连，按各连通分量 $$\mathcal{K}$$ 合并邻域得候选块 $$B_{\mathcal{K}} = \bigcup_{l \in \mathcal{K}} \mathcal{W}(l)$$（若分量内同一性状多个条件独立位点的边际 $$z$$ 异号，则按相连对拆分为 $$\mathcal{W}(a) \cup \mathcal{W}(b)$$），保留满足最小位点数约束的候选块集合 $$\mathcal{B} = \{B_{\mathcal{K}} : \lvert B_{\mathcal{K}}\rvert \ge m_{\min}\}$$。
>   5. **块内共定位检验（事件 $$\mathcal{G}_1$$ 与 $$\mathcal{G}_2$$）**：若 $$\mathcal{B} \neq \varnothing$$（事件 $$\mathcal{G}_1$$），在各块 $$B \in \mathcal{B}$$ 内直接对未经条件调整的原始边际统计量运行 `coloc.abf`，取 $$\max_{B \in \mathcal{B}} \mathrm{PP4}(B)$$；若 $$\mathcal{B} = \varnothing$$（事件 $$\mathcal{G}_2$$），归约为全区域统计量 $$\mathrm{PP4}(\{1, \dots, M\})$$。

**超参数选取依据与统计量性质.** 在真实 LD（$$\widehat{\mathbf{R}} = \mathbf{R}_k$$）、单因果（$$K_k = 1$$）且收缩因子近似相等（$$a_{kj} \approx a_{kc}$$）的理想基准下，与因果变异 $$c$$ 相关为 $$r_{jc}$$ 的代理变异具有边际均值 $$r_{jc}\lambda_c$$ 及对数 ABF 非中心比 $$(a_{kj}/a_{kc})\, r_{jc}^2$$。(i) 取 $$\rho_c = 0.1$$（$$\sqrt{\rho_c} \approx 0.32$$）使边际均值处于 $$0.32\lambda_c$$ 量级，中等信号下多低于 $$\alpha_{\mathrm{lead}} = 10^{-3}$$ 的临界值 $$3.29$$，更强或受参考误差影响的代理交由 Step 2 条件筛选剔除；(ii) 取 $$\rho_b = 0.4$$（$$\sqrt{\rho_b} \approx 0.63$$）将核心代理纳入局部块（实测中位块大小 $$M_B = 42 \ll M$$），模拟中的三档因果 LD（`low` $$[0, 0.2)$$、`mid` $$[0.32, 0.63)$$、`high` $$[0.63, 0.9)$$）以此两刻度作为分档边界；(iii) $$\alpha_{\mathrm{lead}} = 10^{-3}$$ 仅在第一阶段校准窗口上按 5% 误报率下的检出率选定并全程冻结。因候选块经数据驱动选取并取最大值，输出得分属经误报率校准的排序检验统计量；无真值应用时采用与 traceCB 一致的固定三分类规则（$$\mathrm{PP4} > 0.7$$ 判共定位，对应块的 $$\mathrm{PP3} > 0.7$$ 判独立信号，其余为不确定）。

*Remark 6.* 当变异间存在较大 MAF 差异时（例如默认 eQTL 中 $$a_{kj} = 0.681$$ 对 $$a_{kc} = 0.918$$），$$r_{jc}^2 = 0.4$$ 实际对应约 $$29.7\%$$ 而非 $$40\%$$ 的非中心对数 ABF 比例；同时，算法按参考相关性 $$\hat r_{jl}$$ 划定邻域与聚类，其本身并不构成对真实样本内相关性 $$r_{jc}$$ 的严格确定性上界或下界。

**对照方法与统一命名.** 本研究涵盖八种方法配置：`block-coloc`（正式方法：条件 lead 筛选，$$\alpha_{\mathrm{lead}} = 10^{-3}$$，含无合格块时的事件 $$\mathcal{G}_2$$ 全区域归约）、`block-coloc (no fallback)`（无合格块时返回缺失，用于评估全区域归约分支的作用）、`block-coloc (marginal leads)`（仅行 Step 1 边际选点，$$\alpha_{\mathrm{lead}} = 10^{-5}$$，用于评估逐步条件筛选的作用）、全区域 `coloc.abf`、`coloc.susie`、`coloc.susie (fallback)`（Wallace 2021；无 credible set 时归约为全区域 `coloc.abf`）、`coloc.signals (mask)` 与 `coloc.signals (cond)`（Wallace 2020；`maxhits = 3`，`r2thr = 0.01`，`pthr` 与 $$\alpha_{\mathrm{lead}}$$ 相同）。

---

## 4. Reference-LD mismatch and multi-signal lead detection

外部参考面板相关矩阵写为 $$\widehat{\mathbf{R}} = \mathbf{R}_k + \boldsymbol\Delta$$（1000 Genomes EUR，$$n_{\mathrm{ref}} = 489$$）。双元正态近似下非对角样本相关的渐近标准差为 $$(1 - r_{jl}^2)/\sqrt{n_{\mathrm{ref}}}$$（$$r_{jl} \approx 0$$ 时约 $$0.045$$）；在离散二倍体基因型低频强连锁处，高阶峰度显著放大方差——例如在两位点等位基因频率均为 $$\mathrm{MAF} = 0.05$$、$$r = 0.7$$ 且按两位点 Hardy–Weinberg 平衡独立配对抽样单倍型的基准下，$$\sqrt{n_{\mathrm{ref}}}(\operatorname{atanh}\hat r - \operatorname{atanh} r)$$ 的精确渐近方差为 $$5.955$$（二元正态为 $$1.00$$）。在 12 个真实基因组窗口上，参考 LD 与样本内 LD 的非对角元绝对差 $$\lvert r - \hat r\rvert$$ 中位数为 $$0.027$$–$$0.036$$，第 99 百分位为 $$0.12$$–$$0.17$$（典型值 $$0.15$$）。

**固定参考 LD 下条件统计量的完整分布.** 考虑性状 $$k$$ 含单因果 $$c$$（非中心强度 $$\lambda_c$$，$$s_{kc} = +1$$），非因果变异 $$j \neq c$$ 与 $$c$$ 的真实相关为 $$r_{jc}$$，给定**固定**参考相关为 $$\hat r_{jc}$$（忽略微小岭项 $$\eta$$）。在正态近似模型 (1.9) 下，$$(z_{kc}, z_{kj})^\top$$ 服从二元正态分布，均值 $$(\lambda_c, r_{jc}\lambda_c)^\top$$，协方差矩阵 $$\sigma^2_k \begin{pmatrix} 1 & r_{jc} \\ r_{jc} & 1 \end{pmatrix}$$。经线性变换 $$(1 - \hat r_{jc}^2)^{-1/2}(-\hat r_{jc},\ 1)(z_{kc}, z_{kj})^\top$$，统计量 $$z_{kj \mid c}^{\mathrm{ref}} = (z_{kj} - \hat r_{jc} z_{kc}) / \sqrt{1 - \hat r_{jc}^2}$$ 的无条件分布与给定 $$z_{kc}$$ 的条件分布均精确为正态分布：

$$
z_{kj \mid c}^{\mathrm{ref}} \;\sim\; \mathcal{N}\!\left(\frac{(r_{jc} - \hat r_{jc})\,\lambda_c}{\sqrt{1 - \hat r_{jc}^2}},\; \sigma^2_k\left[1 + \frac{2\hat r_{jc}(\hat r_{jc} - r_{jc})}{1 - \hat r_{jc}^2}\right]\right), \qquad z_{kj \mid c}^{\mathrm{ref}} \mid z_{kc} \;\sim\; \mathcal{N}\!\left(\frac{(r_{jc} - \hat r_{jc})\, z_{kc}}{\sqrt{1 - \hat r_{jc}^2}},\; \sigma^2_k\,\frac{1 - r_{jc}^2}{1 - \hat r_{jc}^2}\right) \,. \tag{4.1}
$$

其中无条件方差分子由二次型 $$\begin{pmatrix} -\hat r_{jc} & 1 \end{pmatrix}\begin{pmatrix} 1 & r_{jc} \\ r_{jc} & 1 \end{pmatrix}\begin{pmatrix} -\hat r_{jc} \\ 1 \end{pmatrix} = (1 - \hat r_{jc}^2) + 2\hat r_{jc}(\hat r_{jc} - r_{jc})$$ 直接导出。在一般条件集合 $$\mathcal{S}$$ 下，令 $$\hat{\mathbf{w}} = \widehat{\mathbf{R}}_{j\mathcal{S}}\widehat{\mathbf{R}}_{\mathcal{S}\mathcal{S},\eta}^{-1}$$、$$\hat v_{j \mid \mathcal{S}} = 1 - \hat{\mathbf{w}}\widehat{\mathbf{R}}_{\mathcal{S}j}$$，对应无条件均值与方差为 $$(\mu_{kj} - \hat{\mathbf{w}}\boldsymbol\mu_{k\mathcal{S}})/\sqrt{\hat v_{j \mid \mathcal{S}}}$$ 与 $$\sigma^2_k(1 - 2\hat{\mathbf{w}}\mathbf{R}_{k,\mathcal{S}j} + \hat{\mathbf{w}}\mathbf{R}_{k,\mathcal{S}\mathcal{S}}\hat{\mathbf{w}}^\top)/\hat v_{j \mid \mathcal{S}}$$；当 $$\mathbf{R}_{k,\mathcal{S}\mathcal{S}}$$ 可逆时，给定 $$\mathbf{z}_{k\mathcal{S}}$$ 的条件残差方差为 $$\sigma^2_k(1 - \mathbf{R}_{k,j\mathcal{S}}\mathbf{R}_{k,\mathcal{S}\mathcal{S}}^{-1}\mathbf{R}_{k,\mathcal{S}j})/\hat v_{j \mid \mathcal{S}}$$。（若对随机参考面板抽样平均，无条件分布呈正态混合而非单一正态。）

由式 (4.1)，均值偏移随 $$\lambda_c$$ 线性放大；当 $$\alpha_{\mathrm{lead}} = 10^{-3}$$（$$z_{\alpha/2} = 3.29$$）、$$\hat r_{jc} \approx 0$$ 且取第 99 百分位参考误差 $$0.15$$ 时，均值本身达到显著性临界值的信号强度参考尺度为

$$
\lambda^* = \frac{3.29}{0.15} \approx 22 \,. \tag{4.2}
$$

当 GWAS 遗传度 $$h^2_1 = 2.0\%$$（$$\lambda_1 \in [14.21, 20.20]$$）逼近 $$\lambda^*$$ 时，`coloc.signals (cond)` 将含参考 LD 近似误差的条件统计量 $$z^{\mathrm{ref}}$$ 平方代入式 (2.1) 指数项，抬高阴性尾部并使匹配 5% 误报率下的检出率下跌 $$6.7$$ 个百分点；而 `block-coloc` 在块内直接采用未经条件调整的原始边际统计量，当参考 LD 变化仅引起匹配块与参考块的边界在共通单点 ABF 与尺度下相对其并集偏移时，单个块的后验差受子集与尺度扰动界控制（[Supplementary](supplement.md) Proposition S5），整体跨面板稳健性由 [Figure 4](fig4.md) 配对实验检验。

**零假设有效检验数与区域极值检出率.** 在零假设 $$\mathbf{z}_k \sim \mathcal{N}(\mathbf{0}, \mathbf{R}_k)$$ 下，区域极值假阳性率满足 Bonferroni 上界 $$M\alpha$$。基于相关矩阵特征值计算有效独立检验数（Li & Ji 2005；12 窗口均值 $$M_{\mathrm{eff}} = 111 \ll M$$）及在 $$\widehat{\mathbf{R}}$$ 上的 Monte Carlo 积分，$$\alpha_{\mathrm{lead}} = 10^{-3}$$ 与 $$10^{-4}$$ 下的假 lead 率分别为 $$13.4\%$$ 与 $$1.7\%$$（`null_both` 实测 $$13.9\%$$–$$14.5\%$$ 与 $$1.4\%$$–$$2.1\%$$）。有真信号时，区域极大值汇聚因果与代理信号，检出率显著高于式 (1.10) 单点理论值：

| 参数设定（代表情形） | 非中心参数 $$\lambda$$ | 单点显著率（$$10^{-3}$$） | 区域检出 $$\pi_{\mathrm{lead}}$$（$$10^{-3}$$） | 双因果检出 $$\Pr(\lvert\mathcal{L}\rvert \ge 2)$$（理论 / 实测） |
|---|---:|---:|---:|---|
| 性状 2，$$h^2_2 = 0.75\%, K_2 = 2$$ | 1.94 | 0.09 | 0.28 | 4.0% / 4.0% |
| 性状 2，$$h^2_2 = 1.5\%, K_2 = 2$$（默认） | 2.75 | 0.29 | 0.50 | 14.9% / 15.9% |
| 性状 2，$$h^2_2 = 3.0\%, K_2 = 2$$ | 3.90 | 0.73 | 0.84 | 56.2% / 55.9% |
| 性状 2，$$h^2_2 = 5.0\%, K_2 = 2$$ | 5.06 | 0.96 | 0.98 | 90.2% / 89.3% |
| 性状 1，$$h^2_1 = 0.3\%, K_1 = 2$$（默认） | 5.48 | 0.99 | 0.99 | 93.9% / 93.6% |
| 性状 1，$$h^2_1 = 2.0\%, K_1 = 1$$ | 20.20 | 1.00 | 1.00 | — |

## References

- Benner C, et al. (2017). Prospects of fine-mapping trait-associated genomic regions by using summary statistics from genome-wide association studies. *Am J Hum Genet* 101: 539–551.
- Bulik-Sullivan BK, et al. (2015). LD Score regression distinguishes confounding from polygenicity in genome-wide association studies. *Nat Genet* 47: 291–295.
- Giambartolomei C, et al. (2014). Bayesian test for colocalisation between pairs of genetic association studies using summary statistics. *PLoS Genet* 10: e1004383.
- Hukku A, et al. (2021). Probabilistic colocalization of genetic variants from complex and molecular traits: promise and limitations. *Am J Hum Genet* 108: 25–35.
- Hukku A, Sampson MG, Luca F, Pique-Regi R, Wen X (2022). Analyzing and reconciling colocalization and transcriptome-wide association studies from the perspective of inferential reproducibility. *Am J Hum Genet* 109: 825–837.
- Jiang W, Xiao J, Cai M (2026). traceCB: Trans-ancestry cell-type-specific eQTLs mapping by integrating scRNA-seq and bulk data. *bioRxiv*, doi:10.64898/2026.06.20.733502.
- Li J, Ji L (2005). Adjusting multiple testing in multilocus analyses using the eigenvalues of a correlation matrix. *Heredity* 95: 221–227.
- Wakefield J (2009). Bayes factors for genome-wide association studies: comparison with P-values. *Genet Epidemiol* 33: 79–86.
- Wallace C (2020). Eliciting priors and relaxing the single causal variant assumption in colocalisation analyses. *PLoS Genet* 16: e1008720.
- Wallace C (2021). A more accurate method for colocalisation analysis allowing for multiple causal variants. *PLoS Genet* 17: e1009440.
- Yang J, et al. (2012). Conditional and joint multiple-SNP analysis of GWAS summary statistics identifies additional variants influencing complex traits. *Nat Genet* 44: 369–375.
- Zhang W, et al. (2024). SharePro: an accurate and efficient genetic colocalization method accounting for multiple causal signals. *Bioinformatics* 40: btae295.
- Zou Y, et al. (2022). Fine-mapping from summary data with the "Sum of Single Effects" model. *PLoS Genet* 18: e1010299.
