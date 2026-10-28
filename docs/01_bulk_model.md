# Model

本页从个体水平模型出发，逐项推导模拟中用到的量：效应的随机模型与遗传度、边际效应、summary statistics 的分布、coloc.abf 的贝叶斯因子和后验，以及多因果变异、条件分析和参考 LD 带来的偏差。推导框架与 LD score regression（Bulik-Sullivan et al. 2015）和 IBSEP / traceCB（Jiang et al. 2026）的 summary-statistics 模型一致；共定位部分沿用 Giambartolomei et al. (2014) 与 Wallace (2020, 2021)。

## 1. Individual-level model

一个 cis 区域含 $$M$$ 个 SNP。两个性状 $$k \in \{1, 2\}$$ 在两组不重叠的样本中测量，样本量 $$N_k$$（本项目 $$N_1 = 20{,}000$$，$$N_2 = 1{,}000$$）：

$$
\mathbf{y}_k = \mathbf{X}_k \boldsymbol\beta_k + \boldsymbol\varepsilon_k, \qquad \boldsymbol\varepsilon_k \sim \mathcal{N}\!\left(\mathbf{0},\ \sigma^2_{k}\mathbf{I}_{N_k}\right), \qquad \boldsymbol\varepsilon_1 \perp \boldsymbol\varepsilon_2 . \tag{1.1}
$$

- $$\mathbf{X}_k$$ 的每列 $$\mathbf{x}_{kj}$$ 在样本内标准化：$$\mathbf{1}^\top\mathbf{x}_{kj} = 0$$，$$\mathbf{x}_{kj}^\top\mathbf{x}_{kj} = N_k$$。于是样本内 LD 矩阵为 $$\mathbf{R}_k = \mathbf{X}_k^\top\mathbf{X}_k / N_k$$，元素 $$r_{k,jl}$$。
- $$\boldsymbol\beta_k$$ 是标准化尺度上的联合效应：$$\beta_{kj}$$ 表示基因型增加一个标准差，表型变化多少个表型标准差。
- 两组样本不重叠且噪声独立。这对应 LDSC / traceCB 中 inflation 矩阵 $$\mathbf{C}$$ 的非对角元为 0（§4）。

## 2. Effect-size model and heritability

用指示变量 $$\gamma_{kj} \in \{0, 1\}$$ 表示 SNP $$j$$ 是否为性状 $$k$$ 的因果变异，$$\mathcal{C}_k = \{ j : \gamma_{kj} = 1 \}$$，$$K_k = \lvert \mathcal{C}_k \rvert$$。给定 $$\mathcal{C}_k$$，

$$
\beta_{kj} = \gamma_{kj}\, s_{kj} \sqrt{h^2_k / K_k}, \qquad s_{kj} \overset{\text{iid}}{\sim} \mathrm{Unif}\{-1, +1\}. \tag{2.1}
$$

**Heritability.** 遗传部分 $$g_{ki} = \mathbf{x}_{ki}^\top\boldsymbol\beta_k$$ 在样本内的方差为

$$
\operatorname{Var}(g_k) = \frac{1}{N_k}\boldsymbol\beta_k^\top\mathbf{X}_k^\top\mathbf{X}_k\boldsymbol\beta_k = \boldsymbol\beta_k^\top\mathbf{R}_k\boldsymbol\beta_k
= \sum_{j \in \mathcal{C}_k}\beta_{kj}^2 + \sum_{j \neq l \in \mathcal{C}_k}\beta_{kj}\beta_{kl}\, r_{k,jl}
= \frac{h^2_k}{K_k}\Big(K_k + \sum_{j \neq l} s_{kj}s_{kl}\, r_{k,jl}\Big). \tag{2.2}
$$

对符号求期望，交叉项为 0，$$\mathbb{E}_s[\operatorname{Var}(g_k)] = h^2_k$$。为使每次模拟的局部遗传度严格等于 $$h^2_k$$，代码把 $$\boldsymbol\beta_k$$ 乘以 $$c_k = \left(1 + K_k^{-1}\sum_{j \neq l} s_{kj}s_{kl} r_{k,jl}\right)^{-1/2}$$；当因果变异之间 $$\lvert r\rvert < 0.2$$ 时 $$c_k \in [0.91, 1.12]$$。取 $$\sigma^2_k = 1 - h^2_k$$，则 $$\operatorname{Var}(y_k) = 1$$，$$h^2_k = \operatorname{Var}(g_k)/\operatorname{Var}(y_k)$$ 即通常的（狭义、局部）遗传度。

**Per-SNP heritability.** 若因果集合在 $$M$$ 个 SNP 中均匀随机抽取，$$\Pr(\gamma_{kj} = 1) = K_k / M =: p_k$$，则

$$
\omega_k := \mathbb{E}\big[\beta_{kj}^2\big] = p_k \cdot \frac{h^2_k}{K_k} = \frac{h^2_k}{M}, \qquad \mathbb{E}\big[\beta_{kj}^2 \mid \gamma_{kj} = 1\big] = \frac{h^2_k}{K_k} = \frac{\omega_k}{p_k}. \tag{2.3}
$$

$$\omega_k$$ 是 LDSC 和 traceCB 中 $$\boldsymbol\Omega$$ 的对角元（区域平均的 per-SNP 遗传度，$$h^2/M$$）；$$h^2_k / K_k$$ 是每个因果变异的遗传度。traceCB 的模拟用 $$(h^2, p_{\mathrm{causal}})$$ 参数化，等价于这里的 $$(h^2_k, K_k = M p_k)$$。本项目 $$M \approx 1{,}000$$、$$K_k \le 2$$，即 $$p_k \le 0.002$$，比 traceCB 的 $$p_{\mathrm{causal}} \in [0.005, 0.02]$$ 更稀疏，对应“区域内只有一两个信号”的共定位情形。

**Two traits.** 两个性状的效应结构由支撑集 $$(\mathcal{C}_1, \mathcal{C}_2)$$ 决定。共定位指 $$\mathcal{C}_1 \cap \mathcal{C}_2 \neq \varnothing$$。符号 $$s_{1j}$$、$$s_{2j}$$ 独立，因此即使 $$j$$ 共享，$$\operatorname{Cov}(\beta_{1j}, \beta_{2j}) = 0$$，即 traceCB 意义下的 co-heritability $$\omega_x = 0$$。共定位关心的是两个性状的因果变异是否为同一个，而不是效应是否同向；coloc 的贝叶斯因子只依赖 $$z^2$$（§5），与符号无关。

## 3. Marginal effects

对 SNP $$j$$ 的一元回归系数（标准化尺度）为

$$
\hat b_{kj} = \frac{\mathbf{x}_{kj}^\top\mathbf{y}_k}{\mathbf{x}_{kj}^\top\mathbf{x}_{kj}}
= \frac{1}{N_k}\mathbf{x}_{kj}^\top\mathbf{X}_k\boldsymbol\beta_k + \frac{1}{N_k}\mathbf{x}_{kj}^\top\boldsymbol\varepsilon_k
= \underbrace{\sum_{l=1}^M r_{k,jl}\,\beta_{kl}}_{b_{kj}} + \underbrace{\frac{1}{N_k}\mathbf{x}_{kj}^\top\boldsymbol\varepsilon_k}_{e_{kj}} . \tag{3.1}
$$

$$b_{kj} = (\mathbf{R}_k\boldsymbol\beta_k)_j$$ 是真实的边际效应：非因果 SNP 通过 LD 继承因果变异的效应。对随机效应 (2.1) 求矩：

$$
\mathbb{E}[b_{kj}] = 0, \qquad
\operatorname{Var}(b_{kj}) = \sum_l r_{k,jl}^2\, \mathbb{E}[\beta_{kl}^2] = \omega_k\, \ell_{kj}, \qquad \ell_{kj} := \sum_{l} r_{k,jl}^2 , \tag{3.2}
$$

其中 $$\ell_{kj}$$ 是 SNP $$j$$ 在区域内的 LD score（不同 SNP 的效应符号独立、均值为 0，故交叉项 $$\mathbb{E}[\beta_{kl}\beta_{km}] = 0$$）。给定因果集合，$$b_{kj} = \sqrt{h^2_k/K_k}\sum_{c \in \mathcal{C}_k} s_{kc}\, r_{k,jc}$$（不计缩放 $$c_k$$）。

## 4. Summary-statistic model

**Estimation error.** 由 (3.1)，$$\mathbf{e}_k = \mathbf{X}_k^\top\boldsymbol\varepsilon_k / N_k$$，给定 $$\mathbf{X}_k$$ 精确地有

$$
\mathbf{e}_k \sim \mathcal{N}\!\left(\mathbf{0},\ \frac{\sigma^2_k}{N_k}\mathbf{R}_k\right), \qquad \operatorname{Cov}(\mathbf{e}_1, \mathbf{e}_2) = \mathbf{0}. \tag{4.1}
$$

估计误差在 SNP 之间的相关就是 LD；两性状之间无相关，因为样本不重叠。若两组样本有 $$N_s$$ 人重叠、噪声相关为 $$\rho_\varepsilon$$，则 $$\operatorname{Cov}(e_{1j}, e_{2l}) \approx N_s\rho_\varepsilon\, r_{jl} / (N_1 N_2)$$，这正是 LDSC / traceCB 中 $$\mathbf{C}$$ 的非对角元所吸收的项。

**Standard error.** 残差平方和 $$\mathrm{RSS}_{kj} = \lVert\mathbf{y}_k\rVert^2 - N_k\hat b_{kj}^2$$，且 $$\lVert\mathbf{y}_k\rVert^2 / N_k \approx \operatorname{Var}(y_k) = 1$$，故

$$
\hat s_{kj}^2 = \frac{\mathrm{RSS}_{kj}}{(N_k - 2)\, N_k} \approx \frac{1 - \hat b_{kj}^2}{N_k}. \tag{4.2}
$$

**z-scores.** 当 $$\hat b_{kj}^2 \ll 1$$（每个 SNP 解释的方差很小）时，$$\hat s_{kj} \approx 1/\sqrt{N_k}$$，$$z_{kj} = \hat b_{kj}/\hat s_{kj} \approx \sqrt{N_k}\,(b_{kj} + e_{kj})$$，由 (3.1) 与 (4.1)：

$$
\mathbf{z}_k \mid \boldsymbol\beta_k \sim \mathcal{N}\!\left(\sqrt{N_k}\,\mathbf{R}_k\boldsymbol\beta_k,\ \sigma^2_k\mathbf{R}_k\right) \approx \mathcal{N}\!\left(\sqrt{N_k}\,\mathbf{R}_k\boldsymbol\beta_k,\ \mathbf{R}_k\right), \tag{4.3}
$$

（$$\sigma^2_k = 1 - h^2_k \approx 1$$）。(4.3) 是 SuSiE-RSS（Zou et al. 2022）的似然，也是条件分析（Yang et al. 2012）的基础。对 (2.1) 再取期望，得到 LDSC 的期望 $$\chi^2$$：

$$
\mathbb{E}\big[z_{kj}^2\big] = 1 + N_k\operatorname{Var}(b_{kj}) = 1 + \frac{N_k h^2_k}{M}\,\ell_{kj}. \tag{4.4}
$$

**Non-centrality of a causal variant.** 对因果变异 $$c$$，若它与其他因果变异不连锁，$$b_{kc} = \beta_{kc}$$。精确地，一元回归的 $$z$$ 与样本相关 $$\hat\rho$$ 满足 $$z = \sqrt{N_k - 2}\,\hat\rho / \sqrt{1 - \hat\rho^2}$$，而 $$\rho^2 = \beta_{kc}^2 = h^2_k / K_k$$，故

$$
\lambda_{kc} := \mathbb{E}\lvert z_{kc}\rvert \approx \sqrt{\frac{N_k\, h^2_k / K_k}{1 - h^2_k / K_k}} \approx \sqrt{N_k\, h^2_k / K_k}. \tag{4.5}
$$

本项目主设置：$$\lambda_1 = 7.7$$（$$K_1 = 1$$）/ 5.5（$$K_1 = 2$$）；$$\lambda_2 = 3.9$$（$$K_2 = 1$$）/ 2.7（$$K_2 = 2$$）。一个因果变异以 $$p < \alpha$$ 被检出的概率为 $$\Phi(\lambda - z_{\alpha/2}) + \Phi(-\lambda - z_{\alpha/2})$$；$$\lambda_2 = 2.7$$、$$\alpha = 10^{-3}$$（$$z_{\alpha/2} = 3.29$$）时约 29%。区域内的 lead 是 $$\lvert z\rvert$$ 最大的 SNP，最大值取自因果变异及其 LD 代理，因此区域内找到 lead 的概率更高（$$\lambda = 2.75$$ 时约 50%，见 [Supplementary Note](09_parameters.md) Proposition 5）。

**Per-allele scale.** coloc 的输入是每个等位基因的效应。标准化效应换算为 $$\beta^{\mathrm{allele}}_{kj} = \beta_{kj}\,\sigma_y / \operatorname{sd}(x_j)$$，$$\operatorname{sd}(x_j) = \sqrt{2 f_j(1 - f_j)}$$（$$f_j$$ 为 MAF），对应的方差 $$V_{kj} = \hat s^{\,2,\mathrm{allele}}_{kj} \approx \sigma_y^2 / \{N_k\, 2 f_j(1-f_j)\}$$。

## 5. Colocalisation under the single-causal-variant assumption

**Hypotheses.** coloc.abf（Giambartolomei et al. 2014）在 SCV 假设 $$K_k \le 1$$ 下枚举五类配置：$$H_0$$：$$\mathcal{C}_1 = \mathcal{C}_2 = \varnothing$$；$$H_1$$ / $$H_2$$：只有一个性状有因果；$$H_3$$：$$\mathcal{C}_1 = \{j\}$$，$$\mathcal{C}_2 = \{l\}$$，$$j \neq l$$；$$H_4$$：$$\mathcal{C}_1 = \mathcal{C}_2 = \{j\}$$。SNP 层面的先验为 $$p_1, p_2, p_{12}$$（本项目 $$10^{-4}, 10^{-4}, 10^{-5}$$）。

**Why a single-SNP statistic suffices under SCV.** 设性状 $$k$$ 只有因果 $$j$$。由 (4.3)，$$\mathbb{E}[\mathbf{z}_k] = \sqrt{N_k}\,\beta_{kj}\,\mathbf{r}_{\cdot j}$$，其余 SNP 在给定 $$z_{kj}$$ 时的条件分布为

$$
\mathbf{z}_{k,-j} \mid z_{kj} \sim \mathcal{N}\!\left(\mathbf{r}_{-j,j}\, z_{kj},\ \mathbf{R}_{-j,-j} - \mathbf{r}_{-j,j}\mathbf{r}_{-j,j}^\top\right),
$$

它不含 $$\beta_{kj}$$，且与零假设下的条件分布相同。因此“$$j$$ 是唯一因果”对“无关联”的整区域似然比等于只用 $$z_{kj}$$ 的似然比，不需要 LD。这就是 coloc.abf 只用边际统计量、不用 LD 的依据。若性状 $$k$$ 另有因果 $$l$$，上述条件均值变为 $$\mathbf{r}_{-j,j}z_{kj} + \sqrt{N_k}\,\beta_{kl}(\mathbf{r}_{-j,l} - \mathbf{r}_{-j,j}r_{jl})$$，含 $$\beta_{kl}$$，单 SNP 似然比不再等于整区域似然比，SCV 下的推断失效（§6）。

**Approximate Bayes factor.** 在“$$j$$ 为因果”下取先验 $$\beta^{\mathrm{allele}}_{kj} \sim \mathcal{N}(0, W)$$，$$\hat\beta^{\mathrm{allele}}_{kj} \mid \beta \sim \mathcal{N}(\beta, V_{kj})$$，则边际分布为 $$\mathcal{N}(0, V_{kj} + W)$$；零假设下为 $$\mathcal{N}(0, V_{kj})$$。两者之比（Wakefield 2009）：

$$
\mathrm{ABF}_{kj} = \sqrt{\frac{V_{kj}}{V_{kj} + W}}\exp\!\left\{\frac{\hat\beta^2}{2}\left(\frac{1}{V_{kj}} - \frac{1}{V_{kj} + W}\right)\right\}
= \sqrt{1 - r_{kj}}\,\exp\!\left(\tfrac{1}{2}r_{kj}z_{kj}^2\right), \qquad r_{kj} = \frac{W}{V_{kj} + W}. \tag{5.1}
$$

连续性状取 $$\sqrt{W} = 0.15\,\sigma_y$$。本模型中因果变异的每等位基因效应为 $$\sqrt{h^2_k/K_k}/\sqrt{2f(1-f)}$$：eQTL（$$h^2_2 = 0.015$$，$$K_2 = 1$$）在 MAF 0.05–0.5 时为 0.17–0.40，GWAS（$$h^2_1 = 0.003$$）为 0.08–0.18，与先验标准差 0.15 同一量级。零假设下 $$z \sim \mathcal{N}(0,1)$$，$$\mathbb{E}[e^{r z^2/2}] = (1 - r)^{-1/2}$$，故 $$\mathbb{E}[\mathrm{ABF}_{kj}] = 1$$。

**Posterior.** 每个配置的似然相对 $$H_0$$ 为相应 ABF 的乘积，乘以先验后对配置求和：

$$
\begin{aligned}
\Pr(H_1 \mid D) &\propto p_1\textstyle\sum_j \mathrm{ABF}_{1j}, \qquad
\Pr(H_2 \mid D) \propto p_2\textstyle\sum_j \mathrm{ABF}_{2j},\\
\Pr(H_3 \mid D) &\propto p_1p_2\textstyle\sum_{j \neq l}\mathrm{ABF}_{1j}\mathrm{ABF}_{2l}
= p_1p_2\Big(\textstyle\sum_j\mathrm{ABF}_{1j}\sum_l\mathrm{ABF}_{2l} - \sum_j\mathrm{ABF}_{1j}\mathrm{ABF}_{2j}\Big),\\
\Pr(H_4 \mid D) &\propto p_{12}\textstyle\sum_j \mathrm{ABF}_{1j}\mathrm{ABF}_{2j}, \qquad \Pr(H_0 \mid D) \propto 1 .
\end{aligned} \tag{5.2}
$$

若共享变异为 $$a$$ 且 $$\mathrm{ABF}_{1a}$$ 主导 $$\sum_j\mathrm{ABF}_{1j}$$，

$$
\frac{\Pr(H_4\mid D)}{\Pr(H_3\mid D)} \approx \frac{p_{12}}{p_1p_2}\cdot\frac{\mathrm{ABF}_{2a}}{\sum_{l\neq a}\mathrm{ABF}_{2l}}, \qquad
\frac{\Pr(H_4\mid D)}{\Pr(H_1\mid D)} \approx \frac{p_{12}}{p_1}\,\mathrm{ABF}_{2a}. \tag{5.3}
$$

两个推论：(i) $$p_{12}/(p_1p_2) = 10^3$$，只要性状 2 在 $$a$$ 上有一点证据，$$H_4$$ 就压过 $$H_3$$；(ii) 性状 2 无信号（gwas_only）时 $$\mathbb{E}[\mathrm{ABF}_{2a}] = 1$$，$$\mathrm{PP4} \approx 0.1/(1 + 0.1 + 0.01) \approx 0.09$$，不会接近 0。

## 6. Multiple causal variants

当 $$K_2 = 2$$，$$\mathcal{C}_1 = \{a\}$$，$$\mathcal{C}_2 = \{a, b\}$$（single_double_shared）时，$$\sum_{l\neq a}\mathrm{ABF}_{2l}$$ 包含 $$\mathrm{ABF}_{2b} \approx \sqrt{1-r}\exp(r\lambda_{2b}^2/2)$$。由 (5.3)，$$\mathrm{PP4}/\mathrm{PP3}$$ 随 $$\exp\{r(z_{2a}^2 - z_{2b}^2)/2\}$$ 变化：特有信号 $$b$$ 越强于共享信号 $$a$$，后验越移向 $$H_3$$（Wallace 2020, Fig 7）。处理多因果的现有方法：

| 类别 | 方法 | 如何得到单信号的统计量 | 对 LD 的依赖 |
|---|---|---|---|
| 条件分析 | coloc.signals (cond)（Wallace 2020）；PWCoCo（Robinson et al. 2022） | 用 (6.1) 依次条件化，得到条件独立信号，对每对信号做 (5.2) | 需要与 (4.3) 匹配的 $$\mathbf{R}$$ |
| 屏蔽 | coloc.signals (mask)（Wallace 2020） | 去掉与已找到信号 $$r^2 > r^2_{\mathrm{thr}}$$ 的 SNP，在剩余 SNP 上做 (5.2) | 只用于确定屏蔽范围；假设因果之间连锁平衡 |
| 精细定位 | eCAVIAR（Hormozdiari et al. 2016）；coloc.susie（Wallace 2021）；SharePro（Zhang et al. 2024）；fastENLOC（Hukku et al. 2021） | 用 (4.3) 的多因果似然（SuSiE、DAP-G 等）分解为单效应，对效应对做 (5.2) | 似然本身含 $$\mathbf{R}$$，对 LD 失配最敏感 |

block-coloc 用 (6.1) 计数条件独立的 lead，再把计算限制在 lead 的 LD 邻域内做 (5.2)，使每块内恢复 SCV（[Method](02_block_coloc_pseudocode.md)）。

**Conditional z-scores.** 记已选 lead 集合 $$S$$。由 (4.3) 的多元正态条件分布，

$$
z_{kj} \mid \mathbf{z}_{kS} \sim \mathcal{N}\!\left(\mu_j + \mathbf{R}_{jS}\mathbf{R}_{SS}^{-1}(\mathbf{z}_{kS} - \boldsymbol\mu_S),\ 1 - \mathbf{R}_{jS}\mathbf{R}_{SS}^{-1}\mathbf{R}_{Sj}\right), \qquad \boldsymbol\mu = \sqrt{N_k}\,\mathbf{R}\boldsymbol\beta_k .
$$

若 $$\mathcal{C}_k \subseteq S$$（$$j$$ 的信号全部来自 $$S$$），则 $$\mu_j = \mathbf{R}_{jS}\mathbf{R}_{SS}^{-1}\boldsymbol\mu_S$$，标准化条件统计量

$$
z_{kj\mid S} = \frac{z_{kj} - \mathbf{R}_{jS}\mathbf{R}_{SS}^{-1}\mathbf{z}_{kS}}{\sqrt{1 - \mathbf{R}_{jS}\mathbf{R}_{SS}^{-1}\mathbf{R}_{Sj}}} \sim \mathcal{N}(0, 1); \tag{6.1}
$$

若 $$j$$ 自身为因果，其均值为 $$\sqrt{N_k}\,\beta_{kj}\sqrt{1 - \mathbf{R}_{jS}\mathbf{R}_{SS}^{-1}\mathbf{R}_{Sj}}$$。例：$$S = \{c\}$$，$$\lambda_c = 20$$，非因果 SNP $$j$$ 与 $$c$$ 的 $$r = 0.3$$，则边际 $$\mathbb{E}[z_j] = 6$$（$$p \approx 2\times10^{-9}$$，且 $$r^2 = 0.09 < 0.1$$ 不会被 clump 去掉），而 $$\mathbb{E}[z_{j\mid c}] = 0$$。实现中在 $$\mathbf{R}_{SS}$$ 对角线加 $$10^{-6}$$ 的 ridge。

## 7. Reference LD

方法实际使用参考人群 LD $$\widehat{\mathbf{R}} = \mathbf{R} + \boldsymbol\Delta$$（1000 Genomes EUR，$$n_{\mathrm{ref}} = 489$$）。样本相关系数的渐近方差为 $$(1 - r^2)^2/n_{\mathrm{ref}}$$，故 $$\operatorname{sd}(\Delta_{jl}) \approx (1 - r_{jl}^2)/\sqrt{n_{\mathrm{ref}}}$$，$$r \approx 0$$ 时约 0.045。$$M > n_{\mathrm{ref}}$$ 时 $$\operatorname{rank}(\widehat{\mathbf{R}}) \le n_{\mathrm{ref}} - 1$$，不满秩。

**Bias of the conditional statistic.** 取 $$S = \{c\}$$，真 LD $$r$$，参考 LD $$\hat r$$，非因果 SNP $$j$$：

$$
\mathbb{E}\big[z_{kj\mid c}\big] = \frac{r\lambda_c - \hat r\lambda_c}{\sqrt{1 - \hat r^2}} = \frac{(r - \hat r)\,\lambda_c}{\sqrt{1 - \hat r^2}} . \tag{7.1}
$$

偏差与 $$\lambda_c$$ 成正比：$$r = 0.90$$、$$\hat r = 0.85$$ 时，$$\lambda_c = 20$$ 得 1.9，$$\lambda_c = 40$$ 得 3.8（超过 $$z_{\alpha/2} = 3.29$$）。因此信号越强，参考 LD 越可能制造假的条件独立信号。(5.1)–(5.2) 不含 $$\mathbf{R}$$，coloc.abf 不受影响；SuSiE-RSS 的似然 (4.3) 直接含 $$\widehat{\mathbf{R}}$$，受影响最大（Zou et al. 2022；Kanai et al. 2022）。Benner et al. (2017) 指出参考人群样本量需随 GWAS 样本量增大；Hukku et al. (2021) 发现 LD 失配主要降低共定位功效。

本项目的 $$\widehat{\mathbf{R}}$$：按 rsID、位置和等位基因把每个窗口的 SNP 对齐到 1000 Genomes Phase 3 EUR（GRCh37）。等位基因互换时 $$\hat r_{jl} \to -\hat r_{jl}$$（计数等位基因从 $$a_1$$ 换成 $$a_0$$ 使 $$x \to 2 - x$$）。12 个窗口 11,850 个 SNP 中 11,794 个（99.5%）匹配；$$\widehat{\mathbf{R}}$$ 与样本内 LD 的非对角元相关为 0.95–0.99，中位绝对差 0.03，与 $$(1-r^2)/\sqrt{489}$$ 一致。

## 8. Summary of model quantities

| 量 | 定义 / 推导 | 本项目取值 | 出处 |
|---|---|---|---|
| $$h^2_k$$ | $$\operatorname{Var}(g_k)/\operatorname{Var}(y_k)$$，(2.2) | $$h^2_1$$ 0.1%–2%（主 0.3%），$$h^2_2$$ 0.75%–5%（主 1.5%） | 8 月 bulk；第二段扫描 |
| $$K_k$$ | $$\lvert\mathcal{C}_k\rvert$$ | 0–2，由场景决定 | Wallace (2021)；SharePro 的 $$(K_C, K_S)$$ |
| $$h^2_k/K_k$$ | 每个因果变异的遗传度，(2.3) | 例如 eQTL 双信号 0.75% | — |
| $$\omega_k = h^2_k/M$$ | 区域平均 per-SNP 遗传度，$$\boldsymbol\Omega$$ 的对角元，(2.3) | $$\approx 1.5\times10^{-5}$$（eQTL） | LDSC；traceCB |
| $$\sigma^2_k$$ | $$1 - h^2_k$$，使 $$\operatorname{Var}(y) = 1$$ | — | (A3) |
| $$\mathbf{R}_k$$、$$\ell_{kj}$$ | 样本内 LD 与 LD score，(3.2) | 由 UKBB 基因型计算 | — |
| $$\widehat{\mathbf{R}}$$ | 参考 LD，§7 | 1000G EUR，489 人 | Zou et al. (2022) |
| $$\lambda_{kc}$$ | 因果变异的非中心参数，(4.5) | 7.7 / 5.5；3.9 / 2.7 | — |
| $$W$$ | ABF 先验方差，(5.1) | $$(0.15\,\sigma_y)^2$$ | Giambartolomei et al. (2014) |
| $$p_1, p_2, p_{12}$$ | SNP 层面先验，(5.2) | $$10^{-4}, 10^{-4}, 10^{-5}$$ | coloc 默认值 |

## References

- Benner C, et al. (2017). Prospects of fine-mapping trait-associated genomic regions by using summary statistics from genome-wide association studies. *Am J Hum Genet* 101(4): 539–551.
- Bulik-Sullivan BK, et al. (2015). LD Score regression distinguishes confounding from polygenicity in genome-wide association studies. *Nat Genet* 47(3): 291–295.
- Giambartolomei C, et al. (2014). Bayesian test for colocalisation between pairs of genetic association studies using summary statistics. *PLoS Genet* 10(5): e1004383.
- Hormozdiari F, et al. (2016). Colocalization of GWAS and eQTL signals detects target genes. *Am J Hum Genet* 99(6): 1245–1260.
- Hukku A, et al. (2021). Probabilistic colocalization of genetic variants from complex and molecular traits: promise and limitations. *Am J Hum Genet* 108(1): 25–35.
- Jiang W, Xiao J, Cai M (2026). traceCB: Trans-ancestry cell-type-specific eQTLs mapping by integrating scRNA-seq and bulk data. *bioRxiv*, doi:10.64898/2026.06.20.733502.
- Kanai M, et al. (2022). Meta-analysis fine-mapping is often miscalibrated at single-variant resolution. *Cell Genomics* 2(12): 100210.
- Robinson JW, et al. (2022). An efficient and robust tool for colocalisation: Pair-wise Conditional and Colocalisation (PWCoCo). *bioRxiv* 2022.08.08.503158.
- Wakefield J (2009). Bayes factors for genome-wide association studies: comparison with P-values. *Genet Epidemiol* 33(1): 79–86.
- Wallace C (2020). Eliciting priors and relaxing the single causal variant assumption in colocalisation analyses. *PLoS Genet* 16(4): e1008720.
- Wallace C (2021). A more accurate method for colocalisation analysis allowing for multiple causal variants. *PLoS Genet* 17(9): e1009440.
- Wang G, et al. (2020). A simple new approach to variable selection in regression, with application to genetic fine mapping. *J R Stat Soc B* 82(5): 1273–1300.
- Yang J, et al. (2012). Conditional and joint multiple-SNP analysis of GWAS summary statistics identifies additional variants influencing complex traits. *Nat Genet* 44(4): 369–375.
- Zhang W, et al. (2024). SharePro: an accurate and efficient genetic colocalization method accounting for multiple causal signals. *Bioinformatics* 40(5).
- Zou Y, et al. (2022). Fine-mapping from summary data with the "Sum of Single Effects" model. *PLoS Genet* 18(7): e1010299.
