# Supplementary Note: derivations for the simulation parameters

本补充说明给出 [Simulation](03_simulation.md) §4 中参数取值所依据的推导。记号与 [Model](01_bulk_model.md) 相同。每个命题先给结论，再给证明；无闭式解的量（S4、S5）用 Monte Carlo 数值计算。理论值与第一段模拟的对照见 [Results](04_results.md) §5。

## S1. Assumptions

- **(A1)** 基因型在样本内标准化，$$\mathbf{R}_k = \mathbf{X}_k^\top\mathbf{X}_k/N_k$$。
- **(A2)** 性状 $$k$$ 有 $$K_k$$ 个因果变异，$$\beta_{kj} = c_k\, s_{kj}\sqrt{q_k}$$，$$q_k = h^2_k/K_k$$，$$s_{kj} \overset{\text{iid}}{\sim} \mathrm{Unif}\{-1, +1\}$$，两性状的符号独立；$$c_k$$ 使 $$\boldsymbol\beta_k^\top\mathbf{R}_k\boldsymbol\beta_k = h^2_k$$（Model 式 2.2）。
- **(A3)** $$\boldsymbol\varepsilon_k \sim \mathcal{N}(\mathbf{0}, (1 - h^2_k)\mathbf{I})$$，两组样本不重叠，$$\boldsymbol\varepsilon_1 \perp \boldsymbol\varepsilon_2$$。
- **(A4)** 每个因果变异解释的方差很小：$$q_k \le 0.05$$。
- **(A5)** 除 S8 外，分析用的 LD 等于真实 LD。

## S2. Signal strength of a causal variant

**Proposition 1.** 在 (A1)–(A4) 下，若因果变异 $$c$$ 与同一性状的其他因果变异不连锁，则

$$
z_{kc} \;\dot\sim\; \mathcal{N}\!\left(\lambda_k,\ \frac{1}{1 - q_k}\right), \qquad \lambda_k = \sqrt{\frac{(N_k - 2)\, q_k}{1 - q_k}} \approx \sqrt{\frac{N_k\, q_k}{1 - q_k}} ,
$$

误差为 $$O(N_k^{-1/2})$$。$$\lambda_k$$ 只依赖乘积 $$N_k q_k$$。

*Proof.* 由 (A2)–(A3)，$$x_c$$ 与 $$y_k$$ 的总体相关为 $$\rho = \beta_{kc}$$，$$\rho^2 = q_k$$（$$c$$ 与其他因果不连锁时 $$c_k = 1$$）。一元回归的 $$t$$ 统计量 $$z = \sqrt{N - 2}\,\hat\rho/\sqrt{1 - \hat\rho^2} =: g(\hat\rho)$$。样本相关满足 $$\hat\rho \,\dot\sim\, \mathcal{N}(\rho, (1 - \rho^2)^2/N)$$，且 $$g'(\rho) = \sqrt{N - 2}\,(1 - \rho^2)^{-3/2}$$。由 delta 方法，$$\mathbb{E}[z] \approx g(\rho) = \lambda_k$$，$$\operatorname{Var}(z) \approx g'(\rho)^2 (1 - \rho^2)^2/N \approx (1 - \rho^2)^{-1}$$。$$\square$$

**Corollary 1.** 因果变异本身在水平 $$\alpha$$ 下显著的概率为 $$\Phi(\lambda_k - z_{\alpha/2}) + \Phi(-\lambda_k - z_{\alpha/2})$$，其中 $$z_{\alpha/2}$$ 在 $$5\times10^{-8}$$、$$10^{-4}$$、$$10^{-3}$$ 时分别为 5.45、3.89、3.29。

| 性状 | $$h^2$$ | $$\lambda$$（$$K$$ = 1 / 2） | $$\Pr(p < 5\times10^{-8})$$ | $$\Pr(p < 10^{-4})$$ | $$\Pr(p < 10^{-3})$$ |
|---|---|---|---|---|---|
| 1（$$N_1 = 20{,}000$$） | 0.001 | 4.47 / 3.16 | 0.16 / 0.01 | 0.72 / 0.23 | 0.88 / 0.45 |
| | **0.003** | **7.76 / 5.48** | 0.99 / 0.51 | 1.00 / 0.94 | 1.00 / 0.99 |
| | 0.01 | 14.2 / 10.0 | 1 / 1 | 1 / 1 | 1 / 1 |
| | 0.02 | 20.2 / 14.2 | 1 / 1 | 1 / 1 | 1 / 1 |
| 2（$$N_2 = 1{,}000$$） | 0.0075 | 2.75 / 1.94 | 0.00 / 0.00 | 0.13 / 0.03 | 0.29 / 0.09 |
| | **0.015** | **3.90 / 2.75** | 0.06 / 0.00 | 0.51 / 0.13 | 0.73 / 0.29 |
| | 0.03 | 5.56 / 3.90 | 0.54 / 0.06 | 0.95 / 0.51 | 0.99 / 0.73 |
| | 0.05 | 7.25 / 5.06 | 0.96 / 0.35 | 1.00 / 0.88 | 1.00 / 0.96 |

粗体为默认值。

**Remark (residual variance).** (A3) 中 $$1 - h^2_k$$ 只通过 $$\lambda_k$$ 中的 $$(1 - q_k)^{-1}$$ 起作用；在 (A4) 下修正不到 3%。

## S3. Effect-size distribution

**Proposition 2.** 设 (A2) 中的固定大小效应 $$\lvert\beta\rvert = \sqrt{q}$$ 换成正态效应 $$\beta \sim \mathcal{N}(0, q)$$（traceCB；Wallace 2021）。则

1. 两者的 $$\mathbb{E}[\beta^2] = q$$ 相同，因此 per-SNP 遗传度 $$\omega = h^2/M$$（Model 式 2.3）和 LDSC 意义下的 $$\mathbb{E}[z^2]$$（式 4.4）相同；
2. 正态效应下 $$\mathbb{E}\lvert\beta\rvert = \sqrt{2/\pi}\,\sqrt{q} \approx 0.80\sqrt{q}$$，$$\Pr(\lvert\beta\rvert < \tfrac{1}{2}\sqrt{q}) = 2\Phi(\tfrac12) - 1 \approx 0.38$$；
3. $$K = 1$$ 时，正态效应下 $$\operatorname{Var}(g) = \beta^2 \sim q\,\chi^2_1$$，局部遗传度的变异系数为 $$\sqrt{2}$$；固定大小效应下为 0。

*Proof.* (1) 两者 $$\beta^2$$ 的期望都是 $$q$$，(4.4) 只依赖 $$\mathbb{E}[\beta^2]$$。(2) 由 $$\beta/\sqrt{q} \sim \mathcal{N}(0, 1)$$ 直接计算。(3) $$\operatorname{Var}(\chi^2_1)/\mathbb{E}[\chi^2_1]^2 = 2$$。$$\square$$

因此固定大小效应下，一个场景内所有位点的信号强度由 $$(N, h^2, K)$$ 唯一确定，场景间的差别只来自因果个数与 LD 结构。正态效应作为敏感性分析单独运行（[Results](04_results.md) §7）：由第 2 点，正态效应下有相当比例的因果变异强度不到名义值的一半，第二个信号更难检出，预期 block-coloc 的增益变小而方向不变，与结果一致。

## S4. Two linked causal variants

**Proposition 3.** 设性状有两个因果变异 $$a, b$$，$$r = r_{ab}$$，$$s = s_a s_b$$，$$\lambda = \sqrt{Nq}$$。在 (A1)–(A5) 下，

$$
c = (1 + s r)^{-1/2}, \qquad
\big\lvert\mathbb{E}[z_a]\big\rvert = \lambda\sqrt{1 + s r}, \qquad
\big\lvert\mathbb{E}[z_{b\mid a}]\big\rvert = \lambda\sqrt{1 - s r},
$$

其中 $$z_{b\mid a}$$ 为 Model 式 (6.1) 的条件统计量。

*Proof.* 由式 (2.2)，$$\boldsymbol\beta^\top\mathbf{R}\boldsymbol\beta = c^2 q (2 + 2sr) = h^2 (1 + sr) c^2$$，令其等于 $$h^2$$ 得 $$c$$。由式 (4.3)，$$\mathbb{E}[z_a] = \lambda c (s_a + r s_b)$$，$$\lvert\cdot\rvert = \lambda c\lvert 1 + sr\rvert = \lambda\sqrt{1 + sr}$$。条件统计量 $$z_{b\mid a} = (z_b - r z_a)/\sqrt{1 - r^2}$$ 的均值为 $$\lambda c\, s_b (1 - r^2)/\sqrt{1 - r^2}$$，取绝对值并代入 $$c$$ 得 $$\lambda\sqrt{(1 - r^2)/(1 + sr)} = \lambda\sqrt{1 - sr}$$。$$\square$$

| $$\lvert r\rvert$$ | $$c$$（同号 / 异号） | $$\sqrt{1 + sr}$$（同 / 异） | $$\sqrt{1 - sr}$$（同 / 异） |
|---|---|---|---|
| 0.2 | 0.91 / 1.12 | 1.10 / 0.89 | 0.89 / 1.10 |
| 0.32 | 0.87 / 1.21 | 1.15 / 0.82 | 0.82 / 1.15 |
| 0.63 | 0.78 / 1.64 | 1.28 / 0.61 | 0.61 / 1.28 |
| 0.9 | 0.73 / 3.16 | 1.38 / 0.32 | 0.32 / 1.38 |

**Corollary 2.** (i) $$\lvert r\rvert < 0.2$$（第一段）时，条件信号与两因果独立时相差不超过 11%。(ii) $$\lvert r\rvert \ge 0.63$$ 时，同号位点的条件信号降到 $$0.32\lambda$$–$$0.61\lambda$$；异号位点的边际信号降到同一范围。档位切点 0.32、0.63 即 $$\sqrt{\rho_c}$$、$$\sqrt{\rho_b}$$。

## S5. Leads

记 Algorithm 1（[Method](02_block_coloc_pseudocode.md)）Step 1–2 选出的 lead 集合为 $$\mathcal{L}$$，$$\mathcal{W}(l) = \{j : r^2_{jl} \ge \rho_b\}$$。

**Proposition 4 (leads without signal).** 若性状无因果变异，则 $$\mathbf{z} \sim \mathcal{N}(\mathbf{0}, \mathbf{R})$$，且

$$
\Pr(\lvert\mathcal{L}\rvert \ge 1) = \Pr\Big(\max_j \lvert z_j\rvert > z_{\alpha/2}\Big) \le M\alpha .
$$

若性状只有一个因果变异 $$c$$，且 $$S \ni c$$，则对 $$j \notin S$$，条件统计量 (6.1) 服从 $$\mathcal{N}(0, 1)$$，与 $$\boldsymbol\beta$$ 无关。因此出现第二个（假）lead 的概率近似为同一区域无信号时的 $$\Pr(\lvert\mathcal{L}\rvert \ge 1)$$；由于 Step 1 已去掉与 $$S$$ 中 lead $$r^2 \ge \rho_c$$ 的 SNP，实际略小，且同样不超过 $$M\alpha$$。

*Proof.* 第一式为定义，上界为 Bonferroni。第二部分：Model §6 中 $$\mathcal{C}_k \subseteq S$$ 时 $$\mu_j = \mathbf{R}_{jS}\mathbf{R}_{SS}^{-1}\boldsymbol\mu_S$$，(6.1) 的均值为 0、方差为 1。$$\square$$

$$\max_j\lvert z_j\rvert$$ 的分布依赖 $$\mathbf{R}$$，无闭式解。用 Li & Ji (2005) 的有效检验数近似，$$\Pr \approx 1 - (1 - \alpha)^{M_{\mathrm{eff}}}$$；12 个窗口的 $$M_{\mathrm{eff}}$$ 为 51–169（均值 111），远小于 $$M$$（546–1,528）。在 $$\widehat{\mathbf{R}}$$ 上直接抽样（每窗口 20,000 次）：

| $$\alpha$$ | $$\Pr(\lvert\mathcal{L}\rvert \ge 1)$$，12 窗口均值（范围） | Monte Carlo 标准误 |
|---|---|---|
| $$10^{-3}$$ | 13.4%（6.1%–20.9%） | 0.07 |
| $$10^{-4}$$ | 1.7%（0.8%–2.7%） | 0.03 |
| $$10^{-5}$$ | 0.20%（0.09%–0.32%） | 0.01 |

**Proposition 5 (leads from true signals).** 无闭式解；以下量由 $$\mathbf{z} \sim \mathcal{N}(\boldsymbol\mu, \widehat{\mathbf{R}})$$ 的 Monte Carlo 计算，因果在 MAF ≥ 0.05 的 SNP 中均匀抽取：

- 单因果，$$\boldsymbol\mu = \lambda\,\widehat{\mathbf{R}}_{\cdot c}$$（每窗口 4,000 次，标准误 ≤ 0.3 个百分点）：$$\pi_{\mathrm{lead}}(\lambda) = \Pr(\max_j\lvert z_j\rvert > z_{\alpha/2})$$，以及 lead 邻域包含因果的概率 $$\pi_{\mathcal{W}}(\lambda) = \Pr(c \in \mathcal{W}(\arg\max_j\lvert z_j\rvert))$$；
- 双因果，$$\lvert r_{ab}\rvert < 0.2$$，$$\boldsymbol\mu$$ 按 Proposition 3（每窗口 500 次，标准误 ≤ 0.6 个百分点）：$$\Pr(\lvert\mathcal{L}\rvert \ge 2)$$ 及“两个因果各被一个不同的 lead 覆盖”的概率。

| $$\lambda$$ | $$\pi_{\mathrm{lead}}$$（$$10^{-3}$$ / $$10^{-4}$$ / $$5\times10^{-8}$$） | $$\pi_{\mathcal{W}}$$，给定找到（$$10^{-3}$$） | 双因果：$$\Pr(\lvert\mathcal{L}\rvert \ge 2)$$（$$10^{-3}$$） | 双因果：两个都被覆盖 |
|---|---|---|---|---|
| 1.94 | 0.28 / 0.07 / 0.00 | 0.46 | 0.04 | 0.01 |
| 2.75 | 0.50 / 0.21 / 0.01 | 0.70 | 0.15 | 0.08 |
| 3.16 | 0.63 / 0.33 / 0.02 | 0.78 | 0.28 | 0.19 |
| 3.90 | 0.84 / 0.60 / 0.08 | 0.89 | 0.56 | 0.46 |
| 4.47 | 0.93 / 0.79 / 0.20 | 0.94 | — | — |
| 5.06 | 0.98 / 0.91 / 0.40 | 0.97 | 0.90 | 0.84 |
| 5.48 | 0.99 / 0.96 / 0.56 | 0.98 | 0.94 | 0.89 |
| 7.25 | 1.00 / 1.00 / 0.97 | 1.00 | — | — |
| 10.0 | 1.00 / 1.00 / 1.00 | 1.00 | 0.99 | 0.97 |
| 14.2 | 1.00 / 1.00 / 1.00 | 1.00 | 0.99 | 0.97 |

两个观察：(i) $$\pi_{\mathrm{lead}}(\lambda)$$ 明显大于 Corollary 1 中因果变异本身显著的概率（$$\lambda = 2.75$$ 时 0.50 对 0.29），因为 lead 取自因果变异及其所有 LD 代理的最大值；(ii) 双因果、$$\lambda = 2.75$$ 时，找到两个 lead 的位点中约一半并非两个因果各对应一个 lead（0.08 / 0.15），另一个 lead 来自噪声。

**Corollary 3 (choice of $$\rho_c$$, $$\rho_b$$).** 与因果变异相关为 $$r$$ 的 SNP 边际均值为 $$r\lambda$$。$$r^2 < \rho_c = 0.1$$ 的 SNP 均值小于 $$0.32\lambda$$，$$\lambda \le 10$$ 时达不到 $$10^{-3}$$ 的期望阈值 3.29，因而 clump 一般不会把同一信号计两次；$$\lambda$$ 更大时会，这时需要条件选择（Step 2）。$$\mathcal{W}(l)$$ 包含所有均值至少为 $$\sqrt{\rho_b}\lambda = 0.63\lambda$$ 的代理。

## S6. Approximate Bayes factors

**Proposition 6.** 在 (A1)–(A4) 下，SNP $$j$$（MAF $$f_j$$）的收缩因子与期望对数 ABF 为

$$
r_{kj} = \frac{W}{W + V_{kj}}, \quad V_{kj} \approx \frac{1}{2 f_j (1 - f_j) N_k}; \qquad
\mathbb{E}\big[\log\mathrm{ABF}_{kj}\big] = \tfrac12\log(1 - r_{kj}) + \tfrac12 r_{kj}\big(\mu_{kj}^2 + 1\big),
$$

$$\mu_{kj} = \mathbb{E}[z_{kj}]$$。

*Proof.* $$V$$ 由 Model §4 的每等位基因尺度换算得到。对式 (5.1) 取对数并用 $$\mathbb{E}[z^2] = \mu^2 + 1$$。$$\square$$

| $$N$$ | $$f = 0.05$$ | 0.1 | 0.2 | 0.3 | 0.5 |
|---|---|---|---|---|---|
| 1,000 | 0.68 | 0.80 | 0.88 | 0.90 | 0.92 |
| 20,000 | 0.98 | 0.99 | 0.99 | 0.995 | 0.996 |

$$\sqrt{W} = 0.15$$。$$r = 0.9$$ 时 $$\mathbb{E}[\log\mathrm{ABF}]$$ 在零假设下为 −0.70，$$\mu = 2.75$$ 时 2.6，$$\mu = 3.9$$ 时 6.1。

**Remark (per-allele effects and MAF).** 在 (A2) 下每等位基因效应为 $$\sqrt{q}/\sqrt{2f(1-f)}$$，即效应大小与 $$[2f(1-f)]^{-1/2}$$ 成正比（Speed et al. 2012 中 $$\alpha = -1$$，GCTA 与 LDSC 的默认）。$$\mathrm{MAF} \ge 0.05$$ 使其不超过 $$f = 0.5$$ 时的 2.3 倍。

## S7. Prior mass of the coloc hypotheses

**Proposition 7.** 区域有 $$M$$ 个 SNP 时，(5.2) 中各假设的先验质量与 $$1 : Mp_1 : Mp_2 : M(M-1)p_1p_2 : Mp_{12}$$ 成比例。特别地，

$$
\frac{\Pr(H_4)}{\Pr(H_3)} = \frac{p_{12}}{(M - 1)\,p_1 p_2} .
$$

*Proof.* 对 (5.2) 取所有 ABF = 1（无数据）即得。$$\square$$

| $$M$$ | $$\Pr(H_0)$$ | $$\Pr(H_1) = \Pr(H_2)$$ | $$\Pr(H_3)$$ | $$\Pr(H_4)$$ | $$H_4 : H_3$$ |
|---|---|---|---|---|---|
| 546 | 0.895 | 0.049 | 0.003 | 0.005 | 1.83 |
| 983 | 0.822 | 0.081 | 0.008 | 0.008 | 1.02 |
| 1,528 | 0.744 | 0.114 | 0.017 | 0.011 | 0.65 |

在 coloc 默认先验下，$$p_{12}/(p_1p_2) = 10^3$$ 与 $$M - 1 \approx 10^3$$ 抵消；同一组 ABF 在 SNP 更多的窗口给出更低的 PP4，因此各方法的判定线在校准窗口上按误报率确定（[Simulation](03_simulation.md) §6）。

**Remark (blocks).** 同一结论适用于 block-coloc 的块：块中位 42 个 SNP 时先验比约 24。把块内 $$p_1, p_2, p_{12}$$ 都乘以 $$s = M/M_B$$ 后，块内 $$H_1, H_2, H_4$$ 的先验质量（$$M_B p_k s$$）与整区域相同，$$H_3$$ 为 $$M_B(M_B-1)p_1p_2s^2 \approx M(M-1)p_1p_2$$。但这样每一对 SNP 的 $$H_3$$ 先验被放大 $$s^2$$ 倍，而块内的 SNP 对多为同一信号的代理，ABF 乘积很大，块内 $$H_3$$ 因而被高估。数值结果见 [Results](04_results.md) §8。

## S8. Reference LD

**Proposition 8.** 在 (A1)–(A4) 下，若分析用 $$\hat r$$ 代替真实 $$r$$，非因果 SNP $$j$$ 对强度为 $$\lambda$$ 的因果 $$c$$ 条件化后，条件统计量的期望越过 $$z_{\alpha/2}$$ 当且仅当

$$
\lvert r_{jc} - \hat r_{jc}\rvert > \frac{z_{\alpha/2}\sqrt{1 - \hat r_{jc}^2}}{\lambda} .
$$

*Proof.* 由 Model 式 (7.1)。$$\square$$

$$n_{\mathrm{ref}} = 489$$ 时 $$\operatorname{sd}(\hat r - r) \approx (1 - r^2)/\sqrt{489}$$（$$r = 0$$ 时 0.045）；12 个窗口中 $$\lvert\hat r - r\rvert$$ 的中位数为 0.027–0.036，99% 分位数为 0.12–0.17。$$\alpha = 10^{-3}$$、$$\hat r \approx 0$$ 时，所需偏差在 $$\lambda = 7.8$$（GWAS 默认）时为 0.42，$$\lambda = 20$$（第二段最强）时为 0.16，后者与 99% 分位数相当。

## S9. Precision of the evaluation

**Proposition 9.** 测试窗口上 $$n$$ 个位点的检出率 $$\hat p$$ 的标准误为 $$\sqrt{p(1-p)/n}$$。两方法在同一批位点上的检出率之差 $$\hat\delta = (b_{10} - b_{01})/n$$ 的标准误为 $$\sqrt{\{\pi_d - \delta^2\}/n} \le \sqrt{\pi_d/n}$$，$$\pi_d = (b_{10} + b_{01})/n$$；在 $$\alpha = 0.05$$、功效 80% 下可检出的差异约为 $$2.80\sqrt{\pi_d/n}$$。

*Proof.* 二项分布与配对多项分布的方差。$$\square$$

| 量 | $$n$$ | 标准误（百分点） |
|---|---|---|
| 单场景检出率，第一段 | 600 | $$p = 0.75$$：1.8；$$p = 0.96$$：0.8 |
| 合并检出率，第一段 | 2,400 | 0.9 |
| 单场景检出率，第二段 | 300 | 2.5 |
| 校准误报率（5% 处） | 1,200 | 0.63 |

合并 2,400 个位点、$$\pi_d$$ = 0.03–0.06 时，可检出的差异为 1.0–1.4 个百分点；单场景（600 个位点）约 2.8 个百分点。

## Code

```bash
python3 scripts/param_checks.py                 # S5：Proposition 4 与单因果部分
python3 scripts/param_checks_two.py 1.94 2.75 3.16 3.9 5.06 5.48 10.03 14.21   # S5：双因果部分
```

只用 1000 Genomes EUR 参考 LD（numpy、pandas）；参考 LD 由 `prepare_ref_ld.R` 在本地生成，不在本仓库中。S2–S4、S6、S7、S9 的数值由文中公式直接计算。

## References

见 [Model](01_bulk_model.md#references)，另加：

- Li J, Ji L (2005). Adjusting multiple testing in multilocus analyses using the eigenvalues of a correlation matrix. *Heredity* 95(3): 221–227.
- Speed D, Hemani G, Johnson MR, Balding DJ (2012). Improved heritability estimation from genome-wide SNPs. *Am J Hum Genet* 91(6): 1011–1021.
