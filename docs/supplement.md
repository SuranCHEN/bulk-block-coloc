# Formal Proofs, Sensitivity Analyses, and Extended Ablations

本页给出正文各节引用的形式化命题、引理与分布细节证明（§1），并汇总窗口层面推断、分块规则消融、表型尺度估计、正态效应分布、四类阴性联合校准及归档方案对比（§2）。符号沿用 [Model and Method](method.md)。

---

## 1. Formal propositions, lemmas, and statistical proofs

### 1.1 Fixed-design versus random-design causal $$t$$-statistic distribution

*Remark S1.* 设性状 $$k$$ 仅含单一因果变异 $$c$$（或全体因果相互正交使 $$\xi_k = 1$$）。若 $$(x_{kc}, \tilde y_k)$$ 视为来自真实相关为 $$\rho_{kc} = s_{kc}\sqrt{q_k}$$ 的独立同分布双变量正态总体，样本相关满足 $$\sqrt{N_k}(\hat\rho_{kc} - \rho_{kc}) \xrightarrow{d} \mathcal{N}\big(0, (1 - \rho_{kc}^2)^2\big)$$。对式 (1.7) 变换 $$\psi_{N_k}(\rho) = \sqrt{N_k - 2}\,\rho / \sqrt{1 - \rho^2}$$ 求导得 $$\psi'_{N_k}(\rho) = \sqrt{N_k - 2}\,(1 - \rho^2)^{-3/2}$$，一阶 Delta 方法给出随机设计渐近正态分布

$$
z_{kc} \;\dot\sim\; \mathcal{N}\!\left(\psi_{N_k}(\rho_{kc}),\ \frac{1}{1 - q_k}\right) \,. \tag{S1}
$$

（若存在多个非正交因果变异，一般相关为 $$\rho_{kc} = \xi_k s_{kc}\sqrt{q_k} + \sum_{c' \neq c} r_{k,cc'}\beta_{kc'}$$，且离散基因型下四阶矩修正取决于设计与 MAF。）与之不同，在本研究的固定样本内标准化基因型（$$\mathbf{x}_{kc}^\top \mathbf{x}_{kc} = N_k - 1$$）下，$$z_{kc}$$ 服从式 (1.10) 的精确非中心 $$t_{\nu_k}(s_{kc}\delta_k)$$ 分布（$$\nu_k = N_k - 2$$）。利用大样本展开 $$\sqrt{\nu_k / U} = 1 - \frac{U - \nu_k}{2\nu_k} + O_p(N_k^{-1})$$ 与 $$\operatorname{Var}(U / (2\nu_k)) = (2\nu_k)^{-1}$$，当固定 $$q_k$$ 令 $$N_k \to \infty$$ 时，$$\delta_k^2 / (2\nu_k) \to q_k / \{2(1 - q_k)\}$$，给出固定设计渐近方差

$$
\lim_{N_k \to \infty}\operatorname{Var}(z_{kc} \mid \mathbf{X}_k) = 1 + \frac{q_k}{2(1 - q_k)} < \frac{1}{1 - q_k} \,. \tag{S2}
$$

表 S1 汇总了 $$N = 1{,}000, q = 0.05$$ 下两类模型各自 1,000,000 次独立 Monte Carlo 抽样与精确解析分布及渐近极限的对比：

| 设计模型 | 经验均值（100 万次） | 经验方差（100 万次） | 精确有限样本均值 / 方差 | 对应渐近方差极限 |
|---|---:|---:|---|---:|
| i.i.d. 双变量正态随机设计 | 7.2539 | 1.0545 | 7.2548 / 1.0548 | 1.0526 |
| 固定标准化基因型 OLS 设计 | 7.2579 | 1.0266 | 7.2566 / 1.0284 | 1.0263 |

*Note.* 两种抽样的经验均值与方差相对各自精确有限样本分布的偏差均在 $$1.3$$ 个 Monte Carlo 标准误以内；一阶 Delta 位置近似 $$\psi_{1000}(\sqrt{0.05}) = 7.2475$$ 反映预期的 $$O(N^{-1/2})$$ 有限样本阶数差异。

---

### 1.2 Asymptotic and finite-sample behaviour of whole-region `coloc.abf`

记性状 $$k$$ 的真实边际效应向量为 $$\mathbf{b}_k = \mathbf{R}_k \boldsymbol\beta_k$$，区域最大绝对效应位置为 $$j_k^* = \arg\max_{1 \le j \le M} \lvert b_{kj}\rvert$$。

**Proposition S1（多因果情形下全区域 `coloc.abf` 的渐近极限）.** *在正态近似模型 (1.9)（即 $$z_{kj} = \sqrt{N_k}\,b_{kj} + O_p(1)$$）下，设区域 SNP 数 $$M \ge 2$$ 固定，联合效应 $$\boldsymbol\beta_k$$ 固定且 $$\mathbf{b}_k \neq \mathbf{0}$$，每 SNP 先验 $$(p_1, p_2, p_{12})$$ 为与样本量无关的正常数，各性状最大值点 $$j_1^*, j_2^*$$ 各自唯一，且 $$N_1, N_2 \to \infty$$ 满足 $$N_1 / N_2 \to \kappa \in (0, \infty)$$、$$V_{kj} \asymp 1 / N_k$$、$$W_k$$ 收敛于正常数。则全区域 `coloc.abf` 满足*

$$
\mathrm{PP4} \xrightarrow{p} \begin{cases} 1, & j_1^* = j_2^*, \\ 0, & j_1^* \neq j_2^* . \end{cases} \tag{S3}
$$

*Proof.* 由 $$z_{kj} = \sqrt{N_k}\,b_{kj} + O_p(1)$$ 与 $$V_{kj} \asymp 1 / N_k$$，收缩因子满足 $$a_{kj} = W_k / (W_k + V_{kj}) = 1 - O(1 / N_k)$$ 且 $$\tfrac{1}{2}\log(1 - a_{kj}) = O(\log N_k)$$。代入 [Method](method.md) 式 (2.1) 得 $$\log \mathrm{ABF}_{kj} = \tfrac{1}{2} N_k b_{kj}^2 + O_p(\sqrt{N_k}) + O(\log N_k)$$。在式 (2.2) 中，$$H_4$$ 未归一化权重的对数主项为 $$\tfrac{1}{2}\max_j (N_1 b_{1j}^2 + N_2 b_{2j}^2)$$，而 $$H_3$$ 为 $$\tfrac{1}{2}\max_{i \neq j}(N_1 b_{1i}^2 + N_2 b_{2j}^2)$$。当 $$j_1^* = j_2^*$$ 时，由唯一性知 $$H_4$$ 在 $$j_1^*$$ 处的主项严格高于 $$H_3$$ 至少一个 $$O(N)$$ 正量；当 $$j_1^* \neq j_2^*$$ 时，$$H_3$$ 在 $$(i, j) = (j_1^*, j_2^*)$$ 处严格主导。其余假设 $$H_0, H_1, H_2$$ 至多含单性状线性增长指数项，固定先验仅贡献 $$O(1)$$ 常数项。$$\square$$

*Remark S2.* 在共同样本量 $$N_k$$ 下，由式 (1.6) 有 $$\hat s_{kj}^2 \approx (\operatorname{Var}(y_k) - b_{kj}^2)/N_k$$；因映射 $$b \mapsto b^2 / (\operatorname{Var}(y_k) - b^2)$$ 在 $$\lvert b\rvert < \sqrt{\operatorname{Var}(y_k)}$$ 上关于 $$\lvert b\rvert$$ 严格单调递增，除以残差标准误保持 $$\lvert b_{kj}\rvert$$ 的排序。若外部汇总数据各 SNP 样本量 $$N_{kj}$$ 不等，则主项需保留残差方差分母而成为 $$\tfrac{1}{2} N_{kj} b_{kj}^2 / \{\operatorname{Var}(y_k) - b_{kj}^2\}$$。在本模拟的等幅重标定设计下，未破坏对称性前各因果边缘强度近乎相等，仅因有限样本随机噪声打破并列，故 `single_double_shared` 与 `double_shared1` 中共享变异在含双因果性状内同时最强的比例约为 $$1/2$$ 与 $$1/4$$（实测为 49.0% 与 20%–23%）。

**Lemma S1（五假设后验概率的条件分解恒等式）.** *记 $$S_1 = \sum_{j=1}^M \mathrm{ABF}_{1j}$$、$$S_2 = \sum_{j=1}^M \mathrm{ABF}_{2j}$$、$$S_{12} = \sum_{j=1}^M \mathrm{ABF}_{1j}\mathrm{ABF}_{2j}$$，以性状 1 后验加权均值 $$m_2 = S_{12} / S_1$$，以及性状 1 边缘空信号总后验质量 $$\pi_{0,1} := \mathrm{PP0} + \mathrm{PP2}$$。在式 (2.2) 下恒有*

$$
\mathrm{PP4} = (1 - \pi_{0,1})\,\widetilde{\mathrm{PP4}}, \qquad \widetilde{\mathrm{PP4}} = \frac{(p_{12}/p_1)\, m_2}{1 + (p_{12}/p_1)\, m_2 + p_2(S_2 - m_2)} \,. \tag{S4}
$$

*Proof.* 式 (2.2) 中 $$H_1, H_3, H_4$$ 的未归一化权重之和为 $$p_1 S_1 + p_1 p_2(S_1 S_2 - S_{12}) + p_{12} S_{12} = p_1 S_1\{1 + p_2(S_2 - m_2) + (p_{12}/p_1) m_2\}$$，其在总权重中的占比恰好为 $$1 - \pi_{0,1}$$。将 $$H_4$$ 权重 $$p_{12} S_{12} = p_1 S_1 (p_{12}/p_1) m_2$$ 除以该和并乘以 $$1 - \pi_{0,1}$$ 即得证。$$\square$$

---

### 1.3 Ideal fixed-partition consistency and non-consistency of data-driven blocks

**Proposition S2（固定区域划分下分块后验的渐近性质）.** *设区域先验划分 $$\mathcal{P}$$ 固定。对任一块 $$B \in \mathcal{P}$$ 及满足 $$\mathbf{b}_{k,B} \neq \mathbf{0}$$ 的性状，记唯一块内极大值点为 $$j_k^*(B) = \arg\max_{j \in B}\lvert b_{kj}\rvert$$。在 Proposition S1 渐近条件下：若两性状在 $$B$$ 内均非零，则 $$\mathrm{PP4}(B) \xrightarrow{p} \mathbb{I}\{j_1^*(B) = j_2^*(B)\}$$；若至少一个性状在 $$B$$ 内满足 $$\mathbf{b}_{k,B} = \mathbf{0}$$，则 $$\mathrm{PP4}(B) \xrightarrow{p} 0$$。从而 $$\max_{B \in \mathcal{P}}\mathrm{PP4}(B) \xrightarrow{p} 1$$ 与并集乘积 $$1 - \prod_{B \in \mathcal{P}}\{1 - \mathrm{PP4}(B)\} \xrightarrow{p} 1$$ 均当且仅当存在某块 $$B \in \mathcal{P}$$ 使两性状块内效应均非零且极大值点重合。*

*Proof.* 当 $$\lvert B\rvert \ge 2$$ 且两性状在 $$B$$ 内均非零时，对块 $$B$$ 直接应用 Proposition S1。当单 SNP 块 $$\lvert B\rvert = 1$$（记 $$B = \{j\}$$）时，$$H_3$$ 权重恒为 $$0$$，余下假设中若 $$b_{1j}, b_{2j}$$ 均非零，则 $$p_{12}\mathrm{ABF}_{1j}\mathrm{ABF}_{2j}$$ 在指数主项上严格主导 $$1$$、$$p_1\mathrm{ABF}_{1j}$$ 与 $$p_2\mathrm{ABF}_{2j}$$，故仍有 $$\mathrm{PP4}(B) \xrightarrow{p} 1$$。若任一块 $$B$$（$$\lvert B\rvert \ge 1$$）中存在性状使 $$\mathbf{b}_{k,B} = \mathbf{0}$$，则该性状对一切 $$j \in B$$ 有 $$z_{kj} = O_p(1)$$ 与 $$\mathrm{ABF}_{kj} = O_p(N_k^{-1/2})$$，使 $$H_3$$ 与 $$H_4$$ 权重相对空信号项或单非零性状项依概率趋于 $$0$$，从而 $$\mathrm{PP4}(B) \xrightarrow{p} 0$$。对有限个块取极大或乘积即得等价充要条件。$$\square$$

*Remark S3.* Proposition S2 仅对先验固定划分刻画“块内边际极值是否重合”。欲将该极限提升为真实因果共定位（$$H_4^{\mathrm{R}}$$）的相合估计，还须额外假定真实 $$\mathbf{R}$$ 严格分块对角且每独立块每性状至多一个因果变异。就现行数据驱动的 Algorithm 1 而言，即使已知真实 $$\mathbf{R}$$，在一般情形下仍不具备普遍渐近相合性：

1. **正交基因型与最小块规模的反例**：取 $$\mathbf{R} = \mathbf{I}_M$$ 并考虑 `double_shared1`（性状 1 因果 $$\{A, B\}$$，性状 2 因果 $$\{A, C\}$$）。当 $$N_1, N_2 \to \infty$$ 时全部因果均被 Step 2 检出，但因非对角相关性全零，邻域 $$\mathcal{W}(l) = \{l\}$$ 均仅为单点。Algorithm 1 Step 4 剔除 $$\lvert B\rvert < 2$$ 的候选块后，唯一共享块 $$\{A\}$$ 被丢弃，归约为全区域 `coloc.abf`（事件 $$\mathcal{G}_2$$）。若各因果等强，全区域极限检出率为 $$1/4 < 1$$；若非共享变异 $$B, C$$ 严格强于 $$A$$，由 Proposition S1 有 $$\mathrm{PP4} \xrightarrow{p} 0$$。
2. **固定超参数与有限参考面板**：固定的 lead 阈值 $$\alpha_{\mathrm{lead}} = 10^{-3}$$（不随样本量趋零）、固定岭项 $$\eta = 10^{-6}$$ 以及固定参考面板容量（$$n_{\mathrm{ref}} = 489$$，致条件残余偏移按 $$\sqrt{N_k / n_{\mathrm{ref}}}$$ 放大）均破坏普遍渐近相合性。实现相合需要 $$\alpha_{\mathrm{lead}} \to 0$$、$$\eta \to 0$$、$$\lVert\widehat{\mathbf{R}} - \mathbf{R}_k\rVert = o_p(N_{\max}^{-1/2})$$，并结合严格分块对角且块大小 $$\ge 2$$ 的结构假定。

---

### 1.4 Linked causal variants, out-of-block shift, and LD resolution

**Proposition S3（双连锁因果变异下边际与条件统计量的均值参数）.** *设某性状含两个因果变异 $$a, b$$，样本内相关为 $$r = r_{ab}$$，符号乘积为 $$s = s_a s_b \in \{-1, +1\}$$，基准非中心强度为 $$\lambda = \sqrt{N q}$$（$$q = h^2 / 2$$）。在正态近似模型 (1.9) 下，配合式 (1.4) 精确重标定 $$\xi = (1 + s r)^{-1/2}$$ 及使用真实 $$r$$ 的无岭条件统计量 (3.3)（给出正态近似模型下的均值参数而非有限样本学生化 $$t$$ 统计量的精确期望），有*

$$
\big\lvert\mathbb{E}[z_a]\big\rvert = \lambda\sqrt{1 + s r}, \qquad \big\lvert\mathbb{E}[z_{b \mid a}]\big\rvert = \lambda\sqrt{1 - s r} \,. \tag{S5}
$$

*Proof.* 由 $$\boldsymbol\beta^\top \mathbf{R} \boldsymbol\beta = \xi^2 q(2 + 2 s r) = h^2$$ 得 $$\xi = (1 + s r)^{-1/2}$$。由正态近似模型 (1.9)，$$\mathbb{E}[z_a] = \lambda \xi(s_a + r s_b)$$，取绝对值得 $$\lambda\sqrt{1 + s r}$$。条件统计量 $$z_{b \mid a} = (z_b - r z_a)/\sqrt{1 - r^2}$$ 在式 (1.9) 下的均值为 $$\lambda \xi s_b\sqrt{1 - r^2}$$，代入 $$\xi$$ 即得 $$\lambda\sqrt{1 - s r}$$。当 $$\lvert r\rvert < 0.2$$ 时，$$\xi \in [0.91, 1.12]$$，条件强度因子 $$\sqrt{1 - s r} \in [0.89, 1.10]$$，呈依赖于效应符号一致性的约 $$\pm 10\%$$ 移位（有别于忽略 $$\xi$$ 时 $$\sqrt{1 - r^2}$$ 给出的 $$2\%$$ 衰减）。$$\square$$

**Proposition S4（块外因果变异引起的块内边际均值偏移）.** *设候选块 $$B$$ 不含性状 $$k$$ 的任何真因果变异（$$\mathcal{C}_k \cap B = \varnothing$$）。在正态近似模型 (1.9) 下，记 $$\lambda_{kc} = \sqrt{N_k}\lvert\beta_{kc}\rvert$$，对任意 $$j \in B$$ 有*

$$
\mathbb{E}[z_{kj}] = \sum_{c \in \mathcal{C}_k} s_{kc}\lambda_{kc}\, r_{k,jc}, \qquad \max_{j \in B}\big\lvert\mathbb{E}[z_{kj}]\big\rvert \le \sum_{c \in \mathcal{C}_k}\lambda_{kc}\max_{j \in B}\lvert r_{k,jc}\rvert \,. \tag{S6}
$$

*Proof.* 由正态近似模型 (1.9) 知 $$\mathbb{E}[\mathbf{z}_k] = \sqrt{N_k}\,\mathbf{R}_k \boldsymbol\beta_k$$，因 $$\mathcal{C}_k \cap B = \varnothing$$ 直接展开求和并取三角不等式即得证。$$\square$$

**Lemma S2（五假设模型下 $$\mathrm{PP3} > \mathrm{PP4}$$ 的等价判据）.** *令归一化单点权重为 $$u_j = \mathrm{ABF}_{1j} / S_1$$ 与 $$v_j = \mathrm{ABF}_{2j} / S_2$$。在式 (2.2) 下恒有*

$$
\mathrm{PP3} > \mathrm{PP4} \iff O := \sum_{j=1}^M u_j v_j = \frac{S_{12}}{S_1 S_2} < \frac{1}{1 + p_{12}/(p_1 p_2)} \,, \tag{S7}
$$

*在 `coloc` 默认先验下右端阈值为 $$1 / 1001$$。*

*Proof.* $$H_3$$ 与 $$H_4$$ 未归一化权重之差为 $$p_1 p_2(S_1 S_2 - S_{12}) - p_{12} S_{12} = p_1 p_2 S_1 S_2 - (p_1 p_2 + p_{12}) S_{12}$$，其正负号等价于 $$S_{12}/(S_1 S_2) < p_1 p_2 / (p_1 p_2 + p_{12})$$。$$\square$$

*Remark S4.* 若忽略其余背景 SNP 求和与 ABF 收缩，仅比较单因果构型 $$H_3(A, B)$$（$$\lvert r_{AB}\rvert = r$$）与 $$H_4(A)$$（或 $$H_4(B)$$），期望对数后验比大于零给出定性基准 $$r < r^* = \sqrt{1 - 2\lambda_{\min}^{-2}\log\{p_{12}/(p_1 p_2)\}}$$（$$\lambda_2 = 3.90$$ 时 $$r^* = 0.30$$）。在保留完整 ABF 收缩的理想两 SNP（$$M = 2$$）模型下，$$\lambda_2 = 3.90, r = 0$$ 仅给出 $$\Pr(\mathrm{PP3} > \mathrm{PP4}) = 44.9\%$$；而在 977 SNP 真实窗口上该比例达 95.9%（中位交叉点在 $$r \in [0.6, 0.7]$$）。Lemma S2 阐明了该差异的稀释机制：虽然高度相关的重尾 ABF 分布不能单由零假设均值 $$\mathbb{E}[\mathrm{ABF}_{2j}] = 1$$ 保证分母集中，但实测真实顺式窗口（约 940 个背景及代理 SNP）使 $$S_2$$ 处于 $$10^3$$ 量级，从而将非因果位点归一化权重 $$v_A = \mathrm{ABF}_{2A}/S_2$$ 稀释至 $$1/1001$$ 以下。在补充 LD 扫描（$$h^2_2 = 3.0\%$$，12 窗口）中，`single_distinct` 在第一阶段经验临界值 `PP4 > 0.141` 与固定规则 `PP4 > 0.7` 下的假阳性率由 $$r \in [0.32, 0.45)$$ 时的 5.8% / 0.6% 单调升至 $$[0.55, 0.63)$$ 时的 22.3% / 7.8% 及 $$[0.80, 0.90)$$ 时的 93.3% / 61.7%。

**跨性状中等 LD 情景.** 在多信号场景下人为构造跨性状非共享因果对 $$B$$（性状 1）与 $$C$$（性状 2）处于中等 LD（$$\lvert r_{BC}\rvert \in [0.32, 0.63)$$，中位值 $$0.44$$；其余因果对 $$\lvert r\rvert < 0.2$$；12 窗口 $$\times$$ 50 次重复，沿用第一阶段经验临界值）：

| 场景与设定 | 跨性状 $$B$$–$$C$$ LD | `coloc.abf`（%） | `block-coloc`（%） | 块外条件调整对照（%） | `block-coloc`（$$\mathrm{PP4} > 0.7$$，%） | 分块比例（%） |
|---|---|---:|---:|---:|---:|---:|
| `double_distinct`，$$h^2_1 = 0.3\%$$ | low（$$< 0.2$$） | 3.7 | 4.0 | 3.8 | 0.5 | 1.3 |
| `double_distinct`，$$h^2_1 = 0.3\%$$ | mid（$$0.32$$–$$0.63$$） | 12.8 | **16.3** | 15.8 | **1.8** | 11.5 |
| `double_distinct`，$$h^2_1 = 1.0\%$$ | mid（$$0.32$$–$$0.63$$） | 12.8 | **16.0** | 16.0 | 1.5 | 12.5 |
| `double_shared1`，$$h^2_1 = 0.3\%$$ | low / mid | 50.5 / 52.7 | 59.7 / 59.5 | 59.5 / 60.3 | 31.5 / 28.5 | 37.8 / 42.3 |

*Note.* 第一阶段固定规则 $$\mathrm{PP4} > 0.7$$ 下 `double_distinct` 跨性状 `mid` 档 `coloc.abf` 误报率为 0.8%。就三元相关矩阵半正定性而言，存在共同代理 $$j \in \mathcal{W}(l) \cap \mathcal{W}(l')$$ 的必要条件（符号对齐后）为 $$r_{ll'} \ge 2\rho_b - 1$$（$$\rho_b = 0.4$$ 时必要下界为 $$-0.2$$）；实测中，进入局部块分析的 `double_distinct` 位点有 71.0% 同时包含 $$B$$ 与 $$C$$，导致对块外独立信号作条件调整的对照方案同样无法避免合块干扰，且 `double_shared1` 最佳块覆盖共享变异 $$A$$ 的定位率由 98.0% 降至 78.0%。

---

### 1.5 Subset perturbation bound and empirical quantile bound

**Proposition S5（候选子集变动与表型标准差重估的后验扰动上界）.** *设两个嵌套 SNP 集合满足 $$B \subseteq U \subseteq \{1, \dots, M\}$$，且在各自五假设模型中赋予非零先验质量。*

1. *若对一切 $$j \in B$$，在 $$B$$ 与 $$U$$ 上采用完全相同的单点 $$(\mathrm{ABF}_{1j}, \mathrm{ABF}_{2j})$$ 与每 SNP 先验 $$(p_1, p_2, p_{12})$$，记 $$\Omega_U$$ 为 $$U$$ 上五个假设未归一化权重之和，$$\Omega_B \le \Omega_U$$ 为仅由 $$B$$ 内构型组成的权重之和，则确定性地有*

   $$
   \big\lvert \mathrm{PP4}(B) - \mathrm{PP4}(U) \big\rvert \le 1 - \frac{\Omega_B}{\Omega_U} =: \varepsilon_B \,. \tag{S8}
   $$

2. *进一步设 $$\lvert B\rvert \ge 2$$（以保证所有五个假设未归一化权重严格为正且对数权重有限）。若在 $$B$$ 上重新估计表型标准差 $$\widehat{\mathrm{sdY}}_k$$（从而改变 $$W_k$$），令 $$B$$ 上五个假设的对数未归一化权重相对固定 $$U$$ 尺度时的偏移为 $$\Delta\ell_0, \dots, \Delta\ell_4$$，极差为 $$\delta_B = \max_{0 \le h \le 4}\Delta\ell_h - \min_{0 \le h \le 4}\Delta\ell_h$$，则*

   $$
   \big\lvert \mathrm{PP4}^{\mathrm{est}}(B) - \mathrm{PP4}(U) \big\rvert \le \varepsilon_B + \frac{\delta_B}{4} \,. \tag{S9}
   $$

*Proof.* (1) 记 $$\Omega_U = \Omega_B + R$$（$$R \ge 0$$），$$H_4$$ 未归一化权重满足 $$H_{4,U} = H_{4,B} + R_4$$（$$0 \le R_4 \le R$$）。直接通分得 $$\mathrm{PP4}(U) - \mathrm{PP4}(B) = (R_4 \Omega_B - H_{4,B} R)/(\Omega_B \Omega_U) \in [-(H_{4,B}/\Omega_B)(R/\Omega_U),\ R_4/\Omega_U] \subseteq [-R/\Omega_U,\ R/\Omega_U]$$，即证 (S8)。(2) 在 $$B$$ 上构造函数 $$\ell_h(t) = \ell_h + t\,\Delta\ell_h$$（$$t \in [0, 1]$$），由微积分导数关系，softmax 概率 $$p_h(t) \propto \exp\{\ell_h(t)\}$$ 满足

$$
\frac{\mathrm{d}p_4}{\mathrm{d}t} = p_4(1 - p_4)\left(\Delta\ell_4 - \frac{\sum_{h \neq 4} p_h \Delta\ell_h}{1 - p_4}\right) \,. \tag{S10}
$$

括号内两项均介于 $$[\min_h\Delta\ell_h, \max_h\Delta\ell_h]$$，故差值的绝对值至多为 $$\delta_B$$；结合 $$p_4(1 - p_4) \le 1/4$$ 在 $$t \in [0, 1]$$ 上积分并应用三角不等式即得 (S9)。$$\square$$

**Lemma S3（经验分位数的有限样本扰动界）.** *设在容量为 $$n$$ 的校准阴性样本上，统计量 $$T_{\mathrm{bc}}$$ 与基准统计量 $$X$$ 在事件 $$\mathcal{G}_1$$ 之外逐样本相等，记分块比例 $$\hat\varepsilon = n^{-1}\sum_{i=1}^n \mathbb{I}\{i \in \mathcal{G}_1\}$$。则经验分布函数对一切实数 $$t$$ 满足 $$\lvert\hat F_{\mathrm{bc}}(t) - \hat F_X(t)\rvert \le \hat\varepsilon$$；对任意 $$q \in (\hat\varepsilon, 1 - \hat\varepsilon)$$，经验广义逆满足*

$$
\hat F_X^{-1}(q - \hat\varepsilon) \le \hat F_{\mathrm{bc}}^{-1}(q) \le \hat F_X^{-1}(q + \hat\varepsilon) \,. \tag{S11}
$$

*Proof.* 在 $$\mathcal{G}_1$$ 之外 $$\mathbb{I}\{T_{\mathrm{bc},i} \le t\} = \mathbb{I}\{X_i \le t\}$$，故 $$\lvert\hat F_{\mathrm{bc}}(t) - \hat F_X(t)\rvert \le \hat\varepsilon$$ 对一切 $$t \in \mathbb{R}$$ 成立。令 $$t = \hat F_X^{-1}(q + \hat\varepsilon)$$ 得 $$\hat F_{\mathrm{bc}}(t) \ge \hat F_X(t) - \hat\varepsilon \ge q$$，即得右端上界；对任意 $$t < \hat F_X^{-1}(q - \hat\varepsilon)$$ 有 $$\hat F_{\mathrm{bc}}(t) \le \hat F_X(t) + \hat\varepsilon < q$$，令 $$t \uparrow \hat F_X^{-1}(q - \hat\varepsilon)$$ 得左端下界。$$\square$$

*Remark S5.* Lemma S3 直接限制了经验临界值所处的**经验秩（rank）**区间宽度。阈值数值本身高度接近（例如默认设定下 $$0.145$$ 与 $$0.141$$），还依赖校准阴性经验分布 $$\hat F_X$$ 在 $$q = 0.95$$ 处具有足够的局部概率密度而无明显跳跃。

---

## 2. Extended sensitivity and ablation analyses

### 2.1 Window-level inference accounting for within-window LD clustering

同一基因组窗口内的模拟位点共享真实 LD 矩阵与窗口 SNP 数量 $$M$$。表 S2 对比 `block-coloc` 相对全区域 `coloc.abf` 在 6 个独立测试窗口上的窗口均值自由度 5 的 $$t_5$$ 95% 置信区间、同时对校准窗口（重定 5% 临界值）与测试窗口进行整群重抽样的窗口自举（cluster bootstrap）95% 区间（2,000 次重抽），以及描述性的合并位点层面 Wald 区间（`scripts/window_inference.py`）：

| 模拟设定 | 差值（百分点） | 6 窗口均值 $$t_5$$ 95% CI | 校准—测试窗口 Bootstrap 95% CI | 位点层面 Wald 95% CI（描述性） | 差值 $$> 0$$ 窗口数 |
|---|---:|---|---|---|---:|
| 默认（$$h^2_1 = 0.3\%, h^2_2 = 1.5\%$$） | +3.5 | (+2.2, +4.9) | (+1.6, +4.5) | (+2.7, +4.4) | 6 / 6 |
| $$h^2_2 = 0.75\%$$ | +1.5 | (+0.7, +2.3) | (−1.1, +2.4) | (+0.5, +2.5) | 6 / 6 |
| $$h^2_2 = 3.0\%$$ | +7.1 | (+3.0, +11.2) | (+4.5, +10.2) | (+5.6, +8.6) | 6 / 6 |
| $$h^2_2 = 5.0\%$$ | +13.7 | (+11.0, +16.4) | (+11.5, +15.7) | (+11.7, +15.6) | 6 / 6 |
| $$h^2_1 = 0.1\%$$ | +0.8 | (+0.4, +1.3) | (−1.0, +1.2) | (+0.2, +1.5) | 6 / 6 |
| $$h^2_1 = 1.0\%$$ | +6.2 | (+3.9, +8.4) | (+4.2, +7.9) | (+4.8, +7.6) | 6 / 6 |
| $$h^2_1 = 2.0\%$$ | +6.7 | (+3.8, +9.5) | (+4.4, +8.8) | (+5.2, +8.1) | 6 / 6 |
| 正态效应：默认 / $$h^2_2 = 5.0\%$$ | +1.9 / +5.5 | (+0.8, +3.0) / (+4.3, +6.7) | (+0.8, +2.7) / (+4.0, +6.3) | — | 6 / 6；6 / 6 |

*Note.* 在信号最弱的两组设定（$$h^2_2 = 0.75\%$$ 与 $$h^2_1 = 0.1\%$$），全部 6 个测试窗口差值均为正，但纳入校准窗口临界值重估变异性后的校准—测试窗口 Bootstrap 区间包含零；其余低 LD 设定均严格为正。需说明两点有限窗口性质：(i) 在 $$n = 6$$ 个独立窗口上有放回抽样时，单次漏抽某个极值窗口的概率为 $$(5/6)^6 \approx 33.5\%$$（两端至少漏一端达 $$58.2\%$$），故百分位 bootstrap 区间常略窄于 $$t_5$$ 区间；(ii) 独立的随机单元是窗口而非位点，对于 $$n = 6$$ 个校准窗口，无分布容许界所需次序统计量阶数 $$\lceil (6 + 1)\times 0.95\rceil = 7 > 6$$，因此有限窗口设计不构成对全新染色体区域的无分布（distribution-free）5% 误报保证。

### 2.2 Sensitivity to block-construction rules

在默认（9,600 位点）、$$h^2_2 = 5.0\%$$（4,800 位点）与 $$h^2_1 = 2.0\%$$（4,800 位点）共 19,200 个位点上固定 $$\alpha_{\mathrm{lead}} = 10^{-3}$$，单项改变分块规则并在各自校准窗口重新校准 5% 临界值（测试窗口检出率 / 双性状独立阴性误报率，%）：

| 分块规则消融配置 | 默认（相对基准差） | $$h^2_2 = 5.0\%$$（相对基准差） | $$h^2_1 = 2.0\%$$（相对基准差） | 默认最佳块中位 SNP 数 |
|---|---|---|---|---:|
| **基准规则**（$$\rho_c = 0.1, \rho_b = 0.4$$，含多信号筛选与符号拆分） | 76.8 / 3.8（—） | 97.7 / 6.8（—） | 76.3 / 4.3（—） | 42 |
| 去除边际符号拆分 | 76.8 / 3.8（0.0） | 97.7 / 6.8（0.0） | 76.3 / 4.3（0.0） | 42 |
| 去除多信号筛选 | 76.8 / 3.8（0.0） | 97.7 / 6.8（0.0） | 76.3 / 4.3（0.0） | 46 |
| 邻域阈值 $$\rho_b = 0.2$$ / $$0.3$$ | 77.3 (+0.6) / 77.3 (+0.6) | 97.4 (−0.2) / 97.7 (0.0) | 76.8 (+0.5) / 76.6 (+0.2) | 93 / 59 |
| 邻域阈值 $$\rho_b = 0.5$$ / $$0.6$$ | 76.3 (−0.5) / 75.9 (−0.8) | 97.2 (−0.4) / 97.1 (−0.6) | 76.1 (−0.2) / 75.2 (−1.2) | 34 / 27 |
| 聚类阈值 $$\rho_c = 0.05$$ / $$0.2$$ | 76.8 (0.0) / 76.9 (+0.2) | 97.3 (−0.3) / 97.8 (+0.1) | 76.2 (−0.2) / 76.4 (+0.1) | 42 / 43 |

全部 8 项规则变动的检出率变化幅度均不超过 1.2 个百分点（研究冻结 $$\rho_b = 0.4$$，不作事后调参）。去除多信号筛选条件虽将默认设定阳性分块比例由 29.2% 增至 50.5%，但单共享信号同时主导块内外后验，检出率与误报率均无变化；保留多信号筛选可令单信号位点直接遵循标准 `coloc.abf`。方法学上，按原始边际 $$z$$ 符号拆分连通分量不满足计数等位基因翻转不变性（仅改变符号而保持 LD、$$z^2$$ 与 ABF 不变的等位基因翻转可改变分量拆分）；消融实验表明去除符号拆分在 19,200 个位点中仅影响 8 个位点的连续得分，对检出率与误报率无实质影响（脚本 `scripts/block_ablation.R` 与报告 `results/ablation/ABLATION_*.md` 提供完整数值，基准保持默认参数以与既有运行严格对齐）。

### 2.3 Phenotype scale $$\sigma_{yk}$$ (`sdY`) estimation sensitivity

未传入已知表型标准差 $$\sigma_{yk}$$ 时，`coloc.abf` 在每个候选子集上由 $$(\hat s^2, N, \mathrm{MAF})$$ 线性回归估计 $$\widehat{\mathrm{sdY}}_k$$。在同一 19,200 位点上比较现行块内重估、全区域统一估计一次并固定传入各块、以及固定真实值 $$\sigma_{y1} = \sigma_{y2} = 1$$（`scripts/block_sdy.R`；`results/sdy/SDY_*.md`）：

| $$\sigma_{yk}$$ 估计策略 | 默认：检出率（误报率，%） | $$h^2_2 = 5.0\%$$：检出率（误报率，%） | $$h^2_1 = 2.0\%$$：检出率（误报率，%） | 全位点 $$\lvert\Delta\mathrm{PP4}\rvert$$ 第 99 百分位 |
|---|---|---|---|---:|
| 块内各自重新估计（默认实现） | 76.8 (3.8) | 97.7 (6.8) | 76.3 (4.3) | — |
| 全区域估计一次传入各块 | 76.8 (3.8)，差 0.0 | 97.7 (6.8)，差 0.0 | 76.3 (4.3)，差 0.0 | 0.0011–0.0018 |
| 固定真实总体尺度 $$\sigma_{yk} = 1$$ | 76.7 (3.8)，差 −0.1 | 97.7 (6.8)，差 0.0 | 76.3 (4.3)，差 0.0 | 0.0033–0.0054 |

在所检验的三组网格上，全区域统一估计与块内重估的检出率差精确为 0.0 个百分点（McNemar $$p = 1$$），最佳块变更率 $$\le 0.1\%$$，连续得分第 99 百分位差保持在 $$10^{-3}$$ 量级。以上结果在所检验设定范围内支持块内重估表型尺度并未引入人为功效抬升（因未在各块上显式求值 Proposition S5 的最坏极差 $$\delta_B / 4$$，该经验证据不直接外推至未考察的数据集）。

### 2.4 Normal effect-size distribution and four-class negative calibration

**随机正态效应分布敏感性分析.** 在不作精确 LD 重标定的独立正态随机效应模型 $$\beta_{kj} \sim \mathcal{N}(0, q_k)$$（$$q_k = h^2_k / K_k$$）下，每位点的边缘二阶矩等于名义方差目标 $$\mathbb{E}[\beta_{kj}^2] = q_k$$（与重标定设计下低 LD 双因果诱导二阶矩 $$q_k / (1 - r^2)$$ 仅差 $$O(r^2)$$），但 $$\mathbb{E}\lvert\beta_{kj}\rvert = \sqrt{2/\pi}\sqrt{q_k} \approx 0.80\sqrt{q_k}$$ 且 $$\Pr(\lvert\beta_{kj}\rvert < \tfrac{1}{2}\sqrt{q_k}) = 2\Phi(1/2) - 1 \approx 0.383$$，弱效应位点比例显著增大。在各自校准窗口重新校准临界值后：默认设定下 `block-coloc` 检出率为 39.7%，高于全区域 `coloc.abf` 的 37.8%（+1.9 个百分点；6/6 测试窗口为正，窗口 $$t_5$$ 95% CI: +0.8, +3.0；描述性 McNemar $$p = 2\times 10^{-11}$$），与 `coloc.signals`（39.2%–39.5%）接近；在 $$h^2_2 = 5.0\%$$ 下，`block-coloc` 达 51.7%，高于 `coloc.abf` 的 46.2%（+5.5 个百分点；6/6 测试窗口为正，窗口 $$t_5$$ 95% CI: +4.3, +6.7；描述性 $$p < 10^{-15}$$）。净增益幅度与有效分块分支 $$\mathcal{G}_1$$ 的发生比例（正态默认 11.4%、正态 $$h^2_2 = 5.0\%$$ 22.3%，相应低于重标定默认的 29.2% 与 73.0%）保持一致。

**强信号下 `gwas_only` 的第一类错误膨胀与四类阴性联合校准.** 在 `gwas_only`（单性状强信号、性状 2 全零）下，由 Lemma S1 恒等式 (S4) 知 $$\mathrm{PP4}$$ 中位数稳定在 0.040–0.043（不随 $$h^2_2$$ 下降）。在 $$h^2_2 = 5.0\%$$ 时，仅凭两类双信号独立阴性确定的 5% 临界值降至 0.012–0.014，导致 `coloc.abf` 与 `block-coloc` 在 `gwas_only` 上假阳性率达 100.0%。改用四类阴性各自 95% 分位数的最大值 $$t_* = \max_h \hat F_h^{-1}(0.95)$$ 进行四类阴性联合校准时（测试窗口，%）：

| 设定 | 方法 | 临界值（两类 $$\to$$ 四类） | 测试检出率（变化） | 测试 distinct 误报率 | 测试 `gwas_only` 假阳性率 |
|---|---|---|---:|---:|---:|
| 默认 | `block-coloc` | 0.145 $$\to$$ 0.158 | 76.8 $$\to$$ 75.2 (−1.5) | 3.8 $$\to$$ 3.4 | 4.8 $$\to$$ 4.0 |
| | `coloc.abf` | 0.141 $$\to$$ 0.158 | 73.2 $$\to$$ 71.4 (−1.8) | 3.5 $$\to$$ 3.2 | 5.0 $$\to$$ 4.0 |
| | `coloc.signals (cond)` | 0.210 $$\to$$ 0.210 | 76.2 $$\to$$ 76.2 (0.0) | 3.8 $$\to$$ 3.8 | 2.5 $$\to$$ 2.5 |
| $$h^2_2 = 5.0\%$$ | `block-coloc` | 0.014 $$\to$$ 0.154 | 97.7 $$\to$$ **95.9** (**−1.8**) | 6.8 $$\to$$ 1.3 | 100.0 $$\to$$ **4.3** |
| | `coloc.abf` | 0.012 $$\to$$ 0.154 | 89.8 $$\to$$ **81.9** (**−7.9**) | 5.0 $$\to$$ 0.3 | 100.0 $$\to$$ **4.3** |
| | `coloc.signals (cond)` | 0.152 $$\to$$ 0.182 | 98.2 $$\to$$ 97.8 (−0.3) | 7.5 $$\to$$ 7.0 | 4.3 $$\to$$ 3.7 |
| $$h^2_1 = 2.0\%$$ | `block-coloc` | 0.147 $$\to$$ 0.173 | 76.3 $$\to$$ 74.3 (−2.1) | 4.3 $$\to$$ 4.0 | 5.0 $$\to$$ 4.3 |
| | `coloc.abf` | 0.139 $$\to$$ 0.173 | 70.1 $$\to$$ 67.0 (−3.1) | 3.0 $$\to$$ 2.8 | 5.7 $$\to$$ 4.3 |
| | `coloc.signals (cond)` | 0.270 $$\to$$ 0.290 | 76.3 $$\to$$ 75.1 (−1.2) | 6.5 $$\to$$ 6.3 | 2.0 $$\to$$ 1.7 |

四类阴性联合校准将强 eQTL 下 `gwas_only` 假阳性率压回 4.3%；全区域 `coloc.abf` 因部分非最强共享信号得分滞留于 $$[0.02, 0.15]$$ 损失 7.9 个百分点，而 `block-coloc` 仅下降 1.8 个百分点（相对 `coloc.abf` 优势增至 +14.0 个百分点）。

### 2.5 Archived regional-product formulation and earlier versions

研究亦考察了按互不重叠分块近似合成全局区域后验 $$1 - \prod_b \{1 - \mathrm{PP4}(B_b)\}$$（Hukku et al. 2022；脚本 `scripts/block_partition.R`、`scripts/partition_eval.py`）。在第一阶段测试窗口（5% 匹配误报率）上，含背景余块 `rest` 的乘积后验检出率为 72.5%（较 `block-coloc` 低 4.3 个百分点），不含 `rest` 时为 74.1%（低 2.7 个百分点），仅在联合候选块上求积时与块极大值同为 76.8%。其退化原因在于 $$N_2 = 1{,}000$$ 有限样本下单性状块残留质量（$$\mathrm{ABF}_{2j} \ge \sqrt{1 - a_{2j}} \approx 0.31$$）经多块乘积累积放大，且块间残余 LD 破坏了独立乘积似然的前提，故最终方法选用块极大值统计量。早期探索阶段（固定单点解释方差并输入样本内 LD）存档于 `results/stage1/` 与 `results/stage2/`，均已由本文局部遗传度加外部参考 LD 的主设计取代。
