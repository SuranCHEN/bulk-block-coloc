# Method: block-coloc

block-coloc 的目标是在“一个性状有多个条件独立信号”时恢复 coloc.abf 的单因果假设：只在确有多个信号时，按 LD 把区域分成若干块，使每块内每个性状至多一个信号，再在块内做 coloc.abf。记号见 [Model](01_bulk_model.md)。代码：`R/methods.R`（`run_methods()`），交接包原版在 `R/block_coloc.R`。

## 1. Algorithm

**Input.** 两个性状的 $$(\hat b_{kj}, \hat s_{kj}, \mathrm{MAF}_{kj})_{j=1}^M$$、样本量 $$N_k$$、LD 矩阵 $$\widehat{\mathbf{R}}$$；lead 阈值 $$\alpha_{\mathrm{lead}}$$，clump 阈值 $$\rho_c = 0.1$$，邻域阈值 $$\rho_b = 0.4$$（均为 $$r^2$$）；coloc 先验 $$(p_1, p_2, p_{12})$$。

**Step 1: marginal leads.** 对每个性状 $$k$$，候选集合 $$\mathcal{A}_k = \{ j : p_{kj} < \alpha_{\mathrm{lead}} \}$$。按 $$p$$ 值从小到大贪心 clump：$$j$$ 入选当且仅当它与所有已选 lead 的 $$\hat r^2 < \rho_c$$。记结果为 $$\mathcal{L}^{\mathrm{marg}}_k$$。

**Step 2: conditional leads.** 按 $$p$$ 值顺序遍历 $$\mathcal{L}^{\mathrm{marg}}_k$$，第一个直接保留；之后的 $$j$$ 仅当条件统计量 $$z_{kj\mid S}$$（[Model](01_bulk_model.md) 式 6.1，$$S$$ 为已保留的 lead）满足 $$2\Phi(-\lvert z_{kj\mid S}\rvert) < \alpha_{\mathrm{lead}}$$ 时保留。记结果为 $$\mathcal{L}_k$$。

**Step 3: gate.** 若

$$
\lvert\mathcal{L}_1\rvert \ge 1,\quad \lvert\mathcal{L}_2\rvert \ge 1,\quad \max\!\left(\lvert\mathcal{L}_1\rvert, \lvert\mathcal{L}_2\rvert\right) \ge 2,
$$

进入 Step 4；否则返回整区域 coloc.abf 的 PP4。两个性状都只有一个信号时 SCV 假设成立，整区域 coloc.abf 即为合适的检验（Wallace 2021 的模拟中，真单因果时整区域 coloc 表现最好）。

**Step 4: blocks.** 对每个 lead $$l \in \mathcal{L}_k$$ 定义邻域 $$\mathcal{W}(l) = \{ j : \hat r_{jl}^2 \ge \rho_b \}$$。构造二部图：性状 1 的 lead $$a$$ 与性状 2 的 lead $$b$$ 相连当且仅当 $$\mathcal{W}(a) \cap \mathcal{W}(b) \neq \varnothing$$。每个连通分量 $$\mathcal{K}$$ 给出一个块 $$B_{\mathcal{K}} = \bigcup_{l \in \mathcal{K}} \mathcal{W}(l)$$；若分量内同一性状的 lead 效应方向不一致，则按边拆成 $$\mathcal{W}(a) \cup \mathcal{W}(b)$$。去掉少于 2 个 SNP 的块。

**Step 5: block-level coloc.** 对每个块在 $$B$$ 上计算 coloc.abf，返回 $$\max_B \mathrm{PP4}(B)$$ 及对应的块（作为定位结果）。若 Step 4 没有得到任何块（两性状的 lead 邻域互不相交），称为 *no valid block*，返回整区域 coloc.abf 的 PP4（2026-10-28 定为正式规则）。理由：此时两性状的信号在 LD 上互不相邻，按块无法检验共享；退回整区域与 Step 3 中 SCV 成立时的处理一致。把 no valid block 记为无结果的版本（“无 fallback”）作为敏感性分析报告（[Results](04_results.md) §4）。

```text
Algorithm 1  block-coloc (conditional lead selection)
input : (b̂, ŝ, MAF, N) for traits 1, 2;  LD R̂;  α_lead;  ρ_c = 0.1;  ρ_b = 0.4
1  for k in {1, 2}:
2      L_marg[k] ← GreedyClump({j : p_kj < α_lead}, by p, R̂, ρ_c)
3      L[k] ← first element of L_marg[k]
4      for j in L_marg[k][2:]:
5          if 2Φ(−|z_{kj|L[k]}|) < α_lead:  L[k] ← L[k] ∪ {j}
6  if |L[1]| ≥ 1 and |L[2]| ≥ 1 and max(|L[1]|, |L[2]|) ≥ 2:
7      blocks ← ConnectedComponents(L[1], L[2], W(l) = {j : r̂²_jl ≥ ρ_b})
8      if blocks ≠ ∅:  return max_B coloc.abf(B).PP4
9      else:           return coloc.abf(all SNPs).PP4          (no valid block)
10 return coloc.abf(all SNPs).PP4
```

## 2. Versions compared

| 版本 | 代码名 | lead 计数 | 切块条件（Step 3） | 其他 |
|---|---|---|---|---|
| block-coloc（条件 lead 选择，本文方法） | `V3_gate_conditional`，no valid block 时退回（`V3_gate_conditional_fb`） | 条件独立信号（Step 2） | 有 | 两性状用同一 $$\widehat{\mathbf{R}}$$；无 fallback 的版本作为敏感性分析 |
| block-coloc（边际 lead 选择） | `V2_gate_marginal` | 边际 clump（Step 1） | 有 | 同上 |
| block-coloc（交接包原版，对照） | `V0_handoff` | 边际 clump；无候选时用 $$p$$ 值最小的 5 个 SNP 补位 | 无：两性状都有 lead 即切块；单 lead 时取 $$\max(\mathrm{PP4}_{\mathrm{block}}, \mathrm{PP4}_{\mathrm{whole}})$$ | — |

相对交接包原版的改动及理由：

1. **去掉补位**：在没有信号的性状上，补位用噪声造出 lead，进而建出无意义的块。
2. **去掉 $$\max(\mathrm{PP4}_{\mathrm{block}}, \mathrm{PP4}_{\mathrm{whole}})$$**：对两个分数取最大会系统性抬高 PP4，在相同判定线下增加误报。
3. **切块条件**：SCV 成立时不切块（见 Step 3）。
4. **条件 lead 选择**：由 [Model](01_bulk_model.md) 式 (3.1)，与强信号 $$c$$ 有 LD 的 SNP $$j$$ 的边际均值为 $$r_{jc}\mu_c$$。若 $$r_{jc}^2 < \rho_c$$，它躲过 clump，但仍可能显著，例如 $$\mu_c = 20$$、$$r_{jc} = 0.3$$ 时 $$\mathbb{E}[z_j] = 6$$。边际计数会把它当作第二个信号；条件统计量 (6.1) 的期望为 0，不会误计。

## 3. Relation to existing methods

- 与 **coloc.signals (cond)** 和 **PWCoCo** 相同之处在于都用条件分析分离信号（[Model](01_bulk_model.md) §6）。不同之处在于，它们对每一对条件独立信号分别在条件化后的整区域统计量上做 coloc，block-coloc 则不改动统计量，只把区域限制在 lead 的 LD 邻域内。
- 与 **coloc.signals (mask)** 相似：块外的 SNP 不参与计算，相当于屏蔽。区别是 mask 屏蔽已找到信号的邻域、在剩余区域找下一个信号，block-coloc 保留邻域、在邻域内计算。
- 只在 Step 2 用到 LD 的条件分析；Step 4 的邻域只依赖 $$\hat r^2$$ 是否越过 $$\rho_b$$，对参考 LD 的小误差不敏感。

## 4. Comparison methods

| 方法 | 实现 |
|---|---|
| coloc.abf, whole region | `coloc::coloc.abf`，全部 $$M$$ 个 SNP |
| coloc.susie | 两性状分别 `susieR::susie_rss(z, R̂, n)`（$$L = 10$$，coverage 0.95，purity 0.5），`coloc::coloc.susie` 对所有 credible set 对计算 PP4，取最大；任一性状无 credible set 时记为无结果 |
| coloc.susie + coloc.abf fallback | 同上，无 credible set 时返回整区域 coloc.abf（Wallace 2021 讨论的混合用法） |
| coloc.signals (mask / cond) | `coloc::coloc.signals`，`maxhits = 3`，`r2thr = 0.01`，`pthr` 与 $$\alpha_{\mathrm{lead}}$$ 相同（默认 $$10^{-6}$$ 在 $$N_2 = 1{,}000$$ 时几乎找不到信号），取所有信号对中最大的 PP4 |

所有方法使用相同先验 $$p_1 = p_2 = 10^{-4}$$、$$p_{12} = 10^{-5}$$（coloc 默认值），相同的 $$\widehat{\mathbf{R}}$$。

## 5. Tuning: lead threshold

$$\alpha_{\mathrm{lead}} \in \{10^{-3}, 10^{-4}, 10^{-5}\}$$（两性状相同），对每个 block-coloc 版本分别在**校准窗口**上选择：在误报率固定为 5% 时取校准检出率最高者，并列取 $$10^{-3}$$。8 月旧 bulk 的 $$(5\times10^{-8}, 10^{-5})$$ 作为敏感性分析报告。当前设定下选出的阈值：block-coloc $$10^{-3}$$（校准检出率 76.5%，$$10^{-4}$$ 与 $$10^{-5}$$ 分别为 74.4% 与 73.3%；不含退回规则时同样选 $$10^{-3}$$），边际 lead 选择 $$10^{-5}$$，交接包原版 $$10^{-4}$$。

## References

见 [Model](01_bulk_model.md#references)。
