# 2. block-coloc 流程与伪代码

本页写清 block-coloc 三个版本的具体流程：**block-coloc（条件 lead 选择）**（本文方法）、**block-coloc（边际 lead 选择）**，以及作为对照的 **block-coloc（交接包原版）**。
对应代码：交接包原版在 `R/block_coloc.R`（`block_coloc()`），两个新版本和所有对照方法在 `R/methods.R`（`run_methods()`）。
本页的伪代码与代码逐行对应。

## 名称对照

| 本文名称 | 图表中的英文名 | 代码里的名称 |
|---|---|---|
| block-coloc（条件 lead 选择，本文方法） | block-coloc (conditional lead selection, proposed) | `V3_gate_conditional` |
| block-coloc（边际 lead 选择） | block-coloc (marginal lead selection) | `V2_gate_marginal` |
| block-coloc（交接包原版，对照） | block-coloc (original handoff) | `V0_handoff` |

---

## 2.1 为什么要 block-coloc

coloc.abf 假设**每个性状在区域内至多一个因果变异**。如果区域里有多个独立信号（例如 eQTL 有 A、B 两个信号，GWAS 只和 A 共享），
整区域 coloc.abf 会被最强的信号主导，可能把“共享 A”误判成“不共享”（H3）。

block-coloc 的思路：用 LD 把区域切成几块，**让每块里只剩一个信号**，再在每块里跑 coloc.abf，取 PP4 最大的块。
三个版本的区别在于**什么时候该切块、怎么数信号**。

## 2.2 记号和输入

| 记号 | 含义 |
|---|---|
| 性状 1 / 性状 2 | GWAS / eQTL |
| `beta_t, se_t, p_t` | 性状 t 每个 SNP 的效应、标准误、p 值；`z_t = beta_t / se_t` |
| `R_t` | 与性状 t 统计量匹配的 LD 矩阵（bulk：样本内 LD） |
| `thr_t` | lead 阈值，候选 1e-3 / 1e-4 / 1e-5，在校准窗口上选（见 2.11） |
| `clump_r2 = 0.1` | 两个 lead 的 r² 必须低于它 |
| `block_r2 = 0.4` | lead 的邻域：与 lead 的 r² ≥ 0.4 的 SNP |
| `p1 = p2 = 1e-4, p12 = 1e-5` | coloc.abf 的先验，所有方法相同 |

## 2.3 公共子程序

```text
函数 Clump(候选集合 C, p, R):                       # 贪心 clump
    按 p 从小到大排序 C
    leads ← 空
    对 C 中每个 SNP i:
        若 leads 为空 或 对所有 l ∈ leads 都有 R[i,l]² < clump_r2:
            leads ← leads ∪ {i}
    返回 leads

函数 MarginalLeads(性状 t):                          # 两个新版本用；不补位
    C ← { j : p_t[j] < thr_t }
    若 C 为空: 返回 空
    返回 Clump(C, p_t, R_t)

函数 ConditionalLeads(性状 t):                       # 条件 lead 选择用
    M ← MarginalLeads(t)                              # 已按 p 排序
    若 |M| ≤ 1: 返回 M
    kept ← { M 的第一个 }
    对 M 中其余每个 j（按 p 从小到大）:
        z_cond ← (z_t[j] − R_jS R_SS⁻¹ z_t[S]) / sqrt(1 − R_jS R_SS⁻¹ R_Sj)     # S = kept
        若 2·Φ(−|z_cond|) < thr_t:  kept ← kept ∪ {j}
    返回 kept

函数 BuildBlocks(leads1, leads2):                     # 与 R/block_coloc.R 规则相同
    若 leads1 或 leads2 为空: 返回 空
    对每个 GWAS lead a:  W1[a] ← { j : R_1[j,a]² ≥ block_r2 }        # 两个新版本：用各性状自己的 LD；交接包原版两性状都用 R_1
    对每个 eQTL lead b:  W2[b] ← { j : R_2[j,b]² ≥ block_r2 }
    建二部图：a 与 b 之间连边 ⇔ W1[a] ∩ W2[b] 非空
    对图的每个连通分量 K:
        若 K 内同一性状的 lead 的 z 符号不一致（sign_split）:
            K 的每条边 (a,b) 单独成块：W1[a] ∪ W2[b]
        否则:
            整个分量成块：所有 W1[a] ∪ 所有 W2[b]
    去掉少于 2 个 SNP 的块和重复块
    返回 块列表

函数 BlockBranch(leads1, leads2):
    blocks ← BuildBlocks(leads1, leads2)
    若至少一块 ABF 成功:  返回 各块 PP4 的最大值              # branch = block
    否则:                 返回 无结果（no_valid_block）

函数 ABF(SNP 集合 S):
    在 S 上跑 coloc.abf（两性状的 beta、varbeta、MAF、N，type = quant）
    返回 PP4
```

## 2.4 block-coloc（交接包原版，对照）

```text
对每个性状 t:
    C ← { j : p_t[j] < thr_t }
    若 C 为空:            C ← p 值最小的 5 个 SNP          # top-5 补位
    否则若 |C| < 2:       C ← C ∪ p 值最小的 2 个 SNP      # 至少 2 个候选
    leads_t ← Clump(C, p_t, R_1)                            # 两性状都用 R_1
blocks ← BuildBlocks(leads_1, leads_2)                      # 邻域也都用 R_1
block_pp4 ← 各块 ABF 的最大值（没有有效块时为 0）
whole_pp4 ← ABF(全部 SNP)
若 (|leads_1| ≤ 1 或 |leads_2| ≤ 1) 且 whole_pp4 > block_pp4:
    返回 whole_pp4                                          # 单 lead 时取二者较大者
否则若有有效块:
    返回 block_pp4
否则:
    返回 无结果（no_valid_block）
```

## 2.5 block-coloc（边际 lead 选择）

```text
leads_1 ← MarginalLeads(1);  leads_2 ← MarginalLeads(2)
若 |leads_1| ≥ 1 且 |leads_2| ≥ 1 且 max(|leads_1|, |leads_2|) ≥ 2:     # 切块条件
    返回 BlockBranch(leads_1, leads_2)
否则:
    返回 ABF(全部 SNP)                                       # branch = whole
```

## 2.6 block-coloc（条件 lead 选择，本文方法）

```text
leads_1 ← ConditionalLeads(1);  leads_2 ← ConditionalLeads(2)
若 |leads_1| ≥ 1 且 |leads_2| ≥ 1 且 max(|leads_1|, |leads_2|) ≥ 2:
    返回 BlockBranch(leads_1, leads_2)                       # 只用条件后保留的 lead 建块
否则:
    返回 ABF(全部 SNP)
```

和 block-coloc（边际 lead 选择）只差一处：**第 2 个及以后的 lead，要在“以已选 lead 为条件”后仍然显著才计数**。

## 2.7 改了什么、为什么改

| 版本 | 相对上一版的改动 | 为什么 |
|---|---|---|
| block-coloc（交接包原版） | 交接包原样 | 对照 |
| block-coloc（边际 lead 选择） | ① 去掉 top-5 补位；② 去掉“单 lead 时取 max(block, 整区域)”；③ 两性状各用自己的 LD；④ 加切块条件：至少一个性状有 ≥ 2 个 lead 才切块，否则跑整区域 coloc.abf | ① 补位会在没有信号的性状上强行造出 lead，用噪声建块。② 取两个分数的最大值会系统性抬高 PP4，误报随之上升。③ eQTL 统计量应配 eQTL 样本的 LD。④ 两个性状都只有一个信号时，coloc.abf 的单因果假设成立，没有理由切块；Wallace (2021) 也发现真单因果时整区域 coloc 表现最好 |
| block-coloc（条件 lead 选择） | 第 2 个及以后的 lead 必须在条件分析后仍显著才计数 | 与强信号有 LD 的 SNP 也可能边际显著，又因 r² < 0.1 躲过 clump；它不是条件独立的信号，边际计数却会把它误数成第二个信号（例子见 [bulk 模型 1.7](01_bulk_model.md)） |

**找不到足够 lead 时**：两个新版本都直接跑整区域 coloc.abf；交接包原版则先用 top-5 补位凑出候选。

## 2.8 开发中试过、不再报告的版本

为了记录完整，代码里还有两个版本，第一段主分析中计算过，但从 2026-10-25 起不再报告，第二段也不再计算：

- `V1_block_any`（没有切块条件）：去掉补位和取最大，但只要两个性状都有 lead 就切块。它相当于 block-coloc（边际 lead 选择）去掉第 ④ 条，在第一段里表现最差。
- `V4_cond_block`（块内再条件化）：切块条件同 block-coloc（条件 lead 选择），块内再扣除块外 lead 的影响。在第一段里与 block-coloc（条件 lead 选择）几乎相同，没有带来额外好处。

## 2.9 每个方法输出什么

| 字段 | 含义 |
|---|---|
| `branch` | `block` / `whole`（block-coloc）；`susie` / `whole`（coloc.susie）；`signals`（coloc.signals） |
| `pp4` | 最终的共定位后验概率；无结果时为 NA |
| `pp4_whole` | 整区域 coloc.abf 的 PP4（所有方法都存，便于比较） |
| `n_leads1`, `n_leads2` | 用于判断是否切块的 lead 数（条件 lead 选择时为条件后保留的个数） |
| `n_leads1_marginal`, `n_leads2_marginal` | 边际 clump 后的 lead 数 |
| `n_blocks`, `best_block_snps` | 块数、PP4 最大的块包含哪些 SNP |
| `status` | `ok` / `no_valid_block` / SuSiE 的各种无 credible set 情况等 |

评价时另有 “+ whole-region fallback” 版本：block 分支建不出有效块时，用 `pp4_whole` 代替“无结果”。两种处理方式都报告。

## 2.10 对照方法

| 名称 | 做法 |
|---|---|
| coloc.abf, whole region | 整区域 coloc.abf |
| coloc.susie | 两个性状各跑 `susie_rss(z, R, n)`（L = 10，覆盖率 0.95，纯度 0.5），再用 `coloc.susie` 两两配对，取最大 PP4；任一性状没有 credible set 时记为无结果 |
| coloc.susie + coloc.abf fallback | 同上，但没有 credible set 时退回整区域 coloc.abf（Wallace 2021 中表现最好的混合用法） |
| coloc.signals (mask) | 找到一个信号后屏蔽其附近 SNP，再找下一个（最多 3 个，pthr = 1e-3） |
| coloc.signals (cond) | 逐步条件分析拆信号 |

## 2.11 lead 阈值怎么定

- 三个 block-coloc 版本的 lead 阈值在 1e-3、1e-4、1e-5（两性状相同）中选，**只在 6 个校准窗口上选**：在误报率固定为 5% 时，取校准检出率最高的；并列时取 1e-3。第一段选出的结果：条件 lead 选择 1e-4，边际 lead 选择和交接包原版 1e-5。
- 8 月旧 bulk 用的 GWAS 5e-8、eQTL 1e-5 作为敏感性分析单独报告，不参与选择。

## 参考文献

- Wallace C (2021). A more accurate method for colocalisation analysis allowing for multiple causal variants. *PLoS Genet* 17(9): e1009440.
- Wallace C (2020). Eliciting priors and relaxing the single causal variant assumption in colocalisation analyses. *PLoS Genet* 16(4): e1008720.（coloc.signals）
