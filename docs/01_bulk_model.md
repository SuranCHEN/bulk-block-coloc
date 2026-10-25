# 1. bulk 模型

本页说明 bulk 模拟用的统计模型：数据从哪来、表型怎么生成、统计量怎么算，以及为什么这个模型是所有共定位方法的“标准情形”。
对应代码：`R/simulate_bulk.R`。

---

## 1.1 什么是 bulk

**bulk RNA-seq** 把一块组织里成千上万个细胞一起打碎测序，每个人、每个基因只得到**一个表达量**，相当于所有细胞的平均。
GTEx、eQTLGen 这类 eQTL 资源都是 bulk 数据。

在本项目里，“bulk eQTL”指：**每个人一个连续的表达值，用线性模型做 eQTL**。这与 GWAS 的定量性状完全同构，所以两个性状用同一个模型生成。

## 1.2 什么是真实的，什么是模拟的

| 东西 | 真实还是模拟 | 来源 |
|---|---|---|
| 基因型 | **真实** | 本地 UKBB chr22，12 个约 1 Mb 的窗口 |
| LD | **真实** | 由上面的基因型计算 |
| 因果 SNP、效应大小、方向 | 模拟 | 按场景随机抽取 |
| GWAS 表型、eQTL 表达 | 模拟 | 下面的线性模型 |

因为因果 SNP 是我们自己设的，所以每个位点都知道正确答案（共定位或不共定位），才能计算各方法的检出率和误报率。
结果只能解读为“方法在真实 LD 结构下的表现”，不能解读为任何真实基因或疾病的发现。

## 1.3 样本

| 性状 | 人数 | 文件 |
|---|---|---|
| GWAS（性状 1） | 20,000 | `gwas_genotypes.f32`，LD 为 `ld_gwas.f64` |
| eQTL（性状 2） | 1,000 | `qtl_genotypes_n1000.f32`，LD 为 `ld_qtl_n1000.f64` |

两组人**互不重叠**（已按 `donors.tsv` 核对）。这与真实的 GWAS–eQTL 共定位一致：两个研究来自不同人群。
如果两组人重叠，两个性状的 z 值会因为共享样本而相关，coloc 的“两个数据集独立”这一假设就不成立。

## 1.4 表型怎么生成

对一个性状（GWAS 或 eQTL）、第 k 个人：

```
y_k = Σ_j  b_j · x_kj  +  e_k
```

- `x_kj`：第 j 个因果 SNP 的基因型，在该样本内标准化（均值 0、方差 1）；
- `b_j = s_j · sqrt(q)`：效应大小。`s_j` 随机取 +1 或 −1，`q` 是**这个因果 SNP 解释的表型方差比例**；
- `e_k ~ N(0, 1 − Var(g))`：噪声，使总方差 `Var(y) ≈ 1`。

一个性状有两个因果 SNP 时，两者的 q 相同（1:1）。

**强度 q 的取值：**

| 性状 | N | q 的取值 | 主设置 |
|---|---|---|---|
| GWAS | 20,000 | 0.1%、0.3%、1%、2% | 0.3% |
| eQTL | 1,000 | 0.75%、1.5%、3%、5% | 1.5% |

## 1.5 统计量怎么算

对每个 SNP j 单独做一元线性回归（OLS），得到 summary statistics：

```
beta_j = Σ_k x̃_kj y_k / Σ_k x̃_kj²         (x̃ 为中心化后的剂量)
se_j   = sqrt( RSS_j / (N − 2) / Σ_k x̃_kj² )
z_j    = beta_j / se_j
p_j    = 2 · P(T_{N−2} > |z_j|)
```

这正是真实 GWAS 和 bulk eQTL 分析交出来的 summary statistics。

## 1.6 因果 SNP 的 z 值有多大

一元回归的 z 值约等于 `sqrt(N) × 相关系数(x_j, y)`。因果 SNP 与 y 的相关系数平方就是 q，所以

```
因果 SNP 的 |z| ≈ sqrt(N · q)
```

| | q | N | 期望 \|z\| | 对应 p 值 |
|---|---|---|---|---|
| GWAS | 0.1% / 0.3% / 1% / 2% | 20,000 | 4.5 / 7.7 / 14 / 20 | 8e-6 / 9e-15 / 2e-45 / 6e-89 |
| eQTL | 0.75% / 1.5% / 3% / 5% | 1,000 | 2.7 / 3.9 / 5.5 / 7.1 | 6e-3 / 1e-4 / 4e-8 / 2e-12 |

lead 阈值 1e-3 对应 \|z\| ≈ 3.3。eQTL 最弱一档（\|z\| ≈ 2.7）平均过不了阈值，用来看信号很弱时各方法的表现。

强度分几档的目的：信号太强时所有方法都能做对，比不出差别；从勉强能测到很明显分几档，才能看出各方法在什么条件下开始失效。
（8 月的旧 bulk 设置 h² = 0.2、N = 1,000，因果 \|z\| ≈ 14，属于“很强”的一端。）

## 1.7 关键性质：z ~ N(Rλ, R)

这是理解共定位方法的核心。

把所有 SNP 的 z 值写成向量 `z`，LD 矩阵（SNP 之间的相关系数）写成 `R`，那么在效应较小时：

```
z ~ N( R λ ,  R )
```

- **均值 `Rλ`**：`λ` 是真实效应（只在因果 SNP 上非零，大小约 `sqrt(N q)`）。非因果 SNP j 的期望 z 是 `R_jc · λ_c`，也就是**通过 LD 从因果 SNP “借来”的信号**；
- **协方差 `R`**：即使没有任何信号，不同 SNP 的 z 值也按 LD 相关。

为什么成立：每个 z_j 都是同一个 y 与不同 x_j 的相关，x_j 之间的相关就是 LD，所以 z_j 之间也按 LD 相关。

本项目里的每个方法都依赖它：

| 方法 / 步骤 | 怎么用到 |
|---|---|
| coloc.abf | 假设每个性状只有一个因果变异，其余 SNP 的关联都来自它们与这个因果变异的 LD |
| SuSiE / coloc.susie | `susie_rss` 直接以 `z ~ N(Rλ, R)` 为似然 |
| block-coloc（条件 lead 选择）的条件 z | 就是这个多元正态的条件分布，公式见表下 |
| block 划分 | 用 R 判断哪些 SNP 属于同一个信号 |

给定已选的 lead 集合 S，SNP j 的条件 z 为：

```
z_j|S = (z_j − R_jS R_SS⁻¹ z_S) / sqrt(1 − R_jS R_SS⁻¹ R_Sj)
```

**前提：R 必须和产生 z 的样本匹配。** 线性模型下，匹配的 LD 就是**样本内 LD**：GWAS 用 20,000 人的 LD，eQTL 用 1,000 人的 LD。
（已核对：`ld_qtl_n1000.f64` 与 1,000 人剂量的相关矩阵逐元素相同，最大差 4e-14。）
所以 bulk 是这些方法的假设**全部成立的理想情形**：如果某个方法在 bulk 下表现不好，原因在方法本身，而不是模型设错了。

### 例子：边际显著但不是条件独立的 lead

**条件独立信号**（conditionally independent signal）：以已选的 lead 为条件后仍然显著的信号。这是条件与联合分析（GCTA-COJO，Yang et al. 2012）中的标准说法，block-coloc（条件 lead 选择）只把条件独立的 lead 计为信号。

设 GWAS 的因果 SNP c 很强，`z_c = 20`（q = 2%）。旁边有个非因果 SNP j，与 c 的 `r = 0.3`（`r² = 0.09 < 0.1`）：

- j 的期望 z 是 `0.3 × 20 = 6`，p ≈ 2e-9，远过 lead 阈值 1e-3；
- 它与 c 的 `r² = 0.09 < 0.1`，clump 去不掉它；
- 于是只按边际 p 值数 lead（block-coloc 边际 lead 选择），会把 j 当成第二个信号，数出 **2 个 lead**，走 block 分支。

如果以 c 为条件再看 j（block-coloc 条件 lead 选择）：

```
z_j|c = (6 − 0.3 × 20) / sqrt(1 − 0.09) = 0
```

条件后 j 的 z 为 0，说明它不是条件独立的信号，不计为第二个 lead，于是正确地不切块、直接跑整区域 coloc.abf。
信号弱时这个问题不明显：`z_c = 7.7` 时，`r = 0.3` 的 SNP 期望 z 只有 2.3（p ≈ 0.02），过不了 1e-3。
这就是第二段要变化 GWAS 强度的原因之一。

## 1.8 场景

场景名沿用组里 8 月旧 bulk 的命名（single / double 指每个性状的因果变异个数，shared / distinct 指是否共享）；single_double_shared 是新加的，按同一规则命名。
表中 A、B、C、D 表示不同的因果变异。前 4 个场景与 Wallace (2021) Fig 1 的 4 行一一对应，double_distinct 是补充的多信号阴性对照。

| 场景 | GWAS 因果变异 | eQTL 因果变异 | 真值 | Wallace (2021) Fig 1 | 测什么 |
|---|---|---|---|---|---|
| single_shared | A | A | 共定位 | 第 2 行 | 单信号：不应切块，应跑整区域 |
| single_double_shared | A | A, B | 共定位 | 第 3 行 | eQTL 多一个不共享的信号，整区域 coloc.abf 可能被它干扰 |
| double_shared2 | A, B | A, B | 共定位 | 第 4 行 | 两个信号都共享 |
| single_distinct | A | B | 不共定位 | 第 1 行 | 单信号阴性 |
| double_distinct | A, B | C, D | 不共定位 | （补充） | 多信号阴性：检验切块会不会造成误报 |

single_double_shared 中 eQTL 的两个因果变异（共享的 A、eQTL 特有的 B）强度相同，每次模拟由噪声决定谁更显著；评价时按“观测到的 \|z\| 谁更大”拆开报告。

## 1.9 因果 SNP 之间的 LD

用 `ld_gwas` 中的 \|r\| 控制同一场景里那一对因果 SNP（single_double_shared、double_shared2、single_distinct 的 A 与 B；double_distinct 的 A–B 和 C–D 两对；double_distinct 的交叉对始终 \|r\| < 0.2）：

| 档位 | 因果之间 \|r\| | 对应 r² | 方法会怎样 |
|---|---|---|---|
| 低（主设置） | < 0.2 | < 0.04 | 两个因果成为两个独立 lead，各自成块 |
| 中 | 0.32 – 0.63 | 0.1 – 0.4 | 第二个 lead 会被 clump（r² ≥ 0.1）去掉，但又不在第一个 lead 的邻域（r² < 0.4）里 |
| 高 | 0.63 – 0.9 | 0.4 – 0.8 | 两个因果落进同一个块 |

切点来自方法的两个阈值：clump r² = 0.1、邻域 r² = 0.4 换算成 \|r\| 是 0.316 和 0.632，代码里取两位小数 0.32 和 0.63。所以分档测的正是算法行为发生变化的地方。
主设置 \|r\| < 0.2 沿用 8 月的旧 bulk。上限 0.9 是因为 LD 更高的两个 SNP 基本无法区分，所有方法都会失败。
single_shared 没有成对的因果，不参与分档。所有因果 SNP 都要求在两个样本中 MAF ≥ 0.05。

## 1.10 可比性：公共随机数

同一个（窗口、场景、重复）下：

- 两个性状的噪声 `e` 固定；
- 同一 LD 档位下，因果 SNP 和效应方向也固定。

改变强度 q 时只有 q 在变，其余完全相同。所以不同强度之间的差别来自 q 本身，而不是换了一批随机数。

## 1.11 与 8 月旧 bulk 的区别

| | 8 月旧 bulk | 现在 |
|---|---|---|
| 生成代码 | Python（`generate_coloc.py` 等） | R（`simulate_bulk.R`） |
| 窗口 | 每次随机取 2,000 个连续 SNP | 固定的 12 个约 1 Mb 窗口，6 个校准 + 6 个测试 |
| 样本 | 两性状各随机抽 1,000 人 | GWAS 20,000 人，eQTL 1,000 人，固定且不重叠 |
| 强度 | h² = 0.2 一档 | 每个性状 4 档 |
| 场景 | single_shared、single_distinct、double_shared1、double_shared2、double_distinct | single_shared、single_double_shared、double_shared2、single_distinct、double_distinct |
| 因果间 LD | \|r\| < 0.2 | 主设置 \|r\| < 0.2，另有中、高两档 |
| lead 阈值 | GWAS 5e-8，eQTL 1e-5 | 在 1e-3 / 1e-4 / 1e-5 中于校准窗口上选（旧阈值作敏感性分析） |
| 方法 | block-coloc（交接包原版） | block-coloc（条件 lead 选择 / 边际 lead 选择 / 交接包原版）+ 整区域 coloc.abf、coloc.susie、coloc.signals |

## 参考文献

- Giambartolomei C, et al. (2014). Bayesian test for colocalisation between pairs of genetic association studies using summary statistics. *PLoS Genet* 10(5): e1004383.（`doc/material/coloc_2014.pdf`）
- Wallace C (2021). A more accurate method for colocalisation analysis allowing for multiple causal variants. *PLoS Genet* 17(9): e1009440.（`doc/material/coloc_susie_2021.pdf`）
- Yang J, et al. (2012). Conditional and joint multiple-SNP analysis of GWAS summary statistics identifies additional variants influencing complex traits. *Nat Genet* 44(4): 369–375.
- Wang G, Sarkar A, Carbonetto P, Stephens M (2020). A simple new approach to variable selection in regression, with application to genetic fine mapping. *JRSS B* 82(5): 1273–1300.
