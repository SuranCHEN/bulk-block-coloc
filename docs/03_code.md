# 3. 代码说明与复现

本页说明本仓库的代码结构、每个文件做什么，以及如何复现。

> UKBB 个体级数据受数据使用协议约束，**不在本仓库中**。代码可以在自己有权使用的基因型数据上运行；不需要任何数据也能用合成数据自检。

---

## 3.1 文件一览

| 文件 | 作用 |
|---|---|
| `R/block_coloc.R` | block-coloc（交接包原版）的实现，以及 clump、二部图连通分量等公共函数 |
| `R/methods.R` | 所有方法的定义：block-coloc（条件 lead 选择 / 边际 lead 选择 / 交接包原版）、整区域 coloc.abf、coloc.susie、coloc.signals。输入只是 z 值和 LD |
| `R/simulate_bulk.R` | bulk 模拟函数库：读窗口、抽因果 SNP、线性模型生成表型、OLS 统计量、组装位点 |
| `scripts/run.R` | 在各窗口上模拟并对每个位点跑所有方法，结果按分片保存，可断点续跑 |
| `scripts/evaluate.R` | 读 `run.R` 的结果，定判定线、出表出图、写 RESULTS |
| `scripts/test_synthetic.R` | 用随机生成的假基因型自检 11 项，不需要真实数据 |
| `results/stage1/` | 第一段（主设置）的汇总结果：`RESULTS.md`、汇总表、图 |
| `docs/` | 本网站的文档 |

代码里的方法内部名（如 `V3_gate_conditional`）只是历史编号，报告和图表一律用描述性名称，对照表见 [第 2 页](02_block_coloc_pseudocode.md)。

## 3.2 数据流

```text
<数据目录>/<窗口>/
  gwas_genotypes.f32, qtl_genotypes_n1000.f32, ld_gwas.f64, ld_qtl_n1000.f64, variants.tsv
        │  bulk_load_bank()            读入并中心化；核对 eQTL LD 与样本内相关一致
        ▼
  bulk_sample_causal()                 按场景和 LD 档位抽因果 SNP 与效应方向
  bulk_draws()                         固定两性状的噪声（公共随机数）
        │
  bulk_gwas() / bulk_eqtl()            线性模型生成 y，逐 SNP OLS → beta, se, z, p
        │
  bulk_assemble_locus()                组成方法需要的位点（各性状自己的 MAF 和 LD）
        │
  run_methods()  (R/methods.R)         block-coloc 各版本、coloc.abf、coloc.susie、coloc.signals
        ▼
outputs/bulk_gated_v1/shards_<stage>/  每个（窗口、场景、重复）一个 .rds
        │  run.R 结尾汇总
        ▼
scores_<stage>.tsv.gz, loci_<stage>.tsv.gz
        │  evaluate.R
        ▼
RESULTS.md, main_table.tsv, pooled_table.tsv, gate_table.tsv, thresholds.tsv, *.png
```

### 数据目录格式

每个窗口一个子文件夹（以窗口中心位置命名），包含：

| 文件 | 内容 |
|---|---|
| `variants.tsv` | SNP 表，至少有 `snp`、`pos` 两列 |
| `gwas_genotypes.f32` | GWAS 个体 × SNP 的剂量矩阵，float32、小端、按行存储 |
| `qtl_genotypes_n1000.f32` | eQTL 个体 × SNP 的剂量矩阵，同上 |
| `ld_gwas.f64`、`ld_qtl_n1000.f64` | 两个样本各自的样本内 LD（相关矩阵），float64 |

`scripts/test_synthetic.R --make-bank <目录>` 会按这个格式生成两个合成窗口，可以用来试跑。

## 3.3 `R/simulate_bulk.R` 里的函数

| 函数 | 做什么 |
|---|---|
| `bulk_default_config()` | 所有模拟参数：N、q 主设置、MAF 下限、LD 档位 |
| `bulk_scenarios` | 5 个场景的表：名称、真值、每个性状的因果数 |
| `bulk_load_bank(bank_root, center)` | 读一个窗口 |
| `bulk_sample_causal(bank, scenario, tier, seed)` | 抽因果 SNP；保证档位内的 \|r\| 要求和其他对 \|r\| < 0.2 |
| `bulk_draws(bank, center, scenario, rep)` | 生成两性状的噪声，种子由（窗口、场景、重复）决定 |
| `bulk_gwas()` / `bulk_eqtl()` | 线性模型 + OLS，返回 beta、se、z、p |
| `bulk_assemble_locus()` | 组装位点并附上真值 |

代码里的场景内部编号（`A_A` 等）对应报告中的名称：`A_A` = single_shared，`A_AB` = single_double_shared，`AB_AB` = double_shared2，`A_B` = single_distinct，`AB_CD` = double_distinct。

## 3.4 两段运行

| 阶段 | 设置 | 规模 |
|---|---|---|
| smoke（冒烟测试） | 2 个窗口 × 5 场景 × 2 次重复，所有设置各跑一遍 | 几分钟到十几分钟 |
| main（第一段） | GWAS q 0.3%，eQTL q 1.5%，\|r\| < 0.2；12 窗口 × 5 场景 × 100 次 | 6,000 个位点 |
| grid（第二段） | 从主设置出发每次只变一个：GWAS q、eQTL q、因果间 LD 档位；每格 50 次 | 22,800 个位点 |

第一段里，block-coloc 各版本在 4 组 lead 阈值下各跑一次（1e-3、1e-4、1e-5、旧设置 5e-8 / 1e-5）。第二段只跑报告的三个版本，只用第一段选出的阈值，并沿用第一段的判定线，不重新调。

## 3.5 复现命令

在仓库根目录：

```bash
# 0. 自检（不需要真实数据）
Rscript scripts/test_synthetic.R

# 用合成数据试跑全流程
Rscript scripts/test_synthetic.R --make-bank /tmp/synbank
Rscript scripts/run.R --stage smoke --cores 4 --bank /tmp/synbank
Rscript scripts/evaluate.R --stage smoke

# 1. 第一段（真实数据目录用 --bank 指定）
Rscript scripts/run.R --stage main --cores 4 --bank <数据目录>
Rscript scripts/evaluate.R --stage main

# 2. 第二段（需要第一段的 thresholds.tsv）
Rscript scripts/run.R --stage grid --cores 4 --bank <数据目录>
Rscript scripts/evaluate.R --stage grid
```

- 窗口划分（6 个校准 + 6 个测试）写在 `scripts/run.R` 开头；换数据时需要相应修改。
- 中断后重跑同一条命令即可，已完成的分片会跳过。
- 所有随机数由（窗口、场景、重复、档位）决定，重跑得到完全相同的结果。
- 依赖：R ≥ 4.4，R 包 `coloc`、`susieR`（自检在 coloc 5.2.3、susieR 0.14.2 下通过）。

## 3.6 自检内容（`scripts/test_synthetic.R`）

1. 读窗口：LD 文件等于基因型的相关矩阵；
2. 抽因果：每个场景的因果数、共享数正确，LD 档位和 MAF 要求都满足；
3. OLS：`z = beta / se`，与 `lm()` 结果一致；
4. 没有信号时，z 的方差约为 1，**z 之间的相关等于样本内 LD**（即 `z ~ N(Rλ, R)`）；
5. 因果 SNP 的 |z| 约等于 `sqrt(N q)`；
6. 公共随机数：同样的种子得到同样的因果和噪声。
