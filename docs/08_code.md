# Code

运行脚本在 `scripts/`，模拟函数在 `R/simulate_bulk.R`，方法定义在 `R/methods.R`，交接包原版在 `R/block_coloc.R`。上一版代码见本仓库的 git 历史（2026-10-26 及以前的提交）。

## Files

| 文件 | 作用 |
|---|---|
| `prepare_ref_ld.R` | 读 1000 Genomes EUR chr22 PLINK 文件，把每个窗口的 SNP 按 rsID、位置、等位基因对齐，写出参考 LD（`ld_ref.f64`）和匹配表（`ref_snps.tsv`），并与样本内 LD 比较 |
| `simulate_bulk.R` | 只定义函数：读窗口（限制到参考人群中有的 SNP）、按场景和 LD 档位抽因果变异、按 $$h^2/K$$ 生成表型、逐 SNP OLS、组装位点 |
| `run.R` | 在 12 个窗口上模拟并对每个位点运行所有方法；`--stage smoke / main / grid`，`--ld ref / insample`；按（窗口、场景、重复）分片保存，可断点续跑 |
| `evaluate.R` | 在校准窗口上定判定线和 lead 阈值，在测试窗口上计算功效、误报率、AUC、配对检验，输出表格、图和 `RESULTS*.md` |
| `param_checks.py` | [Supplementary Note](09_parameters.md) Proposition 4–5 的 Monte Carlo（无信号、单因果），并汇总第一段的观察值（Python 3，numpy、pandas；需要本地数据） |
| `grid_paired.py` | 第二段的检出率、误报率（含 fallback 规则）、配对 McNemar 比较和 lead 计数（Python 3，numpy、pandas、scipy；需要本地的评分表） |
| `prior_rule.R` | 块内先验与判定规则（Results §8）：在与第一段相同的位点上计算整区域 coloc.abf、block-coloc（每 SNP / 区域等价先验）、coloc.signals 的 PP3 与 PP4；不改 methods.R；需要本地数据 |
| `prior_rule_eval.py` | 上述结果的评价：误报率匹配与 traceCB 规则（PP4 > 0.7 / PP3 > 0.7），配对比较，并与第一段逐位点核对（Python 3，numpy、pandas、scipy） |
| `param_checks_two.py` | Proposition 5 的双因果部分：按 Algorithm 1 Step 1–2 计数 lead；参数为 $$\lambda$$ |
| `test_synthetic.R` | 合成数据自检 22 项：LD 文件、参考 LD 对齐与符号、因果抽样、OLS、$$\mathbf{z}$$ 的零分布、$$\operatorname{Var}(g) = h^2$$、$$\mathbb{E}\lvert z\rvert \approx \sqrt{N h^2/K}$$、公共随机数 |

## Data flow

```text
data/sc-1mb-pattern-20260909-v1/<window>/          UKBB 基因型（20,000 + 1,000 人）与样本内 LD
data/genotype_reference/1000G_EUR_chr22/            1000 Genomes EUR（489 人）
        │ prepare_ref_ld.R
        ▼
data/ref_ld_1000g_eur_20261027/<window>/            参考 LD + 匹配 SNP
        │ run.R
        ▼
outputs/bulk_h2ref_v1/shards_<stage>/*.rds          每个（窗口、场景、重复）一个分片
outputs/bulk_h2ref_v1/scores_<stage>.tsv.gz         每个位点 × 方法的 PP4、分支、lead 数
outputs/bulk_h2ref_v1/loci_<stage>.tsv.gz           每个位点的输入特征（因果 |z|、因果间 r 等）
        │ evaluate.R
        ▼
outputs/bulk_h2ref_v1/RESULTS*.md, *.tsv, *.png
```

## Reproduction

```bash
Rscript scripts/test_synthetic.R
Rscript scripts/prepare_ref_ld.R
Rscript scripts/run.R --stage smoke --cores 4 && Rscript scripts/evaluate.R --stage smoke
Rscript scripts/run.R --stage main  --cores 4 && Rscript scripts/evaluate.R --stage main
Rscript scripts/run.R --stage grid  --cores 3 && Rscript scripts/evaluate.R --stage grid
python3 scripts/grid_paired.py
Rscript scripts/prior_rule.R --cores 4 && python3 scripts/prior_rule_eval.py --out outputs/bulk_h2ref_v1_prior --stage1 outputs/bulk_h2ref_v1
python3 scripts/param_checks.py
python3 scripts/param_checks_two.py 1.94 2.75 3.16 3.9 5.06 5.48 10.03 14.21
```

UK Biobank 个体级数据只保存在本地，不在任何仓库中；公开仓库只含代码、文档和汇总结果。
