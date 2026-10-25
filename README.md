# bulk-block-coloc

在 GWAS × bulk eQTL 共定位中，一个区域里可能有多个独立信号，这时 coloc.abf 的“每个性状至多一个因果变异”假设不成立。
本仓库评估 **block-coloc**：只在确有多个独立信号时，用 LD 把区域切块、每块跑 coloc.abf；其余情况直接跑整区域 coloc.abf。
我们在真实 UKBB chr22 基因型上做模拟，与整区域 coloc.abf、coloc.susie、coloc.signals 在相同误报率下比较。

**文档网站**：https://SuranCHEN.github.io/bulk-block-coloc/

## 内容

| 位置 | 内容 |
|---|---|
| `docs/` | 文档（即网站内容）：bulk 模型、block-coloc 流程与伪代码、代码说明、结果 |
| `R/` | 方法与模拟代码：`block_coloc.R`、`methods.R`、`simulate_bulk.R` |
| `scripts/` | 运行脚本：`run.R`、`evaluate.R`、`test_synthetic.R` |
| `results/stage1/` | 第一段（主设置）的汇总结果和图 |

报告的方法：

| 名称 | 代码里的内部名 |
|---|---|
| block-coloc（条件 lead 选择，本文方法） | `V3_gate_conditional` |
| block-coloc（边际 lead 选择） | `V2_gate_marginal` |
| block-coloc（交接包原版，对照） | `V0_handoff` |

对照：整区域 coloc.abf、coloc.susie（及无 credible set 时退回 coloc.abf 的用法）、coloc.signals（mask / cond）。

## 快速开始（不需要真实数据）

```bash
# 需要 R ≥ 4.4 和 R 包 coloc、susieR
Rscript scripts/test_synthetic.R                         # 自检 11 项
Rscript scripts/test_synthetic.R --make-bank /tmp/synbank  # 生成两个合成窗口
Rscript scripts/run.R --stage smoke --cores 4 --bank /tmp/synbank
Rscript scripts/evaluate.R --stage smoke
```

完整的运行步骤和数据格式见 [docs/03_code.md](docs/03_code.md)。

## 数据

- 本仓库**不含任何 UKBB 数据**（个体级基因型、样本编号、LD 矩阵）。UKBB 数据受数据使用协议约束，只保存在有权使用者的本地。
- 只有基因型和 LD 来自真实数据；因果变异、效应和表型均为模拟，结果不对应任何真实基因或疾病。

## 许可证

GPL-3.0，见 `LICENSE`。
