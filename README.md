# bulk-block-coloc

Colocalisation when a genomic region contains multiple causal variants: generative model, theoretical properties, and simulation evaluation of **block-coloc** on UK Biobank chromosome 22 genotypes.

**block-coloc** defaults to whole-region `coloc.abf` unless stepwise conditional analysis detects multiple independent signals ($\max(\lvert\mathcal{L}_1\rvert, \lvert\mathcal{L}_2\rvert) \ge 2$) in at least one trait. When multiple signals are detected, it forms joint linkage-disequilibrium (LD) neighbourhood blocks ($\hat r^2 \ge \rho_b$) around conditional leads, evaluates `coloc.abf` on unadjusted marginal summary statistics within each joint block, and falls back to the whole region when no admissible joint block is formed. Across two continuous traits simulated on non-overlapping cohorts (trait 1: $N_1 = 20{,}000$, GWAS-like; trait 2: $N_2 = 1{,}000$, bulk-eQTL-like) with external 1000 Genomes EUR reference LD ($n_{\mathrm{ref}} = 489$), we benchmark **block-coloc** against `block-coloc (no fallback)`, `block-coloc (marginal leads)`, whole-region `coloc.abf`, `coloc.susie`, `coloc.susie (fallback)`, `coloc.signals (mask)`, and `coloc.signals (cond)`.

**Documentation site**: https://SuranCHEN.github.io/bulk-block-coloc/ (Chinese)

- Model, Algorithm, and Analytical Properties: [docs/method.md](docs/method.md)
- Simulation Design, Scenarios, and Evaluation Protocol: [docs/simulation.md](docs/simulation.md)
- Figure 1 (Workflow and Branch Decomposition): [docs/fig1.md](docs/fig1.md)
- Figure 2 (Power, Type I Error, and Signal-Strength Grids): [docs/fig2.md](docs/fig2.md)
- Figure 3 (Multi-Signal Behaviour of Whole-Region `coloc.abf`): [docs/fig3.md](docs/fig3.md)
- Figure 4 (Matched In-Sample LD vs 1000G Reference LD): [docs/fig4.md](docs/fig4.md)
- Figure 5 (Fixed Three-Way Calling Rule and Block-Prior Sensitivity): [docs/fig5.md](docs/fig5.md)
- Supplementary Proofs, Sensitivity Analyses, and Extended Ablations: [docs/supplement.md](docs/supplement.md)

## Summary of findings

**Reference-LD performance and mechanism.** Calibrated to a $5\%$ false-positive rate (FPR) on independent genomic windows, **block-coloc** detects $76.8\%$ of colocalised test loci versus $73.2\%$ for whole-region `coloc.abf` ($+3.5$ percentage points; $6/6$ test windows positive, window-level $95\%$ $t_5$ CI $[+2.2, +4.9]$; descriptive locus-level McNemar $p < 10^{-15}$) and is comparable to `coloc.signals (cond)` ($76.2\%$; paired difference $+0.6$ percentage points, $95\%$ $t_5$ CI $[-1.6, +2.8]$). As eQTL regional heritability increases ($h^2_2 \in \{1.5\%, 3.0\%, 5.0\%\}$), the advantage over `coloc.abf` widens to $+4.0$, $+7.1$, and $+13.7$ percentage points ($95.9\%$ vs $82.2\%$ at $h^2_2 = 5.0\%$), reflecting that whole-region `coloc.abf` asymptotically tests coincidence of the strongest marginal signals ($\mathbb{1}\{j_1^* = j_2^*\}$) rather than shared causality ($\mathcal{C}_1 \cap \mathcal{C}_2 \neq \varnothing$; Supplement Proposition S1). Law-of-total-probability branch decomposition shows that the gain arises from loci entering valid joint blocks ($\mathcal{G}_1$), while replacing matched in-sample LD ($\mathbf{R}_1, \mathbf{R}_2$) with 1000 Genomes EUR reference LD ($\widehat{\mathbf{R}}$) changes **block-coloc** pooled detection by at most $0.5$ percentage points (window means $-0.4$ to $+0.2$) compared with losses of $3.3$ to $6.7$ percentage points for `coloc.signals (cond)` (a deterministic subset-perturbation bound clarifies why minor block-boundary shifts have negligible impact when marginal inputs stay fixed; Supplement Proposition S5).

**Scope and limitations.** When matched in-sample LD is available, `coloc.signals (cond)` outperforms **block-coloc** ($79.5\%$ vs $77.0\%$ at default; $82.9\%$ vs $76.8\%$ at $h^2_1 = 2.0\%$); **block-coloc** is specifically motivated for summary-statistic analyses relying on external reference LD. At the evaluated sample sizes, effect strengths, and calibrated $5\%$ FPR cut-offs, discrimination of distinct causal variants degrades sharply once $\lvert r\rvert \ge 0.32$ (a finite-sample power boundary rather than asymptotic non-identifiability when $\lvert r\rvert < 1$), and under moderate cross-trait LD ($\lvert r\rvert \in [0.32, 0.63)$) at multi-signal loci **block-coloc** shows a higher FPR than `coloc.abf` ($16.3\%$ vs $12.8\%$ calibrated; $1.8\%$ vs $0.8\%$ under $\mathrm{PP4} > 0.7$). Under fixed-threshold calling ($\mathrm{PP4} > 0.7$), part of the higher call rate reflects elevated prior odds inside smaller blocks ($\Pr(H_4)/\Pr(H_3) = p_{12}/\{(M_B - 1)p_1 p_2\} \approx 24.4$ at $M_B = 42$ vs $\approx 1.1$ whole-region), whereas matched-FPR gains persist under region-equivalent block priors ($+3.0$ vs $+3.5$ percentage points, where the finite-set factor $(M_B - 1)/M_B$ yields an approximate odds match). Under $h^2_2 = 5.0\%$, calibrating solely on two-trait distinct negatives lowers the threshold to $0.012$–$0.014$, which single-trait (`gwas_only`) loci cross ($0.0\%$ under $\mathrm{PP4} > 0.7$). Finally, with fixed $\alpha_{\mathrm{lead}}$, finite $n_{\mathrm{ref}}$, or orthogonal genotypes ($\mathbf{R} = \mathbf{I}_M$) under a minimum block size of two, the data-driven algorithm is not universally consistent. Empirical evaluations use $6$ calibration and $6$ test windows on chromosome 22.

## Repository layout

| Path | Content |
|---|---|
| `docs/` | Streamlined documentation site (`method`, `simulation`, `fig1`–`fig5`, `supplement`) |
| `R/` | `block_coloc.R` (original implementation), `methods.R` (compared methods), `simulate_bulk.R` (data-generating process) |
| `scripts/` | Simulation runners, evaluation scripts, sensitivity analyses, and figure generation scripts |
| `results/` | Saved aggregate tables and markdown reports |

## Quick start (synthetic genotypes, no individual-level data required)

```bash
Rscript scripts/test_synthetic.R
Rscript scripts/test_synthetic.R --make-bank /tmp/synbank
Rscript scripts/run.R --stage smoke --cores 4 --bank /tmp/synbank --ref /tmp/synbank_ref --windows 31263287,21059630
Rscript scripts/evaluate.R --stage smoke
Rscript scripts/plot_main_figures.R
```

## Data and license

This repository contains no individual-level UK Biobank genotypes, sample identifiers, or in-sample LD matrices. Only genotype LD structure is real; causal variants, effects, and phenotypes are simulated. Licensed under GPL-3.0 (`LICENSE`).
