# bulk-block-coloc

Colocalisation when a region contains more than one causal variant: a simulation study of **block-coloc** on UK Biobank chr22 genotypes.

block-coloc runs whole-region coloc.abf unless conditional analysis finds multiple independent signals; in that case it splits the region into LD blocks and runs coloc.abf within blocks (and falls back to the whole region if no block contains leads of both traits). We simulate two continuous traits on real genotypes (trait 1: N = 20,000, GWAS-like; trait 2: N = 1,000, bulk-eQTL-like; non-overlapping samples), give every method the 1000 Genomes EUR reference LD, and compare with whole-region coloc.abf, coloc.susie and coloc.signals at a matched 5% false-positive rate.

**Documentation**: https://SuranCHEN.github.io/bulk-block-coloc/ （中文）

- Model: [docs/01_bulk_model.md](docs/01_bulk_model.md)
- Method: [docs/02_block_coloc_pseudocode.md](docs/02_block_coloc_pseudocode.md)
- Simulation: [docs/03_simulation.md](docs/03_simulation.md)
- Results: [docs/04_results.md](docs/04_results.md)
- Supplementary Note (derivations): [docs/09_parameters.md](docs/09_parameters.md)

## Main result

**Stage 1 (default setting).** At a matched 5% false-positive rate, block-coloc detects 76.8% of colocalised loci, 3.5 points more than whole-region coloc.abf (73.2%, McNemar p < 1e-15) and similar to coloc.signals (cond) (76.2%). The gain is limited because, with the regional eQTL heritability split over two causal variants, each eQTL signal has expected |z| of about 2.75 and a second eQTL signal is found in only about 15% of loci.

**Stage 2 (one factor at a time).** The gain over coloc.abf grows with eQTL strength: +4.0, +7.1 and +13.7 points at eQTL h2 = 1.5%, 3% and 5% (95.9% vs 82.2%). The fraction of two-causal eQTL loci with two conditional leads (4%, 16%, 56%, 89%) matched the values predicted before the run. coloc.signals (cond) detects more when either trait is strong, but its false-positive rate rises to 6-9% with the stage-1 cut-offs. When the causal variants of a trait are in LD (|r| >= 0.32), no method separates distinct causal variants.

block-coloc returns the whole-region coloc.abf PP4 when no block contains leads of both traits; the version without this rule is reported as a sensitivity analysis.

## Layout

| Path | Content |
|---|---|
| `docs/` | Documentation site |
| `R/` | `block_coloc.R` (original implementation), `methods.R` (all compared methods), `simulate_bulk.R` (data-generating process) |
| `scripts/` | `prepare_ref_ld.R`, `run.R`, `evaluate.R`, `test_synthetic.R` |
| `results/` | Aggregate tables and figures |

## Quick start (no individual-level data needed)

```bash
# R >= 4.4 with coloc and susieR
Rscript scripts/test_synthetic.R                                   # 22 checks on synthetic genotypes
Rscript scripts/test_synthetic.R --make-bank /tmp/synbank          # synthetic windows + reference LD
Rscript scripts/run.R --stage smoke --cores 4 --bank /tmp/synbank --ref /tmp/synbank_ref --windows 31263287,21059630
Rscript scripts/evaluate.R --stage smoke
```

## Data

This repository contains no UK Biobank data (individual-level genotypes, sample identifiers or in-sample LD matrices). Only genotypes and LD are real; causal variants, effects and phenotypes are simulated and do not correspond to any real gene or disease. The 1000 Genomes reference panel is public.

## License

GPL-3.0, see `LICENSE`.
