#!/usr/bin/env Rscript
# Block priors and decision rules (user request 2026-10-29).
#
# Author: Claude (written for the user's project; executed by Codex).
# Bulk only. methods.R and R/block_coloc.R are sourced read-only; nothing in
# them is changed.
#
# Question
#   coloc.abf priors are per SNP (p1, p2, p12). A block with M_B SNPs then has
#   prior odds H4 : H3 = p12 / ((M_B - 1) p1 p2), about 22 for a 46-SNP block,
#   versus about 1 for the whole region (M about 1,000) (Supplementary Note,
#   Proposition 7). Part of block-coloc's PP4 may therefore come from the
#   smaller region rather than from the data. This run recomputes, on exactly
#   the same simulated loci as stage 1 (same seeds), for every locus:
#     whole      whole-region coloc.abf, PP0-PP4
#     bc_orig    block-coloc as in stage 1 (conditional leads, alpha 1e-3,
#                per-SNP priors inside each block), PP3 / PP4 of the best block
#     bc_adj     the same blocks, priors scaled by s = M / M_B inside a block:
#                p1 s, p2 s, p12 s. Then the prior probabilities of H0-H4 in
#                the block equal those of the whole region (M_B p1 s = M p1,
#                M_B (M_B - 1) p1 p2 s^2 ~ M (M - 1) p1 p2, M_B p12 s = M p12).
#     signals_mask, signals_cond   coloc.signals, PP3 / PP4 of the signal pair
#                with the largest PP4 (same settings as methods.R)
#   prior_rule_eval.py then compares (i) matched 5% FPR cut-offs and (ii) the
#   fixed rule used in traceCB (Jiang et al. 2026): colocalised if PP4 > 0.7,
#   independent if PP3 > 0.7, otherwise undetermined.
#
# Usage (from the repository root):
#   Rscript scripts/bulk_h2ref_v1/prior_rule.R --cores 4                     # default cell, 100 reps
#   Rscript scripts/bulk_h2ref_v1/prior_rule.R --cores 4 --h2_eqtl 0.05 --reps 50
# Options: --bank --ref --out --windows --reps --h2_gwas --h2_eqtl --cores
# Default output: outputs/bulk_h2ref_v1_prior (suffix _g<h2_gwas>_e<h2_eqtl> if
# either h2 differs from 0.003 / 0.015). Resumable: one shard per
# (window, scenario, replicate).

suppressPackageStartupMessages(library(parallel))
.args <- commandArgs(trailingOnly = TRUE)
.opt <- function(name, default) {
  i <- match(paste0("--", name), .args)
  if (is.na(i) || i == length(.args)) default else .args[i + 1]
}
script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
.find_repo <- function(start) { d <- normalizePath(start)
  repeat { if (file.exists(file.path(d, "R", "block_coloc.R"))) return(d)
    p <- dirname(d); if (p == d) stop("cannot find the repository root (folder with R/block_coloc.R)"); d <- p } }
.first_file <- function(...) { f <- c(...); f <- f[file.exists(f)]; if (!length(f)) stop("missing: ", paste(c(...), collapse = " or ")); f[1] }
repo <- .find_repo(dirname(script_file))
source(file.path(repo, "R", "block_coloc.R"))
source(.first_file(file.path(repo, "R", "methods.R"), file.path(repo, "scripts", "sc_gated_v2", "task01", "methods.R")))
source(.first_file(file.path(repo, "R", "simulate_bulk.R"), file.path(repo, "scripts", "bulk_h2ref_v1", "simulate_bulk.R")))

bank_root <- .opt("bank", file.path(repo, "data", "sc-1mb-pattern-20260909-v1"))
ref_root <- .opt("ref", file.path(repo, "data", "ref_ld_1000g_eur_20261027"))
H2G <- as.numeric(.opt("h2_gwas", "0.003")); H2E <- as.numeric(.opt("h2_eqtl", "0.015"))
suffix <- if (H2G != 0.003 || H2E != 0.015) sprintf("_g%s_e%s", format(H2G), format(H2E)) else ""
out <- .opt("out", file.path(repo, "outputs", paste0("bulk_h2ref_v1_prior", suffix)))
cores <- as.integer(.opt("cores", "4"))
n_rep <- as.integer(.opt("reps", if (suffix == "") "100" else "50"))
dir.create(out, recursive = TRUE, showWarnings = FALSE)

CALIBRATION <- c(31263287, 37199593, 39694225, 42313324, 44324730, 46518082)
TEST <- c(21059630, 24009444, 26904965, 29033189, 33343871, 34894283)
.w <- as.character(.opt("windows", "all"))
windows <- if (!identical(.w, "all")) as.numeric(strsplit(.w, ",")[[1]]) else c(CALIBRATION, TEST)
CFG <- bulk_default_config()
GATED <- gated_default_config()          # lead threshold 1e-3 (stage-1 choice), priors 1e-4 / 1e-4 / 1e-5

quiet <- function(expr) { val <- NULL; invisible(utils::capture.output(val <- suppressWarnings(suppressMessages(expr)))); val }

# coloc.abf on SNP subset idx with given priors; returns PP0-PP4.
abf_pp <- function(loc, idx, p1, p2, p12) {
  f <- tryCatch(quiet(coloc::coloc.abf(.coloc_dataset(loc, idx, 1), .coloc_dataset(loc, idx, 2), p1 = p1, p2 = p2, p12 = p12)),
                error = function(e) NULL)
  if (is.null(f)) return(rep(NA_real_, 5))
  unname(f$summary[c("PP.H0.abf", "PP.H1.abf", "PP.H2.abf", "PP.H3.abf", "PP.H4.abf")])
}

signals_pp <- function(loc, mode) {
  idx <- seq_len(nrow(loc$sumstats))
  r <- tryCatch(quiet(coloc::coloc.signals(
    .coloc_dataset(loc, idx, 1, with_ld = TRUE), .coloc_dataset(loc, idx, 2, with_ld = TRUE),
    method = mode, p1 = GATED$p1, p2 = GATED$p2, p12 = GATED$p12,
    maxhits = GATED$signals_maxhits, r2thr = GATED$signals_r2thr, pthr = GATED$signals_pthr)), error = function(e) NULL)
  if (is.null(r) || !nrow(r$summary)) return(c(NA_real_, NA_real_))
  s <- as.data.frame(r$summary); k <- which.max(s$PP.H4.abf)
  c(s$PP.H3.abf[k], s$PP.H4.abf[k])
}

one_locus <- function(bank, cz, gw, eq) {
  loc <- bulk_assemble_locus(bank, cz, gw, eq, ld = "ref")
  M <- nrow(loc$sumstats)
  whole <- abf_pp(loc, seq_len(M), GATED$p1, GATED$p2, GATED$p12)
  m1 <- .leads_marginal(loc, 1, GATED); m2 <- .leads_marginal(loc, 2, GATED)
  l1 <- .leads_conditional(loc, 1, GATED, m1)$kept; l2 <- .leads_conditional(loc, 2, GATED, m2)$kept
  go <- length(l1) >= 1L && length(l2) >= 1L && max(length(l1), length(l2)) >= 2L
  branch <- "whole"; nb <- 0L
  orig <- adj <- c(NA_real_, NA_real_); size_orig <- size_adj <- NA_integer_
  if (go) {
    blocks <- .build_blocks(loc, l1, l2, GATED); nb <- length(blocks)
    if (nb) {
      po <- t(vapply(blocks, function(b) abf_pp(loc, b$indices, GATED$p1, GATED$p2, GATED$p12)[4:5], numeric(2)))
      pa <- t(vapply(blocks, function(b) { s <- M / length(b$indices)
        abf_pp(loc, b$indices, GATED$p1 * s, GATED$p2 * s, GATED$p12 * s)[4:5] }, numeric(2)))
      ko <- if (any(is.finite(po[, 2]))) which.max(po[, 2]) else NA
      ka <- if (any(is.finite(pa[, 2]))) which.max(pa[, 2]) else NA
      if (!is.na(ko)) { orig <- po[ko, ]; size_orig <- length(blocks[[ko]]$indices) }
      if (!is.na(ka)) { adj <- pa[ka, ]; size_adj <- length(blocks[[ka]]$indices) }
      branch <- if (is.na(ko)) "no_valid_block" else "block"
    } else branch <- "no_valid_block"
  }
  sm <- signals_pp(loc, "mask"); sc <- signals_pp(loc, "cond")
  data.frame(n_snps = M, n_leads1 = length(l1), n_leads2 = length(l2), branch = branch, n_blocks = nb,
             whole_pp0 = whole[1], whole_pp1 = whole[2], whole_pp2 = whole[3], whole_pp3 = whole[4], whole_pp4 = whole[5],
             bc_orig_pp3 = orig[1], bc_orig_pp4 = orig[2], bc_orig_size = size_orig,
             bc_adj_pp3 = adj[1], bc_adj_pp4 = adj[2], bc_adj_size = size_adj,
             sig_mask_pp3 = sm[1], sig_mask_pp4 = sm[2], sig_cond_pp3 = sc[1], sig_cond_pp4 = sc[2])
}

run_unit <- function(bank, center, scenario, rep) {
  sc <- bulk_scenarios[bulk_scenarios$scenario == scenario, ]
  dr <- bulk_draws(bank, center, scenario, rep)
  cz <- bulk_sample_causal(bank, scenario, "low", .bulk_seed(center, scenario, rep, "causal", "low"), CFG)
  gw <- bulk_gwas(bank, cz, H2G, dr); eq <- bulk_eqtl(bank, cz, H2E, dr)
  r <- one_locus(bank, cz, gw, eq)
  cbind(data.frame(window = center, set = if (center %in% CALIBRATION) "calibration" else "test",
                   scenario = scenario, truth = sc$truth, rep = rep, h2_gwas = H2G, h2_eqtl = H2E), r)
}

shard_dir <- file.path(out, "shards_prior"); dir.create(shard_dir, showWarnings = FALSE)
t0 <- Sys.time()
for (center in windows) {
  ids <- expand.grid(rep = seq_len(n_rep), scenario = bulk_scenarios$scenario, stringsAsFactors = FALSE)
  ids$file <- file.path(shard_dir, sprintf("%s_%s_%03d.rds", format(center, scientific = FALSE), ids$scenario, ids$rep))
  todo <- ids[!file.exists(ids$file), ]
  message(sprintf("[%s] window %s: %d of %d units to run", format(Sys.time(), "%H:%M:%S"),
                  format(center, scientific = FALSE), nrow(todo), nrow(ids)))
  if (!nrow(todo)) next
  bank <- bulk_load_bank(bank_root, center, CFG, ref_root = ref_root)
  res <- mclapply(seq_len(nrow(todo)), function(i) {
    x <- tryCatch(run_unit(bank, center, todo$scenario[i], todo$rep[i]), error = function(e) e)
    if (inherits(x, "error")) return(sprintf("%s %s %d: %s", center, todo$scenario[i], todo$rep[i], conditionMessage(x)))
    tmp <- paste0(todo$file[i], ".tmp"); saveRDS(x, tmp); file.rename(tmp, todo$file[i]); NA_character_
  }, mc.cores = cores, mc.preschedule = FALSE)
  err <- unlist(res); err <- err[!is.na(err)]
  if (length(err)) cat(err, file = file.path(out, "errors_prior.txt"), sep = "\n", append = TRUE)
  rm(bank); invisible(gc())
  message(sprintf("  elapsed %.1f min", as.numeric(Sys.time() - t0, units = "mins")))
}
files <- list.files(shard_dir, pattern = "\\.rds$", full.names = TRUE)
d <- do.call(rbind, lapply(files, readRDS))
write.table(d, gzfile(file.path(out, "prior_rule_loci.tsv.gz")), sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo_prior.txt"))
message(sprintf("done: %d loci; total %.1f min; output %s", nrow(d), as.numeric(Sys.time() - t0, units = "mins"), out))
