#!/usr/bin/env Rscript
# Phenotype-scale (sdY) estimation sensitivity of block-coloc and whole-region coloc.abf.
#
# Background: bulk_assemble_locus() sets no sdY, so coloc.abf estimates sdY on every SNP
# subset it is given (coloc:::sdY.est: regress 2*N*MAF*(1-MAF) on 1/varbeta through the origin). The ABF
# prior SD is 0.15*sdY, so each block gets its own prior scale. The simulated phenotypes have Var(y) = 1 in
# expectation (simulate_bulk.R header), so sdY = 1 is the true value.
#
# Same loci as stage 1 / stage 2 (simulation copied unchanged from block_ablation.R / block_partition.R).
# methods.R, simulate_bulk.R and R/block_coloc.R are sourced read-only. Stage-1 block-coloc rule only
# (clump 0.1, neighbourhood 0.4, sign split, gate, lead threshold 1e-3). Leads and blocks do not depend on
# sdY, so the three arms use the same blocks; only the PP values (and so the best block) can change:
#   est    current implementation: sdY estimated by coloc on each subset   (must reproduce stage 1 / 2)
#   whole  sdY estimated once per trait on the whole region, passed to every subset (fixed scale per locus)
#   one    sdY = 1 for both traits (true value), whole region and blocks
# Output: one row per (locus, arm). The adopted rule returns pp4_whole for status whole / no_valid_block.
# Extra columns: sdY estimates on the whole region and on the best block of arm est (own implementation
# of coloc:::sdY.est; sdy_check compares it with coloc:::sdY.est on the whole region).
#
# Usage (from the hand-off root):
#   Rscript scripts/bulk_h2ref_v1/block_sdy.R --cores 4                         # default cell, 100 reps
#   Rscript scripts/bulk_h2ref_v1/block_sdy.R --cores 4 --h2_eqtl 0.05 --reps 50
#   Rscript scripts/bulk_h2ref_v1/block_sdy.R --cores 4 --h2_gwas 0.02 --reps 50
# Output: outputs/bulk_h2ref_v1_sdy[_g<h2_gwas>_e<h2_eqtl>]/sdy_loci.tsv.gz. Resumable (one shard per unit).

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
out <- .opt("out", file.path(repo, "outputs", paste0("bulk_h2ref_v1_sdy", suffix)))
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

# coloc.abf on SNP subset idx; returns PP0-PP4 (NA on error).
abf_pp <- function(loc, idx, p1 = GATED$p1, p2 = GATED$p2, p12 = GATED$p12) {
  f <- tryCatch(quiet(coloc::coloc.abf(.coloc_dataset(loc, idx, 1), .coloc_dataset(loc, idx, 2), p1 = p1, p2 = p2, p12 = p12)),
                error = function(e) NULL)
  if (is.null(f)) return(rep(NA_real_, 5))
  unname(f$summary[c("PP.H0.abf", "PP.H1.abf", "PP.H2.abf", "PP.H3.abf", "PP.H4.abf")])
}

# Own copy of coloc:::sdY.est (coloc 5.x); NA instead of an error when the coefficient is negative.
sdy_est <- function(vbeta, maf, n) {
  oneover <- 1 / vbeta; nvx <- 2 * n * maf * (1 - maf)
  cf <- tryCatch(unname(coef(lm(nvx ~ oneover - 1))[["oneover"]]), error = function(e) NA_real_)
  if (!is.finite(cf) || cf < 0) NA_real_ else sqrt(cf)
}
sdy_of <- function(loc, idx, trait) {
  ss <- loc$sumstats
  sdy_est(ss[[paste0("se_trait", trait)]][idx]^2, ss[[paste0("maf_trait", trait)]][idx], loc[[paste0("n", trait)]])
}
ARMS <- c("est", "whole", "one")

one_locus <- function(bank, cz, gw, eq) {
  loc0 <- bulk_assemble_locus(bank, cz, gw, eq, ld = "ref")
  M <- nrow(loc0$sumstats); snp <- loc0$sumstats$snp; tr <- loc0$truth
  has <- function(idx, nm) if (!length(idx)) NA else length(nm) > 0 && any(nm %in% snp[idx])
  all_idx <- seq_len(M)
  sw <- c(sdy_of(loc0, all_idx, 1), sdy_of(loc0, all_idx, 2))
  chk <- tryCatch(max(abs(sw - c(quiet(coloc:::sdY.est(loc0$sumstats$se_trait1^2, loc0$sumstats$maf_trait1, loc0$n1)),
                                 quiet(coloc:::sdY.est(loc0$sumstats$se_trait2^2, loc0$sumstats$maf_trait2, loc0$n2))))),
                  error = function(e) NA_real_)
  # leads and blocks: sdY-free, computed once
  cfg <- GATED
  m1 <- .leads_marginal(loc0, 1, cfg); m2 <- .leads_marginal(loc0, 2, cfg)
  l1 <- .leads_conditional(loc0, 1, cfg, m1)$kept; l2 <- .leads_conditional(loc0, 2, cfg, m2)$kept
  go <- length(l1) >= 1L && length(l2) >= 1L && max(length(l1), length(l2)) >= 2L
  blocks <- if (go) .build_blocks(loc0, l1, l2, cfg) else list()
  rows <- vector("list", length(ARMS)); best_est <- integer()
  for (ai in seq_along(ARMS)) {
    arm <- ARMS[ai]; loc <- loc0
    if (arm == "whole") { loc$sdY1 <- sw[1]; loc$sdY2 <- sw[2] }
    if (arm == "one") { loc$sdY1 <- 1; loc$sdY2 <- 1 }
    whole <- abf_pp(loc, all_idx)
    status <- "whole"; pp <- if (go) c(NA_real_, NA_real_) else whole[4:5]; nb <- length(blocks); best <- integer()
    if (go) {
      status <- "no_valid_block"
      if (nb) {
        p <- t(vapply(blocks, function(b) abf_pp(loc, b$indices), numeric(5)))
        if (any(is.finite(p[, 5]))) {
          k <- which.max(p[, 5]); pp <- p[k, 4:5]; best <- blocks[[k]]$indices; status <- "block"
        }
      }
    }
    if (arm == "est") best_est <- best
    rows[[ai]] <- data.frame(arm = arm, n_leads1 = length(l1), n_leads2 = length(l2), gate = go, n_blocks = nb,
                             status = status, pp3 = pp[1], pp4 = pp[2], pp3_whole = whole[4], pp4_whole = whole[5],
                             best_size = length(best), best_has_shared = has(best, tr$shared),
                             best_same_as_est = if (length(best) && length(best_est)) identical(best, best_est) else NA,
                             stringsAsFactors = FALSE)
  }
  d <- do.call(rbind, rows)
  d$sdy1_whole <- sw[1]; d$sdy2_whole <- sw[2]; d$sdy_check <- chk
  d$sdy1_best_est <- if (length(best_est)) sdy_of(loc0, best_est, 1) else NA_real_
  d$sdy2_best_est <- if (length(best_est)) sdy_of(loc0, best_est, 2) else NA_real_
  d
}

run_unit <- function(bank, center, scenario, rep) {
  sc <- bulk_scenarios[bulk_scenarios$scenario == scenario, ]
  dr <- bulk_draws(bank, center, scenario, rep)
  cz <- bulk_sample_causal(bank, scenario, "low", .bulk_seed(center, scenario, rep, "causal", "low"), CFG)
  gw <- bulk_gwas(bank, cz, H2G, dr); eq <- bulk_eqtl(bank, cz, H2E, dr)
  r <- one_locus(bank, cz, gw, eq)
  data.frame(window = center, set = if (center %in% CALIBRATION) "calibration" else "test",
             scenario = scenario, truth = sc$truth, rep = rep, h2_gwas = H2G, h2_eqtl = H2E, r,
             stringsAsFactors = FALSE)
}

shard_dir <- file.path(out, "shards_sdy"); dir.create(shard_dir, showWarnings = FALSE)
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
  if (length(err)) cat(err, file = file.path(out, "errors_sdy.txt"), sep = "\n", append = TRUE)
  rm(bank); invisible(gc())
  message(sprintf("  elapsed %.1f min", as.numeric(Sys.time() - t0, units = "mins")))
}
files <- list.files(shard_dir, pattern = "\\.rds$", full.names = TRUE)
d <- do.call(rbind, lapply(files, readRDS))
write.table(d, gzfile(file.path(out, "sdy_loci.tsv.gz")), sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo_sdy.txt"))
message(sprintf("done: %d loci; total %.1f min; output %s", nrow(d), as.numeric(Sys.time() - t0, units = "mins"), out))
