#!/usr/bin/env Rscript
# Ablation of block-coloc's construction rules (user decision 2026-11-16).
#
# Same loci as stage 1 / stage 2 (simulation copied unchanged from block_partition.R
# and prior_rule.R: same seeds, causal draws, phenotypes and summary statistics).
# methods.R, simulate_bulk.R and R/block_coloc.R are sourced read-only.
#
# For every locus, block-coloc (conditional leads, lead threshold 1e-3, coloc priors
# 1e-4 / 1e-4 / 1e-5) is recomputed under each configuration below; one output row
# per (locus, configuration):
#   base     clump r2 0.1, neighbourhood r2 0.4, sign split on, gate on  (= stage-1 method)
#   nosplit  sign split off (a connected component with mixed lead signs is one block)
#   rb0.2 / rb0.3 / rb0.5 / rb0.6   neighbourhood threshold rho_b (r2) changed
#   rc0.05 / rc0.2                  clump threshold rho_c (r2) changed (changes the leads)
#   nogate   blocks whenever both traits have >= 1 conditional lead (no ">= 2 leads" gate)
# status: whole (gate not passed) / block / no_valid_block. The adopted rule returns
# pp4_whole for whole and no_valid_block. Reproduction check: for config base, pp4 and
# pp4_whole must equal V3_gate_conditional (1e-3) pp4 and pp4_whole in the stage-1 /
# stage-2 score files (ablation_eval.py does this).
#
# Usage (from the hand-off root):
#   Rscript scripts/bulk_h2ref_v1/block_ablation.R --cores 4                         # default cell, 100 reps
#   Rscript scripts/bulk_h2ref_v1/block_ablation.R --cores 4 --h2_eqtl 0.05 --reps 50
#   Rscript scripts/bulk_h2ref_v1/block_ablation.R --cores 4 --h2_gwas 0.02 --reps 50
# Output: outputs/bulk_h2ref_v1_ablation[_g<h2_gwas>_e<h2_eqtl>]/ablation_loci.tsv.gz.
# Resumable: one shard per (window, scenario, replicate).

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
out <- .opt("out", file.path(repo, "outputs", paste0("bulk_h2ref_v1_ablation", suffix)))
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

CONFIGS <- list(
  base = list(),
  nosplit = list(sign_split = FALSE),
  rb0.2 = list(block_r2 = 0.2), rb0.3 = list(block_r2 = 0.3), rb0.5 = list(block_r2 = 0.5), rb0.6 = list(block_r2 = 0.6),
  rc0.05 = list(clump_r2 = 0.05), rc0.2 = list(clump_r2 = 0.2),
  nogate = list(nogate = TRUE))

one_locus <- function(bank, cz, gw, eq) {
  loc <- bulk_assemble_locus(bank, cz, gw, eq, ld = "ref")
  M <- nrow(loc$sumstats); snp <- loc$sumstats$snp; tr <- loc$truth
  has <- function(idx, nm) if (!length(idx)) NA else length(nm) > 0 && any(nm %in% snp[idx])
  whole <- abf_pp(loc, seq_len(M))
  cache <- new.env(hash = TRUE)
  blockpp <- function(idx) {
    key <- paste(idx, collapse = ",")
    v <- cache[[key]]
    if (is.null(v)) { v <- abf_pp(loc, idx); assign(key, v, envir = cache) }
    v
  }
  leadcache <- list()
  rows <- vector("list", length(CONFIGS))
  for (ci in seq_along(CONFIGS)) {
    cn <- names(CONFIGS)[ci]
    cfg <- utils::modifyList(GATED, CONFIGS[[cn]])
    lk <- format(cfg$clump_r2)
    if (is.null(leadcache[[lk]])) {
      m1 <- .leads_marginal(loc, 1, cfg); m2 <- .leads_marginal(loc, 2, cfg)
      leadcache[[lk]] <- list(l1 = .leads_conditional(loc, 1, cfg, m1)$kept, l2 = .leads_conditional(loc, 2, cfg, m2)$kept)
    }
    l1 <- leadcache[[lk]]$l1; l2 <- leadcache[[lk]]$l2
    go <- length(l1) >= 1L && length(l2) >= 1L && (isTRUE(cfg$nogate) || max(length(l1), length(l2)) >= 2L)
    status <- "whole"; pp <- if (go) c(NA_real_, NA_real_) else whole[4:5]; nb <- 0L; best <- integer()
    if (go) {
      blocks <- .build_blocks(loc, l1, l2, cfg); nb <- length(blocks)
      status <- "no_valid_block"
      if (nb) {
        p <- t(vapply(blocks, function(b) blockpp(b$indices), numeric(5)))
        if (any(is.finite(p[, 5]))) {
          k <- which.max(p[, 5]); pp <- p[k, 4:5]; best <- blocks[[k]]$indices; status <- "block"
        }
      }
    }
    rows[[ci]] <- data.frame(config = cn, n_leads1 = length(l1), n_leads2 = length(l2), gate = go, n_blocks = nb,
                             status = status, pp3 = pp[1], pp4 = pp[2], pp3_whole = whole[4], pp4_whole = whole[5],
                             best_size = length(best), best_has_shared = has(best, tr$shared),
                             best_has_c1 = has(best, tr$causal_gwas), best_has_c2 = has(best, tr$causal_eqtl),
                             stringsAsFactors = FALSE)
  }
  do.call(rbind, rows)
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

shard_dir <- file.path(out, "shards_ablation"); dir.create(shard_dir, showWarnings = FALSE)
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
  if (length(err)) cat(err, file = file.path(out, "errors_ablation.txt"), sep = "\n", append = TRUE)
  rm(bank); invisible(gc())
  message(sprintf("  elapsed %.1f min", as.numeric(Sys.time() - t0, units = "mins")))
}
files <- list.files(shard_dir, pattern = "\\.rds$", full.names = TRUE)
d <- do.call(rbind, lapply(files, readRDS))
write.table(d, gzfile(file.path(out, "ablation_loci.tsv.gz")), sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo_ablation.txt"))
message(sprintf("done: %d loci; total %.1f min; output %s", nrow(d), as.numeric(Sys.time() - t0, units = "mins"), out))
