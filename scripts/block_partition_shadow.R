#!/usr/bin/env Rscript
# Shadow-signal scenarios for block-coloc (user decision 2026-11-13).
#
# Copy of block_partition.R with (i) its own causal sampler and (ii) the
# V4_cond_block score (block SNPs conditioned on the same trait's leads outside
# the block, methods.R .condition_outside) next to the stage-1 score.
# simulate_bulk.R and methods.R are sourced read-only.
#
# Purpose: test Proposition E (docs/11_theory.md). In the main design every pair
# of causal variants of DIFFERENT traits has |r| < 0.2, so a signal of one trait
# casts almost no shadow on a causal variant of the other. Here one cross-trait
# pair is put in a chosen LD range (--cross mid: 0.32 <= |r| < 0.63, default;
# --cross low: |r| < 0.2, the control drawn by the same sampler):
#   AB_AC  GWAS A, B; eQTL A, C; cross pair B-C; every other pair |r| < 0.2. Truth: colocalised.
#   AB_CD  GWAS A, B; eQTL C, D; cross pair B-C; every other pair |r| < 0.2. Truth: not colocalised.
# Only these two scenarios are run.
#   Rscript scripts/bulk_h2ref_v1/block_partition_shadow.R --cores 4 --reps 50 --cross mid
#   Rscript scripts/bulk_h2ref_v1/block_partition_shadow.R --cores 4 --reps 50 --cross low
# (add --h2_gwas 0.01 for a stronger GWAS shadow). Output:
# outputs/bulk_h2ref_v1_shadow_<cross>[_g<h2_gwas>_e<h2_eqtl>]/.
#
# Below: header of block_partition.R, unchanged.
#
# Bulk only. methods.R, simulate_bulk.R and R/block_coloc.R are sourced
# read-only; nothing in them is changed. The simulation part (seeds, causal
# draws, phenotypes, summary statistics, leads, blocks) is copied unchanged
# from prior_rule.R, so the loci are exactly those of stage 1 / Results §8.
#
# What is new: for every locus that passes the gate (both traits >= 1
# conditional lead, one trait >= 2), the region is split into DISJOINT parts
#   joint   the blocks of stage 1 (.build_blocks)
#   trait1  r2 >= 0.4 neighbourhood of a GWAS lead that is in no joint block
#   trait2  r2 >= 0.4 neighbourhood of an eQTL lead that is in no joint block
#   rest    all SNPs in none of the above
# A SNP in more than one neighbourhood goes to the part whose lead has the
# largest r2 with it (ties: first part). coloc.abf (per-SNP priors) is run on
# every part and PP0-PP4 are saved. Loci that do not pass the gate have one
# part, the whole region. The region-level posterior is computed afterwards
# (partition_eval.py):
#   A  = prod_b (1 - PP4_b),  C1 = prod_b (PP0_b + PP2_b),
#   C2 = prod_b (PP0_b + PP1_b),  C0 = prod_b PP0_b
#   H0 = C0, H1 = C2 - C0, H2 = C1 - C0, H3 = A - C1 - C2 + C0, H4 = 1 - A.
# stage1_pp4 (max over the original, possibly overlapping blocks) must equal
# bc_orig_pp4 of prior_rule.R locus by locus; this is the reproduction check.
#
# Usage (from the repository root, i.e. the folder containing R/block_coloc.R):
#   Rscript scripts/bulk_h2ref_v1/block_partition.R --cores 4                       # default cell, 100 reps
#   Rscript scripts/bulk_h2ref_v1/block_partition.R --cores 4 --h2_eqtl 0.05 --reps 50
# Options: --bank --ref --out --windows --reps --h2_gwas --h2_eqtl --cores
# Default output: outputs/bulk_h2ref_v1_partition (suffix _g<h2_gwas>_e<h2_eqtl>
# if either h2 differs from 0.003 / 0.015). Resumable: one shard per
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
CROSS <- .opt("cross", "mid"); stopifnot(CROSS %in% c("low", "mid"))
suffix <- paste0("_", CROSS, if (H2G != 0.003 || H2E != 0.015) sprintf("_g%s_e%s", format(H2G), format(H2E)) else "")
SCEN <- c("AB_AC", "AB_CD")
out <- .opt("out", file.path(repo, "outputs", paste0("bulk_h2ref_v1_shadow", suffix)))
cores <- as.integer(.opt("cores", "4"))
n_rep <- as.integer(.opt("reps", "50"))
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

# Causal SNPs with one cross-trait pair (B-C) in the chosen LD range.
shadow_sample_causal <- function(bank, scenario, cross, seed, cfg = CFG) {
  set.seed(seed)
  R <- abs(bank$ld_gwas)
  elig <- which(bank$gwas$maf >= cfg$maf_min & bank$qtl$maf >= cfg$maf_min)
  rng <- if (cross == "low") c(0, cfg$r_other_max) else cfg$ld_tiers$mid
  far <- function(cand, chosen) apply(R[cand, chosen, drop = FALSE] < cfg$r_other_max, 1, all)
  pick <- function(cand) if (length(cand)) cand[sample.int(length(cand), 1L)] else NA_integer_
  g <- NULL
  for (attempt in seq_len(cfg$max_draw_attempts)) {
    A <- pick(elig)
    rest <- setdiff(elig, A)
    B <- pick(rest[far(rest, A)]); if (is.na(B)) next
    rest <- setdiff(rest, B)
    C <- pick(rest[far(rest, A) & R[rest, B] >= rng[1] & R[rest, B] < rng[2]]); if (is.na(C)) next
    if (scenario == "AB_AC") { g <- c(A, B); e <- c(A, C); break }
    rest <- setdiff(rest, C)
    D <- pick(rest[far(rest, c(A, B, C))]); if (is.na(D)) next
    g <- c(A, B); e <- c(C, D); break
  }
  if (is.null(g)) stop("could not draw causal SNPs for ", scenario, " cross ", cross)
  list(scenario = scenario, tier = "low", gwas = g, eqtl = e,
       sign_gwas = sample(c(-1, 1), length(g), replace = TRUE),
       sign_eqtl = sample(c(-1, 1), length(e), replace = TRUE),
       shared = intersect(g, e), B = B, C = C, cross_r = unname(bank$ld_gwas[B, C]))
}

one_locus <- function(bank, cz, gw, eq) {
  loc <- bulk_assemble_locus(bank, cz, gw, eq, ld = "ref")
  M <- nrow(loc$sumstats); all_idx <- seq_len(M); snp <- loc$sumstats$snp; tr <- loc$truth
  flags <- function(idx) c(length(tr$shared) > 0 && any(tr$shared %in% snp[idx]),
                           any(tr$causal_gwas %in% snp[idx]), any(tr$causal_eqtl %in% snp[idx]))
  has <- function(idx, nm) if (!length(idx)) NA else length(nm) > 0 && any(nm %in% snp[idx])
  whole <- abf_pp(loc, all_idx)
  m1 <- .leads_marginal(loc, 1, GATED); m2 <- .leads_marginal(loc, 2, GATED)
  l1 <- .leads_conditional(loc, 1, GATED, m1)$kept; l2 <- .leads_conditional(loc, 2, GATED, m2)$kept
  go <- length(l1) >= 1L && length(l2) >= 1L && max(length(l1), length(l2)) >= 2L
  stage1_pp3 <- stage1_pp4 <- v4_pp3 <- v4_pp4 <- NA_real_; n_orig <- 0L; overlap_snps <- 0L
  s1_idx <- v4_idx <- integer()
  parts <- list()
  if (go) {
    blocks <- .build_blocks(loc, l1, l2, GATED); n_orig <- length(blocks)
    if (n_orig) {
      po <- t(vapply(blocks, function(b) abf_pp(loc, b$indices), numeric(5)))
      if (any(is.finite(po[, 5]))) { k <- which.max(po[, 5]); stage1_pp3 <- po[k, 4]; stage1_pp4 <- po[k, 5]
        s1_idx <- blocks[[k]]$indices }
      pv <- t(vapply(blocks, function(b) {
        target <- tryCatch(.condition_outside(loc, b, l1, l2, GATED), error = function(e) NULL)
        if (is.null(target)) rep(NA_real_, 5) else abf_pp(target, b$indices) }, numeric(5)))
      if (any(is.finite(pv[, 5]))) { k <- which.max(pv[, 5]); v4_pp3 <- pv[k, 4]; v4_pp4 <- pv[k, 5]
        v4_idx <- blocks[[k]]$indices }
    }
    cand <- lapply(blocks, function(b) list(type = "joint", leads = c(b$leads1, b$leads2), idx = b$indices))
    used1 <- unique(unlist(lapply(blocks, `[[`, "leads1"))); used2 <- unique(unlist(lapply(blocks, `[[`, "leads2")))
    for (l in setdiff(l1, used1))
      cand[[length(cand) + 1L]] <- list(type = "trait1", leads = l, idx = which(loc$ld1[, l]^2 >= GATED$block_r2))
    for (l in setdiff(l2, used2))
      cand[[length(cand) + 1L]] <- list(type = "trait2", leads = l, idx = which(loc$ld2[, l]^2 >= GATED$block_r2))
    LD <- loc$ld1
    score <- matrix(-Inf, M, length(cand))
    for (p in seq_along(cand)) {
      idx <- cand[[p]]$idx
      if (length(idx)) score[idx, p] <- apply(LD[idx, cand[[p]]$leads, drop = FALSE]^2, 1, max)
    }
    n_in <- rowSums(is.finite(score))
    overlap_snps <- sum(n_in > 1L)
    owner <- rep(NA_integer_, M)
    if (any(n_in > 0L)) owner[n_in > 0L] <- max.col(score[n_in > 0L, , drop = FALSE], ties.method = "first")
    for (p in seq_along(cand)) {
      idx <- which(owner == p)
      if (length(idx)) parts[[length(parts) + 1L]] <- list(type = cand[[p]]$type, idx = idx)
    }
    rest <- which(n_in == 0L)
    if (length(rest)) parts[[length(parts) + 1L]] <- list(type = "rest", idx = rest)
  }
  if (!length(parts)) parts <- list(list(type = "whole", idx = all_idx))
  pp <- t(vapply(parts, function(p) if (p$type == "whole") whole else abf_pp(loc, p$idx), numeric(5)))
  fl <- t(vapply(parts, function(p) flags(p$idx), logical(3)))
  part_str <- paste(vapply(seq_along(parts), function(i)
    paste(parts[[i]]$type, length(parts[[i]]$idx), paste(sprintf("%.17g", pp[i, ]), collapse = ","),
          paste(as.integer(fl[i, ]), collapse = ","), sep = "|"), ""), collapse = ";")
  data.frame(n_snps = M, n_leads1 = length(l1), n_leads2 = length(l2), gate = go, n_blocks_orig = n_orig,
             stage1_pp3 = stage1_pp3, stage1_pp4 = stage1_pp4, overlap_snps = overlap_snps,
             v4_pp3 = v4_pp3, v4_pp4 = v4_pp4,
             stage1_block_has_shared = has(s1_idx, tr$shared), stage1_block_has_B = has(s1_idx, bank$snp[cz$B]),
             stage1_block_has_C = has(s1_idx, bank$snp[cz$C]),
             v4_block_has_shared = has(v4_idx, tr$shared), v4_block_has_C = has(v4_idx, bank$snp[cz$C]),
             whole_pp0 = whole[1], whole_pp1 = whole[2], whole_pp2 = whole[3], whole_pp3 = whole[4], whole_pp4 = whole[5],
             n_parts = length(parts), parts = part_str, stringsAsFactors = FALSE)
}

run_unit <- function(bank, center, scenario, rep) {
  sc <- bulk_scenarios[bulk_scenarios$scenario == scenario, ]
  dr <- bulk_draws(bank, center, scenario, rep)
  cz <- shadow_sample_causal(bank, scenario, CROSS, .bulk_seed(center, scenario, rep, "causal", "shadow", CROSS), CFG)
  gw <- bulk_gwas(bank, cz, H2G, dr); eq <- bulk_eqtl(bank, cz, H2E, dr)
  r <- one_locus(bank, cz, gw, eq)
  cbind(data.frame(window = center, set = if (center %in% CALIBRATION) "calibration" else "test",
                   scenario = scenario, truth = sc$truth, rep = rep, h2_gwas = H2G, h2_eqtl = H2E,
                   cross = CROSS, cross_abs_r = round(abs(cz$cross_r), 3)), r)
}

shard_dir <- file.path(out, "shards_partition"); dir.create(shard_dir, showWarnings = FALSE)
t0 <- Sys.time()
for (center in windows) {
  ids <- expand.grid(rep = seq_len(n_rep), scenario = SCEN, stringsAsFactors = FALSE)
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
  if (length(err)) cat(err, file = file.path(out, "errors_partition.txt"), sep = "\n", append = TRUE)
  rm(bank); invisible(gc())
  message(sprintf("  elapsed %.1f min", as.numeric(Sys.time() - t0, units = "mins")))
}
files <- list.files(shard_dir, pattern = "\\.rds$", full.names = TRUE)
d <- do.call(rbind, lapply(files, readRDS))
write.table(d, gzfile(file.path(out, "partition_loci.tsv.gz")), sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo_partition.txt"))
message(sprintf("done: %d loci; total %.1f min; output %s", nrow(d), as.numeric(Sys.time() - t0, units = "mins"), out))
