#!/usr/bin/env Rscript
# Bulk block-coloc (h2 split over causal SNPs; 1000 Genomes reference LD):
# simulate loci on the 12 windows and run every method on the same inputs.
# Evaluation is done separately by evaluate.R.
#
# Author: Claude (written for the user's project; executed by Codex).
# Bulk only: no single-cell code is sourced or run.
#
# Usage (from the repository root):
#   Once: Rscript scripts/bulk_h2ref_v1/prepare_ref_ld.R   (1000G EUR reference LD)
#   Smoke test (1 calibration + 1 test window, 2 replicates, every cell):
#     Rscript scripts/bulk_h2ref_v1/run.R --stage smoke --cores 4
#   Stage 1 (main setting, 100 replicates):
#     Rscript scripts/bulk_h2ref_v1/run.R --stage main --cores 4
#   Stage 2 (one factor at a time, 50 replicates; needs thresholds.tsv from
#   evaluate.R --stage main, and runs only after the user has seen stage 1):
#     Rscript scripts/bulk_h2ref_v1/run.R --stage grid --cores 4
# Options: --bank <dir> --ref <dir> --ld ref|insample --out <dir>
#          --windows <comma list> --reps <n>
#   --ld ref (default): every method gets the 1000G EUR reference LD.
#   --ld insample: each trait gets its own in-sample LD (comparison run; the
#     default output folder is then outputs/bulk_h2ref_v1_insample).
#   Both choices use the same SNPs (those present in the reference).
# Resumable: one shard per (window, scenario, replicate) in <out>/shards_<stage>.
#
# Design
#   Windows: 6 calibration + 6 test (same split as the earlier 12-window work).
#   Scenarios: A_A, A_AB, AB_AB, AB_AC (shared) and A_B, AB_CD (not shared),
#     plus the extra negatives gwas_only and null_both; see
#     simulate_bulk.R.
#   Strength: local h2 of each trait, split equally over its K causal SNPs.
#   Stage 1 cell: h2_gwas 0.3%, h2_eqtl 1.5%, LD tier low.
#   Stage 2 cells (one factor at a time from stage 1):
#     h2_gwas 0.1%, 1%, 2%; h2_eqtl 0.75%, 3%, 5%; LD tier mid, high.
#   Methods: all 10 in methods.R at lead threshold 1e-3; V0-V4 again at 1e-4,
#     1e-5 and at the earlier bulk setting (GWAS 5e-8, eQTL 1e-5). In stage 2
#     V0-V4 run only at the lead threshold selected in stage 1.
#   LD given to the methods: 1000 Genomes EUR (489 people) for both traits
#     (--ld ref), or each trait's in-sample LD (--ld insample). The causal SNPs
#     and LD tiers are always defined by the true (in-sample GWAS) LD.
#   coloc.abf priors p1 = p2 = 1e-4, p12 = 1e-5 (methods.R defaults).

suppressPackageStartupMessages(library(parallel))

.args <- commandArgs(trailingOnly = TRUE)
.opt <- function(name, default) {
  i <- match(paste0("--", name), .args)
  if (is.na(i)) return(default)
  if (i == length(.args) || startsWith(.args[i + 1], "--")) return(TRUE)
  .args[i + 1]
}
script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
# Repository root = the nearest parent folder that contains R/block_coloc.R.
# Works in the project layout (scripts/bulk_h2ref_v1/) and in the public
# repository layout (scripts/ with R/methods.R and R/simulate_bulk.R).
.find_repo <- function(start) { d <- normalizePath(start)
  repeat { if (file.exists(file.path(d, "R", "block_coloc.R"))) return(d)
    p <- dirname(d); if (p == d) stop("cannot find the repository root (folder with R/block_coloc.R)"); d <- p } }
.first_file <- function(...) { f <- c(...); f <- f[file.exists(f)]; if (!length(f)) stop("missing: ", paste(c(...), collapse = " or ")); f[1] }
repo <- .find_repo(dirname(script_file))
stage <- as.character(.opt("stage", "main"))
if (!stage %in% c("smoke", "main", "grid")) stop("--stage must be smoke, main or grid")
bank_root <- .opt("bank", file.path(repo, "data", "sc-1mb-pattern-20260909-v1"))
ref_root <- .opt("ref", file.path(repo, "data", "ref_ld_1000g_eur_20261027"))
LDMODE <- as.character(.opt("ld", "ref"))
if (!LDMODE %in% c("ref", "insample")) stop("--ld must be ref or insample")
if (!dir.exists(ref_root)) stop("reference LD folder not found: ", ref_root, " (run prepare_ref_ld.R first)")
out <- .opt("out", file.path(repo, "outputs", if (LDMODE == "ref") "bulk_h2ref_v1" else "bulk_h2ref_v1_insample"))
dir.create(out, recursive = TRUE, showWarnings = FALSE)
cores <- as.integer(.opt("cores", "4"))
n_rep <- as.integer(.opt("reps", switch(stage, smoke = "2", main = "100", grid = "50")))

source(file.path(repo, "R", "block_coloc.R"))
source(.first_file(file.path(repo, "R", "methods.R"),                                   # public repository
                   file.path(repo, "scripts", "sc_gated_v2", "task01", "methods.R")))  # project (generic method definitions)
source(.first_file(file.path(repo, "R", "simulate_bulk.R"),
                   file.path(repo, "scripts", "bulk_h2ref_v1", "simulate_bulk.R")))

CALIBRATION <- c(31263287, 37199593, 39694225, 42313324, 44324730, 46518082)
TEST <- c(21059630, 24009444, 26904965, 29033189, 33343871, 34894283)
.w <- as.character(.opt("windows", "all"))
windows <- if (!identical(.w, "all")) as.numeric(strsplit(.w, ",")[[1]]) else
  if (stage == "smoke") c(CALIBRATION[1], TEST[1]) else c(CALIBRATION, TEST)

CFG <- bulk_default_config()
GATED <- gated_default_config()
MAIN <- data.frame(h2_gwas = CFG$h2_gwas, h2_eqtl = CFG$h2_eqtl, tier = "low", stringsAsFactors = FALSE)
GRID <- rbind(data.frame(h2_gwas = c(0.001, 0.01, 0.02), h2_eqtl = CFG$h2_eqtl, tier = "low"),
              data.frame(h2_gwas = CFG$h2_gwas, h2_eqtl = c(0.0075, 0.03, 0.05), tier = "low"),
              data.frame(h2_gwas = CFG$h2_gwas, h2_eqtl = CFG$h2_eqtl, tier = c("mid", "high")),
              stringsAsFactors = FALSE)
CELLS <- switch(stage, main = MAIN, grid = GRID, smoke = rbind(MAIN, GRID))
BLOCK_METHODS <- c("V0_handoff", "V1_block_any", "V2_gate_marginal", "V3_gate_conditional", "V4_cond_block")
OTHER_METHODS <- setdiff(gated_methods, BLOCK_METHODS)

# Lead thresholds (GWAS, eQTL) for V0-V4. Stage 2 uses the stage-1 choice.
if (stage == "grid") {
  thf <- file.path(out, "thresholds.tsv")
  if (!file.exists(thf)) stop("stage 2 needs ", thf, " (run evaluate.R --stage main first)")
  th <- utils::read.delim(thf, stringsAsFactors = FALSE)
  # only the block-coloc versions that are reported (user decision 2026-10-25)
  th <- th[th$selected & th$method %in% c("V3_gate_conditional", "V2_gate_marginal", "V0_handoff"), ]
  LEAD <- unique(data.frame(t1 = th$lead_threshold1, t2 = th$lead_threshold2))
  LEAD_OF <- setNames(paste(th$lead_threshold1, th$lead_threshold2), th$method)
} else {
  LEAD <- data.frame(t1 = c(1e-3, 1e-4, 1e-5, 5e-8), t2 = c(1e-3, 1e-4, 1e-5, 1e-5))
  LEAD_OF <- NULL
}

shard_dir <- file.path(out, paste0("shards_", stage))
dir.create(shard_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------- helpers ----

quiet <- function(expr) { val <- NULL; invisible(utils::capture.output(val <- suppressWarnings(suppressMessages(expr)))); val }

# SuSiE localisation: union of the two credible sets of the signal pair with the largest PP4.
susie_unit <- function(f1, f2) {
  if (is.null(f1) || is.null(f2) || inherits(f1, "error") || inherits(f2, "error")) return(character())
  if (!length(f1$sets$cs) || !length(f2$sets$cs)) return(character())
  cs <- tryCatch(quiet(coloc::coloc.susie(f1, f2, p1 = GATED$p1, p2 = GATED$p2, p12 = GATED$p12)),
                 error = function(e) NULL)
  if (is.null(cs) || !nrow(cs$summary)) return(character())
  s <- as.data.frame(cs$summary); k <- which.max(s$PP.H4.abf)
  nm1 <- names(f1$pip); nm2 <- names(f2$pip)
  unique(c(nm1[f1$sets$cs[[paste0("L", s$idx1[k])]]], nm2[f2$sets$cs[[paste0("L", s$idx2[k])]]]))
}

one_locus <- function(bank, cz, gw, eq, cell) {
  loc <- bulk_assemble_locus(bank, cz, gw, eq, ld = LDMODE)
  names_ok <- function(f) { if (!inherits(f, "error") && is.null(names(f$pip))) names(f$pip) <- loc$sumstats$snp; f }
  f1 <- tryCatch(names_ok(.susie_fit(loc, 1, GATED)), error = function(e) e)
  f2 <- tryCatch(names_ok(.susie_fit(loc, 2, GATED)), error = function(e) e)
  if (!inherits(f1, "error")) loc$susie1 <- f1
  if (!inherits(f2, "error")) loc$susie2 <- f2
  res <- quiet(run_methods(loc, OTHER_METHODS, GATED))
  res$lead_threshold1 <- GATED$p_threshold1; res$lead_threshold2 <- GATED$p_threshold2
  blk <- lapply(seq_len(nrow(LEAD)), function(i) {
    m <- if (is.null(LEAD_OF)) BLOCK_METHODS else names(LEAD_OF)[LEAD_OF == paste(LEAD$t1[i], LEAD$t2[i])]
    r <- quiet(run_methods(loc, m, utils::modifyList(GATED, list(p_threshold1 = LEAD$t1[i], p_threshold2 = LEAD$t2[i]))))
    r$lead_threshold1 <- LEAD$t1[i]; r$lead_threshold2 <- LEAD$t2[i]; r })
  res <- do.call(rbind, c(list(res), blk))
  su <- susie_unit(loc$susie1, loc$susie2)
  res$unit_snps <- ifelse(res$branch == "block", res$best_block_snps,
                          ifelse(res$branch == "susie", paste(su, collapse = ";"), ""))
  z1 <- loc$sumstats$beta_trait1 / loc$sumstats$se_trait1
  z2 <- loc$sumstats$beta_trait2 / loc$sumstats$se_trait2
  zg <- abs(gw$z[cz$gwas]); ze <- abs(eq$z[cz$eqtl])
  lc <- data.frame(n_snps = nrow(loc$sumstats),
    max_abs_z_gwas = max(abs(z1)), max_abs_z_eqtl = max(abs(z2)),
    causal_abs_z_gwas = paste(round(zg, 2), collapse = ";"),
    causal_abs_z_eqtl = paste(round(ze, 2), collapse = ";"),
    pair_abs_r = paste(round(abs(cz$pair_r), 3), collapse = ";"),
    cross_abs_r = paste(round(abs(cz$cross_r), 3), collapse = ";"),
    eqtl_shared_stronger = if (cz$scenario %in% c("A_AB", "AB_AC")) ze[1] > ze[2] else NA,
    gwas_shared_stronger = if (cz$scenario == "AB_AC") zg[1] > zg[2] else NA,
    var_y_gwas = gw$var_y, var_y_eqtl = eq$var_y, var_g_gwas = gw$var_g, var_g_eqtl = eq$var_g,
    n_cs_gwas = if (inherits(f1, "error")) NA else length(f1$sets$cs),
    n_cs_eqtl = if (inherits(f2, "error")) NA else length(f2$sets$cs))
  res <- cbind(cell, res, row.names = NULL)
  attr(res, "locus") <- cbind(cell, lc, row.names = NULL)
  res
}

run_unit <- function(bank, center, scenario, rep) {
  sc <- bulk_scenarios[bulk_scenarios$scenario == scenario, ]
  cells <- CELLS[sc$has_tier_pair | CELLS$tier == "low", , drop = FALSE]
  dr <- bulk_draws(bank, center, scenario, rep)
  czs <- list(); gws <- list(); eqs <- list(); rows <- list(); loci <- list(); skipped <- character()
  for (i in seq_len(nrow(cells))) {
    cell <- cells[i, , drop = FALSE]; tr <- cell$tier
    # a tier whose causal SNPs cannot be drawn in this window is skipped and
    # logged; the other tiers of the unit still run
    if (is.null(czs[[tr]])) czs[[tr]] <- tryCatch(bulk_sample_causal(bank, scenario, tr, .bulk_seed(center, scenario, rep, "causal", tr), CFG),
                                                  error = function(e) e)
    cz <- czs[[tr]]
    if (inherits(cz, "error")) { skipped <- unique(c(skipped, sprintf("%s %s %d tier %s skipped: %s", center, scenario, rep, tr, conditionMessage(cz)))); next }
    kg <- paste(tr, cell$h2_gwas); ke <- paste(tr, cell$h2_eqtl)
    if (is.null(gws[[kg]])) gws[[kg]] <- bulk_gwas(bank, cz, cell$h2_gwas, dr)
    if (is.null(eqs[[ke]])) eqs[[ke]] <- bulk_eqtl(bank, cz, cell$h2_eqtl, dr)
    r <- one_locus(bank, cz, gws[[kg]], eqs[[ke]], cell)
    lc <- attr(r, "locus"); attr(r, "locus") <- NULL
    shared <- bank$snp[cz$shared]; gc_snp <- bank$snp[cz$gwas]; ec_snp <- bank$snp[cz$eqtl]
    units <- strsplit(r$unit_snps, ";")
    r$unit_size <- vapply(units, function(u) sum(u != ""), 0L)
    r$unit_has_shared <- vapply(units, function(u) if (!length(shared) || !length(u)) NA else any(shared %in% u), NA)
    r$unit_has_gwas_causal <- vapply(units, function(u) if (!length(u)) NA else any(gc_snp %in% u), NA)
    r$unit_has_eqtl_causal <- vapply(units, function(u) if (!length(u)) NA else any(ec_snp %in% u), NA)
    r$best_block_snps <- NULL; r$unit_snps <- NULL
    rows[[i]] <- r
    loci[[i]] <- cbind(lc, causal_gwas = paste(gc_snp, collapse = ";"), causal_eqtl = paste(ec_snp, collapse = ";"))
  }
  set <- if (center %in% CALIBRATION) "calibration" else "test"
  id <- data.frame(window = center, set = set, scenario = scenario, truth = sc$truth, rep = rep,
                   n_causal_gwas = sc$n_causal_gwas, n_causal_eqtl = sc$n_causal_eqtl, ld_panel = LDMODE)
  if (!length(rows)) stop(paste(skipped, collapse = "; "))
  list(scores = cbind(id, do.call(rbind, rows), row.names = NULL),
       loci = cbind(id, do.call(rbind, loci), row.names = NULL), skipped = skipped)
}

# ------------------------------------------------------------------ run ----

t_start <- Sys.time()
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
    tmp <- paste0(todo$file[i], ".tmp"); saveRDS(x, tmp); file.rename(tmp, todo$file[i])
    if (length(x$skipped)) paste(x$skipped, collapse = "\n") else NA_character_
  }, mc.cores = cores, mc.preschedule = FALSE)
  err <- unlist(res); err <- err[!is.na(err)]
  if (length(err)) {
    cat(err, file = file.path(out, paste0("errors_", stage, ".txt")), sep = "\n", append = TRUE)
    message(length(err), " unit errors or skipped LD tiers (see errors_", stage, ".txt)")
  }
  rm(bank); invisible(gc())
  message(sprintf("  elapsed %.1f min", as.numeric(Sys.time() - t_start, units = "mins")))
}

# collect (only the windows of this run's window set)
files <- list.files(shard_dir, pattern = "\\.rds$", full.names = TRUE)
all <- lapply(files, readRDS)
scores <- do.call(rbind, lapply(all, `[[`, "scores"))
loci <- do.call(rbind, lapply(all, `[[`, "loci"))
write.table(scores, gzfile(file.path(out, paste0("scores_", stage, ".tsv.gz"))), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(loci, gzfile(file.path(out, paste0("loci_", stage, ".tsv.gz"))), sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(out, paste0("sessionInfo_", stage, ".txt")))
message(sprintf("LD given to the methods: %s; output %s", LDMODE, out))
message(sprintf("done: %d score rows from %d units; total %.1f min", nrow(scores), length(files),
                as.numeric(Sys.time() - t_start, units = "mins")))
