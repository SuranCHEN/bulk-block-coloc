#!/usr/bin/env Rscript
# Checks of simulate_bulk.R on SYNTHETIC genotypes (no UKBB data needed).
#
# Author: Claude (written for the user's project).
#
# Usage (from the repository root):
#   Rscript scripts/bulk_gated_v1/test_synthetic.R [--make-bank <dir>]
# Without --make-bank it builds a small synthetic bank in a temporary folder,
# runs the checks below and prints PASS / FAIL for each. With --make-bank it
# also writes two synthetic windows (named like one calibration and one test
# window) to <dir>, so run.R can be exercised end to end with --bank <dir>.
#
# Checks
#   1. Bank loader: in-sample LD files equal cor() of the genotypes.
#   2. Causal sampling: right numbers of causal SNPs per trait, shared sets,
#      LD tier respected for the tier pair(s), all other pairs |r| < 0.2.
#   3. OLS statistics: z = beta / se; null z has variance about 1.
#   4. Null z across replicates are correlated as the in-sample LD (the
#      property z ~ N(R lambda, R) that every method relies on).
#   5. Causal |z| is about sqrt(N q) for a single causal SNP.
#   6. Common random numbers: same (window, scenario, rep) gives the same
#      causal SNPs and residuals; a different rep gives different ones.

.args <- commandArgs(trailingOnly = TRUE)
script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
# Repository root = the nearest parent folder that contains R/block_coloc.R.
# Works in the project layout (scripts/bulk_gated_v1/) and in the public
# repository layout (scripts/ with R/methods.R and R/simulate_bulk.R).
.find_repo <- function(start) { d <- normalizePath(start)
  repeat { if (file.exists(file.path(d, "R", "block_coloc.R"))) return(d)
    p <- dirname(d); if (p == d) stop("cannot find the repository root (folder with R/block_coloc.R)"); d <- p } }
.first_file <- function(...) { f <- c(...); f <- f[file.exists(f)]; if (!length(f)) stop("missing: ", paste(c(...), collapse = " or ")); f[1] }
repo <- .find_repo(dirname(script_file))
source(.first_file(file.path(repo, "R", "simulate_bulk.R"), file.path(repo, "scripts", "bulk_gated_v1", "simulate_bulk.R")))
i <- match("--make-bank", .args); bank_out <- if (is.na(i)) NULL else .args[i + 1]

results <- list()
check <- function(name, ok, detail = "") {
  results[[length(results) + 1]] <<- data.frame(check = name, result = if (isTRUE(ok)) "PASS" else "FAIL", detail = detail)
  message(sprintf("[%s] %s  %s", if (isTRUE(ok)) "PASS" else "FAIL", name, detail))
}

# ---------------------------------------------------- synthetic window ----

# Haplotypes from a thresholded Gaussian with block-wise AR(1) correlation,
# giving a mix of weak, moderate and strong LD.
make_window <- function(dir, center, p = 300L, n_gwas = 20000L, n_qtl = 1000L, seed = 1) {
  set.seed(seed)
  rho <- rep(c(0.5, 0.8, 0.95, 0.7, 0.9), length.out = p)
  af <- stats::runif(p, 0.06, 0.5)
  hap <- function(n) {
    Z <- matrix(stats::rnorm(n * p), n, p)
    for (j in 2:p) Z[, j] <- rho[j] * Z[, j - 1] + sqrt(1 - rho[j]^2) * Z[, j]
    sweep(Z, 2, stats::qnorm(1 - af), ">") * 1
  }
  geno <- function(n) hap(n) + hap(n)
  Gg <- geno(n_gwas); Gq <- geno(n_qtl)
  dir.create(file.path(dir, center), recursive = TRUE, showWarnings = FALSE)
  d <- file.path(dir, center)
  snp <- sprintf("rs%07d", seq_len(p))
  utils::write.table(data.frame(chrom = 22, snp = snp, cm = 0, pos = center - 500000 + 3000 * seq_len(p), a1 = "A", a0 = "G",
                                genotype_variant_index = seq_len(p)), file.path(d, "variants.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  wb <- function(M, f, size) { con <- file(file.path(d, f), "wb"); writeBin(as.vector(t(M)), con, size = size, endian = "little"); close(con) }
  wb(Gg, "gwas_genotypes.f32", 4L); wb(Gq, sprintf("qtl_genotypes_n%d.f32", n_qtl), 4L)
  wb(stats::cor(Gg), "ld_gwas.f64", 8L); wb(stats::cor(Gq), sprintf("ld_qtl_n%d.f64", n_qtl), 8L)
  invisible(d)
}

root <- if (is.null(bank_out)) file.path(tempdir(), "synbank") else bank_out
make_window(root, 31263287, seed = 11)   # same name as a calibration window
make_window(root, 21059630, seed = 22)   # same name as a test window
cfg <- bulk_default_config()

# 1. loader
bank <- tryCatch(bulk_load_bank(root, 31263287, cfg, check_ld = TRUE), error = function(e) e)
check("1 loader + QTL LD check", !inherits(bank, "error"), if (inherits(bank, "error")) conditionMessage(bank) else sprintf("p = %d", bank$p))
dg <- max(abs(stats::cor(bank$gwas$Gc) - bank$ld_gwas))
check("1 GWAS LD file equals cor(genotypes)", dg < 1e-6, sprintf("max diff %.1e", dg))

# 2. causal sampling
R <- abs(bank$ld_gwas)
tiers <- cfg$ld_tiers
bad <- character()
for (s in bulk_scenarios$scenario) for (t in names(tiers)) {
  sc <- bulk_scenarios[bulk_scenarios$scenario == s, ]
  if (!sc$has_tier_pair && t != "low") next
  for (r in 1:30) {
    cz <- tryCatch(bulk_sample_causal(bank, s, t, .bulk_seed("t", s, t, r), cfg), error = function(e) e)
    if (inherits(cz, "error")) { bad <- c(bad, paste(s, t, conditionMessage(cz))); next }
    if (length(cz$gwas) != sc$n_causal_gwas || length(cz$eqtl) != sc$n_causal_eqtl) bad <- c(bad, paste(s, t, "count"))
    exp_shared <- switch(s, A_A = 1, A_AB = 1, AB_AB = 2, A_B = 0, AB_CD = 0)
    if (length(cz$shared) != exp_shared) bad <- c(bad, paste(s, t, "shared"))
    pr <- abs(cz$pair_r)
    if (length(pr) && any(pr < tiers[[t]][1] | pr >= tiers[[t]][2])) bad <- c(bad, paste(s, t, "tier", paste(round(pr, 2), collapse = ",")))
    if (s == "AB_CD" && any(R[cz$gwas, cz$eqtl] >= cfg$r_other_max)) bad <- c(bad, paste(s, t, "cross pair"))
    if (any(c(bank$gwas$maf[unique(c(cz$gwas, cz$eqtl))], bank$qtl$maf[unique(c(cz$gwas, cz$eqtl))]) < cfg$maf_min)) bad <- c(bad, paste(s, t, "maf"))
  }
}
check("2 causal sampling (5 scenarios x tiers x 30 draws)", !length(bad), if (length(bad)) paste(head(unique(bad), 5), collapse = "; ") else "")

# 3. OLS z and null variance
dr <- bulk_draws(bank, 31263287, "A_B", 1)
cz0 <- list(gwas = integer(), eqtl = integer(), sign_gwas = numeric(), sign_eqtl = numeric())
gw0 <- bulk_gwas(bank, cz0, 0, dr); eq0 <- bulk_eqtl(bank, cz0, 0, dr)
check("3 z = beta / se", max(abs(gw0$z - gw0$beta / gw0$se)) < 1e-10)
fit <- stats::lm(I(dr$eps_eqtl - mean(dr$eps_eqtl)) ~ bank$qtl$Gc[, 5])
z_lm <- summary(fit)$coefficients[2, 3]
check("3 OLS z matches lm()", abs(z_lm - eq0$z[5]) < 1e-8, sprintf("lm %.4f vs %.4f", z_lm, eq0$z[5]))

# 4. null z correlation equals in-sample LD (eQTL sample, 2,000 null replicates)
Z <- t(vapply(1:2000, function(r) { set.seed(1000 + r)
  bulk_eqtl(bank, cz0, 0, list(eps_eqtl = stats::rnorm(bank$qtl$n)))$z }, numeric(bank$p)))
vz <- stats::median(apply(Z, 2, stats::var))
dz <- abs(stats::cor(Z) - bank$ld_qtl)
check("4 null z: median variance about 1", abs(vz - 1) < 0.08, sprintf("median var %.3f", vz))
check("4 null cor(z) equals in-sample LD", stats::quantile(dz, .99) < 0.08,
      sprintf("99th pct |diff| %.3f (Monte Carlo sd about %.3f)", stats::quantile(dz, .99), 1 / sqrt(2000)))

# 5. causal |z| about sqrt(N q), single causal SNP (A-A), 200 replicates
zs <- t(vapply(1:200, function(r) {
  cz <- bulk_sample_causal(bank, "A_A", "low", .bulk_seed("z", r), cfg)
  d <- bulk_draws(bank, 31263287, "A_A", r)
  c(abs(bulk_gwas(bank, cz, 0.003, d)$z[cz$gwas]), abs(bulk_eqtl(bank, cz, 0.015, d)$z[cz$eqtl]))
}, numeric(2)))
eg <- sqrt(20000 * 0.003); ee <- sqrt(1000 * 0.015)
check("5 causal |z| GWAS about sqrt(N q)", abs(mean(zs[, 1]) - eg) < 0.4, sprintf("mean %.2f vs %.2f", mean(zs[, 1]), eg))
check("5 causal |z| eQTL about sqrt(N q)", abs(mean(zs[, 2]) - ee) < 0.4, sprintf("mean %.2f vs %.2f", mean(zs[, 2]), ee))

# 6. common random numbers
a <- bulk_sample_causal(bank, "AB_CD", "low", .bulk_seed(1, "AB_CD", 3, "causal", "low"), cfg)
b <- bulk_sample_causal(bank, "AB_CD", "low", .bulk_seed(1, "AB_CD", 3, "causal", "low"), cfg)
c4 <- bulk_sample_causal(bank, "AB_CD", "low", .bulk_seed(1, "AB_CD", 4, "causal", "low"), cfg)
d1 <- bulk_draws(bank, 1, "AB_CD", 3); d2 <- bulk_draws(bank, 1, "AB_CD", 3); d3 <- bulk_draws(bank, 1, "AB_CD", 4)
check("6 same seed -> same causal SNPs and residuals", identical(a, b) && identical(d1, d2))
check("6 different rep -> different draws", !identical(d1$eps_gwas, d3$eps_gwas) && !identical(c(a$gwas, a$eqtl), c(c4$gwas, c4$eqtl)))

res <- do.call(rbind, results)
cat("\n"); print(res, row.names = FALSE)
if (any(res$result == "FAIL")) quit(status = 1)
