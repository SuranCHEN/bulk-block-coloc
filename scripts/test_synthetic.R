#!/usr/bin/env Rscript
# Checks of simulate_bulk.R and prepare_ref_ld.R on SYNTHETIC genotypes
# (no UKBB data and no 1000 Genomes data needed).
#
# Author: Claude (written for the user's project).
#
# Usage (from the repository root):
#   Rscript scripts/bulk_h2ref_v1/test_synthetic.R [--make-bank <dir>]
# Without --make-bank it builds a small synthetic bank and a synthetic
# reference panel in a temporary folder, runs the checks below and prints
# PASS / FAIL for each. With --make-bank it writes two synthetic windows to
# <dir> and their reference LD to <dir>_ref, so run.R can be exercised end to
# end with --bank <dir> --ref <dir>_ref.
#
# Checks
#   1. Bank loader: in-sample LD files equal cor() of the genotypes.
#   2. prepare_ref_ld.R: reads a PLINK .bed correctly, drops SNPs missing from
#      the reference, flips the sign of LD for swapped alleles (the reference
#      LD equals cor() of the reference dosages of the bank's counted allele).
#   3. Loader with the reference: SNP set restricted to matched SNPs; every
#      matrix and vector has the same SNPs in the same order.
#   4. Causal sampling: right numbers of causal SNPs per trait, shared sets,
#      LD tier respected for the tier pair(s), all other pairs |r| < 0.2.
#   5. OLS statistics: z = beta / se; null z has variance about 1; null z
#      across replicates are correlated as the in-sample LD.
#   6. Heritability: Var(g) equals h2 exactly; causal |z| is about
#      sqrt(N h2 / K) for K = 1 and K = 2 causal SNPs.
#   7. Locus assembly: ld = "ref" gives both traits the reference LD,
#      ld = "insample" gives each trait its own in-sample LD.
#   8. Common random numbers.
#   9. effect = "normal" (sensitivity analysis): E[Var(g)] about h2 and Var(g)
#      variable across replicates (coefficient of variation about sqrt(2) for
#      K = 1); same causal SNPs, signs and residuals as the fixed design;
#      magnitudes reproducible from the seed; the fixed design is unchanged.

.args <- commandArgs(trailingOnly = TRUE)
script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
# Repository root = the nearest parent folder that contains R/block_coloc.R.
# Works in the project layout (scripts/bulk_h2ref_v1/) and in the public
# repository layout (scripts/ with R/methods.R and R/simulate_bulk.R).
.find_repo <- function(start) { d <- normalizePath(start)
  repeat { if (file.exists(file.path(d, "R", "block_coloc.R"))) return(d)
    p <- dirname(d); if (p == d) stop("cannot find the repository root (folder with R/block_coloc.R)"); d <- p } }
.first_file <- function(...) { f <- c(...); f <- f[file.exists(f)]; if (!length(f)) stop("missing: ", paste(c(...), collapse = " or ")); f[1] }
repo <- .find_repo(dirname(script_file))
source(.first_file(file.path(repo, "R", "simulate_bulk.R"), file.path(repo, "scripts", "bulk_h2ref_v1", "simulate_bulk.R")))
source(.first_file(file.path(repo, "R", "methods.R"), file.path(repo, "scripts", "sc_gated_v2", "task01", "methods.R")))
prep_script <- .first_file(file.path(dirname(normalizePath(script_file)), "prepare_ref_ld.R"))
i <- match("--make-bank", .args); bank_out <- if (is.na(i)) NULL else .args[i + 1]

results <- list()
check <- function(name, ok, detail = "") {
  results[[length(results) + 1]] <<- data.frame(check = name, result = if (isTRUE(ok)) "PASS" else "FAIL", detail = detail)
  message(sprintf("[%s] %s  %s", if (isTRUE(ok)) "PASS" else "FAIL", name, detail))
}

# ---------------------------------------------------- synthetic window ----

# Haplotypes from a thresholded Gaussian with block-wise AR(1) correlation,
# giving a mix of weak, moderate and strong LD. The reference panel (489
# people) comes from the same model, so its LD differs from the in-sample LD
# by sampling noise only, as for a well-matched real reference.
make_window <- function(dir, center, p = 300L, n_gwas = 20000L, n_qtl = 1000L, n_ref = 489L, seed = 1) {
  set.seed(seed)
  rho <- rep(c(0.5, 0.8, 0.95, 0.7, 0.9), length.out = p)
  af <- stats::runif(p, 0.06, 0.5)
  hap <- function(n) {
    Z <- matrix(stats::rnorm(n * p), n, p)
    for (j in 2:p) Z[, j] <- rho[j] * Z[, j - 1] + sqrt(1 - rho[j]^2) * Z[, j]
    sweep(Z, 2, stats::qnorm(1 - af), ">") * 1
  }
  geno <- function(n) hap(n) + hap(n)
  Gg <- geno(n_gwas); Gq <- geno(n_qtl); Gr <- geno(n_ref)
  d <- file.path(dir, center); dir.create(d, recursive = TRUE, showWarnings = FALSE)
  snp <- sprintf("rs%d%04d", center %% 1000, seq_len(p))
  pos <- center - 500000 + 3000 * seq_len(p)
  utils::write.table(data.frame(chrom = 22, snp = snp, cm = 0, pos = pos, a1 = "A", a0 = "G",
                                genotype_variant_index = seq_len(p)), file.path(d, "variants.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  wb <- function(M, f, size) { con <- file(file.path(d, f), "wb"); writeBin(as.vector(t(M)), con, size = size, endian = "little"); close(con) }
  wb(Gg, "gwas_genotypes.f32", 4L); wb(Gq, sprintf("qtl_genotypes_n%d.f32", n_qtl), 4L)
  wb(stats::cor(Gg), "ld_gwas.f64", 8L); wb(stats::cor(Gq), sprintf("ld_qtl_n%d.f64", n_qtl), 8L)
  # reference: every 20th SNP is absent; every 3rd SNP has its alleles swapped
  # in the .bim (so the .bed stores the dosage of a0 = 2 - dosage of a1);
  # a few calls are missing.
  inref <- seq_len(p) %% 20 != 0
  swap <- seq_len(p) %% 3 == 0
  Gr[sample(length(Gr), 200)] <- NA
  list(snp = snp, pos = pos, Gr = Gr, inref = inref, swap = swap)
}

# Write a SNP-major PLINK 1 .bed/.bim/.fam from dosages of the .bim A1 allele.
write_plink <- function(prefix, chr, snp, pos, A1, A2, D) {
  n <- nrow(D); nb <- (n + 3L) %/% 4L
  code <- function(x) ifelse(is.na(x), 1L, ifelse(x == 2, 0L, ifelse(x == 1, 2L, 3L)))
  con <- file(paste0(prefix, ".bed"), "wb")
  writeBin(as.raw(c(0x6c, 0x1b, 0x01)), con)
  for (j in seq_len(ncol(D))) {
    cc <- c(code(D[, j]), rep(0L, 4L * nb - n))
    m <- matrix(cc, 4)
    writeBin(as.raw(m[1, ] + 4L * m[2, ] + 16L * m[3, ] + 64L * m[4, ]), con)
  }
  close(con)
  utils::write.table(data.frame(chr, snp, 0, pos, A1, A2), paste0(prefix, ".bim"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  utils::write.table(data.frame(seq_len(n), seq_len(n), 0, 0, 0, -9), paste0(prefix, ".fam"), sep = " ", quote = FALSE, row.names = FALSE, col.names = FALSE)
}

root <- if (is.null(bank_out)) file.path(tempdir(), "synbank") else bank_out
ref_root <- paste0(root, "_ref")
w1 <- make_window(root, 31263287, seed = 11)   # same name as a calibration window
w2 <- make_window(root, 21059630, seed = 22)   # same name as a test window
# one reference PLINK file for both windows (like a chromosome-wide panel)
dref <- list(); for (w in list(w1, w2)) {
  k <- which(w$inref)
  D <- w$Gr[, k]; D[, w$swap[k]] <- 2 - D[, w$swap[k]]
  dref[[length(dref) + 1]] <- list(snp = w$snp[k], pos = w$pos[k], A1 = ifelse(w$swap[k], "G", "A"), A2 = ifelse(w$swap[k], "A", "G"), D = D)
}
plink_dir <- file.path(dirname(root), paste0(basename(root), "_plink")); dir.create(plink_dir, showWarnings = FALSE)
prefix <- file.path(plink_dir, "synref")
write_plink(prefix, 22, unlist(lapply(dref, `[[`, "snp")), unlist(lapply(dref, `[[`, "pos")),
            unlist(lapply(dref, `[[`, "A1")), unlist(lapply(dref, `[[`, "A2")), do.call(cbind, lapply(dref, `[[`, "D")))
rs <- system2(file.path(R.home("bin"), "Rscript"), c(shQuote(prep_script), "--bank", shQuote(root), "--plink", shQuote(prefix), "--out", shQuote(ref_root)),
              stdout = TRUE, stderr = TRUE)
cfg <- bulk_default_config()

# 1. loader (no reference)
bank0 <- tryCatch(bulk_load_bank(root, 31263287, cfg, check_ld = TRUE), error = function(e) e)
check("1 loader + QTL LD check", !inherits(bank0, "error"), if (inherits(bank0, "error")) conditionMessage(bank0) else sprintf("p = %d", bank0$p))
dg <- max(abs(stats::cor(bank0$gwas$Gc) - bank0$ld_gwas))
check("1 GWAS LD file equals cor(genotypes)", dg < 1e-6, sprintf("max diff %.1e", dg))

# 2. prepare_ref_ld.R
ok_prep <- is.null(attr(rs, "status")) && file.exists(file.path(ref_root, "31263287", "ld_ref.f64"))
check("2 prepare_ref_ld.R ran", ok_prep, if (!ok_prep) paste(tail(rs, 5), collapse = " | ") else "")
info <- utils::read.delim(file.path(ref_root, "31263287", "ref_snps.tsv"), stringsAsFactors = FALSE)
check("2 SNPs absent from the reference are dropped", identical(info$bank_index, which(w1$inref)),
      sprintf("%d of %d kept", nrow(info), length(w1$inref)))
check("2 swapped alleles get sign -1", identical(info$sign == -1, w1$swap[w1$inref]))
Gr <- w1$Gr[, w1$inref]; for (j in seq_len(ncol(Gr))) Gr[is.na(Gr[, j]), j] <- mean(Gr[, j], na.rm = TRUE)
m <- nrow(info)
con <- file(file.path(ref_root, "31263287", "ld_ref.f64"), "rb"); Rr <- matrix(readBin(con, "double", m * m, 8L), m, m, byrow = TRUE); close(con)
dr <- max(abs(Rr - stats::cor(Gr)))
check("2 reference LD equals cor() of reference dosages of a1", dr < 1e-10, sprintf("max diff %.1e", dr))
sm <- utils::read.delim(file.path(ref_root, "summary.tsv"))
check("2 summary: strongly linked pairs agree in sign with the in-sample LD", all(sm$sign_agree_r_gt_0.3 == 1),
      sprintf("cor %s", paste(round(sm$cor_with_insample, 3), collapse = ", ")))

# 3. loader with the reference
bank <- tryCatch(bulk_load_bank(root, 31263287, cfg, check_ld = TRUE, ref_root = ref_root), error = function(e) e)
check("3 loader with reference", !inherits(bank, "error"), if (inherits(bank, "error")) conditionMessage(bank) else sprintf("p = %d of %d", bank$p, bank$p_bank))
same <- identical(rownames(bank$ld_ref), bank$snp) && identical(rownames(bank$ld_gwas), bank$snp) &&
  ncol(bank$gwas$Gc) == bank$p && ncol(bank$qtl$Gc) == bank$p && length(bank$gwas$maf) == bank$p &&
  max(abs(bank$ld_gwas - bank0$ld_gwas[w1$inref, w1$inref])) == 0
check("3 every matrix restricted to the same matched SNPs", isTRUE(same))

# 4. causal sampling
R <- abs(bank$ld_gwas)
tiers <- cfg$ld_tiers
bad <- character(); skipped <- character()
# double_shared1 needs an A with two partners in the tier. The synthetic
# genotypes have weaker LD than real data and may have none in the high tier
# (on the 12 real windows there are 400-1,300 such SNPs per window), so that
# combination is skipped here when it is impossible.
ab_ac_ok <- function(t) { T <- R >= tiers[[t]][1] & R < tiers[[t]][2]; diag(T) <- FALSE
  bc <- if (t == "low") cfg$r_other_max else tiers[[t]][2]
  any(vapply(seq_len(nrow(T)), function(a) { B <- which(T[a, ]); if (length(B) < 2) return(FALSE)
    S <- R[B, B] < bc; diag(S) <- FALSE; any(S) }, NA)) }
for (s in bulk_scenarios$scenario) for (t in names(tiers)) {
  sc <- bulk_scenarios[bulk_scenarios$scenario == s, ]
  if (!sc$has_tier_pair && t != "low") next
  if (s == "AB_AC" && !ab_ac_ok(t)) { skipped <- c(skipped, paste(s, t)); next }
  for (r in 1:30) {
    cz <- tryCatch(bulk_sample_causal(bank, s, t, .bulk_seed("t", s, t, r), cfg), error = function(e) e)
    if (inherits(cz, "error")) { bad <- c(bad, paste(s, t, conditionMessage(cz))); next }
    if (length(cz$gwas) != sc$n_causal_gwas || length(cz$eqtl) != sc$n_causal_eqtl) bad <- c(bad, paste(s, t, "count"))
    exp_shared <- switch(s, A_A = 1, A_AB = 1, AB_AB = 2, AB_AC = 1, A_B = 0, AB_CD = 0, gwas_only = 0, null_both = 0)
    if (length(cz$shared) != exp_shared) bad <- c(bad, paste(s, t, "shared"))
    pr <- abs(cz$pair_r)
    if (length(pr) && any(pr < tiers[[t]][1] | pr >= tiers[[t]][2])) bad <- c(bad, paste(s, t, "tier", paste(round(pr, 2), collapse = ",")))
    if (s == "AB_CD" && any(R[cz$gwas, cz$eqtl] >= cfg$r_other_max)) bad <- c(bad, paste(s, t, "cross pair"))
    if (s == "AB_AC") {
      if (cz$gwas[1] != cz$eqtl[1] || cz$gwas[2] == cz$eqtl[2]) bad <- c(bad, paste(s, t, "structure"))
      bc_max <- if (t == "low") cfg$r_other_max else tiers[[t]][2]
      if (abs(cz$cross_r) >= bc_max || abs(abs(cz$cross_r) - R[cz$gwas[2], cz$eqtl[2]]) > 1e-12) bad <- c(bad, paste(s, t, "B-C"))
    }
    if (any(c(bank$gwas$maf[unique(c(cz$gwas, cz$eqtl))], bank$qtl$maf[unique(c(cz$gwas, cz$eqtl))]) < cfg$maf_min)) bad <- c(bad, paste(s, t, "maf"))
  }
}
check("4 causal sampling (8 scenarios x tiers x 30 draws)", !length(bad),
      paste(c(if (length(bad)) paste(head(unique(bad), 5), collapse = "; "),
              if (length(skipped)) paste("skipped (impossible in synthetic LD):", paste(skipped, collapse = ", "))), collapse = " "))

# 5. OLS z, null variance, null correlation = in-sample LD
dr0 <- bulk_draws(bank, 31263287, "A_B", 1)
cz0 <- list(gwas = integer(), eqtl = integer(), sign_gwas = numeric(), sign_eqtl = numeric())
gw0 <- bulk_gwas(bank, cz0, 0.003, dr0); eq0 <- bulk_eqtl(bank, cz0, 0.015, dr0)
check("5 z = beta / se", max(abs(gw0$z - gw0$beta / gw0$se)) < 1e-10)
check("5 no causal SNP -> Var(g) = 0", gw0$var_g == 0 && eq0$var_g == 0)
fit <- stats::lm(I(dr0$eps_eqtl - mean(dr0$eps_eqtl)) ~ bank$qtl$Gc[, 5])
z_lm <- summary(fit)$coefficients[2, 3]
check("5 OLS z matches lm()", abs(z_lm - eq0$z[5]) < 1e-8, sprintf("lm %.4f vs %.4f", z_lm, eq0$z[5]))
Z <- t(vapply(1:2000, function(r) { set.seed(1000 + r)
  bulk_eqtl(bank, cz0, 0, list(eps_eqtl = stats::rnorm(bank$qtl$n)))$z }, numeric(bank$p)))
vz <- stats::median(apply(Z, 2, stats::var))
dz <- abs(stats::cor(Z) - bank$ld_qtl)
check("5 null z: median variance about 1", abs(vz - 1) < 0.08, sprintf("median var %.3f", vz))
check("5 null cor(z) equals in-sample LD", stats::quantile(dz, .99) < 0.08,
      sprintf("99th pct |diff| %.3f (Monte Carlo sd about %.3f)", stats::quantile(dz, .99), 1 / sqrt(2000)))

# 6. heritability: Var(g) = h2; causal |z| about sqrt(N h2 / K)
zs <- t(vapply(1:200, function(r) {
  c1 <- bulk_sample_causal(bank, "A_A", "low", .bulk_seed("z1", r), cfg)
  c2 <- bulk_sample_causal(bank, "AB_AB", "low", .bulk_seed("z2", r), cfg)
  d <- bulk_draws(bank, 31263287, "A_A", r)
  g1 <- bulk_gwas(bank, c1, 0.003, d); e1 <- bulk_eqtl(bank, c1, 0.015, d)
  g2 <- bulk_gwas(bank, c2, 0.003, d); e2 <- bulk_eqtl(bank, c2, 0.015, d)
  c(abs(g1$z[c1$gwas]), abs(e1$z[c1$eqtl]), mean(abs(g2$z[c2$gwas])), mean(abs(e2$z[c2$eqtl])),
    g1$var_g, e1$var_g, g2$var_g, e2$var_g, g2$var_y, e2$var_y)
}, numeric(10)))
vg_err <- max(abs(zs[, 5:8] - rep(c(0.003, 0.015, 0.003, 0.015), each = nrow(zs))))
check("6 Var(g) equals h2 exactly (K = 1 and K = 2)", vg_err < 1e-12, sprintf("max |Var(g) - h2| %.1e", vg_err))
check("6 Var(y) about 1", abs(mean(zs[, 9:10]) - 1) < 0.02, sprintf("mean %.3f", mean(zs[, 9:10])))
ex <- c(sqrt(20000 * 0.003), sqrt(1000 * 0.015), sqrt(20000 * 0.003 / 2), sqrt(1000 * 0.015 / 2))
obs <- colMeans(zs[, 1:4])
check("6 causal |z| about sqrt(N h2 / K): GWAS K=1, eQTL K=1, GWAS K=2, eQTL K=2", all(abs(obs - ex) < 0.45),
      sprintf("observed %s vs expected %s", paste(round(obs, 2), collapse = "/"), paste(round(ex, 2), collapse = "/")))

# 7. locus assembly
cz <- bulk_sample_causal(bank, "A_AB", "low", .bulk_seed("loc", 1), cfg)
d <- bulk_draws(bank, 31263287, "A_AB", 1)
gw <- bulk_gwas(bank, cz, 0.003, d); eq <- bulk_eqtl(bank, cz, 0.015, d)
lr <- bulk_assemble_locus(bank, cz, gw, eq, ld = "ref")
li <- bulk_assemble_locus(bank, cz, gw, eq, ld = "insample")
k <- lr$sumstats$snp
check("7 ld = ref: both traits get the reference LD",
      max(abs(lr$ld1 - bank$ld_ref[k, k])) < 1e-12 && max(abs(lr$ld2 - bank$ld_ref[k, k])) < 1e-12)
check("7 ld = insample: each trait gets its own in-sample LD",
      max(abs(li$ld1 - bank$ld_gwas[k, k])) < 1e-12 && max(abs(li$ld2 - bank$ld_qtl[k, k])) < 1e-12)

# 8. common random numbers
a <- bulk_sample_causal(bank, "AB_CD", "low", .bulk_seed(1, "AB_CD", 3, "causal", "low"), cfg)
b <- bulk_sample_causal(bank, "AB_CD", "low", .bulk_seed(1, "AB_CD", 3, "causal", "low"), cfg)
c4 <- bulk_sample_causal(bank, "AB_CD", "low", .bulk_seed(1, "AB_CD", 4, "causal", "low"), cfg)
d1 <- bulk_draws(bank, 1, "AB_CD", 3); d2 <- bulk_draws(bank, 1, "AB_CD", 3); d3 <- bulk_draws(bank, 1, "AB_CD", 4)
check("8 same seed -> same causal SNPs and residuals", identical(a, b) && identical(d1, d2))
check("8 different rep -> different draws", !identical(d1$eps_gwas, d3$eps_gwas) && !identical(c(a$gwas, a$eqtl), c(c4$gwas, c4$eqtl)))

# 9. effect = "normal"
nv <- t(vapply(1:2000, function(r) {
  c1 <- bulk_sample_causal(bank, "A_A", "low", .bulk_seed("n1", r), cfg)
  c2 <- bulk_sample_causal(bank, "AB_AB", "low", .bulk_seed("n2", r), cfg)
  c1 <- c(c1, bulk_effect_magnitudes(1, "A_A", r, "low", c1)); c2 <- c(c2, bulk_effect_magnitudes(1, "AB_AB", r, "low", c2))
  d <- bulk_draws(bank, 1, "A_A", r)
  c(bulk_eqtl(bank, c1, 0.015, d)$var_g, bulk_eqtl(bank, c2, 0.015, d)$var_g, bulk_gwas(bank, c2, 0.003, d)$var_g)
}, numeric(3)))
mv <- colMeans(nv) / c(0.015, 0.015, 0.003); cv1 <- stats::sd(nv[, 1]) / mean(nv[, 1])
check("9 normal effects: mean Var(g) / h2 about 1 (eQTL K=1, eQTL K=2, GWAS K=2)", all(abs(mv - 1) < 0.1),
      sprintf("%s", paste(round(mv, 3), collapse = "/")))
check("9 normal effects: Var(g) varies, CV about sqrt(2) for K = 1", abs(cv1 - sqrt(2)) < 0.15, sprintf("CV %.2f", cv1))
cz9 <- bulk_sample_causal(bank, "AB_AB", "low", .bulk_seed("n9", 1), cfg)
m1 <- bulk_effect_magnitudes(1, "AB_AB", 1, "low", cz9); m2 <- bulk_effect_magnitudes(1, "AB_AB", 1, "low", cz9)
m3 <- bulk_effect_magnitudes(1, "AB_AB", 2, "low", cz9)
check("9 normal effects: magnitudes reproducible, differ across reps", identical(m1, m2) && !identical(m1, m3))
d9 <- bulk_draws(bank, 1, "AB_AB", 1)
f9 <- bulk_eqtl(bank, cz9, 0.015, d9); n9 <- bulk_eqtl(bank, c(cz9, m1), 0.015, d9)
f9b <- .bulk_trait(bank$qtl, cz9$eqtl, cz9$sign_eqtl, 0.015, d9$eps_eqtl)
check("9 fixed design unchanged; normal design differs only through effect sizes",
      identical(f9, f9b) && !identical(f9$z, n9$z) && abs(n9$var_y - (n9$var_g + 1 - 0.015)) < 0.1,
      sprintf("Var(y) %.3f", n9$var_y))

res <- do.call(rbind, results)
cat("\n"); print(res, row.names = FALSE)
if (!is.null(bank_out)) message("synthetic bank: ", root, "\nreference LD: ", ref_root,
                                "\nrun.R: --bank ", root, " --ref ", ref_root)
if (any(res$result == "FAIL")) quit(status = 1)
