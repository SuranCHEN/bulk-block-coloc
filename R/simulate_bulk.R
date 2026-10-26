# Bulk simulation library: one GWAS trait and one bulk eQTL trait, both from
# a LINEAR model, on real UKBB chr22 genotypes.
#
# Author: Claude (written for the user's project; executed by Codex).
# This file only DEFINES functions. Nothing runs when it is sourced.
# It does not use any single-cell code.
#
# ---------------------------------------------------------------------------
# Data (read-only, never copied): data/sc-1mb-pattern-20260909-v1/<window>/
#   gwas_genotypes.f32      20,000 GWAS individuals x p SNPs (dosage)
#   qtl_genotypes_n1000.f32  1,000 QTL individuals x p SNPs (dosage)
#   ld_gwas.f64             in-sample LD of the 20,000 GWAS individuals
#   ld_qtl_n1000.f64        in-sample LD of the 1,000 QTL individuals
#   variants.tsv, reference_scale.tsv
# The folder name contains "sc" for historical reasons only; here we use
# nothing but the genotypes and LD. GWAS and QTL individuals do not overlap
# (checked in donors.tsv).
#
# Generative model, the same for both traits (trait 1 = GWAS, trait 2 = eQTL):
#   y_k = sum_j b_j x_kj + e_k,  x_kj = genotype standardised in that sample,
#   b_j = s_j * sqrt(q)  for each causal SNP j of that trait, s_j = +1/-1,
#   e_k ~ N(0, 1 - Var(g)), so Var(y) is about 1 and EACH causal SNP explains
#   a fraction q of the trait variance. Two causal SNPs in one trait get the
#   same q (1:1). Expected causal |z| is about sqrt(N q).
#   Summary statistics: per-allele OLS of y on each SNP's dosage.
#   Matched LD for each trait is the in-sample LD of that trait's individuals
#   (exact for OLS statistics; ld_qtl_n1000 equals cor() of the QTL dosages).
#
# Scenarios (Wallace 2021, PLoS Genet 17:e1009440, Fig 1, plus AB-CD).
# Letters are causal SNPs; left of "-" = GWAS, right = eQTL.
#   A_A    A-A    GWAS A,    eQTL A      shared       truth 1
#   A_AB   A-AB   GWAS A,    eQTL A,B    A shared     truth 1
#   AB_AB  AB-AB  GWAS A,B,  eQTL A,B    both shared  truth 1
#   A_B    A-B    GWAS A,    eQTL B      not shared   truth 0
#   AB_CD  AB-CD  GWAS A,B,  eQTL C,D    not shared   truth 0
# Extra negatives (added 2026-10-25; reported for the false-positive rate only,
# never used to calibrate the PP4 cut-offs):
#   gwas_only     GWAS A,    eQTL none                truth 0
#   null_both     GWAS none, eQTL none                truth 0
# Report names (evaluate.R): A_A single_shared, A_AB single_double_shared,
# AB_AB double_shared2, A_B single_distinct, AB_CD double_distinct.
#
# LD between causal SNPs (absolute r in the GWAS in-sample LD):
#   "tier pair(s)": A-B in A_AB, AB_AB, A_B; A-B and C-D in AB_CD.
#     low  : |r| <  0.2            (main setting, as in the earlier bulk work)
#     mid  : 0.32 <= |r| < 0.63    (r^2 between the clump 0.1 and block 0.4 cut-offs)
#     high : 0.63 <= |r| < 0.9     (r^2 >= 0.4: the two SNPs fall in one block)
#   every other pair of different causal SNPs: |r| < 0.2.
#   A_A has no pair and is not run for the mid / high tiers.
#
# Common random numbers: for a given (window, scenario, rep) the residual
# draws of each trait are fixed, and for a given tier so are the causal SNPs
# and signs; only q changes across strength levels.
# ---------------------------------------------------------------------------

bulk_default_config <- function() {
  list(
    n_gwas = 20000L,          # GWAS individuals in the bank
    n_qtl = 1000L,            # QTL individuals (bank files exist for 100, 200, 1000)
    q_gwas = 0.003,           # variance explained per GWAS causal SNP (main 0.3%)
    q_eqtl = 0.015,           # variance explained per eQTL causal SNP (main 1.5%)
    maf_min = 0.05,           # causal SNPs: in-sample MAF >= 0.05 in both samples
    r_other_max = 0.2,        # |r| for causal pairs outside the tier pair(s)
    ld_tiers = list(low = c(0, 0.2), mid = c(0.32, 0.63), high = c(0.63, 0.9)),
    max_draw_attempts = 500L
  )
}

bulk_scenarios <- data.frame(
  scenario = c("A_A", "A_AB", "AB_AB", "A_B", "AB_CD", "gwas_only", "null_both"),
  label = c("A-A", "A-AB", "AB-AB", "A-B", "AB-CD", "gwas_only", "null_both"),
  truth = c(1L, 1L, 1L, 0L, 0L, 0L, 0L),
  n_causal_gwas = c(1L, 1L, 2L, 1L, 2L, 1L, 0L),
  n_causal_eqtl = c(1L, 2L, 2L, 1L, 2L, 0L, 0L),
  has_tier_pair = c(FALSE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE),
  calibration_negative = c(FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE),
  stringsAsFactors = FALSE
)

# ------------------------------------------------------------- bank I/O ----

.bulk_read_matrix <- function(path, n, p, size) {
  expected <- as.double(n) * as.double(p) * as.double(size)
  if (!file.exists(path) || file.info(path)$size != expected)
    stop("matrix size mismatch: ", path, " (expected ", expected, " bytes)")
  con <- file(path, "rb"); on.exit(close(con))
  v <- readBin(con, "double", n = as.double(n) * p, size = size, endian = "little")
  matrix(v, nrow = n, byrow = TRUE)
}

# Centred genotypes and their column sums of squares for one sample.
.bulk_sample <- function(G) {
  af <- colMeans(G) / 2
  Gc <- sweep(G, 2, 2 * af, "-")
  ss <- colSums(Gc^2)
  list(Gc = Gc, ss = ss, sd = sqrt(ss / (nrow(G) - 1)), maf = pmin(af, 1 - af), n = nrow(G))
}

#' Load one 1-Mb window. Genotypes stay on this computer; nothing is copied.
#' check_ld = TRUE recomputes the QTL in-sample LD and stops if the stored
#' file differs (cheap; the GWAS check is left to the smoke test).
bulk_load_bank <- function(bank_root, center, cfg = bulk_default_config(), check_ld = TRUE) {
  dir <- file.path(bank_root, format(center, scientific = FALSE, trim = TRUE))
  rd <- function(f) utils::read.delim(file.path(dir, f), stringsAsFactors = FALSE)
  variants <- rd("variants.tsv")
  if (anyDuplicated(variants$snp)) stop("duplicated SNP ids in ", dir)
  snp <- variants$snp; p <- length(snp)
  ld <- function(f) {
    R <- .bulk_read_matrix(file.path(dir, f), p, p, 8L)
    if (max(abs(R - t(R))) > 1e-8 || max(abs(diag(R) - 1)) > 1e-6) stop("invalid LD: ", f)
    dimnames(R) <- list(snp, snp); R
  }
  gw <- .bulk_sample(.bulk_read_matrix(file.path(dir, "gwas_genotypes.f32"), cfg$n_gwas, p, 4L))
  qt <- .bulk_sample(.bulk_read_matrix(file.path(dir, sprintf("qtl_genotypes_n%d.f32", cfg$n_qtl)), cfg$n_qtl, p, 4L))
  ld_gwas <- ld("ld_gwas.f64")
  ld_qtl <- ld(sprintf("ld_qtl_n%d.f64", cfg$n_qtl))
  if (check_ld) {
    ok <- qt$ss > 1e-8
    C <- stats::cor(qt$Gc[, ok, drop = FALSE])
    d <- max(abs(C - ld_qtl[ok, ok]))
    if (d > 1e-6) stop("ld_qtl file differs from the in-sample QTL LD by ", d)
  }
  list(center = center, snp = snp, p = p, pos = variants$pos,
       gwas = gw, qtl = qt, ld_gwas = ld_gwas, ld_qtl = ld_qtl)
}

# ------------------------------------------------------------------ seeds ----

.bulk_seed <- function(...) {
  v <- utf8ToInt(paste(..., sep = "|")); h <- 7
  for (x in v) h <- (h * 131 + x) %% 2147483629
  as.integer(h)
}

# ------------------------------------------------------- causal sampling ----

#' Draw causal SNPs for one scenario and LD tier. Returns indices for
#' GWAS / eQTL, signs, the shared set and the realised |r| of the tier pair(s).
bulk_sample_causal <- function(bank, scenario, tier = "low", seed, cfg = bulk_default_config()) {
  sc <- bulk_scenarios[bulk_scenarios$scenario == scenario, ]
  if (!nrow(sc)) stop("unknown scenario: ", scenario)
  if (!tier %in% names(cfg$ld_tiers)) stop("unknown LD tier: ", tier)
  if (!sc$has_tier_pair && tier != "low") stop(scenario, " has no causal pair; use tier 'low'")
  set.seed(seed)
  R <- abs(bank$ld_gwas)
  elig <- which(bank$gwas$maf >= cfg$maf_min & bank$qtl$maf >= cfg$maf_min)
  lo <- cfg$ld_tiers[[tier]][1]; hi <- cfg$ld_tiers[[tier]][2]
  in_tier <- function(cand, a) R[cand, a] >= lo & R[cand, a] < hi
  far <- function(cand, chosen) if (!length(chosen)) rep(TRUE, length(cand)) else
    apply(R[cand, chosen, drop = FALSE] < cfg$r_other_max, 1, all)
  pick <- function(cand) if (length(cand)) cand[sample.int(length(cand), 1L)] else NA_integer_
  if (scenario == "null_both") { g <- integer(); e <- integer() } else
  for (attempt in seq_len(cfg$max_draw_attempts)) {
    A <- pick(elig)
    if (scenario == "A_A") { g <- A; e <- A; break }
    if (scenario == "gwas_only") { g <- A; e <- integer(); break }
    others <- setdiff(elig, A)
    B <- pick(others[in_tier(others, A)])
    if (is.na(B)) next
    if (scenario %in% c("A_AB", "AB_AB", "A_B")) {
      g <- switch(scenario, A_AB = A, AB_AB = c(A, B), A_B = A)
      e <- switch(scenario, A_AB = c(A, B), AB_AB = c(A, B), A_B = B)
      break
    }
    # AB_CD: C and D form a tier pair; every cross pair has |r| < r_other_max
    rest <- setdiff(elig, c(A, B))
    C <- pick(rest[far(rest, c(A, B))])
    if (is.na(C)) next
    rest <- setdiff(rest, C)
    D <- pick(rest[in_tier(rest, C) & far(rest, c(A, B))])
    if (is.na(D)) next
    g <- c(A, B); e <- c(C, D)
    break
  }
  if (!exists("g", inherits = FALSE))
    stop("could not draw causal SNPs for ", scenario, " tier ", tier)
  pair_r <- switch(scenario, A_A = numeric(), gwas_only = numeric(), null_both = numeric(), A_AB = bank$ld_gwas[e[1], e[2]],
                   AB_AB = bank$ld_gwas[g[1], g[2]], A_B = bank$ld_gwas[g, e],
                   AB_CD = c(bank$ld_gwas[g[1], g[2]], bank$ld_gwas[e[1], e[2]]))
  list(scenario = scenario, tier = tier, gwas = g, eqtl = e,
       sign_gwas = sample(c(-1, 1), length(g), replace = TRUE),
       sign_eqtl = sample(c(-1, 1), length(e), replace = TRUE),
       shared = intersect(g, e), pair_r = unname(pair_r))
}

# ----------------------------------------------- common random numbers ----

bulk_draws <- function(bank, center, scenario, rep) {
  set.seed(.bulk_seed(center, scenario, rep, "gwas")); e1 <- rnorm(bank$gwas$n)
  set.seed(.bulk_seed(center, scenario, rep, "eqtl")); e2 <- rnorm(bank$qtl$n)
  list(eps_gwas = e1, eps_eqtl = e2)
}

# ------------------------------------------------------ linear trait ----

# One linear trait on one sample. Returns OLS summary statistics for every SNP.
.bulk_trait <- function(S, causal, signs, q, eps) {
  n <- S$n
  g <- numeric(n)
  if (length(causal)) {
    X <- sweep(S$Gc[, causal, drop = FALSE], 2, S$sd[causal], "/")
    g <- drop(X %*% (signs * sqrt(q)))
  }
  ve <- max(1 - stats::var(g), 0.05)
  y <- g + sqrt(ve) * eps
  y <- y - mean(y)
  Xty <- drop(crossprod(S$Gc, y))
  beta <- Xty / S$ss
  rss <- pmax(sum(y^2) - beta * Xty, 1e-12)
  se <- sqrt(rss / (n - 2) / S$ss)
  z <- beta / se
  valid <- is.finite(z) & S$ss > 1e-8
  list(beta = beta, se = se, z = z, p = 2 * stats::pt(-abs(z), df = n - 2),
       valid = valid, n = n, q = q, var_y = stats::var(y), var_g = stats::var(g))
}

bulk_gwas <- function(bank, causal, q, draws)
  .bulk_trait(bank$gwas, causal$gwas, causal$sign_gwas, q, draws$eps_gwas)

bulk_eqtl <- function(bank, causal, q, draws)
  .bulk_trait(bank$qtl, causal$eqtl, causal$sign_eqtl, q, draws$eps_eqtl)

# ----------------------------------------------------------- assemble ----

#' Build a methods.R locus (trait 1 = GWAS, trait 2 = eQTL). Each trait gets
#' its own in-sample MAF and in-sample LD. Truth labels are attached.
bulk_assemble_locus <- function(bank, causal, gw, eq) {
  keep <- which(gw$valid & eq$valid)
  ss <- data.frame(snp = bank$snp[keep],
                   beta_trait1 = gw$beta[keep], se_trait1 = gw$se[keep], pval_trait1 = gw$p[keep],
                   beta_trait2 = eq$beta[keep], se_trait2 = eq$se[keep], pval_trait2 = eq$p[keep],
                   maf_trait1 = bank$gwas$maf[keep], maf_trait2 = bank$qtl$maf[keep],
                   stringsAsFactors = FALSE)
  loc <- new_locus(ss, n1 = gw$n, n2 = eq$n,
                   ld1 = bank$ld_gwas[keep, keep], ld2 = bank$ld_qtl[keep, keep])
  loc$truth <- list(scenario = causal$scenario, tier = causal$tier,
                    truth = bulk_scenarios$truth[bulk_scenarios$scenario == causal$scenario],
                    causal_gwas = bank$snp[causal$gwas], causal_eqtl = bank$snp[causal$eqtl],
                    shared = bank$snp[causal$shared])
  loc
}
