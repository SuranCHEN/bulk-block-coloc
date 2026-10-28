# Bulk simulation library: two linear traits (trait 1 = large sample, called
# "GWAS"; trait 2 = small sample, called "eQTL") on real UKBB chr22 genotypes.
# Methods receive the 1000 Genomes EUR reference LD (default) or in-sample LD.
#
# Author: Claude (written for the user's project; executed by Codex).
# This file only DEFINES functions. Nothing runs when it is sourced.
# It does not use any single-cell code.
#
# Changes from scripts/bulk_gated_v1/simulate_bulk.R (user decisions 2026-10-27):
#   1. Strength is the trait's local heritability h2 in the window, split
#      equally over its K causal SNPs (per-SNP h2 = h2 / K), as in the
#      earlier (August) bulk work. Before, every causal SNP explained a fixed q.
#   2. The genetic value is rescaled so that its in-sample variance is exactly
#      h2 (it differs from K * h2/K only through LD between causal SNPs).
#   3. Methods get the 1000 Genomes EUR LD (489 people) instead of the in-sample
#      LD, as in a real summary-statistics analysis. The SNP set of every window
#      is restricted to SNPs present in the reference, for both LD choices, so a
#      comparison of the two LD choices differs only in the LD matrix.
#   4. (2026-10-28, sensitivity analysis requested by the user) effect = "normal":
#      the effect of each causal SNP is drawn from N(0, h2 / K), independently
#      for the two traits, as in traceCB (Jiang et al. 2026). The genetic value
#      is NOT rescaled, so Var(g) equals h2 only in expectation, and
#      e ~ N(0, 1 - h2). The draw is b_j = s_j * |u_j| * sqrt(h2 / K) with the
#      same causal SNPs, signs s_j and residuals as the "fixed" design and
#      u_j ~ N(0, 1) from a separate seed, so the two designs are paired.
#      effect = "fixed" (default) is unchanged: the results of stages 1 and 2
#      are reproduced exactly.
#
# ---------------------------------------------------------------------------
# Data (read-only, never copied): data/sc-1mb-pattern-20260909-v1/<window>/
#   gwas_genotypes.f32      20,000 GWAS individuals x p SNPs (dosage of a1)
#   qtl_genotypes_n1000.f32  1,000 QTL individuals x p SNPs (dosage of a1)
#   ld_gwas.f64             in-sample LD of the 20,000 GWAS individuals
#   ld_qtl_n1000.f64        in-sample LD of the 1,000 QTL individuals
#   variants.tsv            snp, pos, a1 (counted allele), a0
# Reference LD (written by prepare_ref_ld.R from public 1000 Genomes data):
#   data/ref_ld_1000g_eur_20261027/<window>/ref_snps.tsv   matched SNPs
#   data/ref_ld_1000g_eur_20261027/<window>/ld_ref.f64     LD, aligned to a1
# GWAS and QTL individuals do not overlap (checked in donors.tsv).
#
# Generative model, the same for both traits:
#   y_k = sum_j b_j x_kj + e_k,  x_kj = genotype standardised in that sample,
#   b_j = s_j * sqrt(h2 / K) for the K causal SNPs of that trait, s_j = +1/-1,
#   g = sum_j b_j x_j is then rescaled so that Var(g) = h2 in the sample,
#   e_k ~ N(0, 1 - h2), so Var(y) is about 1.
#   (effect = "normal": b_j = s_j |u_j| sqrt(h2 / K), u_j ~ N(0, 1), no rescaling.)
#   Expected causal |z| is about sqrt(N h2 / K) when the causal SNPs are in
#   low LD. Summary statistics: per-allele OLS of y on each SNP's dosage.
#
# Scenarios (Wallace 2021, PLoS Genet 17:e1009440, Fig 1, plus AB-CD).
# Letters are causal SNPs; left of "-" = GWAS, right = eQTL.
#   A_A    GWAS A,    eQTL A      shared       truth 1   (report: single_shared)
#   A_AB   GWAS A,    eQTL A,B    A shared     truth 1   (single_double_shared)
#   AB_AB  GWAS A,B,  eQTL A,B    both shared  truth 1   (double_shared2)
#   AB_AC  GWAS A,B,  eQTL A,C    A shared     truth 1   (double_shared1; added 2026-10-27)
#   A_B    GWAS A,    eQTL B      not shared   truth 0   (single_distinct)
#   AB_CD  GWAS A,B,  eQTL C,D    not shared   truth 0   (double_distinct)
# Extra negatives (false-positive rate only, never used for calibration):
#   gwas_only     GWAS A,    eQTL none                truth 0
#   null_both     GWAS none, eQTL none                truth 0
#
# LD between causal SNPs (absolute r in the GWAS in-sample LD, i.e. the true LD):
#   "tier pair(s)": A-B in A_AB, AB_AB, A_B; A-B and C-D in AB_CD;
#     A-B (GWAS pair) and A-C (eQTL pair) in AB_AC.
#     low  : |r| <  0.2            (main setting)
#     mid  : 0.32 <= |r| < 0.63
#     high : 0.63 <= |r| < 0.9
#   every other pair of different causal SNPs: |r| < 0.2, except B-C in AB_AC
#   in the mid / high tiers: B and C are both linked to A there, so |r(B,C)|
#   < 0.2 is often impossible; B-C is then only required to be below the
#   tier's upper bound (realised value reported as cross_r).
#
# Common random numbers: for a given (window, scenario, rep) the residual
# draws of each trait are fixed, and for a given tier so are the causal SNPs
# and signs; only h2 changes across strength levels.
# ---------------------------------------------------------------------------

bulk_default_config <- function() {
  list(
    n_gwas = 20000L,          # GWAS individuals in the bank
    n_qtl = 1000L,            # QTL individuals (bank files exist for 100, 200, 1000)
    h2_gwas = 0.003,          # local h2 of the GWAS trait, split over its causal SNPs (main 0.3%)
    h2_eqtl = 0.015,          # local h2 of the eQTL trait, split over its causal SNPs (main 1.5%)
    maf_min = 0.05,           # causal SNPs: in-sample MAF >= 0.05 in both samples
    r_other_max = 0.2,        # |r| for causal pairs outside the tier pair(s)
    ld_tiers = list(low = c(0, 0.2), mid = c(0.32, 0.63), high = c(0.63, 0.9)),
    max_draw_attempts = 500L,
    effect = "fixed"          # "fixed" (main design) or "normal" (sensitivity analysis, see header)
  )
}

bulk_scenarios <- data.frame(
  scenario = c("A_A", "A_AB", "AB_AB", "AB_AC", "A_B", "AB_CD", "gwas_only", "null_both"),
  label = c("A-A", "A-AB", "AB-AB", "AB-AC", "A-B", "AB-CD", "gwas_only", "null_both"),
  truth = c(1L, 1L, 1L, 1L, 0L, 0L, 0L, 0L),
  n_causal_gwas = c(1L, 1L, 2L, 2L, 1L, 2L, 1L, 0L),
  n_causal_eqtl = c(1L, 2L, 2L, 2L, 1L, 2L, 0L, 0L),
  has_tier_pair = c(FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE),
  calibration_negative = c(FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE),
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

.bulk_subset <- function(S, keep) list(Gc = S$Gc[, keep, drop = FALSE], ss = S$ss[keep], sd = S$sd[keep],
                                       maf = S$maf[keep], n = S$n)

#' Reference LD of one window (written by prepare_ref_ld.R). Returns the bank
#' column indices of the matched SNPs and the LD among them, already aligned
#' to the bank's counted allele a1.
bulk_read_ref <- function(ref_root, center, variants) {
  dir <- file.path(ref_root, format(center, scientific = FALSE, trim = TRUE))
  rs <- utils::read.delim(file.path(dir, "ref_snps.tsv"), stringsAsFactors = FALSE)
  if (!all(variants$snp[rs$bank_index] == rs$snp)) stop("ref_snps.tsv does not match variants.tsv in ", dir)
  if (is.unsorted(rs$bank_index, strictly = TRUE)) stop("ref_snps.tsv must follow the bank order in ", dir)
  m <- nrow(rs)
  R <- .bulk_read_matrix(file.path(dir, "ld_ref.f64"), m, m, 8L)
  if (max(abs(R - t(R))) > 1e-8 || max(abs(diag(R) - 1)) > 1e-6) stop("invalid reference LD in ", dir)
  dimnames(R) <- list(rs$snp, rs$snp)
  list(index = rs$bank_index, ld = R, info = rs)
}

#' Load one 1-Mb window. Genotypes stay on this computer; nothing is copied.
#' With ref_root, every SNP not in the reference panel is dropped (from the
#' genotypes, the in-sample LD and the reference LD alike).
#' check_ld = TRUE recomputes the QTL in-sample LD and stops if the stored
#' file differs.
bulk_load_bank <- function(bank_root, center, cfg = bulk_default_config(), check_ld = TRUE, ref_root = NULL) {
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
  ld_ref <- NULL; keep <- seq_len(p)
  if (!is.null(ref_root)) {
    ref <- bulk_read_ref(ref_root, center, variants)
    keep <- ref$index; ld_ref <- ref$ld
    gw <- .bulk_subset(gw, keep); qt <- .bulk_subset(qt, keep)
    ld_gwas <- ld_gwas[keep, keep]; ld_qtl <- ld_qtl[keep, keep]
  }
  list(center = center, snp = snp[keep], p = length(keep), p_bank = p, pos = variants$pos[keep],
       gwas = gw, qtl = qt, ld_gwas = ld_gwas, ld_qtl = ld_qtl, ld_ref = ld_ref)
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
    if (scenario == "AB_AC") {
      # B (GWAS pair) and C (eQTL pair) both in the tier with A; B-C below
      # r_other_max (low tier) or below the tier's upper bound (mid / high),
      # see the header. Every B in the tier is tried before re-drawing A.
      bc_max <- if (tier == "low") cfg$r_other_max else hi
      Bs <- others[in_tier(others, A)]
      if (length(Bs) < 2) next
      Bs <- Bs[sample.int(length(Bs))]
      C <- NA_integer_
      for (B in Bs) {
        cand <- setdiff(Bs, B)
        C <- pick(cand[R[cand, B] < bc_max])
        if (!is.na(C)) break
      }
      if (is.na(C)) next
      g <- c(A, B); e <- c(A, C)
      break
    }
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
                   AB_AC = c(bank$ld_gwas[g[1], g[2]], bank$ld_gwas[e[1], e[2]]),
                   AB_CD = c(bank$ld_gwas[g[1], g[2]], bank$ld_gwas[e[1], e[2]]))
  list(scenario = scenario, tier = tier, gwas = g, eqtl = e,
       sign_gwas = sample(c(-1, 1), length(g), replace = TRUE),
       sign_eqtl = sample(c(-1, 1), length(e), replace = TRUE),
       shared = intersect(g, e), pair_r = unname(pair_r),
       cross_r = if (scenario == "AB_AC") unname(bank$ld_gwas[g[2], e[2]]) else numeric())
}

# ----------------------------------------------- common random numbers ----

bulk_draws <- function(bank, center, scenario, rep) {
  set.seed(.bulk_seed(center, scenario, rep, "gwas")); e1 <- rnorm(bank$gwas$n)
  set.seed(.bulk_seed(center, scenario, rep, "eqtl")); e2 <- rnorm(bank$qtl$n)
  list(eps_gwas = e1, eps_eqtl = e2)
}

#' Effect magnitudes for effect = "normal": |u_j|, u_j ~ N(0, 1), one per
#' causal SNP of each trait, drawn independently for the two traits (a shared
#' SNP gets two independent magnitudes). Seeded by (window, scenario, rep,
#' tier) only, so they do not change with h2 (common random numbers).
bulk_effect_magnitudes <- function(center, scenario, rep, tier, causal) {
  set.seed(.bulk_seed(center, scenario, rep, "effect", tier))
  list(mag_gwas = abs(stats::rnorm(length(causal$gwas))), mag_eqtl = abs(stats::rnorm(length(causal$eqtl))))
}

# ------------------------------------------------------ linear trait ----

# One linear trait on one sample with local heritability h2 split equally over
# its K causal SNPs. Returns OLS summary statistics for every SNP.
# mags = NULL: fixed effects sqrt(h2 / K), Var(g) rescaled to h2 exactly.
# mags given (effect = "normal"): effects signs * mags * sqrt(h2 / K), no
# rescaling, residual variance 1 - h2 (traceCB).
.bulk_trait <- function(S, causal, signs, h2, eps, mags = NULL) {
  n <- S$n; K <- length(causal)
  if (h2 < 0 || h2 >= 0.95) stop("h2 must be in [0, 0.95)")
  if (!is.null(mags) && length(mags) != K) stop("mags must have one value per causal SNP")
  g <- numeric(n)
  if (K && h2 > 0) {
    X <- sweep(S$Gc[, causal, drop = FALSE], 2, S$sd[causal], "/")
    if (is.null(mags)) {
      g <- drop(X %*% (signs * sqrt(h2 / K)))
      vg <- stats::var(g)
      if (vg > 1e-12) g <- g * sqrt(h2 / vg)     # Var(g) = h2 exactly in this sample
    } else {
      g <- drop(X %*% (signs * mags * sqrt(h2 / K)))   # Var(g) = h2 in expectation only
    }
  }
  ve <- if (is.null(mags)) 1 - stats::var(g) else 1 - h2
  y <- g + sqrt(ve) * eps
  y <- y - mean(y)
  Xty <- drop(crossprod(S$Gc, y))
  beta <- Xty / S$ss
  rss <- pmax(sum(y^2) - beta * Xty, 1e-12)
  se <- sqrt(rss / (n - 2) / S$ss)
  z <- beta / se
  valid <- is.finite(z) & S$ss > 1e-8
  list(beta = beta, se = se, z = z, p = 2 * stats::pt(-abs(z), df = n - 2),
       valid = valid, n = n, h2 = h2, h2_per_snp = if (K) h2 / K else 0,
       var_y = stats::var(y), var_g = stats::var(g))
}

# causal$mag_gwas / causal$mag_eqtl are set only for effect = "normal"
# (bulk_effect_magnitudes); otherwise they are NULL and the fixed design is used.
bulk_gwas <- function(bank, causal, h2, draws)
  .bulk_trait(bank$gwas, causal$gwas, causal$sign_gwas, h2, draws$eps_gwas, causal$mag_gwas)

bulk_eqtl <- function(bank, causal, h2, draws)
  .bulk_trait(bank$qtl, causal$eqtl, causal$sign_eqtl, h2, draws$eps_eqtl, causal$mag_eqtl)

# ----------------------------------------------------------- assemble ----

#' Build a methods.R locus (trait 1 = GWAS, trait 2 = eQTL). Each trait gets
#' its own in-sample MAF. LD given to the methods:
#'   ld = "ref"      both traits use the 1000 Genomes EUR reference LD (default)
#'   ld = "insample" each trait uses its own in-sample LD
#' Truth labels are attached.
bulk_assemble_locus <- function(bank, causal, gw, eq, ld = c("ref", "insample")) {
  ld <- match.arg(ld)
  keep <- which(gw$valid & eq$valid)
  ss <- data.frame(snp = bank$snp[keep],
                   beta_trait1 = gw$beta[keep], se_trait1 = gw$se[keep], pval_trait1 = gw$p[keep],
                   beta_trait2 = eq$beta[keep], se_trait2 = eq$se[keep], pval_trait2 = eq$p[keep],
                   maf_trait1 = bank$gwas$maf[keep], maf_trait2 = bank$qtl$maf[keep],
                   stringsAsFactors = FALSE)
  if (ld == "ref") {
    if (is.null(bank$ld_ref)) stop("reference LD not loaded (bulk_load_bank(..., ref_root = ...))")
    ld1 <- ld2 <- bank$ld_ref[keep, keep]
  } else {
    ld1 <- bank$ld_gwas[keep, keep]; ld2 <- bank$ld_qtl[keep, keep]
  }
  loc <- new_locus(ss, n1 = gw$n, n2 = eq$n, ld1 = ld1, ld2 = ld2)
  loc$truth <- list(scenario = causal$scenario, tier = causal$tier,
                    truth = bulk_scenarios$truth[bulk_scenarios$scenario == causal$scenario],
                    causal_gwas = bank$snp[causal$gwas], causal_eqtl = bank$snp[causal$eqtl],
                    shared = bank$snp[causal$shared], ld = ld)
  loc
}
