#!/usr/bin/env Rscript
# Reference LD from 1000 Genomes EUR (public data) for the 12 bulk windows.
#
#
# Usage (from the repository root; run once, takes about a minute):
#   Rscript scripts/bulk_h2ref_v1/prepare_ref_ld.R
# Options:
#   --bank  <dir>     windows (default data/sc-1mb-pattern-20260909-v1)
#   --plink <prefix>  1000G EUR chr22 PLINK files (default
#                     data/genotype_reference/1000G_EUR_chr22/1000G.EUR.QC.maf.22)
#   --out   <dir>     output (default data/ref_ld_1000g_eur_20261027)
#
# For every window:
#   1. Match each bank SNP (variants.tsv) to the reference by rsID and position.
#      Alleles must be identical (sign +1) or swapped (sign -1); anything else
#      is dropped. Palindromic SNPs were already removed from the bank.
#   2. Read the reference genotypes from the PLINK .bed file (dosage of the
#      .bim A1 allele), mean-impute missing calls, drop SNPs that are
#      monomorphic in the reference.
#   3. LD = Pearson correlation of the reference dosages, multiplied by
#      sign_i * sign_j so that it refers to the bank's counted allele a1
#      (the same allele as the simulated summary statistics).
#   4. Write <out>/<window>/ref_snps.tsv and <out>/<window>/ld_ref.f64
#      (float64, little-endian, row-major, m x m), and compare with the
#      in-sample GWAS LD of the same SNPs (summary.tsv). The comparison is a
#      check of allele alignment: the signs of |r| > 0.3 pairs must agree; on
#      the real windows the off-diagonal entries correlate at about 0.95-0.99.
# Only public 1000 Genomes genotypes are read here, plus the bank's SNP list
# and in-sample LD for the check. Nothing is uploaded.

.args <- commandArgs(trailingOnly = TRUE)
.opt <- function(name, default) { i <- match(paste0("--", name), .args)
  if (is.na(i) || i == length(.args)) default else .args[i + 1] }
script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
.find_repo <- function(start) { d <- normalizePath(start)
  repeat { if (file.exists(file.path(d, "R", "block_coloc.R"))) return(d)
    p <- dirname(d); if (p == d) stop("cannot find the repository root (folder with R/block_coloc.R)"); d <- p } }
repo <- .find_repo(dirname(script_file))
bank_root <- .opt("bank", file.path(repo, "data", "sc-1mb-pattern-20260909-v1"))
plink <- .opt("plink", file.path(repo, "data", "genotype_reference", "1000G_EUR_chr22", "1000G.EUR.QC.maf.22"))
out <- .opt("out", file.path(repo, "data", "ref_ld_1000g_eur_20261027"))
dir.create(out, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------ PLINK .bed ----

bim <- utils::read.table(paste0(plink, ".bim"), colClasses = c("character", "character", "numeric", "numeric", "character", "character"),
                         col.names = c("chr", "snp", "cm", "pos", "A1", "A2"))
n_ref <- length(readLines(paste0(plink, ".fam")))
nb <- (n_ref + 3L) %/% 4L
bed <- paste0(plink, ".bed")
con <- file(bed, "rb"); magic <- readBin(con, "raw", 3L); close(con)
if (!identical(as.integer(magic), c(0x6cL, 0x1bL, 0x01L))) stop(bed, " is not a SNP-major PLINK .bed file")
if (file.info(bed)$size != 3 + as.double(nb) * nrow(bim)) stop(bed, " size does not match .bim / .fam")

# Dosage of the .bim A1 allele: code 00 -> 2, 01 -> missing, 10 -> 1, 11 -> 0.
read_bed <- function(rows) {
  con <- file(bed, "rb"); on.exit(close(con))
  lut <- c(2, NA, 1, 0)
  vapply(rows, function(i) {
    seek(con, 3 + (i - 1) * as.double(nb))
    b <- as.integer(readBin(con, "raw", nb))
    codes <- as.vector(rbind(bitwAnd(b, 3L), bitwAnd(bitwShiftR(b, 2L), 3L),
                             bitwAnd(bitwShiftR(b, 4L), 3L), bitwAnd(bitwShiftR(b, 6L), 3L)))[seq_len(n_ref)]
    lut[codes + 1L]
  }, numeric(n_ref))
}

# ---------------------------------------------------------------- windows ----

windows <- sort(list.dirs(bank_root, full.names = FALSE, recursive = FALSE))
windows <- windows[grepl("^[0-9]+$", windows) & file.exists(file.path(bank_root, windows, "variants.tsv"))]
if (!length(windows)) stop("no windows found in ", bank_root)

summ <- list()
for (w in windows) {
  v <- utils::read.delim(file.path(bank_root, w, "variants.tsv"), stringsAsFactors = FALSE)
  p <- nrow(v)
  i <- match(v$snp, bim$snp)
  found <- !is.na(i) & bim$pos[ifelse(is.na(i), 1L, i)] == v$pos
  sgn <- rep(NA_real_, p)
  same <- found & bim$A1[ifelse(is.na(i), 1L, i)] == v$a1 & bim$A2[ifelse(is.na(i), 1L, i)] == v$a0
  swap <- found & bim$A1[ifelse(is.na(i), 1L, i)] == v$a0 & bim$A2[ifelse(is.na(i), 1L, i)] == v$a1
  sgn[same] <- 1; sgn[swap] <- -1
  k <- which(!is.na(sgn))
  X <- read_bed(i[k])
  miss <- colMeans(is.na(X))
  mu <- colMeans(X, na.rm = TRUE)
  for (j in seq_len(ncol(X))) X[is.na(X[, j]), j] <- mu[j]
  poly <- apply(X, 2, stats::sd) > 0
  k <- k[poly]; X <- X[, poly, drop = FALSE]; mu <- mu[poly]; miss <- miss[poly]
  R <- stats::cor(X) * outer(sgn[k], sgn[k])
  R <- (R + t(R)) / 2; diag(R) <- 1
  ref_af1 <- mu / 2                                # frequency of the .bim A1 allele
  info <- data.frame(snp = v$snp[k], bank_index = k, pos = v$pos[k], a1 = v$a1[k], a0 = v$a0[k],
                     ref_A1 = bim$A1[i[k]], ref_A2 = bim$A2[i[k]], sign = sgn[k],
                     ref_maf = round(pmin(ref_af1, 1 - ref_af1), 5), ref_missing = round(miss, 5))
  d <- file.path(out, w); dir.create(d, showWarnings = FALSE)
  utils::write.table(info, file.path(d, "ref_snps.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  con <- file(file.path(d, "ld_ref.f64"), "wb"); writeBin(as.vector(t(R)), con, size = 8L, endian = "little"); close(con)

  # check against the in-sample GWAS LD of the same SNPs
  chk <- c(cor_with_insample = NA, median_abs_diff = NA, p99_abs_diff = NA, sign_agree_r_gt_0.3 = NA)
  f <- file.path(bank_root, w, "ld_gwas.f64")
  if (file.exists(f)) {
    con <- file(f, "rb"); G <- matrix(readBin(con, "double", n = as.double(p) * p, size = 8L, endian = "little"), p, p, byrow = TRUE); close(con)
    G <- G[k, k]; ut <- upper.tri(R); a <- R[ut]; b <- G[ut]; big <- abs(b) > 0.3
    chk <- c(cor_with_insample = stats::cor(a, b), median_abs_diff = stats::median(abs(a - b)),
             p99_abs_diff = unname(stats::quantile(abs(a - b), 0.99)), sign_agree_r_gt_0.3 = mean(sign(a[big]) == sign(b[big])))
  }
  summ[[w]] <- data.frame(window = w, bank_snps = p, matched = length(k), not_in_reference = sum(!found),
                          allele_mismatch = sum(found & is.na(sgn)), monomorphic_in_reference = sum(!poly),
                          swapped_alleles = sum(sgn[k] < 0), n_reference = n_ref,
                          t(round(chk, 4)), check.names = FALSE)
  message(sprintf("window %s: %d of %d SNPs matched; cor with in-sample LD %.3f", w, length(k), p, chk[1]))
}
summ <- do.call(rbind, summ)
utils::write.table(summ, file.path(out, "summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(c("# 1000 Genomes EUR reference LD for the bulk windows",
             "",
             sprintf("Written by scripts/bulk_h2ref_v1/prepare_ref_ld.R on %s.", format(Sys.time(), "%Y-%m-%d %H:%M")),
             sprintf("Reference: %s (%d individuals, public 1000 Genomes Phase 3 EUR, GRCh37).", basename(plink), n_ref),
             "Per window: ref_snps.tsv (matched SNPs, bank order) and ld_ref.f64 (float64, little-endian, row-major, m x m),",
             "LD aligned to the bank's counted allele a1. summary.tsv compares it with the in-sample GWAS LD."),
           file.path(out, "README.md"))
print(summ, row.names = FALSE)
# A wrong allele alignment flips the sign of r for every pair involving that SNP,
# so strongly linked pairs must agree in sign. (On the real windows the
# correlation with the in-sample LD is about 0.95-0.99; it is lower on
# synthetic data, whose LD has less structure, so it is only a gross check here.)
bad <- summ$sign_agree_r_gt_0.3 < 0.99 | summ$cor_with_insample < 0.5
if (any(bad, na.rm = TRUE)) stop("allele alignment looks wrong in window(s): ", paste(summ$window[which(bad)], collapse = ", "))
message("done: ", out)
