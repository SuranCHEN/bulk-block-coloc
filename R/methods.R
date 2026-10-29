# Gated block-coloc: method definitions for the v2 comparison.
#
# This file only DEFINES methods. It does not simulate data and does not
# modify R/block_coloc.R, whose internal helpers are reused unchanged.
#
# One locus is passed as a list (see new_locus()):
#   sumstats : data.frame, one row per SNP, columns
#              snp, beta_trait1, se_trait1, pval_trait1,
#              beta_trait2, se_trait2, pval_trait2, and maf
#              (or maf_trait1 + maf_trait2). Trait 1 = GWAS, trait 2 = eQTL.
#   n1, n2   : sample sizes (GWAS individuals; eQTL donors).
#   ld1, ld2 : signed LD correlation matrices matched to trait 1 / trait 2
#              statistics, rows and columns named by SNP.
#   susie1, susie2 : optional pre-fitted susieR objects. When absent, SuSiE is
#              fitted with susie_rss(z, R, n). The count-model simulation will
#              pass its own fits, because NB score statistics need the
#              weighted-likelihood SuSiE fit rather than plain z + LD.
#   sdY1, sdY2 : optional; passed to coloc as the trait standard deviation.
#
# Methods (names are fixed and used in every report):
#   V0_handoff          original block_coloc(): top-2/top-5 lead fill-in and
#                       max(block, whole) when either trait has <= 1 lead.
#   V1_block_any        block branch if BOTH traits have >= 1 lead passing the
#                       threshold; otherwise whole-locus coloc.abf.
#   V2_gate_marginal    block branch only if both traits have >= 1 lead AND at
#                       least one trait has >= 2 leads (marginal clumping).
#   V3_gate_conditional as V2, but a second or later lead counts only if it is
#                       still significant after conditioning on the leads
#                       already kept; blocks are built from the kept leads.
#   V4_cond_block       PROPOSED (not yet agreed): V3 gating and blocks, but
#                       inside each block every SNP's z is first conditioned
#                       on the same trait's kept leads OUTSIDE that block, so
#                       shadows of out-of-block signals are removed before
#                       coloc.abf. Close in spirit to coloc.signals(cond).
#   whole_abf           coloc.abf on all SNPs.
#   susie_raw           coloc.susie; NA when either trait has no credible set.
#   susie_fallback      coloc.susie; whole-locus coloc.abf when either trait
#                       has no credible set.
#   signals_mask        coloc.signals(method = "mask").
#   signals_cond        coloc.signals(method = "cond").
#
# Every method returns one row with the same columns (see .gated_row()).
# For the block branch with no valid block, pp4 is NA (status
# "no_valid_block"); the whole-locus PP4 is always stored in pp4_whole so the
# two policies ("NA" vs "fall back to whole") can be compared later without
# rerunning anything.

gated_default_config <- function() {
  list(
    p_threshold1 = 1e-3,       # lead threshold, GWAS (user decision 2026-10-10)
    p_threshold2 = 1e-3,       # lead threshold, eQTL
    clump_r2 = 0.10,           # handoff default
    block_r2 = 0.40,           # handoff default
    sign_split = TRUE,         # handoff default
    min_component_snps = 2L,   # handoff default
    p1 = 1e-4, p2 = 1e-4, p12 = 1e-5,  # coloc priors, identical for all methods
    cond_ridge = 1e-6,         # numerical ridge for R_SS inverse in conditioning
    susie_L = 10L, susie_coverage = 0.95, susie_purity = 0.5,  # susieR defaults
    signals_maxhits = 3L,      # coloc.signals default
    signals_r2thr = 0.01,      # coloc.signals default
    signals_pthr = 1e-3        # set equal to the lead threshold (coloc default 1e-6)
  )
}

gated_methods <- c(
  "V0_handoff", "V1_block_any", "V2_gate_marginal", "V3_gate_conditional",
  "V4_cond_block", "whole_abf", "susie_raw", "susie_fallback", "signals_mask", "signals_cond"
)

# ---------------------------------------------------------------- inputs ----

new_locus <- function(sumstats, n1, n2, ld1, ld2 = ld1,
                      susie1 = NULL, susie2 = NULL, sdY1 = NULL, sdY2 = NULL) {
  sumstats <- as.data.frame(sumstats, stringsAsFactors = FALSE)
  need <- c("snp", "beta_trait1", "se_trait1", "pval_trait1",
            "beta_trait2", "se_trait2", "pval_trait2")
  miss <- setdiff(need, names(sumstats))
  if (length(miss)) stop("sumstats is missing: ", paste(miss, collapse = ", "))
  if (!"maf_trait1" %in% names(sumstats)) sumstats$maf_trait1 <- sumstats$maf
  if (!"maf_trait2" %in% names(sumstats)) sumstats$maf_trait2 <- sumstats$maf
  if (is.null(sumstats$maf_trait1) || is.null(sumstats$maf_trait2))
    stop("sumstats needs maf, or maf_trait1 and maf_trait2")
  sumstats$snp <- as.character(sumstats$snp)
  if (anyDuplicated(sumstats$snp)) stop("duplicated SNP ids")
  for (k in c("se_trait1", "se_trait2"))
    if (any(!is.finite(sumstats[[k]]) | sumstats[[k]] <= 0)) stop(k, " must be finite and > 0")
  align <- function(LD, label) {
    LD <- as.matrix(LD); storage.mode(LD) <- "double"
    if (is.null(rownames(LD)) || is.null(colnames(LD))) stop(label, " needs SNP row and column names")
    if (!all(sumstats$snp %in% rownames(LD)) || !all(sumstats$snp %in% colnames(LD)))
      stop(label, " does not contain every SNP")
    LD <- LD[sumstats$snp, sumstats$snp, drop = FALSE]
    if (any(!is.finite(LD))) stop(label, " has non-finite values")
    LD <- (LD + t(LD)) / 2; diag(LD) <- 1
    LD[LD > 1] <- 1; LD[LD < -1] <- -1
    LD
  }
  n1 <- as.numeric(n1); n2 <- as.numeric(n2)
  if (!is.finite(n1) || n1 <= 0 || !is.finite(n2) || n2 <= 0) stop("n1, n2 must be > 0")
  list(sumstats = sumstats, n1 = n1, n2 = n2,
       ld1 = align(ld1, "ld1"), ld2 = align(ld2, "ld2"),
       susie1 = susie1, susie2 = susie2, sdY1 = sdY1, sdY2 = sdY2)
}

.z <- function(locus, trait) {
  ss <- locus$sumstats
  ss[[paste0("beta_trait", trait)]] / ss[[paste0("se_trait", trait)]]
}

# ------------------------------------------------------- lead selection ----

# Marginal leads: SNPs with p < threshold, greedily clumped at r2 < clump_r2
# using the trait's own LD. No fill-in: a trait with no SNP below the
# threshold has zero leads.
.leads_marginal <- function(locus, trait, config) {
  p <- locus$sumstats[[paste0("pval_trait", trait)]]
  LD <- locus[[paste0("ld", trait)]]
  thr <- config[[paste0("p_threshold", trait)]]
  cand <- which(is.finite(p) & p < thr)
  if (!length(cand)) return(integer())
  .block_coloc_greedy_clump(cand, p, LD, config$clump_r2)
}

# Conditional z of SNP j given the kept set S:
#   z_j|S = (z_j - R_jS R_SS^-1 z_S) / sqrt(1 - R_jS R_SS^-1 R_Sj)
.conditional_z <- function(z, LD, j, S, ridge) {
  if (!length(S)) return(z[j])
  RSS <- LD[S, S, drop = FALSE] + diag(ridge, length(S))
  RjS <- LD[j, S, drop = FALSE]
  w <- RjS %*% solve(RSS)
  num <- z[j] - drop(w %*% z[S])
  den2 <- 1 - drop(w %*% t(RjS))
  if (!is.finite(den2) || den2 <= 1e-8) return(0)
  num / sqrt(den2)
}

# Conditional leads: walk the marginal leads in p-value order; the first is
# kept; each later one is kept only if its conditional p-value (given all
# leads kept so far) is still below the threshold.
.leads_conditional <- function(locus, trait, config, marginal = NULL) {
  if (is.null(marginal)) marginal <- .leads_marginal(locus, trait, config)
  if (length(marginal) <= 1L) return(list(kept = marginal, z_cond = .z(locus, trait)[marginal]))
  z <- .z(locus, trait)
  LD <- locus[[paste0("ld", trait)]]
  thr <- config[[paste0("p_threshold", trait)]]
  kept <- marginal[1L]; zc <- z[marginal[1L]]
  for (j in marginal[-1L]) {
    zj <- .conditional_z(z, LD, j, kept, config$cond_ridge)
    if (2 * pnorm(-abs(zj)) < thr) { kept <- c(kept, j); zc <- c(zc, zj) }
  }
  list(kept = kept, z_cond = zc)
}

# Conditional z for many SNPs at once, given a set S (same formula as above).
.conditional_z_vec <- function(z, LD, idx, S, ridge) {
  if (!length(S)) return(z[idx])
  RSSi <- solve(LD[S, S, drop = FALSE] + diag(ridge, length(S)))
  RiS <- LD[idx, S, drop = FALSE]
  W <- RiS %*% RSSi
  num <- z[idx] - drop(W %*% z[S])
  den2 <- 1 - rowSums(W * RiS)
  out <- num / sqrt(pmax(den2, 1e-8))
  out[den2 <= 1e-8] <- 0   # SNP is (almost) a copy of a conditioning lead
  out
}

# --------------------------------------------------- block construction ----

# Same rules as .block_coloc_build_components() in R/block_coloc.R, except that
# each trait's neighbourhoods use that trait's own LD (ld1 for GWAS leads, ld2
# for eQTL leads). With ld1 == ld2 the result is identical to the original.
.build_blocks <- function(locus, leads1, leads2, config) {
  if (!length(leads1) || !length(leads2)) return(list())
  win1 <- lapply(leads1, function(i) which(locus$ld1[, i]^2 >= config$block_r2))
  win2 <- lapply(leads2, function(i) which(locus$ld2[, i]^2 >= config$block_r2))
  overlap <- matrix(FALSE, length(leads1), length(leads2))
  for (a in seq_along(leads1)) for (b in seq_along(leads2))
    overlap[a, b] <- length(intersect(win1[[a]], win2[[b]])) > 0L
  comps <- .block_coloc_graph_components(overlap)
  s1 <- sign(.z(locus, 1)[leads1]); s2 <- sign(.z(locus, 2)[leads2])
  blocks <- list()
  add <- function(rows, cols, reason) {
    idx <- sort(unique(c(unlist(win1[rows]), unlist(win2[cols]))))
    if (length(idx) >= config$min_component_snps)
      blocks[[length(blocks) + 1L]] <<- list(indices = idx, leads1 = leads1[rows],
                                              leads2 = leads2[cols], reason = reason)
  }
  for (cp in comps) {
    mixed <- length(unique(s1[cp$trait1])) > 1L || length(unique(s2[cp$trait2])) > 1L
    if (isTRUE(config$sign_split) && mixed) {
      edges <- which(overlap[cp$trait1, cp$trait2, drop = FALSE], arr.ind = TRUE)
      for (e in seq_len(nrow(edges)))
        add(cp$trait1[edges[e, 1]], cp$trait2[edges[e, 2]], "mixed_sign_pair_split")
    } else add(cp$trait1, cp$trait2, "connected_component")
  }
  if (length(blocks)) {
    keys <- vapply(blocks, function(b) paste(b$indices, collapse = ","), "")
    blocks <- blocks[!duplicated(keys)]
  }
  blocks
}

# ------------------------------------------------------------- coloc.abf ----

.coloc_dataset <- function(locus, idx, trait, with_ld = FALSE) {
  ss <- locus$sumstats
  d <- list(snp = ss$snp[idx],
            beta = ss[[paste0("beta_trait", trait)]][idx],
            varbeta = ss[[paste0("se_trait", trait)]][idx]^2,
            MAF = ss[[paste0("maf_trait", trait)]][idx],
            N = locus[[paste0("n", trait)]], type = "quant")
  sdY <- locus[[paste0("sdY", trait)]]
  if (!is.null(sdY)) d$sdY <- sdY
  if (with_ld) d$LD <- locus[[paste0("ld", trait)]][idx, idx, drop = FALSE]
  d
}

.abf <- function(locus, idx, config) {
  fit <- NULL
  invisible(utils::capture.output(fit <- suppressWarnings(suppressMessages(
    coloc::coloc.abf(.coloc_dataset(locus, idx, 1), .coloc_dataset(locus, idx, 2),
                     p1 = config$p1, p2 = config$p2, p12 = config$p12)))))
  s <- fit$summary
  list(pp3 = unname(s["PP.H3.abf"]), pp4 = unname(s["PP.H4.abf"]))
}

.whole <- function(locus, config) .abf(locus, seq_len(nrow(locus$sumstats)), config)

# ------------------------------------------------------------ result row ----

.gated_row <- function(method, branch = NA_character_, status = "ok",
                       pp4 = NA_real_, pp4_whole = NA_real_,
                       leads1 = integer(), leads2 = integer(),
                       n_leads1_marginal = NA_integer_, n_leads2_marginal = NA_integer_,
                       n_blocks = NA_integer_, best_block = integer(), locus = NULL,
                       detail = NA_character_) {
  snp <- if (is.null(locus)) character() else locus$sumstats$snp
  data.frame(
    method = method, branch = branch, status = status,
    pp4 = pp4, pp4_whole = pp4_whole,
    n_leads1 = length(leads1), n_leads2 = length(leads2),
    n_leads1_marginal = n_leads1_marginal, n_leads2_marginal = n_leads2_marginal,
    leads1 = paste(snp[leads1], collapse = ";"),
    leads2 = paste(snp[leads2], collapse = ";"),
    n_blocks = n_blocks,
    best_block_snps = paste(snp[best_block], collapse = ";"),
    detail = detail, stringsAsFactors = FALSE
  )
}

# V4: replace each block SNP's beta by (conditional z) x se, conditioning each
# trait on its kept leads that are not part of this block.
.condition_outside <- function(locus, block, l1, l2, config) {
  idx <- block$indices
  for (t in 1:2) {
    out <- setdiff(if (t == 1) l1 else l2, if (t == 1) block$leads1 else block$leads2)
    out <- setdiff(out, idx)
    if (!length(out)) next
    zc <- .conditional_z_vec(.z(locus, t), locus[[paste0("ld", t)]], idx, out, config$cond_ridge)
    se <- locus$sumstats[[paste0("se_trait", t)]][idx]
    locus$sumstats[[paste0("beta_trait", t)]][idx] <- zc * se
  }
  locus
}

# -------------------------------------------------------- block versions ----

.run_block_version <- function(method, locus, config, whole) {
  m1 <- .leads_marginal(locus, 1, config)
  m2 <- .leads_marginal(locus, 2, config)
  if (method %in% c("V3_gate_conditional", "V4_cond_block")) {
    l1 <- .leads_conditional(locus, 1, config, m1)$kept
    l2 <- .leads_conditional(locus, 2, config, m2)$kept
  } else { l1 <- m1; l2 <- m2 }
  go_block <- length(l1) >= 1L && length(l2) >= 1L &&
    (method == "V1_block_any" || max(length(l1), length(l2)) >= 2L)
  common <- list(leads1 = l1, leads2 = l2, n_leads1_marginal = length(m1),
                 n_leads2_marginal = length(m2), locus = locus, pp4_whole = whole$pp4)
  if (!go_block)
    return(do.call(.gated_row, c(list(method = method, branch = "whole",
                                      pp4 = whole$pp4, n_blocks = 0L), common)))
  blocks <- .build_blocks(locus, l1, l2, config)
  pp <- vapply(blocks, function(b) {
    target <- if (method == "V4_cond_block") .condition_outside(locus, b, l1, l2, config) else locus
    f <- tryCatch(.abf(target, b$indices, config), error = function(e) NULL)
    if (is.null(f)) NA_real_ else f$pp4
  }, numeric(1))
  if (!length(pp) || all(!is.finite(pp)))
    return(do.call(.gated_row, c(list(method = method, branch = "block",
                                      status = "no_valid_block", n_blocks = length(blocks)), common)))
  k <- which.max(pp)
  do.call(.gated_row, c(list(method = method, branch = "block", pp4 = pp[k],
                             n_blocks = length(blocks), best_block = blocks[[k]]$indices,
                             detail = blocks[[k]]$reason), common))
}

.run_v0 <- function(locus, config, whole) {
  cfg <- utils::modifyList(block_coloc_default_config(), list(
    trait1_p_threshold = config$p_threshold1, trait2_p_threshold = config$p_threshold2,
    lead_clump_r2 = config$clump_r2, block_r2 = config$block_r2,
    sign_split = config$sign_split, hybrid_whole_locus = TRUE,
    min_component_snps = config$min_component_snps,
    p1 = config$p1, p2 = config$p2, p12 = config$p12))
  fit <- block_coloc(locus$sumstats, locus$ld1, locus$n1, locus$n2, config = cfg)
  s <- fit$summary
  snp <- locus$sumstats$snp
  l1 <- match(fit$leads$snp[fit$leads$trait == "trait1"], snp)
  l2 <- match(fit$leads$snp[fit$leads$trait == "trait2"], snp)
  best <- if (is.finite(s$selected_component_id))
    match(fit$membership$snp[fit$membership$component_id == s$selected_component_id], snp) else integer()
  .gated_row("V0_handoff",
             branch = switch(s$selected_source, whole_locus_single_lead_hybrid = "whole",
                             maximum_block_pp4 = "block", no_overlapping_valid_block = "block"),
             status = if (s$selected_source == "no_overlapping_valid_block") "no_valid_block" else "ok",
             pp4 = if (s$selected_source == "no_overlapping_valid_block") NA_real_ else s$pp4,
             pp4_whole = whole$pp4, leads1 = l1, leads2 = l2,
             n_blocks = s$n_components, best_block = best, locus = locus,
             detail = paste0("fallback1=", s$trait1_used_fallback, ";fallback2=", s$trait2_used_fallback,
                             ";source=", s$selected_source, ";uses_ld1_only"))
}

# ----------------------------------------------------------------- SuSiE ----

.susie_fit <- function(locus, trait, config) {
  pre <- locus[[paste0("susie", trait)]]
  if (!is.null(pre)) return(pre)
  susieR::susie_rss(z = .z(locus, trait), R = locus[[paste0("ld", trait)]],
                    n = locus[[paste0("n", trait)]], L = config$susie_L,
                    coverage = config$susie_coverage, min_abs_corr = config$susie_purity)
}

.run_susie <- function(locus, config, whole) {
  f1 <- tryCatch(.susie_fit(locus, 1, config), error = function(e) e)
  f2 <- tryCatch(.susie_fit(locus, 2, config), error = function(e) e)
  bad <- inherits(f1, "error") || inherits(f2, "error")
  ncs <- function(f) if (inherits(f, "error")) NA_integer_ else length(f$sets$cs)
  n1 <- ncs(f1); n2 <- ncs(f2)
  status <- if (bad) "susie_error" else if (n1 == 0L && n2 == 0L) "no_cs_both" else
    if (n1 == 0L) "no_cs_trait1" else if (n2 == 0L) "no_cs_trait2" else "ok"
  pp <- NA_real_
  if (status == "ok") {
    cs <- tryCatch(suppressWarnings(suppressMessages(coloc::coloc.susie(
      f1, f2, p1 = config$p1, p2 = config$p2, p12 = config$p12))), error = function(e) e)
    if (inherits(cs, "error")) status <- "coloc_susie_error" else
      pp <- max(cs$summary$PP.H4.abf, na.rm = TRUE)
  }
  det <- paste0("n_cs1=", n1, ";n_cs2=", n2)
  raw <- .gated_row("susie_raw", branch = "susie", status = status, pp4 = pp,
                    pp4_whole = whole$pp4, detail = det)
  fb_from_whole <- !(status %in% c("ok"))
  fb <- .gated_row("susie_fallback", branch = if (fb_from_whole) "whole" else "susie",
                   status = if (fb_from_whole) paste0("fallback_whole:", status) else "ok",
                   pp4 = if (fb_from_whole) whole$pp4 else pp, pp4_whole = whole$pp4, detail = det)
  rbind(raw, fb)
}

# --------------------------------------------------------- coloc.signals ----

.run_signals <- function(locus, config, whole, mode) {
  idx <- seq_len(nrow(locus$sumstats))
  res <- tryCatch(suppressWarnings(suppressMessages(coloc::coloc.signals(
    .coloc_dataset(locus, idx, 1, with_ld = TRUE), .coloc_dataset(locus, idx, 2, with_ld = TRUE),
    method = mode, p1 = config$p1, p2 = config$p2, p12 = config$p12,
    maxhits = config$signals_maxhits, r2thr = config$signals_r2thr,
    pthr = config$signals_pthr))), error = function(e) e)
  name <- paste0("signals_", mode)
  if (inherits(res, "error"))
    return(.gated_row(name, branch = "signals", status = paste0("error:", conditionMessage(res)),
                      pp4_whole = whole$pp4))
  s <- as.data.frame(res$summary)
  .gated_row(name, branch = "signals", pp4 = max(s$PP.H4.abf, na.rm = TRUE),
             pp4_whole = whole$pp4, detail = paste0("n_signal_pairs=", nrow(s)))
}

# ------------------------------------------------------------ public API ----

#' Run one or more methods on one locus. Returns a data.frame, one row per method.
run_methods <- function(locus, methods = gated_methods, config = gated_default_config()) {
  config <- utils::modifyList(gated_default_config(), config)
  bad <- setdiff(methods, gated_methods)
  if (length(bad)) stop("unknown method: ", paste(bad, collapse = ", "))
  whole <- .whole(locus, config)
  out <- list()
  for (m in methods) {
    if (m %in% c("V1_block_any", "V2_gate_marginal", "V3_gate_conditional", "V4_cond_block"))
      out[[m]] <- .run_block_version(m, locus, config, whole)
    else if (m == "V0_handoff") out[[m]] <- .run_v0(locus, config, whole)
    else if (m == "whole_abf")
      out[[m]] <- .gated_row("whole_abf", branch = "whole", pp4 = whole$pp4, pp4_whole = whole$pp4)
    else if (m == "signals_mask") out[[m]] <- .run_signals(locus, config, whole, "mask")
    else if (m == "signals_cond") out[[m]] <- .run_signals(locus, config, whole, "cond")
  }
  if (any(c("susie_raw", "susie_fallback") %in% methods)) {
    su <- .run_susie(locus, config, whole)
    out$susie <- su[su$method %in% methods, , drop = FALSE]
  }
  res <- do.call(rbind, out)
  rownames(res) <- NULL
  res[match(intersect(gated_methods, methods), res$method), , drop = FALSE]
}
