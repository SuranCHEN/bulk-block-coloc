# LD-block colocalization for loci that may contain more than one causal signal.
#
# The method intentionally contains no cosine-similarity, optimal-transport, or
# cross-population code. It uses one aligned LD matrix, finds approximately
# independent leads in each trait, builds overlapping LD neighborhoods, and
# runs coloc.abf separately in each resulting component.

block_coloc_default_config <- function() {
  list(
    trait1_p_threshold = 5e-8,#trait进入到候选集合的pval阈值，trait1用的是GWAS，使用严格的5e-8
    trait2_p_threshold = 1e-5,#trait2通常是eQTL/QTL样本量或者功效较低，选值比GWAS宽松，pval只负责SNP进入到候选的资格，不决定PP4
    lead_min_k = 2L,#如果已经有显著 SNP，但候选数量不足 2 个，代码会补充 p-value 最小的 SNP，不保证最后一定有2个lead，后续LD clumping可能合并
    fallback_top_k = 5L,#如果没有任何 SNP 通过 p-value 阈值，就取 p-value 最小的前 5 个，fallback 不代表“这些 SNP 已经被认定为真实信号”，只是避免弱信号性状被硬阈值完全排除。
    lead_clump_r2 = 0.10,#独立lead参数
    block_r2 = 0.40,#每个lead周围的LD neighborhood定义为r² >= block_r2
    sign_split = TRUE,#lead的z- score方向可能不同，按照实际相交的lead pair拆开
    hybrid_whole_locus = TRUE,
    min_component_snps = 2L,#一个 component 至少需要多少个有效 SNP 才运行 coloc.abf
    p1 = 1e-4,#Trait 1 每个 SNP 成为 causal 的先验概率。
    p2 = 1e-4,#Trait 2 每个 SNP 成为 causal 的先验概率。
    p12 = 1e-5#同一个 SNP 同时影响两个性状的先验概率。这个值大pp4大
  )
}

.block_coloc_scalar <- function(x, fallback = NA_real_) {
  value <- suppressWarnings(as.numeric(x))
  if (length(value) == 1L && is.finite(value)) value else fallback
}

.block_coloc_validate_config <- function(config) {
  required <- names(block_coloc_default_config())
  missing <- setdiff(required, names(config))
  if (length(missing)) {
    stop("Configuration is missing: ", paste(missing, collapse = ", "))
  }
  probability_names <- c(
    "trait1_p_threshold", "trait2_p_threshold", "lead_clump_r2",
    "block_r2", "p1", "p2", "p12"
  )
  for (name in probability_names) {
    value <- .block_coloc_scalar(config[[name]])
    if (!is.finite(value) || value <= 0 || value > 1) {
      stop(name, " must be in (0, 1]")
    }
  }
  config$lead_min_k <- as.integer(config$lead_min_k)
  config$fallback_top_k <- as.integer(config$fallback_top_k)
  config$min_component_snps <- as.integer(config$min_component_snps)
  if (config$lead_min_k < 1L || config$fallback_top_k < 1L ||
      config$min_component_snps < 1L) {
    stop("lead_min_k, fallback_top_k, and min_component_snps must be positive")
  }
  config
}

.block_coloc_validate_inputs <- function(sumstats, LD, n1, n2, type1, type2, s1, s2) {
  needed <- c(
    "snp", "beta_trait1", "se_trait1", "pval_trait1",
    "beta_trait2", "se_trait2", "pval_trait2"
  )
  missing <- setdiff(needed, names(sumstats))
  if (length(missing)) stop("sumstats is missing: ", paste(missing, collapse = ", "))
  if (!"maf" %in% names(sumstats) &&
      !all(c("maf_trait1", "maf_trait2") %in% names(sumstats))) {
    stop("sumstats needs maf, or both maf_trait1 and maf_trait2")
  }
  if (!"maf_trait1" %in% names(sumstats)) sumstats$maf_trait1 <- sumstats$maf
  if (!"maf_trait2" %in% names(sumstats)) sumstats$maf_trait2 <- sumstats$maf
  if (anyDuplicated(sumstats$snp)) stop("sumstats$snp must be unique")
  if (!is.matrix(LD)) LD <- as.matrix(LD)
  storage.mode(LD) <- "double"
  if (nrow(LD) != ncol(LD)) stop("LD must be square")

  snps <- as.character(sumstats$snp)
  if (!is.null(rownames(LD)) && !is.null(colnames(LD))) {
    if (!all(snps %in% rownames(LD)) || !all(snps %in% colnames(LD))) {
      stop("Every sumstats SNP must occur in both LD row and column names")
    }
    LD <- LD[snps, snps, drop = FALSE]
  } else if (nrow(LD) != nrow(sumstats)) {
    stop("Unnamed LD must have the same number and order of SNPs as sumstats")
  } else {
    rownames(LD) <- colnames(LD) <- snps
  }
  if (any(!is.finite(LD))) stop("LD contains non-finite values")
  if (max(abs(LD - t(LD))) > 1e-6) {
    warning("LD was not exactly symmetric; replacing it with (LD + t(LD))/2")
  }
  LD <- (LD + t(LD)) / 2
  LD[LD > 1] <- 1
  LD[LD < -1] <- -1
  diag(LD) <- 1

  n1 <- .block_coloc_scalar(n1)
  n2 <- .block_coloc_scalar(n2)
  if (!is.finite(n1) || n1 <= 0 || !is.finite(n2) || n2 <= 0) {
    stop("n1 and n2 must be positive")
  }
  type1 <- match.arg(type1, c("quant", "cc"))
  type2 <- match.arg(type2, c("quant", "cc"))
  if (type1 == "cc" && (!is.finite(.block_coloc_scalar(s1)) || s1 <= 0 || s1 >= 1)) {
    stop("s1 must be in (0, 1) for a case-control trait")
  }
  if (type2 == "cc" && (!is.finite(.block_coloc_scalar(s2)) || s2 <= 0 || s2 >= 1)) {
    stop("s2 must be in (0, 1) for a case-control trait")
  }

  list(
    sumstats = as.data.frame(sumstats, stringsAsFactors = FALSE),
    LD = LD, n1 = n1, n2 = n2, type1 = type1, type2 = type2,
    s1 = s1, s2 = s2
  )
}

.block_coloc_greedy_clump <- function(indices, pvalues, LD, r2_cut) {
  indices <- unique(as.integer(indices))
  if (!length(indices)) return(integer())
  ordered <- indices[order(pvalues[indices], indices, na.last = TRUE)]
  selected <- integer()
  for (index in ordered) {
    if (!length(selected) || all(LD[index, selected]^2 < r2_cut)) {
      selected <- c(selected, index)
    }
  }
  selected
}#候选 SNP 按 p-value 从小到大处理：最小 p-value 的 SNP 先成为 lead；后续 SNP 如果和已有 lead 的 r² >= 0.10，就作为 proxy 丢弃；只有与所有已有 lead 都满足 r² < 0.10 的 SNP，才成为新的 lead。

.block_coloc_select_trait <- function(pvalues, betas, ses, LD, threshold, config) {
  pvalues <- suppressWarnings(as.numeric(pvalues))
  finite <- which(
    is.finite(pvalues) & is.finite(betas) & is.finite(ses) & ses > 0
  )
  if (!length(finite)) {
    return(list(candidates = integer(), leads = integer(), used_fallback = FALSE))
  }

  candidates <- finite[pvalues[finite] < threshold]#严格阈值会减少候选SNP，可能漏掉真实信号，更可能触发fallback；宽松阈值会增加候选SNP，可能产生更多lead，更容易把噪声引入到block
  used_fallback <- !length(candidates)
  ordered <- finite[order(pvalues[finite], finite)]
  if (used_fallback) {
    candidates <- head(ordered, config$fallback_top_k)
  } else if (length(candidates) < config$lead_min_k) {
    candidates <- union(candidates, head(ordered, config$lead_min_k))
  }
  candidates <- candidates[order(pvalues[candidates], candidates)]
  leads <- .block_coloc_greedy_clump(
    indices = candidates,
    pvalues = pvalues,
    LD = LD,
    r2_cut = config$lead_clump_r2
  )
  list(
    candidates = as.integer(candidates),
    leads = as.integer(leads),
    used_fallback = used_fallback
  )
}

.block_coloc_graph_components <- function(overlap) {
  if (!length(overlap) || !any(overlap)) return(list())
  n1 <- nrow(overlap)
  n2 <- ncol(overlap)
  seen1 <- rep(FALSE, n1)
  seen2 <- rep(FALSE, n2)
  output <- list()

  for (start in seq_len(n1)) {
    if (seen1[start] || !any(overlap[start, ])) next
    queue1 <- start
    queue2 <- integer()
    component1 <- integer()
    component2 <- integer()
    seen1[start] <- TRUE
    repeat {
      while (length(queue1)) {
        index <- queue1[[1L]]
        queue1 <- queue1[-1L]
        component1 <- union(component1, index)
        neighbors <- which(overlap[index, ])
        new <- neighbors[!seen2[neighbors]]
        if (length(new)) {
          seen2[new] <- TRUE
          queue2 <- union(queue2, new)
        }
      }
      while (length(queue2)) {
        index <- queue2[[1L]]
        queue2 <- queue2[-1L]
        component2 <- union(component2, index)
        neighbors <- which(overlap[, index])
        new <- neighbors[!seen1[neighbors]]
        if (length(new)) {
          seen1[new] <- TRUE
          queue1 <- union(queue1, new)
        }
      }
      if (!length(queue1) && !length(queue2)) break
    }
    output[[length(output) + 1L]] <- list(
      trait1 = sort(component1), trait2 = sort(component2)
    )
  }
  output
}

.block_coloc_build_components <- function(selection1, selection2, sumstats, LD, config) {
  leads1 <- selection1$leads
  leads2 <- selection2$leads
  if (!length(leads1) || !length(leads2)) return(list())

  windows1 <- lapply(leads1, function(index) which(LD[, index]^2 >= config$block_r2))
  windows2 <- lapply(leads2, function(index) which(LD[, index]^2 >= config$block_r2))
  overlap <- matrix(FALSE, nrow = length(leads1), ncol = length(leads2))
  for (i in seq_along(leads1)) {
    for (j in seq_along(leads2)) {
      overlap[i, j] <- length(intersect(windows1[[i]], windows2[[j]])) > 0L
    }
  }

  graph_components <- .block_coloc_graph_components(overlap)
  z1 <- sumstats$beta_trait1 / sumstats$se_trait1
  z2 <- sumstats$beta_trait2 / sumstats$se_trait2
  signs1 <- sign(z1[leads1])
  signs2 <- sign(z2[leads2])
  blocks <- list()

  add_block <- function(rows, cols, graph_id, split_reason) {
    indices <- sort(unique(c(
      unlist(windows1[rows], use.names = FALSE),
      unlist(windows2[cols], use.names = FALSE)
    )))
    if (length(indices) >= config$min_component_snps) {
      blocks[[length(blocks) + 1L]] <<- list(
        indices = indices,
        trait1_lead_positions = rows,
        trait2_lead_positions = cols,
        graph_component = graph_id,
        split_reason = split_reason
      )
    }
  }

  for (graph_id in seq_along(graph_components)) {
    component <- graph_components[[graph_id]]
    mixed <- length(unique(signs1[component$trait1])) > 1L ||
      length(unique(signs2[component$trait2])) > 1L
    if (isTRUE(config$sign_split) && mixed) {
      local <- which(overlap[component$trait1, component$trait2, drop = FALSE], arr.ind = TRUE)
      if (nrow(local)) {
        for (edge in seq_len(nrow(local))) {
          add_block(
            component$trait1[local[edge, "row"]],
            component$trait2[local[edge, "col"]],
            graph_id,
            "mixed_lead_sign_pair_split"
          )
        }
      }
    } else {
      add_block(component$trait1, component$trait2, graph_id, "connected_component")
    }
  }

  if (length(blocks)) {
    keys <- vapply(blocks, function(block) paste(block$indices, collapse = ","), character(1L))
    blocks <- blocks[!duplicated(keys)]
  }
  attr(blocks, "overlap") <- overlap
  attr(blocks, "windows1") <- windows1
  attr(blocks, "windows2") <- windows2
  blocks
}

.block_coloc_dataset <- function(sumstats, indices, trait, n, type, case_fraction) {
  beta_name <- paste0("beta_trait", trait)
  se_name <- paste0("se_trait", trait)
  dataset <- list(
    snp = as.character(sumstats$snp[indices]),
    beta = as.numeric(sumstats[[beta_name]][indices]),
    varbeta = as.numeric(sumstats[[se_name]][indices])^2,
    MAF = as.numeric(sumstats[[paste0("maf_trait", trait)]][indices]),
    N = as.numeric(n),
    type = type
  )
  if (type == "cc") dataset$s <- as.numeric(case_fraction)
  dataset
}

.block_coloc_run_abf <- function(inputs, indices, config) {
  ss <- inputs$sumstats
  valid <- indices[
    is.finite(ss$beta_trait1[indices]) & is.finite(ss$se_trait1[indices]) &
      ss$se_trait1[indices] > 0 & is.finite(ss$beta_trait2[indices]) &
      is.finite(ss$se_trait2[indices]) & ss$se_trait2[indices] > 0 &
      is.finite(ss$maf_trait1[indices]) & ss$maf_trait1[indices] > 0 &
      ss$maf_trait1[indices] < 1 & is.finite(ss$maf_trait2[indices]) &
      ss$maf_trait2[indices] > 0 & ss$maf_trait2[indices] < 1
  ]
  if (length(valid) < config$min_component_snps) {
    stop("Too few finite SNPs remain for coloc.abf")
  }

  fit <- NULL
  invisible(capture.output(
    fit <- suppressWarnings(suppressMessages(coloc::coloc.abf(
      dataset1 = .block_coloc_dataset(
        ss, valid, 1L, inputs$n1, inputs$type1, inputs$s1
      ),
      dataset2 = .block_coloc_dataset(
        ss, valid, 2L, inputs$n2, inputs$type2, inputs$s2
      ),
      p1 = config$p1,
      p2 = config$p2,
      p12 = config$p12
    ))),
    file = NULL
  ))
  summary <- fit$summary
  snp_pp4 <- suppressWarnings(as.numeric(fit$results$SNP.PP.H4))
  best_index <- if (length(snp_pp4) && any(is.finite(snp_pp4))) which.max(snp_pp4) else NA_integer_
  list(
    pp0 = .block_coloc_scalar(summary[["PP.H0.abf"]]),
    pp1 = .block_coloc_scalar(summary[["PP.H1.abf"]]),
    pp2 = .block_coloc_scalar(summary[["PP.H2.abf"]]),
    pp3 = .block_coloc_scalar(summary[["PP.H3.abf"]]),
    pp4 = .block_coloc_scalar(summary[["PP.H4.abf"]]),
    n_snps = length(valid),
    best_snp = if (is.finite(best_index)) as.character(fit$results$snp[[best_index]]) else NA_character_,
    best_snp_pp4 = if (is.finite(best_index)) snp_pp4[[best_index]] else NA_real_,
    indices = valid
  )
}

#' Run LD-block + coloc.abf on one aligned locus.
#'
#' @param sumstats Data frame with one row per SNP and the standard columns
#'   documented in README.md.
#' @param LD Signed LD-correlation matrix aligned to sumstats$snp.
#' @param n1,n2 Trait sample sizes.
#' @param type1,type2 "quant" or "cc".
#' @param s1,s2 Case fractions for case-control traits.
#' @param config Configuration returned by block_coloc_default_config().
#' @return A list with locus summary, per-component coloc results, block SNP
#'   membership, and selected lead SNPs.
block_coloc <- function(sumstats,
                        LD,
                        n1,
                        n2,
                        type1 = "quant",
                        type2 = "quant",
                        s1 = NA_real_,
                        s2 = NA_real_,
                        config = block_coloc_default_config()) {
  if (!requireNamespace("coloc", quietly = TRUE)) {
    stop("R package 'coloc' is required")
  }
  config <- .block_coloc_validate_config(utils::modifyList(
    block_coloc_default_config(), config
  ))
  inputs <- .block_coloc_validate_inputs(
    sumstats, LD, n1, n2, type1, type2, s1, s2
  )
  ss <- inputs$sumstats

  selection1 <- .block_coloc_select_trait(
    ss$pval_trait1, ss$beta_trait1, ss$se_trait1, inputs$LD,
    config$trait1_p_threshold, config
  )
  selection2 <- .block_coloc_select_trait(
    ss$pval_trait2, ss$beta_trait2, ss$se_trait2, inputs$LD,
    config$trait2_p_threshold, config
  )
  blocks <- .block_coloc_build_components(
    selection1, selection2, ss, inputs$LD, config
  )

  component_rows <- list()
  membership_rows <- list()
  fits <- list()
  if (length(blocks)) {
    for (component_id in seq_along(blocks)) {
      block <- blocks[[component_id]]
      fit <- tryCatch(
        .block_coloc_run_abf(inputs, block$indices, config),
        error = function(error) error
      )
      if (inherits(fit, "error")) {
        component_rows[[component_id]] <- data.frame(
          component_id = component_id,
          graph_component = block$graph_component,
          split_reason = block$split_reason,
          trait1_leads = paste(ss$snp[selection1$leads[block$trait1_lead_positions]], collapse = ";"),
          trait2_leads = paste(ss$snp[selection2$leads[block$trait2_lead_positions]], collapse = ";"),
          n_snps = length(block$indices),
          pp0 = NA_real_, pp1 = NA_real_, pp2 = NA_real_, pp3 = NA_real_, pp4 = NA_real_,
          best_snp = NA_character_, best_snp_pp4 = NA_real_,
          status = "coloc_error", error_message = conditionMessage(fit),
          stringsAsFactors = FALSE
        )
      } else {
        fits[[component_id]] <- fit
        component_rows[[component_id]] <- data.frame(
          component_id = component_id,
          graph_component = block$graph_component,
          split_reason = block$split_reason,
          trait1_leads = paste(ss$snp[selection1$leads[block$trait1_lead_positions]], collapse = ";"),
          trait2_leads = paste(ss$snp[selection2$leads[block$trait2_lead_positions]], collapse = ";"),
          n_snps = fit$n_snps,
          pp0 = fit$pp0, pp1 = fit$pp1, pp2 = fit$pp2, pp3 = fit$pp3, pp4 = fit$pp4,
          best_snp = fit$best_snp, best_snp_pp4 = fit$best_snp_pp4,
          status = "ok", error_message = NA_character_,
          stringsAsFactors = FALSE
        )
      }
      membership_rows[[component_id]] <- data.frame(
        component_id = component_id,
        snp = ss$snp[block$indices],
        stringsAsFactors = FALSE
      )
    }
  }
  components <- if (length(component_rows)) {
    do.call(rbind, component_rows)
  } else {
    data.frame(
      component_id = integer(), graph_component = integer(), split_reason = character(),
      trait1_leads = character(), trait2_leads = character(), n_snps = integer(),
      pp0 = numeric(), pp1 = numeric(), pp2 = numeric(), pp3 = numeric(), pp4 = numeric(),
      best_snp = character(), best_snp_pp4 = numeric(), status = character(),
      error_message = character(), stringsAsFactors = FALSE
    )
  }
  membership <- if (length(membership_rows)) {
    do.call(rbind, membership_rows)
  } else {
    data.frame(component_id = integer(), snp = character(), stringsAsFactors = FALSE)
  }

  ok_components <- which(components$status == "ok" & is.finite(components$pp4))
  best_component <- if (length(ok_components)) {
    ok_components[[which.max(components$pp4[ok_components])]]
  } else {
    NA_integer_
  }
  block_pp4 <- if (is.finite(best_component)) components$pp4[[best_component]] else 0

  whole <- tryCatch(
    .block_coloc_run_abf(inputs, seq_len(nrow(ss)), config),
    error = function(error) error
  )
  whole_pp4 <- if (inherits(whole, "error")) NA_real_ else whole$pp4
  use_whole <- isTRUE(config$hybrid_whole_locus) &&
    (length(selection1$leads) <= 1L || length(selection2$leads) <= 1L) &&
    is.finite(whole_pp4) && whole_pp4 > block_pp4

  if (use_whole) {
    selected <- whole
    selected_source <- "whole_locus_single_lead_hybrid"
    selected_component <- NA_integer_
  } else if (is.finite(best_component)) {
    selected <- fits[[best_component]]
    selected_source <- "maximum_block_pp4"
    selected_component <- components$component_id[[best_component]]
  } else {
    selected <- list(pp0 = NA_real_, pp1 = NA_real_, pp2 = NA_real_, pp3 = 0, pp4 = 0)
    selected_source <- "no_overlapping_valid_block"
    selected_component <- NA_integer_
  }

  summary <- data.frame(
    method = "ldblock_coloc_abf_hybrid",
    selected_source = selected_source,
    selected_component_id = selected_component,
    score = selected$pp4,
    pp0 = selected$pp0,
    pp1 = selected$pp1,
    pp2 = selected$pp2,
    pp3 = selected$pp3,
    pp4 = selected$pp4,
    maximum_block_pp4 = block_pp4,
    whole_locus_pp4 = whole_pp4,
    n_snps = nrow(ss),
    n_candidates_trait1 = length(selection1$candidates),
    n_candidates_trait2 = length(selection2$candidates),
    n_leads_trait1 = length(selection1$leads),
    n_leads_trait2 = length(selection2$leads),
    n_components = nrow(components),
    trait1_used_fallback = selection1$used_fallback,
    trait2_used_fallback = selection2$used_fallback,
    stringsAsFactors = FALSE
  )
  leads <- rbind(
    data.frame(
      trait = "trait1", snp = ss$snp[selection1$leads],
      pvalue = ss$pval_trait1[selection1$leads],
      used_fallback = selection1$used_fallback,
      stringsAsFactors = FALSE
    ),
    data.frame(
      trait = "trait2", snp = ss$snp[selection2$leads],
      pvalue = ss$pval_trait2[selection2$leads],
      used_fallback = selection2$used_fallback,
      stringsAsFactors = FALSE
    )
  )

  list(
    summary = summary,
    components = components,
    membership = membership,
    leads = leads,
    config = config
  )
}
