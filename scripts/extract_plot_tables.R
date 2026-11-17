#!/usr/bin/env Rscript
# Read-only extraction of exact summary tables (with Wilson 95% CIs, window-level
# t_5 95% CIs, ROC step coordinates, and PP4 distribution quantiles) from saved
# per-locus score files for publication-grade main figures (Figures 1-5 & Supp Fig 1).
#
# No simulation and no colocalisation method is run here; existing files are never modified.
#
# Usage (from clike_ldblock_coloc_handoff_20260811 root):
#   Rscript scripts/bulk_h2ref_v1/extract_plot_tables.R \
#     --root outputs \
#     --out  outputs/bulk_h2ref_v1_plot_tables

.args <- commandArgs(trailingOnly = TRUE)
.opt <- function(name, default) {
  i <- match(paste0("--", name), .args)
  if (is.na(i) || i == length(.args)) default else .args[i + 1L]
}

out_root <- .opt("root", "outputs")
dest_dir <- .opt("out", file.path(out_root, "bulk_h2ref_v1_plot_tables"))
dir.create(dest_dir, recursive = TRUE, showWarnings = FALSE)

BLOCK0 <- c("V0_handoff", "V1_block_any", "V2_gate_marginal", "V3_gate_conditional", "V4_cond_block")
SCEN_ALL  <- c("A_A", "A_AB", "AB_AB", "AB_AC", "A_B", "AB_CD", "gwas_only", "null_both")
SCEN_POS  <- c("A_A", "A_AB", "AB_AB", "AB_AC")
SCEN_NEG  <- c("A_B", "AB_CD")
SCEN_CORE <- c(SCEN_POS, SCEN_NEG)
SCEN_LABEL <- c(
  A_A = "single_shared", A_AB = "single_double_shared",
  AB_AB = "double_shared2", AB_AC = "double_shared1",
  A_B = "single_distinct", AB_CD = "double_distinct",
  gwas_only = "gwas_only", null_both = "null_both"
)

# Clean method set (6 external/core methods + 2 internal ablations; excludes deprecated V0_handoff)
KEEP_KEYS <- c(
  "V3_gate_conditional_fb@1e-03",
  "whole_abf",
  "susie_raw",
  "susie_fallback",
  "signals_mask",
  "signals_cond",
  "V3_gate_conditional@1e-03",
  "V2_gate_marginal@1e-05"
)

fe <- function(x) formatC(x, format = "e", digits = 0)
wilson <- function(k, n) {
  if (n == 0L) return(c(NA_real_, NA_real_))
  p <- k / n; z <- 1.95996398454005
  den <- 1 + z^2 / n
  ctr <- (p + z^2 / (2 * n)) / den
  hlf <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / den
  c(ctr - hlf, ctr + hlf)
}
auc_calc <- function(s, y) {
  r <- rank(s); n1 <- sum(y == 1L); n0 <- sum(y == 0L)
  (sum(r[y == 1L]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}
win_t_ci <- function(v) {
  m <- mean(v); s <- stats::sd(v); n <- length(v)
  hw <- stats::qt(0.975, df = n - 1L) * s / sqrt(n)
  c(mean = m, lo = m - hw, hi = m + hw)
}

rd_gz <- function(path) {
  if (!file.exists(path)) stop("Missing input file: ", path)
  utils::read.delim(gzfile(path), stringsAsFactors = FALSE)
}

load_scores <- function(dir_path, stage_name) {
  d <- rd_gz(file.path(dir_path, paste0("scores_", stage_name, ".tsv.gz")))
  fb <- d[d$method %in% BLOCK0, ]
  fb$method <- paste0(fb$method, "_fb")
  fb$pp4 <- ifelse(is.finite(fb$pp4), fb$pp4, fb$pp4_whole)
  d <- rbind(d, fb)
  is_blk <- d$method %in% c(BLOCK0, paste0(BLOCK0, "_fb"))
  d$key <- ifelse(
    is_blk,
    ifelse(d$lead_threshold1 == d$lead_threshold2,
           paste0(d$method, "@", fe(d$lead_threshold1)),
           paste0(d$method, "@", fe(d$lead_threshold1), "/", fe(d$lead_threshold2))),
    d$method
  )
  d$score <- ifelse(is.finite(d$pp4), d$pp4, -Inf)
  d$no_call <- !is.finite(d$pp4)
  d
}

main_dir <- file.path(out_root, "bulk_h2ref_v1")
th <- utils::read.delim(file.path(main_dir, "thresholds.tsv"), stringsAsFactors = FALSE)
dm <- load_scores(main_dir, "main")
dg <- load_scores(main_dir, "grid")

# ----------------------------------------------------------------------
# 1. Stage 1 scenario rates + Wilson 95% CIs
# ----------------------------------------------------------------------
tm <- dm[dm$set == "test" & dm$key %in% KEEP_KEYS, ]
tm$cutoff <- th$cutoff[match(tm$key, th$key)]
tm$call <- tm$score > tm$cutoff

s1_rows <- list()
for (sc in SCEN_ALL) {
  for (ky in KEEP_KEYS) {
    sub <- tm[tm$scenario == sc & tm$key == ky, ]
    n <- nrow(sub); k <- sum(sub$call)
    ci <- wilson(k, n)
    s1_rows[[length(s1_rows) + 1L]] <- data.frame(
      scenario_id = sc, scenario = unname(SCEN_LABEL[sc]),
      truth = sub$truth[1L], key = ky,
      n = n, k_call = k, rate = k / n,
      wilson_lo = ci[1L], wilson_hi = ci[2L],
      no_call_rate = mean(sub$no_call),
      cutoff = sub$cutoff[1L],
      stringsAsFactors = FALSE
    )
  }
}
s1_tab <- do.call(rbind, s1_rows)

# Hard internal verification against saved main_table.tsv
mt_ref <- utils::read.delim(file.path(main_dir, "main_table.tsv"), stringsAsFactors = FALSE)
chk_bc <- s1_tab[s1_tab$key == "V3_gate_conditional_fb@1e-03" & s1_tab$scenario %in% mt_ref$scenario, ]
mt_bc  <- mt_ref[mt_ref$method == "block-coloc (conditional lead selection, proposed)", ]
mt_bc  <- mt_bc[match(chk_bc$scenario, mt_bc$scenario), ]
if (any(chk_bc$k_call != mt_bc$n_positive) || any(chk_bc$n != mt_bc$n_total)) {
  stop("Verification failed: stage1_scenarios_ci.tsv disagrees with main_table.tsv")
}
utils::write.table(s1_tab, file.path(dest_dir, "stage1_scenarios_ci.tsv"),
                   sep = "\t", quote = FALSE, row.names = FALSE)

# ----------------------------------------------------------------------
# 2. Stage 1 ROC coordinates & exact AUC on core test scenarios
# ----------------------------------------------------------------------
tm_core <- tm[tm$scenario %in% SCEN_CORE, ]
roc_rows <- list()
for (ky in KEEP_KEYS) {
  sub <- tm_core[tm_core$key == ky, ]
  s1 <- sub$score[sub$truth == 1L]
  s0 <- sub$score[sub$truth == 0L]
  a  <- auc_calc(sub$score, sub$truth)
  t_grid <- sort(unique(c(Inf, sub$score, -Inf)), decreasing = TRUE)
  fpr_v <- vapply(t_grid, function(t) mean(s0 >= t), 0.0)
  tpr_v <- vapply(t_grid, function(t) mean(s1 >= t), 0.0)
  keep_idx <- !duplicated(paste(fpr_v, tpr_v, sep = "_"))
  roc_rows[[length(roc_rows) + 1L]] <- data.frame(
    key = ky, auc = a, fpr = fpr_v[keep_idx], tpr = tpr_v[keep_idx],
    stringsAsFactors = FALSE
  )
}
roc_tab <- do.call(rbind, roc_rows)
old_plain_roc <- file.path(dest_dir, "stage1_roc_points.tsv")
if (file.exists(old_plain_roc)) unlink(old_plain_roc)
roc_gz <- gzfile(file.path(dest_dir, "stage1_roc_points.tsv.gz"), "wt")
utils::write.table(roc_tab, roc_gz, sep = "\t", quote = FALSE, row.names = FALSE)
close(roc_gz)

# ----------------------------------------------------------------------
# 3. Stage 2 one-factor scans + Wilson 95% CIs (verified against grid_table.tsv)
# ----------------------------------------------------------------------
dm_rep50 <- dm[dm$rep <= max(dg$rep), names(dg)]
d_all_grid <- rbind(dg, dm_rep50)
d_all_grid <- d_all_grid[d_all_grid$set == "test" & d_all_grid$key %in% KEEP_KEYS, ]
d_all_grid$cutoff <- th$cutoff[match(d_all_grid$key, th$key)]
d_all_grid$call <- d_all_grid$score > d_all_grid$cutoff

build_scan_rows <- function(df_sub, scan_name, vary_col) {
  vals <- if (vary_col == "tier") c("low", "mid", "high") else as.character(sort(unique(df_sub[[vary_col]])))
  out <- list()
  for (v in vals) {
    for (ky in KEEP_KEYS) {
      sub <- df_sub[as.character(df_sub[[vary_col]]) == v & df_sub$key == ky, ]
      pos <- sub[sub$truth == 1L, ]; neg <- sub[sub$truth == 0L, ]
      kp <- sum(pos$call); np <- nrow(pos); cip <- wilson(kp, np)
      kn <- sum(neg$call); nn <- nrow(neg); cin <- wilson(kn, nn)
      out[[length(out) + 1L]] <- data.frame(
        scan = scan_name, level = v, key = ky,
        n_pos = np, k_pos = kp, tpr = kp / np, tpr_lo = cip[1L], tpr_hi = cip[2L],
        n_neg = nn, k_neg = kn, fpr = kn / nn, fpr_lo = cin[1L], fpr_hi = cin[2L],
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, out)
}

d_core_g <- d_all_grid[d_all_grid$scenario %in% SCEN_CORE, ]
s2_gwas <- build_scan_rows(d_core_g[d_core_g$h2_eqtl == 0.015 & d_core_g$tier == "low", ], "h2_gwas", "h2_gwas")
s2_eqtl <- build_scan_rows(d_core_g[d_core_g$h2_gwas == 0.003 & d_core_g$tier == "low", ], "h2_eqtl", "h2_eqtl")
d_tier  <- d_all_grid[d_all_grid$h2_gwas == 0.003 & d_all_grid$h2_eqtl == 0.015 &
                        d_all_grid$scenario %in% c("A_AB", "AB_AB", "AB_AC", "A_B", "AB_CD"), ]
s2_tier <- build_scan_rows(d_tier, "tier", "tier")
s2_tab  <- rbind(s2_gwas, s2_eqtl, s2_tier)

# Verify against grid_table.tsv on h2_gwas and h2_eqtl scans
gt_ref <- utils::read.delim(file.path(main_dir, "grid_table.tsv"), stringsAsFactors = FALSE)
for (i in seq_len(nrow(s2_gwas))) {
  r <- s2_gwas[i, ]
  m <- gt_ref[gt_ref$h2_gwas == as.numeric(r$level) & gt_ref$h2_eqtl == 0.015 &
                gt_ref$tier == "low" & gt_ref$method == r$key, ]
  if (nrow(m) != 1L || abs(r$tpr - m$TPR) > 1e-12 || abs(r$fpr - m$FPR) > 1e-12) {
    stop("Verification failed: stage2_grid_ci.tsv disagrees with grid_table.tsv")
  }
}
utils::write.table(s2_tab, file.path(dest_dir, "stage2_grid_ci.tsv"),
                   sep = "\t", quote = FALSE, row.names = FALSE)

# ----------------------------------------------------------------------
# 4. Figure 1 exact branch decomposition + 6-window t_5 95% CI on net gain
# ----------------------------------------------------------------------
branch_one_setting <- function(raw_df, t_abf, t_bc, setting_id, setting_label) {
  pos <- raw_df[raw_df$set == "test" & raw_df$scenario %in% SCEN_POS, ]
  id_cols <- c("window", "scenario", "rep", "h2_gwas", "h2_eqtl", "tier")
  w <- pos[pos$method == "whole_abf", c(id_cols, "pp4")]
  names(w)[names(w) == "pp4"] <- "abf"
  b <- pos[pos$method == "V3_gate_conditional" &
             abs(pos$lead_threshold1 - 1e-3) < 1e-12 &
             abs(pos$lead_threshold2 - 1e-3) < 1e-12,
           c(id_cols, "pp4", "pp4_whole", "branch", "status")]
  m <- merge(w, b, by = id_cols)
  m$bc <- ifelse(is.finite(m$pp4), m$pp4, m$pp4_whole)
  m$g  <- ifelse(m$branch == "whole", "G0", ifelse(m$status == "no_valid_block", "G2", "G1"))
  if (max(abs(m$bc[m$g != "G1"] - m$abf[m$g != "G1"])) > 1e-12) {
    stop("Branch identity violated in ", setting_id)
  }
  c_abf <- m$abf > t_abf; c_bc <- m$bc > t_bc
  n <- nrow(m)
  g1_idx <- m$g == "G1"; g0_idx <- m$g == "G0"; g2_idx <- m$g == "G2"
  g1_pp  <- 100 * sum(c_bc[g1_idx] - c_abf[g1_idx]) / n
  g02_pp <- 100 * sum(c_bc[!g1_idx] - c_abf[!g1_idx]) / n
  net_pp <- 100 * mean(c_bc - c_abf)
  per_win <- tapply(100 * (c_bc - c_abf), m$window, mean)
  wci <- win_t_ci(as.numeric(per_win))
  data.frame(
    setting_id = setting_id, setting_label = setting_label, n_pos = n,
    share_G0 = 100 * mean(g0_idx), share_G1 = 100 * mean(g1_idx), share_G2 = 100 * mean(g2_idx),
    g1_pp = g1_pp, g02_pp = g02_pp, net_pp = net_pp,
    win_mean_pp = unname(wci["mean"]), win_lo_pp = unname(wci["lo"]), win_hi_pp = unname(wci["hi"]),
    wins_pos = sum(per_win > 0),
    stringsAsFactors = FALSE
  )
}

t_abf_main <- th$cutoff[th$key == "whole_abf"]
t_bc_main  <- th$cutoff[th$key == "V3_gate_conditional_fb@1e-03"]
dm_raw <- rd_gz(file.path(main_dir, "scores_main.tsv.gz"))
dg_raw <- rd_gz(file.path(main_dir, "scores_grid.tsv.gz"))
win_set <- unique(dm_raw[, c("window", "set")])
if (!"set" %in% names(dg_raw)) dg_raw$set <- win_set$set[match(dg_raw$window, win_set$window)]

b_list <- list(
  branch_one_setting(dm_raw[dm_raw$h2_gwas == 0.003 & dm_raw$h2_eqtl == 0.015 & dm_raw$tier == "low", ],
                     t_abf_main, t_bc_main, "default", "Default (0.3%, 1.5%)"),
  branch_one_setting(dg_raw[dg_raw$h2_gwas == 0.001 & dg_raw$h2_eqtl == 0.015 & dg_raw$tier == "low", ],
                     t_abf_main, t_bc_main, "gwas_0001", "GWAS h2 = 0.1%"),
  branch_one_setting(dg_raw[dg_raw$h2_gwas == 0.003 & dg_raw$h2_eqtl == 0.0075 & dg_raw$tier == "low", ],
                     t_abf_main, t_bc_main, "eqtl_0075", "eQTL h2 = 0.75%"),
  branch_one_setting(dg_raw[dg_raw$h2_gwas == 0.003 & dg_raw$h2_eqtl == 0.03 & dg_raw$tier == "low", ],
                     t_abf_main, t_bc_main, "eqtl_03", "eQTL h2 = 3.0%"),
  branch_one_setting(dg_raw[dg_raw$h2_gwas == 0.003 & dg_raw$h2_eqtl == 0.05 & dg_raw$tier == "low", ],
                     t_abf_main, t_bc_main, "eqtl_05", "eQTL h2 = 5.0%"),
  branch_one_setting(dg_raw[dg_raw$h2_gwas == 0.01 & dg_raw$h2_eqtl == 0.015 & dg_raw$tier == "low", ],
                     t_abf_main, t_bc_main, "gwas_01", "GWAS h2 = 1.0%"),
  branch_one_setting(dg_raw[dg_raw$h2_gwas == 0.02 & dg_raw$h2_eqtl == 0.015 & dg_raw$tier == "low", ],
                     t_abf_main, t_bc_main, "gwas_02", "GWAS h2 = 2.0%")
)

for (norm_spec in list(
  list(dir = "bulk_h2ref_v1_normal", id = "normal_default", lbl = "Normal (default)"),
  list(dir = "bulk_h2ref_v1_normal_g0.003_e0.05", id = "normal_e05", lbl = "Normal (eQTL 5.0%)")
)) {
  ndir <- file.path(out_root, norm_spec$dir)
  if (file.exists(file.path(ndir, "scores_main.tsv.gz"))) {
    nth <- utils::read.delim(file.path(ndir, "thresholds.tsv"), stringsAsFactors = FALSE)
    ndf <- rd_gz(file.path(ndir, "scores_main.tsv.gz"))
    b_list[[length(b_list) + 1L]] <- branch_one_setting(
      ndf,
      nth$cutoff[nth$key == "whole_abf"],
      nth$cutoff[nth$key == "V3_gate_conditional_fb@1e-03"],
      norm_spec$id, norm_spec$lbl
    )
  }
}
b_tab <- do.call(rbind, b_list)
utils::write.table(b_tab, file.path(dest_dir, "fig1_branch_ci.tsv"),
                   sep = "\t", quote = FALSE, row.names = FALSE)

# ----------------------------------------------------------------------
# 5. Supplementary Fig 1 PP4 quantiles + test-window score table (no genotypes)
# ----------------------------------------------------------------------
q_rows <- list()
for (ky in KEEP_KEYS) {
  for (sc in SCEN_ALL) {
    sub <- tm[tm$key == ky & tm$scenario == sc, ]
    ok  <- sub$pp4[is.finite(sub$pp4)]
    qs  <- if (length(ok) >= 5L) stats::quantile(ok, c(0.05, 0.25, 0.50, 0.75, 0.95)) else rep(NA_real_, 5L)
    q_rows[[length(q_rows) + 1L]] <- data.frame(
      key = ky, scenario_id = sc, scenario = unname(SCEN_LABEL[sc]),
      n = nrow(sub), n_finite = length(ok), no_call_rate = mean(!is.finite(sub$pp4)),
      cutoff = sub$cutoff[1L],
      q05 = qs[1L], q25 = qs[2L], q50 = qs[3L], q75 = qs[4L], q95 = qs[5L],
      stringsAsFactors = FALSE
    )
  }
}
q_tab <- do.call(rbind, q_rows)
utils::write.table(q_tab, file.path(dest_dir, "supp_pp4_quantiles.tsv"),
                   sep = "\t", quote = FALSE, row.names = FALSE)

pts <- tm[, c("key", "scenario", "pp4")]
pts$scenario <- unname(SCEN_LABEL[pts$scenario])
gz_con <- gzfile(file.path(dest_dir, "supp_pp4_points.tsv.gz"), "wt")
utils::write.table(pts, gz_con, sep = "\t", quote = FALSE, row.names = FALSE)
close(gz_con)

cat("SUCCESS: wrote verified plot tables to", dest_dir, "\n")
