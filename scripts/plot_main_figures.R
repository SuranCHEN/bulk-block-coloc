#!/usr/bin/env Rscript
# Generate unified publication-grade figures for bulk-block-coloc.
# Reads strictly from verified aggregate tables under results/plot_tables/ and results/.

suppressPackageStartupMessages({
  library(grDevices)
  library(graphics)
  library(stats)
  library(utils)
})

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
if (length(file_arg) > 0L) {
  script_path <- normalizePath(sub("^--file=", "", file_arg[1L]))
  repo_root <- normalizePath(file.path(dirname(script_path), ".."))
} else {
  repo_root <- normalizePath(".")
}

out_dir <- file.path(repo_root, "docs")
pt_dir  <- file.path(repo_root, "results", "plot_tables")
res_dir <- file.path(repo_root, "results")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# Unified Statistical Palette & Plotmath Typography Helpers
# ---------------------------------------------------------------------------
COL_BC      <- "#1D4ED8"  # Deep Royal Blue: block-coloc
COL_BC_PALE <- "#EFF6FF"
COL_ABF     <- "#D97706"  # Warm Amber: coloc.abf
COL_COND    <- "#0F766E"  # Deep Teal: coloc.signals (cond)
COL_MASK    <- "#7C3AED"  # Muted Violet: coloc.signals (mask)
COL_SUSIE   <- "#64748B"  # Slate Gray: coloc.susie (fallback)
COL_RAW     <- "#94A3B8"  # Light Slate: coloc.susie
COL_NOFB    <- "#0284C7"  # Cyan-Blue: block-coloc (no fallback)
COL_MARG    <- "#38BDF8"  # Sky Blue: block-coloc (marginal leads)
COL_EQPR    <- "#059669"  # Emerald: region-equivalent block prior

COL_POS     <- "#059669"
COL_NEG     <- "#DC2626"
COL_TEXT    <- "#0F172A"
COL_MUTED   <- "#475569"
COL_GRID    <- "#E2E8F0"
COL_ZERO    <- "#334155"

open_png <- function(filename, width_in, height_in, res = 300L) {
  path <- file.path(out_dir, filename)
  png_type <- if (identical(Sys.info()[["sysname"]], "Darwin")) "quartz" else "cairo"
  png(filename = path, width = width_in, height = height_in, units = "in",
      res = res, bg = "white", type = png_type, family = "sans")
  path
}

wilson_ci <- function(k, n, z = 1.959963984540054) {
  p <- k / n
  z2 <- z * z
  denom <- 1 + z2 / n
  center <- (p + z2 / (2 * n)) / denom
  half <- (z / denom) * sqrt(p * (1 - p) / n + z2 / (4 * n * n))
  list(rate = p, lo = pmax(0, center - half), hi = pmin(1, center + half))
}

draw_grid_y <- function(yticks) {
  for (y in yticks) abline(h = y, col = COL_GRID, lwd = 0.7, lty = 1)
}

draw_grid_x <- function(xticks) {
  for (x in xticks) abline(v = x, col = COL_GRID, lwd = 0.7, lty = 1)
}

draw_panel_title <- function(tag, title_text, line = 0.8) {
  if (is.expression(title_text) || is.call(title_text)) {
    mtext(title_text, side = 3, adj = 0, line = line, cex = 0.94, col = COL_TEXT)
  } else {
    mtext(paste0(tag, "   ", title_text), side = 3, adj = 0, line = line,
          cex = 0.94, font = 2, col = COL_TEXT)
  }
}

draw_box <- function(x, y, w, h, title, body = NULL, fill = "#FFFFFF",
                     border = "#94A3B8", lwd = 1.3, title_col = COL_TEXT,
                     body_col = COL_MUTED, title_cex = 0.82, body_cex = 0.72) {
  rect(x - w / 2, y - h / 2, x + w / 2, y + h / 2, col = fill, border = border, lwd = lwd)
  if (is.null(body)) {
    text(x, y, labels = title, col = title_col, cex = title_cex, font = 2)
  } else {
    text(x, y + h * 0.17, labels = title, col = title_col, cex = title_cex, font = 2)
    text(x, y - h * 0.18, labels = body, col = body_col, cex = body_cex)
  }
}

draw_arrow <- function(x0, y0, x1, y1, col = COL_MUTED, lwd = 1.35, len = 0.08) {
  arrows(x0, y0, x1, y1, col = col, lwd = lwd, length = len, angle = 22)
}

# ---------------------------------------------------------------------------
# Figure 1: Workflow (A) + Law-of-Total-Probability Branch Decomposition (B)
# ---------------------------------------------------------------------------
make_fig1 <- function() {
  df <- read.delim(file.path(pt_dir, "fig1_branch_ci.tsv"), stringsAsFactors = FALSE)
  ord_ids <- c("normal_e05", "normal_default",
               "gwas_02", "gwas_01", "gwas_0001",
               "eqtl_05", "eqtl_03", "eqtl_0075",
               "default")
  df <- df[match(ord_ids, df$setting_id), ]

  row_expr <- c(
    normal_e05     = expression(paste("Normal effect, eQTL ", italic(h)[2]^2 == "5.0%")),
    normal_default = expression(paste("Normal effect, default")),
    gwas_02        = expression(paste("GWAS ", italic(h)[1]^2 == "2.0%")),
    gwas_01        = expression(paste("GWAS ", italic(h)[1]^2 == "1.0%")),
    gwas_0001      = expression(paste("GWAS ", italic(h)[1]^2 == "0.1%")),
    eqtl_05        = expression(paste("eQTL ", italic(h)[2]^2 == "5.0%")),
    eqtl_03        = expression(paste("eQTL ", italic(h)[2]^2 == "3.0%")),
    eqtl_0075      = expression(paste("eQTL ", italic(h)[2]^2 == "0.75%")),
    default        = expression(paste("Default (", italic(h)[1]^2 == "0.3%, ", italic(h)[2]^2 == "1.5%)"))
  )

  path <- open_png("fig1_overview.png", 12.2, 5.5)
  on.exit(dev.off(), add = TRUE)

  layout(matrix(c(1, 2), nrow = 1), widths = c(0.98, 1.25))

  # Panel A: Algorithmic Workflow (strict plotmath inequalities, no ASCII >=)
  par(mar = c(1.0, 0.8, 2.2, 0.8), xpd = TRUE)
  plot.new()
  plot.window(xlim = c(0, 100), ylim = c(0, 100))
  draw_panel_title("A", "block-coloc statistical procedure", line = 0.6)

  draw_box(50, 91, 92, 11.5,
           "Marginal summary stats + reference LD",
           expression(paste("Greedy clump (", italic(p) < 10^{-3}, ", ", hat(italic(r))^2 < 0.1,
                            ") + conditional leads ", italic(L)[1], ", ", italic(L)[2])),
           fill = "#F8FAFC", border = COL_MUTED)
  draw_arrow(50, 85.2, 50, 79.8)

  draw_box(50, 73.5, 80, 12.0,
           expression(paste("Multi-signal screen: ", "|", italic(L)[1], "|" >= 1, ",  |",
                            italic(L)[2], "|" >= 1, ",  max" >= 2)),
           expression(paste("Require ", {} >= 2, " conditional leads in at least one trait")),
           fill = "#F8FAFC", border = COL_TEXT, lwd = 1.4)

  # Left branch: G0
  lines(c(19, 19), c(67.5, 58.2), col = COL_MUTED, lwd = 1.35)
  draw_arrow(19, 58.2, 19, 54.0)
  text(13.0, 62.5, "No", cex = 0.82, font = 2, col = COL_MUTED)
  draw_box(19, 47.0, 36, 13.0,
           expression(paste("Branch ", italic(G)[0], ": whole region")),
           "coloc.abf (56.8% of positives)",
           fill = "#F1F5F9", border = "#64748B", title_col = COL_TEXT, title_cex = 0.80, body_cex = 0.70)

  # Right branch -> neighbourhood blocks
  lines(c(69, 69), c(67.5, 62.0), col = COL_BC, lwd = 1.35)
  draw_arrow(69, 62.0, 69, 56.8)
  text(75.5, 62.5, "Yes", cex = 0.82, font = 2, col = COL_BC)
  draw_box(69, 49.5, 56, 14.0,
           "Joint LD-neighbourhood blocks",
           expression(paste("Neighbourhoods {", italic(j) : hat(italic(r))[italic(jl)]^2 >= 0.4,
                            "}, ", italic(M)[italic(B)] >= 2)),
           fill = COL_BC_PALE, border = COL_BC, title_col = COL_BC, title_cex = 0.80, body_cex = 0.69)

  # Sub-branches: G2 fallback & G1 main
  lines(c(52, 23), c(42.5, 33.0), col = COL_MUTED, lwd = 1.35)
  draw_arrow(23, 33.0, 21, 29.5)
  text(44.5, 34.4, "No admissible block", cex = 0.71, font = 2, col = COL_MUTED)
  draw_box(20, 21.0, 39, 15.5,
           expression(paste("Branch ", italic(G)[2], ": fallback")),
           "coloc.abf (14.1% of positives)",
           fill = "#FFF7ED", border = COL_ABF, title_col = "#9A3412", title_cex = 0.80, body_cex = 0.70)

  lines(c(72, 72), c(42.5, 31.5), col = COL_BC, lwd = 1.5)
  draw_arrow(72, 31.5, 72, 29.8)
  text(83.5, 36.2, expression(paste({} >= 1, " block")), cex = 0.77, font = 2, col = COL_BC)
  rect(44, 9.7, 98, 29.3, col = COL_BC_PALE, border = COL_BC, lwd = 2.0)
  text(71, 24.4, expression(paste("Branch ", italic(G)[1], ": block coloc.abf")),
       col = COL_BC, cex = 0.83, font = 2)
  text(71, 19.2, expression(paste("Unadjusted marginal stats inside ", italic(B))),
       col = COL_TEXT, cex = 0.71)
  text(71, 14.0, expression(paste("Return ", max[italic(B)], " PP4(", italic(B), ") (29.2% of positives)")),
       col = COL_TEXT, cex = 0.70)

  text(50, 2.8,
       expression(paste("Identical scores on ", italic(G)[0] %union% italic(G)[2],
                        "; net gain arises from ", italic(G)[1], ".")),
       cex = 0.73, col = COL_MUTED, font = 3)

  # Panel B: Law-of-Total-Probability Decomposition + Window t5 95% CI (1-decimal uniform)
  par(mar = c(4.2, 12.6, 2.2, 2.4), xpd = FALSE)
  n <- nrow(df)
  y <- seq_len(n)
  xlim <- c(-2.0, 24.5)
  plot(NA, xlim = xlim, ylim = c(0.4, n + 1.45), yaxt = "n", xaxt = "n",
       xlab = expression(paste("Rate difference: block-coloc ", -{}, " coloc.abf (percentage points)")),
       ylab = "", bty = "n", cex.lab = 0.92)
  draw_panel_title("B", "Law-of-total-probability decomposition by branch", line = 0.6)

  xticks <- seq(0, 16, by = 4)
  draw_grid_x(xticks)
  segments(0, 0.4, 0, n + 0.55, col = COL_ZERO, lwd = 1.1, lty = 2)
  abline(h = c(2.5, 5.5, 8.5), col = "#CBD5E1", lwd = 0.8, lty = 3)

  bh <- 0.56
  for (i in seq_len(n)) {
    g1  <- df$g1_pp[i]
    g02 <- df$g02_pp[i]
    if (g1 != 0) rect(0, y[i] - bh / 2, g1, y[i] + bh / 2, col = COL_BC, border = NA)
    if (g02 < 0) rect(g02, y[i] - bh / 2, 0, y[i] + bh / 2, col = "#F97316", border = NA)
    segments(df$win_lo_pp[i], y[i], df$win_hi_pp[i], y[i], col = COL_TEXT, lwd = 1.6)
    segments(df$win_lo_pp[i], y[i] - 0.14, df$win_lo_pp[i], y[i] + 0.14, col = COL_TEXT, lwd = 1.3)
    segments(df$win_hi_pp[i], y[i] - 0.14, df$win_hi_pp[i], y[i] + 0.14, col = COL_TEXT, lwd = 1.3)
    points(df$win_mean_pp[i], y[i], pch = 19, col = COL_TEXT, cex = 1.05)

    lbl <- gsub("-", "\u2212", sprintf("%+.1f (%+.1f, %+.1f)", df$win_mean_pp[i], df$win_lo_pp[i], df$win_hi_pp[i]))
    text(max(df$win_hi_pp[i], df$g1_pp[i]) + 0.35, y[i], labels = lbl,
         adj = c(0, 0.5), cex = 0.76, col = COL_TEXT)
  }

  axis(1, at = c(-2, xticks), labels = c("-2", xticks), cex.axis = 0.88, col = COL_MUTED)
  par(xpd = TRUE)
  for (i in seq_len(n)) {
    text(-2.2, y[i], labels = row_expr[df$setting_id[i]], adj = c(1, 0.5),
         cex = 0.84, col = COL_TEXT)
  }
  legend(x = 10.2, y = n + 1.48, xjust = 0.5, yjust = 1,
         legend = c(expression(paste(italic(G)[1], " block")),
                    expression(paste(italic(G)[0] %union% italic(G)[2], " shift")),
                    expression(paste("Net (95% ", italic(t)[5], " CI)"))),
         fill = c(COL_BC, "#F97316", NA),
         border = c(NA, NA, NA),
         pch = c(NA, NA, 19),
         lty = c(NA, NA, 1),
         col = c(NA, NA, COL_TEXT),
         ncol = 3, bty = "n", cex = 0.75, x.intersp = 0.45)
  par(xpd = FALSE)
  message("Wrote: ", path)
}

# ---------------------------------------------------------------------------
# Figure 2 Stage 1: Scenario Bars with Wilson 95% CI (A) + ROC with AUC CIs (B)
# ---------------------------------------------------------------------------
make_fig2_stage1 <- function() {
  sc  <- read.delim(file.path(pt_dir, "stage1_scenarios_ci.tsv"), stringsAsFactors = FALSE)
  roc <- read.delim(gzfile(file.path(pt_dir, "stage1_roc_points.tsv.gz")), stringsAsFactors = FALSE)
  pt  <- read.delim(file.path(res_dir, "h2ref_stage1", "pooled_table.tsv"), stringsAsFactors = FALSE)

  scen_order <- c("single_shared", "single_double_shared", "double_shared2", "double_shared1",
                  "single_distinct", "double_distinct", "gwas_only", "null_both")
  scen_labs <- c("single\nshared", "single\ndouble\nshared", "double\nshared2", "double\nshared1",
                 "single\ndistinct", "double\ndistinct", "gwas\nonly", "null\nboth")

  keys_bar <- c("V3_gate_conditional_fb@1e-03", "whole_abf", "susie_fallback",
                "signals_mask", "signals_cond")
  labs_bar <- c("block-coloc", "coloc.abf", "coloc.susie (fallback)",
                "coloc.signals (mask)", "coloc.signals (cond)")
  cols_bar <- c(COL_BC, COL_ABF, COL_SUSIE, COL_MASK, COL_COND)

  # Map pooled_table AUC [95% CI] strings verbatim
  auc_map <- list(
    "V3_gate_conditional_fb@1e-03" = pt$AUC[pt$method == "block-coloc (conditional lead selection, proposed)"],
    "V2_gate_marginal@1e-05"       = pt$AUC[pt$method == "block-coloc (marginal lead selection)"],
    "V3_gate_conditional@1e-03"    = pt$AUC[pt$method == "block-coloc (no fallback)"],
    "signals_cond"                 = pt$AUC[pt$method == "coloc.signals (cond)"],
    "signals_mask"                 = pt$AUC[pt$method == "coloc.signals (mask)"],
    "whole_abf"                    = pt$AUC[pt$method == "coloc.abf, whole region"],
    "susie_fallback"               = pt$AUC[pt$method == "coloc.susie + coloc.abf fallback"],
    "susie_raw"                    = pt$AUC[pt$method == "coloc.susie"]
  )

  path <- open_png("fig2_stage1.png", 12.4, 5.2)
  on.exit(dev.off(), add = TRUE)

  layout(matrix(c(1, 2), nrow = 1), widths = c(1.38, 1.0))

  # Panel A: Grouped Bar Chart + Wilson 95% CIs
  par(mar = c(4.8, 4.1, 2.8, 0.8))
  n_sc <- length(scen_order)
  n_m  <- length(keys_bar)
  centers <- seq_len(n_sc)
  bw <- 0.145
  offsets <- (seq_len(n_m) - (n_m + 1) / 2) * bw

  plot(NA, xlim = c(0.45, n_sc + 0.55), ylim = c(0, 108), xaxt = "n", yaxt = "n",
       xlab = "", ylab = "Colocalisation call rate (%)", bty = "n")
  draw_panel_title("A", "Matched 5% FPR across 8 scenarios (test windows)", line = 1.2)

  rect(0.48, 0, 4.46, 102, col = "#F0FDF4", border = NA)
  rect(4.54, 0, 6.46, 102, col = "#FEF2F2", border = NA)
  rect(6.54, 0, 8.52, 102, col = "#FFFBEB", border = NA)
  draw_grid_y(seq(0, 100, by = 20))
  abline(h = 5, col = COL_NEG, lty = 2, lwd = 1.0)
  abline(v = c(4.5, 6.5), col = COL_MUTED, lty = 3, lwd = 1.0)

  for (si in seq_len(n_sc)) {
    all_zero <- TRUE
    for (mi in seq_len(n_m)) {
      row <- sc[sc$scenario == scen_order[si] & sc$key == keys_bar[mi], ]
      x0 <- centers[si] + offsets[mi]
      r  <- row$rate * 100
      lo <- row$wilson_lo * 100
      hi <- row$wilson_hi * 100
      if (r > 0) {
        all_zero <- FALSE
        rect(x0 - bw * 0.45, 0, x0 + bw * 0.45, r, col = cols_bar[mi], border = NA)
        segments(x0, lo, x0, hi, col = COL_TEXT, lwd = 1.1)
      }
    }
    if (all_zero) {
      text(centers[si], 2.0, "0.0%", cex = 0.72, col = COL_MUTED)
    }
  }

  axis(2, at = seq(0, 100, by = 20), las = 1, cex.axis = 0.88, col = COL_MUTED)
  axis(1, at = centers, labels = scen_labs, tick = FALSE, line = 0.5, cex.axis = 0.73)
  text(2.5, 104.5, "Colocalised (power)", cex = 0.80, font = 2, col = "#166534")
  text(5.5, 104.5, "2-trait distinct", cex = 0.80, font = 2, col = "#991B1B")
  text(7.5, 104.5, "Extra negatives", cex = 0.80, font = 2, col = "#92400E")
  legend("topright", legend = labs_bar, fill = cols_bar, border = NA, bty = "n",
         cex = 0.74, inset = c(0.01, 0.08))

  # Panel B: Test-Window ROC Curves with Bootstrap AUC 95% CIs
  par(mar = c(4.8, 4.4, 2.8, 1.4))
  plot(NA, xlim = c(0, 1), ylim = c(0, 1.02), xaxs = "i", yaxs = "i",
       xaxt = "n", yaxt = "n",
       xlab = "False-positive rate (2-trait distinct, n = 1,200)",
       ylab = "True-positive rate (colocalised, n = 2,400)", bty = "n")
  draw_panel_title("B", "Test-window ROC curves (AUC [95% CI])", line = 1.2)

  draw_grid_x(seq(0.2, 0.8, by = 0.2))
  draw_grid_y(seq(0.2, 1.0, by = 0.2))
  abline(0, 1, col = "#CBD5E1", lty = 3, lwd = 1.0)
  abline(v = 0.05, col = COL_NEG, lty = 2, lwd = 1.0)
  text(0.07, 0.06, "5% FPR", adj = c(0, 0), cex = 0.74, col = COL_NEG)

  roc_specs <- list(
    list(key = "V3_gate_conditional_fb@1e-03", lab = "block-coloc",                 col = COL_BC,    lty = 1, lwd = 2.6),
    list(key = "V2_gate_marginal@1e-05",       lab = "block-coloc (marginal)",      col = COL_MARG,  lty = 1, lwd = 1.6),
    list(key = "V3_gate_conditional@1e-03",    lab = "block-coloc (no fallback)",   col = COL_NOFB,  lty = 1, lwd = 1.5),
    list(key = "signals_cond",                 lab = "coloc.signals (cond)",        col = COL_COND,  lty = 2, lwd = 2.0),
    list(key = "signals_mask",                 lab = "coloc.signals (mask)",        col = COL_MASK,  lty = 2, lwd = 1.8),
    list(key = "whole_abf",                    lab = "coloc.abf",                   col = COL_ABF,   lty = 2, lwd = 2.0),
    list(key = "susie_fallback",               lab = "coloc.susie (fallback)",      col = COL_SUSIE, lty = 4, lwd = 1.7),
    list(key = "susie_raw",                    lab = "coloc.susie",                 col = COL_RAW,   lty = 3, lwd = 1.6)
  )

  for (idx in rev(seq_along(roc_specs))) {
    sp <- roc_specs[[idx]]
    sub <- roc[roc$key == sp$key, ]
    sub <- sub[order(sub$fpr, sub$tpr), ]
    lines(sub$fpr, sub$tpr, col = sp$col, lty = sp$lty, lwd = sp$lwd)
  }
  leg_txt <- vapply(roc_specs, function(sp) sprintf("%s: %s", sp$lab, auc_map[[sp$key]]), "")
  legend("bottomright", legend = leg_txt,
         col = vapply(roc_specs, `[[`, "", "col"),
         lty = vapply(roc_specs, `[[`, 1, "lty"),
         lwd = vapply(roc_specs, `[[`, 1, "lwd"),
         bty = "n", cex = 0.65, inset = c(0.01, 0.02))
  axis(1, cex.axis = 0.86, col = COL_MUTED)
  axis(2, las = 1, cex.axis = 0.86, col = COL_MUTED)
  message("Wrote: ", path)
}

# ---------------------------------------------------------------------------
# Figure 2 Stage 2: Log-Scale Parameter Scans with Wilson 95% CIs
# ---------------------------------------------------------------------------
make_fig2_stage2 <- function() {
  gd <- read.delim(file.path(pt_dir, "stage2_grid_ci.tsv"), stringsAsFactors = FALSE)

  methods <- c("V3_gate_conditional_fb@1e-03", "whole_abf", "susie_fallback",
               "signals_mask", "signals_cond")
  m_labs  <- c("block-coloc", "coloc.abf", "coloc.susie (fallback)",
               "coloc.signals (mask)", "coloc.signals (cond)")
  m_cols  <- c(COL_BC, COL_ABF, COL_SUSIE, COL_MASK, COL_COND)
  m_pchs  <- c(19, 15, 18, 17, 4)
  m_ltys  <- c(1, 2, 4, 2, 2)
  m_lwds  <- c(2.4, 1.8, 1.6, 1.8, 1.9)

  path <- open_png("fig2_stage2.png", 12.0, 6.8)
  on.exit(dev.off(), add = TRUE)

  layout(matrix(c(1, 2, 3,
                  4, 5, 6,
                  7, 7, 7), nrow = 3, byrow = TRUE),
         heights = c(1.0, 1.0, 0.22))
  par(oma = c(0, 0, 1.8, 0))

  # Both h2_1 and h2_2 on log10 scale so 0.3% default tick is cleanly labelled
  scans <- list(
    list(scan = "h2_gwas", levels = c("0.001", "0.003", "0.01", "0.02"),
         xvals = log10(c(0.1, 0.3, 1.0, 2.0)), xlabs = c("0.1%", "0.3%*", "1.0%", "2.0%"),
         xlab = expression(paste("GWAS ", italic(h)[1]^2, " (%, log scale)")),
         title = expression(bold(paste("GWAS ", italic(h)[1]^2, " scan"))), fpr_ylim = c(0, 12.2)),
    list(scan = "h2_eqtl", levels = c("0.0075", "0.015", "0.03", "0.05"),
         xvals = log10(c(0.75, 1.5, 3.0, 5.0)), xlabs = c("0.75%", "1.5%*", "3.0%", "5.0%"),
         xlab = expression(paste("eQTL ", italic(h)[2]^2, " (%, log scale)")),
         title = expression(bold(paste("eQTL ", italic(h)[2]^2, " scan"))), fpr_ylim = c(0, 12.2)),
    list(scan = "tier", levels = c("low", "mid", "high"),
         xvals = c(1, 2, 3), xlabs = c("low* (<0.2)", "mid (0.32-0.63)", "high (0.63-0.9)"),
         xlab = expression(paste("Within-trait causal |", italic(r), "| tier")),
         title = expression(bold("Causal-pair LD scan")), fpr_ylim = c(0, 46))
  )

  # Top row: Detection rate (%)
  for (ci in seq_along(scans)) {
    sp <- scans[[ci]]
    par(mar = c(2.4, if (ci == 1) 4.2 else 3.0, 2.2, 1.0))
    xr <- range(sp$xvals)
    pad <- diff(xr) * 0.08
    plot(NA, xlim = c(xr[1] - pad, xr[2] + pad), ylim = c(38, 102),
         xaxt = "n", yaxt = "n", xlab = "",
         ylab = if (ci == 1) "Detection rate (%)" else "", bty = "n")
    mtext(sp$title, side = 3, line = 0.4, cex = 0.90, col = COL_TEXT)
    draw_grid_y(seq(40, 100, by = 20))
    axis(2, at = seq(40, 100, by = 20), las = 1, cex.axis = 0.84, col = COL_MUTED)
    axis(1, at = sp$xvals, labels = sp$xlabs, cex.axis = 0.80, col = COL_MUTED)

    jit <- seq(-0.02, 0.02, length.out = length(methods)) * diff(xr)
    for (mi in rev(seq_along(methods))) {
      sub <- gd[gd$scan == sp$scan & gd$key == methods[mi], ]
      sub <- sub[match(sp$levels, sub$level), ]
      xv <- sp$xvals + jit[mi]
      lines(xv, sub$tpr * 100, col = m_cols[mi], lty = m_ltys[mi], lwd = m_lwds[mi])
      segments(xv, sub$tpr_lo * 100, xv, sub$tpr_hi * 100, col = m_cols[mi], lwd = 1.1)
      points(xv, sub$tpr * 100, col = m_cols[mi], pch = m_pchs[mi], cex = 1.05)
    }
  }

  # Bottom row: FPR (%)
  for (ci in seq_along(scans)) {
    sp <- scans[[ci]]
    par(mar = c(4.0, if (ci == 1) 4.2 else 3.0, 1.4, 1.0))
    xr <- range(sp$xvals)
    pad <- diff(xr) * 0.08
    plot(NA, xlim = c(xr[1] - pad, xr[2] + pad), ylim = sp$fpr_ylim,
         xaxt = "n", yaxt = "n", xlab = sp$xlab,
         ylab = if (ci == 1) "False-positive rate (%)" else "", bty = "n")
    yticks <- if (max(sp$fpr_ylim) > 15) seq(0, 40, by = 10) else seq(0, 10, by = 2.5)
    draw_grid_y(yticks)
    abline(h = 5, col = COL_NEG, lty = 2, lwd = 1.0)
    axis(2, at = yticks, las = 1, cex.axis = 0.84, col = COL_MUTED)
    axis(1, at = sp$xvals, labels = sp$xlabs, cex.axis = 0.80, col = COL_MUTED)

    jit <- seq(-0.02, 0.02, length.out = length(methods)) * diff(xr)
    for (mi in rev(seq_along(methods))) {
      sub <- gd[gd$scan == sp$scan & gd$key == methods[mi], ]
      sub <- sub[match(sp$levels, sub$level), ]
      xv <- sp$xvals + jit[mi]
      lines(xv, sub$fpr * 100, col = m_cols[mi], lty = m_ltys[mi], lwd = m_lwds[mi])
      segments(xv, sub$fpr_lo * 100, xv, sub$fpr_hi * 100, col = m_cols[mi], lwd = 1.1)
      points(xv, sub$fpr * 100, col = m_cols[mi], pch = m_pchs[mi], cex = 1.05)
    }
  }

  par(mar = c(0, 0, 0, 0))
  plot.new()
  plot.window(xlim = c(0, 1), ylim = c(0, 1), xaxs = "i", yaxs = "i")
  seg_w <- 0.032
  gap_in <- 0.007
  gap_between <- 0.032
  lab_w <- strwidth(m_labs, units = "user", cex = 0.85)
  total_w <- sum(seg_w + gap_in + lab_w) + (length(m_labs) - 1) * gap_between
  x_cur <- 0.5 - total_w / 2
  for (mi in seq_along(m_labs)) {
    segments(x_cur, 0.5, x_cur + seg_w, 0.5, col = m_cols[mi], lty = m_ltys[mi], lwd = m_lwds[mi])
    points(x_cur + seg_w / 2, 0.5, col = m_cols[mi], pch = m_pchs[mi], cex = 1.05)
    text(x_cur + seg_w + gap_in, 0.5, labels = m_labs[mi], adj = c(0, 0.5), cex = 0.85, col = COL_TEXT)
    x_cur <- x_cur + seg_w + gap_in + lab_w[mi] + gap_between
  }
  mtext("C   Stage 2 single-factor parameter scans (test windows, * marks Stage 1 default)",
        side = 3, outer = TRUE, adj = 0.02, line = 0.2, cex = 0.96, font = 2, col = COL_TEXT)
  message("Wrote: ", path)
}

# ---------------------------------------------------------------------------
# Figure 3: Prop S1 Stratified Power (A, B) + Transition Gap Curve (C)
# ---------------------------------------------------------------------------
make_fig3 <- function() {
  st <- read.delim(file.path(res_dir, "h2ref_stage2", "propC_strata.tsv"), stringsAsFactors = FALSE)
  gp <- read.delim(file.path(pt_dir, "propC_gap_iqr.tsv"), stringsAsFactors = FALSE)

  path <- open_png("fig3_propC.png", 12.2, 4.4)
  on.exit(dev.off(), add = TRUE)

  layout(matrix(c(1, 2, 3), nrow = 1), widths = c(1.0, 1.0, 1.08))

  h2_vals <- c(0.75, 1.5, 3.0, 5.0)
  h2_log  <- log10(h2_vals)
  h2_raw  <- c(0.0075, 0.015, 0.03, 0.05)
  scens   <- c("A_AB", "AB_AC")
  titles  <- c("Trait 1: A; Trait 2: A, B", "Trait 1: A, B; Trait 2: A, C")
  tags    <- c("A", "B")

  for (si in 1:2) {
    par(mar = c(4.3, 4.1, 2.6, 0.9))
    plot(NA, xlim = range(h2_log) + c(-0.08, 0.08), ylim = c(0, 104),
         xaxt = "n", yaxt = "n",
         xlab = expression(paste("eQTL ", italic(h)[2]^2, " (%, log scale)")),
         ylab = if (si == 1) "Detection rate (%)" else "", bty = "n")
    draw_panel_title(tags[si], titles[si], line = 0.9)
    draw_grid_y(seq(0, 100, by = 20))
    axis(2, at = seq(0, 100, by = 20), las = 1, cex.axis = 0.86, col = COL_MUTED)
    axis(1, at = h2_log, labels = c("0.75%", "1.5%", "3.0%", "5.0%"), cex.axis = 0.84, col = COL_MUTED)

    curves <- list(
      list(method = "whole", ss = "True",  col = COL_ABF, lty = 1, pch = 15, j = -0.012, lab = "coloc.abf (A strongest)"),
      list(method = "bc",    ss = "True",  col = COL_BC,  lty = 1, pch = 19, j =  0.012, lab = "block-coloc (A strongest)"),
      list(method = "whole", ss = "False", col = COL_ABF, lty = 2, pch = 0,  j = -0.012, lab = "coloc.abf (A not strongest)"),
      list(method = "bc",    ss = "False", col = COL_BC,  lty = 2, pch = 1,  j =  0.012, lab = "block-coloc (A not strongest)")
    )
    for (cv in curves) {
      sub <- st[st$scenario == scens[si] & st$method == cv$method & st$shared_strongest == cv$ss, ]
      sub <- sub[match(h2_raw, sub$h2_eqtl), ]
      xv <- h2_log + cv$j
      lines(xv, sub$rate, col = cv$col, lty = cv$lty, lwd = 2.1)
      segments(xv, sub$lo, xv, sub$hi, col = cv$col, lwd = 1.1)
      points(xv, sub$rate, col = cv$col, pch = cv$pch, cex = 1.05, bg = "white")
    }
    if (si == 1) {
      legend("bottomright", legend = vapply(curves, `[[`, "", "lab"),
             col = vapply(curves, `[[`, "", "col"),
             lty = vapply(curves, `[[`, 1, "lty"),
             pch = vapply(curves, `[[`, 1, "pch"),
             lwd = 2.0, bty = "n", cex = 0.73, inset = c(0.01, 0.02))
    }
  }

  # Panel C: Transition Gap Curve (Median + IQR)
  par(mar = c(4.3, 4.1, 2.6, 1.2))
  xr <- range(gp$gap_med) + c(-1.5, 1.5)
  plot(NA, xlim = xr, ylim = c(0, 1.03), xaxt = "n", yaxt = "n",
       xlab = expression(paste("Signal-strength gap ", italic(g) == (italic(z)[paste(2, italic(B))]^2 - italic(z)[paste(2, italic(A))]^2) / 2)),
       ylab = expression(paste("Whole-region coloc.abf ", PP4)), bty = "n")
  draw_panel_title("C", "PP4 across signal gap (n = 3,000)", line = 0.9)
  draw_grid_y(seq(0, 1, by = 0.25))
  abline(v = 0, col = COL_ZERO, lty = 2, lwd = 1.1)
  abline(h = 0.5, col = COL_MUTED, lty = 3, lwd = 1.0)

  lines(gp$gap_med, gp$q50, col = COL_ABF, lwd = 2.1)
  segments(gp$gap_med, gp$q25, gp$gap_med, gp$q75, col = COL_ABF, lwd = 2.0)
  points(gp$gap_med, gp$q50, pch = 19, col = COL_TEXT, cex = 1.05)

  text(-15.0, 0.55, expression(paste(italic(g) < 0, ": shared ", italic(A), " stronger")),
       adj = c(0, 0), cex = 0.74, col = COL_MUTED)
  text(15.0, 0.55, expression(paste(italic(g) > 0, ": private ", italic(B), " stronger")),
       adj = c(1, 0), cex = 0.74, col = COL_MUTED)
  axis(1, cex.axis = 0.84, col = COL_MUTED)
  axis(2, at = seq(0, 1, by = 0.25), las = 1, cex.axis = 0.84, col = COL_MUTED)
  legend("bottomleft", legend = c("Bin median PP4", "Interquartile range (IQR)"),
         pch = c(19, NA), lty = c(1, 1), lwd = c(2.0, 2.2),
         col = c(COL_TEXT, COL_ABF), bty = "n", cex = 0.74, inset = c(0.01, 0.02))
  message("Wrote: ", path)
}

# ---------------------------------------------------------------------------
# Figure 4: Balanced 3-Setting Dumbbell Rates (A) + Complete Forest Plot (B)
# ---------------------------------------------------------------------------
make_fig4 <- function() {
  # Complete verified 3-setting paired results from results/h2ref_insample_*/LD_COMPARE.md
  rates_df <- data.frame(
    setting = rep(c("default", "g02", "e05"), each = 5),
    method  = rep(c("block-coloc", "coloc.abf", "coloc.susie (fallback)",
                    "coloc.signals (mask)", "coloc.signals (cond)"), times = 3),
    ins = c(77.0, 73.2, 72.5, 78.5, 79.5,
            76.8, 70.1, 72.2, 78.7, 82.9,
            97.5, 89.8, 87.5, 95.8, 97.9),
    ref = c(76.8, 73.2, 72.6, 78.2, 76.2,
            76.3, 70.1, 70.9, 75.2, 76.2,
            97.7, 89.8, 87.2, 94.8, 98.2),
    col = rep(c(COL_BC, COL_ABF, COL_SUSIE, COL_MASK, COL_COND), times = 3),
    stringsAsFactors = FALSE
  )

  # Complete 3-setting window-level paired CIs for all 4 LD-dependent methods
  forest_df <- data.frame(
    setting = rep(c("default", "g02", "e05"), each = 4),
    method  = rep(c("block-coloc", "coloc.susie (fallback)",
                    "coloc.signals (mask)", "coloc.signals (cond)"), times = 3),
    mean = c(-0.3,  0.0, -0.3, -3.3,
             -0.4, -1.3, -3.5, -6.7,
              0.2, -0.3, -1.1,  0.3),
    lo   = c(-0.8, -0.3, -0.8, -4.4,
             -1.3, -2.2, -4.9, -7.8,
             -0.5, -2.0, -2.3, -0.6),
    hi   = c( 0.3,  0.3,  0.3, -2.2,
              0.4, -0.3, -2.1, -5.5,
              0.8,  1.5,  0.2,  1.1),
    col  = rep(c(COL_BC, COL_SUSIE, COL_MASK, COL_COND), times = 3),
    stringsAsFactors = FALSE
  )

  path <- open_png("fig4_ld_comparison.png", 12.2, 5.8)
  on.exit(dev.off(), add = TRUE)

  layout(matrix(c(1, 2), nrow = 1), widths = c(1.08, 1.0))

  # Panel A: High-resolution paired dumbbell plot across all 3 settings x 5 methods
  # Y coordinates: eQTL 5% (1..5), gap, GWAS 2% (7..11), gap, Default (13..17)
  ya <- c(17:13, 11:7, 5:1)
  par(mar = c(4.6, 10.8, 2.8, 1.0))
  plot(NA, xlim = c(68, 101), ylim = c(0.4, 18.4), xaxt = "n", yaxt = "n",
       xlab = "Detection rate at matched 5% FPR (%)", ylab = "", bty = "n")
  draw_panel_title("A", "Matched in-sample LD vs 1000G reference LD", line = 1.1)
  draw_grid_x(seq(70, 100, by = 5))
  abline(h = c(6, 12), col = "#CBD5E1", lwd = 0.9, lty = 2)

  for (i in seq_len(nrow(rates_df))) {
    yy <- ya[i]
    segments(rates_df$ins[i], yy, rates_df$ref[i], yy, col = rates_df$col[i], lwd = 2.0)
    points(rates_df$ins[i], yy, pch = 21, bg = "white", col = rates_df$col[i], cex = 1.25, lwd = 1.6)
    points(rates_df$ref[i], yy, pch = 19, col = rates_df$col[i], cex = 1.25)
  }

  axis(1, at = seq(70, 100, by = 5), cex.axis = 0.85, col = COL_MUTED)
  par(xpd = TRUE)
  for (i in seq_len(nrow(rates_df))) {
    text(67.2, ya[i], labels = rates_df$method[i], adj = c(1, 0.5),
         cex = 0.78, font = if (rates_df$method[i] == "block-coloc") 2 else 1, col = rates_df$col[i])
  }
  text(68.2, 17.85, "Default", adj = c(0, 0.5), cex = 0.80, font = 2, col = COL_TEXT)
  text(68.2, 11.85, expression(bold(paste("GWAS ", italic(h)[1]^2 == "2.0%"))),
       adj = c(0, 0.5), cex = 0.80, col = COL_TEXT)
  text(68.2,  5.85, expression(bold(paste("eQTL ", italic(h)[2]^2 == "5.0%"))),
       adj = c(0, 0.5), cex = 0.80, col = COL_TEXT)
  legend(x = 84, y = 18.8, legend = c("Matched in-sample LD", "1000G EUR reference LD"),
         pch = c(21, 19), pt.bg = c("white", COL_TEXT), col = c(COL_TEXT, COL_TEXT),
         bty = "n", cex = 0.76, pt.cex = 1.15)
  par(xpd = FALSE)

  # Panel B: Complete 3-setting Forest Plot (12 rows across 4 LD-dependent methods)
  yb <- c(14:11, 9.5:6.5, 5:2)
  par(mar = c(4.6, 10.6, 2.8, 2.4))
  plot(NA, xlim = c(-9.0, 9.0), ylim = c(0.4, 15.6), yaxt = "n", xaxt = "n",
       xlab = expression(paste("Paired change: Reference ", -{}, " Matched LD (percentage points)")),
       ylab = "", bty = "n", cex.lab = 0.90)
  draw_panel_title("B", expression(bold(paste("B   6-test-window mean & 95% ", italic(t)[5], " CI"))), line = 1.1)
  draw_grid_x(seq(-8, 2, by = 2))
  segments(0, 1.2, 0, 15.1, col = COL_ZERO, lty = 2, lwd = 1.1)
  abline(h = c(5.75, 10.25), col = "#CBD5E1", lwd = 0.9, lty = 2)

  for (i in seq_len(nrow(forest_df))) {
    yy <- yb[i]
    segments(forest_df$lo[i], yy, forest_df$hi[i], yy, col = forest_df$col[i], lwd = 2.0)
    segments(forest_df$lo[i], yy - 0.14, forest_df$lo[i], yy + 0.14, col = forest_df$col[i], lwd = 1.3)
    segments(forest_df$hi[i], yy - 0.14, forest_df$hi[i], yy + 0.14, col = forest_df$col[i], lwd = 1.3)
    points(forest_df$mean[i], yy, pch = 19, col = forest_df$col[i], cex = 1.18)
    m_fmt  <- if (abs(forest_df$mean[i]) < 1e-9) "0.0" else sprintf("%+.1f", forest_df$mean[i])
    ci_txt <- gsub("-", "\u2212", sprintf("%s (%+.1f, %+.1f)", m_fmt, forest_df$lo[i], forest_df$hi[i]))
    text(2.1, yy, labels = ci_txt, adj = c(0, 0.5), cex = 0.74, col = COL_TEXT)
  }

  axis(1, at = seq(-8, 2, by = 2), cex.axis = 0.85, col = COL_MUTED)
  par(xpd = TRUE)
  for (i in seq_len(nrow(forest_df))) {
    text(-9.3, yb[i], labels = forest_df$method[i], adj = c(1, 0.5),
         cex = 0.77, font = if (forest_df$method[i] == "block-coloc") 2 else 1, col = forest_df$col[i])
  }
  text(-8.8, 14.95, "Default", adj = c(0, 0.5), cex = 0.78, font = 2, col = COL_TEXT)
  text(-8.8, 10.05, expression(bold(paste("GWAS ", italic(h)[1]^2 == "2.0%"))),
       adj = c(0, 0.5), cex = 0.78, col = COL_TEXT)
  text(-8.8,  5.55, expression(bold(paste("eQTL ", italic(h)[2]^2 == "5.0%"))),
       adj = c(0, 0.5), cex = 0.78, col = COL_TEXT)
  text(0, 0.8, "coloc.abf (LD-free): 0.0 (0.0, 0.0) across all settings",
       cex = 0.72, font = 3, col = COL_MUTED)
  par(xpd = FALSE)
  message("Wrote: ", path)
}

# ---------------------------------------------------------------------------
# Figure 5: Fixed Three-Way Calling Rule (A) + Data vs Block-Prior Leverage (B)
# ---------------------------------------------------------------------------
make_fig5 <- function() {
  tc_def <- read.delim(file.path(res_dir, "h2ref_prior_default", "prior_rule_traceCB.tsv"), stringsAsFactors = FALSE)
  tc_e05 <- read.delim(file.path(res_dir, "h2ref_prior_e005",    "prior_rule_traceCB.tsv"), stringsAsFactors = FALSE)
  mt_def <- read.delim(file.path(res_dir, "h2ref_prior_default", "prior_rule_matched.tsv"), stringsAsFactors = FALSE)
  mt_e05 <- read.delim(file.path(res_dir, "h2ref_prior_e005",    "prior_rule_matched.tsv"), stringsAsFactors = FALSE)

  m_keys <- c("coloc.abf (whole region)",
              "block-coloc, per-SNP priors (stages 1-2)",
              "block-coloc, region-equivalent block priors",
              "coloc.signals (mask)",
              "coloc.signals (cond)")
  m_short <- c("coloc.abf", "block-coloc\n(per-SNP)", "block-coloc\n(region-equiv)",
               "coloc.signals\n(mask)", "coloc.signals\n(cond)")

  path <- open_png("fig5_fixed_rule.png", 12.2, 5.4)
  on.exit(dev.off(), add = TRUE)

  layout(matrix(c(1, 2), nrow = 1), widths = c(1.12, 1.0))

  # Panel A: Stacked Fixed-Rule Calling Composition
  par(mar = c(5.0, 4.0, 2.8, 0.9))
  n_m <- length(m_keys)
  xc <- seq_len(n_m)
  plot(NA, xlim = c(0.45, n_m + 0.55), ylim = c(0, 115), xaxt = "n", yaxt = "n",
       xlab = "", ylab = "Share of positive test loci (%)", bty = "n")
  draw_panel_title("A", expression(bold(paste("A   Fixed rule (PP4 > 0.7 / PP3 > 0.7): Default (L) vs ", italic(h)[2]^2 == "5% (R)"))), line = 1.0)
  draw_grid_y(seq(0, 100, by = 20))

  bw <- 0.30
  for (mi in seq_len(n_m)) {
    for (si in 1:2) {
      src <- if (si == 1) tc_def else tc_e05
      row <- src[src$method == m_keys[mi] & src$group == "shared (4 scenarios)", ]
      c_pp4 <- row$colocalised * 100
      c_pp3 <- row$independent * 100
      c_unc <- row$undetermined * 100
      x0 <- xc[mi] + (if (si == 1) -0.20 else 0.20)
      rect(x0 - bw / 2, 0, x0 + bw / 2, c_pp4, col = COL_BC, border = "white", lwd = 0.5)
      rect(x0 - bw / 2, c_pp4, x0 + bw / 2, c_pp4 + c_pp3, col = COL_NEG, border = "white", lwd = 0.5)
      rect(x0 - bw / 2, c_pp4 + c_pp3, x0 + bw / 2, 100, col = "#CBD5E1", border = "white", lwd = 0.5)
      if (c_pp4 >= 9) text(x0, c_pp4 / 2, sprintf("%.0f%%", c_pp4), col = "white", cex = 0.68, font = 2)
      if (c_unc >= 9) text(x0, c_pp4 + c_pp3 + c_unc / 2, sprintf("%.0f%%", c_unc), col = COL_TEXT, cex = 0.66)
      if (si == 1) {
        text(x0, 103.0, "Def.", cex = 0.62, col = COL_MUTED)
      } else {
        text(x0, 103.0, expression(paste(italic(h)[2]^2 == "5%")), cex = 0.62, col = COL_MUTED)
      }
    }
  }
  axis(2, at = seq(0, 100, by = 20), las = 1, cex.axis = 0.85, col = COL_MUTED)
  axis(1, at = xc, labels = m_short, tick = FALSE, line = 0.5, cex.axis = 0.71)
  legend("top",
         legend = c("Colocalised (PP4 > 0.7)", "Independent (PP3 > 0.7)", "Inconclusive"),
         fill = c(COL_BC, COL_NEG, "#CBD5E1"), border = NA, horiz = TRUE,
         bty = "n", cex = 0.74, inset = c(0, -0.02))

  # Panel B: Matched 5% FPR vs Fixed Rule PP4 > 0.7 with FPR annotations
  par(mar = c(5.4, 4.2, 2.8, 0.9))
  p_keys <- m_keys[1:3]
  p_labs <- c("coloc.abf",
              "block-coloc (per-SNP)",
              expression(paste("block-coloc (equiv ", italic(M)/italic(M)[italic(B)], ")")))
  p_cols <- c(COL_ABF, COL_BC, COL_EQPR)

  xg <- 1:4
  plot(NA, xlim = c(0.4, 4.6), ylim = c(0, 124), xaxt = "n", yaxt = "n",
       xlab = "", ylab = "Positive detection rate (%)", bty = "n")
  draw_panel_title("B", "Separating ranking gain from small-block prior odds", line = 1.0)
  draw_grid_y(seq(0, 100, by = 20))
  abline(v = 2.5, col = COL_MUTED, lty = 3, lwd = 1.0)

  bw_b <- 0.255
  offs_b <- c(-1, 0, 1) * bw_b
  for (gi in 1:4) {
    n_pos <- if (gi %in% c(1, 3)) 2400L else 1200L
    for (pi in 1:3) {
      if (gi == 1) p_val <- mt_def$TPR[mt_def$method == p_keys[pi]]
      if (gi == 2) p_val <- mt_e05$TPR[mt_e05$method == p_keys[pi]]
      if (gi == 3) p_val <- tc_def$colocalised[tc_def$method == p_keys[pi] & tc_def$group == "shared (4 scenarios)"]
      if (gi == 4) p_val <- tc_e05$colocalised[tc_e05$method == p_keys[pi] & tc_e05$group == "shared (4 scenarios)"]
      ci <- wilson_ci(round(p_val * n_pos), n_pos)
      x0 <- xg[gi] + offs_b[pi]
      rect(x0 - bw_b * 0.42, 0, x0 + bw_b * 0.42, ci$rate * 100, col = p_cols[pi], border = NA)
      segments(x0, ci$lo * 100, x0, ci$hi * 100, col = COL_TEXT, lwd = 1.2)
      y_bump <- if (gi <= 2 && pi == 2) 4.6 else 1.8
      text(x0, ci$hi * 100 + y_bump, sprintf("%.1f", ci$rate * 100), cex = 0.65, col = COL_TEXT)
    }
  }
  axis(2, at = seq(0, 100, by = 20), las = 1, cex.axis = 0.85, col = COL_MUTED)
  par(xpd = TRUE)
  text(1, -7.5,  "Matched 5%", cex = 0.75, col = COL_TEXT)
  text(1, -13.5, "(Default)",  cex = 0.75, col = COL_TEXT)
  text(2, -7.5,  "Matched 5%", cex = 0.75, col = COL_TEXT)
  text(2, -13.5, expression(paste("(", italic(h)[2]^2 == "5%)")), cex = 0.75, col = COL_TEXT)
  text(3, -7.5,  "PP4 > 0.7",  cex = 0.75, col = COL_TEXT)
  text(3, -13.5, "(Default)",  cex = 0.75, col = COL_TEXT)
  text(4, -7.5,  "PP4 > 0.7",  cex = 0.75, col = COL_TEXT)
  text(4, -13.5, expression(paste("(", italic(h)[2]^2 == "5%)")), cex = 0.75, col = COL_TEXT)

  text(1, -21.5, "FPR: 3.5-3.8%", cex = 0.66, col = COL_MUTED, font = 3)
  text(2, -21.5, "FPR: 5.0-6.8%", cex = 0.66, col = COL_MUTED, font = 3)
  text(3, -21.5, expression(italic(paste("FPR " <= "0.2%"))), cex = 0.66, col = COL_MUTED)
  text(4, -21.5, expression(italic(paste("FPR " <= "0.3%"))), cex = 0.66, col = COL_MUTED)
  par(xpd = FALSE)
  legend("topleft", legend = p_labs, fill = p_cols, border = NA,
         bty = "n", cex = 0.74, inset = c(0.01, -0.01))
  message("Wrote: ", path)
}

# ---------------------------------------------------------------------------
# Supplementary Figure 1: Canonical Names & Plotmath Quantile Distributions
# ---------------------------------------------------------------------------
make_supp_fig1 <- function() {
  qdf <- read.delim(file.path(pt_dir, "supp_pp4_quantiles.tsv"), stringsAsFactors = FALSE)
  pts <- read.delim(gzfile(file.path(pt_dir, "supp_pp4_points.tsv.gz")), stringsAsFactors = FALSE)

  m_order <- c("V3_gate_conditional_fb@1e-03", "whole_abf", "signals_cond", "signals_mask",
               "susie_fallback", "V3_gate_conditional@1e-03", "V2_gate_marginal@1e-05", "susie_raw")
  m_titles <- list(
    expression(bold(paste("block-coloc (", italic(t)[italic(m)] == "0.145)"))),
    expression(bold(paste("coloc.abf (", italic(t)[italic(m)] == "0.141)"))),
    expression(bold(paste("coloc.signals (cond) (", italic(t)[italic(m)] == "0.210)"))),
    expression(bold(paste("coloc.signals (mask) (", italic(t)[italic(m)] == "0.171)"))),
    expression(bold(paste("coloc.susie (fallback) (", italic(t)[italic(m)] == "0.139)"))),
    expression(bold(paste("block-coloc (no fallback) (", italic(t)[italic(m)] == "0.107)"))),
    expression(bold(paste("block-coloc (marginal leads) (", italic(t)[italic(m)] == "0.140)"))),
    expression(bold(paste("coloc.susie (", italic(t)[italic(m)] == "0.011)")))
  )

  sc_names <- c("single_shared", "single_double_shared", "double_shared2", "double_shared1",
                "single_distinct", "double_distinct", "gwas_only", "null_both")
  sc_cols  <- c(rep(COL_POS, 4), rep(COL_NEG, 2), rep(COL_ABF, 2))
  sc_fills <- c(rep("#BBF7D0", 4), rep("#FECACA", 2), rep("#FDE68A", 2))

  set.seed(20261118)
  path <- open_png("supp_fig1_pp4_dist.png", 12.6, 8.2)
  on.exit(dev.off(), add = TRUE)

  par(mfrow = c(2, 4), mar = c(7.2, 3.9, 2.3, 0.7), oma = c(0.6, 0.8, 0.4, 0.4), mgp = c(2.4, 0.7, 0))

  for (mi in seq_along(m_order)) {
    key <- m_order[mi]
    plot(NA, xlim = c(0.4, 8.6), ylim = c(-0.16, 1.03), xaxt = "n", yaxt = "n",
         xlab = "", ylab = if (mi %in% c(1, 5)) "PP4" else "", bty = "n")
    draw_grid_y(seq(0, 1, by = 0.25))
    abline(h = -0.04, col = "#94A3B8", lwd = 0.9)
    abline(v = c(4.5, 6.5), col = COL_GRID, lty = 2, lwd = 0.9)

    qsub <- qdf[qdf$key == key, ]
    cut_v <- qsub$cutoff[1L]
    abline(h = cut_v, col = COL_TEXT, lty = 2, lwd = 1.1)

    for (si in seq_along(sc_names)) {
      qr <- qsub[qsub$scenario == sc_names[si], ]
      pvals <- pts$pp4[pts$key == key & pts$scenario == sc_names[si] & is.finite(pts$pp4)]
      if (length(pvals) > 120L) pvals <- sample(pvals, 120L)
      if (length(pvals) > 0L) {
        points(si + runif(length(pvals), -0.18, 0.18), pvals,
               pch = 16, cex = 0.36, col = adjustcolor(sc_cols[si], alpha.f = 0.20))
      }
      if (is.finite(qr$q50)) {
        segments(si, qr$q05, si, qr$q95, col = sc_cols[si], lwd = 1.2)
        rect(si - 0.24, qr$q25, si + 0.24, qr$q75, col = adjustcolor(sc_fills[si], alpha.f = 0.8),
             border = sc_cols[si], lwd = 1.1)
        segments(si - 0.24, qr$q50, si + 0.24, qr$q50, col = COL_TEXT, lwd = 2.0)
      }
      nc_pct <- round(qr$no_call_rate * 100)
      if (nc_pct > 0) {
        text(si, -0.105, sprintf("%d%%", nc_pct), cex = 0.64, font = 2, col = COL_TEXT)
      } else {
        text(si, -0.105, "0%", cex = 0.60, col = COL_MUTED)
      }
    }
    mtext(m_titles[[mi]], side = 3, line = 0.5, cex = 0.76, col = COL_TEXT)
    axis(2, at = c(-0.105, 0, 0.25, 0.5, 0.75, 1.0),
         labels = c("No call", "0", "0.25", "0.5", "0.75", "1"),
         las = 1, cex.axis = 0.70, col = COL_MUTED)
    if (mi >= 5) {
      axis(1, at = seq_along(sc_names), labels = sc_names, las = 2, cex.axis = 0.70, tick = FALSE, col = COL_TEXT)
    } else {
      axis(1, at = seq_along(sc_names), labels = rep("", length(sc_names)), tick = FALSE)
    }
  }
  message("Wrote: ", path)
}

make_fig1()
make_fig2_stage1()
make_fig2_stage2()
make_fig3()
make_fig4()
make_fig5()
make_supp_fig1()
message("All 7 figures generated successfully.")
