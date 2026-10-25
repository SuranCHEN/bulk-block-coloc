#!/usr/bin/env Rscript
# Bulk block-coloc: evaluation. Reads scores written by run.R; writes
# tables, figures and RESULTS.md. No simulation and no method is run here.
#
# Author: Claude (written for the user's project; executed by Codex).
#
# Usage (from the repository root):
#   Rscript scripts/bulk_gated_v1/evaluate.R --stage smoke
#   Rscript scripts/bulk_gated_v1/evaluate.R --stage main
#   Rscript scripts/bulk_gated_v1/evaluate.R --stage grid    (needs the main stage first)
#
# Statistical rules (fixed before seeing results)
#   Score = PP4 of the method. A missing PP4 (no result) is a negative call
#     and is reported separately as the no-call rate; denominators always
#     include every simulated locus.
#   Matched false-positive rate: on CALIBRATION windows, each method gets the
#     smallest PP4 cut-off t such that the share of negative loci
#     (single_distinct and double_distinct pooled) with PP4 > t is <= 5%.
#   Lead threshold for V0-V4: chosen on calibration windows among 1e-3, 1e-4,
#     1e-5 (same for both traits) as the one with the highest calibration
#     detection rate at matched 5% FPR (ties -> 1e-3). The "_fb" version of a
#     method uses the same lead threshold as the method itself. The earlier
#     bulk setting (GWAS 5e-8, eQTL 1e-5) is reported as a sensitivity only.
#   All reported performance comes from TEST windows with the frozen cut-offs.
#   Intervals: Wilson 95% CI for rates; paired differences use McNemar's exact
#     test and a Wald CI for paired proportions; AUC CI from 1,000 bootstrap
#     resamples of loci (stratified by truth).
#   Localisation (only for methods that return a sub-region): among detected
#     positives, the share whose selected block / SuSiE credible-set pair
#     contains a shared causal SNP, with the median size of that sub-region.

.args <- commandArgs(trailingOnly = TRUE)
.opt <- function(name, default) { i <- match(paste0("--", name), .args)
  if (is.na(i) || i == length(.args)) default else .args[i + 1] }
script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
# Repository root = the nearest parent folder that contains R/block_coloc.R.
# Works in the project layout (scripts/bulk_gated_v1/) and in the public
# repository layout (scripts/ with R/methods.R and R/simulate_bulk.R).
.find_repo <- function(start) { d <- normalizePath(start)
  repeat { if (file.exists(file.path(d, "R", "block_coloc.R"))) return(d)
    p <- dirname(d); if (p == d) stop("cannot find the repository root (folder with R/block_coloc.R)"); d <- p } }
.first_file <- function(...) { f <- c(...); f <- f[file.exists(f)]; if (!length(f)) stop("missing: ", paste(c(...), collapse = " or ")); f[1] }
repo <- .find_repo(dirname(script_file))
stage <- .opt("stage", "main")
out <- .opt("out", file.path(repo, "outputs", "bulk_gated_v1"))
FPR <- as.numeric(.opt("fpr", "0.05"))
B <- as.integer(.opt("boot", "1000"))
set.seed(20261025)

BLOCK0 <- c("V0_handoff", "V1_block_any", "V2_gate_marginal", "V3_gate_conditional", "V4_cond_block")
# "_fb" = same method, but a block branch with no valid block falls back to the
# whole-locus PP4 instead of giving no result.
BLOCK <- c(BLOCK0, paste0(BLOCK0, "_fb"))
OTHER <- c("whole_abf", "susie_raw", "susie_fallback", "signals_mask", "signals_cond")
SCEN <- c("A_A", "A_AB", "AB_AB", "A_B", "AB_CD")
# Scenario names shown in reports (the group's earlier bulk naming; user decision
# 2026-10-25). Internal ids (A_A, ...) stay in the score files.
SCEN_LABEL <- c(A_A = "single_shared", A_AB = "single_double_shared", AB_AB = "double_shared2",
                A_B = "single_distinct", AB_CD = "double_distinct")
SCEN_POS <- c("A_A", "A_AB", "AB_AB"); SCEN_NEG <- c("A_B", "AB_CD")
SYM <- c(1e-3, 1e-4, 1e-5)

# Reported block-coloc versions (user decision 2026-10-25): conditional lead
# selection (proposed method), marginal lead selection, and the original handoff version as a
# control. V1 and V4 are still computed by run.R stage 1 but not reported.
REPORT_BLOCK <- c("V3_gate_conditional", "V2_gate_marginal", "V0_handoff")

# Names shown in tables and figures (internal names stay in thresholds.tsv).
LABEL <- c(V3_gate_conditional = "block-coloc (conditional lead selection, proposed)",
           V2_gate_marginal = "block-coloc (marginal lead selection)",
           V0_handoff = "block-coloc (original handoff)",
           V1_block_any = "block-coloc (always split)", V4_cond_block = "block-coloc (conditioned blocks)",
           whole_abf = "coloc.abf, whole region", susie_raw = "coloc.susie",
           susie_fallback = "coloc.susie + coloc.abf fallback",
           signals_mask = "coloc.signals (mask)", signals_cond = "coloc.signals (cond)")
SHORT <- c(V3_gate_conditional = "block-coloc, conditional (proposed)", V2_gate_marginal = "block-coloc, marginal",
           V0_handoff = "block-coloc, handoff", V1_block_any = "block-coloc, always split", V4_cond_block = "block-coloc, conditioned blocks",
           whole_abf = "coloc.abf", susie_raw = "coloc.susie", susie_fallback = "coloc.susie + fallback",
           signals_mask = "signals (mask)", signals_cond = "signals (cond)")
base_of <- function(k) sub("_fb$", "", sub("@.*", "", k))
legacy_of <- function(k) grepl("@.*/", k)
name_of <- function(k, short = FALSE) { fb <- grepl("_fb$", sub("@.*", "", k))
  paste0((if (short) SHORT else LABEL)[base_of(k)], ifelse(fb, " + whole-region fallback", ""),
         ifelse(legacy_of(k), " [earlier bulk thresholds]", "")) }
lead_of <- function(k) ifelse(grepl("@", k), sub("/", " / ", sub(".*@", "", k)), "-")
label_of <- function(k) name_of(k)
rd <- function(f) utils::read.delim(gzfile(file.path(out, f)), stringsAsFactors = FALSE)
fmt <- function(x, d = 3) ifelse(is.na(x), "NA", formatC(x, format = "f", digits = d))
pct <- function(x) ifelse(is.na(x), "NA", sprintf("%.1f%%", 100 * x))
fe <- function(x) formatC(x, format = "e", digits = 0)
md_table <- function(d) c(paste0("| ", paste(names(d), collapse = " | "), " |"),
                          paste0("|", paste(rep("---", ncol(d)), collapse = "|"), "|"),
                          apply(d, 1, function(r) paste0("| ", paste(r, collapse = " | "), " |")), "")
wilson <- function(k, n) { if (n == 0) return(c(NA, NA)); p <- k / n; z <- 1.959964
  c((p + z^2 / (2 * n) - z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2))) / (1 + z^2 / n),
    (p + z^2 / (2 * n) + z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2))) / (1 + z^2 / n)) }
ci_txt <- function(x) { x <- as.logical(x); ci <- wilson(sum(x), length(x)); sprintf("[%s, %s]", pct(ci[1]), pct(ci[2])) }
rate_ci <- function(x) { x <- as.logical(x); sprintf("%s %s", pct(mean(x)), ci_txt(x)) }
auc <- function(s, y) { r <- rank(s); n1 <- sum(y == 1); n0 <- sum(y == 0)
  (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0) }
score_of <- function(pp4) ifelse(is.finite(pp4), pp4, -Inf)
key_of <- function(d) ifelse(d$method %in% BLOCK,
  ifelse(d$lead_threshold1 == d$lead_threshold2, paste0(d$method, "@", fe(d$lead_threshold1)),
         paste0(d$method, "@", fe(d$lead_threshold1), "/", fe(d$lead_threshold2))), d$method)
# The three-line conclusion is written by Claude into conclusions_<stage>.md
# (in the output folder) and inserted here, so re-running this script keeps it.
conclusions <- function(st) { f <- file.path(out, paste0("conclusions_", st, ".md"))
  if (file.exists(f)) readLines(f, warn = FALSE) else "(to be written by Claude after reviewing the tables; Codex leaves this section empty)" }
cutoff <- function(s_neg, fpr) { s <- sort(s_neg, decreasing = TRUE); k <- floor(fpr * length(s)); if (k >= length(s)) -Inf else s[k + 1] }

load_scores <- function(st) {
  d <- rd(paste0("scores_", st, ".tsv.gz"))
  fb <- d[d$method %in% BLOCK0, ]
  fb$method <- paste0(fb$method, "_fb"); fb$pp4 <- ifelse(is.finite(fb$pp4), fb$pp4, fb$pp4_whole)
  d <- rbind(d, fb)
  d$key <- key_of(d); d$score <- score_of(d$pp4); d$no_call <- !is.finite(d$pp4)
  d$id <- paste(d$window, d$scenario, d$rep, d$q_gwas, d$q_eqtl, d$tier, sep = "|")
  d
}

pooled <- function(x, boot = TRUE) {
  pos <- x[x$truth == 1, ]; neg <- x[x$truth == 0, ]
  yy <- c(rep(1, nrow(pos)), rep(0, nrow(neg))); ss <- c(pos$score, neg$score)
  a <- auc(ss, yy)
  ci <- if (boot) quantile(replicate(B, { i1 <- sample(which(yy == 1), replace = TRUE); i0 <- sample(which(yy == 0), replace = TRUE)
    auc(ss[c(i1, i0)], yy[c(i1, i0)]) }), c(.025, .975)) else c(NA, NA)
  data.frame(n_pos = nrow(pos), TPR = rate_ci(pos$call), n_neg = nrow(neg), FPR = rate_ci(neg$call),
             AUC = if (boot) sprintf("%s [%s, %s]", fmt(a), fmt(ci[1]), fmt(ci[2])) else fmt(a),
             no_call = pct(mean(x$no_call)), branch_block = pct(mean(x$branch == "block", na.rm = TRUE)))
}

paired <- function(ts, keys, ref) do.call(rbind, lapply(setdiff(keys, ref), function(k) {
  one <- function(truth) {
    a <- ts[ts$key == k & ts$truth == truth, c("id", "call")]; b <- ts[ts$key == ref & ts$truth == truth, c("id", "call")]
    m <- merge(a, b, by = "id"); n <- nrow(m)
    b10 <- sum(m$call.x & !m$call.y); b01 <- sum(!m$call.x & m$call.y)
    dlt <- (b10 - b01) / n; se <- sqrt(b10 + b01 - (b10 - b01)^2 / n) / n
    p <- if (b10 + b01 == 0) 1 else stats::binom.test(b10, b10 + b01)$p.value
    sprintf("%+.1f pts [%+.1f, %+.1f], p=%s", 100 * dlt, 100 * (dlt - 1.96 * se), 100 * (dlt + 1.96 * se), format.pval(p, digits = 2))
  }
  data.frame(method = name_of(k), vs = name_of(ref), TPR_difference = one(1), FPR_difference = one(0))
}))

gate_table <- function(ts, keys) do.call(rbind, lapply(keys, function(k) do.call(rbind, lapply(SCEN, function(s) {
  x <- ts[ts$key == k & ts$scenario == s, ]
  if (!nrow(x)) return(NULL)
  multi <- max(x$n_causal_gwas[1], x$n_causal_eqtl[1]) >= 2
  data.frame(method = name_of(k), scenario = SCEN_LABEL[s], true_signals = paste0(x$n_causal_gwas[1], "+", x$n_causal_eqtl[1]),
             should_go = if (multi) "block" else "whole", n = nrow(x),
             branch_block = pct(mean(x$branch == "block")),
             correct_branch = pct(mean(x$branch == if (multi) "block" else "whole")),
             leads_match_truth = pct(mean(x$n_leads1 == x$n_causal_gwas & x$n_leads2 == x$n_causal_eqtl)),
             median_leads = sprintf("%s+%s", fmt(stats::median(x$n_leads1), 0), fmt(stats::median(x$n_leads2), 0)))
}))))

# ------------------------------------------------------------ smoke ----

evaluate_smoke <- function() {
  d <- rd("scores_smoke.tsv.gz"); lc <- rd("loci_smoke.tsv.gz")
  med <- function(v) stats::median(as.numeric(unlist(strsplit(as.character(v), ";"))), na.rm = TRUE)
  rng <- function(v) { x <- as.numeric(unlist(strsplit(as.character(v), ";"))); if (!length(x) || all(is.na(x))) "-" else sprintf("%.2f-%.2f", min(x, na.rm = TRUE), max(x, na.rm = TRUE)) }
  cells <- unique(lc[, c("scenario", "q_gwas", "q_eqtl", "tier")])
  san <- do.call(rbind, lapply(seq_len(nrow(cells)), function(i) {
    x <- merge(lc, cells[i, ])
    data.frame(scenario = SCEN_LABEL[cells$scenario[i]], q_gwas = cells$q_gwas[i], q_eqtl = cells$q_eqtl[i], tier = cells$tier[i], n = nrow(x),
               causal_z_gwas = fmt(med(x$causal_abs_z_gwas), 1), expected_gwas = fmt(sqrt(20000 * cells$q_gwas[i]), 1),
               causal_z_eqtl = fmt(med(x$causal_abs_z_eqtl), 1), expected_eqtl = fmt(sqrt(1000 * cells$q_eqtl[i]), 1),
               pair_abs_r = rng(x$pair_abs_r), var_y = sprintf("%s / %s", fmt(stats::median(x$var_y_gwas), 2), fmt(stats::median(x$var_y_eqtl), 2)))
  }))
  san <- san[order(san$tier, san$scenario, san$q_gwas, san$q_eqtl), ]
  st <- as.data.frame(table(method = d$method, status = d$status))
  st <- st[st$Freq > 0, ]
  errf <- file.path(out, "errors_smoke.txt")
  L <- c("# Smoke test (2 windows x 5 scenarios x 2 replicates, every cell)", "",
         sprintf("Generated %s. Score rows: %d. Loci: %d.", format(Sys.time(), "%Y-%m-%d %H:%M"), nrow(d), nrow(lc)), "",
         "## Simulation sanity", "",
         "Causal |z| (median over loci and causal SNPs) should be close to sqrt(N q). pair_abs_r must lie inside the tier",
         "(low < 0.2; mid 0.32-0.63; high 0.63-0.9). var_y is the median trait variance (GWAS / eQTL), close to 1.", "",
         md_table(san), "## Method status counts", "", md_table(st),
         "## Errors", "", if (file.exists(errf)) readLines(errf) else "None.")
  writeLines(L, file.path(out, "RESULTS_smoke.md"))
  message("wrote ", file.path(out, "RESULTS_smoke.md"))
}

# ------------------------------------------------------------- main ----

evaluate_main <- function() {
  d <- load_scores("main")
  loci <- rd("loci_main.tsv.gz")
  loci$id <- paste(loci$window, loci$scenario, loci$rep, loci$q_gwas, loci$q_eqtl, loci$tier, sep = "|")
  cal <- d[d$set == "calibration", ]; tst <- d[d$set == "test", ]
  keys <- unique(d$key)
  th <- do.call(rbind, lapply(keys, function(k) {
    x <- cal[cal$key == k, ]
    t <- cutoff(x$score[x$truth == 0], FPR)
    data.frame(key = k, method = x$method[1], lead_threshold1 = x$lead_threshold1[1], lead_threshold2 = x$lead_threshold2[1],
               cutoff = t, calib_fpr = mean(x$score[x$truth == 0] > t), calib_tpr = mean(x$score[x$truth == 1] > t),
               calib_n_neg = sum(x$truth == 0), calib_n_pos = sum(x$truth == 1))
  }))
  th$legacy <- th$method %in% BLOCK & th$lead_threshold1 != th$lead_threshold2
  # lead threshold choice for V0-V4 (symmetric thresholds only); _fb follows its method
  th$selected <- !th$method %in% BLOCK
  for (m in BLOCK0) {
    x <- th[th$method == m & !th$legacy, ]; best <- max(x$calib_tpr)
    pick <- x[x$calib_tpr == best, ]
    pick <- if (any(pick$lead_threshold1 == 1e-3)) pick[pick$lead_threshold1 == 1e-3, ][1, ] else pick[1, ]
    th$selected[th$method %in% c(m, paste0(m, "_fb")) & th$lead_threshold1 == pick$lead_threshold1 & !th$legacy] <- TRUE
  }
  write.table(th, file.path(out, "thresholds.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  sel <- th$key[th$selected]
  MAIN_KEYS <- c(sel[match(REPORT_BLOCK, sub("@.*", "", sel))], OTHER)
  FB_KEYS <- sel[match(paste0(REPORT_BLOCK, "_fb"), sub("@.*", "", sel))]
  LEG_KEYS <- th$key[th$legacy & th$method %in% REPORT_BLOCK]
  tst$call <- tst$score > th$cutoff[match(tst$key, th$key)]
  ts <- tst[tst$key %in% c(MAIN_KEYS, FB_KEYS, LEG_KEYS), ]
  tm <- ts[ts$key %in% MAIN_KEYS, ]

  L <- c("# Bulk block-coloc: main comparison (stage 1)", "",
         sprintf("Generated %s. Setting: GWAS N 20,000, q 0.3%% per causal SNP; eQTL N 1,000, q 1.5%% per causal SNP; causal |r| < 0.2.",
                 format(Sys.time(), "%Y-%m-%d %H:%M")),
         sprintf("Calibration windows: %s. Test windows: %s.", paste(sort(unique(cal$window)), collapse = ", "), paste(sort(unique(tst$window)), collapse = ", ")),
         sprintf("Loci: %d calibration, %d test (per method). Target FPR %.0f%% on calibration windows (negatives: single_distinct and double_distinct).",
                 length(unique(cal$id)), length(unique(tst$id)), 100 * FPR), "",
         "## Three-line conclusion", "", conclusions("main"), "")

  # 1. main table in the AGENTS.md format: scenario x method, test windows
  mt <- do.call(rbind, lapply(SCEN, function(s) do.call(rbind, lapply(MAIN_KEYS, function(k) {
    x <- tm[tm$scenario == s & tm$key == k, ]
    data.frame(scenario = SCEN_LABEL[s], truth = if (s %in% SCEN_POS) "shared" else "not shared", method = name_of(k), lead_p = lead_of(k),
               n_total = nrow(x), n_positive = sum(x$call), rate = pct(mean(x$call)), CI95 = ci_txt(x$call),
               no_call_rate = pct(mean(x$no_call)), branch_block = pct(mean(x$branch == "block", na.rm = TRUE)))
  }))))
  names(mt)[names(mt) == "CI95"] <- "95% CI"
  write.table(mt, file.path(out, "main_table.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  L <- c(L, "## 1. Main table (test windows, frozen cut-offs)", "",
         "n_positive = loci called colocalised (PP4 above the method's cut-off). For shared scenarios rate is the detection rate;",
         "for not-shared scenarios it is the false-positive rate. Denominators include every locus; no result counts as not called.", "",
         md_table(mt))

  # 2. pooled
  pt <- do.call(rbind, lapply(MAIN_KEYS, function(k) cbind(data.frame(method = name_of(k), lead_p = lead_of(k)), pooled(tm[tm$key == k, ]))))
  write.table(pt, file.path(out, "pooled_table.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  L <- c(L, "## 2. Pooled over scenarios (test windows)", "",
         "TPR pools single_shared, single_double_shared, double_shared2; FPR pools single_distinct and double_distinct. AUC does not depend on any cut-off.", "", md_table(pt))

  # 3. cut-offs
  tt <- th[base_of(th$key) %in% c(REPORT_BLOCK, OTHER), ]
  tt <- tt[order(match(base_of(tt$key), c(REPORT_BLOCK, OTHER)), grepl("_fb@", tt$key), tt$legacy, -tt$lead_threshold1), ]
  L <- c(L, "## 3. Cut-offs chosen on calibration windows", "",
         "lead_p = lead threshold (GWAS / eQTL when they differ). '*' marks the lead threshold kept for each block-coloc version.", "",
         md_table(data.frame(method = name_of(tt$key), lead_p = paste0(lead_of(tt$key), ifelse(tt$selected & tt$method %in% BLOCK, " *", "")),
                             cutoff = ifelse(is.finite(tt$cutoff), fmt(tt$cutoff, 4), "-Inf"),
                             calib_FPR = pct(tt$calib_fpr), calib_TPR = pct(tt$calib_tpr))))

  # 4. paired comparisons
  pc <- rbind(paired(tm, MAIN_KEYS, "susie_fallback"), paired(tm, MAIN_KEYS, "whole_abf"))
  L <- c(L, "## 4. Paired differences on the same test loci (McNemar)", "",
         "Difference = method minus reference, in percentage points, with 95% CI and exact McNemar p-value.", "", md_table(pc))

  # 5. branch choice
  gk <- MAIN_KEYS[base_of(MAIN_KEYS) %in% BLOCK0]
  gt <- gate_table(tm, gk)
  write.table(gt, file.path(out, "gate_table.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  L <- c(L, "## 5. Branch choice: split into blocks or not (test windows)", "",
         "should_go: block when either trait has 2 causal SNPs, whole otherwise. leads_match_truth: the method's lead counts",
         "(GWAS, eQTL) equal the true numbers of causal SNPs.", "",
         md_table(gt))

  # 6. single_double_shared split by which eQTL signal is stronger
  ab <- merge(tm[tm$scenario == "A_AB", ], loci[, c("id", "eqtl_shared_stronger")], by = "id")
  sp <- do.call(rbind, lapply(MAIN_KEYS, function(k) { x <- ab[ab$key == k, ]
    data.frame(method = name_of(k), n_shared_stronger = sum(x$eqtl_shared_stronger), TPR_shared_stronger = rate_ci(x$call[x$eqtl_shared_stronger]),
               n_other_stronger = sum(!x$eqtl_shared_stronger), TPR_other_stronger = rate_ci(x$call[!x$eqtl_shared_stronger])) }))
  L <- c(L, "## 6. single_double_shared split by which eQTL signal is stronger (test windows)", "",
         "The eQTL has two causal variants: one shared with the GWAS and one eQTL-only. 'stronger' compares their observed eQTL |z|.", "", md_table(sp))

  # 7. localisation
  loc_tab <- do.call(rbind, lapply(MAIN_KEYS, function(k) {
    x <- tm[tm$key == k & tm$truth == 1 & tm$call, ]
    u <- x[!is.na(x$unit_has_shared), ]
    data.frame(method = name_of(k), detected = nrow(x), with_subregion = nrow(u),
               contains_shared = if (nrow(u)) rate_ci(u$unit_has_shared) else "NA",
               median_size_snps = if (nrow(u)) fmt(stats::median(u$unit_size), 0) else "NA")
  }))
  L <- c(L, "## 7. Localisation among detected test positives", "",
         "Sub-region = best block (block branch) or the two credible sets of the best SuSiE signal pair.",
         "Whole-locus calls and coloc.signals calls have no sub-region and are not counted here.", "", md_table(loc_tab))

  # 8. no-valid-block policy and legacy lead thresholds
  ft <- do.call(rbind, lapply(c(MAIN_KEYS[base_of(MAIN_KEYS) %in% BLOCK0], FB_KEYS, LEG_KEYS), function(k)
    cbind(data.frame(method = name_of(k), lead_p = lead_of(k), cutoff = fmt(th$cutoff[th$key == k], 4)), pooled(ts[ts$key == k, ], boot = FALSE))))
  L <- c(L, "## 8. Sensitivity: fallback policy and the earlier bulk lead thresholds (test windows)", "",
         "'+ whole-region fallback' = when the block branch builds no valid block, use the whole-region PP4 instead of no result.",
         "'[earlier bulk thresholds]' = lead thresholds GWAS 5e-8, eQTL 1e-5 (August bulk setting), each with its own calibrated cut-off.", "", md_table(ft))

  # 9. inputs
  lt <- loci[loci$set == "test", ]
  med <- function(v) stats::median(as.numeric(unlist(strsplit(as.character(v), ";"))), na.rm = TRUE)
  inp <- do.call(rbind, lapply(SCEN, function(s) { x <- lt[lt$scenario == s, ]
    data.frame(scenario = SCEN_LABEL[s], n = nrow(x), causal_z_gwas = fmt(med(x$causal_abs_z_gwas), 1),
               causal_z_eqtl = fmt(med(x$causal_abs_z_eqtl), 1), pair_abs_r_median = fmt(med(x$pair_abs_r), 2),
               gwas_no_cs = pct(mean(x$n_cs_gwas == 0, na.rm = TRUE)), eqtl_no_cs = pct(mean(x$n_cs_eqtl == 0, na.rm = TRUE))) }))
  L <- c(L, "## 9. Inputs (test windows)", "",
         sprintf("Expected causal |z| about sqrt(N q): GWAS %.1f, eQTL %.1f.", sqrt(20000 * 0.003), sqrt(1000 * 0.015)), "", md_table(inp))

  # figures
  cols <- grDevices::hcl.colors(length(MAIN_KEYS), "Dark 3")
  png(file.path(out, "matched_fpr_by_scenario.png"), width = 1600, height = 1000, res = 130)
  m <- sapply(MAIN_KEYS, function(k) sapply(SCEN, function(s) mean(tm$call[tm$key == k & tm$scenario == s])))
  par(mar = c(4, 4, 3, 1))
  barplot(t(m), beside = TRUE, col = cols, names.arg = SCEN_LABEL[SCEN], ylim = c(0, 1.32), axes = FALSE,
          ylab = "Share called colocalised", main = "Test windows at matched 5% FPR (left 3: shared; right 2: not shared)")
  axis(2, at = seq(0, 1, 0.2))
  abline(v = 3 * (length(MAIN_KEYS) + 1) + 0.5, lty = 2); abline(h = FPR, lty = 3)
  legend("top", legend = label_of(MAIN_KEYS), fill = cols, bty = "n", cex = 0.7, ncol = 2)
  dev.off()
  png(file.path(out, "roc_test.png"), width = 1100, height = 1000, res = 140)
  plot(NA, xlim = c(0, 1), ylim = c(0, 1), xlab = "False-positive rate", ylab = "Detection rate",
       main = "ROC on test windows (main setting)"); abline(0, 1, lty = 3, col = "grey60"); abline(v = FPR, lty = 2, col = "grey60")
  for (i in seq_along(MAIN_KEYS)) { x <- tm[tm$key == MAIN_KEYS[i], ]
    th_grid <- sort(unique(c(Inf, x$score, -Inf)), decreasing = TRUE)
    s1 <- x$score[x$truth == 1]; s0 <- x$score[x$truth == 0]
    lines(vapply(th_grid, function(t) mean(s0 >= t), 0), vapply(th_grid, function(t) mean(s1 >= t), 0),
          col = cols[i], lwd = 2, lty = if (base_of(MAIN_KEYS[i]) %in% BLOCK0) 1 else 2) }
  legend("bottomright", legend = label_of(MAIN_KEYS), col = cols, lwd = 2,
         lty = ifelse(base_of(MAIN_KEYS) %in% BLOCK0, 1, 2), bty = "n", cex = 0.6)
  dev.off()
  L <- c(L, "## Figures", "", "- Main figure: `matched_fpr_by_scenario.png`", "- `roc_test.png`", "",
         "## Reproduce", "", "```",
         "Rscript scripts/bulk_gated_v1/run.R --stage main --cores 4 && Rscript scripts/bulk_gated_v1/evaluate.R --stage main",
         "```", "", "## Files changed (git diff --stat)", "", "(Codex pastes the output of `git diff --stat HEAD~1` here after committing)")
  writeLines(L, file.path(out, "RESULTS.md"))
  message("wrote ", file.path(out, "RESULTS.md"))
}

# ------------------------------------------------------------- grid ----

evaluate_grid <- function() {
  th <- utils::read.delim(file.path(out, "thresholds.tsv"), stringsAsFactors = FALSE)
  sel <- th$key[th$selected & !grepl("_fb@", th$key)]
  g <- load_scores("grid"); m <- load_scores("main"); m <- m[m$rep <= max(g$rep), ]
  d <- rbind(g, m[, names(g)])
  d <- d[d$set == "test" & d$key %in% sel, ]
  d$call <- d$score > th$cutoff[match(d$key, th$key)]
  keys <- c(sel[match(REPORT_BLOCK, sub("@.*", "", sel))], OTHER); keys <- keys[!is.na(keys)]
  one <- function(x) data.frame(TPR = mean(x$call[x$truth == 1]), FPR = mean(x$call[x$truth == 0]),
                                AUC = if (length(unique(x$truth)) == 2) auc(x$score, x$truth) else NA,
                                n_pos = sum(x$truth == 1), n_neg = sum(x$truth == 0))
  cells <- unique(d[, c("q_gwas", "q_eqtl", "tier")])
  tab <- do.call(rbind, lapply(seq_len(nrow(cells)), function(i) do.call(rbind, lapply(keys, function(k) {
    x <- d[d$q_gwas == cells$q_gwas[i] & d$q_eqtl == cells$q_eqtl[i] & d$tier == cells$tier[i] & d$key == k, ]
    cbind(cells[i, ], data.frame(method = k, name = name_of(k)), one(x), row.names = NULL) }))))
  write.table(tab, file.path(out, "grid_table.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  wide <- function(sub, vary, what) { w <- reshape(sub[, c(vary, "method", what)], idvar = vary, timevar = "method", direction = "wide")
    nm <- sub(paste0("^", what, "\\."), "", names(w)); names(w) <- ifelse(nm == vary, vary, name_of(nm, short = TRUE))
    for (j in 2:ncol(w)) w[[j]] <- pct(w[[j]]); w[order(w[[vary]]), ] }
  qg <- tab[tab$q_eqtl == 0.015 & tab$tier == "low", ]
  qe <- tab[tab$q_gwas == 0.003 & tab$tier == "low", ]
  # LD tiers: pooled over the scenarios that have a causal pair (single_shared is not run for mid / high)
  dt <- d[d$q_gwas == 0.003 & d$q_eqtl == 0.015 & d$scenario != "A_A", ]
  tier_tab <- do.call(rbind, lapply(c("low", "mid", "high"), function(t) do.call(rbind, lapply(keys, function(k)
    cbind(data.frame(tier = t, method = k), one(dt[dt$tier == t & dt$key == k, ]))))))
  tier_tab$tier <- factor(tier_tab$tier, c("low", "mid", "high"))
  by_sc <- do.call(rbind, lapply(c("low", "mid", "high"), function(t) do.call(rbind, lapply(setdiff(SCEN, "A_A"), function(s) {
    x <- dt[dt$tier == t & dt$scenario == s, ]
    cbind(data.frame(tier = t, scenario = SCEN_LABEL[s], n = sum(x$key == keys[1])),
          data.frame(setNames(lapply(keys, function(k) pct(mean(x$call[x$key == k]))), name_of(keys, short = TRUE)), check.names = FALSE)) }))))
  gk <- keys[base_of(keys) %in% REPORT_BLOCK]
  gt <- do.call(rbind, lapply(c("low", "mid", "high"), function(t) cbind(data.frame(tier = t), gate_table(dt[dt$tier == t, ], gk))))
  L <- c("# Bulk block-coloc: stage 2, one factor at a time (test windows, frozen cut-offs)", "",
         sprintf("Generated %s. Replicates per cell: %d. Cut-offs and lead thresholds from stage 1 (thresholds.tsv), not re-tuned.",
                 format(Sys.time(), "%Y-%m-%d %H:%M"), max(d$rep)), "",
         "## Three-line conclusion", "", conclusions("grid"), "",
         "## Detection rate vs GWAS q (eQTL q 1.5%, low LD)", "", md_table(wide(qg, "q_gwas", "TPR")),
         "## False-positive rate vs GWAS q", "", md_table(wide(qg, "q_gwas", "FPR")),
         "## Detection rate vs eQTL q (GWAS q 0.3%, low LD)", "", md_table(wide(qe, "q_eqtl", "TPR")),
         "## False-positive rate vs eQTL q", "", md_table(wide(qe, "q_eqtl", "FPR")),
         "## LD between causal SNPs (main strengths; single_double_shared, double_shared2 pooled for TPR; single_distinct, double_distinct for FPR)", "",
         "low |r| < 0.2; mid 0.32-0.63; high 0.63-0.9.", "",
         md_table(wide(tier_tab, "tier", "TPR")), md_table(wide(tier_tab, "tier", "FPR")),
         "### By scenario", "", md_table(by_sc), "### Branch choice by LD tier", "", md_table(gt))
  cols <- grDevices::hcl.colors(length(keys), "Dark 3")
  png(file.path(out, "grid_lines.png"), width = 1800, height = 1100, res = 130)
  par(mfrow = c(2, 3))
  specs <- list(list(qg, "q_gwas", "TPR"), list(qe, "q_eqtl", "TPR"), list(tier_tab, "tier", "TPR"),
                list(qg, "q_gwas", "FPR"), list(qe, "q_eqtl", "FPR"), list(tier_tab, "tier", "FPR"))
  for (sp in specs) {
    x <- sp[[1]]; v <- sp[[2]]; y <- sp[[3]]
    xv <- if (v == "tier") as.integer(x$tier) else x[[v]]
    plot(NA, xlim = range(xv), ylim = c(0, if (y == "TPR") 1 else max(0.2, max(x[[y]], na.rm = TRUE))),
         log = if (v != "tier") "x" else "", xaxt = if (v == "tier") "n" else "s", xlab = v, ylab = y, main = paste(y, "vs", v))
    if (v == "tier") axis(1, 1:3, c("low", "mid", "high"))
    for (i in seq_along(keys)) { z <- x[x$method == keys[i], ]; zx <- if (v == "tier") as.integer(z$tier) else z[[v]]
      o <- order(zx); lines(zx[o], z[[y]][o], type = "b", col = cols[i], pch = 16, lty = if (base_of(keys[i]) %in% BLOCK0) 1 else 2) }
    if (y == "FPR") abline(h = FPR, lty = 3)
  }
  legend("topright", legend = label_of(keys), col = cols, lty = ifelse(base_of(keys) %in% BLOCK0, 1, 2), bty = "n", cex = .55)
  dev.off()
  writeLines(c(L, "Figure: `grid_lines.png`"), file.path(out, "RESULTS_grid.md"))
  message("wrote ", file.path(out, "RESULTS_grid.md"))
}

switch(stage, smoke = evaluate_smoke(), main = evaluate_main(), grid = evaluate_grid(),
       stop("--stage must be smoke, main or grid"))
