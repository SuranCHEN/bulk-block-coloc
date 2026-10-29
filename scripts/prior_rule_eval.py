#!/usr/bin/env python3
"""Evaluate prior_rule.R output: block priors x decision rules.

Methods (block-coloc uses the adopted whole-region fallback when no valid block):
  whole_abf     whole-region coloc.abf
  bc_orig       block-coloc, per-SNP priors inside blocks (as in stages 1 and 2)
  bc_adj        block-coloc, block priors scaled so that H0-H4 have the whole-region prior mass
  signals_mask, signals_cond
Decision rules:
  matched   cut-off on PP4 chosen on calibration windows (single_distinct + double_distinct)
            so that the false-positive rate is 5%; reported on test windows
  traceCB   fixed rule of Jiang et al. (2026): colocalised PP4 > 0.7, independent PP3 > 0.7,
            otherwise undetermined; reported on test windows (no calibration)
Optional check (--stage1 <dir>): PP4 of whole_abf, bc_orig and coloc.signals must equal
the stage-1 (scores_main) or stage-2 (scores_grid) scores of the same cell (same loci, same seeds).

Usage: python3 scripts/bulk_h2ref_v1/prior_rule_eval.py --out outputs/bulk_h2ref_v1_prior [--stage1 outputs/bulk_h2ref_v1]
Writes prior_rule_matched.tsv, prior_rule_traceCB.tsv, prior_rule_paired.tsv and PRIOR_RULE.md in --out.
Needs numpy, pandas, scipy.
"""
import argparse, os
import numpy as np
import pandas as pd
from scipy.stats import binomtest

ap = argparse.ArgumentParser()
ap.add_argument("--out", required=True)
ap.add_argument("--stage1", default=None)
a = ap.parse_args()
d = pd.read_csv(os.path.join(a.out, "prior_rule_loci.tsv.gz"), sep="\t")

POS = ["A_A", "A_AB", "AB_AB", "AB_AC"]; NEG = ["A_B", "AB_CD"]
LABEL = {"A_A": "single_shared", "A_AB": "single_double_shared", "AB_AB": "double_shared2", "AB_AC": "double_shared1",
         "A_B": "single_distinct", "AB_CD": "double_distinct", "gwas_only": "gwas_only", "null_both": "null_both"}
NAME = {"whole_abf": "coloc.abf (whole region)", "bc_orig": "block-coloc, per-SNP priors (stages 1-2)",
        "bc_adj": "block-coloc, region-equivalent block priors", "signals_mask": "coloc.signals (mask)",
        "signals_cond": "coloc.signals (cond)"}
M = list(NAME)

def pp(meth):
    if meth == "whole_abf": return d.whole_pp3, d.whole_pp4
    if meth in ("bc_orig", "bc_adj"):   # fallback: whole-region values when there is no block result
        p3, p4 = d[f"{meth}_pp3"], d[f"{meth}_pp4"]
        use = p4.isna()
        return p3.where(~use, d.whole_pp3), p4.where(~use, d.whole_pp4)
    k = meth.replace("signals_", "sig_")
    return d[f"{k}_pp3"], d[f"{k}_pp4"]

for m in M:
    p3, p4 = pp(m); d[f"{m}_p3"] = p3; d[f"{m}_p4"] = p4
key = ["window", "scenario", "rep"]
cal, tst = d[d.set == "calibration"], d[d.set == "test"]

def cutoff(neg_scores, fpr=0.05):
    s = np.sort(neg_scores)[::-1]; k = int(np.floor(fpr * len(s)))
    return -np.inf if k >= len(s) else s[k]

def auc(s, y):
    r = pd.Series(s).rank().values; n1 = (y == 1).sum(); n0 = (y == 0).sum()
    return (r[y == 1].sum() - n1 * (n1 + 1) / 2) / (n1 * n0)

def mcnemar(A, B):
    b10, b01, n = int((A & ~B).sum()), int((~A & B).sum()), len(A)
    dl = (b10 - b01) / n; se = np.sqrt(max((b10 + b01) / n - dl ** 2, 0) / n)
    return dl, dl - 1.96 * se, dl + 1.96 * se, (binomtest(b10, b10 + b01, 0.5).pvalue if b10 + b01 else 1.0)

# ---- rule 1: matched 5% FPR
rows, calls = [], {}
for m in M:
    sc = lambda x: x[f"{m}_p4"].fillna(-1).values
    cn = cal[cal.scenario.isin(NEG)]
    t = cutoff(sc(cn))
    c = tst[f"{m}_p4"].fillna(-1) > t
    calls[m] = pd.Series(c.values, index=pd.MultiIndex.from_frame(tst[key]))
    tc = tst[tst.scenario.isin(POS + NEG)]
    r = dict(method=NAME[m], cutoff=round(t, 4),
             TPR=c[tst.scenario.isin(POS)].mean(), FPR=c[tst.scenario.isin(NEG)].mean(),
             gwas_only=c[tst.scenario == "gwas_only"].mean(), null_both=c[tst.scenario == "null_both"].mean(),
             AUC=auc(tc[f"{m}_p4"].fillna(-1).values, tc.scenario.isin(POS).astype(int).values))
    for s in POS + NEG: r[LABEL[s]] = c[tst.scenario == s].mean()
    rows.append(r)
matched = pd.DataFrame(rows)

# ---- rule 2: traceCB fixed rule
rows2 = []
for m in M:
    p3, p4 = tst[f"{m}_p3"], tst[f"{m}_p4"]
    cls = np.where(p4 > 0.7, "colocalised", np.where(p3 > 0.7, "independent", "undetermined"))
    for grp, scen in [("shared (4 scenarios)", POS), ("single_distinct", ["A_B"]), ("double_distinct", ["AB_CD"]),
                      ("gwas_only", ["gwas_only"]), ("null_both", ["null_both"])]:
        sel = tst.scenario.isin(scen).values
        rows2.append(dict(method=NAME[m], group=grp, n=int(sel.sum()),
                          colocalised=(cls[sel] == "colocalised").mean(), independent=(cls[sel] == "independent").mean(),
                          undetermined=(cls[sel] == "undetermined").mean()))
trace = pd.DataFrame(rows2)

# ---- paired comparisons (positives) under both rules
fixed = {m: pd.Series((tst[f"{m}_p4"] > 0.7).values, index=pd.MultiIndex.from_frame(tst[key])) for m in M}
pos_idx = calls["whole_abf"].index[calls["whole_abf"].index.get_level_values("scenario").isin(POS)]
neg_idx = calls["whole_abf"].index[calls["whole_abf"].index.get_level_values("scenario").isin(NEG)]
rows3 = []
for rule, C in [("matched 5% FPR", calls), ("traceCB PP4 > 0.7", fixed)]:
    for x, y in [("bc_adj", "whole_abf"), ("bc_orig", "whole_abf"), ("bc_adj", "bc_orig"), ("bc_adj", "signals_cond"), ("bc_adj", "signals_mask")]:
        for part, idx in [("TPR", pos_idx), ("FPR", neg_idx)]:
            dl, lo, hi, p = mcnemar(C[x].loc[idx].values, C[y].loc[idx].values)
            rows3.append(dict(rule=rule, method=NAME[x], vs=NAME[y], part=part, diff=dl, lo=lo, hi=hi, p=p))
paired = pd.DataFrame(rows3)

# ---- prior diagnostics on loci that used a block
blk = d[d.branch == "block"]
diag = dict(block_loci=len(blk), median_block_snps=float(blk.bc_orig_size.median()), median_region_snps=float(d.n_snps.median()),
            median_pp4_orig=float(blk.bc_orig_pp4.median()), median_pp4_adj=float(blk.bc_adj_pp4.median()),
            median_pp4_whole=float(blk.whole_pp4.median()))

# ---- optional consistency check with stage 1
check = ""
if a.stage1:
    s1 = None
    for f in ("scores_main.tsv.gz", "scores_grid.tsv.gz"):   # the cell is in stage 1 (main) or stage 2 (grid)
        fp = os.path.join(a.stage1, f)
        if not os.path.exists(fp): continue
        x = pd.read_csv(fp, sep="\t")
        x = x[(x.tier == "low") & np.isclose(x.h2_gwas, d.h2_gwas.iloc[0]) & np.isclose(x.h2_eqtl, d.h2_eqtl.iloc[0])]
        if len(x): s1 = x; break
    if s1 is None: raise SystemExit("no stage-1/2 scores for this cell in " + a.stage1)
    out = []
    # V3_gate_conditional (no fallback) in stage 1: whole-region PP4 when it does not split, the best
    # block's PP4 when it splits, NA when no valid block. Rebuild the same quantity from this run.
    d["v3_equiv_pp4"] = np.where(d.branch == "whole", d.whole_pp4, d.bc_orig_pp4)
    v3 = s1[(s1.method == "V3_gate_conditional") & (s1.lead_threshold1 == 1e-3)][key + ["branch"]].merge(d[key + ["branch"]], on=key, suffixes=("_s1", ""))
    n_branch = int((v3.branch_s1 != v3.branch.replace({"no_valid_block": "block"})).sum())
    out.append(f"branch (whole / block) mismatches for V3_gate_conditional: {n_branch} of {len(v3)}")
    for meth, col in [("whole_abf", "whole_pp4"), ("V3_gate_conditional", "v3_equiv_pp4"), ("signals_mask", "sig_mask_pp4"), ("signals_cond", "sig_cond_pp4")]:
        x = s1[(s1.method == meth) & ((s1.lead_threshold1 == 1e-3) | (meth != "V3_gate_conditional"))][key + ["pp4"]]
        mg = x.merge(d[key + [col]], on=key)
        diff = (mg.pp4 - mg[col]).abs()
        out.append(f"{meth}: {len(mg)} loci, max |diff| {np.nanmax(diff.values) if len(mg) else float('nan'):.2e}, NA mismatches {int((mg.pp4.isna() != mg[col].isna()).sum())}")
    check = "; ".join(out)

matched.to_csv(os.path.join(a.out, "prior_rule_matched.tsv"), sep="\t", index=False)
trace.to_csv(os.path.join(a.out, "prior_rule_traceCB.tsv"), sep="\t", index=False)
paired.to_csv(os.path.join(a.out, "prior_rule_paired.tsv"), sep="\t", index=False)

pct = lambda v: f"{100 * v:.1f}%"
def md(df, fmt):
    cols = list(df.columns)
    L = ["| " + " | ".join(cols) + " |", "|" + "---|" * len(cols)]
    for _, r in df.iterrows(): L.append("| " + " | ".join(fmt(c, r[c]) for c in cols) + " |")
    return "\n".join(L)
f1 = lambda c, v: (f"{v:.3f}" if c == "AUC" else (str(v) if c in ("method", "cutoff") else pct(v)))
f2 = lambda c, v: (str(v) if c in ("method", "group", "n") else pct(v))
f3 = lambda c, v: (str(v) if c in ("rule", "method", "vs", "part") else (f"{v:.2g}" if c == "p" else f"{100 * v:+.1f}"))
h = d.iloc[0]
# three-line conclusion kept in conclusions_prior.md in --out (written separately), so re-running keeps it
cf = os.path.join(a.out, "conclusions_prior.md")
concl = open(cf).read().strip() if os.path.exists(cf) else "(to be written separately)"
txt = [f"# Block priors and decision rules (GWAS h2 {h.h2_gwas:g}, eQTL h2 {h.h2_eqtl:g}, fixed effects, low LD)", "",
       f"Loci: {len(d)} ({(d.set == 'test').sum()} test). Same simulated loci as stage 1 (same seeds). block-coloc uses the whole-region fallback.",
       "block-coloc, region-equivalent block priors: p1, p2, p12 multiplied by M / M_B inside a block of M_B SNPs (region of M SNPs).", "",
       "## Three-line conclusion", "", concl, "",
       "## 1. Matched 5% FPR (cut-off on calibration windows; test windows)", "", md(matched, f1), "",
       "## 2. traceCB rule (PP4 > 0.7 colocalised, PP3 > 0.7 independent; test windows)", "", md(trace, f2), "",
       "## 3. Paired differences on the same test loci (percentage points, McNemar)", "", md(paired, f3), "",
       "## 4. Prior diagnostics (loci that used a block)", "", "; ".join(f"{k} {v:.3g}" if isinstance(v, float) else f"{k} {v}" for k, v in diag.items()), ""]
if check: txt += ["## 5. Consistency with stage 1 (PP4 should be identical)", "", check, ""]
with open(os.path.join(a.out, "PRIOR_RULE.md"), "w") as fh: fh.write("\n".join(txt))
print("\n".join(txt))
