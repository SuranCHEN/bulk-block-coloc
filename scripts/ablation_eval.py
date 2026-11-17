#!/usr/bin/env python3
"""Evaluate the block-coloc construction ablation (block_ablation.R).

Usage: python3 ablation_eval.py <ablation_loci.tsv.gz> [<scores file with the same loci>] [--out FILE.md]
The optional scores file (scores_main.tsv.gz or scores_grid.tsv.gz) is used for the reproduction check:
config 'base' must reproduce V3_gate_conditional (lead threshold 1e-3) pp4 / pp4_whole exactly.

For every configuration: adopted rule (PP4 of the best block; whole-region PP4 when the gate is not passed
or no valid block), cut-off at 5% FPR on calibration windows (single_distinct + double_distinct), detection
and false-positive rates on test windows, by scenario; share of positives taking the block branch;
localisation (best block contains the shared causal variant); fixed rule PP4 > 0.7 / PP3 > 0.7; paired
comparison with 'base' on the same positive test loci (locus level: McNemar; window level: t interval over
the 6 test windows).
"""
import sys, math
import numpy as np, pandas as pd
from scipy.stats import binomtest, t as tdist

POS = ["A_A", "A_AB", "AB_AB", "AB_AC"]; NEG = ["A_B", "AB_CD"]
LAB = {"A_A": "single_shared", "A_AB": "single_double_shared", "AB_AB": "double_shared2", "AB_AC": "double_shared1",
       "A_B": "single_distinct", "AB_CD": "double_distinct", "gwas_only": "gwas_only", "null_both": "null_both"}
ID = ["window", "scenario", "rep"]
DESC = {"base": "stage-1 rule (clump 0.1, neighbourhood 0.4, sign split, gate)", "nosplit": "no sign split",
        "rb0.2": "neighbourhood r2 >= 0.2", "rb0.3": "neighbourhood r2 >= 0.3", "rb0.5": "neighbourhood r2 >= 0.5",
        "rb0.6": "neighbourhood r2 >= 0.6", "rc0.05": "clump r2 < 0.05", "rc0.2": "clump r2 < 0.2",
        "nogate": "no gate (block whenever both traits have a lead)"}


def cutoff(v):
    v = np.sort(np.nan_to_num(np.asarray(v, float), nan=-1.0))[::-1]
    return v[int(math.floor(0.05 * len(v)))]


def tobool(s):
    return s.astype(str).str.upper().map({"TRUE": True, "FALSE": False})


def main(argv):
    out = argv[argv.index("--out") + 1] if "--out" in argv else None
    args = [a for a in argv if not a.startswith("--") and a != out]
    d = pd.read_csv(args[0], sep="\t", low_memory=False)
    d["score"] = np.where(d.status == "block", d.pp4, d.pp4_whole)
    d["score3"] = np.where(d.status == "block", d.pp3, d.pp3_whole)
    cfgs = [c for c in DESC if c in set(d.config)] + sorted(set(d.config) - set(DESC))
    hg, he = d.h2_gwas.iloc[0], d.h2_eqtl.iloc[0]
    L = [f"# block-coloc construction ablation: h2_gwas = {hg}, h2_eqtl = {he}", "",
         f"{d[ID].drop_duplicates().shape[0]} loci, {len(cfgs)} configurations. Same loci for every configuration; "
         "each configuration gets its own cut-off at 5% FPR on the calibration windows.", ""]
    # reproduction check
    if len(args) > 1:
        s = pd.read_csv(args[1], sep="\t", low_memory=False)
        s = s[(s.method == "V3_gate_conditional") & np.isclose(s.lead_threshold1, 1e-3) & np.isclose(s.lead_threshold2, 1e-3)
              & np.isclose(s.h2_gwas, hg) & np.isclose(s.h2_eqtl, he) & (s.tier == "low")]
        b = d[d.config == "base"].merge(s[ID + ["pp4", "pp4_whole"]], on=ID, suffixes=("", "_s1"))
        fin = np.isfinite(b.pp4) & np.isfinite(b.pp4_s1); nan_same = bool(np.all(np.isfinite(b.pp4) == np.isfinite(b.pp4_s1)))
        md = float(np.max(np.abs(b.pp4[fin] - b.pp4_s1[fin]))) if fin.any() else float("nan")
        mw = float(np.max(np.abs(b.pp4_whole - b.pp4_whole_s1)))
        L += [f"**Reproduction check** (config base vs `{args[1]}`): {len(b)} loci matched; block PP4 max |diff| = {md:.3g}; "
              f"NA pattern identical: {nan_same}; whole-region PP4 max |diff| = {mw:.3g}.", ""]
    piv = {c: d[d.config == c].set_index(ID).sort_index() for c in cfgs}
    idx = piv["base"].index
    for c in cfgs:
        piv[c] = piv[c].reindex(idx)
    meta = piv["base"]
    st = meta.set.values; sc = idx.get_level_values("scenario"); win = idx.get_level_values("window")
    cal = (st == "calibration") & np.isin(sc, NEG); pos = (st == "test") & np.isin(sc, POS)
    neg = (st == "test") & np.isin(sc, NEG); go = (st == "test") & (sc == "gwas_only")
    thr = {c: cutoff(piv[c].score.values[cal]) for c in cfgs}
    call = {c: np.nan_to_num(piv[c].score.values, nan=-1) > thr[c] for c in cfgs}
    L += ["## Matched 5% FPR, test windows", "",
          "| configuration | cut-off | detection | FPR (distinct) | gwas_only | block branch (positives) | best block has shared variant | diff vs base (pp), McNemar p | window mean (t 95% CI), windows > 0 |",
          "|---|---|---|---|---|---|---|---|---|"]
    for c in cfgs:
        x = piv[c]; cb, cc = call["base"][pos], call[c][pos]
        blk = (x.status.values == "block") & pos
        hs = tobool(x.best_has_shared)[blk]
        if c == "base":
            cmp = "—"; wcmp = "—"
        else:
            b10 = int((cc & ~cb).sum()); b01 = int((~cc & cb).sum())
            p = binomtest(b10, b10 + b01, 0.5).pvalue if b10 + b01 else 1.0
            per = np.array([100 * (cc[win[pos] == w].mean() - cb[win[pos] == w].mean()) for w in np.unique(win[pos])])
            q = tdist.ppf(0.975, len(per) - 1); se = per.std(ddof=1) / math.sqrt(len(per))
            cmp = f"{100 * (b10 - b01) / len(cc):+.1f}, p = {p:.2g}"
            wcmp = f"{per.mean():+.1f} ({per.mean() - q * se:+.1f}, {per.mean() + q * se:+.1f}), {int((per > 0).sum())}/{len(per)}"
        L.append(f"| {c}: {DESC.get(c, c)} | {thr[c]:.3f} | {100 * cc.mean():.1f} | {100 * call[c][neg].mean():.1f} | {100 * call[c][go].mean():.1f} | "
                 f"{100 * blk.sum() / pos.sum():.1f}% | {100 * hs.mean():.1f}% | {cmp} | {wcmp} |")
    L += ["", "## By scenario (detection or false-positive rate, %)", "",
          "| configuration | " + " | ".join(LAB[s] for s in POS + NEG) + " |", "|---|" + "---|" * (len(POS) + len(NEG))]
    for c in cfgs:
        L.append(f"| {c} | " + " | ".join(f"{100 * call[c][(st == 'test') & (sc == s)].mean():.1f}" for s in POS + NEG) + " |")
    L += ["", "## Fixed rule (PP4 > 0.7 colocalised; PP3 > 0.7 independent), test windows", "",
          "| configuration | colocalised: called coloc / indep | single_distinct: indep (coloc) | double_distinct: indep (coloc) |", "|---|---|---|---|"]
    for c in cfgs:
        x = piv[c]; p4 = np.nan_to_num(x.score.values, nan=0); p3 = np.nan_to_num(x.score3.values, nan=0)
        f = lambda m: (100 * np.mean(p4[m] > 0.7), 100 * np.mean(p3[m] > 0.7))
        a = f(pos); s1 = f((st == "test") & (sc == "A_B")); s2 = f((st == "test") & (sc == "AB_CD"))
        L.append(f"| {c} | {a[0]:.1f} / {a[1]:.1f} | {s1[1]:.1f} ({s1[0]:.1f}) | {s2[1]:.1f} ({s2[0]:.1f}) |")
    L += ["", "## Leads and blocks (all windows)", "",
          "| configuration | single-signal loci with >= 2 leads in a trait | two-causal trait with >= 2 leads (trait 1 / trait 2) | median blocks (block branch) | median best-block size |",
          "|---|---|---|---|---|"]
    one = np.isin(sc, ["A_A", "A_B"])
    for c in cfgs:
        x = piv[c]
        m1 = np.isin(sc, ["AB_AB", "AB_AC", "AB_CD"]); m2 = np.isin(sc, ["A_AB", "AB_AB", "AB_AC", "AB_CD"])
        blk = x.status.values == "block"
        L.append(f"| {c} | {100 * np.mean((x.n_leads1.values[one] >= 2) | (x.n_leads2.values[one] >= 2)):.1f}% | "
                 f"{100 * np.mean(x.n_leads1.values[m1] >= 2):.1f}% / {100 * np.mean(x.n_leads2.values[m2] >= 2):.1f}% | "
                 f"{np.median(x.n_blocks.values[blk]) if blk.any() else float('nan'):.0f} | {np.median(x.best_size.values[blk]) if blk.any() else float('nan'):.0f} |")
    txt = "\n".join(L) + "\n"
    if out:
        open(out, "w").write(txt)
    print(txt)


if __name__ == "__main__":
    main(sys.argv[1:])
