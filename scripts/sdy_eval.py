#!/usr/bin/env python3
"""Evaluate the sdY sensitivity run (block_sdy.R).

Usage: python3 sdy_eval.py <sdy_loci.tsv.gz> [<scores file with the same loci>] [--out FILE.md]
The optional scores file (scores_main.tsv.gz or scores_grid.tsv.gz) is used for the reproduction check:
arm 'est' must reproduce V3_gate_conditional (lead threshold 1e-3) pp4 / pp4_whole exactly.

Per arm (est / whole / one): block-coloc (adopted rule) and whole-region coloc.abf, each with its own cut-off
at 5% FPR on calibration windows (single_distinct + double_distinct); test-window detection, FPR, gwas_only;
paired block-coloc minus coloc.abf (McNemar; window t interval); fixed rule PP4 > 0.7 / PP3 > 0.7.
Diagnostics: sdY estimates (whole region vs best block), |PP4 change| between arms, best-block changes.
"""
import sys, math
import numpy as np, pandas as pd
from scipy.stats import binomtest, t as tdist

POS = ["A_A", "A_AB", "AB_AB", "AB_AC"]; NEG = ["A_B", "AB_CD"]
LAB = {"A_A": "single_shared", "A_AB": "single_double_shared", "AB_AB": "double_shared2", "AB_AC": "double_shared1",
       "A_B": "single_distinct", "AB_CD": "double_distinct"}
ID = ["window", "scenario", "rep"]
ARMS = ["est", "whole", "one"]
DESC = {"est": "sdY estimated on each subset (current)", "whole": "sdY estimated once on the whole region",
        "one": "sdY = 1 (true value)"}


def cutoff(v):
    v = np.sort(np.nan_to_num(np.asarray(v, float), nan=-1.0))[::-1]
    return v[int(math.floor(0.05 * len(v)))]


def paired(ca, cb, win):
    b10 = int((ca & ~cb).sum()); b01 = int((~ca & cb).sum())
    p = binomtest(b10, b10 + b01, 0.5).pvalue if b10 + b01 else 1.0
    per = np.array([100 * (ca[win == w].mean() - cb[win == w].mean()) for w in np.unique(win)])
    q = tdist.ppf(0.975, len(per) - 1); se = per.std(ddof=1) / math.sqrt(len(per))
    return (f"{100 * (b10 - b01) / len(ca):+.1f}, p = {p:.2g}",
            f"{per.mean():+.1f} ({per.mean() - q * se:+.1f}, {per.mean() + q * se:+.1f}), {int((per > 0).sum())}/{len(per)}")


def qs(x, ps=(0.5, 0.9, 0.99)):
    x = np.asarray(x, float); x = x[np.isfinite(x)]
    return " / ".join(f"{np.quantile(x, p):.4g}" for p in ps) if len(x) else "—"


def main(argv):
    out = argv[argv.index("--out") + 1] if "--out" in argv else None
    args = [a for a in argv if not a.startswith("--") and a != out]
    d = pd.read_csv(args[0], sep="\t", low_memory=False)
    d["score"] = np.where(d.status == "block", d.pp4, d.pp4_whole)
    d["score3"] = np.where(d.status == "block", d.pp3, d.pp3_whole)
    hg, he = d.h2_gwas.iloc[0], d.h2_eqtl.iloc[0]
    L = [f"# sdY sensitivity: h2_gwas = {hg}, h2_eqtl = {he}", "",
         f"{d[ID].drop_duplicates().shape[0]} loci, arms: " + "; ".join(f"`{a}` = {DESC[a]}" for a in ARMS) + ".",
         "Leads and blocks are the same in every arm; each method and arm gets its own cut-off at 5% FPR on the calibration windows.", ""]
    if len(args) > 1:
        s = pd.read_csv(args[1], sep="\t", low_memory=False)
        s = s[(s.method == "V3_gate_conditional") & np.isclose(s.lead_threshold1, 1e-3) & np.isclose(s.lead_threshold2, 1e-3)
              & np.isclose(s.h2_gwas, hg) & np.isclose(s.h2_eqtl, he) & (s.tier == "low")]
        b = d[d.arm == "est"].merge(s[ID + ["pp4", "pp4_whole"]], on=ID, suffixes=("", "_s1"))
        fin = np.isfinite(b.pp4) & np.isfinite(b.pp4_s1); nan_same = bool(np.all(np.isfinite(b.pp4) == np.isfinite(b.pp4_s1)))
        md = float(np.max(np.abs(b.pp4[fin] - b.pp4_s1[fin]))) if fin.any() else float("nan")
        mw = float(np.max(np.abs(b.pp4_whole - b.pp4_whole_s1)))
        L += [f"**Reproduction check** (arm est vs `{args[1]}`): {len(b)} loci matched; block PP4 max |diff| = {md:.3g}; "
              f"NA pattern identical: {nan_same}; whole-region PP4 max |diff| = {mw:.3g}.", ""]
    piv = {a: d[d.arm == a].set_index(ID).sort_index() for a in ARMS}
    idx = piv["est"].index
    for a in ARMS:
        piv[a] = piv[a].reindex(idx)
    meta = piv["est"]
    wdiff = float(np.nanmax(np.abs(piv["whole"].pp4_whole.values - piv["est"].pp4_whole.values)))
    L += [f"**Internal check**: whole-region PP4, arm whole vs arm est, max |diff| = {wdiff:.3g} (should be ~0: same sdY); "
          f"own sdY estimate vs coloc:::sdY.est, max |diff| = {np.nanmax(meta.sdy_check.values):.3g}; "
          f"loci with NA sdy_check: {int(np.isnan(meta.sdy_check.values).sum())}.", ""]
    st = meta.set.values; sc = idx.get_level_values("scenario"); win = idx.get_level_values("window").values
    cal = (st == "calibration") & np.isin(sc, NEG); pos = (st == "test") & np.isin(sc, POS)
    neg = (st == "test") & np.isin(sc, NEG); go = (st == "test") & (sc == "gwas_only")
    L += ["## Matched 5% FPR, test windows", "",
          "| arm | method | cut-off | detection | FPR (distinct) | gwas_only | block-coloc − coloc.abf (pp), McNemar p | window mean (t 95% CI), windows > 0 |",
          "|---|---|---|---|---|---|---|---|"]
    calls = {}
    for a in ARMS:
        x = piv[a]
        for m, col in (("block-coloc", "score"), ("coloc.abf", "pp4_whole")):
            v = x[col].values; t = cutoff(v[cal]); c = np.nan_to_num(v, nan=-1) > t; calls[(a, m)] = c
            cmp, wcmp = paired(calls[(a, "block-coloc")][pos], c[pos], win[pos]) if m == "coloc.abf" else ("", "")
            L.append(f"| {a} | {m} | {t:.4f} | {100 * c[pos].mean():.1f} | {100 * c[neg].mean():.1f} | {100 * c[go].mean():.1f} | {cmp} | {wcmp} |")
    L += ["", "Block-coloc, arm vs arm est (same positive test loci):", "", "| arm | diff vs est (pp), McNemar p | window mean (t 95% CI) |", "|---|---|---|"]
    for a in ["whole", "one"]:
        cmp, wcmp = paired(calls[(a, "block-coloc")][pos], calls[("est", "block-coloc")][pos], win[pos])
        L.append(f"| {a} | {cmp} | {wcmp} |")
    L += ["", "## By scenario (block-coloc / coloc.abf, %)", "",
          "| arm | " + " | ".join(LAB[s] for s in POS + NEG) + " |", "|---|" + "---|" * (len(POS) + len(NEG))]
    for a in ARMS:
        L.append(f"| {a} | " + " | ".join(
            f"{100 * calls[(a, 'block-coloc')][(st == 'test') & (sc == s)].mean():.1f} / {100 * calls[(a, 'coloc.abf')][(st == 'test') & (sc == s)].mean():.1f}"
            for s in POS + NEG) + " |")
    L += ["", "## Fixed rule (PP4 > 0.7 colocalised; PP3 > 0.7 independent), test windows, block-coloc", "",
          "| arm | colocalised: called coloc / indep | single_distinct: indep (coloc) | double_distinct: indep (coloc) | gwas_only: coloc |", "|---|---|---|---|---|"]
    for a in ARMS:
        x = piv[a]; p4 = np.nan_to_num(x.score.values, nan=0); p3 = np.nan_to_num(x.score3.values, nan=0)
        f = lambda m: (100 * np.mean(p4[m] > 0.7), 100 * np.mean(p3[m] > 0.7))
        A = f(pos); s1 = f((st == "test") & (sc == "A_B")); s2 = f((st == "test") & (sc == "AB_CD")); g = f(go)
        L.append(f"| {a} | {A[0]:.1f} / {A[1]:.1f} | {s1[1]:.1f} ({s1[0]:.1f}) | {s2[1]:.1f} ({s2[0]:.1f}) | {g[0]:.1f} |")
    blk = meta.status.values == "block"
    L += ["", "## Diagnostics (all windows)", "",
          "| quantity | median / p90 / p99 |", "|---|---|",
          f"| sdY trait 1, whole region | {qs(meta.sdy1_whole)} |",
          f"| sdY trait 2, whole region | {qs(meta.sdy2_whole)} |",
          f"| sdY trait 1, best block / whole region (block branch, arm est) | {qs(meta.sdy1_best_est[blk] / meta.sdy1_whole[blk])} |",
          f"| sdY trait 2, best block / whole region (block branch, arm est) | {qs(meta.sdy2_best_est[blk] / meta.sdy2_whole[blk])} |",
          f"| abs(log sdY2 ratio), block branch | {qs(np.abs(np.log(meta.sdy2_best_est[blk] / meta.sdy2_whole[blk])))} |"]
    for a in ["whole", "one"]:
        L.append(f"| block-coloc abs(PP4 {a} − est), all loci | {qs(np.abs(piv[a].score.values - meta.score.values))} |")
        L.append(f"| coloc.abf abs(PP4 {a} − est), all loci | {qs(np.abs(piv[a].pp4_whole.values - meta.pp4_whole.values))} |")
        sb = piv[a].best_same_as_est.astype(str).str.upper().values[blk]
        L.append(f"| best block differs from arm est ({a}), block branch | {100 * np.mean(sb == 'FALSE'):.1f}% |")
    L.append(f"| block subsets where the sdY estimate is NA (negative coefficient), block branch | "
             f"{100 * np.mean(~np.isfinite(meta.sdy2_best_est.values[blk]) | ~np.isfinite(meta.sdy1_best_est.values[blk])):.1f}% |")
    txt = "\n".join(L) + "\n"
    if out:
        open(out, "w").write(txt)
    print(txt)


if __name__ == "__main__":
    main(sys.argv[1:])
