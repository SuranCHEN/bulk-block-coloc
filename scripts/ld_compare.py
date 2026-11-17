#!/usr/bin/env python3
"""Matched (in-sample) LD vs 1000 Genomes reference LD, on the same simulated loci.

Both runs use the same windows, SNPs, causal variants, phenotypes and summary
statistics (run.R --ld ref / --ld insample); only the LD matrix given to the
methods differs. coloc.abf uses no LD, so its PP4 must be identical in the two
runs: this is checked first and the script stops if it fails.

Usage (from the bulk-block-coloc root or the hand-off root):
  python3 ld_compare.py <scores_ref.tsv.gz> <scores_insample.tsv.gz> [--h2_gwas X] [--h2_eqtl X] [--out FILE.md]
<scores_ref> may be scores_main.tsv.gz or scores_grid.tsv.gz; the cell is selected by --h2_gwas / --h2_eqtl
(default 0.003 / 0.015), LD tier low. Replicates present in both files are used.

Reported for every method available in both runs:
  A  detection and false-positive rate on test windows at a matched 5% FPR, cut-offs set on the calibration
     windows separately for each LD (each method at its own 5%)
  B  the same with the fixed rule PP4 > 0.7 (no calibration)
  C  paired change, reference minus in-sample, on the same loci (locus-level McNemar; window-level t interval)
  D  difference in differences: (block-coloc - method) under reference LD minus the same under in-sample LD
  E  lead counts (conditional selection) against the true number of causal variants, and SuSiE no-call rates
"""
import sys, math
import numpy as np, pandas as pd
from scipy.stats import binomtest, t as tdist

POS = ["A_A", "A_AB", "AB_AB", "AB_AC"]; NEG = ["A_B", "AB_CD"]
LAB = {"A_A": "single_shared", "A_AB": "single_double_shared", "AB_AB": "double_shared2", "AB_AC": "double_shared1",
       "A_B": "single_distinct", "AB_CD": "double_distinct", "gwas_only": "gwas_only", "null_both": "null_both"}
ID = ["window", "scenario", "rep"]
# name -> (method, lead threshold or None, fallback to whole-region PP4 when no valid block)
KEYS = {
    "coloc.abf": ("whole_abf", None, False),
    "block-coloc (conditional leads, 1e-3)": ("V3_gate_conditional", 1e-3, True),
    "block-coloc (marginal leads, 1e-3)": ("V2_gate_marginal", 1e-3, True),
    "block-coloc (marginal leads, 1e-5)": ("V2_gate_marginal", 1e-5, True),
    "block, no gate (marginal leads, 1e-3)": ("V1_block_any", 1e-3, True),
    "block-coloc (original handoff, 1e-4)": ("V0_handoff", 1e-4, True),
    "coloc.susie": ("susie_raw", None, False),
    "coloc.susie + fallback": ("susie_fallback", None, False),
    "coloc.signals (mask)": ("signals_mask", None, False),
    "coloc.signals (cond)": ("signals_cond", None, False),
}
MAINKEY = "block-coloc (conditional leads, 1e-3)"


def opt(argv, name, default):
    return argv[argv.index(name) + 1] if name in argv else default


def load(path, hg, he):
    d = pd.read_csv(path, sep="\t", low_memory=False)
    d = d[np.isclose(d.h2_gwas, hg) & np.isclose(d.h2_eqtl, he) & (d.tier == "low")]
    return d


def wide(d):
    cols = {}
    for name, (m, t, fb) in KEYS.items():
        x = d[d.method == m]
        if t is not None:
            x = x[np.isclose(x.lead_threshold1, t) & np.isclose(x.lead_threshold2, t)]
        if not len(x):
            continue
        pp = x.pp4.where(np.isfinite(x.pp4), x.pp4_whole) if fb else x.pp4
        cols[name] = pd.Series(pp.values.astype(float), index=pd.MultiIndex.from_frame(x[ID]))
    w = pd.DataFrame(cols)
    return w


def cutoff(v):
    v = np.sort(np.nan_to_num(v, nan=-1.0))[::-1]
    return v[int(math.floor(0.05 * len(v)))]


def win_t(a, b, win):
    per = np.array([100 * (a[win == w].mean() - b[win == w].mean()) for w in np.unique(win)])
    k = len(per); m = per.mean()
    if k < 2:
        return m, np.nan, np.nan, per
    q = tdist.ppf(0.975, k - 1); se = per.std(ddof=1) / math.sqrt(k)
    return m, m - q * se, m + q * se, per


def main(argv):
    ref_path, ins_path = [a for a in argv if not a.startswith("--") and a not in
                          (opt(argv, "--h2_gwas", None), opt(argv, "--h2_eqtl", None), opt(argv, "--out", None))][:2]
    hg = float(opt(argv, "--h2_gwas", "0.003")); he = float(opt(argv, "--h2_eqtl", "0.015"))
    dr, di = load(ref_path, hg, he), load(ins_path, hg, he)
    if not len(dr) or not len(di):
        sys.exit(f"cell h2_gwas={hg}, h2_eqtl={he} not found in both files")
    wr, wi = wide(dr), wide(di)
    names = [n for n in KEYS if n in wr.columns and n in wi.columns]
    idx = wr.index.intersection(wi.index)
    wr, wi = wr.loc[idx, names], wi.loc[idx, names]
    meta = pd.DataFrame(index=idx).reset_index()
    sets = dr.drop_duplicates("window").set_index("window")["set"]
    meta["set"] = meta.window.map(sets).values
    L = [f"# Matched (in-sample) LD vs 1000 Genomes reference LD: h2_gwas = {hg}, h2_eqtl = {he}", "",
         f"Loci in both runs: {len(idx)} ({meta.rep.nunique()} replicates x {meta.window.nunique()} windows x {meta.scenario.nunique()} scenarios). "
         f"Reference: `{ref_path}`; in-sample: `{ins_path}`.", ""]
    # ---- identity check
    a, b = wr["coloc.abf"].values, wi["coloc.abf"].values
    md = float(np.nanmax(np.abs(a - b))); same_na = bool(np.all(np.isnan(a) == np.isnan(b)))
    L += [f"**Check.** Whole-region coloc.abf PP4, reference vs in-sample: max |diff| = {md:.3g}; NA pattern identical: {same_na}.", ""]
    if md > 1e-8 or not same_na:
        L += ["STOP: the two runs are not on the same loci (coloc.abf uses no LD and must be identical)."]
        print("\n".join(L)); sys.exit(1)
    cal = (meta.set == "calibration").values; test = (meta.set == "test").values
    sc = meta.scenario.values; win = meta.window.values
    pos = test & np.isin(sc, POS); neg = test & np.isin(sc, NEG); cneg = cal & np.isin(sc, NEG); go = test & (sc == "gwas_only")
    S = {"ref": wr, "ins": wi}
    thr = {k: {n: cutoff(S[k][n].values[cneg]) for n in names} for k in S}
    call = {k: {n: np.nan_to_num(S[k][n].values, nan=-1.0) > thr[k][n] for n in names} for k in S}
    fix = {k: {n: np.nan_to_num(S[k][n].values, nan=-1.0) > 0.7 for n in names} for k in S}
    pct = lambda v: f"{100 * np.mean(v):.1f}" if len(v) else "NA"
    # ---- A
    L += ["## A. Matched 5% FPR (cut-offs from the calibration windows of each LD), test windows", "",
          "Detection rate / false-positive rate (single_distinct + double_distinct) / gwas_only, %; cut-off in brackets.", "",
          "| method | in-sample LD | reference LD | detection: reference − in-sample |", "|---|---|---|---|"]
    for n in names:
        ci, cr = call["ins"][n], call["ref"][n]
        L.append(f"| {n} | {pct(ci[pos])} / {pct(ci[neg])} / {pct(ci[go])} [{thr['ins'][n]:.3f}] | "
                 f"{pct(cr[pos])} / {pct(cr[neg])} / {pct(cr[go])} [{thr['ref'][n]:.3f}] | {100 * (cr[pos].mean() - ci[pos].mean()):+.1f} |")
    L += ["", "By scenario (detection or false-positive rate, in-sample → reference):", "",
          "| method | " + " | ".join(LAB[s] for s in POS + NEG) + " |", "|---|" + "---|" * (len(POS) + len(NEG))]
    for n in names:
        L.append(f"| {n} | " + " | ".join(f"{pct(call['ins'][n][test & (sc == s)])} → {pct(call['ref'][n][test & (sc == s)])}" for s in POS + NEG) + " |")
    # ---- B
    L += ["", "## B. Fixed rule PP4 > 0.7, test windows", "",
          "| method | colocalised loci called, in-sample → reference | not-colocalised loci called, in-sample → reference |", "|---|---|---|"]
    for n in names:
        L.append(f"| {n} | {pct(fix['ins'][n][pos])} → {pct(fix['ref'][n][pos])} | {pct(fix['ins'][n][neg])} → {pct(fix['ref'][n][neg])} |")
    # ---- C
    L += ["", "## C. Paired change in detection, reference − in-sample (same positive test loci, matched 5% FPR)", "",
          "| method | locus-level diff, McNemar p | window mean (t 95% CI) | windows < 0 | per-window |", "|---|---|---|---|---|"]
    for n in names:
        cr, ci = call["ref"][n][pos], call["ins"][n][pos]
        b10 = int((cr & ~ci).sum()); b01 = int((~cr & ci).sum())
        p = binomtest(b10, b10 + b01, 0.5).pvalue if b10 + b01 else 1.0
        m, lo, hi, per = win_t(cr, ci, win[pos])
        L.append(f"| {n} | {100 * (b10 - b01) / len(cr):+.1f}, p = {p:.2g} | {m:+.1f} ({lo:+.1f}, {hi:+.1f}) | {int((per < 0).sum())}/{len(per)} | "
                 + ", ".join(f"{x:+.1f}" for x in per) + " |")
    # ---- D
    L += ["", f"## D. Difference in differences: ({MAINKEY} − method) under reference LD minus the same under in-sample LD", "",
          "Positive: the method loses more than block-coloc when the matched LD is replaced by the reference LD.", "",
          "| method | block-coloc − method, in-sample | block-coloc − method, reference | difference (window t 95% CI) | windows > 0 |", "|---|---|---|---|---|"]
    if MAINKEY in names:
        for n in names:
            if n in (MAINKEY, "coloc.abf"):
                continue
            w = win[pos]; ws = np.unique(w)
            dd = np.array([100 * ((call["ref"][MAINKEY][pos][w == x].mean() - call["ref"][n][pos][w == x].mean())
                                  - (call["ins"][MAINKEY][pos][w == x].mean() - call["ins"][n][pos][w == x].mean())) for x in ws])
            k = len(dd); m = dd.mean(); q = tdist.ppf(0.975, k - 1) if k > 1 else np.nan; se = dd.std(ddof=1) / math.sqrt(k) if k > 1 else np.nan
            di_ = 100 * (call["ins"][MAINKEY][pos].mean() - call["ins"][n][pos].mean()); dr_ = 100 * (call["ref"][MAINKEY][pos].mean() - call["ref"][n][pos].mean())
            L.append(f"| {n} | {di_:+.1f} | {dr_:+.1f} | {m:+.1f} ({m - q * se:+.1f}, {m + q * se:+.1f}) | {int((dd > 0).sum())}/{k} |")
    # ---- E
    L += ["", "## E. Signal counts and no-call rates (all windows)", ""]
    def v3(d):
        x = d[(d.method == "V3_gate_conditional") & np.isclose(d.lead_threshold1, 1e-3) & np.isclose(d.lead_threshold2, 1e-3)]
        return x.set_index(ID).loc[idx]
    xr, xi = v3(dr), v3(di)
    L += ["Conditional lead selection (lead threshold 1e-3): share of loci where the number of leads equals the number of causal variants, in-sample → reference; and share taking the block branch.", "",
          "| scenario | trait 1: leads = causal | trait 2: leads = causal | block branch | no valid block |", "|---|---|---|---|---|"]
    for s in POS + NEG + ["gwas_only", "null_both"]:
        a_, b_ = xi[xi.index.get_level_values("scenario") == s], xr[xr.index.get_level_values("scenario") == s]
        if not len(a_):
            continue
        f = lambda x, c, k: pct((x[c] == x[k]).values)
        L.append(f"| {LAB[s]} | {f(a_, 'n_leads1', 'n_causal_gwas')} → {f(b_, 'n_leads1', 'n_causal_gwas')} | {f(a_, 'n_leads2', 'n_causal_eqtl')} → {f(b_, 'n_leads2', 'n_causal_eqtl')} | "
                 f"{pct((a_.branch == 'block').values)} → {pct((b_.branch == 'block').values)} | {pct((a_.status == 'no_valid_block').values)} → {pct((b_.status == 'no_valid_block').values)} |")
    L += ["", "No result (PP4 = NA), positives on test windows, in-sample → reference: " +
          "; ".join(f"{n} {pct(np.isnan(wi[n].values[pos]))} → {pct(np.isnan(wr[n].values[pos]))}" for n in names if n.startswith("coloc.s")), ""]
    txt = "\n".join(L) + "\n"
    out = opt(argv, "--out", None)
    if out:
        open(out, "w").write(txt)
    print(txt)


if __name__ == "__main__":
    main(sys.argv[1:])
