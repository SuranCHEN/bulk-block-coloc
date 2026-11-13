#!/usr/bin/env python3
"""Window-level inference for the bulk block-coloc comparisons.

Loci within a window share the LD structure and the window-specific prior-mass
effect (Proposition 7), so they are not independent replicates of "a new region".
For each comparison we report, on the same positive test loci:
  locus     McNemar exact test and paired Wald CI (as in Results; loci independent)
  win_t     mean of the 6 per-test-window differences, t interval (df = 5)
  win_boot  cluster bootstrap over test windows of the pooled difference,
            cut-offs fixed at their reported values
  two_level cluster bootstrap over calibration windows (cut-offs re-derived at 5%
            FPR in each resample) and over test windows; this also carries the
            sampling error of the cut-offs
Usage: python3 window_inference.py <out.md>
"""
import sys, math
import numpy as np, pandas as pd
from scipy.stats import binomtest, t as tdist

import os
# Run from the bulk-block-coloc repository root; override with BULK_OUT / BULK_RES.
ROOT = os.environ.get("BULK_OUT", "../clike_ldblock_coloc_handoff_20260811/outputs/").rstrip("/") + "/"
RES = os.environ.get("BULK_RES", "results/").rstrip("/") + "/"
POS = ["A_A", "A_AB", "AB_AB", "AB_AC"]; NEG = ["A_B", "AB_CD"]
METHODS = {"block-coloc": ("V3_gate_conditional", True), "coloc.abf": ("whole_abf", False),
           "signals (cond)": ("signals_cond", False), "signals (mask)": ("signals_mask", False),
           "coloc.susie + fb": ("susie_fallback", False)}
KEY = {"block-coloc": "V3_gate_conditional_fb@1e-03", "coloc.abf": "whole_abf", "signals (cond)": "signals_cond",
       "signals (mask)": "signals_mask", "coloc.susie + fb": "susie_fallback"}
B = int(os.environ.get("BULK_BOOT", "2000"))
rng = np.random.default_rng(20261113)


def load(path):
    d = pd.read_csv(path, sep="\t", low_memory=False)
    out = {}
    for name, (m, fb) in METHODS.items():
        x = d[d.method == m]
        if m.startswith("V3"):
            x = x[(x.lead_threshold1 == 1e-3) & (x.lead_threshold2 == 1e-3)]
        x = x.copy()
        pp = x.pp4.where(np.isfinite(x.pp4), x.pp4_whole) if fb else x.pp4
        x["score"] = np.nan_to_num(pp.values.astype(float), nan=-1.0)
        out[name] = x.set_index(["window", "scenario", "rep", "h2_gwas", "h2_eqtl", "tier"])["score"]
    df = pd.DataFrame(out).reset_index()
    sets = d.drop_duplicates("window").set_index("window")["set"]
    df["set"] = df.window.map(sets)
    return df


def cutoff(scores):
    v = np.sort(scores)[::-1]
    return v[int(math.floor(0.05 * len(v)))]


def mcnemar(a, b):
    n = len(a); b10 = int(np.sum(a & ~b)); b01 = int(np.sum(~a & b)); d = (b10 - b01) / n
    se = math.sqrt(max((b10 + b01) / n - d * d, 0) / n)
    p = binomtest(b10, b10 + b01, 0.5).pvalue if b10 + b01 else 1.0
    return 100 * d, 100 * (d - 1.96 * se), 100 * (d + 1.96 * se), p


def analyse(test_pos, cal_neg, thr, comps):
    """test_pos: positive test loci (all methods' scores); cal_neg: calibration negatives
    used to (re)derive cut-offs in the two-level bootstrap; thr: reported cut-offs."""
    wins = np.sort(test_pos.window.unique()); cwins = np.sort(cal_neg.window.unique())
    det = {m: test_pos[m].values > thr[m] for m in METHODS}
    # per-window arrays for bootstraps
    tw = {w: test_pos.window.values == w for w in wins}
    cw = {w: cal_neg[cal_neg.window == w] for w in cwins}
    rows = []
    for a, b in comps:
        dd, lo, hi, p = mcnemar(det[a], det[b])
        per = np.array([100 * (det[a][tw[w]].mean() - det[b][tw[w]].mean()) for w in wins])
        k = len(per); m_ = per.mean(); se = per.std(ddof=1) / math.sqrt(k); q = tdist.ppf(0.975, k - 1)
        # fixed-threshold window bootstrap
        sa = np.array([det[a][tw[w]].sum() for w in wins]); sb = np.array([det[b][tw[w]].sum() for w in wins])
        nn = np.array([tw[w].sum() for w in wins])
        bs = []
        for _ in range(B):
            i = rng.integers(0, k, k); bs.append(100 * (sa[i].sum() - sb[i].sum()) / nn[i].sum())
        # two-level bootstrap
        sc_a = {w: test_pos[a].values[tw[w]] for w in wins}; sc_b = {w: test_pos[b].values[tw[w]] for w in wins}
        ca = {w: cw[w][a].values for w in cwins}; cb = {w: cw[w][b].values for w in cwins}
        b2 = []
        for _ in range(B):
            ic = rng.choice(cwins, len(cwins)); it = rng.choice(wins, k)
            ta = cutoff(np.concatenate([ca[w] for w in ic])); tb = cutoff(np.concatenate([cb[w] for w in ic]))
            xa = np.concatenate([sc_a[w] for w in it]); xb = np.concatenate([sc_b[w] for w in it])
            b2.append(100 * ((xa > ta).mean() - (xb > tb).mean()))
        rows.append(dict(comp=f"{a} − {b}", locus=f"{dd:+.1f} ({lo:+.1f}, {hi:+.1f}), p = {p:.1g}",
                         win_t=f"{m_:+.1f} ({m_-q*se:+.1f}, {m_+q*se:+.1f})",
                         win_boot="({:+.1f}, {:+.1f})".format(*np.percentile(bs, [2.5, 97.5])),
                         two_level="({:+.1f}, {:+.1f})".format(*np.percentile(b2, [2.5, 97.5])),
                         pos=f"{int((per > 0).sum())}/{k}", per=", ".join(f"{x:+.1f}" for x in per)))
    return rows


def table(rows):
    L = ["| comparison | locus-level (McNemar) | window mean, t 95% CI | window bootstrap | two-level bootstrap | windows > 0 | per-window differences |",
         "|---|---|---|---|---|---|---|"]
    for r in rows:
        L.append(f"| {r['comp']} | {r['locus']} | {r['win_t']} | {r['win_boot']} | {r['two_level']} | {r['pos']} | {r['per']} |")
    return "\n".join(L)


def power_line(test_pos, thr):
    return ", ".join(f"{m} {100*np.mean(test_pos[m].values > thr[m]):.1f}" for m in METHODS)


def thresholds(path):
    th = pd.read_csv(path, sep="\t").set_index("key")["cutoff"]
    return {m: th[KEY[m]] for m in METHODS}


COMPS = [("block-coloc", "coloc.abf"), ("block-coloc", "signals (cond)"), ("block-coloc", "signals (mask)"),
         ("block-coloc", "coloc.susie + fb")]
out = ["# Window-level inference (draft)", "",
       "Detection-rate differences in percentage points on the same positive test loci. Test windows: 6. "
       "Locus-level intervals treat loci as independent; window-level ones treat windows as the unit; the two-level bootstrap "
       f"also re-derives the 5% FPR cut-offs from resampled calibration windows ({B} resamples each).", ""]

main = load(ROOT + "bulk_h2ref_v1/scores_main.tsv.gz")
thr1 = thresholds(RES + "h2ref_stage1/thresholds.tsv")
calneg1 = main[(main.set == "calibration") & main.scenario.isin(NEG)]
tp = main[(main.set == "test") & main.scenario.isin(POS)]
out += ["## Stage 1 (default setting)", "", "Reproduced detection rates (%): " + power_line(tp, thr1), "",
        table(analyse(tp, calneg1, thr1, COMPS)), ""]

grid = load(ROOT + "bulk_h2ref_v1/scores_grid.tsv.gz")
out += ["## Stage 2 (cut-offs from stage 1; two-level bootstrap re-derives them from stage-1 calibration windows)", ""]
cells = [(0.003, 0.0075, "low"), (0.003, 0.03, "low"), (0.003, 0.05, "low"), (0.001, 0.015, "low"),
         (0.01, 0.015, "low"), (0.02, 0.015, "low"), (0.003, 0.015, "mid"), (0.003, 0.015, "high")]
for g, e, tr in cells:
    x = grid[(grid.h2_gwas == g) & (grid.h2_eqtl == e) & (grid.tier == tr)]
    tpx = x[(x.set == "test") & x.scenario.isin(POS)]
    out += [f"### h2_gwas = {g}, h2_eqtl = {e}, LD tier {tr} (n positive test loci = {len(tpx)})", "",
            "Detection rates (%): " + power_line(tpx, thr1), "",
            table(analyse(tpx, calneg1, thr1, COMPS[:2])), ""]

for lab, f, t in [("Normal effects, default strength", "bulk_h2ref_v1_normal/scores_main.tsv.gz", "h2ref_normal_default/thresholds.tsv"),
                  ("Normal effects, h2_eqtl = 5%", "bulk_h2ref_v1_normal_g0.003_e0.05/scores_main.tsv.gz", "h2ref_normal_e005/thresholds.tsv")]:
    d = load(ROOT + f); th = thresholds(RES + t)
    cn = d[(d.set == "calibration") & d.scenario.isin(NEG)]; tpx = d[(d.set == "test") & d.scenario.isin(POS)]
    out += [f"## {lab}", "", "Detection rates (%): " + power_line(tpx, th), "",
            table(analyse(tpx, cn, th, COMPS[:3])), ""]

open(sys.argv[1], "w").write("\n".join(out) + "\n")
print("\n".join(out))
