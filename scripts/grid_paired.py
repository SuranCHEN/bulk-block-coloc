#!/usr/bin/env python3
"""Stage 2 (grid): detection / false-positive rates with the frozen stage-1
cut-offs, paired comparisons (McNemar exact test), and lead counts used in
doc/bulk/04_results.md section 7.

The default cell (h2_gwas 0.003, h2_eqtl 0.015, low LD) is not re-simulated in
stage 2; as in evaluate.R it is taken from stage 1 (replicates 1-50).
"block-coloc + fallback" = block-coloc (conditional lead selection) whose NA
result (no valid block) is replaced by the whole-region PP4, with its own
stage-1 cut-off (thresholds.tsv key V3_gate_conditional_fb@1e-03).

Run from the project root:  python3 scripts/bulk_h2ref_v1/grid_paired.py
Writes outputs/bulk_h2ref_v1/grid_paired.tsv and grid_leads.tsv.
Needs numpy, pandas, scipy.
"""
import numpy as np
import pandas as pd
from scipy.stats import binomtest

OUT = "outputs/bulk_h2ref_v1"
KEY = {"whole_abf": "whole_abf", "susie_fallback": "susie_fallback", "signals_mask": "signals_mask",
       "signals_cond": "signals_cond", "V3_gate_conditional": "V3_gate_conditional@1e-03"}
POS, NEG = ["A_A", "A_AB", "AB_AB", "AB_AC"], ["A_B", "AB_CD"]
CELLS = [(0.003, 0.0075, "low"), (0.003, 0.015, "low"), (0.003, 0.03, "low"), (0.003, 0.05, "low"),
         (0.001, 0.015, "low"), (0.01, 0.015, "low"), (0.02, 0.015, "low"),
         (0.003, 0.015, "mid"), (0.003, 0.015, "high")]
PAIRS = [("block_coloc", "whole_abf"), ("block_coloc_fb", "whole_abf"), ("block_coloc", "signals_cond"),
         ("block_coloc_fb", "signals_cond"), ("block_coloc_fb", "signals_mask")]

g = pd.read_csv(f"{OUT}/scores_grid.tsv.gz", sep="\t")
m = pd.read_csv(f"{OUT}/scores_main.tsv.gz", sep="\t")
m = m[(m.rep <= 50) & m.method.isin(list(KEY)) & (m.lead_threshold1 == 1e-3)]
s = pd.concat([g[g.method.isin(list(KEY))], m], ignore_index=True)
cut = pd.read_csv(f"{OUT}/thresholds.tsv", sep="\t").set_index("key").cutoff
s["call"] = False
for meth, k in KEY.items():
    i = s.method == meth
    s.loc[i, "call"] = s.loc[i, "pp4"].fillna(-1) > cut[k]
fb = s[s.method == "V3_gate_conditional"].copy()
fb["method"] = "block_coloc_fb"
fb["pp4"] = fb.pp4.fillna(fb.pp4_whole)
fb["call"] = fb.pp4 > cut["V3_gate_conditional_fb@1e-03"]
s.loc[s.method == "V3_gate_conditional", "method"] = "block_coloc"
s = pd.concat([s, fb], ignore_index=True)
s["call"] = s["call"].astype(bool)

def cell(d, h1, h2, tier):
    return d[np.isclose(d.h2_gwas, h1) & np.isclose(d.h2_eqtl, h2) & (d.tier == tier)]

rows = []
t = s[s.set == "test"]
for h1, h2, tier in CELLS:
    d = cell(t, h1, h2, tier)
    p, n = d[d.scenario.isin(POS)], d[d.scenario.isin(NEG)]
    base = dict(h2_gwas=h1, h2_eqtl=h2, tier=tier)
    for meth in ["block_coloc", "block_coloc_fb", "whole_abf", "susie_fallback", "signals_mask", "signals_cond"]:
        pm, nm = p[p.method == meth], n[n.method == meth]
        rows.append(dict(base, kind="rate", a=meth, b="", n=len(pm), tpr=pm.call.mean(), fpr=nm.call.mean(),
                         no_result=pm.pp4.isna().mean()))
    for a, b in PAIRS:
        A = p[p.method == a].set_index(["window", "scenario", "rep"]).call
        B = p[p.method == b].set_index(["window", "scenario", "rep"]).call.reindex(A.index)
        b10, b01, N = int((A & ~B).sum()), int((~A & B).sum()), len(A)
        dl = (b10 - b01) / N
        se = np.sqrt(((b10 + b01) / N - dl ** 2) / N)
        rows.append(dict(base, kind="paired", a=a, b=b, n=N, diff=dl, lo=dl - 1.96 * se, hi=dl + 1.96 * se,
                         b10=b10, b01=b01, p=binomtest(b10, b10 + b01, 0.5).pvalue if b10 + b01 else 1.0))
pd.DataFrame(rows).to_csv(f"{OUT}/grid_paired.tsv", sep="\t", index=False)

v = s[(s.method == "block_coloc") & (s.tier == "low")]
lead = []
for h2 in [0.0075, 0.015, 0.03, 0.05]:
    q = v[np.isclose(v.h2_gwas, 0.003) & np.isclose(v.h2_eqtl, h2) & (v.n_causal_eqtl == 2)]
    lead.append(dict(trait=2, h2=h2, n=len(q), p_ge2=(q.n_leads2 >= 2).mean()))
for h1 in [0.001, 0.003, 0.01, 0.02]:
    q = v[np.isclose(v.h2_gwas, h1) & np.isclose(v.h2_eqtl, 0.015) & (v.n_causal_gwas == 2)]
    lead.append(dict(trait=1, h2=h1, n=len(q), p_ge2=(q.n_leads1 >= 2).mean()))
pd.DataFrame(lead).to_csv(f"{OUT}/grid_leads.tsv", sep="\t", index=False)
print(pd.DataFrame(lead).round(3).to_string(index=False))
