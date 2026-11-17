#!/usr/bin/env python3
"""Why is block-coloc insensitive to the LD source? Same loci under in-sample and reference LD.
Branch transitions (G0 gate not passed / G1 valid block / G2 no valid block) and the resulting PP4 changes,
compared with coloc.signals (cond). Usage: python3 ld_mechanism.py  (env BULK_OUT)"""
import os, math, numpy as np, pandas as pd
ROOT = os.environ.get("BULK_OUT", "../clike_ldblock_coloc_handoff_20260811/outputs").rstrip("/") + "/"
POS = ["A_A", "A_AB", "AB_AB", "AB_AC"]; NEG = ["A_B", "AB_CD"]; ID = ["window", "scenario", "rep"]
def cutoff(v):
    v = np.sort(np.nan_to_num(np.asarray(v, float), nan=-1.0))[::-1]; return v[int(math.floor(0.05 * len(v)))]
def frame(path, hg, he):
    d = pd.read_csv(path, sep="\t", low_memory=False)
    d = d[np.isclose(d.h2_gwas, hg) & np.isclose(d.h2_eqtl, he) & (d.tier == "low")]
    sets = d.drop_duplicates("window").set_index("window")["set"]
    b = d[(d.method == "V3_gate_conditional") & np.isclose(d.lead_threshold1, 1e-3) & np.isclose(d.lead_threshold2, 1e-3)].set_index(ID)
    s = d[d.method == "signals_cond"].set_index(ID)["pp4"].rename("cond")
    x = b[["pp4", "pp4_whole", "branch", "status", "n_leads1", "n_leads2", "unit_has_shared"]].join(s).reset_index()
    x["set"] = x.window.map(sets); x["bc"] = np.where(np.isfinite(x.pp4), x.pp4, x.pp4_whole)
    x["g"] = np.where(x.branch == "whole", "G0", np.where(x.status == "no_valid_block", "G2", "G1"))
    return x.set_index(ID)
for lab, ref, ins, hg in (("default", "bulk_h2ref_v1/scores_main.tsv.gz", "bulk_h2ref_v1_insample/scores_main.tsv.gz", 0.003),
                          ("h2_gwas 2%", "bulk_h2ref_v1/scores_grid.tsv.gz", "bulk_h2ref_v1_insample_g0.02_e0.015/scores_main.tsv.gz", 0.02)):
    r, i = frame(ROOT + ref, hg, 0.015), frame(ROOT + ins, hg, 0.015)
    idx = r.index.intersection(i.index); r, i = r.loc[idx], i.loc[idx]
    sc = idx.get_level_values("scenario"); st = r.set.values
    tr = {k: cutoff(x.bc.values[(st == "calibration") & np.isin(sc, NEG)]) for k, x in (("r", r), ("i", i))}
    tc = {k: cutoff(x.cond.values[(st == "calibration") & np.isin(sc, NEG)]) for k, x in (("r", r), ("i", i))}
    print(f"\n######## {lab}: {len(idx)} loci; block-coloc cut-off in-sample {tr['i']:.3f} / ref {tr['r']:.3f}; signals_cond {tc['i']:.3f} / {tc['r']:.3f}")
    for grp, m in (("positives (test)", (st == "test") & np.isin(sc, POS)), ("negatives (test)", (st == "test") & np.isin(sc, NEG)), ("all loci", np.ones(len(idx), bool))):
        T = pd.crosstab(pd.Series(i.g.values[m], name="in-sample"), pd.Series(r.g.values[m], name="reference"), normalize="all").round(3) * 100
        print(f"--- {grp}: branch transition (% of loci)\n{T.to_string()}")
    m = (st == "test") & np.isin(sc, POS)
    same = (i.g.values == r.g.values)
    d_bc = np.abs(r.bc.values - i.bc.values); d_cd = np.abs(np.nan_to_num(r.cond.values, nan=0) - np.nan_to_num(i.cond.values, nan=0))
    for name, mm in (("branch unchanged, G0/G2 (PP4 = whole in both)", same & (i.g.values != "G1")), ("branch unchanged, G1", same & (i.g.values == "G1")), ("branch changed", ~same)):
        k = m & mm
        if k.sum() == 0: continue
        ci = i.bc.values[k] > tr["i"]; cr = r.bc.values[k] > tr["r"]
        print(f"  positives, {name}: {100*k.sum()/m.sum():.1f}% of positives; median |dPP4| {np.median(d_bc[k]):.3g}, 90th pct {np.quantile(d_bc[k],.9):.3g}; "
              f"call changed {100*np.mean(ci!=cr):.1f}% (gained {int((cr&~ci).sum())}, lost {int((ci&~cr).sum())})")
    ci = i.cond.values[m] > tc["i"]; cr = r.cond.values[m] > tc["r"]
    print(f"  signals_cond, positives: median |dPP4| {np.median(d_cd[m]):.3g}, 90th pct {np.quantile(d_cd[m],.9):.3g}; call changed {100*np.mean(ci!=cr):.1f}% (gained {int((cr&~ci).sum())}, lost {int((ci&~cr).sum())})")
    ch = m & ~same
    if ch.sum():
        tt = pd.Series([f"{a}->{b}" for a, b in zip(i.g.values[ch], r.g.values[ch])]).value_counts()
        print("  positives with branch change, by transition:", dict(tt))
        for t_ in tt.index[:4]:
            a, b = t_.split("->"); k = m & (i.g.values == a) & (r.g.values == b)
            print(f"    {t_}: n={k.sum()}, median PP4 in-sample {np.median(i.bc.values[k]):.3f} -> ref {np.median(r.bc.values[k]):.3f}; called in-sample {100*np.mean(i.bc.values[k]>tr['i']):.0f}% -> ref {100*np.mean(r.bc.values[k]>tr['r']):.0f}%")
    # negatives
    n_ = (st == "test") & np.isin(sc, NEG)
    print(f"  negatives: block-coloc called in-sample {100*np.mean(i.bc.values[n_]>tr['i']):.1f}% / ref {100*np.mean(r.bc.values[n_]>tr['r']):.1f}%; "
          f"signals_cond {100*np.mean(i.cond.values[n_]>tc['i']):.1f}% / {100*np.mean(r.cond.values[n_]>tc['r']):.1f}%; "
          f"median negative PP4 signals_cond in-sample {np.nanmedian(i.cond.values[n_]):.3f} / ref {np.nanmedian(r.cond.values[n_]):.3f}; 95th pct {np.nanquantile(i.cond.values[n_],.95):.3f} / {np.nanquantile(r.cond.values[n_],.95):.3f}")
