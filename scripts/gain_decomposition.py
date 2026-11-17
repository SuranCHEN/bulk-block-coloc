#!/usr/bin/env python3
"""Where does block-coloc's gain over whole-region coloc.abf come from?

Same loci, same reference LD, existing per-locus scores (no new simulation).
block-coloc = V3_gate_conditional at lead threshold 1e-3, no valid block -> whole-region PP4.
Every positive test locus falls in exactly one branch of block-coloc:
  G0  gate not passed            -> block-coloc PP4 == whole-region PP4
  G1  block branch, valid block  -> block-coloc PP4 = max block PP4
  G2  block branch, no valid block -> block-coloc PP4 == whole-region PP4 (fallback)
The detection-rate difference is decomposed exactly as
  sum_g w_g [P(bc > t_bc | g) - P(abf > t_abf | g)]
and within G1 further into a PP4 effect at the common cut-off t_abf and a cut-off effect.
Usage: python3 gain_decomposition.py [--out FILE.md]   (env BULK_OUT = hand-off outputs folder)
"""
import os, sys, math
import numpy as np, pandas as pd

ROOT = os.environ.get("BULK_OUT", "../clike_ldblock_coloc_handoff_20260811/outputs").rstrip("/") + "/"
POS = ["A_A", "A_AB", "AB_AB", "AB_AC"]; NEG = ["A_B", "AB_CD"]
LAB = {"A_A": "single_shared", "A_AB": "single_double_shared", "AB_AB": "double_shared2", "AB_AC": "double_shared1",
       "A_B": "single_distinct", "AB_CD": "double_distinct"}
K = ["window", "scenario", "rep", "h2_gwas", "h2_eqtl", "tier"]


def cutoff(v):
    v = np.sort(np.nan_to_num(np.asarray(v, float), nan=-1.0))[::-1]
    return v[int(math.floor(0.05 * len(v)))]


def frame(path):
    d = pd.read_csv(path, sep="\t", low_memory=False)
    sets = d.drop_duplicates("window").set_index("window")["set"]
    w = d[d.method == "whole_abf"].set_index(K)[["pp4"]].rename(columns={"pp4": "abf"})
    b = d[(d.method == "V3_gate_conditional") & np.isclose(d.lead_threshold1, 1e-3) & np.isclose(d.lead_threshold2, 1e-3)].set_index(K)
    x = w.join(b[["pp4", "pp4_whole", "branch", "status", "n_leads1", "n_leads2", "unit_has_shared", "unit_size"]], how="inner")
    x = x.reset_index()
    x["set"] = x.window.map(sets)
    x["bc"] = np.where(np.isfinite(x.pp4), x.pp4, x.pp4_whole)          # adopted rule (fallback)
    x["nofb"] = x.pp4                                                    # no-fallback sensitivity version
    x["g"] = np.where(x.branch == "whole", "G0", np.where(x.status == "no_valid_block", "G2", "G1"))
    chk = np.nanmax(np.abs(x.loc[x.g != "G1", "bc"] - x.loc[x.g != "G1", "abf"]))
    return x, chk


def decompose(x, t_abf, t_bc, t_nofb, title, L):
    tp = x[(x.set == "test") & x.scenario.isin(POS)]; tn = x[(x.set == "test") & x.scenario.isin(NEG)]
    n = len(tp)
    abf = tp.abf.values > t_abf; bc = tp.bc.values > t_bc; bcc = tp.bc.values > t_abf
    nofb = np.nan_to_num(tp.nofb.values, nan=-1) > t_nofb
    L += [f"### {title}", "",
          f"Cut-offs (5% FPR on calibration windows): coloc.abf {t_abf:.4f}, block-coloc {t_bc:.4f}, block-coloc without fallback {t_nofb:.4f}. "
          f"Positive test loci: {n}. Detection: coloc.abf {100*abf.mean():.1f}%, block-coloc {100*bc.mean():.1f}% (difference {100*(bc.mean()-abf.mean()):+.1f}), "
          f"without fallback {100*nofb.mean():.1f}%.", "",
          "| branch | share of positives | coloc.abf | block-coloc | contribution to the difference (pp) | of which: PP4 change at common cut-off | of which: cut-off change |",
          "|---|---|---|---|---|---|---|"]
    tot = 0
    for g, name in (("G0", "G0 gate not passed (PP4 identical)"), ("G1", "G1 valid block (PP4 = max block)"), ("G2", "G2 no valid block (PP4 identical, fallback)")):
        m = (tp.g == g).values
        if not m.any():
            continue
        w = m.mean(); c = 100 * w * (bc[m].mean() - abf[m].mean()); tot += c
        pe = 100 * w * (bcc[m].mean() - abf[m].mean()); ce = 100 * w * (bc[m].mean() - bcc[m].mean())
        L.append(f"| {name} | {100*w:.1f}% | {100*abf[m].mean():.1f} | {100*bc[m].mean():.1f} | {c:+.2f} | {pe:+.2f} | {ce:+.2f} |")
    L.append(f"| total | 100% | {100*abf.mean():.1f} | {100*bc.mean():.1f} | {tot:+.2f} | | |")
    # why the cut-off moves: negatives at common cut-off
    L += ["", "Negatives (calibration windows) above coloc.abf's cut-off, by branch — this sets block-coloc's own cut-off:", ""]
    cn = x[(x.set == "calibration") & x.scenario.isin(NEG)]
    parts = []
    for g in ("G0", "G1", "G2"):
        m = (cn.g == g).values
        parts.append(f"{g}: share {100*m.mean():.1f}%, above t_abf coloc.abf {100*np.mean(cn.abf.values[m] > t_abf):.1f}% / block-coloc {100*np.mean(cn.bc.values[m] > t_abf):.1f}%")
    L += ["- " + "; ".join(parts), ""]
    # G1 by what triggered the gate and by scenario
    g1 = tp[tp.g == "G1"]
    if len(g1):
        trig2 = (g1.n_leads2 >= 2).values
        a1 = g1.abf.values > t_abf; b1 = g1.bc.values > t_bc
        L += [f"G1 by trigger: trait 2 has >= 2 leads: {100*trig2.mean():.1f}% of G1, gain {100*(b1[trig2].mean()-a1[trig2].mean()):+.1f} pp within; "
              f"only trait 1 has >= 2 leads: {100*(~trig2).mean():.1f}% of G1, gain {100*(b1[~trig2].mean()-a1[~trig2].mean()):+.1f} pp within.", ""]
        gain = (~a1) & b1; loss = a1 & (~b1)
        hs = g1.unit_has_shared.astype(str).str.upper().map({"TRUE": True, "FALSE": False}).values
        L += [f"G1 flips: {int(gain.sum())} loci missed by coloc.abf and called by block-coloc "
              f"(selected block contains the shared causal variant in {100*np.mean(hs[gain] == True):.1f}%); "
              f"{int(loss.sum())} called by coloc.abf and missed by block-coloc. Median block size {int(np.median(g1.unit_size))} SNPs.", ""]
    # by scenario
    L += ["| scenario | G0 / G1 / G2 share | contribution G1 (pp, within scenario) | contribution G0+G2 (pp) | total |", "|---|---|---|---|---|"]
    for s in POS:
        m = (tp.scenario == s).values
        sh = [100 * np.mean(tp.g.values[m] == g) for g in ("G0", "G1", "G2")]
        c1 = 100 * np.mean(((tp.g.values == "G1") & m)[m] * (bc[m].astype(int) - abf[m].astype(int)))
        c02 = 100 * np.mean(((tp.g.values != "G1") & m)[m] * (bc[m].astype(int) - abf[m].astype(int)))
        L.append(f"| {LAB[s]} | {sh[0]:.0f} / {sh[1]:.0f} / {sh[2]:.0f} | {c1:+.1f} | {c02:+.1f} | {c1+c02:+.1f} |")
    # no-fallback version, decomposed against block-coloc
    lostG2 = 100 * np.mean((tp.g.values == "G2") & bc)
    L += ["", f"Without fallback: G2 positives get no result; {lostG2:.1f}% of positives are called by block-coloc in G2 and lost. "
          f"Its lower cut-off ({t_nofb:.3f}) comes from the 0.0 / NA negatives in G2 (no-valid-block negatives are removed from the top of the ranking).", ""]
    return L


def main(argv):
    L = ["# Ablation: where block-coloc's gain over coloc.abf comes from", "",
         "Existing per-locus scores, reference LD. Exact decomposition of the detection-rate difference on the same positive test loci (6 test windows) by the branch block-coloc takes. "
         "In G0 and G2 block-coloc returns the whole-region PP4, so any difference there comes only from the two methods' different cut-offs.", ""]
    main_, chk = frame(ROOT + "bulk_h2ref_v1/scores_main.tsv.gz")
    L.append(f"Check: block-coloc PP4 equals whole-region PP4 outside G1, max |diff| = {chk:.2g}.\n")
    cal = main_[(main_.set == "calibration") & main_.scenario.isin(NEG)]
    t = (cutoff(cal.abf), cutoff(cal.bc), cutoff(cal.nofb))
    L += ["## Stage 1 (default setting)", ""]; decompose(main_, *t, "Default: h2_gwas 0.3%, h2_eqtl 1.5%", L)
    grid, _ = frame(ROOT + "bulk_h2ref_v1/scores_grid.tsv.gz")
    L += ["## Stage 2 (stage-1 cut-offs)", ""]
    for (hg, he, tr), x in grid.groupby(["h2_gwas", "h2_eqtl", "tier"]):
        if tr != "low":
            continue
        decompose(x, *t, f"h2_gwas {hg}, h2_eqtl {he}", L)
    for lab, f in (("Normal effects, default", "bulk_h2ref_v1_normal/scores_main.tsv.gz"), ("Normal effects, h2_eqtl 5%", "bulk_h2ref_v1_normal_g0.003_e0.05/scores_main.tsv.gz")):
        x, _ = frame(ROOT + f); c = x[(x.set == "calibration") & x.scenario.isin(NEG)]
        L += [f"## {lab} (own cut-offs)", ""]; decompose(x, cutoff(c.abf), cutoff(c.bc), cutoff(c.nofb), lab, L)
    txt = "\n".join(L) + "\n"
    if "--out" in argv:
        open(argv[argv.index("--out") + 1], "w").write(txt)
    print(txt)


if __name__ == "__main__":
    main(sys.argv[1:])
