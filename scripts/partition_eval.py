#!/usr/bin/env python3
"""Evaluate the partitioned (factorised) coloc posterior against stage-1 block-coloc.

Input: partition_loci.tsv.gz written by scripts/bulk_h2ref_v1/block_partition.R
(optionally prior_rule_loci.tsv.gz from the same cell, for the reproduction check).

Statistics compared (all on the same loci):
  whole        whole-region coloc.abf PP4
  bc_stage1    block-coloc as in the paper: max PP4 over joint blocks; whole-region
               PP4 when the gate fails or there is no valid block
  joint_prod   1 - prod (1 - PP4_b) over joint blocks only (else as bc_stage1):
               isolates "max vs combine"
  region       region-level posterior from ALL disjoint parts (joint blocks,
               single-trait neighbourhoods, rest of region); loci that fail the
               gate have one part (the whole region), so region = whole there
  region_nr    as region but without the "rest" part (sensitivity)

For region / region_nr the full H0-H4 are available:
  A = prod(1-PP4_b), C1 = prod(PP0_b+PP2_b), C2 = prod(PP0_b+PP1_b), C0 = prod(PP0_b)
  H0 = C0, H1 = C2-C0, H2 = C1-C0, H3 = A-C1-C2+C0, H4 = 1-A

Usage: python3 partition_eval.py <partition_loci.tsv.gz> [prior_rule_loci.tsv.gz] [--out DIR]
"""
import sys, os, math
import numpy as np
import pandas as pd
from scipy.stats import binomtest, norm, t as tdist

POS = ["A_A", "A_AB", "AB_AB", "AB_AC"]
CAL_NEG = ["A_B", "AB_CD"]
EXTRA_NEG = ["gwas_only", "null_both"]
LABEL = {"A_A": "single_shared", "A_AB": "single_double_shared", "AB_AB": "double_shared2",
         "AB_AC": "double_shared1", "A_B": "single_distinct", "AB_CD": "double_distinct",
         "gwas_only": "gwas_only", "null_both": "null_both"}
STATS = ["whole", "bc_stage1", "joint_prod", "region", "region_nr"]
RULE_STATS = {"whole": ("whole_pp3", "whole"), "bc_stage1": ("bc_stage1_pp3", "bc_stage1"),
              "region": ("region_h3", "region"), "region_nr": ("region_nr_h3", "region_nr")}
B_BOOT = 2000


def parse_parts(s):
    out = []
    for p in s.split(";"):
        typ, size, pp, fl = p.split("|")
        v = np.array([float("nan") if x == "NA" else float(x) for x in pp.split(",")])
        f = [int(x) for x in fl.split(",")]
        out.append((typ, int(size), v, f))
    return out


def region_post(parts):
    """Region-level H0..H4 from disjoint parts; parts with NA PPs are skipped."""
    A = C1 = C2 = C0 = 1.0
    used = 0
    for _, _, v, _ in parts:
        if not np.all(np.isfinite(v)):
            continue
        used += 1
        A *= v[0] + v[1] + v[2] + v[3]          # = 1 - PP4, without cancellation
        C1 *= v[0] + v[2]
        C2 *= v[0] + v[1]
        C0 *= v[0]
    if not used:
        return [np.nan] * 5
    h = [C0, C2 - C0, C1 - C0, A - C1 - C2 + C0, 1.0 - A]
    return [min(max(x, 0.0), 1.0) for x in h]


def add_stats(d):
    rows = []
    for _, r in d.iterrows():
        parts = parse_parts(r["parts"])
        joint = [p for p in parts if p[0] == "joint" and np.all(np.isfinite(p[2]))]
        gated_block = bool(r["gate"]) and np.isfinite(r["stage1_pp4"])
        bc = r["stage1_pp4"] if gated_block else r["whole_pp4"]
        bc3 = r["stage1_pp3"] if gated_block else r["whole_pp3"]
        if gated_block and joint:
            jp = 1.0 - np.prod([p[2][:4].sum() for p in joint])
        else:
            jp = bc
        reg = region_post(parts)
        reg_nr = region_post([p for p in parts if p[0] != "rest"]) if len(parts) > 1 else reg
        types = [p[0] for p in parts]
        rows.append(dict(bc_stage1=bc, bc_stage1_pp3=bc3, joint_prod=jp,
                         region=reg[4], region_h3=reg[3], region_h0=reg[0], region_h1=reg[1], region_h2=reg[2],
                         region_nr=reg_nr[4], region_nr_h3=reg_nr[3],
                         n_joint=types.count("joint"), n_t1=types.count("trait1"), n_t2=types.count("trait2"),
                         has_rest=int("rest" in types),
                         n_na_parts=sum(not np.all(np.isfinite(p[2])) for p in parts)))
    s = pd.DataFrame(rows, index=d.index)
    d = pd.concat([d, s], axis=1)
    d["whole"] = d["whole_pp4"]
    return d


def threshold(x):
    """5% matched FPR: the floor(0.05 n) + 1-th largest score (NA -> -1)."""
    v = np.sort(np.nan_to_num(np.asarray(x, float), nan=-1.0))[::-1]
    return v[int(math.floor(0.05 * len(v)))]


def wilson(k, n, z=1.96):
    if n == 0:
        return (np.nan, np.nan)
    p = k / n; den = 1 + z * z / n
    c = (p + z * z / (2 * n)) / den; h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / den
    return (c - h, c + h)


def auc(pos, neg):
    pos = np.nan_to_num(np.asarray(pos, float), nan=-1.0); neg = np.nan_to_num(np.asarray(neg, float), nan=-1.0)
    allv = np.concatenate([pos, neg]); ranks = pd.Series(allv).rank().values
    return (ranks[:len(pos)].sum() - len(pos) * (len(pos) + 1) / 2) / (len(pos) * len(neg))


def mcnemar(a, b):
    """a, b boolean arrays on the same loci. Returns diff (pp), Wald CI, exact p."""
    a = np.asarray(a, bool); b = np.asarray(b, bool); n = len(a)
    b10 = int(np.sum(a & ~b)); b01 = int(np.sum(~a & b))
    d = (b10 - b01) / n
    se = math.sqrt(max((b10 + b01) / n - d * d, 0) / n)
    p = binomtest(b10, b10 + b01, 0.5).pvalue if b10 + b01 > 0 else 1.0
    return 100 * d, 100 * (d - 1.96 * se), 100 * (d + 1.96 * se), p, b10, b01


def window_level(a, b, win):
    """Per-window differences in detection rate; t interval and window bootstrap."""
    df = pd.DataFrame({"a": np.asarray(a, float), "b": np.asarray(b, float), "w": np.asarray(win)})
    g = df.groupby("w")[["a", "b"]].mean()
    diffs = 100 * (g["a"] - g["b"]).values
    k = len(diffs)
    m = diffs.mean(); se = diffs.std(ddof=1) / math.sqrt(k) if k > 1 else np.nan
    q = tdist.ppf(0.975, k - 1) if k > 1 else np.nan
    rng = np.random.default_rng(20261113)
    # cluster bootstrap over windows of the pooled difference (windows weighted by n)
    wins = g.index.values; n_w = df.groupby("w").size().reindex(wins).values
    sa = df.groupby("w")["a"].sum().reindex(wins).values; sb = df.groupby("w")["b"].sum().reindex(wins).values
    boot = []
    for _ in range(B_BOOT):
        i = rng.integers(0, k, k)
        boot.append(100 * (sa[i].sum() - sb[i].sum()) / n_w[i].sum())
    lo, hi = np.percentile(boot, [2.5, 97.5])
    return dict(window_diffs=diffs, mean=m, t_lo=m - q * se, t_hi=m + q * se, boot_lo=lo, boot_hi=hi,
                n_pos=int(np.sum(diffs > 0)), k=k)


def fmt(x, nd=1):
    return "NA" if x is None or (isinstance(x, float) and not np.isfinite(x)) else f"{x:.{nd}f}"


def evaluate(d, title):
    L = [f"# {title}", ""]
    cal = d[d.set == "calibration"]; test = d[d.set == "test"]
    # diagnostics
    g = d[d.scenario.isin(POS + CAL_NEG + EXTRA_NEG)]
    L += ["## Diagnostics", "",
          "| scenario | n | gate passed | >=2 joint blocks | >=1 single-trait part | mean parts (gated) | overlap SNPs > 0 (gated) | NA parts |",
          "|---|---|---|---|---|---|---|---|"]
    for sc in POS + CAL_NEG + EXTRA_NEG:
        x = g[g.scenario == sc]
        if not len(x):
            continue
        gx = x[x.gate.astype(bool)]
        L.append(f"| {LABEL[sc]} | {len(x)} | {100*x.gate.astype(bool).mean():.1f}% | {100*(x.n_joint>=2).mean():.1f}% | "
                 f"{100*((x.n_t1+x.n_t2)>=1).mean():.1f}% | {fmt(gx.n_parts.mean() if len(gx) else np.nan,2)} | "
                 f"{fmt(100*(gx.overlap_snps>0).mean() if len(gx) else np.nan)}% | {int(x.n_na_parts.sum())} |")
    L.append("")
    # thresholds
    thr = {s: threshold(cal[cal.scenario.isin(CAL_NEG)][s]) for s in STATS}
    L += ["## Matched 5% FPR (thresholds from calibration windows)", "",
          "Thresholds: " + ", ".join(f"{s} {thr[s]:.4g}" for s in STATS), "",
          "Detection rate (%) on test windows; type I error (%) for negatives.", "",
          "| scenario | n | " + " | ".join(STATS) + " |", "|---|---|" + "---|" * len(STATS)]
    for sc in POS + ["pooled"] + CAL_NEG + EXTRA_NEG:
        x = test[test.scenario.isin(POS)] if sc == "pooled" else test[test.scenario == sc]
        if not len(x):
            continue
        cells = []
        for s in STATS:
            k = int((np.nan_to_num(x[s].values, nan=-1) > thr[s]).sum()); n = len(x)
            lo, hi = wilson(k, n)
            cells.append(f"{100*k/n:.1f} ({100*lo:.1f}-{100*hi:.1f})")
        L.append(f"| {LABEL.get(sc, sc)} | {len(x)} | " + " | ".join(cells) + " |")
    pos = test[test.scenario.isin(POS)]; neg = test[test.scenario.isin(CAL_NEG)]
    L += ["", "AUC (test; positives vs single/double_distinct): " +
          ", ".join(f"{s} {auc(pos[s], neg[s]):.3f}" for s in STATS), ""]
    # paired comparisons
    L += ["## Paired comparisons on the same positive test loci", "",
          "Locus-level: McNemar exact test, paired Wald CI (treats loci as independent).",
          "Window-level: mean of the 6 per-window differences with a t interval, and a cluster bootstrap over windows "
          f"({B_BOOT} resamples) of the pooled difference. With 6 windows these intervals are wide by construction; "
          "they describe generalisation to new regions, the locus-level ones only to these 6 windows.", "",
          "| comparison | locus-level diff (95% CI), p | window-level mean (t 95% CI) | window bootstrap 95% CI | windows > 0 | per-window diffs |",
          "|---|---|---|---|---|---|"]
    det = {s: np.nan_to_num(pos[s].values, nan=-1) > thr[s] for s in STATS}
    comps = [("region", "whole"), ("region", "bc_stage1"), ("joint_prod", "bc_stage1"),
             ("region_nr", "region"), ("bc_stage1", "whole")]
    for a, b in comps:
        dd, lo, hi, p, b10, b01 = mcnemar(det[a], det[b])
        w = window_level(det[a], det[b], pos.window.values)
        L.append(f"| {a} − {b} | {dd:+.1f} ({lo:+.1f}, {hi:+.1f}), p = {p:.2g} | {w['mean']:+.1f} ({w['t_lo']:+.1f}, {w['t_hi']:+.1f}) | "
                 f"({w['boot_lo']:+.1f}, {w['boot_hi']:+.1f}) | {w['n_pos']}/{w['k']} | " +
                 ", ".join(f"{x:+.1f}" for x in w["window_diffs"]) + " |")
    L.append("")
    # traceCB rule
    L += ["## traceCB rule (PP4 > 0.7 colocalised, PP3 > 0.7 independent, else undetermined), test windows", "",
          "| statistic | colocalised loci: coloc / indep / undet | single_distinct: indep (coloc) | double_distinct: indep (coloc) | gwas_only: coloc |",
          "|---|---|---|---|---|"]
    for s, (c3, c4) in RULE_STATS.items():
        def rule(x):
            p4 = np.nan_to_num(x[c4].values, nan=0); p3 = np.nan_to_num(x[c3].values, nan=0)
            return 100 * np.mean(p4 > 0.7), 100 * np.mean(p3 > 0.7), 100 * np.mean((p4 <= 0.7) & (p3 <= 0.7))
        a = rule(pos); sd = rule(test[test.scenario == "A_B"]); dd_ = rule(test[test.scenario == "AB_CD"])
        go = rule(test[test.scenario == "gwas_only"]) if (test.scenario == "gwas_only").any() else (np.nan,) * 3
        L.append(f"| {s} | {a[0]:.1f} / {a[1]:.1f} / {a[2]:.1f} | {sd[1]:.1f} ({sd[0]:.1f}) | {dd_[1]:.1f} ({dd_[0]:.1f}) | {fmt(go[0])} |")
    L.append("")
    # how the region posterior treats no-valid-block loci
    nv = test[test.gate.astype(bool) & ~np.isfinite(test.stage1_pp4)]
    if len(nv):
        L += ["## Loci with no valid block (test windows)", "",
              "| scenario | n | whole PP4 median | region PP4 median | region PP3 median | region PP3 > 0.7 |", "|---|---|---|---|---|---|"]
        for sc in POS + CAL_NEG:
            x = nv[nv.scenario == sc]
            if len(x):
                L.append(f"| {LABEL[sc]} | {len(x)} | {x.whole.median():.3f} | {x.region.median():.3f} | {x.region_h3.median():.3f} | {100*(x.region_h3>0.7).mean():.1f}% |")
        L.append("")
    return "\n".join(L), thr


def main(argv):
    args = [a for a in argv if not a.startswith("--")]
    out = argv[argv.index("--out") + 1] if "--out" in argv else os.path.dirname(os.path.abspath(args[0]))
    if "--out" in argv:
        args = [a for a in args if a != out]
    d = pd.read_csv(args[0], sep="\t")
    d = add_stats(d)
    txt, thr = evaluate(d, f"Partitioned coloc posterior: h2_gwas = {d.h2_gwas.iloc[0]}, h2_eqtl = {d.h2_eqtl.iloc[0]}")
    if len(args) > 1:
        pr = pd.read_csv(args[1], sep="\t")
        m = d.merge(pr[["window", "scenario", "rep", "bc_orig_pp4", "whole_pp4"]], on=["window", "scenario", "rep"],
                    suffixes=("", "_prior"))
        a = m.stage1_pp4.values; b = m.bc_orig_pp4.values
        both_na = np.isnan(a) & np.isnan(b); fin = np.isfinite(a) & np.isfinite(b)
        txt = (f"Reproduction check against prior_rule.R: {len(m)} loci matched; block PP4 max |diff| = "
               f"{np.max(np.abs(a[fin]-b[fin])) if fin.any() else float('nan'):.3g}; NA pattern identical: "
               f"{bool(np.all(both_na | fin))}; whole PP4 max |diff| = {np.max(np.abs(m.whole_pp4 - m.whole_pp4_prior)):.3g}\n\n") + txt
    os.makedirs(out, exist_ok=True)
    with open(os.path.join(out, "PARTITION.md"), "w") as f:
        f.write(txt + "\n")
    cols = ["window", "set", "scenario", "truth", "rep", "gate", "n_parts", "n_joint", "n_t1", "n_t2", "has_rest",
            "whole", "bc_stage1", "joint_prod", "region", "region_h3", "region_nr", "region_nr_h3", "bc_stage1_pp3", "whole_pp3"]
    d[cols].to_csv(os.path.join(out, "partition_stats.tsv.gz"), sep="\t", index=False)
    print(txt)


if __name__ == "__main__":
    main(sys.argv[1:])
