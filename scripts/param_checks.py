#!/usr/bin/env python3
"""Numerical checks for doc/bulk/05_parameters.md.

Uses only the 1000 Genomes EUR reference LD (public) and aggregate stage-1
output; writes aggregate tables to outputs/bulk_h2ref_v1/param_checks_*.tsv.
Approximation: z ~ N(lambda * R[:, c], R) with R = reference LD (the simulation
itself uses in-sample UKBB LD; the two correlate 0.95-0.99, see Model section 7).

Run from the project root:  python3 scripts/bulk_h2ref_v1/param_checks.py
Needs numpy and pandas only.
"""
import os
import numpy as np
import pandas as pd
from statistics import NormalDist

REF = "data/ref_ld_1000g_eur_20261027"
OUT = "outputs/bulk_h2ref_v1"
ALPHA = {"1e-3": 1e-3, "1e-4": 1e-4, "1e-5": 1e-5, "5e-8": 5e-8}
LAMS = [1.94, 2.75, 3.16, 3.9, 4.47, 5.06, 5.48, 5.56, 7.25, 7.76]
N_NULL, N_ALT = 20000, 4000
zcut = {a: NormalDist().inv_cdf(1 - v / 2) for a, v in ALPHA.items()}
rng = np.random.default_rng(20261028)

null_rows, alt_rows, ld = [], [], {}
for w in sorted(int(x) for x in os.listdir(REF) if x.isdigit()):
    snp = pd.read_csv(f"{REF}/{w}/ref_snps.tsv", sep="\t")
    M = len(snp)
    R = np.fromfile(f"{REF}/{w}/ld_ref.f64").reshape(M, M)
    ld[w] = (snp, R)
    ev, V = np.linalg.eigh(R)
    A = V * np.sqrt(np.clip(ev, 0, None))          # z = A u  ~  N(0, R)
    k = np.abs(ev)
    meff = float(np.sum((k >= 1) + (k - np.floor(k))))  # Li & Ji (2005)
    mx = np.abs(A @ rng.standard_normal((M, N_NULL))).max(0)
    row = dict(window=w, M=M, Meff_LiJi=round(meff, 1))
    for a in ALPHA:
        row["P_any_lead_null_" + a] = float(np.mean(mx > zcut[a]))
    null_rows.append(row)
    elig = np.where(snp.ref_maf >= 0.05)[0]
    Z0 = A @ rng.standard_normal((M, N_ALT))
    cs = rng.choice(elig, N_ALT)
    for lam in LAMS:
        Z = Z0 + lam * R[:, cs]
        lead = np.abs(Z).argmax(0)
        top = np.abs(Z).max(0)
        r2 = R[lead, cs] ** 2
        for a in ["1e-3", "1e-4", "5e-8"]:
            det = top > zcut[a]
            alt_rows.append(dict(
                window=w, lam=lam, alpha=a, P_lead_found=float(det.mean()),
                P_causal_in_lead_nbhd=float(np.mean(r2 >= 0.4)),
                P_causal_in_lead_nbhd_given_found=float(np.mean(r2[det] >= 0.4)) if det.any() else np.nan))

T1 = pd.DataFrame(null_rows)
T2 = pd.DataFrame(alt_rows).groupby(["lam", "alpha"]).mean(numeric_only=True).drop(columns="window").reset_index()
T1.to_csv(f"{OUT}/param_checks_null.tsv", sep="\t", index=False)
T2.to_csv(f"{OUT}/param_checks_lead.tsv", sep="\t", index=False)

# ---- stage-1 observations (aggregate)
d = pd.read_csv(f"{OUT}/loci_main.tsv.gz", sep="\t")
s = pd.read_csv(f"{OUT}/scores_main.tsv.gz", sep="\t")
zs = lambda x: [float(t) for t in str(x).replace(",", ";").split(";") if t not in ("nan", "")]
obs = []
for sc in ["A_A", "A_AB", "AB_AB"]:
    q = d[d.scenario == sc]
    for tr, col in [("gwas", "causal_abs_z_gwas"), ("eqtl", "causal_abs_z_eqtl")]:
        z = np.concatenate([zs(x) for x in q[col]])
        obs.append(dict(quantity=f"mean |z| causal {tr}", scenario=sc, value=round(z.mean(), 3), n=len(z)))
v = s[s.method == "V3_gate_conditional"]
for (t1, t2), q in v.groupby(["lead_threshold1", "lead_threshold2"]):
    if t1 != t2:
        continue
    for sc in ["null_both", "gwas_only", "A_A"]:
        x = q[q.scenario == sc]
        for nm, val in [("P(n_leads1>=1)", np.mean(x.n_leads1 >= 1)), ("P(n_leads2>=1)", np.mean(x.n_leads2 >= 1)),
                        ("P(n_leads1>=2)", np.mean(x.n_leads1 >= 2)), ("P(n_leads2>=2)", np.mean(x.n_leads2 >= 2)),
                        ("P(block branch)", np.mean(x.branch != "whole"))]:
            obs.append(dict(quantity=f"{nm} alpha={t1:g}", scenario=sc, value=round(float(val), 4), n=len(x)))
x = v[(v.lead_threshold1 == 1e-3) & (v.lead_threshold2 == 1e-3) & (v.scenario == "A_A")]
x = x.merge(d[["window", "scenario", "rep", "tier", "causal_eqtl"]], on=["window", "scenario", "rep", "tier"])
r2 = []
for _, r in x.iterrows():
    if pd.isna(r.leads2):
        continue
    snp, R = ld[r.window]
    idx = {s_: i for i, s_ in enumerate(snp.snp)}
    c, l = idx.get(r.causal_eqtl), idx.get(str(r.leads2).replace(",", ";").split(";")[0])
    if c is not None and l is not None:
        r2.append(R[l, c] ** 2)
r2 = np.array(r2)
obs.append(dict(quantity="P(r2(first trait-2 lead, causal) >= 0.4) alpha=0.001", scenario="A_A", value=round(float(np.mean(r2 >= 0.4)), 4), n=len(r2)))
O = pd.DataFrame(obs)
O.to_csv(f"{OUT}/param_checks_observed.tsv", sep="\t", index=False)

pd.set_option("display.width", 200)
print(T1.to_string(index=False))
print(T1.drop(columns="window").mean().round(4).to_string())
print(T2.round(3).to_string(index=False))
print(O.to_string(index=False))
