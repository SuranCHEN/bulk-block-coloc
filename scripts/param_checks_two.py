#!/usr/bin/env python3
"""Monte Carlo for Proposition 5 in doc/bulk/09_parameters.md: lead counts when a
trait has two causal variants (|r| < 0.2, random signs, Var(g) rescaled), using
Steps 1-2 of block-coloc (greedy clump r2 < 0.1, then conditional p < alpha).
Truth and analysis both use the 1000G EUR reference LD (approximation).
Usage: python3 scripts/bulk_h2ref_v1/param_checks_two.py LAMBDA [LAMBDA ...]
Appends to outputs/bulk_h2ref_v1/param_checks_two.tsv. numpy + pandas only.
"""
import os, sys
import numpy as np, pandas as pd
from statistics import NormalDist
REF, OUT = "data/ref_ld_1000g_eur_20261027", "outputs/bulk_h2ref_v1/param_checks_two.tsv"
ALPHA, RHO_C, RHO_B, NDRAW = 1e-3, 0.1, 0.4, 500
zc = NormalDist().inv_cdf(1 - ALPHA / 2)
lams = [float(x) for x in sys.argv[1:]]
rng = np.random.default_rng(int(1000 * sum(lams)) + 7)

def leads(z, R):
    cand = np.where(np.abs(z) > zc)[0]
    cand = cand[np.argsort(-np.abs(z[cand]))]
    marg = []
    for j in cand:
        if all(R[j, l] ** 2 < RHO_C for l in marg):
            marg.append(j)
    if not marg:
        return []
    L = [marg[0]]
    for j in marg[1:]:
        S = np.array(L)
        RSS = R[np.ix_(S, S)] + 1e-6 * np.eye(len(S))
        w = np.linalg.solve(RSS, R[S, j])
        v = 1 - R[j, S] @ w
        if v <= 1e-8:
            continue
        if abs((z[j] - w @ z[S]) / np.sqrt(v)) > zc:
            L.append(j)
    return L

rows = []
for w in sorted(int(x) for x in os.listdir(REF) if x.isdigit()):
    snp = pd.read_csv(f"{REF}/{w}/ref_snps.tsv", sep="\t"); M = len(snp)
    R = np.fromfile(f"{REF}/{w}/ld_ref.f64").reshape(M, M)
    ev, V = np.linalg.eigh(R); A = V * np.sqrt(np.clip(ev, 0, None))
    elig = np.where(snp.ref_maf >= 0.05)[0]
    for lam in lams:
        n1 = n2 = both = 0
        for _ in range(NDRAW):
            while True:
                a, b = rng.choice(elig, 2, replace=False)
                if abs(R[a, b]) < 0.2:
                    break
            s = rng.choice([-1, 1], 2)
            c = (1 + s[0] * s[1] * R[a, b]) ** -0.5
            z = lam * c * (s[0] * R[:, a] + s[1] * R[:, b]) + A @ rng.standard_normal(M)
            L = leads(z, R)
            n1 += len(L) >= 1; n2 += len(L) >= 2
            if len(L) >= 2:
                ta = [l for l in L if R[l, a] ** 2 >= RHO_B]; tb = [l for l in L if R[l, b] ** 2 >= RHO_B]
                both += any(x != y for x in ta for y in tb)
        rows.append(dict(window=w, lam=lam, n=NDRAW, P_ge1=n1 / NDRAW, P_ge2=n2 / NDRAW, P_both_tagged=both / NDRAW))
T = pd.DataFrame(rows)
T.to_csv(OUT, sep="\t", index=False, mode="a", header=not os.path.exists(OUT))
for lam, q in T.groupby("lam"):
    out = {k: (q[k].mean(), np.sqrt((q[k] * (1 - q[k])).mean() / (NDRAW * len(q)))) for k in ["P_ge1", "P_ge2", "P_both_tagged"]}
    print(lam, {k: f"{m:.3f} ({se:.3f})" for k, (m, se) in out.items()})
