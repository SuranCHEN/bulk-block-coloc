import numpy as np, pandas as pd, math
"""Checks of Propositions C-F (docs/11_theory.md) on stage-1/2 per-locus outputs.
Run from the bulk-block-coloc repository root; override the output folder with BULK_OUT."""
import os
R=os.environ.get("BULK_OUT","../clike_ldblock_coloc_handoff_20260811/outputs/").rstrip("/")+"/bulk_h2ref_v1/"
K=["window","scenario","rep","h2_gwas","h2_eqtl","tier"]
def scores(f):
    d=pd.read_csv(R+f,sep="\t",low_memory=False)
    w=d[d.method=="whole_abf"].set_index(K)["pp4"].rename("whole")
    b=d[(d.method=="V3_gate_conditional")&(d.lead_threshold1==1e-3)&(d.lead_threshold2==1e-3)].copy()
    b["bc"]=b.pp4.where(np.isfinite(b.pp4),b.pp4_whole); b=b.set_index(K)["bc"]
    return pd.concat([w,b],axis=1).reset_index()
loc=pd.concat([pd.read_csv(R+"loci_main.tsv.gz",sep="\t"),pd.read_csv(R+"loci_grid.tsv.gz",sep="\t")])
sc=pd.concat([scores("scores_main.tsv.gz"),scores("scores_grid.tsv.gz")])
d=loc.merge(sc,on=K).drop_duplicates(K)
T_W, T_B = 0.141194270215242, 0.144807117782326
lam=lambda N,h2,Kc: math.sqrt(N*(h2/Kc)/(1-h2/Kc))
# ---- 1. resolution limit: single_distinct (A in trait 1, B in trait 2), all windows
x=d[d.scenario=="A_B"].copy()
x["lam2"]=[lam(1000,h,1) for h in x.h2_eqtl]; x["lam1"]=[lam(20000,h,1) for h in x.h2_gwas]
x["r"]=pd.to_numeric(x.pair_abs_r,errors="coerce")
x["xmin"]=np.minimum(x.lam1**2,x.lam2**2)*(1-x.r**2)
thr=2*math.log(1e-5/(1e-4*1e-4))
print("2 log(p12/(p1 p2)) =",round(thr,2))
bins=[0,6,10,thr,18,25,40,1e3]
x["bin"]=pd.cut(x.xmin,bins)
g=x.groupby("bin",observed=True).agg(n=("whole","size"),r_med=("r","median"),pp4_med=("whole","median"),
    fp_whole=("whole",lambda v:100*np.mean(v>T_W)),pp4_gt_half=("whole",lambda v:100*np.mean(v>0.5)),fp_bc=("bc",lambda v:100*np.mean(v>T_B)))
print("\n1. single_distinct: false colocalisation vs x = min_k lambda_k^2 (1-r^2)\n",g.round(3).to_string())
# r* boundary for each lambda2 cell
for h in [0.0075,0.015,0.03,0.05]:
    l=lam(1000,h,1); print(f"h2_eqtl {h}: lambda2 {l:.2f}, r* = {math.sqrt(max(0,1-thr/l**2)):.2f}")
# ---- 2. consistency: whole coloc detects iff shared is strongest in each trait
y=d[d.scenario.isin(["A_AB","AB_AC"])&(d.tier=="low")&(d.h2_gwas==0.003)].copy()
def dom(r):
    e = r.eqtl_shared_stronger; g_ = r.gwas_shared_stronger
    ok_e = True if pd.isna(e) else bool(e); ok_g = True if pd.isna(g_) else bool(g_)
    return ok_e and ok_g
y["shared_dominant"]=y.apply(dom,axis=1)
z=y.groupby(["scenario","h2_eqtl","shared_dominant"]).agg(n=("whole","size"),
   whole=("whole",lambda v:100*np.mean(v>T_W)),bc=("bc",lambda v:100*np.mean(v>T_B)),
   whole_pp4_med=("whole","median")).round(3)
print("\n2. detection (%) by whether the shared variant has the largest |z| in every trait with two causals\n",z.to_string())
print("\nshare dominant:",y.groupby(["scenario","h2_eqtl"]).shared_dominant.mean().round(3).to_string())
print("\n--- check: single_distinct by cell (all 12 windows) ---")
x["cell"]=x.h2_gwas.astype(str)+"/"+x.h2_eqtl.astype(str)+"/"+x.tier
print(x.groupby(["cell","bin"],observed=True).agg(n=("whole","size"),r_med=("r","median"),fp_whole=("whole",lambda v:100*np.mean(v>T_W))).round(2).to_string())
print("\n--- single_distinct default lambda2, bc vs whole by tier/bin ---")
xx=x[(x.h2_gwas==0.003)&(x.h2_eqtl==0.015)]
print(xx.groupby(["tier","bin"],observed=True).agg(n=("whole","size"),r=("r","median"),whole=("whole",lambda v:100*np.mean(v>T_W)),bc=("bc",lambda v:100*np.mean(v>T_B))).round(1).to_string())
print(d[(d.scenario=="A_AB")&(d.h2_eqtl==0.015)&(d.tier=="low")&(d.h2_gwas==0.003)].shape)
