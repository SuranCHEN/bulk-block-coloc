"""Conditional vs marginal lead selection inside block-coloc, same loci, reference LD (stage 1, stage 2, normal effects).\nRun from the bulk-block-coloc root; override the per-locus output folder with BULK_OUT."""
import os, math, numpy as np, pandas as pd
from scipy.stats import binomtest, t as tdist
ROOT=os.environ.get("BULK_OUT","../clike_ldblock_coloc_handoff_20260811/outputs").rstrip("/")+"/"
POS=["A_A","A_AB","AB_AB","AB_AC"]; NEG=["A_B","AB_CD"]
K=["window","scenario","rep","h2_gwas","h2_eqtl","tier"]
def cut(v): v=np.sort(v)[::-1]; return v[int(math.floor(0.05*len(v)))]
def load(f):
    d=pd.read_csv(ROOT+f,sep="\t",low_memory=False); return d
def table(d, keys, label):
    """keys: dict name -> (method, thr1, thr2, fallback)"""
    out={}
    for name,(m,t1,t2,fb) in keys.items():
        x=d[d.method==m]
        if t1 is not None: x=x[(x.lead_threshold1==t1)&(x.lead_threshold2==t2)]
        pp=x.pp4.where(np.isfinite(x.pp4),x.pp4_whole) if fb else x.pp4
        s=pd.Series(np.nan_to_num(pp.values.astype(float),nan=-1.0),index=pd.MultiIndex.from_frame(x[K]))
        out[name]=s
    df=pd.DataFrame(out).reset_index(); df["set"]=df.window.map(d.drop_duplicates("window").set_index("window")["set"])
    cal=df[(df.set=="calibration")&df.scenario.isin(NEG)]; tp=df[(df.set=="test")&df.scenario.isin(POS)]; tn=df[(df.set=="test")&df.scenario.isin(NEG)]
    thr={n:cut(cal[n].values) for n in keys}
    print(f"\n== {label}: test power / FPR at matched 5% (cut-off)")
    for n in keys:
        byS=" ".join(f"{sc}:{100*np.mean(tp[tp.scenario==sc][n]>thr[n]):.1f}" for sc in POS)
        go=df[(df.set=='test')&(df.scenario=='gwas_only')]
        print(f"  {n:34s} power {100*np.mean(tp[n]>thr[n]):5.1f}  FPR {100*np.mean(tn[n]>thr[n]):4.1f}  gwas_only {100*np.mean(go[n]>thr[n]):4.1f}  cut {thr[n]:.3f} | {byS}")
    return df,thr,tp
def paired(tp,thr,a,b):
    da=tp[a].values>thr[a]; db=tp[b].values>thr[b]; n=len(da)
    b10=int((da&~db).sum()); b01=int((~da&db).sum()); p=binomtest(b10,b10+b01,.5).pvalue if b10+b01 else 1
    per=np.array([100*(da[tp.window.values==w].mean()-db[tp.window.values==w].mean()) for w in np.sort(tp.window.unique())])
    q=tdist.ppf(.975,len(per)-1); se=per.std(ddof=1)/math.sqrt(len(per))
    print(f"  {a} − {b}: {100*(b10-b01)/n:+.1f} (McNemar p={p:.2g}); window mean {per.mean():+.1f} ({per.mean()-q*se:+.1f}, {per.mean()+q*se:+.1f}); windows>0 {int((per>0).sum())}/{len(per)}; per-window {', '.join(f'{x:+.1f}' for x in per)}")
for lab,f in [("FIXED effects, default (stage 1)","bulk_h2ref_v1/scores_main.tsv.gz"),("NORMAL effects, default","bulk_h2ref_v1_normal/scores_main.tsv.gz"),("NORMAL effects, h2_eqtl 5%","bulk_h2ref_v1_normal_g0.003_e0.05/scores_main.tsv.gz")]:
    d=load(f); keys={"whole_abf":("whole_abf",None,None,False),"susie+fb":("susie_fallback",None,None,False),"signals_cond":("signals_cond",None,None,False)}
    for t in (1e-3,1e-4,1e-5):
        keys[f"cond_fb@{t:g}"]=("V3_gate_conditional",t,t,True); keys[f"marg_fb@{t:g}"]=("V2_gate_marginal",t,t,True)
    keys["blockany_fb@0.001"]=("V1_block_any",1e-3,1e-3,True)
    df,thr,tp=table(d,keys,lab)
    for t in ("0.001","0.0001","1e-05"): paired(tp,thr,f"cond_fb@{t}",f"marg_fb@{t}")
    paired(tp,thr,"cond_fb@0.001","marg_fb@1e-05"); paired(tp,thr,"marg_fb@0.001","whole_abf"); paired(tp,thr,"cond_fb@0.001","susie+fb")
    if "FIXED" in lab:
        x=d[(d.method=="V3_gate_conditional")&(d.lead_threshold1==1e-3)&(d.lead_threshold2==1e-3)]
        print("  lead counts (all windows): scenario | trait1 marginal>=2 / cond>=2 | trait2 marginal>=2 / cond>=2 | true K1,K2")
        for sc,g in x.groupby("scenario"):
            print(f"    {sc:10s} {100*np.mean(g.n_leads1_marginal>=2):5.1f} / {100*np.mean(g.n_leads1>=2):5.1f} | {100*np.mean(g.n_leads2_marginal>=2):5.1f} / {100*np.mean(g.n_leads2>=2):5.1f} | {g.n_causal_gwas.iloc[0]},{g.n_causal_eqtl.iloc[0]}")
# stage 2: selected versions only (cond@1e-3, marginal@1e-5), stage-1 cut-offs
d1=load("bulk_h2ref_v1/scores_main.tsv.gz"); g=load("bulk_h2ref_v1/scores_grid.tsv.gz")
keys={"whole_abf":("whole_abf",None,None,False),"susie+fb":("susie_fallback",None,None,False),"susie_raw":("susie_raw",None,None,False),"signals_cond":("signals_cond",None,None,False),"cond_fb@0.001":("V3_gate_conditional",1e-3,1e-3,True),"marg_fb@1e-05":("V2_gate_marginal",1e-5,1e-5,True)}
df1,thr,_=table(d1,keys,"stage-1 cut-offs for stage 2")
print("\n== Stage 2 (stage-1 cut-offs): power (FPR)  [cell: whole | susie+fb | susie_raw power, no-call% | signals_cond | cond_fb | marg_fb]")
out={}
for name,(m,t1,t2,fb) in keys.items():
    x=g[g.method==m]
    if t1 is not None: x=x[(x.lead_threshold1==t1)&(x.lead_threshold2==t2)]
    pp=x.pp4.where(np.isfinite(x.pp4),x.pp4_whole) if fb else x.pp4
    out[name]=pd.Series(np.nan_to_num(pp.values.astype(float),nan=-1.0),index=pd.MultiIndex.from_frame(x[K]))
G=pd.DataFrame(out).reset_index(); G["set"]=G.window.map(g.drop_duplicates("window").set_index("window")["set"])
for (hg,he,tr),x in G[G.set=="test"].groupby(["h2_gwas","h2_eqtl","tier"]):
    tp=x[x.scenario.isin(POS)]; tn=x[x.scenario.isin(NEG)]
    s=" | ".join(f"{100*np.mean(tp[n]>thr[n]):.1f} ({100*np.mean(tn[n]>thr[n]):.1f})" for n in ["whole_abf","susie+fb","signals_cond","cond_fb@0.001","marg_fb@1e-05"])
    print(f"  g{hg} e{he} {tr}: {s} | susie_raw no-call {100*np.mean(tp['susie_raw']<0):.0f}%")
    if tr=="low": paired(tp,thr,"cond_fb@0.001","marg_fb@1e-05"); paired(tp,thr,"cond_fb@0.001","susie+fb")
