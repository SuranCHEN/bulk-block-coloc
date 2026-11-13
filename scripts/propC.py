"""Proposition C checks and figures (docs/11_theory.md §1). Run from the bulk-block-coloc root;
override the per-locus output folder with BULK_OUT and the figure folder with FIG_OUT."""
import os, math, numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
R=os.environ.get("BULK_OUT","../clike_ldblock_coloc_handoff_20260811/outputs/").rstrip("/")+"/bulk_h2ref_v1/"
FIG=os.environ.get("FIG_OUT","docs")
K=["window","scenario","rep","h2_gwas","h2_eqtl","tier"]
T_W,T_B=0.141194270215242,0.144807117782326
BLUE,ORANGE,INK,MUTED,SURF="#2a78d6","#eb6834","#0b0b0b","#52514e","#fcfcfb"
def scores(f):
    d=pd.read_csv(R+f,sep="\t",low_memory=False)
    w=d[d.method=="whole_abf"].set_index(K)["pp4"].rename("whole")
    b=d[(d.method=="V3_gate_conditional")&(d.lead_threshold1==1e-3)&(d.lead_threshold2==1e-3)].copy()
    b["bc"]=b.pp4.where(np.isfinite(b.pp4),b.pp4_whole)
    return pd.concat([w,b.set_index(K)["bc"]],axis=1).reset_index()
loc=pd.concat([pd.read_csv(R+"loci_main.tsv.gz",sep="\t"),pd.read_csv(R+"loci_grid.tsv.gz",sep="\t")])
d=loc.merge(pd.concat([scores("scores_main.tsv.gz"),scores("scores_grid.tsv.gz")]),on=K).drop_duplicates(K)
d=d[(d.tier=="low")&(d.h2_gwas==0.003)]
def dom(r):
    return all(True if pd.isna(v) else bool(v) for v in (r.eqtl_shared_stronger,r.gwas_shared_stronger))
def wilson(k,n,z=1.96):
    p=k/n; den=1+z*z/n; c=(p+z*z/(2*n))/den; h=z*math.sqrt(p*(1-p)/n+z*z/(4*n*n))/den; return c-h,c+h
# ---------- figure 1
y=d[d.scenario.isin(["A_AB","AB_AC"])].copy(); y["dom"]=y.apply(dom,axis=1)
rows=[]
for (sc,h,dm),g in y.groupby(["scenario","h2_eqtl","dom"]):
    for m,t in (("whole",T_W),("bc",T_B)):
        k=int((g[m]>t).sum()); lo,hi=wilson(k,len(g))
        pw=g.groupby("window").apply(lambda x:(x[m]>t).mean())
        rows.append(dict(scenario=sc,h2_eqtl=h,shared_strongest=dm,method=m,n=len(g),rate=100*k/len(g),lo=100*lo,hi=100*hi,win_min=100*pw.min(),win_max=100*pw.max()))
tab=pd.DataFrame(rows); tab.to_csv(os.path.join(FIG,"propC_strata.tsv"),sep="\t",index=False); print(tab.round(1).to_string())
plt.rcParams.update({"font.size":10,"axes.edgecolor":MUTED,"axes.labelcolor":MUTED,"xtick.color":MUTED,"ytick.color":MUTED,"font.family":"DejaVu Sans"})
fig,axs=plt.subplots(1,2,figsize=(10,4.8),sharey=True,facecolor=SURF)
names={"A_AB":"single_double_shared\n(trait 1: A;  trait 2: A, B)","AB_AC":"double_shared1\n(trait 1: A, B;  trait 2: A, C)"}
for ax,sc in zip(axs,["A_AB","AB_AC"]):
    ax.set_facecolor(SURF)
    for m,col,lab in (("whole",BLUE,"coloc.abf"),("bc",ORANGE,"block-coloc")):
        for dm,ls,mk in ((True,"-","o"),(False,"--","s")):
            t=tab[(tab.scenario==sc)&(tab.method==m)&(tab.shared_strongest==dm)].sort_values("h2_eqtl")
            x=np.arange(len(t))
            ax.errorbar(x,t.rate,yerr=[t.rate-t.lo,t.hi-t.rate],color=col,ls=ls,marker=mk,ms=6,lw=2,capsize=0,elinewidth=1,mec=SURF,mew=1.5,zorder=3)
            if dm is False:
                ax.annotate(lab,(x[-1],t.rate.iloc[-1]),xytext=(9,0),textcoords="offset points",ha="left",va="center",color=INK,fontsize=9)
            elif m=="bc":
                ax.annotate("both methods",(x[1],t.rate.iloc[1]),xytext=(8,-12),textcoords="offset points",ha="left",va="center",color=INK,fontsize=9)
    ax.set_xlim(-0.2,3.75); ax.set_xticks(range(4)); ax.set_xticklabels(["0.75%","1.5%","3%","5%"]); ax.set_xlabel("local heritability of trait 2 (eQTL, N = 1,000)")
    ax.set_title(names[sc],fontsize=10,color=INK,loc="left"); ax.set_ylim(0,104); ax.grid(axis="y",color="#e6e5e0",lw=0.8); ax.set_axisbelow(True)
    for s in ("top","right","left"): ax.spines[s].set_visible(False)
    ax.tick_params(length=0)
axs[0].set_ylabel("loci called colocalised (%)")
axs[0].text(0.98,0.04,"solid: the shared variant has the largest |z|\namong the causal variants of each trait\ndashed: it does not",transform=axs[0].transAxes,va="bottom",ha="right",fontsize=8.5,color=MUTED)
from matplotlib.lines import Line2D
fig.legend(handles=[Line2D([],[],color=BLUE,lw=2,marker="o",mec=SURF,label="whole-region coloc.abf"),Line2D([],[],color=ORANGE,lw=2,marker="o",mec=SURF,label="block-coloc")],loc="upper left",frameon=False,ncol=2,fontsize=9,bbox_to_anchor=(0.0,0.955))
fig.suptitle("When the shared variant is not the strongest signal, whole-region coloc.abf stops improving",x=0.01,ha="left",fontsize=11,color=INK,y=1.0)
fig.tight_layout(rect=(0,0,1,0.9)); fig.savefig(os.path.join(FIG,"propC_strata.png"),dpi=160,facecolor=SURF)
# ---------- finite-sample check: single_double_shared, PP4 vs z-gap in trait 2
a=d[d.scenario=="A_AB"].copy()
z=a.causal_abs_z_eqtl.str.split(";",expand=True).astype(float); a["zA"]=z[0]; a["zB"]=z[1]
a["gap"]=0.5*(a.zB**2-a.zA**2)          # > 0: distinct variant B stronger than shared A in trait 2
a["bin"]=pd.cut(a.gap,[-60,-12,-8,-5,-3,-1,1,3,5,8,12,60])
g=a.groupby("bin",observed=True).agg(n=("whole","size"),gap=("gap","median"),pp4_med=("whole","median"),pp4_mean=("whole","mean"))
g["logistic_0.9"]=1/(1+np.exp(0.9*g.gap-math.log(1000)))
print("\nsingle_double_shared, all trait-2 strengths, low LD: whole-region PP4 vs gap = (zB^2 - zA^2)/2 in trait 2\n",g.round(3).to_string())
print("corr(sign(gap<0), whole>T):", np.mean((a.gap<0)==(a.whole>T_W)).round(3))
for h,gg in a.groupby("h2_eqtl"): print(h, "P(A stronger) %.3f"%np.mean(gg.gap<0), " detect|A stronger %.3f  detect|B stronger %.3f"%((gg.whole[gg.gap<0]>T_W).mean(),(gg.whole[gg.gap>0]>T_W).mean()))
g.to_csv(os.path.join(FIG,"propC_gap.tsv"),sep="\t")

q=a.groupby("bin",observed=True)["whole"].quantile([.25,.75]).unstack()
fig,ax=plt.subplots(figsize=(6.4,4.2),facecolor=SURF); ax.set_facecolor(SURF)
ax.axvline(0,color="#d0cfc8",lw=1); ax.axhline(0.5,color="#d0cfc8",lw=1)
ax.vlines(g.gap,q[.25],q[.75],color=BLUE,lw=1,alpha=.6)
ax.plot(g.gap,g.pp4_med,color=BLUE,lw=2,marker="o",ms=6,mec=SURF,mew=1.5)
ax.set_xlabel("(z_B² − z_A²) / 2 in trait 2   (A shared, B specific to trait 2)"); ax.set_ylabel("whole-region coloc.abf PP4 (median, IQR)")
ax.text(-14,0.06,"shared variant stronger",color=MUTED,fontsize=8.5); ax.text(14.5,0.9,"specific variant stronger",color=MUTED,fontsize=8.5,ha="right")
ax.set_title("PP4 falls as the trait-2-specific signal overtakes the shared one",fontsize=10,color=INK,loc="left")
ax.set_ylim(-0.02,1.02); ax.grid(axis="y",color="#e6e5e0",lw=0.8); ax.set_axisbelow(True); ax.tick_params(length=0)
for sp in ("top","right","left"): ax.spines[sp].set_visible(False)
fig.tight_layout(); fig.savefig(os.path.join(FIG,"propC_gap.png"),dpi=160,facecolor=SURF)
