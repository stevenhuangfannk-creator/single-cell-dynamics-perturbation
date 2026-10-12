"""Independent redraws of saved Knk results; no network reconstruction or model fitting."""
from pathlib import Path
import argparse
import json
import textwrap
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
import networkx as nx

plt.rcParams.update({"font.family":"sans-serif", "font.sans-serif":["Arial","DejaVu Sans"],
    "font.size":8, "axes.titlesize":10, "axes.labelsize":8, "xtick.labelsize":7,
    "ytick.labelsize":7, "legend.fontsize":7, "svg.fonttype":"none", "pdf.fonttype":42,
    "axes.spines.top":False, "axes.spines.right":False, "figure.facecolor":"white",
    "axes.facecolor":"white", "savefig.facecolor":"white"})
BLUE="#2584C0"; RED="#C85B5B"; GRAY="#838D96"; ORANGE="#E69F00"
TITLES = {1:"WT调控网络（显示子图）",2:"SPI1虚拟敲除网络",3:"WT与KO网络比较",4:"差异调控扰动排序",5:"扰动统计显著性",6:"SPI1中心调控网络",7:"GO/KEGG/Reactome富集",8:"候选KO网络扰动热图",9:"单基因与联合KO响应比较",10:"随机种子和细胞抽样稳定性",11:"候选扰动图（内部参考）"}
MEANINGS = {1:"显示WT的强权重PCR网络关联，不能据边符号认定真实激活或抑制。",2:"与WT共用布局，SPI1行出边归零，其他边保持一致。",3:"核对SPI1出边被移除且其他行不变；网络结构变化由算法定义。",4:"非target基因按WT–KO manifold distance排序，颜色仅为官方另算的方向预测。",5:"距离与BH调整p值的关系；按构建必然变化的SPI1 target已排除。",6:"SPI1出边中绝对权重最大的16个邻居；边权是模型关联而非实验证实调控。",7:"非target BH<0.05响应基因在600个测量网络基因背景上的超几何富集。",8:"9个候选分别KO的官方heat distance，颜色是行内标准化log1p距离，无表达方向。",9:"在相同WT上比较SPI1单KO与SPI1+CEBPB联合KO的manifold distance。",10:"600gene的6组种子/抽样与1200gene/10nets的3组敏感性均排名不稳定，不能据主seed宣称可靠机制。",11:"官方perturbationMap坐标；与同一模型SPI1 signed heat-KO内部参考的相似性，非疾病签名，只有9候选。"}

def network_view(WT):
    neighbors = WT.loc["SPI1"].abs().drop("SPI1").nlargest(15).index.tolist()
    strength = WT.abs().sum(axis=0)+WT.abs().sum(axis=1)
    hubs = [g for g in strength.nlargest(40).index if g not in neighbors and g!="SPI1"][:16]
    nodes = ["SPI1"]+neighbors+hubs
    sub=WT.loc[nodes,nodes]; pairs=[]
    for a in nodes:
        for b in nodes:
            if a!=b and sub.loc[a,b]!=0: pairs.append((a,b,float(sub.loc[a,b])))
    best=sorted(pairs,key=lambda v:abs(v[2]),reverse=True)[:80]
    chosen={(a,b):(a,b,w) for a,b,w in best}
    for b in neighbors:
        chosen[("SPI1",b)]=("SPI1",b,float(WT.loc["SPI1",b]))
    graph=nx.DiGraph(); graph.add_nodes_from(nodes)
    for a,b,w in chosen.values(): graph.add_edge(a,b,weight=w)
    positions=nx.circular_layout(graph)
    return graph,positions

def draw_network(ax, graph, positions, ko=False, radial_labels=False):
    nodes=list(graph.nodes)
    nx.draw_networkx_nodes(graph,positions,ax=ax,node_size=[180 if g=="SPI1" else 70 for g in nodes],node_color=[ORANGE if g=="SPI1" else "#DCECF6" for g in nodes],edgecolors="white",linewidths=0.5)
    edges=[(a,b) for a,b in graph.edges if not (ko and a=="SPI1")]
    nx.draw_networkx_edges(graph,positions,ax=ax,edgelist=edges,edge_color=[RED if graph[a][b]["weight"]>0 else BLUE for a,b in edges],width=0.55,alpha=0.55,arrows=True,arrowsize=5,node_size=85)
    if radial_labels:
        for g,p in positions.items():
            ax.text(p[0]*1.1,p[1]*1.1,g,ha="left" if p[0]>0.1 else "right" if p[0]<-0.1 else "center",va="center",fontsize=7)
    else:
        nx.draw_networkx_labels(graph,positions,ax=ax,font_size=7,font_family="Arial")
    ax.margins(0.20); ax.axis("off")

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--module",type=Path,default=Path("modules/01_scTenifoldKnk"))
    ap.add_argument("--figure",default="all",help="all or 01 through 11")
    ap.add_argument("--dpi",type=int,default=600)
    args=ap.parse_args(); module=args.module; results=module/"results"; figures=module/"figures"
    figures.mkdir(parents=True,exist_ok=True)
    WT=pd.read_csv(module/"objects/WT.csv",index_col=0)
    KO=pd.read_csv(module/"objects/KO.csv",index_col=0)
    dr=pd.read_csv(results/"single_SPI1.csv").set_index("gene")
    multi=pd.read_csv(results/"multi_SPI1_CEBPB.csv").set_index("gene")
    distances=pd.read_csv(results/"candidate_screen_distance.csv",index_col=0)
    graph,pos=network_view(WT)
    requested=list(range(1,12)) if args.figure=="all" else [int(args.figure)]
    entries=[]
    for i in requested:
        if i in (1,2):
            fig,ax=plt.subplots(figsize=(5.5,4.5)); draw_network(ax,graph,pos,ko=i==2,radial_labels=True)
            ax.set_title("WT GRN" if i==1 else "Virtual SPI1 KO GRN",loc="left")
            fig.text(0.04,0.02,"32 nodes; 80 strongest displayed edges + 15 SPI1 edges | Positive / negative PCR weights",fontsize=7)
        elif i==3:
            fig,axes=plt.subplots(1,2,figsize=(7.2,3.4),layout="constrained")
            values=[int((WT.loc["SPI1"]!=0).sum()),int((KO.loc["SPI1"]!=0).sum())]
            axes[0].bar(["WT","SPI1 KO"],values,color=[BLUE,ORANGE]); axes[0].set_ylabel("SPI1 outgoing nonzero edges")
            for j,v in enumerate(values): axes[0].text(j,v+8,str(v),ha="center",fontsize=8)
            axes[0].set_ylim(0,max(values)*1.15)
            delta=WT-KO; removed=delta.loc["SPI1"].abs().nlargest(15)
            axes[1].barh(removed.index[::-1],removed.values[::-1],color=ORANGE); axes[1].set_xlabel("Removed absolute edge weight")
            unchanged=bool(np.array_equal(WT.drop(index="SPI1").values,KO.drop(index="SPI1").values))
            fig.suptitle(f"WT vs KO | Other regulator rows identical: {unchanged}")
        elif i==4:
            top=dr.drop(index="SPI1").sort_values("distance",ascending=False).head(20).iloc[::-1]
            fig,ax=plt.subplots(figsize=(5.5,4.9),layout="constrained")
            ax.barh(top.index,top.distance,color=[RED if d=="up" else BLUE for d in top.direction]); ax.set_xlabel("Manifold perturbation distance")
            ax.set_title("Differential regulation ranking | SPI1 excluded",loc="left")
            ax.plot([],[],color=RED,label="Predicted up"); ax.plot([],[],color=BLUE,label="Predicted down"); ax.legend(loc="lower right")
        elif i==5:
            fig,ax=plt.subplots(figsize=(5.5,4.2),layout="constrained")
            data=dr.drop(index="SPI1"); y=-np.log10(np.maximum(data["p.adj"].to_numpy(),1e-300))
            ax.scatter(data.distance,y,s=12,c=np.where(data["p.adj"]<0.05,RED,GRAY),alpha=0.75)
            ax.axhline(-np.log10(0.05),linestyle="--",color="black",lw=0.6)
            for gene in data.nlargest(4,"distance").index: ax.annotate(gene,(data.loc[gene,"distance"],-np.log10(max(data.loc[gene,"p.adj"],1e-300))),xytext=(3,3),textcoords="offset points",fontsize=7)
            ax.set_xlabel("Manifold perturbation distance"); ax.set_ylabel("-log10(BH-adjusted p)")
            ax.set_title("Perturbation significance | SPI1 target excluded",loc="left")
        elif i==6:
            fig,ax=plt.subplots(figsize=(5.7,4.8)); neighbors=WT.loc["SPI1"].abs().drop("SPI1").nlargest(16).index
            center=nx.DiGraph(); center.add_node("SPI1")
            for g in neighbors: center.add_edge("SPI1",g,weight=float(WT.loc["SPI1",g]))
            p={"SPI1":np.array([0.,0.])}; p.update({g:np.array([np.cos(t),np.sin(t)]) for g,t in zip(neighbors,np.linspace(0,2*np.pi,len(neighbors),endpoint=False))})
            draw_network(ax,center,p); ax.set_title("SPI1-centered WT outgoing network",loc="left")
            fig.text(0.05,0.02,"Top 16 absolute outgoing weights | Signed PCR association, not causal edge direction",fontsize=7)
        elif i==7:
            enrich=pd.read_csv(results/"enrichment_all_terms.csv"); enrich=enrich[enrich.contrast=="single_SPI1"]
            fig,axes=plt.subplots(1,3,figsize=(10.8,4.6),layout="constrained")
            for ax,(lib,data) in zip(axes,enrich.groupby("library",sort=False)):
                top=data.sort_values(["p_adj","p_value"]).head(6).iloc[::-1]
                labels=[textwrap.fill(str(t),31) for t in top.term]
                x=-np.log10(np.maximum(top.p_adj.to_numpy(),1e-300))
                ax.barh(np.arange(len(top)),x,color=np.where(top.p_adj<0.05,BLUE,GRAY)); ax.set_yticks(np.arange(len(top)),labels)
                ax.axvline(-np.log10(0.05),color="black",ls="--",lw=0.6); ax.set_xlabel("-log10(BH-adjusted p)")
                ax.set_title(lib.replace("_"," "),fontsize=9)
            fig.suptitle("SPI1 KO enrichment | 14 non-target hits | 600-gene measured background",fontsize=10)
        elif i==8:
            # Row standardization is explicitly a display transform of unsigned distances.
            selected=distances.mean(axis=0).sort_values(ascending=False).head(30).index
            x=np.log1p(distances.loc[:,selected].to_numpy()); x=(x-x.mean(axis=1,keepdims=True))/np.maximum(x.std(axis=1,keepdims=True),1e-12)
            fig,ax=plt.subplots(figsize=(8.3,3.8),layout="constrained")
            im=ax.imshow(x,aspect="auto",cmap="viridis"); ax.set_xticks(range(len(selected)),selected,rotation=70,ha="right",rotation_mode="anchor"); ax.set_yticks(range(len(distances)),distances.index)
            ax.set_title("Candidate KO heat distance | 30 highest mean-distance genes",loc="left")
            fig.colorbar(im,ax=ax,label="Row-standardized log1p(unsigned distance)",shrink=0.8)
        elif i==9:
            genes=dr.index.difference(["SPI1","CEBPB"]); x=dr.loc[genes,"distance"]; y=multi.loc[genes,"distance"]
            fig,axes=plt.subplots(1,2,figsize=(7.6,3.8),layout="constrained")
            axes[0].scatter(x,y,s=9,color=BLUE,alpha=0.6); axes[0].set_xlabel("SPI1 KO distance"); axes[0].set_ylabel("SPI1 + CEBPB KO distance")
            for gene,offset in zip(y.nlargest(3).index,[(-55,4),(6,8),(6,-12)]):
                axes[0].annotate(gene,(x[gene],y[gene]),xytext=offset,textcoords="offset points",fontsize=7)
            sel=(x+y).nlargest(12).index; a=np.arange(len(sel))
            axes[1].barh(a-0.18,x[sel],height=0.35,color=BLUE,label="SPI1 KO"); axes[1].barh(a+0.18,y[sel],height=0.35,color=ORANGE,label="Joint KO")
            axes[1].set_yticks(a,sel); axes[1].invert_yaxis(); axes[1].set_xlabel("Manifold distance"); axes[1].legend(loc="lower right")
            fig.suptitle("Single vs joint KO | Same WT; both target genes excluded")
        elif i==10:
            stability=pd.read_csv(results/"stability_summary.csv")
            fig,grid=plt.subplots(2,2,figsize=(7.3,6.7),layout="constrained")
            axes=grid[0]
            for n,d in stability.groupby("cells"):
                axes[0].plot(d.seed,d.spearman,"o-",label=f"{n} cells/network",color=BLUE if n==400 else ORANGE)
                axes[1].plot(d.seed,d.top30_jaccard,"o-",label=f"{n} cells/network",color=BLUE if n==400 else ORANGE)
            for ax in axes: ax.set_xticks([1,2,3]); ax.set_xlabel("Random seed"); ax.set_ylim(-0.05,1.05)
            axes[0].set_ylabel("Spearman to seed 1, 400 cells"); axes[1].set_ylabel("Top-30 Jaccard to reference")
            axes[0].legend(loc="upper right")
            sensitivity=pd.read_csv(results/"sensitivity1200/parameter_sensitivity_summary.csv")
            grid[1,0].bar(sensitivity.seed,sensitivity.spearman_to_600,color=GRAY)
            grid[1,0].set_xticks([1,2,3]); grid[1,0].set_xlabel("1200-gene / 10-network seed")
            grid[1,0].set_ylabel("Spearman to 600-gene reference"); grid[1,0].set_ylim(0,1)
            corr=pd.read_csv(results/"sensitivity1200/within_1200_seed_spearman.csv",index_col=0)
            im=grid[1,1].imshow(corr,cmap="viridis",vmin=0,vmax=1)
            grid[1,1].set_xticks(range(3),corr.columns); grid[1,1].set_yticks(range(3),corr.index)
            for a in range(3):
                for b in range(3): grid[1,1].text(b,a,f"{corr.iloc[a,b]:.2f}",ha="center",va="center",color="black" if a==b else "white",fontsize=8)
            grid[1,1].set_title("Within 1200-gene panel: Spearman",fontsize=9)
            fig.colorbar(im,ax=grid[1,1],shrink=0.75)
            fig.suptitle("Ranking stability fails | 6 primary runs + 3 broader-panel runs retained")
        elif i==11:
            pm=pd.read_csv(results/"candidate_perturbation_map.csv")
            fig,ax=plt.subplots(figsize=(5.8,4.3),layout="constrained")
            im=ax.scatter(pm.x,pm.y,c=pm.similarity,cmap="viridis",s=65,edgecolor="white",lw=0.6)
            offsets={"CEBPD":(-6,13),"FOS":(-22,-12),"JUN":(10,-7)}
            for row in pm.itertuples(): ax.annotate(row.gene,(row.x,row.y),xytext=offsets.get(row.gene,(4,4)),textcoords="offset points",fontsize=8)
            ax.margins(0.23); ax.set_xlabel("Cosine similarity to internal SPI1 heat-KO profile"); ax.set_ylabel("Residual profile axis (official SVD)")
            ax.set_title("Candidate perturbation map | 9 targets; internal reference",loc="left")
            fig.colorbar(im,ax=ax,label="Internal-reference similarity",shrink=0.75)
        else: raise ValueError("Figure must be 01 through 11")
        base=f"Knk_Figure{i:02d}"
        for ext in ("png","pdf","svg"):
            fig.savefig(figures/f"{base}.{ext}",dpi=args.dpi,bbox_inches="tight",facecolor="white")
        plt.close(fig)
        root="modules/01_scTenifoldKnk"
        entries.append(dict(id=f"Knk-{i:02d}",title_cn=TITLES[i],meaning=MEANINGS[i],dataset="10x PBMC3k, 480 CD14+ monocytes, 600 selected genes"+("; 1200-gene sensitivity" if i==10 else ""),reference="https://doi.org/10.1016/j.patter.2022.100434",analysis_script=f"{root}/scripts/"+("enrich.py" if i==7 else "run_pbmc.R"),plot_script=f"{root}/scripts/plot_figures.py",redraw_command=f"python {root}/scripts/plot_figures.py --figure {i:02d}",png=f"{root}/figures/{base}.png",pdf=f"{root}/figures/{base}.pdf",svg=f"{root}/figures/{base}.svg",status="RUN_OK / PARTIAL",validity="真实计算输出；seed/subsample稳定性差，探索性候选结果，未实验验证。"))
    manifest_path=figures/"figure_manifest.json"
    if manifest_path.exists() and args.figure!="all":
        old=json.loads(manifest_path.read_text(encoding="utf-8")); entries=[e for e in old if e["id"] not in {v["id"] for v in entries}]+entries
    manifest_path.write_text(json.dumps(sorted(entries,key=lambda v:v["id"]),ensure_ascii=False,indent=2),encoding="utf-8")
    print(f"Saved {len(requested)} actual figures at {args.dpi} DPI plus PDF/SVG.")

if __name__=="__main__": main()
