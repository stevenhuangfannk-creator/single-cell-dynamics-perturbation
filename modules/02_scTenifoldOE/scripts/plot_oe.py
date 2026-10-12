"""Redraw individual experimental OE figures from saved network/score tables."""
import argparse
import json
import re
import textwrap
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

mpl.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Arial", "DejaVu Sans"],
                     "font.size": 7, "axes.titlesize": 9, "axes.labelsize": 7,
                     "svg.fonttype": "none", "pdf.fonttype": 42,
                     "figure.facecolor": "white", "axes.facecolor": "white",
                     "axes.spines.top": False, "axes.spines.right": False,
                     "legend.frameon": False})
MODULE = Path(__file__).resolve().parents[1]
SOURCE = "Experimental: virt_overexpression 0.1.0 (virchow-zhang) + official R manifold/DR"
COLORS = {"WT": "#989898", "KO": "#627CBB", "OE_1.5": "#83B6A4", "OE_2": "#347F70", "OE_3": "#174F43"}


def read_dr(path, target):
    table = pd.read_csv(path).set_index("gene")
    return table.loc[table.index != target].copy()


def draw_network(ax, matrix, genes, positions, target, title, limit=70, edge_scale=None):
    sub = matrix.loc[genes, genes]
    edges = [(a, b, float(sub.loc[a, b])) for a in genes for b in genes if a != b and sub.loc[a, b] != 0]
    ranked = sorted(edges, key=lambda x: abs(x[2]), reverse=True)
    centered = [e for e in ranked if e[0] == target][:min(18, limit)]
    edges = centered + [e for e in ranked if e not in centered][:limit - len(centered)]
    graph = nx.DiGraph()
    graph.add_nodes_from(genes)
    graph.add_weighted_edges_from(edges)
    maxw = edge_scale or max((abs(w) for _, _, w in edges), default=1)
    nx.draw_networkx_nodes(graph, positions, ax=ax, node_size=[160 if g == target else 75 for g in genes],
                           node_color=["#C47C37" if g == target else "#D7E7E2" for g in genes], linewidths=.3, edgecolors="#456257")
    nx.draw_networkx_edges(graph, positions, ax=ax, edge_color=["#568C7A" if w > 0 else "#876BA2" for _, _, w in edges],
                           width=[.4 + 1.2 * abs(w) / maxw for _, _, w in edges], alpha=.6, arrows=True,
                           arrowsize=8, node_size=120, connectionstyle="arc3,rad=0.08")
    label_positions = {g: (positions[g] * 1.13 if g != target else positions[g]) for g in genes}
    nx.draw_networkx_labels(graph, label_positions, ax=ax, font_size=7, font_family="Arial")
    ax.set_xlim(-1.32, 1.32)
    ax.set_ylim(-1.32, 1.32)
    ax.set_title(title)
    ax.axis("off")


def prepare_layout(wt, target, output):
    neighbors = wt.loc[target].abs().sort_values(ascending=False)
    genes = [target] + [g for g in neighbors.index if g != target and neighbors[g] > 0][:18]
    # Target centered circular graph layout is a display convention, not cellular coordinates.
    angles = np.linspace(0, 2 * np.pi, len(genes) - 1, endpoint=False)
    pos = {target: np.array([0., 0.])}
    pos.update({g: np.array([np.cos(a), np.sin(a)]) for g, a in zip(genes[1:], angles)})
    pd.DataFrame.from_dict(pos, orient="index", columns=["x", "y"]).to_csv(output / "results" / "network_layout.csv")
    return genes, pos


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--module", type=Path, default=MODULE)
    parser.add_argument("--figure", default="all", help="all or a two-digit figure number")
    parser.add_argument("--dpi", type=int, default=600)
    args = parser.parse_args()
    root = args.module
    meta = json.loads((root / "reports" / "run_metadata.json").read_text(encoding="utf-8-sig"))
    target = meta["target"]
    wt = pd.read_csv(root / "objects" / "WT.csv", index_col=0)
    oe = pd.read_csv(root / "objects" / "OE_2.csv", index_col=0)
    ko = pd.read_csv(root / "objects" / "KO.csv", index_col=0)
    scores = {key: read_dr(root / "results" / f"{key}_dr.csv", target) for key in ["KO", "OE_1", "OE_1.5", "OE_2", "OE_3"]}
    metrics = pd.read_csv(root / "results" / "sensitivity_metrics.csv")
    genes, positions = prepare_layout(wt, target, root)
    manifest = []

    def save(fig, number, title, files, meaning):
        fig.subplots_adjust(bottom=.32 if number == "02" else .16)
        fig.text(.5, .02, SOURCE, ha="center", fontsize=6)
        stem = f"OE_Figure{number}_{title}"
        for ext in ["png", "pdf", "svg"]:
            fig.savefig(root / "figures" / f"{stem}.{ext}", dpi=args.dpi, facecolor="white", bbox_inches="tight")
        plt.close(fig)
        titles_cn = {"01": "虚拟过表达网络", "02": "WT与OE网络差异", "03": "OE扰动响应排序", "04": "SPI1中心调控网络",
                     "05": "过表达强度敏感性", "06": "KO与OE扰动比较", "07": "KO/OE共享与特异响应", "08": "KO/OE通路富集",
                     "09": "不同强度响应热图", "10": "三状态调控网络景观"}
        rel = "modules/02_scTenifoldOE/"
        manifest.append({"number": number, "id": "OE_" + number, "title_cn": titles_cn[number], "stem": stem, "status": "RUN_OK", "method": SOURCE,
                         "dataset": "10x PBMC3k; CD14+ Monocytes; 600 genes × 480 cells; shared Knk WT; SPI1",
                         "source_tables": files, "interpretation": meaning,
                         "meaning": meaning, "reference": "https://github.com/virchow-zhang/virt_overexpression/tree/860d0ea59dc0aeb564a286efae00f562ce84e138",
                         "analysis_script": rel + ("scripts/enrich_compare.py" if number == "08" else "scripts/run_oe.py"),
                         "plot_script": rel + "scripts/plot_oe.py",
                         "png": rel + "figures/" + stem + ".png", "pdf": rel + "figures/" + stem + ".pdf", "svg": rel + "figures/" + stem + ".svg",
                         "validity": "Experimental网络模拟；同一WT，无实验生物学重复；distance不含表达方向，非真实mRNA/CRISPR预测。"})

    requested = lambda num: args.figure in ["all", num]
    if requested("01"):
        fig, ax = plt.subplots(figsize=(7.2, 4.6))
        draw_network(ax, oe, genes, positions, target, f"Virtual OE network: {target}, edge factor 2", 70)
        save(fig, "01", "Virtual_OE_Network", ["objects/OE_2.csv", "results/network_layout.csv"], "固定布局，靶点边优先并补最强关联边，共最多70边；颜色表示PC回归权重符号，非真实因果激活/抑制。")
    if requested("02"):
        fig, ax = plt.subplots(figsize=(7.2, 3.8))
        values = (oe - wt).loc[genes, genes].to_numpy()
        limit = max(np.abs(values).max(), 1e-12)
        im = ax.imshow(values, cmap="PuOr_r", vmin=-limit, vmax=limit, aspect="auto")
        ax.set_xticks(range(len(genes)), genes, rotation=90, rotation_mode="anchor", ha="right", fontsize=7)
        ax.set_yticks(range(len(genes)), genes, fontsize=7)
        ax.set(xlabel="Target gene", ylabel="Source gene", title="WT to OE: adjacency weight difference")
        fig.colorbar(im, ax=ax, label="OE minus WT edge weight", shrink=.75)
        save(fig, "02", "Network_Difference", ["objects/WT.csv", "objects/OE_2.csv"], "仅SPI1行的正出边增强；负边及其他基因边权原样保留。")
    if requested("03"):
        top = scores["OE_2"].sort_values("distance", ascending=False).head(20).iloc[::-1]
        fig, ax = plt.subplots(figsize=(7.2, 4.5))
        ax.barh(top.index, top.distance, color=np.where(top["p.adj"] < .05, "#347F70", "#AAB6B2"))
        ax.set(xlabel="Manifold distance (unsigned)", title=f"OE response ranking: {target}, edge factor 2")
        save(fig, "03", "Perturbation_Ranking", ["results/OE_2_dr.csv"], "按无方向manifold距离排序，排除靶点；绿色FDR<0.05，灰色非显著。")
    if requested("04"):
        fig, ax = plt.subplots(figsize=(7.2, 4.6))
        target_only = oe.copy()
        target_only.loc[target_only.index != target, :] = 0
        draw_network(ax, target_only, genes, positions, target, f"{target}-centered OE outgoing regulation", 40)
        save(fig, "04", "Regulatory_Network", ["objects/OE_2.csv", "results/network_layout.csv"], "显示预先指定靶点的最强18个出边邻居；模型权重是回归关联。")
    if requested("05"):
        sub = metrics[metrics.condition.str.startswith("OE")].copy()
        sub["factor"] = sub.condition.str.replace("OE_", "", regex=False).astype(float)
        fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.1))
        axes[0].plot(sub.factor, sub.delta_frobenius, "o-", color="#347F70")
        axes[1].plot(sub.factor, sub.median_distance, "o-", color="#347F70")
        axes[0].set(ylabel="Adjacency difference Frobenius norm", title="Applied network dose")
        axes[1].set(ylabel="Median manifold distance", title="Downstream network response")
        for ax in axes:
            ax.set(xlabel="Positive target-edge multiplier", xticks=[1, 1.5, 2, 3])
        fig.tight_layout()
        save(fig, "05", "Strength_Sensitivity", ["results/sensitivity_metrics.csv"], "同一WT下网络剂量与响应的敏感性；单次WT，无重复误差条，不等于mRNA倍数。")
    if requested("06"):
        pair = scores["KO"][["distance"]].rename(columns={"distance": "KO"}).join(scores["OE_2"][["distance"]].rename(columns={"distance": "OE"}))
        rho = spearmanr(pair.KO, pair.OE).statistic
        pair.to_csv(root / "results" / "KO_OE_comparison.csv")
        fig, ax = plt.subplots(figsize=(5.2, 4.2))
        ax.scatter(pair.KO, pair.OE, s=9, alpha=.7, color="#347F70", linewidths=0)
        for gene, offset in zip(pair.OE.nlargest(3).index, [(4, 4), (4, 12), (4, -8)]):
            ax.annotate(gene, (pair.loc[gene, "KO"], pair.loc[gene, "OE"]), xytext=offset, textcoords="offset points", fontsize=7)
        ax.set(xlabel="KO manifold distance", ylabel="OE manifold distance", title=f"KO versus OE 2: Spearman rho = {rho:.3f}")
        save(fig, "06", "KO_OE_Comparison", ["results/KO_OE_comparison.csv"], "相同WT、靶点、MA维数和统计规则；比较无方向距离，相关性不证明相反调控。")
    if requested("07"):
        sig = {k: set(v.index[v["p.adj"] < .05]) for k, v in scores.items() if k in ["KO", "OE_2"]}
        groups = {"Shared": sig["KO"] & sig["OE_2"], "KO specific": sig["KO"] - sig["OE_2"], "OE specific": sig["OE_2"] - sig["KO"]}
        pd.DataFrame([{"group": k, "gene": g} for k, v in groups.items() for g in sorted(v)], columns=["group", "gene"]).to_csv(root / "results" / "shared_specific_responses.csv", index=False)
        fig, ax = plt.subplots(figsize=(5.5, 3.3))
        ax.bar(list(groups), [len(v) for v in groups.values()], color=["#7C8C86", "#627CBB", "#347F70"])
        ax.set(ylabel="Genes (target excluded)", title="FDR < 0.05: shared and specific responses")
        save(fig, "07", "Shared_Specific_Responses", ["results/shared_specific_responses.csv", "results/KO_dr.csv", "results/OE_2_dr.csv"], "共享/特异取决于当前WT与FDR阈值，空集合如实显示，未显著不等于无生物学作用。")
    if requested("08"):
        enrichment = root / "results" / "pathway_enrichment.csv"
        if enrichment.exists():
            table = pd.read_csv(enrichment)
            fig, ax = plt.subplots(figsize=(7.2, 9.4))
            valid = table[(table["overlap"] > 0) & table.group.isin(["KO", "OE 2"])].sort_values("p_adj")
            valid = valid.groupby(["group", "library"], sort=False).head(2).iloc[::-1]
            if len(valid):
                # Wrap full term names; omit only terminal Reactome R-HSA identifiers for display.
                short_libraries = {"GO_Biological_Process_2025": "GO BP", "KEGG_2021_Human": "KEGG", "Reactome_2022": "Reactome"}
                labels = [f"{r.group} | {short_libraries.get(r.library, r.library)}\n" + textwrap.fill(re.sub(r"\s+R-HSA-\d+$", "", r.term), width=43)
                          for r in valid.itertuples()]
                ax.barh(labels, -np.log10(valid.p_adj.clip(lower=np.finfo(float).tiny)), color="#347F70")
                ax.axvline(-np.log10(.05), color="#777777", ls="--", lw=.7)
                ax.set(xlabel="-log10 BH-adjusted ORA p", title="KO/OE ORA: top 2 per method and database")
                fig.subplots_adjust(left=.47)
            else:
                ax.axis("off")
                ax.text(.5, .5, "No pathway overlaps in the significant response sets\nNo enriched pathway claim", ha="center", va="center", fontsize=9)
            save(fig, "08", "Pathway_Enrichment", ["results/pathway_enrichment.csv"], "ORA背景为当前网络所有非靶基因；按集合与数据库内BH校正，虚线q=0.05；名称换行保留完整术语，仅省略末尾Reactome R-HSA编号，完整ID见结果表。")
    if requested("09"):
        topgenes = scores["OE_2"].distance.nlargest(20).index
        heat = pd.DataFrame({key: scores[key].distance for key in ["OE_1", "OE_1.5", "OE_2", "OE_3"]}).loc[topgenes]
        heat.to_csv(root / "results" / "strength_response_heatmap.csv")
        fig, ax = plt.subplots(figsize=(5.8, 4.6))
        im = ax.imshow(heat, cmap="viridis", aspect="auto")
        ax.set_xticks(range(4), ["1 (null)", "1.5", "2", "3"])
        ax.set_yticks(range(len(heat)), heat.index, fontsize=7)
        ax.set(xlabel="Positive target-edge multiplier", title="OE response across strengths (same WT)")
        fig.colorbar(im, ax=ax, label="Unsigned manifold distance", shrink=.8)
        save(fig, "09", "Response_Heatmap", ["results/strength_response_heatmap.csv"], "相同基因在强度1/1.5/2/3下的实际无方向距离；未逐行标准化。")
    if requested("10"):
        fig, axes = plt.subplots(1, 3, figsize=(7.2, 3.5))
        overview_genes = genes[:9]
        angles = np.linspace(0, 2 * np.pi, len(overview_genes) - 1, endpoint=False)
        overview_pos = {target: np.array([0., 0.])}
        overview_pos.update({g: np.array([np.cos(a), np.sin(a)]) for g, a in zip(overview_genes[1:], angles)})
        pd.DataFrame.from_dict(overview_pos, orient="index", columns=["x", "y"]).to_csv(root / "results/landscape_layout.csv")
        for ax, name, matrix in zip(axes, ["WT", "KO", "OE 2"], [wt, ko, oe]):
            sub = matrix.copy()
            sub.loc[sub.index != target] = 0
            draw_network(ax, sub, overview_genes, overview_pos, target, name, 20, edge_scale=oe.loc[target, genes].abs().max())
            ax.set_xlim(-1.7, 1.7)
        save(fig, "10", "Regulatory_Landscape", ["objects/WT.csv", "objects/KO.csv", "objects/OE_2.csv", "results/landscape_layout.csv"], "固定同一布局及边宽尺度的三状态SPI1最强8个出边邻居；展示网络模型，不是表达或细胞动力学景观。")
    index = root / "figures" / "figure_manifest.json"
    if args.figure != "all" and index.exists():
        existing = json.loads(index.read_text(encoding="utf-8"))
        manifest = [r for r in existing if r["number"] != args.figure] + manifest
    index.write_text(json.dumps(sorted(manifest, key=lambda x: x["number"]), indent=2, ensure_ascii=False), encoding="utf-8")


if __name__ == "__main__":
    main()
