# scTenifoldKnk：单细胞网络虚拟敲除 / Virtual knockout

本模块已真实执行官方 R `scTenifoldKnk 2.0.0` 与 `scTenifoldNet 1.4.3`。官方模拟示例、PBMC 单基因 KO、联合 KO、9 候选筛查、方向敏感性、GO/KEGG/Reactome 富集、负对照以及 9 组种子/参数运行均有结果。**软件执行状态 `RUN_OK`；科学复现状态 `PARTIAL`：PBMC 排名稳定性差，结果只能用于探索，不能声称可靠机制或真实 CRISPR 效应。**

| 入口 | 内容 |
|---|---|
| [复现报告](reports/REPRODUCTION_CN.md) | 参数、实际数值、不稳定结果与未完成事项 |
| [运行配置](config/pbmc.json) | 已执行 600 genes / 480 cells 的参数 |
| [图件索引](figures/figure_manifest.json) | 11 图和重绘脚本一一对应 |
| [数据登记](../../data_registry/knk_pbmc_CN.md) | 原始 counts、标签、来源、校验与排除规则 |
| [版本表](../../environments/knk_R_versions.csv) | 实际 R 软件版本 |
| [新数据配置蓝图](config/own_data.template.json) | PI / Healthy 未来接入要求，尚未运行用户数据 |

## 方法与解释 / Method

问题是“在一个特定细胞类型的 WT 调控网络中，去除某个候选基因的网络作用后，哪些基因的网络表征发生较大变化”。输入必须是 genes × cells 的非负原始 counts，基因和细胞 ID 唯一，候选存在且具有足够表达/网络信息。

流程：原始 counts → 细胞/基因 QC → CPM normalization → 对多个细胞子样本进行 principal component regression → CP tensor decomposition 降噪 → WT GRN → 将 KO 基因对应行出边归零 → WT/KO manifold alignment → 基因扰动 distance → χ²(df=1) 与 BH → 功能富集。

官方最终 `tensorNetworks$WT[i,j]` 表示 regulator i → target j；源码在 PCR / tensor 后执行了转置。本仓库直接保存官方最终矩阵，未再次转置。KO 是对网络的操作，**不产生 KO 后真实单细胞表达矩阵**。网络权重符号是 PCR 关联，不足以证明激活/抑制或因果边。

`distance` 是 WT/KO 低维表征的欧氏距离；`FC` 是 distance² 相对于非 target distance² 均值的比值，**不是表达 fold change，也不应标为传统差异表达 log2FC**。2.0 的 `direction` / `directionScore` 是从 WT log1p(CPM) 相关 heat kernel 单独推断的方向，与 MA distance 不同。本分析主方向使用 `dr_directionRegressLibSize=TRUE`，同时保存未回归 library size 的方向以供敏感性检查；方向没有独立实验验证。

2.0 新增 `heatKernel`、`hkManifoldAlignment`、`knockoutDirection`、`perturbationMap`。Single / joint KO 默认 `ma_method="manifold"`；`transcriptomeWide=TRUE` 默认 heat。后者根据 WT heat kernel 与平均 logCPM 扩散计算 distance，不能把 heat 距离数值与 MA 距离直接混合比较。这里真实数据只筛查 9 个预先指定候选；全 100 基因 heat screen 仅在官方模拟输入上执行。

## 已验证运行 / Quick start

以下命令从**仓库根目录**运行，`python` 应指向 [shared-requirements.txt](../../environments/shared-requirements.txt) 安装的共享环境；本次完整版本记录见 [shared.lock.txt](../../environments/shared.lock.txt)（Python3.13）。Windows 本次使用 R 4.3.1；分析全程 `nCores=1`，避免影响 scTour GPU 训练。

```powershell
# 独立 R library；不覆盖其他科研项目的环境。
$env:KNK_R_LIB = "work/knk/Rlib"
Rscript environments/knk_install.R

# A：按官方 README 的负二项示例生成 counts，保留 QC / MA / DR / heat screen。
Rscript modules/01_scTenifoldKnk/scripts/run_official_example.R

# B: Download raw 10x counts, transfer published labels by barcode, audit targets.
python modules/01_scTenifoldKnk/scripts/prepare_pbmc.py --data work/knk/data --n-genes 600
Rscript modules/01_scTenifoldKnk/scripts/run_pbmc.R

# 富集在本地计算，自动下载公开 GMT；背景限定600个分析基因。
python modules/01_scTenifoldKnk/scripts/enrich.py --gmt work/knk/data/gmt

# 保存对象后可独立绘图，无需再次网络建模。
python modules/01_scTenifoldKnk/scripts/plot_figures.py --figure all --dpi 600
python modules/01_scTenifoldKnk/scripts/plot_figures.py --figure 10 --dpi 600

# 新环境中验证关键科学不变量。
python -m unittest discover -s tests -p "test_knk*.py"
```

`Figure10` 含额外1200-gene敏感性结果，首次从零重建前应执行下面一段，再执行绘图；其他10张可单独重绘，不依赖1200结果。

```powershell
python modules/01_scTenifoldKnk/scripts/prepare_pbmc.py --module modules/01_scTenifoldKnk/objects/sensitivity1200 --data work/knk/data --n-genes 1200
Rscript modules/01_scTenifoldKnk/scripts/run_sensitivity.R
```

本次安装过程中 `rjson` 源码缺少 Windows compiler、`reshape2` 缺少 `plyr` 引发失败，已通过对应 R 4.3 binary 安装补齐，随后两个官方包从 pinned CRAN source 成功安装。日志保留。`environments/knk_install.R` 对官方包固定版本及源 tarball MD5，其他依赖版本以实际版本表为准；这是已验证环境记录，尚非容器级完整依赖锁。CRAN未来可能撤下旧R二进制；安装器未实现所有依赖的archive-binary回退。若旧R依赖下载受阻，需升级R或安装匹配Rtools后按版本表源码重建，不能宣称当前安装器已在完全空白R环境通过。

## 参数与输出 / Parameters and outputs

`config/pbmc.json` 的数值模型/稳定性参数被 `run_pbmc.R` 实际读取。数据选择规则由 `prepare_pbmc.py` 定义；`run_pbmc.R` 和 `plot_figures.py` 是固定PBMC演示，目标SPI1、联合SPI1+CEBPB及对应文件名保持固定，JSON的目标说明字段不能自动切换这套示例。9个候选均先通过≥5%检出门槛，NFKB1因4.79%检出按规则排除。新数据请使用下面的generic入口，不能仅修改JSON目标字段便宣称已替换分析。

正式参数：`nc_nNet=8`、`nc_nCells=400`、`nc_nComp=3`、`nc_q=0.9`、`nc_lambda=0`、`td_K=3`、`td_maxIter=1000`、`td_maxError=1e-5`、`ma_nDim=2`、`seed=1`。原始 counts 已在 Python 按完整转录组检查 library/mito，再选择600基因，所以 R 设 `qc=FALSE`，避免对子集错误重复应用全转录组 library 阈值。

| 路径 | 真实内容与解释 |
|---|---|
| `objects/pbmc_cd14_counts.mtx` + genes/cells TSV | 原始 UMI 子集，600 × 480；不是 processed h5ad `.X` |
| `objects/WT.csv` / `KO.csv` / `multi_KO.csv` | 完整网络；RDS counterpart 可快速读取 |
| `results/WT_edges.csv` / `KO_edges.csv` | regulator、target、weight 的边表 |
| `results/single_SPI1.csv` / `multi_SPI1_CEBPB.csv` | 完整非筛选 DR 表，含 target；下游解释排除 target |
| `results/single_SPI1_MA.csv` | 1200 × 2 的 WT/KO 表征（X_ / Y_ 前缀） |
| `results/candidate_screen_distance.csv` / `_direction.csv` | 9 × 600 heat-screen，非完整转录组 |
| `results/candidate_perturbation_map.csv` | 官方 `perturbationMap(plot=FALSE)` 输出，参考为同模型 SPI1 heat profile |
| `results/enrichment_all_terms.csv` | 所有符合背景大小门槛的 ORA term，含不显著及零 overlap |
| `results/stability_summary.csv` / `sensitivity1200/` | 全部 seed / sample / panel 的差异，不能选择性删除 |
| `results/negative_control_WT_vs_WT.csv` | 无变化网络负控，全部 p=1 |
| `figures/Knk_Figure01`–`11` | PNG 600dpi、PDF、SVG，白底与可编辑字体 |

对象中的 `rds`、原始下载、counts / h5ad 不上传 GitHub，由脚本重建；网络 CSV 和全部结果表已保留，支持从 Git clone 直接进行大多数重绘。图04/05距离和统计、图06局部网络、图07富集、图08heat结果、图09single/joint、图10稳定性、图11候选map表达不同问题；不应据网络图美观程度评判方法有效性。

## 新数据适配与限制 / Reuse and limitations

未来 Seurat 应导出原始 counts，AnnData 应明确 `layers['counts']` 而不是 scaled `.X`，10x MEX/h5ad可通过 [统一输入脚本](../../shared/prepare_input.py) 检查轴、整数性、唯一 ID、目标、细胞数；[统一 Seurat 导出](../../shared/export_seurat.R) 和 [PI 数据配置](../../shared/pi_input.example.yaml) 提供入口。`own_data.template.json` 是要求蓝图，尚未实现从该 JSON 自动运行分析；不得把蓝图状态写成已运行。

通过审计、限定单细胞类型及合理gene panel后，通用 [run_custom.R](scripts/run_custom.R) 可直接读取genes × cells的MEX与对应唯一ID文件，单目标或逗号分隔联合目标均可。此入口要求新/空output目录，不覆盖历史；默认最多3000genes、20networks、4cores，先验证输入与目标≥5%检出。输入QC需此前在全转录组完成，默认不再重复QC。

```powershell
Rscript modules/01_scTenifoldKnk/scripts/run_custom.R audited_counts.mtx genes.tsv cells.tsv work/my_KO SPI1,CEBPB
# 可选第6个参数只覆盖共有模型参数，例如已有pbmc.json。
Rscript modules/01_scTenifoldKnk/scripts/run_custom.R audited_counts.mtx genes.tsv cells.tsv work/my_KO_single SPI1 modules/01_scTenifoldKnk/config/pbmc.json
```

输出完整`output.rds`、WT/KO、DR、MA、参数/版本/检测率/边数；状态`UNVALIDATED`提醒新样本还需seed及供者稳定性。当前11图脚本专用于PBMC固定文件/目标，新数据可读取generic输出CSV自行改绘图目标/label；本阶段未实现一条命令生成任意输入的完整11图，也未在PI数据运行。

分别对 Fibroblast、Macrophage、Monocyte、Neutrophil、Osteogenic lineage、Endothelial cells 形成有足够细胞和供者支持的网络。先审查 donor、batch、ambient RNA、目标检出、不同供者重复与seed稳健性；混合细胞类型网络容易将细胞身份差异当作调控。官方建议约500 cells或更多，本示例480是小规模验证；新研究应扩大细胞/基因范围，调整网络数并重新验证，不直接照搬不稳定结果。

仅静态共表达及网络操作不足以支持因果推断，必须以真实 perturb-seq/CRISPR、独立样本或相关实验证据校验。对极低表达目标、极少细胞、无出边目标、高批次混杂、跨细胞类型混合、没有独立重复的数据，不建议直接给出机制结论。PI vs Healthy 是疾病状态比较，不能直接作为发育时间。

与 scTenifoldNet 相比，Knk只需WT网络，通过人工网络KO比较；Net原本比较两个观测条件网络。与OE模块相比，KO是目标行归零，OE使用第三方实际边增强算法且被标记 Experimental；二者均非真实表达矩阵生成器。本报告保留不稳定结果和统计假设局限，未对用户PI/Healthy数据开展正式分析。

## 参考 / References

- Osorio D et al. (2022). *scTenifoldKnk: An efficient virtual knockout tool for gene function predictions via single-cell gene regulatory network perturbation*. Patterns 3, 100434. [DOI](https://doi.org/10.1016/j.patter.2022.100434) · [开放全文](https://pmc.ncbi.nlm.nih.gov/articles/PMC9058914/)
- [scTenifoldKnk 官方源码/README/NEWS](https://github.com/cailab-tamu/scTenifoldKnk) · [CRAN](https://cran.r-project.org/package=scTenifoldKnk)
- [scTenifoldNet 官方源码](https://github.com/cailab-tamu/scTenifoldNet)，原方法 DOI [10.1016/j.patter.2020.100139](https://doi.org/10.1016/j.patter.2020.100139)
- Enrichr annotation resource: Xie Z et al. (2021), [10.1002/cpz1.90](https://doi.org/10.1002/cpz1.90)。本实现使用独立本地超几何检验，非Enrichr web默认background。
