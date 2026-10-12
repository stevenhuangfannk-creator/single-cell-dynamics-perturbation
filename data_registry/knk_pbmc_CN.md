# PBMC 数据登记 / Data registry

本模块下载日期：2026-10-12。物种 Homo sapiens；10x PBMC3k 为健康供者外周血单核细胞。此数据是官方公开技术示例，未提供 GEO accession，不能杜撰 GEO 编号。

| 资源 | 下载及用途 | 文件大小 | SHA-256 |
|---|---|---:|---|
| 原始 UMI counts | [10x MEX](https://cf.10xgenomics.com/samples/cell-exp/1.1.0/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz)；所有网络计算的唯一表达来源 | 7,621,991 bytes | `847d6ebd9a1ec9a768f2be7e40ca42cbfe75ebeb6d76a4c24167041699dc28b5` |
| 处理对象 | [cellxgene h5ad](https://raw.githubusercontent.com/chanzuckerberg/cellxgene/main/example-dataset/pbmc3k.h5ad)；仅转移细胞 barcode 对应的已发表 `louvain` 细胞类型标签 | 24,653,425 bytes | `0db367b991dd95809732b218539ede489bea99113807f62ebd7ccc970025fe38` |

两个校验值均与 [Scanpy 官方当前 registry](https://github.com/scverse/scanpy/blob/main/src/scanpy/datasets/registry.yaml) 一致。标签处理过程见 [Scanpy 官方 PBMC 教程](https://scanpy.readthedocs.io/en/latest/tutorials/basics/clustering-2017.html)。**处理 h5ad 的 scaled / regressed `.X` 未用于网络输入**。

原始 MEX：32,738 features × 2,700 cells。95 个重复 symbol 通过 counts 求和合并，得到 32,643 唯一 symbol。注释对象 2,638 cells，其中 480 个 CD14+ Monocytes；barcode 全部在原始 counts 中匹配。原始 counts 复查 library ≥500、mitochondrial fraction <0.1 后仍有 480 cells。删除核糖体/线粒体基因并要求检出率 ≥5%，3,363 个基因合格；按 log1p(CPM) 方差选取 600 genes，强制保留符合相同表达质检的预先指定调控候选。完整排除与选择记录在 [gene_selection_audit.csv](../modules/01_scTenifoldKnk/results/gene_selection_audit.csv)。

SPI1 检出率 63.33%，均值 1.3604 UMI/cell。NFKB1 检出率 4.79%，按预登记规则排除；不以显著性决定候选。正式 screen 是 9 个候选分别 KO，不代表完成真实全转录组筛查。1200-gene sensitivity 按同一规则扩大 panel，保留所有 seed 结果。

下载与重建：`python modules/01_scTenifoldKnk/scripts/prepare_pbmc.py --data work/knk/data`。原始下载、counts、RDS、h5ad 均在 Git ignore 范围，不能仅从 Git clone 即读取大对象；先执行此脚本和 R 分析入口。正式结果 CSV 与重绘图件可直接浏览。

富集注释从 [Enrichr 官方 gene-set libraries](https://maayanlab.cloud/Enrichr/#libraries) 下载 `GO_Biological_Process_2025`、`KEGG_2021_Human`、`Reactome_2022`。下载 URL、大小、SHA-256、term 数量存于 [enrichment_summary.json](../modules/01_scTenifoldKnk/results/enrichment_summary.json)。本项目在本地计算 ORA，背景是已分析的 600 genes；并未向网络服务上传用户数据。

软件来源与 commit SHA 见 [knk_software_provenance.json](knk_software_provenance.json)，R 依赖见 [knk_R_versions.csv](../environments/knk_R_versions.csv)。
