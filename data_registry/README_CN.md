# 公开数据索引 / Public data registry

本次只分析公开 counts。文件 SHA-256、大小及处理记录均来自实际下载；模型和原始文件不直接上传 GitHub，可用下列脚本重建。不同模块对同一输入的用途应分别核查。

| 数据 / Dataset | 物种与细胞 | 用途 | 原始与分析规模 | 下载、检查与记录 |
|---|---|---|---|---|
| 10x PBMC3k | 人，健康供者 PBMC；分析 CD14+ Monocytes | KO 与 OE 共用同一官方 WT | 原始2,700细胞、32,738 features；480单核细胞、600基因；1200基因敏感性另存 | [counts与标签来源、SHA](knk_pbmc_CN.md)；[下载及预处理](../modules/01_scTenifoldKnk/scripts/prepare_pbmc.py) |
| GSE95315 | 小鼠，齿状回颗粒细胞谱系 | scTour Figure 2a | 原始5,454细胞；选4,007细胞；检出≥20细胞后11,819基因，500 HVGs | [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE95315)；[实际下载清单](sctour_A_preprocessing.json)；[下载及预处理](../modules/03_scTour/scripts/prepare_geo.py) |
| GSE104323 | 小鼠，颗粒细胞与未成熟锥体神经元谱系 | scTour Figure 2b–i | 原始24,185细胞；选15,174细胞；过滤后15,015基因，2,000 HVGs | [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE104323)；[实际下载清单](sctour_B_preprocessing.json)；[下载及预处理](../modules/03_scTour/scripts/prepare_geo.py) |

PBMC 的处理 h5ad 仅提供 barcode 对应注释，scaled `.X` 未输入模型。GSE95315 的标签内嵌在原始表中；GSE104323 的 metadata 按 cell ID 对齐。完整筛选与阶段数见模块复现报告，而不是由 UMAP 外观反推。

在仓库根目录、已安装对应环境后运行：

```powershell
work/shared-env/Scripts/python.exe modules/01_scTenifoldKnk/scripts/prepare_pbmc.py --data work/knk/data --n-genes 600
work/sctour/.venv/Scripts/python.exe modules/03_scTour/scripts/prepare_geo.py --dataset A --raw-dir work/sctour/raw --out-dir work/sctour/prepared --download
work/sctour/.venv/Scripts/python.exe modules/03_scTour/scripts/prepare_geo.py --dataset B --raw-dir work/sctour/raw --out-dir work/sctour/prepared --download
```

GEO脚本对已下载的固定文件核对大小、SHA及gzip CRC；网络不可用时可以保留已校验的本地文件运行。PBMC重建方法及可能随上游变化的下载地址见其专门登记。

GO/KEGG/Reactome 的公开 GMT 注释来自 Enrichr，库版本、完整 URL 和 SHA 保存在 [KO富集摘要](../modules/01_scTenifoldKnk/results/enrichment_summary.json) 与 [OE注释来源](../modules/02_scTenifoldOE/reports/pathway_provenance.json)。背景限定为实际网络基因；疾病名称通路不是对供者的诊断。

软件核查另见 [官方KO来源](knk_software_provenance.json)、[第三方OE来源](oe_sources.json) 和 [参考文献](../docs/REFERENCES.md)。未来私有 PI/Healthy 数据要求见 [输入与设计审查](../docs/INPUT_AND_PI_CN.md)；本次未运行私有数据正式结论。
