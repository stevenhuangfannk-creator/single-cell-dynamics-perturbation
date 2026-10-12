# scTenifoldOE 虚拟过表达 / Experimental virtual overexpression

本模块完成了真实 PBMC CD14+ Monocytes 的 SPI1 网络过表达模拟、1.5/2/3 强度分析以及同一 WT 上的 KO/OE 比较。计算状态为 **RUN_OK**；方法验证状态为 **Experimental / PARTIAL**，不能称为正式官方 scTenifoldOE，也没有真实过表达实验验证。

实际 OE 算法来自第三方 [`virt_overexpression 0.1.0`](https://github.com/virchow-zhang/virt_overexpression/tree/860d0ea59dc0aeb564a286efae00f562ce84e138)。调用源码增强靶点行的正出边，再调用官方 R manifold alignment 与 differential regulation。详细许可、维护、方向矛盾、API 和兼容性审查见 [SOURCE_AUDIT_CN.md](reports/SOURCE_AUDIT_CN.md)。所有图页脚注明算法来源。

## 方法及输入 / Method and inputs

它回答“在由这一类细胞拟合的网络中，增强候选基因的正调控关联边，哪些基因的网络表示变化较大”。这是计算模拟产生的功能假设。网络差异不能直接变成真实表达矩阵、表达 log2FC 或实验因果结论。

必须先在模块 01 完成 counts 检查、细胞类型选择、候选靶点表达审查与官方 WT 构建。这里输入 **相同基因顺序的 WT.csv 与 KO.csv**，WT 行为 source、列为 target；不因过表达重新抽样或重新拟合基础网络。要求 ID 唯一、有限值、对角线为零、靶点存在且目标行含正边。若无正出边，默认停止并报告信息量问题，不自动复制入边或填造相关边。

真实验证数据是 10x PBMC3k 原始 counts 与 Scanpy 公开处理对象的细胞注释匹配得到的 480 个 CD14+ Monocytes，网络含 600 个基因；SPI1 检出率为 63.33%。基因面板与 WT 来源见 [模块 01](../01_scTenifoldKnk/README_CN.md)。这是一个小规模真实数据验证，尚不是全转录组预测或供体重复的生物学验证。

OE 操作是 `W[g,j] *= factor`，只对 `W[g,j]>0` 执行。负出边不变，因此它与将全部目标出边归零的 KO 在符号处理上不同。边强度 2 不是 mRNA 表达 2 倍。随后将 WT 与扰动网络进行同构 manifold alignment，按基因比较表示距离，使用官方 Knk 的卡方和 BH 规则。主分析不使用方向预测。

## 安装 / Installation

建议 Python 3.11 专属环境。pandas 3 的同质 DataFrame 会触发上游只读 `.values` 错误；本模块固定 pandas 2.2.3。以下 PowerShell 命令从仓库根执行，所有分析始终使用 `$Py` 指定的同一环境。若没有 `py` launcher，可用已安装的 Python 3.11 可执行文件代替 `py -3.11`。

```powershell
py -3.11 -m venv work/oe-env
$Py = (Resolve-Path work/oe-env/Scripts/python.exe).Path
& $Py -m pip install -r environments/oe-requirements.txt
& $Py modules/02_scTenifoldOE/scripts/install_upstream.py --cache work/oe-upstream
```

安装脚本锁定两个 GitHub revision，记录下载 ZIP 的 SHA256，以 `--no-deps` 安装，避免上游最低版本声明自动升级已经验证的依赖。torch CPU 版本足够用于本模块。R 依赖为 `scTenifoldNet 1.4.3`、`scTenifoldKnk 2.0.0`，安装步骤沿用模块 01。上游 OE 许可声明材料不完整，源码不复制进本仓库。

## 主比较运行 / Shared-WT analysis

用实际 R 安装路径替换示例中的路径，`--r-library` 指向安装模块 01 R 包的 library；下面命令不含用户私有数据。

```powershell
& $Py modules/02_scTenifoldOE/scripts/run_oe.py --wt modules/01_scTenifoldKnk/objects/WT.csv --ko modules/01_scTenifoldKnk/objects/KO.csv --target SPI1 --folds 1 1.5 2 3 --ma-dim 2 --seed 1 --rscript "C:/Program Files/R/R-4.3.1/bin/Rscript.exe" --r-library work/knk/Rlib
& $Py modules/02_scTenifoldOE/scripts/download_annotations.py --output work/oe-gmt
& $Py modules/02_scTenifoldOE/scripts/enrich_compare.py --gmt work/oe-gmt
& $Py modules/02_scTenifoldOE/scripts/plot_oe.py --dpi 600
& $Py -m unittest discover -s tests -p "test_oe*.py" -v
```

`--folds 1` 是无扰动负对照，1.5/2/3 是网络剂量敏感性。MA 固定 `d=2,nCores=1`、seed 1，统计为 `empiricalNull=FALSE,direction=NULL`。`--output` 可指定另一个完整模块目录，避免覆盖既有结果；`plot_oe.py --module` 应读取同一路径。`config/pbmc_cd14.json` 保存本次参数记录；脚本运行参数为真实入口，配置文件不伪称被自动读取。

完整10图与富集套件固定读取1/1.5/2/3，运行整套时保留这四个强度。`run_oe.py`可单独计算其他强度，但缺少这些默认表时不能直接调用现有完整绘图/富集套件。自定义输出后，富集也要使用`enrich_compare.py --module 同一路径`。图件manifest和部分中文图注为本次PBMC/SPI1示例；更换数据或靶点时须同步修订来源和注释，不能沿用PBMC登记。

运行导出 `objects/{WT,KO,OE_1,OE_1.5,OE_2,OE_3}.csv`、每个条件的 `_alignment.csv` 与 `_dr.csv`、RDS 对象、运行版本/哈希/耗时、强度汇总、KO/OE 共同与特异基因以及 GO/KEGG/Reactome ORA 全检验表。RDS 与原始矩阵按项目规则不上传，可由脚本重建；CSV 和图件保留供检查和重绘。

## 第三方公开 API 验证 / Published API execution

```powershell
& $Py modules/02_scTenifoldOE/scripts/verify_public_api.py --counts modules/01_scTenifoldKnk/objects/pbmc_cd14_counts.mtx --genes modules/01_scTenifoldKnk/objects/pbmc_cd14_genes.tsv --cells modules/01_scTenifoldKnk/objects/pbmc_cd14_cells.tsv --target SPI1 --output modules/02_scTenifoldOE/results/native_api
```

实际执行 README 的 `scTenifoldOECUDA` 高层 API，CPU、3 个 bootstrap 网络、每网 300 个细胞、seed 1、2 倍正边增强。输入由模块 01 先审核，演示关闭重复细胞/基因过滤。上游没有提供可下载示例数据或预期结果，所以状态为 API **RUN_OK**，不能称为原论文 OE 数值复现。其 WT 使用 bootstrap mean 近似，并与主官方 CP-WT 比较完全分开保存。

## 图件及单图修改 / Figures and redraw

10 张真实图均输出白底 600 DPI PNG、PDF 与可编辑文本 SVG；[figure_manifest.json](figures/figure_manifest.json) 记录输入表、脚本、数据、科学边界与文件路径。

| 编号 | 解读 |
|---|---|
| 01 | SPI1 邻域 OE 网络：目标边优先及最强边，示图布局无生物学坐标含义 |
| 02 | OE−WT 邻接差值：检查只有目标行正边发生预定变化 |
| 03 | 排除 SPI1 的无方向距离排序：绿色 FDR<0.05，灰色非显著 |
| 04 | SPI1 中心出边：绿色正关联、紫色负关联，非验证的激活/抑制 |
| 05 | 1/1.5/2/3 网络差异与距离：单一 WT 的剂量敏感性，无重复误差条 |
| 06 | 统一 WT/MA/DR 的 KO 和 OE 无方向分数对比 |
| 07 | FDR<0.05 的共享、KO 特异、OE 特异集合，严格排除靶点 |
| 08 | GO/KEGG/Reactome、每方法每库最前 2 项 ORA；全通路表保留零命中 |
| 09 | 相同基因在不同网络强度下的实际距离，未逐行标准化 |
| 10 | 同一布局/绝对边宽尺度的 WT、KO、OE 靶点网络，不是细胞动力学轨迹 |

只改图 05 或图 09 时，例如：

```powershell
& $Py modules/02_scTenifoldOE/scripts/plot_oe.py --figure 05 --dpi 600
& $Py modules/02_scTenifoldOE/scripts/plot_oe.py --figure 09 --dpi 600
```

修改 `plot_oe.py` 的 figsize、点大小、字体或色图即可，无需重新构建 WT 和重新比较网络。网络固定布局写在 `results/network_layout.csv`；所有网络图按相同确定性规则生成，不能把图的 node 位置解释为单细胞表达空间。脚本缺少富集结果时不会生成图 08 来冒充已完成分析。

## 本次结果与适用边界 / Results and limits

SPI1 共有 567 条非零出边，其中 284 条正边被 OE 放大。2 倍 OE 有 9 个非靶显著响应，KO 有 14 个；7 个共同、7 个 KO 特异、2 个 OE 特异。1.5/2/3 的前 20 排序 Jaccard=1，距离 Spearman≈0.999999，表示这一 WT 下剂量排行稳定。未证明不同 WT、不同供体、细胞数或算法下 OE 预测稳定。

共同响应中的多个 HLA-II 基因支持抗原呈递相关假设。KEGG 会重复命中多种带 HLA 模块的疾病通路，不能把这些条目直接解释成 PBMC 供体患病或 SPI1 对具体疾病有因果作用。KO/OE 特异集合在三数据库中未得到 q<0.05 的通路，保留这个阴性结果。

不建议在靶点低检出/无正出边、混合不同细胞类型、极小细胞数、批次驱动网络或未经 QC 的归一化值上直接运行。未来 PI / Healthy 可先经 shared 接口读取 Seurat 导出 counts、AnnData counts layer、10x MEX，按 Fibroblast/Macrophage/Monocyte 等细胞类型构建独立 WT，预先登记候选基因与检出阈值。不同疾病组是比较条件，不是时间顺序。当前没有运行用户 PI 数据。

仍未完成：独立实验 benchmark、跨供体 OE 生物学验证、跨随机 WT 的 OE 稳定性，以及上游缺失 LICENSE 的澄清。本次运行成功不解除这些限制。
