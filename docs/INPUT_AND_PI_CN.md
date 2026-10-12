# 数据接入与 PI/Healthy 审查 / Input and future adaptation

本阶段只运行公开数据。以下入口用于未来数据格式审查，不输出 PI/Healthy 的正式生物学结论。

## 原始计数要求 / Raw counts

Python/AnnData 按细胞×基因管理；R/scTenifold 输入按基因×细胞管理。必须追溯到原始非负整数 counts；数值为整数是必要检查，不能单凭整数性证明数据确为原始 counts。归一化、log、scale 数据不能交给 NB 模型冒充原始 counts。

基因和细胞 ID 必须唯一，metadata 必须按 ID 对齐。工具拒绝重复名称，避免静默修改靶基因。10x 的重复 symbol 可显式选择 `gene_id: id` 使用唯一 feature ID，但靶点也必须改成同一 ID 系统。多组学 MEX 只保留 `Gene Expression`；此筛选须结合原始 feature 表审查。

## 安装与运行 / Run

在仓库根目录执行。以下通用 Python 环境只用于输入和表格审查；三个模型使用各自模块环境。

```powershell
python -m venv .venv
.venv\Scripts\python -m pip install -r environments/shared-requirements.txt
.venv\Scripts\python shared/prepare_input.py --config shared/pi_input.example.yaml --output data/pi_checked --assess-pi
```

Linux/macOS 将 `.venv\Scripts\python` 替换为 `.venv/bin/python`。先复制并修改配置中的真实文件与 metadata 字段。h5ad 默认读取 `layers['counts']`；只有检查原对象后才可显式设置 `counts_layer: X`。原始对象保留不变，导出 `validated_counts.h5ad` 和 `input_audit.json`。

10x MEX 的 `input` 指向含 matrix.mtx、features.tsv/genes.tsv、barcodes.tsv 的目录，均可 gzip。metadata 为首列细胞 ID 的 TSV，必须覆盖全部输入细胞；已有注释列遇到重名会报错，以避免覆写。

Seurat 可先转换为 MEX，并保持 metadata 行顺序：

```powershell
Rscript shared/export_seurat.R data/input.rds data/my_counts RNA
```

导出依赖 SeuratObject 和 Matrix。Seurat v5 有多个 counts layer 时须先审核样本并明确 JoinLayers；脚本拒绝擅自合并。这一 Seurat 导出入口的版本兼容性须在实际用户对象上再检查，不能视为已分析其私有数据。

## 细胞类型与设计 / Study design

重点评估 Fibroblast、Macrophage、Monocyte、Neutrophil、Osteogenic lineage 和 Endothelial cells。每种类型分别检查细胞数、供体数、批次、目标检出率与基因信息量；不要把不同类型合成网络后解释为类型特异调控。候选基因须在分析前依据生物学和检出率选择，不能依据显著结果反选。

`--assess-pi` 要求 cell_type、condition、batch、donor 等字段，报告条件内的供体数、细胞数和批次交叉表。条件与批次/供体完全混杂时，应先修正设计或限定结论。正式网络通常需要足够同类型细胞；具体最小阈值与抽样规模以模块配置及稳定性结果为准。

scTour 的使用还需外部发育阶段、连续状态、标记基因与合理组织结构等证据。PI 与 Healthy 是条件标签，不能自动设为起点与终点。必须把计算拟时序和真实发育时间区分开；现有审查脚本不自动宣称轨迹可信。
