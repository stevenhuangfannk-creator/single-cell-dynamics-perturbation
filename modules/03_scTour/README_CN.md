# scTour：表达矩阵中的细胞动力学 / Cellular dynamics

本模块使用 [LiQian-XC/sctour](https://github.com/LiQian-XC/sctour) 的真实 `sctour==1.0.0` API，复算 [Li 2023, Genome Biology](https://doi.org/10.1186/s13059-023-02988-9) Figure 2。原始 GEO 表直接读取，没有把不明来源的二次 h5ad 当成原始 counts。

运行状态与实际数值以 [复现报告](reports/REPRODUCTION_CN.md)、`results/*/run_summary.json` 和 [图件索引](../../gallery/FIGURE_INDEX_CN.md) 为准。`RUN_OK` 表示软件执行及对象验证通过；`PARTIAL` 表示真实计算已完成但严格原图逐点复现条件尚不足。Smoke Test 仅验 API，不支持发育结论。

本轮已完成11个正式400-epoch模型、完整八训练比例及Figure 2 a–i/组合图，另保存独立Smoke结果。全套10个主图通过PNG600、PDF≥5pt、SVG文字及实际布局/目视检查；A模型非单调结果、B1/B5终末阶段错序与未核实S2参数保留，模块仍为`PARTIAL`。

## 方法与证据边界 / Method

输入是细胞 × 基因的表达矩阵。NB/ZINB 模式要求 `.X` 为非负整数 counts；官方软件还提供 MSE 模式。本仓库的训练入口限定原始 counts 的 NB/ZINB，正式数据使用 NB，并保留 `layers['counts']`；归一化 MSE 流程没有在这里执行验证。scTour 不需要 spliced/unspliced；不要人为伪造这些矩阵。

编码器从表达学习时间 `t∈[0,1]` 和五维潜变量 `z`，神经 ODE 学习 `dz/dt=f(z,t)`。最终潜空间为 `alpha_z*z + alpha_predz*z_t`，`z_t` 来自 ODE 积分。时间是转录组排序估计，VF 是潜空间导数。静态单细胞表达能够支持动态假说，不能单独证明细胞发生真实谱系转化，也不等于实验 RNA velocity。

官方实现让编码器/解码器使用 CUDA，ODE 求解部分使用 CPU。GPU 可用时自动使用 GPU，否则 CPU；不把 CPU 中的 ODE 运算写成全 GPU。保存环境、实际训练比例、种子、loss、耗时、原时间和方向校正，使用官方 `save_model` 与 `predict.load_model` 检查时间推断往返误差。

`ptime_before_manual`已经包含官方基于检出基因数的automatic方向校正，仅表示没有做额外manual reverse；同理，抽样指标`rho_before_alignment`也不是未经校正encoder时间的相关。`ptime`与`rho_aligned`保留进一步生物学方向审查后的结果；阶段并非严格反向时不会为了相关或外观强制翻转。

## 原始数据与论文配置 / Source data

| 数据 | 原始细胞 | 论文谱系细胞 | 基因筛选 | HVG | 原始空间 | 潜空间混合 |
|---|---:|---:|---|---:|---|---|
| GSE95315 | 5,454 | 4,007 | 检出 ≥20 个谱系细胞 | 500 | PCA40 UMAP | 0.2 z + 0.8 z_t |
| GSE104323 | 24,185 | 15,174 | 检出 ≥20 个谱系细胞 | 2,000 | PCA30 UMAP | 0.4 z + 0.6 z_t |

A 内嵌 age/cluster 两行元数据；B 用独立 metadata 文件按 cell ID 对齐。B 原始 `nIPC-perin` 与 `nIPC` 合并为原图的 `nIPC`，原标签完整保留。A 神经母细胞的两个类别在 B 原注释中只有 `Neuroblast`，不能照搬七色模板拆出未经证实的类别。

原始文件 URL、压缩大小、SHA256、gzip CRC、完整筛选数见 [A 清单](../../data_registry/sctour_A_preprocessing.json)、[B 清单](../../data_registry/sctour_B_preprocessing.json)。正式输入为基因筛选后的 counts，用 `seurat_v3` 在原始 counts 上选择 HVGs。每个文件都读取到 gzip EOF 校验 CRC，再生成 QC h5ad 与 HVG 表。

`config/paper_A.json` 使用 Methods 的 500 HVGs；官方后处理教程使用 2,000 HVGs，独立保存在 `config/tutorial_A.json`。教程配置不混入论文配置。论文 Results 使用20%训练表述，Methods 没有逐主面板给出全部 UMAP 参数、种子或固定比例。这里 seed 0、邻居 15、UMAP min_dist 0.1；A/B 主展示预登记按20%训练，B95%作为抽样稳定性参考，全部400 epochs。原先按软件默认的A90%两种子及B95%结果完整保留在诊断报告，不能因为其好坏选择性删去。主展示比例按论文依据事先登记，不能按视觉或发育趋势选择模型。

## 安装与运行 / Quick start

在仓库根目录运行以下 PowerShell 命令。需要 `uv`，可用 `python -m pip install uv` 安装；所有环境与数据都放在工作目录，避免影响其他项目。

```powershell
./environments/sctour-setup.ps1
$Py = './work/sctour/.venv/Scripts/python.exe'
$Mod = './modules/03_scTour'
& $Py "$Mod/scripts/prepare_geo.py" --dataset A --raw-dir work/sctour/raw --out-dir work/sctour/prepared --download
& $Py "$Mod/scripts/prepare_geo.py" --dataset B --raw-dir work/sctour/raw --out-dir work/sctour/prepared --download
& $Py "$Mod/scripts/run_sctour.py" --config "$Mod/config/paper_A.json" --input-path work/sctour/prepared/A_train_counts.h5ad --output work/sctour/smoke_A --smoke
& $Py "$Mod/scripts/run_sctour.py" --config "$Mod/config/paper_A.json" --input-path work/sctour/prepared/A_train_counts.h5ad --output work/sctour/full_A20
& $Py "$Mod/scripts/run_sctour.py" --config "$Mod/config/paper_B.json" --input-path work/sctour/prepared/B_train_counts.h5ad --output work/sctour/subsamples/B_95
& $Py "$Mod/scripts/run_subsampling.py" --config "$Mod/config/paper_B.json" --input-path work/sctour/prepared/B_train_counts.h5ad --output work/sctour/subsamples
& $Py "$Mod/scripts/plot_figure2.py" --a-path work/sctour/full_A20/inferred.h5ad --b-path work/sctour/subsamples/B_20/inferred.h5ad --subsample-dir work/sctour/subsamples --config "$Mod/config/plot.json" --output "$Mod/figures"
& $Py "$Mod/scripts/export_run_evidence.py" --work-root work/sctour --module $Mod
```

服务器可按 `environments/sctour-lock.txt` 创建等价 Python 3.11 环境。CPU 在安装时加 `-CPU`。注意现代 PyTorch 改变了 `torch.load` 默认行为；这里固定 2.5.1，与原软件检查点加载方式兼容。只加载本仓库可信来源的本地模型检查点。

正式运行输出 `inferred.h5ad`、`model.pth`、`loss.tsv`、`train_cells.tsv`、`validation_cells.tsv`、`inference.tsv.gz`、配置、环境与摘要。模型、QC/full对象及大结果留在 `work/`，GitHub 只保存登记、代码、图件及精简结果。`objects/README.md` 解释重建与保存路径。

`--load model.pth` 用于原输入的恢复：必须提供同序细胞、同序基因、相同 counts、相同 QC 和训练配置；入口在推断前核对它们。不同数据的预测应使用官方 `sct.predict` 新数据 API，不把恢复训练对象的时间直接赋给新输入。`run_subsampling.py` 的自定义 `--fractions` 必须包含95，或输出目录已有完成的95%参考。缓存会核对配置和输入；换数据、seed、HVG或超参数时使用新输出目录。

训练/验证划分沿用官方`split_index`：训练数量为`ceil(n_cells*percent)`，验证从其余细胞独立抽取，数量为`min(floor(n_train*0.1), n_cells-n_train)`。因此95%模型的验证比例受剩余细胞数限制，其他比例通常有相当于训练数量10%的独立验证细胞。实际ID及人数保存于`train_cells.tsv`、`validation_cells.tsv`和`results/run_comparison.tsv`，不把验证细胞计入训练百分比。

## 单图修改 / Redraw one panel

诊断运行也能从fresh clone重建：B95已在Quick start建立，A90需要以下两个独立运行。它们用于保留失败和随机种子证据，不替换A20主图。

```powershell
& $Py "$Mod/scripts/run_sctour.py" --config "$Mod/config/diagnostic_A90.json" --input-path work/sctour/prepared/A_train_counts.h5ad --output work/sctour/full_A
& $Py "$Mod/scripts/run_sctour.py" --config "$Mod/config/diagnostic_A90.json" --input-path work/sctour/prepared/A_train_counts.h5ad --output work/sctour/full_A_seed1 --seed 1
& $Py "$Mod/scripts/plot_figure2.py" --a-path work/sctour/full_A/inferred.h5ad --panel a --prefix A90_seed0 --config "$Mod/config/plot.json" --output "$Mod/figures/diagnostics"
& $Py "$Mod/scripts/plot_figure2.py" --b-path work/sctour/subsamples/B_95/inferred.h5ad --panel g --prefix B95_scTour --config "$Mod/config/plot.json" --output "$Mod/figures/diagnostics"
```

B95的b–f可把最后一条的 `--panel g` 分别改为b、c、d、e、f。`--prefix`精确产生现有诊断文件名，诊断重绘不修改主图索引；每张诊断图在manifest都有完整重绘命令。首次导出A90结果需再次运行 `export_run_evidence.py`。

例如只修改 `plot.json` 的 `stream_density`，重画 g，无需训练：

```powershell
& $Py "$Mod/scripts/plot_figure2.py" --b-path work/sctour/subsamples/B_20/inferred.h5ad --panel g --config "$Mod/config/plot.json" --output "$Mod/figures"
```

支持 `point_size/point_alpha/stream_density/stream_linewidth/stream_arrowsize/font_size/legend_size/output_dpi`。`colorbar_position` 保存default、i和composite的 `[left,bottom,width,height]` 画布相对坐标，可单独调色条位置。每个图导出 PNG600、PDF、SVG；组合图文件名为 `scTour_Figure2_Reproduction`。UMAP 纵横比例固定相等，禁止用拉伸或搬移细胞去拟合原图。官方高维邻域投影同时使用 `zs_key='X_TNODE'`、`vf_key='X_VF'`、`use_rep_neigh='X_TNODE'`、`E_key`、`t_key='ptime'`；`X_VF[:,:2]` 从未被当作 UMAP 位移。

在上述单图命令中改 `--panel`，并使用下表输入参数。推断对象旁必须有 `config_used.json` 和完成的400-epoch `run_summary.json`；Smoke对象不能用于Figure2。h需要全部八个比例的摘要及推断表，核对真实比例、全部同序细胞及一致标签；i必须是完成的20%正式模型。没有完成的面板会报错，不生成占位图或清除已有索引。

| 面板 | 输入参数 |
|---|---|
| a | `--a-path work/sctour/full_A20/inferred.h5ad` |
| b–g | `--b-path work/sctour/subsamples/B_20/inferred.h5ad` |
| h | `--subsample-dir work/sctour/subsamples` |
| i | `--subsample-dir work/sctour/subsamples`，读取其中 `B_20/inferred.h5ad` |
| composite | 同时提供a、b–g和h的输入参数；只重绘组合图，要求完整a–i结果 |

图 a 的时间和方向来自 A 模型，b/c/d 使用同一 B PCA UMAP，f/g 使用同一 B 潜空间 UMAP，i 使用独立20%模型的完整细胞推断。e/h 是真实细胞排序条带。inset 不重新计算坐标。

## 如何解读和替换自己的数据 / Interpretation

类别 UMAP 检查连续转录组结构；时间渐变检验已知发育阶段是否顺序合理；流线检验邻域间方向；排序条带显示批次是否挤占时间结构；抽样相关衡量不同模型一致性。三者都不能代替谱系追踪或真实干预实验。

对于 PI/Healthy 或其他自己的数据，先使用仓库共同输入适配器读取 h5ad、10x MEX、Seurat 转换出的 counts，验证基因/细胞 ID 唯一、整数、非负、有限。应先选同一可信连续谱系，保留样本批次和原注释；单纯 PI 与 Healthy 两组差异不构成时间顺序。本阶段没有运行任何用户 PI 数据结论。

共享入口产生 `validated_counts.h5ad` 后，以下准备命令可把常见的 `cell_type` 字段转为 runner 所需的分类字段 `celltype`，复制真实批次，计算全输入基因的 `n_genes_by_counts`，过滤低检出基因并选择HVG。`--celltypes` 必须填写经过生物学审查的实际类别；下面的 Fibroblast 只是替换示例，并不默认证明该类具有连续轨迹。

```powershell
# 先按项目首页建立 shared-env（包含PyYAML），共享检查使用它。
Copy-Item shared/pi_input.example.yaml work/my_input.yaml
# 先修改 work/my_input.yaml 的真实路径、counts_layer 和元数据字段，再运行：
$SharedPy = './work/shared-env/Scripts/python.exe'
& $SharedPy shared/prepare_input.py --config work/my_input.yaml --output work/my_input --assess-pi --export-mex
& $Py "$Mod/scripts/prepare_custom.py" --input-path work/my_input/mex --output work/my_sctour/prepared --celltype-key cell_type --batch-key batch --celltypes Fibroblast --hvg 2000
```

这里使用共享入口导出的原始MEX和元数据作为跨环境桥，避免新AnnData/Pandas的nullable-string编码无法被scTour旧版本读取。矩阵、基因、细胞及注释按字符串ID核对，保留`001`等前导零及字面基因名`NA`；整数counts仅在数值完全不变时转为模型需要的float32。公开A数据已实际通过这一跨环境桥。已有兼容的旧格式counts h5ad也可直接传给 `prepare_custom.py`。

复制论文配置至自己的配置文件，令 `hvg` 与上一步一致，并检查 `pca_components` 小于输入维度。`stage_order` 只填写有独立生物学证据的阶段；也可设 `[]`，此时阶段方向输出 `NOT_ASSESSED`、相关和单调性为null，不做手动反转。随后先用 `--smoke` 检查接口，再决定是否开展正式分析。Figure2绘图脚本的A/B配色对应论文类别；新数据需据实际类别修改 `PALETTES` 和排序条带，不把 PI 类别套进神经发生标签。未覆盖的类别会显式报错。

如果没有可靠阶段注释，不应自动用人为顺序校正时间；先审查 marker、样本年龄、独立证据再决定方向。仅终末细胞、强周期过程、多分支根状态含混或离散炎症亚群时要谨慎使用：论文承认循环过程、不同终末状态之间的伪动力学及复杂分支串接风险。与 RNA velocity 相比，本方法不需转录剪接信息，但独立生物学证据仍不可省略。

## 验证 / Validation

```powershell
$env:SCTOUR_RESULT_DIR = 'work/sctour/full_A20'
$env:SCTOUR_SMOKE_DIR = 'work/sctour/smoke_A'
& $Py -m unittest discover -s tests -p 'test_sctour*.py' -v
```

检查整数/非有限/负值/重复 ID，实际 h5ad、时间范围、五维潜变量/VF及模型往返。绘图源检查、PDF 字体检查和逐面板目视审查记录在报告。只有运行成功不代表严格重现，也不代表因果正确。
