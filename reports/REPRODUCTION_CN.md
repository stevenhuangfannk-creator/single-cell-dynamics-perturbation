# 总复现报告与未完成清单 / Reproduction and remaining work

最终核查日期：2026-10-12。三个方法的已执行范围、38张真实图件、单图重绘、中文教程及数据/环境记录已整理交付。**总体科学状态为PARTIAL**：软件运行和结果追溯通过；KO排名不稳定、OE为第三方Experimental、scTour A发育趋势失败，且S2超参数未核实，因此不宣称严格全面复现。状态定义见 [方法比较](../docs/METHOD_COMPARISON_CN.md)。

## 已执行证据 / Executed evidence

| 模块 | 真实输入与执行范围 | 结果与科学状态 |
|---|---|---|
| scTenifoldKnk | 官方负二项示例100基因×2,000细胞；PBMC3k原始counts中480个CD14+ Monocytes、600基因；SPI1、SPI1+CEBPB、9候选heat KO、负对照、三库ORA；1,200基因三seed扩展 | `RUN_OK / PARTIAL`。单KO14个非靶显著基因，联合KO11个；6组主稳定性与3组扩展均保存，排名不稳定 |
| Experimental OE | 同一官方CP-WT；SPI1正目标行出边1/1.5/2/3倍；官方R MA/DR比较；第三方原生API另行验证 | `RUN_OK / PARTIAL`。1倍负控无显著响应；1.5/2/3倍各9个非靶响应；2倍与KO共享7、KO特异7、OE特异2；倍数不是mRNA表达倍数 |
| scTour A | GSE95315原始5,454→4,007颗粒谱系细胞；11,819检出筛选基因→500HVG；PCA40、latent .2z+.8zt；A20 seed0及A90 seed0/1各400epochs | 计算`RUN_OK`、趋势`PARTIAL`。A20阶段中位数Spearman=.1、非单调；两次A90也失败。原始counts、标签及官方API审查未发现已证实调用错误，根因UNKNOWN |
| scTour B | GSE104323原始24,185→15,174目标细胞；15,015筛选基因→2,000HVG；PCA30、latent .4z+.6zt；八比例各400epochs、seed0、全部细胞推断 | 计算`RUN_OK`、科学`PARTIAL`。B20/B95相关0.981548；B1/B5阶段非单调结果保留；高全局相关不能代替逐阶段趋势或实验验证 |

PBMC标签按barcode取自处理对象，网络输入为原始UMI counts。GEO直接读取原始表达表及metadata。URL、SHA256、gzip CRC、筛选/排除规则见 [数据登记](../data_registry/README_CN.md)；方法论文和官方来源见 [参考文献](../docs/REFERENCES.md)。官方模拟和tiny接口fixture没有冒充真实生物学结果。

详见 [KO报告](../modules/01_scTenifoldKnk/reports/REPRODUCTION_CN.md)、[OE报告](../modules/02_scTenifoldOE/reports/REPRODUCTION_CN.md)、[scTour逐面板报告](../modules/03_scTour/reports/REPRODUCTION_CN.md)、[A20独立审计](../modules/03_scTour/reports/A20_independent_api_audit.json)及[B八比例独立审计](../modules/03_scTour/reports/B_subsampling_independent_audit.json)。

## scTour训练与抽样 / Training and subsampling

预登记A20为a来源，B20为b–g/i来源，B95仅作h的相关参考；没有按结果好坏挑模型替换A失败。每个正式模型均400epochs，模型保存重载的ptime最大误差0，完整参数、loss、训练/验证ID和全细胞推断见 [运行表](../modules/03_scTour/results/run_comparison.tsv)及各运行目录。

| B训练比例 | Epochs | 训练秒数 | 相对95%的同序全细胞Spearman | 官方反向/手动反向；阶段 |
|---|---:|---:|---:|---|
| 1% | 400 | 349.27 | 0.901819 | False/False；非单调 |
| 5% | 400 | 655.47 | 0.926148 | False/False；非单调 |
| 10% | 400 | 690.24 | 0.958580 | True/False；严格递增 |
| 20% | 400 | 2032.95 | 0.981548 | False/False；严格递增 |
| 40% | 400 | 1731.80 | 0.992942 | False/False；严格递增 |
| 60% | 400 | 2147.98 | 0.996251 | False/False；严格递增 |
| 80% | 400 | 3420.00 | 0.996549 | False/False；严格递增 |
| 95% | 400 | 1518.83 | 1.000000 | True/False；严格递增 |

最小相关系数为0.901819。每个比例仅一个seed，细胞和抽样模型不充当生物学重复，没有伪造误差条。当前八比例网格覆盖1–95%，不能称已核对原文每个抽样间隔。

方向先保留官方基于基因检出数的自动校正；只有可靠真实阶段严格递减且相关为负时才额外使用官方reverse_time和VF反向显示。本轮正式B网格无额外手动反转；A的非单调关系不能用整体反向修复。`ptime_before_manual`与`rho_before_alignment`均已经包含官方自动定向，仅早于本仓库额外手动校正。VF使用完整5D与官方邻域投影，未取前两维冒充UMAP运动。

## 图件与复建 / Figures and reconstruction

| 图件组 | 已交付真实图件 |
|---|---|
| KO | 11张：网络、单/多KO响应、三库富集、候选扰动图及稳定性 |
| Experimental OE | 10张：网络增强、剂量、KO/OE共享/特异与通路比较；各图注明第三方来源 |
| scTour主图 | a–i九张独立图及完整Figure2组合，共10张；a六幅横图，下部b/c、d/e、f/g、h/i四列 |
| scTour诊断 | A90 seed0的a、B95的b–g，共7张；保留参考及异常，不替代20%主展示 |

合计38张，均为PNG600dpi/PDF/SVG。图件索引、含义、数据、分析/绘图脚本及独立重绘命令见 [画廊](../gallery/FIGURE_INDEX_CN.md)。完整组合见 [scTour Figure 2](../modules/03_scTour/figures/scTour_Figure2_Reproduction.png)。原图比例和配色作为参照，实际UMAP几何、异常时间和分支未经人为移动或手画修正。

KO/OE读取已提交网络CSV与结果表进行单图重绘；scTour读取本地inferred.h5ad和抽样结果，不重训。原始counts、RDS、AnnData及模型权重由.gitignore排除；首次克隆须按教程下载/训练重建，不能从图像恢复模型。每个模块教程含顺序、参数、输出、环境与限制，不是任意数据一键生成固定论文图的承诺。

## 验证与新数据 / Validation and new inputs

本机共同CI检查26项通过，见 [tests_delivery.log](tests_delivery.log)。scTour独立环境12项必要回归全部通过、无skip，见 [tests_final.log](../modules/03_scTour/reports/tests_final.log)：真实A20保存对象、Smoke checkpoint恢复身份边界、MEX字符串ID/metadata对齐、缓存counts/注释、h八行及正负相关、未完成h/i拒绝与真实比例校验。软件fixture仅验证接口与拒绝边界。

新增独立组合重绘入口`--panel composite`，缺少完整a–i正式输入会被拒绝；该边界的额外1项回归通过，见 [test_composite_boundary.log](../modules/03_scTour/reports/test_composite_boundary.log)，未重训模型。B原始2,000HVG × 15,174细胞的30,348,000个count项与prepared逐值一致；八模型的counts、身份和注释独立分块核验通过，证据见B八比例审计。

最终 [仓库完整性审计](repository_integrity_audit.json)核对网络不变量、三格式文件、PNG600DPI和本地文档链接；[绘图预检](figure_preflight.json)包含三个绘图源码及全部38个PDF的≥5pt实际文字检查。源码WARN的原因与逐面板/组合目视审查记录见各模块报告；自动检查不替代科学准确性评价。

[共享输入验证](INPUT_VALIDATION_CN.md)覆盖h5ad/10x MEX的counts、ID、轴向和metadata；现代AnnData nullable-string HDF5不能被旧scTour环境读取的问题已实际复现，并由共享导出MEX再旧环境准备的桥解决。公开A的4,007细胞/500HVG跨环境准备成功；`001`前导零、字面`NA`、逆序metadata及counts逐值回归通过。通用KO的80×150真实PBMC tiny单/联合KO入口通过，仅作为接口验证。真实Seurat对象尚未验证。

OE的WT/KO CSV两侧也发现并修复了数字ID/字面NA改写，DR靶点排除同步保持字符串身份。三个回归在两套Python环境通过，真实R及Python→R读取保持身份/权重；600×600既有公开矩阵逐元素不变。修复前失败与实际验证见 [OE网络CSV审计](../modules/02_scTenifoldOE/reports/NETWORK_CSV_INPUT_QA_CN.md)。没有重拟合网络或修改MA/DR算法。

PI/Healthy只提供 [格式与设计审查](../docs/INPUT_AND_PI_CN.md)，没有执行用户数据正式分析。未来必须先评估同细胞类型/可信谱系、全转录组QC、目标表达检出、供者与批次设计及资源；PI与Healthy不能直接解释为时间先后。

## 环境与资源 / Environments and resources

| 用途 | 实际版本与执行 |
|---|---|
| 共享输入 | Python3.13.15、NumPy2.5.4、pandas3.0.6；[shared-requirements.txt](../environments/shared-requirements.txt)用于新安装，shared.lock.txt含本机路径仅供追溯 |
| KO与OE官方R比较 | R4.3.1、scTenifoldKnk2.0.0、scTenifoldNet1.4.3；固定包版本/来源，单核运行；全部R依赖尚非容器级锁定 |
| OE正式Python | Python3.11.11、NumPy1.26.4、pandas2.2.3、SciPy1.13.1；固定virt_overexpression与virt_knock源码revision；CPU计算；pandas3上游只读数组错误已真实复现，wrapper明确拒绝 |
| scTour | Python3.11.11、sctour1.0.0、Scanpy1.9.8、AnnData0.10.8、torch2.5.1+cu121；CUDA编码/解码、官方CPU ODE；安装脚本依赖检查通过 |

硬件Windows11、i9-13900HX、32逻辑CPU、15.69GiB RAM、RTX4070 Laptop GPU8GB，见 [hardware.json](hardware.json)。11个正式模型共4,400epochs，训练计时合计3.84小时。训练时间不含全部数据下载、预处理及绘图；系统并行负载和内存压力影响明显，不能作为独占硬件benchmark。

KO主流程34.16秒，1,200基因三seed55.61秒，未记录主KO独立进程峰值RAM。OE含导入主比较54.56秒，Python与R子进程抽样RSS加总峰750.6MiB；共享页可重复计数、短峰可漏。scTour早期没有进程/CUDA峰字段的运行保持null，后期分别记录Windows peak working set与Torch peak allocated。整张GPU监控包含其他应用，不称模型独占显存峰。逐模型、抽样口径和原始记录见 [资源报告](RUN_RESOURCES_CN.md)。

## GitHub交付 / GitHub delivery

仓库为 [single-cell-dynamics-perturbation](https://github.com/stevenhuangfannk-creator/single-cell-dynamics-perturbation)，Private。代码、真实结果、图件和日志分阶段提交。早期网络失败时使用REST Git objects；普通Git已恢复，main跟踪origin/main，原本地验证历史保留在local-validation-history-20261012，没有删除或覆盖其他项目。

[交付快照](github_delivery.json)证明其中明确记录的SHA、Private设置、远端/本地相同文件树及该SHA的Actions。快照提交本身会产生新SHA；最新HEAD再以只读--no-save核查，避免审计文件递归修改自身。CI不下载GEO、不重训模型，也不能证明生物学有效性。

## 未完成与失败证据 / Remaining limits and failures

| 项目 | 实际边界与后续所需证据 |
|---|---|
| KO生物学稳定性 | 主非靶距离跨seed/抽样Spearman约.073–.251；1,200基因面板跨seed约.186–.219。结果仅探索性，需独立供者、充分网络信息与真实实验 |
| OE可靠性 | 不存在已核实的官方scTenifoldOE；第三方Experimental。固定WT剂量排名稳定，跨WT/供者和真实OE未验证；缺独立LICENSE文件，许可材料待澄清，未复制其源码入库 |
| scTour A趋势 | A20及两个A90正式模型均非单调；排除counts/标签错位，根因未知，完整失败保留 |
| 论文严格参数 | S2多来源下载未获得有效文件，实际尝试日志保留；batch_size1024和重建权重各.5为官方API默认，seurat_v3与原文严格一致性未确认；原UMAP参数/权重/全部抽样ID未获得 |
| 泛化与更大规模 | 大规模全转录组KO、跨供者benchmark、真实Seurat对象、用户PI/Healthy正式分析未执行；通用接口通过不等于全部新数据科学泛化通过 |
| 环境范围 | 完整R容器锁定、全部正式CPU重训未验证；已验证Windows/GPU步骤，不扩称其他平台全流程验证 |

KO的FC和无方向网络距离不是表达上/下调；ORA按实际测量背景校正，疾病通路名不是供者诊断。scTour可用counts而不要求spliced/unspliced，但静态转录组时间/VF预测不能证明真实谱系转化。seed/细胞抽样不是生物学重复。

自检：真实counts与软件fixture分开；所有失败模型保留；同WT比较一致；负控、测量背景与多重校正明确；38图均可追溯并独立重绘；八比例模型均实际400epochs；未验证工作和科学状态没有标为REPRODUCED。
