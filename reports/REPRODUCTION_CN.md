# 总复现报告与未完成清单 / Reproduction and remaining work

核查日期：2026-10-12。本报告是交付过程中的实际结果快照，scTour 后续完成时须同步更新下表及模块报告。总体状态为 **PARTIAL**：已有真实软件运行、公开数据结果和可追踪图件，但没有达到独立实验验证或原论文逐点数值复现。`RUN_OK` 只说明指定计算与输出检查通过；状态定义见 [方法比较](../docs/METHOD_COMPARISON_CN.md)。

## 已执行证据 / Executed evidence

| 模块 | 真实输入与已执行范围 | 当前结果与科学限制 |
|---|---|---|
| scTenifoldKnk | 官方负二项示例100基因×2,000细胞；真实 PBMC3k 中480个 CD14+ Monocytes、600基因；SPI1、SPI1+CEBPB、9候选 heat KO、负对照、三库 ORA；另有1,200基因敏感性 | 软件 `RUN_OK`，科学 `PARTIAL`。SPI1 非靶显著基因14个，联合KO排除两靶点后11个；6组主稳定性与3组更宽面板运行均保存，排名不稳定 |
| Experimental OE | 复用同一官方 CP-WT，600基因×480细胞；SPI1 正目标行出边增强1/1.5/2/3倍；同WT KO/OE比较；真实第三方原生API另存 | 软件 `RUN_OK`，科学 `Experimental / PARTIAL`。1倍负控无显著响应；1.5/2/3倍各9个非靶响应；2倍与KO共同7个、KO特异7个、OE特异2个 |
| scTour A | 原始 GSE95315：5,454细胞中4,007个论文颗粒谱系细胞，11,819个检出筛选基因中选500 HVGs | A20(seed0)与A90(seed0/1)均完成400 epochs，模型往返时间误差0；A20阶段中位数Spearman=0.1且非单调，生物学趋势未通过。原始counts、标签及API调用审查没有发现已证实的调用错误，根因未知 |
| scTour B | 原始 GSE104323：24,185细胞中15,174个论文谱系细胞，15,015个筛选基因中选2,000 HVGs | B20%主模型和B95%参考模型各400 epochs已完成，全部细胞推断及模型往返检查通过；颗粒谱系阶段中位数 Spearman≈1。其他抽样比例与 h 仍在运行，整体保持 `PARTIAL` |

PBMC 网络输入来自原始 UMI counts；处理 h5ad 只提供按 barcode 对应的标签。GEO 数据直接恢复原始表及其元数据，下载 SHA256、gzip CRC、筛选规则和排除记录见 [数据登记](../data_registry/README_CN.md)。官方模拟示例单独登记，未冒充真实样本或原论文 SERGIO 图复现。

详细证据见 [KO报告](../modules/01_scTenifoldKnk/reports/REPRODUCTION_CN.md)、[OE报告](../modules/02_scTenifoldOE/reports/REPRODUCTION_CN.md)、[scTour报告](../modules/03_scTour/reports/REPRODUCTION_CN.md)、[A20独立API审查](../modules/03_scTour/reports/A20_independent_api_audit.json) 和 [scTour已完成运行表](../modules/03_scTour/results/run_comparison.tsv)。

## 图件与复建 / Figures and reconstruction

| 图件组 | 当前有实际三格式文件的图件 | 仍未完成 |
|---|---|---|
| KO | 11张；PNG600dpi、PDF、SVG，完整网络与统计表支持单图重绘 | 没有用未执行的全转录组结果补图 |
| OE | 10张；PNG600dpi、PDF、SVG，同WT剂量与KO/OE比较，页脚保留第三方来源 | 跨WT/供者稳定性与真实OE实验验证图未执行 |
| scTour主面板 | A20的a与B20的b–g、i，共8个；各有 PNG/PDF/SVG | 抽样面板 h、完整 Figure2 组合图及最终整体QA |
| scTour诊断 | A90 seed0的 a、B95的 b–g，共7个；保留实际异常和参考模型 | 诊断图不替代预登记20%主展示 |

当前合计21张KO/OE图，加15张scTour主/诊断图；已按三个 manifest 检查三格式文件存在。当前36张图的600DPI与三格式检查通过，84项仓库完整性检查无失败，全部主图与诊断PDF通过≥5pt字体检查。图库及逐图代码入口见 [图件索引](../gallery/FIGURE_INDEX_CN.md)，最终须加入h及组合图后重建，不能把当前数量当成最终数量。

KO/OE单图脚本读取已提交CSV，不重新构建网络或执行MA。scTour单图读取本地 `inferred.h5ad`、同目录配置及必要抽样结果，不重新训练。原始下载、counts、RDS、AnnData、模型权重留在被Git排除的本地工作目录；首次克隆须按模块教程下载/训练重建这些对象，不能从图像恢复模型。

## 环境与资源口径 / Environments and resources

| 用途 | 已实际使用的关键版本 |
|---|---|
| 共同输入/表格环境 | Python3.13.15、NumPy2.5.4、pandas3.0.6；新建环境用 [shared-requirements.txt](../environments/shared-requirements.txt)；`shared.lock.txt`含本机安装路径，仅作版本记录 |
| KO及OE官方比较 | R4.3.1、scTenifoldKnk2.0.0、scTenifoldNet1.4.3，R计算单核；官方包版本与源校验固定，全部R依赖尚非容器级完整锁定 |
| OE正式Python环境 | Python3.11.11、NumPy1.26.4、pandas2.2.3、SciPy1.13.1、torch2.5.1+cu121；OE计算使用CPU；实际 `virt_overexpression0.1.0`及`virt_knock0.2.0`源码revision已固定 |
| scTour | Python3.11.11、sctour1.0.0、Scanpy1.9.8、AnnData0.10.8、torch2.5.1+cu121；模型编码/解码使用CUDA，官方ODE仍在CPU |

硬件为Windows11、i9-13900HX、32逻辑CPU、约15.69 GiB RAM、RTX4070 Laptop GPU（8,188 MiB），见 [硬件登记](hardware.json)。不同环境独立，pandas3导致第三方OE只读数组错误已真实复现并由wrapper拒绝；正式OE使用pandas2.2.3，未修补第三方源码。

| 运行 | 实际记录 | 不能据此声称 |
|---|---|---|
| KO主流程与扩展 | 主单/多KO、筛查、6组稳定性及负控34.16秒；1,200基因三seed合计55.61秒 | 没有独立进程峰值RAM记录；23.04/115.2 MB原始张量估算不等于实际峰值 |
| OE主比较 | 含导入复跑54.56秒；Python及当前R子进程抽样RSS加总峰值750.6 MiB，0.1秒采样 | RSS可重复计共享页、遗漏短峰，不能当作精确独占内存基准；见 [资源记录](../modules/02_scTenifoldOE/reports/resource_usage.json) |
| scTour正式训练 | A90两seed训练407.98/359.71秒；A20训练518.47秒；B20训练2,032.95秒；B95训练1,518.83秒 | 这些为各模型记录的训练时间，不含全部下载/预处理/绘图；并行应用与系统压力影响耗时，正式CPU重训未验证 |

## 科学与新数据边界 / Scientific and reuse limits

KO的6组主运行除参考自比较外，非靶距离Spearman约0.073–0.251；1,200基因面板内部不同seed约0.186–0.219。其单seed显著基因和通路只能作为探索性假设。OE强度排名Spearman>0.999999、top20 Jaccard=1仅成立于一个固定WT，不能抵消KO的跨拟合不稳定性。没有匹配真实KO/OE实验、独立供者重复或外部表型benchmark；随机seed和细胞抽样不是生物学重复。

KO的`FC`不是表达倍数，OE边权倍数不是mRNA倍数；无方向网络距离不能推断真实上/下调。ORA按实际测量背景及所有eligible terms校正：KO背景600基因，OE背景599个非靶基因。疾病名称通路仅表示基因集重叠，不是供者诊断或因果机制。

scTour时间和潜空间VF是静态转录组模型估计，不能证明谱系转化或当作实验RNA velocity。保留官方方向与人工方向审查；缺少可靠阶段或少于两个有效阶段时应记为未评估，不人为设定发育顺序。PI/Healthy分组不能直接当作时间起点/终点，跨分支、周期和终末细胞风险仍需外部证据。

共同输入工具已验证counts、唯一ID、轴方向、metadata对齐及审计后MEX导出往返，见 [输入验证](INPUT_VALIDATION_CN.md)；[新数据说明](../docs/INPUT_AND_PI_CN.md)提供h5ad/10x/Seurat接入。generic KO及scTour准备入口是复用路径，不表示已在用户私有数据完成模型或泛化验证。Seurat脚本尚未读取真实用户对象；用户数据须先完成同细胞类型/可信谱系、全转录组QC、供者/批次设计、目标检出和资源检查。

## GitHub交付与验证 / Delivery verification

仓库为 [single-cell-dynamics-perturbation](https://github.com/stevenhuangfannk-creator/single-cell-dynamics-perturbation)，登记为Private。已按输入骨架、KO、OE、保存结果检查、审计MEX桥接及通用KO验证分阶段提交上传。早期Git传输失败时使用GitHub REST Git objects发布；原本地提交保留在 `local-validation-history-20261012`。当前Git传输已恢复，`main`跟踪`origin/main`，之后使用普通commit/push。

[交付审计](github_delivery.json)记录已同步阶段的`local_tree == remote_tree`和轻量Actions成功，说明该阶段文件树一致；这与commit SHA相同是不同检查。该记录目前对应较早已同步阶段，未覆盖仍在生成的全部scTour结果及本报告。最终须再次提交、发布并核对远端tree与Actions，保留实际SHA映射，不能把阶段检查当作最终完整交付证明。轻量CI不下载完整GEO、不重训模型，也不证明生物学有效性。

## 未完成清单 / Remaining work

| 项目 | 当前状态与完成依据 |
|---|---|
| scTour B抽样网格及主Figure2 | A20主a已完成并保留非单调结果；其余B比例仍在运行，只接受全部真实400epoch摘要、推断对象、h/组合三格式图及QA作为完成证据 |
| 最终图库、总报告、GitHub同步 | 待scTour真实结果确定后更新；重建索引、运行必要检查、核对最终文件树与Actions |
| 科学准确性与稳健性 | 尚无匹配KO/OE实验、跨供者/WT验证、可靠KO目标确认或scTour独立谱系证据 |
| 更大规模与新用户数据 | PBMC全转录组、大规模多供者分析及PI/Healthy正式结论未执行；不能直接套用当前不稳定小面板结论 |
| 环境与上游材料 | 完整容器级R锁定、旧R依赖撤下时的所有回退、全部正式CPU重训未验证；第三方OE缺少独立LICENSE文件，许可材料待澄清 |
| 泛化接入 | generic模型和真实Seurat对象仍需独立验证；固定论文/PBMC绘图套件不等于任意数据一键生成全部图 |

自检：真实counts和模拟测试分开、失败seed保留、同WT比较一致、负控及测量背景明确、已完成图可追踪、未完成计算未标成功。当前交付支持方法运行与探索性结果复算，科学状态保持`PARTIAL`。
