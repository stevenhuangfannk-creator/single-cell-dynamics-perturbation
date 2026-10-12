# scTour Figure 2 复现报告 / Reproduction report

状态：`PARTIAL`。2026-10-12已完成真实数据预处理、官方模型训练、对象保存加载及Figure 2 a–i和组合图；严格原论文所有配置及生物学趋势尚不能全部确认。正式训练共11个模型、4,400 epochs，另有200细胞/3epochs的Smoke Test，没有用演示结果替代论文复现。

## 来源与实现 / Sources and implementation

已读取 [原始论文与 Methods](https://doi.org/10.1186/s13059-023-02988-9)、[PMC Figure 2](https://pmc.ncbi.nlm.nih.gov/articles/PMC10290357/#Fig2)、[官方基本教程](https://sctour.readthedocs.io/en/latest/notebook/scTour_inference_basic.html)、[官方方向修正教程](https://sctour.readthedocs.io/en/latest/notebook/scTour_inference_PostInference_adjustment.html) 及 [向量场源码](https://github.com/LiQian-XC/sctour/blob/main/sctour/vector_field.py)。原图顶部 a 是六幅横图；下部依次为 b/c、d/e、f/g、h/i 四列。原图参考只存工作缓存。

实际使用 PyPI `sctour==1.0.0`；安装的 `train.py` 和 `model.py` 与本轮下载的官方 main（commit `4d26abb936b7872706e6540f38bab6251159a475`）逐文件一致。Python 3.11.11、torch 2.5.1+cu121、Scanpy 1.9.8、AnnData 0.10.8、NumPy 1.26.4、SciPy 1.13.1。固定 torch 2.5.1 保留官方模型加载语义，依赖兼容性检查通过；完整锁定见 `environments/sctour-lock.txt`。

Supplementary Table S2 的初次下载返回HTML，未当作有效XLSX读取。`batch_size=1024`、`alpha_recon_lec=alpha_recon_lode=0.5`采用当前官方API默认值；不能把它们称为已经逐字段核对的论文S2参数。HVG算法使用 `seurat_v3`，其严格原文一致性同样未确认。

NVIDIA RTX 4070 Laptop GPU（8,188 MiB），i9-13900HX、32 logical CPUs、15.7 GiB RAM。编码器/解码器按官方实现使用 CUDA；ODE 求解按原软件在 CPU 运行。每个实际模型保存耗时与 loss；后续模型另存 Windows 进程 peak working set、CUDA peak allocated 和实际 device。系统有并行应用，耗时不应视为独占硬件基准。

## 原始数据审计 / Raw-data audit

Supplementary Table S2 的追加获取仍失败，PMC、EuropePMC及出版社附件的实际尝试见 [supplement_access_attempts.json](supplement_access_attempts.json)。没有将HTML或无效XLSX当作已读附件。

A 的原始 5,454 细胞恢复 4,007 个颗粒谱系细胞；14,545 基因经检出 ≥20 细胞规则剩 11,819，再选500 HVGs。B 的原始24,185细胞恢复15,174个目标细胞；27,933基因过滤后15,015，再选2,000 HVGs。两组与论文细胞数完全吻合，元数据按ID匹配，不以行号猜测。Gzip CRC、SHA256、大小和URL见 `data_registry/sctour_*_preprocessing.json`。

| 数据 | 类别与数量 |
|---|---|
| A | nIPC 88；Neuroblast_1 97；Neuroblast_2 777；Granule-immature 1,333；Granule-mature 1,712 |
| B | nIPC 821（含原始 nIPC-perin）；Neuroblast 1,381；Immature-Pyr 4,520；Immature-GC 2,419；GC-juv 3,420；GC-adult 2,613 |

在原始counts上用 `seurat_v3` 选HVG，QC对象与训练对象独立保存；训练 `.X` 不被归一化/PCA绘图覆盖。坐标对象保留原细胞次序。批次定义为细胞ID中的原始10x文库前缀（如`10X79_1`），原始年龄及细胞类别同时保留。

训练队列结束后，独立审计用固定小块读取原始B gzip和HDF5，没有载入完整矩阵或模型。2,000HVG×15,174细胞的30,348,000个count项（5,081,261个非零项）的有序逐基因摘要与prepared完全相同；八个模型的`.X`及`layers['counts']`的CSR `data/indices/indptr`直接逐块比较相同。细胞/基因ID及顺序、原始标签/年龄、QC基因检出数和模型方向标记均一致；导出TSV与HDF5的时间差小于1e-7，见 [B_subsampling_independent_audit.json](B_subsampling_independent_audit.json)。原始文件SHA256独立重算且gzip读至末尾无CRC错误；审计不重新训练或加载checkpoint。

## 已完成运行与异常 / Executed runs

精确运行清单见 `results/run_comparison.tsv`，每个目录同时保存 `config_used.json`、环境、loss、实际训练/验证细胞ID及全细胞推断表。耗时、实际模型device及已记录的内存峰值汇总于 `results/resource_comparison.tsv`；早期模型未采集的字段留空。只导出完成了 `run_summary.json` 的目录。尚未完成的训练不会被列为成功。

最初 A 使用软件默认90%训练，seed0与seed1各400epochs，分别耗时407.98与359.71秒。训练/验证最终loss为123.72/121.03与123.06/125.64；模型保存重新载入的ptime最大误差均0。

两次 A90 模型均未恢复正确发育顺序：seed0阶段中位数为 nIPC 0.0028、Neuroblast_1 0.9750、Neuroblast_2 0.8443、Granule-immature 0.5829、Granule-mature 0.4583；与已知顺序相关为0。nIPC和后续阶段呈非单调关系，不能用整体反转修复。seed1重复了同一现象。失败证据保留在结果目录与 `figures/diagnostics/A90_seed0_Figure2a.*`。

针对异常，执行独立 source-to-input 审计：原gzip中的4,007个细胞标签及500个HVG全部2,003,500个count逐值与训练h5ad完全一致，见 `results/A_source_identity_audit.json`；排除细胞标签与矩阵列错位。源码、QC/HVG输入和官方API对照没有发现调用错误。没有证实收敛、软件版本或HVG算法中的哪项导致该偏差，不给未经验证的根因。

论文Results使用20%训练的表述，Methods没有给出每个主面板的完整训练比例；因此预登记追加 A20/500HVG/400epochs/seed0，并独立记录。A20作为贴近论文训练设计的分析，评价取决于实际结果，A90失败不删除。`paper_A.json` 默认20%，原90%配置保存为 `diagnostic_A90.json`。官方后处理教程的2,000HVG另存配置，不替代论文500HVG。

A20 已完成400epochs，训练518.47秒，loss132.65/129.20，checkpoint往返ptime误差0。实际训练/验证张量逐值与种子重放的细胞ID吻合。其阶段中位数 nIPC .0541 → Neuroblast_1 .8760 → Neuroblast_2 .9097 → Granule-immature .5163 → Granule-mature .4302，阶段Spearman=.1，仍非单调。官方gene-count自动方向标记为True，但没有额外手动反转；原结果完整用于主a面板，不选择另外一个好看的seed或HVG配置来替换失败。

B95正式400epochs耗时1518.83秒，loss600.29/594.85，模型往返ptime误差0。颗粒谱系中位数为 nIPC .0595 → Neuroblast .2789 → Immature-GC .4533 → GC-juv .8544 → GC-adult .9471，严格递增（阶段中位数Spearman=1）。Immature-Pyr .5373，没有把锥体与颗粒细胞人为排成单一成熟序列。B95是Figure2h参考模型；其b–g图保留在diagnostics而非主展示。B20按正文设计预登记为主b–g/i来源，不依赖是否比B95更好看或更吻合预期。

B20 已完成400epochs，训练2032.95秒，loss655.35/650.18，checkpoint往返ptime误差0。阶段中位数 nIPC .0404 → Neuroblast .2295 → Immature-GC .4821 → GC-juv .9237 → GC-adult .9786 严格递增；Immature-Pyr .4632为另一分支。官方与手动反向标记均False。与B95同序15,174个细胞的ptime Spearman为 .981548。该数值支持跨训练比例排序一致性，但不是独立实验验证；本轮B只有单seed。

B的1/5/10/20/40/60/80/95%八个模型均已实际完成400epochs，用全部15,174细胞推断。原队列已退出，最新版缓存核对与八行相关表重新生成成功，见 [subsampling_cache_final.log](subsampling_cache_final.log)；12个完成目录（11个正式模型和Smoke）已导出，见 [export_final.log](export_final.log)。每个比例一个seed；该网格覆盖1–95%范围，但原图h未逐行给出完整比例列表，因此不声称与原文每个抽样间隔完全一致。

| 训练比例 | 训练/验证细胞 | 训练秒数 | 相对95%的Spearman | 官方自动反向 | 额外手动反向 | 颗粒五阶段中位数 |
|---|---:|---:|---:|---|---|---|
| 1% | 152/15 | 349.27 | 0.9018189393514768 | False | False | 未严格递增，幼年高于成年 |
| 5% | 759/75 | 655.47 | 0.9261481997790174 | False | False | 未严格递增，幼年略高于成年 |
| 10% | 1,518/151 | 690.24 | 0.9585804421806908 | True | False | 严格递增 |
| 20% | 3,035/303 | 2,032.95 | 0.9815482870312411 | False | False | 严格递增 |
| 40% | 6,070/607 | 1,731.80 | 0.9929416446473271 | False | False | 严格递增 |
| 60% | 9,105/910 | 2,147.98 | 0.9962505361822267 | False | False | 严格递增 |
| 80% | 12,140/1,214 | 3,420.00 | 0.9965489058639511 | False | False | 严格递增 |
| 95% | 14,416/758 | 1,518.83 | 0.9999999999999998 | True | False | 严格递增 |

95%行是参考模型与自身的比较，软件浮点结果接近1，不是独立验证。八个模型本轮都没有额外manual reverse，因此`rho_before_alignment`与`rho_aligned`相同；10%与95%的官方automatic反向已经包含在两列时间中。

## 图件与统计解释 / Figures and interpretation

已完成的B1模型400epochs/349.27秒，loss685.22/700.14，阶段中位数Spearman=.9，但GC-juv .9913高于GC-adult .9800，终末两阶段次序未通过。该反例保留在B1摘要；不能把全局模型相关或外观相似解释为每个阶段都准确。

B5与B10也已完成400epochs，分别耗时655.47与690.24秒，模型往返误差均0。B5的GC-juv .973078略高于GC-adult .972200，仍未通过严格阶段单调性；B10的五个颗粒阶段严格递增。B40完成400epochs，1731.80秒，loss624.02/624.51；B60完成400epochs，2147.98秒，loss609.44/606.50。B40/B60的颗粒阶段均严格递增、模型往返误差0。

B80在2026-10-12实际完成400epochs，训练3419.9980254999828秒、总运行3606.848881999962秒，loss602.167525697108/609.6748388752317。进程peak working set为1,121,587,200 bytes，CUDA peak allocated为172,468,224 bytes；官方与额外手动反向均False，checkpoint往返ptime误差0。颗粒阶段中位数为 nIPC .077074 → Neuroblast .222472 → Immature-GC .334009 → GC-juv .791954 → GC-adult .971440，严格递增；Immature-Pyr .413975保留为另一分支。全部八比例已经完成，没有以中途400/400日志代替最终摘要。

`figures/figure_manifest.json` 把每个真实面板与数据、分析入口、绘图入口及三格式文件一一对应。a 的五维latent与VF来自A；b/c/d使用同一个B PCA UMAP；f/g使用同一个B latent UMAP；i使用独立20%训练模型自己的latent UMAP。所有时间色条范围0–1。e/h包含每个观测细胞，不为布局抽样删除细胞。

时间方向先保留官方基于基因检出数的自动校正；仅当真实阶段中位数严格递减且与已知顺序相关为负时，使用官方 `reverse_time`，同时官方VF绘图`reverse=True`。保存原时间、校正时间和方向标记。高维VF从未截取前两维充当UMAP位移；使用官方潜空间余弦相似度、递归邻域、时间邻域和嵌入位移进行投影。

h以95%训练为参考计算全部同序细胞的Spearman；`ptime_before_manual`和`rho_before_alignment`已经包含官方基于检出基因数的automatic方向校正，仅表示尚未执行本仓库额外的manual reverse，不是未经校正的encoder原始时间。`ptime`与`rho_aligned`是在这一基础上按可靠生物学阶段核查并在满足严格反向条件时手动校正的结果。相关系数是模型排序一致性，不是独立实验因果验证。单seed没有人为添加误差条，也不从不同细胞冒充生物学重复。

`pvalue_aligned`保留SciPy按细胞计算的软件返回值，用于追踪计算；细胞不是独立模型或生物学重复，不据此宣称统计显著性。大样本高相关时返回0可能是浮点下溢，不能将其表述为精确的`p=0`；本模块的稳健性解释以排序相关与阶段检查为主。

计算结果、发育趋势、视觉风格分别评价；严格复现不由外观决定。现代环境与未发表的UMAP参数、原始模型权重及全套训练抽样ID限制逐点数值一致；UMAP随机几何差异如实保留，不能移动细胞或手绘分化路径。

| 面板 | 原始→当前细胞数；输入 | 嵌入；主模型 | 计算结果 | 生物学趋势 | 视觉与差异 |
|---|---|---|---|---|---|
| a | A 5,454→4,007；500HVG | PCA40 UMAP + 5D latent UMAP；20%、400epochs、seed0、z/zt=.2/.8 | 真实模型与官方VF | 未恢复单调成熟顺序 | 六幅横向图、A独立palette；真实异常保留 |
| b | B 24,185→15,174；2,000HVG | PCA30 UMAP；B20、400epochs、seed0、z/zt=.4/.6 | 完整细胞注释 | 支持两类神经元谱系标签 | 原始UMAP部分亚群分离更明显，未经搬移连接 |
| c | 同b | 与b逐点相同PCA坐标 | 完整11个原始文库批次 | 批次分布供检查，未宣称已校正批次 | tab20类别编码，批次原值保留 |
| d | 同b | 与b相同PCA坐标；B20 | 全部细胞ptime范围0–1 | 颗粒阶段中位数递增 | viridis连续色条，几何差异来自真实UMAP |
| e | 同b | B20全体细胞按ptime稳定排序 | 一列一个细胞、批次及逐类型条带 | 分支与成熟状态时序吻合 | 与原图同类层叠结构；横轴标明rank而非线性time |
| f | 同b | B20五维X_TNODE UMAP | 类型及同坐标批次inset | 潜空间形成颗粒/锥体两分支 | inset置右上空白，未遮盖主轨迹 |
| g | 同b | 与f相同坐标；完整5D VF官方投影 | 黑色流线及箭头 | 主要方向符合nIPC向两个神经元分支 | 合法密度/宽度/箭头参数；无手画轨迹 |
| h | B全部15,174；2,000HVG | 1/5/10/20/40/60/80/95%各400epochs，seed0 | 完整八行相关表已生成 | 10–95%阶段递增；1/5%终末阶段错序保留 | 仅八个实际比例行，右侧显示95%参考相关 |
| i | 同b | B20自己的latent UMAP；全体细胞 | 真正20%模型，与b–g同一预登记模型 | 方向及完整分支吻合B阶段关系 | 细胞类型inset、单独留白与0–1色条 |

A/B各自类别和颜色完整定义在绘图脚本 `PALETTES`，未按截图七色硬拆不存在的类别。九个独立面板a–i及完整组合均已输出PNG600/PDF/SVG；manifest共17项，其中10个主图、7个保留的A90/B95诊断图，全部有独立重绘入口。

## 检查与未完成项 / Validation and limits

最终10个主图的文件名全集、主manifest10项及三格式文件均通过核对；全部PDF可审计文字≥5pt，PNG元数据600DPI、SVG保留文字，见 [main_figure_file_qa.json](main_figure_file_qa.json) 与 [main_figure_file_qa.log](main_figure_file_qa.log)。h及组合图的实际renderer标签边界均通过，见 [figure_layout_qa.json](figure_layout_qa.json)：e最长标签距c轴区域最小5.74 pixels，h的ylabel距f轴区域27.81 pixels，h的xlabel距i轴区域13.02 pixels，均在100DPI检查画布中测量，实际600DPI导出相同比例。

模块绘图代理与项目主代理独立目视h、完整组合预览及原分辨率c/e、f/h裁片通过；11批次共享图例、h八行及相关列、f/i的inset和i色条可读。A的非单调偏差完整显示。最终组合绘图exit0，见 [redraw_final.log](redraw_final.log)；h独立边界检查见 [redraw_h_layout.log](redraw_h_layout.log)。首次QA包装的NumPy bool JSON序列化失败已修正为原生bool，失败日志保留在 [redraw_final_layout_serialization_failure.log](redraw_final_layout_serialization_failure.log)，没有改动推断数值来解决绘图检查。

A20进程peak working set为1,324,638,208 bytes，CUDA peak allocated为48,850,432 bytes；总运行598.33秒包含推断和嵌入，训练518.47秒。B20早期运行未记录相同的peak字段，不补造资源数字。

输入验证包含负值/非有限/非整数/重复ID的拒绝测试。A90实际对象3项测试通过；B20对象及真实Smoke checkpoint边界回归合计5项通过，验证同细胞数但ID不同或counts不同的恢复被拒绝，空stage_order输出NOT_ASSESSED/null且不反转。时间范围、五维latent/VF、两个二维嵌入及保存载入检查通过，但A的生物学排序未通过。新数据准备桥在公开A counts上实际执行，生成4,007细胞/500HVG输入，没有训练PI模型。

现代共享AnnData/Pandas写出的nullable-string HDF5编码已实际复现不能被scTour旧环境读取；改用共享入口的原始MEX与metadata TSV后，公开A数据跨环境准备成功（4,007细胞/500HVG），见 `shared_mex_bridge.log` 与 `shared_mex_bridge_legacy.log`。额外软件回归发现Pandas2的`index_col=0`会在`dtype=str`下仍把`001`索引转为整数；scTour现先读取字符串再设置索引，保留前导零及字面基因名`NA`，且逆序元数据可按ID恢复，见 `mex_string_id_regression.json`。小型软件fixture不是生物学结果；完整MEX准备回归已实际通过（60×200软件counts、50HVG，逐值与字符串ID校验）。

最新必要回归共12项全部通过，无skip，见 `tests_final.log`：实际A20保存对象、真实Smoke检查点恢复边界、MEX字符串身份及按ID对齐、旧缓存注释/QC、抽样指标中的`001`/字面`NA`、h的完整八行及正负相关侧标签、未完成i/部分h拒绝、比例与细胞次序核对。没有在这些回归检查中重新训练模型。

组合图另支持`--panel composite`，从完整a–i结果只重绘组合并保留已有单图索引，避免为了修改组合布局重复导出a–g。其完整结果拒绝边界另外执行1项回归并通过（0.117秒、无skip），见 [test_composite_boundary.log](test_composite_boundary.log)；原12项日志没有改成13项。

最新图形源检查17 PASS、3 WARN、0 FAIL；WARN为配置变量的字体大小无法从静态字面值审计（已限制下限5pt并已检查全部10个主PDF实际文字）、无TIFF（用户要求PNG/PDF/SVG）及12英寸独立a宽度不同于常见183mm提交尺寸（本库用于Figure2横图阅读，非期刊已验收版式）。官方streamplot的少量masked-grid NaN线宽导致PDF初次导出失败；只将这些无定义的样式线宽置0，未改VF、有限几何或细胞坐标，保留真实流线计算。

本轮计划的正式模型和图件均已执行；未验证的精确原论文超参数、不满足发育趋势的A结果、B1/B5局部错序及单seed稳健性限制仍保留，因此保持 `PARTIAL`，不标记 `REPRODUCED`。本模块没有运行PI/Healthy用户数据的正式结论。模型从静态转录组估计连续排序，不能证明真实谱系转化；周期、终末状态和复杂根分支的限制详见论文Discussion及模块教程。
