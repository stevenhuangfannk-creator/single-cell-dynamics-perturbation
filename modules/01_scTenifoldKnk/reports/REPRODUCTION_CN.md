# scTenifoldKnk 复现记录 / Reproduction report

执行日期 2026-10-12。软件运行 `RUN_OK`，方法展示与公开数据验证 `PARTIAL`。全部结果为本次真实执行；未将官方负二项模拟 counts 称为真实样本，未用图件替代未完成实验验证。

## 完成标准与结果 / Acceptance and evidence

| 验收项 | 证据 | 状态 |
|---|---|---|
| 官方源码与2.0 API核实 | CRAN 2.0.0，Net1.4.3；源码 SHA / tarball校验已登记；读取 README、NEWS、R源及man page | RUN_OK |
| 原始论文核实 | Europe PMC full-text XML，PMC9058914；检查WT网络、行KO、quasi-MA、χ²及富集过程 | RUN_OK |
| 官方示例 | README相同seed1负二项100genes×2000cells、ng10KO；10nets、500cells/net、K3、d2 | RUN_OK |
| 公开counts来源与注释对应 | 原始10x MEX32738×2700，独立processed对象仅供barcode标签，两个SHA均官方registry一致 | RUN_OK |
| Single / Multi KO | 共享WT；SPI1行和SPI1+CEBPB行归零；完整MA/DR/网络/RDS输出 | RUN_OK |
| 候选逐个KO | 2.0官方API `transcriptomeWide=TRUE, gKO=9 qualified candidates`，heat9×600；非真实全转录组 | RUN_OK / PARTIAL |
| 方向/扰动图 | 官方方向函数；library regression sensitivity；官方perturbationMap内部参考，非独立验证 | RUN_OK / PARTIAL |
| ORA | GO/KEGG/Reactome GMT来源校验；600测量背景；全部term表 | RUN_OK |
| 稳定性 | 600genes的6组、1200genes/10nets的3组，全部保存；排名一致性差 | PARTIAL |
| 软件负控 | WT vs WT，600个基因全部p=1 | RUN_OK |
| 图件 | 11张白底600dpiPNG及PDF/SVG，独立重绘、manifest、字体及视觉QA | RUN_OK |
| 真实KO实验预测准确性 | 没有匹配SPI1真实KO实验/独立供者验证，不声称已证实 | PARTIAL |

原论文 2022 使用 SERGIO 模块模拟和 Nkx2-1/Trem2/Hnf4a/Hnf4g 等数据、以及约5271microglia的系统KO。本仓库官方负二项README示例**不是SERGIO原论文图2复现**，PBMC也不是原文原始数据；称为“官方示例与公开数据方法验证”。2.0方向与heat扩展在2026新增，不能归为2022论文已验证全部功能。

## 数据与运行 / Data and execution

PBMC3k annotations：2638cells中480CD14+Monocytes，原始counts重新QC后480保留；95重复symbol求和合并。基因≥5%检出、排除ribo/mito、按logCPM variance取600并保留符合相同质检的先验调控候选；合格背景3363gene。该600panel及1200扩展不构成全转录组网络。因library门槛在原始完整counts上已检查，R分析不再对缩减panel重新QC。

主目标SPI1检测63.33%（均值1.3604UMI），CEBPB66.46%、IRF816.46%。NFKB1仅4.79%检出，按5%门槛排除。其他9个候选全部写入审计表，不根据响应统计选择。SPI1WT非零出边567，WT总非零边319727；CP重构使原来PCR阈值后的网络重新变稠密，不能把最终边数误称为仍然只留top10%。

官方示例QC：2000cells → library>=30保留1914 → outlier保留1907 → mito ratio<=0.1保留1000；100genes保留。ng10KO得到6/100 FDR<0.05（含ng10），同seed全部100heat KO得到相同WT；初次0.84s，增加QC维度记录后重复执行0.92s。参考README未给固定数值金标准，因此核对的是参数、维度、网络归零逻辑、可重复seed及官方统计公式，不宣称与未发布数值逐位一致。官方模拟与全部真实运行均单核。

PBMC主配置：8nets×400cells、nComp3、q0.9、lambda0、CP K3、maxIter1000、tol1e-5、MA2、seed1、方向回归library size。结果SPI1之外14个FDR<0.05，联合KO排除SPI1/CEBPB后11个。主排名包括HLA-DRB1、HLA-DPB1、HLA-DPA1、HLA-DRB5，反映同一seed下的网络扰动；不能据此推定SPI1独立控制特定通路。

主CP日志 `norm explained: 10.8%`，官方示例8.5%；这说明rank3重构保留的网络结构有限，**不是10.8%的生物学方差或预测准确率**。源码指标是CP分解的norm解释比例，不能将“程序complete”写成已达到1e-5重构误差。保留原始日志。扩大1200panel/10nets仍约9%且排名不稳定，不以更改参数美化主结论。

## 稳定性与科学边界 / Robustness

基准为600gene/400cells/seed1，比较非target距离Spearman与top30 Jaccard；不是供者级生物学重复。基准自比较1仅为技术检验，不能算作独立稳定性证据。

| cells/network | seed | Spearman | Top30 Jaccard | 显著non-target |
|---:|---:|---:|---:|---:|
|240|1|0.157|0.034|12|
|240|2|0.073|0.000|11|
|240|3|0.090|0.000|23|
|400|1|1.000|1.000|14|
|400|2|0.152|0.034|18|
|400|3|0.251|0.277|10|

预登记更宽panel敏感性：同一480cells、1200genes、10nets、400cells/net、K3、d2、seed1/2/3。三次耗时19.17/16.46/19.98s；与600panel基准共有599非target基因，Spearman0.352/0.148/0.098、top30Jaccard0.034/0/0。1200panel内部种子两两Spearman0.186–0.219。扩大panel与恢复默认10net并未使排名稳定；未删除失败seed，未换主分析。

因此，当前图04/05/07的显著基因与通路是探索性可复算输出，不是稳健生物学发现。缺少真实KO、供者重复及匹配表型；`FC` 不是表达倍数，direction仅为WT相关heat扩散预测。共表达与单细胞测序深度、细胞状态、ambient RNA可混杂；理论χ²及BH不保证模型错误/相关基因下的生物学FDR。还未评估empirical null、更多供者、约3000gene或充分细胞网络，不应把这些尚未运行项目写成成功。

## 富集与Perturbation Map / Enrichment and map

三个Enrichr GMT分别为GO_BP2025、KEGG2021Human、Reactome2022，term与背景交集大小3–300。每contrast query是BH<0.05排除target；600gene背景，不用全人类基因背景。hypergeometric upper tail；在每library及contrast中对**所有**eligible terms（包括zero overlap）BH。两contrast共3626term记录，82行BH<0.05；一条term出现多个contrast不代表独立证据。疾病名称例如Asthma仅表示共享HLA基因集富集，不能解释为供者患此病或KO造成此病。

Figure11确实读取官方`perturbationMap(plot=FALSE)`输出。reference是相同模型SPI1的signed heat profile（自身target清零），只测试API与内部轮廓相似性；**不是疾病logFC signature、独立阳性对照或泛化验证**。9candidate map不等于全转录组筛查，不能据相似性宣称疾病因果靶点。

## 资源、错误与后续 / Resources and limits

R4.3.1 Windows11，单核BLAS/OMP。600gene×600gene×8net原始double tensor约23.04MB；1200×1200×10约115.2MB，不含副本/MA中间对象。主全部单/多/screen/6stability/负控耗时34.16s；额外1200三seed55.61s。无独立进程级peak RAM精确记录，因此不能把raw tensor大小称为实际峰值；资源由项目总资源日志另记录。并行scTour训练期间监测available RAM，敏感性执行完成后释放R进程；后续script每seed前设置0.8GiB可用RAM门槛。

安装初次源编译rjson失败，补R4.3 binary；Knk加载时缺plyr，补binary后源码安装成功。`git clone` / codeload偶发网络连接失败，官方CRAN tarball、GitHub REST/raw与EuropePMC可访问；已用校验登记而未使用未知镜像。R启动提示部分locale=C.UTF-8设置失败，运行实际locale=C，未影响counts/数值；保留日志不隐藏。

复核config数值读取时重复运行600gene流程，在完成single/multi/screen及部分seed后因系统availableRAM降至约0.3GiB终止。初次完整运行的summary、pooled distance与完整日志保留；检查六个单seed DR中599非target distance与初次pooled matrix逐值相同（最大差0），没有用部分重复替换完整稳定性结论。

通用`run_custom.R`独立tiny API smoke现已完成：真实PBMC保存panel前80gene、前150CD14cells，6370nonzero；目标SPI1/CEBPB检测66%/64%，3net×100cells、seed7、MA2。单SPI1与联合KO分别成功，80行DR距离finite、完整WT/KO/MA/RDS保存，非target网络行不变。联合计算0.45s。此前统一0.8GiB门槛两次正确阻止运行；经规模评估，将≤100gene且≤200cell tiny入口门槛设0.25GiB，其余仍0.8GiB，真实包加载前后内存数值记录于`custom_smoke_summary.json`，未stub资源检查。联合after-packages available0.406GiB，rawtensor仅153600bytes；直接Rterm单KO进程每50ms采样观察峰值RSS270700544bytes（约258MiB，采样峰值而非无限精度瞬时峰值）。全部tiny运行位于独立work目录，仅用于验证API，不输出额外生物学图，不替换正式结论，也未运行PI数据。

未完成的是：真实PBMC全转录组KO、对匹配真实KO的生物学准确性评估、独立供者重复与稳定目标确认、完整容器级R dependency lock、用户PI/Healthy正式分析。当前480细胞少于官方约500建议，600gene/8nets是资源受控验证，不能直接用于未审查的新数据研究结论。

自检：原始counts与标签分离、target存在/表达/边信息已审查、single/multi矩阵行操作正确、WT复用一致、负控allp1、结果可读取、directionScore target NA按官方定义而保留、全部seed未筛选、富集background正确、图件数据可追踪、11图可单图修改。尚未达到“稳健可实验验证的SPI1机制结论”，故状态保持PARTIAL。
