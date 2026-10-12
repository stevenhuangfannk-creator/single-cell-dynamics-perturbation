# 真实运行报告 / Execution report

核查日期：2026-10-12。模块计算 **RUN_OK**，科学验证 **Experimental / PARTIAL**。主数据与全部图均来自真实 PBMC counts 及计算输出，没有模拟数据代替生物学结果。小型三基因矩阵仅用于单元测试，不进入画廊或主结果。

## 已实际执行 / Executed

1. 核查官方包/API、第三方源码、版本、许可、提交、维护与测试情况，记录为 [SOURCE_AUDIT_CN.md](SOURCE_AUDIT_CN.md)。未确认官方独立 OE 软件。
2. 安装并实际调用 `virt_overexpression 0.1.0`。README 高层 API 使用真实 600 genes × 480 cells，3 个网络/每网 300 个细胞，在 CPU 上完成；分析函数耗时约 1.84 秒，输出见 `../results/native_api/`。这是公开 API 执行验证，不是无参考数值的“严格论文复现”。
3. 从模块 01 复用官方 CP tensor WT 与 SPI1 KO。SPI1 在 480 个 CD14+ Monocytes 中检出率 63.33%，567 条非零出边中有 284 条正边。本模块没有额外翻转 WT。
4. 执行强度 1/1.5/2/3；始终使用真实上游 Level 1 正目标行增强，无 Level 2/3。KO 与各 OE 使用同一官方 R MA/DR 函数、d=2、seed=1、nCores=1、empiricalNull=FALSE。KO 重算的 distance 与模块 01 正式单 KO 表在数值容差内一致。
5. 完成背景限定 GO、KEGG、Reactome 超几何 ORA，生成全部 10 张真实图的 PNG/PDF/SVG、图件 manifest 和单图重绘入口。

## 结果 / Results

| 条件 | 修改边数 | 邻接差异 Frobenius norm | 非靶 FDR<0.05 基因数 |
|---|---:|---:|---:|
| KO | 567 | 0.43454 | 14 |
| OE 1，负对照 | 0 | 0 | 0 |
| OE 1.5 | 284 | 0.13992 | 9 |
| OE 2 | 284 | 0.27984 | 9 |
| OE 3 | 284 | 0.55967 | 9 |

负对照所有 p.adj=1，manifold 距离仅浮点残差（<1e-12）。OE 1.5 与 OE 2、OE 3 与 OE 2 的无方向距离 Spearman 均 >0.999999，前 20 Jaccard=1。这个结论只覆盖同一 WT 的强度变化，不覆盖重新拟合 WT 或供体变异。

2 倍 OE 与 KO 的距离 Spearman≈0.487。7 个共同响应为 HLA-DPA1、HLA-DPB1、HLA-DQB1、HLA-DRB1、HLA-DRB5、IFITM3、PRELID1；7 个 KO 特异为 CTSD、DYNLT1、HNRNPC、NAP1L1、PPT1、RBP7、S100A12；2 个 OE 特异为 CALM2、TIMP1。差异集合仅定义为当前阈值下是否显著，不能推出两个实验中的基因作用相反。

富集背景是 599 个非靶网络基因。每个集合与数据库分别对所有背景内大小 3–300 的通路 BH 校正，包括零命中（p=1）通路：GO BP 959 项、KEGG 201 项、Reactome 647 项。共同集合有 GO/KEGG/Reactome 显著项 14/24/7，KO 特异与 OE 特异在所有库中都无显著项。HLA-II 模块造成多个免疫/疾病名称通路重复命中，不能把通路名称当供体诊断或已证明的机制。原始表 `pathway_enrichment.csv` 保留完整名称和命中基因。

## 环境、耗时与资源 / Runtime and resources

最终分析/绘图环境：Python 3.11.11、numpy 1.26.4、pandas 2.2.3、scipy 1.13.1、torch 2.5.1+cu121（本模块明确使用 CPU）、matplotlib 3.8.4、networkx 3.6.1；R 4.3.1、scTenifoldNet 1.4.3、scTenifoldKnk 2.0.0。精确依赖见 `../../../environments/oe-analysis.lock.txt`；跨机器 CPU 安装入口见 `oe-requirements.txt`。

在 32 逻辑 CPU、15.69 GiB RAM 主机上串行调用 R，主比较记录分析时间与外部监控时间两种口径：包含解释器导入的复跑耗时 54.56 秒，Python + 当前 R 子进程的抽样 RSS 总峰值约 750.6 MiB。采样间隔 0.1 秒，RSS 加总可能重复计算共享页，也可能漏掉短峰，具体见 `resource_usage.json`。默认限制最多 1000 网络基因，需要资源评估后显式修改 `--max-genes`。

## 校验 / Validation

`unit_tests.log` 记录 10 项检查，9 项通过，1 项只在 pandas 3 下适用的 guard 测试按环境跳过。涵盖：输入轴/ID/非有限值/靶点、拒绝 fallback、真实上游正边修改、负对照、所有条件有限结果、KO 与官方结果一致、富集背景与零命中通路、10 图和数据文件路径。

在 pandas 3.0.6 中独立复现了上游 `np.fill_diagonal(oe.values,0)` 的只读错误，完整错误记录见 `pandas3_upstream_failure.log`；wrapper 对 pandas 3 明确拒绝，另有兼容性 guard 测试。最终正式分析使用 pandas 2.2.3，未修补第三方源码。README 从建环境到分析、富集、测试、重绘均显式使用同一 `$Py`。

## 尚未完成 / Remaining limitations

未进行独立真实 OE 实验 benchmark、跨供体或跨随机 WT 的 OE 可靠性验证；第三方 LICENSE 材料仍不完整。第三方原生 bootstrap mean 网络与官方 CP tensor 网络不等价；原生 API 与主比较分别保存。OE 没有方向预测，边权倍数不是 mRNA 倍数，网络显著性不是传统差异表达显著性。这些限制保留为 PARTIAL，不影响已完成的软件运行事实。
