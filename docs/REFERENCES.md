# 来源与参考文献 / Sources

| 方法/数据 | 原始来源 | 用途 |
|---|---|---|
| scTenifoldKnk | Osorio et al. (2022), Patterns, [DOI 10.1016/j.patter.2022.100434](https://doi.org/10.1016/j.patter.2022.100434)；[官方源码](https://github.com/cailab-tamu/scTenifoldKnk) | KO 方法及当前 API 核查 |
| scTenifoldNet | [官方源码](https://github.com/cailab-tamu/scTenifoldNet) | 网络、tensor 与 alignment 组件 |
| Experimental OE | [virt_overexpression](https://github.com/virchow-zhang/virt_overexpression) | 第三方算法审查与运行；不称为官方 scTenifoldOE |
| scTour | Li Q. (2023), Genome Biology, [DOI 10.1186/s13059-023-02988-9](https://doi.org/10.1186/s13059-023-02988-9)；[官方源码](https://github.com/LiQian-XC/sctour) | Figure 2、Methods、真实训练与推断 |
| scTour 教程 | [基础推断](https://sctour.readthedocs.io/en/latest/notebook/scTour_inference_basic.html)、[方向调整](https://sctour.readthedocs.io/en/latest/notebook/scTour_inference_PostInference_adjustment.html)、[向量场 API](https://sctour.readthedocs.io/en/latest/api_vf.html) | API 和绘图流程；教程配置与论文配置分开 |
| GSE95315 | [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE95315) | scTour Figure 2a 的齿状回数据 |
| GSE104323 | [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE104323) | scTour Figure 2b–i 的齿状回数据 |
| 两套齿状回数据原研究 | Hochgerner H, Zeisel A, Lönnerberg P, Linnarsson S. (2018), Nature Neuroscience 21:290–299. [DOI 10.1038/s41593-017-0056-2](https://doi.org/10.1038/s41593-017-0056-2) | 原始研究的细胞注释与发育背景；由 scTour 正文参考文献核对 |

实际下载 URL、字节数、SHA-256、筛选数量、软件 SHA 和 API 审查结果见 data_registry 与模块报告。文献链接与方法参数来源应和真实数据登记一起读取，不能仅依据本文的名称推断完成状态。
