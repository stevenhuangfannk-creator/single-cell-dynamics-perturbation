# 方法来源与源码审查 / Source audit

核查日期：2026-10-12。结论：本模块保留项目名称 `scTenifoldOE`，实现必须标为 **Experimental**。公开检索、当前 CRAN `PACKAGES` 名单和已安装 `scTenifoldKnk 2.0.0`、`scTenifoldNet 1.4.3` 的导出函数中，未找到可确认的独立官方 `scTenifoldOE` 或 `scTenifoldOverexpression` 包/API。搜索到论文中同名的定制流程，不能据此认定正式软件存在。这里的“未找到”有明确检索范围，不声称证明全世界不存在该软件。

实际使用 [virt_overexpression](https://github.com/virchow-zhang/virt_overexpression/tree/860d0ea59dc0aeb564a286efae00f562ce84e138) 的公开函数 `virtual_overexpression`，作者为第三方 `virchow-zhang`，不是 scTenifoldKnk 原作者团队发布的官方过表达方法。

| 项目 | 核查结果 |
|---|---|
| OE 版本 / 固定 revision | `0.1.0` / `860d0ea59dc0aeb564a286efae00f562ce84e138` |
| 必需 companion | `virt_knock 0.2.0` / `53abe39616bce36aa0fb3db8d15cc0e689cb1b54` |
| 声明的成熟度 | pyproject 的 Development Status 为 Alpha |
| 创建 / 最近推送 | 2026-06-27；截至核查日仅 2 个提交，均在同日 |
| 发布 / 测试 | 无 tags、GitHub releases、tests 目录或科学回归测试 |
| CI | 仅自动 Dependency Graph；未见分析/单元测试工作流 |
| 依赖 | numpy、pandas、scipy、torch、tqdm、virt_knock；富集可选 gseapy |
| OE 许可 | pyproject / README 声明 MIT，但 LICENSE 链接缺失、GitHub license 字段 null；许可材料不完整 |
| companion 许可 | virt_knock 仓库存在 MIT LICENSE |
| 本仓库处理 | 锁定源码版本，通过下载/安装调用；不复制分发第三方源码，不把许可不完整描述成已审结 |

具体源码：[OE core.py](https://github.com/virchow-zhang/virt_overexpression/blob/860d0ea59dc0aeb564a286efae00f562ce84e138/virt_overexpression/core.py)、[pyproject.toml](https://github.com/virchow-zhang/virt_overexpression/blob/860d0ea59dc0aeb564a286efae00f562ce84e138/pyproject.toml)、[companion core.py](https://github.com/virchow-zhang/virt_knock/blob/53abe39616bce36aa0fb3db8d15cc0e689cb1b54/virt_knock/core.py)。GitHub 元数据可由 [OE API](https://api.github.com/repos/virchow-zhang/virt_overexpression) 与 [commits API](https://api.github.com/repos/virchow-zhang/virt_overexpression/commits) 重新查询。上游未来维护状态可能改变。

## 实际 OE 动作 / Actual perturbation

对靶基因 `g`，源码首先取目标行，正权重乘以 `oe_fold`，负权重保持不变，最终对角线归零。这是 Level 1。本次 SPI1 在共享 WT 的目标行有 284 条正边，真实执行始终走 Level 1。没有利用后备路径凑出扰动。

上游称目标行为“incoming edges”，但其 companion 的 `_make_one_pc_network` 将预测目标的回归系数写入行，随后高层流水线 `.T` 转置。[官方 Knk 源码](https://github.com/cailab-tamu/scTenifoldKnk/blob/master/R/scTenifoldKnk.R) 同样先转置 WT 再将靶点行归零。因此在最终 WT 中 `WT[i,j]` 是 **i 到 j 的出边**。本仓库共享官方最终 WT，不额外转置，图与报告按实际出边解释；原日志仍保留上游的 incoming 文字，不能照日志误解方向。

Level 2 在没有正行边时复制列到目标行，并放大该列的正值。这不是简单“增强原有出边”，可能改变方向与其他符号权重。Level 3 将绝对 Spearman 相关映射成正边，会把负相关也转成正权重，且没有候选相关时 `max()` 可失败。本项目默认拒绝 Level 2/3，要求先报告靶点网络信息量，再另行设计分析。调用前严格检查缺失靶点、重复 ID、轴顺序、非有限值、非零对角线和无正边，而不是接受上游跳过缺失靶点的行为。

## 与官方 Knk 的兼容性 / Compatibility

`virt_knock::tensor_decompose_mean` 实际是 bootstrap 网络逐元素平均和缩放，没有执行 CANDECOMP/PARAFAC；不能把它自动当作官方 `tensorDecomposition` 的等价实现。其 manifold 用 float32 全矩阵 `torch.linalg.eigh`，官方 R 实现也不保证产生同一个数值嵌入。第三方原生 API 演示因此单独保存到 `results/native_api/`，不与官方 WT 的主 KO/OE 比较混合。

主比较直接复用模块 01 官方构建的 WT，调用真实第三方 OE 边修改函数，再对 KO 与全部 OE 使用相同 `scTenifoldNet::manifoldAlignment(d=2,nCores=1)` 和 `scTenifoldKnk::dRegulation(gKO="SPI1",empiricalNull=FALSE,direction=NULL)`，随机种子为 1。使用官方函数并不使 OE 生物学模型获得官方或实验验证。

已复现 **pandas 3 兼容性错误**：同质 dtype DataFrame 的 `.values` 是只读，源码 `np.fill_diagonal(oe.values,0)` 报 `ValueError: underlying array is read-only`。混合 dtype 的 CSV 列可能偶然产生可写副本，所以一次真实数据成功不能证明 pandas 3 可靠。模块固定兼容的 pandas 2.2.3，wrapper 对 pandas >=3 明确拒绝；没有修改第三方源码来隐藏错误。

上游 `examples/batch_oe.py` 从 `virt_overexpression` 顶层导入 `run_qc`、`build_pc_networks` 等未被顶层导出的名字，不能直接按该脚本宣称成功。实际验证的是 README 中公开的高层 `scTenifoldOECUDA` API；上游没有附参考数据或预期数值，本次用真实 PBMC counts 补齐输入，参数明确记录为小规模执行验证。

## 科学解释边界 / Interpretation boundary

边权乘 2、3 并不代表靶基因 mRNA 增加 2、3 倍。正/负回归边权也不是已验证的激活/抑制因果关系。输出 distance、FC、Z 是网络差异统计量，不是表达 log2FC。未调用 KO 专用方向预测来推断 OE 上调/下调。拟合的网络来自一个公开供体和有限 600 基因面板，结果为功能假设；需要独立 perturb-seq、CRISPR 或过表达实验验证。

方法背景：[scTenifoldKnk 论文](https://doi.org/10.1016/j.patter.2022.100434)、[scTenifoldNet 论文](https://doi.org/10.1016/j.patter.2020.100139)、[官方 scTenifoldNet](https://github.com/cailab-tamu/scTenifoldNet)。这些论文不构成对本第三方 OE 操作的独立验证。
