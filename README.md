# 单细胞动力学与虚拟基因扰动方法库 / Single-cell dynamics and virtual perturbation

本仓库以真实公开数据和官方软件为基础，建设 scTenifoldKnk、Experimental OE 和 scTour 三个独立方法模块。当前正在执行环境核查、数据下载和计算复现；运行状态以各模块报告为准，不能把目录建立或 Smoke Test 当作论文复现。

An evidence-grounded collection of reproducible single-cell regulatory perturbation and dynamics workflows. Experimental outputs are computational hypotheses, not validated perturbation or lineage experiments.

## 方法入口 / Modules

| 模块 | 方法来源 | 入口 |
|---|---|---|
| 虚拟 KO | [官方 scTenifoldKnk](https://github.com/cailab-tamu/scTenifoldKnk) | [01_scTenifoldKnk](modules/01_scTenifoldKnk/README_CN.md) |
| 虚拟 OE（Experimental） | [第三方 virt_overexpression](https://github.com/virchow-zhang/virt_overexpression) | [02_scTenifoldOE](modules/02_scTenifoldOE/README_CN.md) |
| 动力学推断 | [官方 scTour](https://github.com/LiQian-XC/sctour) | [03_scTour](modules/03_scTour/README_CN.md) |

## 科学边界 / Interpretation

scTenifoldKnk 的网络扰动距离、FC 与统计显著性不等于真实表达差异。OE 网络增强倍数不等于 mRNA 过表达倍数。scTour 可使用表达矩阵而不要求 spliced/unspliced，但静态表达中的动力学预测不能直接证明真实谱系转化。PI 与 Healthy 的分组不是时间顺序。

大型原始数据、AnnData、模型及中间对象不上传 Git；通过数据登记、固定配置和下载/运行脚本重建。小型结果表、真实图件及日志可纳入版本管理。
