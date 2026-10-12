# 实际运行与资源 / Executed resource records

硬件与软件版本见 [hardware.json](hardware.json) 和各模块环境文件。下表仅列已完成的真实运行；未记录的峰值不补造。

| scTour运行 | Epochs | 训练秒数 | 全流程秒数 | Windows进程 peak working set / MiB | Torch CUDA peak allocated / MiB |
|---|---:|---:|---:|---:|---:|
| A20_seed0 | 400 | 518.47 | 598.33 | 1263.27 | 46.59 |
| A90_seed0 | 400 | 407.98 | 未记录 | 未记录 | 未记录 |
| A90_seed1 | 400 | 359.71 | 未记录 | 未记录 | 未记录 |
| B_1 | 400 | 349.27 | 512.70 | 1028.59 | 156.96 |
| B_10 | 400 | 690.24 | 728.39 | 1259.65 | 164.48 |
| B_20 | 400 | 2032.95 | 未记录 | 未记录 | 未记录 |
| B_40 | 400 | 1731.80 | 1768.56 | 1468.22 | 164.48 |
| B_5 | 400 | 655.47 | 818.08 | 1348.44 | 157.25 |
| B_60 | 400 | 2147.98 | 2331.31 | 1512.41 | 164.48 |
| B_80 | 400 | 3420.00 | 3606.85 | 1069.63 | 164.48 |
| B_95 | 400 | 1518.83 | 未记录 | 未记录 | 未记录 |
| smoke_A（Smoke，仅API检查） | 3 | 14.27 | 未记录 | 未记录 | 未记录 |

训练时间按runner计时；正式模型使用CUDA编码/解码和官方CPU ODE。全流程计时包含训练后的推断及嵌入，不含原始数据下载/预处理或独立绘图。加载对象验证不计作重新训练。Windows进程peak working set与Torch CUDA allocated是不同口径；早期没有峰值字段的运行保留“未记录”。

KO主流程34.16秒；1200基因三seed扩展55.61秒。主KO没有独立进程峰值RAM记录；通用80×150 tiny入口单/联合KO各0.45秒，单KO Rterm以50ms采样的RSS峰约258MiB，详见 [tiny验证](../modules/01_scTenifoldKnk/reports/custom_smoke_summary.json)。tiny仅验证接口。

OE主比较含导入54.56秒；Python与当时R子进程RSS抽样加总峰750.6MiB（0.1秒采样），共享页可重复计数、短峰可遗漏，详见 [OE资源证据](../modules/02_scTenifoldOE/reports/resource_usage.json)。

系统资源共582次抽样，范围2026-10-12T18:08:30+08:00至2026-10-12T23:26:06+08:00；最低可用RAM 5.0MiB。整卡显存抽样最大7125.0MiB，**包含其他应用，不能称为模型独占显存峰值**。

原始样本见 [resource_samples.csv](resource_samples.csv)，机器可读摘要见 [resource_summary.json](resource_summary.json)。系统压力与并行应用影响耗时，本次数据不构成独占硬件benchmark。
