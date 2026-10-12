# 网络 CSV 身份回归 / Network CSV identity regression

核查日期：2026-10-12。状态：`RUN_OK`，范围为共享 WT/KO 主入口的 CSV 读入、导出与响应计数。没有重新拟合网络，没有执行新的扰动、MA 或 DR 分析，也没有运行用户数据。

## 修复前证据 / Before fix

原 Python `pd.read_csv(..., index_col=0)` 会把全数字字符串基因 ID `001/002` 转为整数 `1/2`，列名却保留 `001/002`，因而错误拒绝合法方阵；字面 ID `NA` 被转成缺失值。原 R `read.csv(..., row.names=1)` 同样丢失前导零，并对 `NA` 行名报缺失错误。四个实际探针已保存在 [network_csv_before_fix.json](network_csv_before_fix.json)。

后续 DR 表的默认解析也会把 `gene='NA'` 变成缺失，使靶点错误计入非靶显著响应。小型 CLI 回归在修复前实际得到 `[2,2]` 而非预期的 `[1,1]`，失败日志见 [network_csv_DR_before_fix.log](network_csv_DR_before_fix.log)。该测试用临时表替代 R 输出，只验证输入身份和计数逻辑，不代表计算了新的生物学结果。

## 修复与实际验证 / Fix and checks

Python 显式对首列使用字符串 converter 后设置 index，保留权重的原数值解析；R 首列指定 `character`，其余列指定 `numeric`，再设置矩阵行名；DR 表只对 `gene` 列指定字符串 converter。扰动和官方 MA/DR 调用保持原算法。

| 检查 | 结果与证据 |
|---|---|
| Python 3.11 / pandas 2.2.3 | 3 项必要测试全部通过：[日志](network_csv_python_regression.log) |
| shared Python 3.13 / pandas 3.0.6 | 同 3 项轻量测试全部通过；不导入第三方 OE/torch：[日志](network_csv_shared_regression.log) |
| 实际 R 4.3.1 读入 | `001/002`、`001/002/NA` 的行列身份一致，正/负权重精确相等：[日志](network_csv_R_regression.log) |
| Python 实际导出 → R 实际读入 | 三个 ID 和全部小型矩阵权重精确相等：[日志](network_csv_cross_language.log) |
| 已保存公开 WT/KO | 600×600 两个矩阵的顺序 ID 和数值逐元素完全不变：[检查](network_csv_existing_matrix_check.json) |

这些测试使用可精确表示的正/负小数和零，检查数据类型与逐元素相等；没有将 ID 汉化或重新编号。R 日志保留本机 locale 启动警告，读取断言实际通过。

从仓库根目录可重新执行必要检查：

```powershell
& $Py -m unittest tests.test_oe_csv -v
Rscript tests/check_oe_network_csv.R
```

`$Py` 使用模块教程的已安装 Python 环境。R 检查只提取生产脚本的 CSV reader，不载入分析包或调用 MA/DR。Python 的 CLI 计数测试隔离第三方计算，只检验输入/结果表读入和靶点排除；轻量 GitHub CI 已加入这 3 项测试。

## 新数据边界 / Reuse boundary

这一修复验证共享 WT/KO 网络 CSV 的身份保持。第三方原生 counts API 演示、10 图绘制和完整富集仍按本次 PBMC/SPI1 的明确范围使用；替换数据时仍需审查网络信息、版本、目标、来源注释和科学稳定性。CSV 接口通过不构成任意新输入、跨供体或真实扰动实验的生物学验证。
