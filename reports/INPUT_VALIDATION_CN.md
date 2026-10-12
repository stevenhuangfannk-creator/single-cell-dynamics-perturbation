# 输入工具验证 / Input validation

2026-10-12 在 Windows、Python 3.13.15 独立项目环境中执行：

```text
python -m unittest tests.test_input -v
Ran 6 tests ... OK
```

覆盖非法/非整数 counts、重复基因、空细胞、缺失靶点、MEX 基因×细胞转换、metadata ID 对齐、明确 counts layer 的 h5ad 读取和 PI 审查输出。测试使用小型人工矩阵，目的仅为验证软件拒绝非法输入和保持轴/身份，不是生物学模拟结果。环境主依赖精确版本见 environments/shared-requirements.txt，完整本机包清单见 environments/shared.lock.txt。

新增AnnData→MEX导出往返测试，确认counts值、基因×细胞轴方向、唯一基因ID及细胞顺序均保留。Seurat 导出脚本提供版本分支但没有读取真实 Seurat 对象，后续仍须以实际对象核验。PI 审查不会评估真实谱系证据或自动生成疾病结论。模型运行与科学有效性由各模块单独检查。
