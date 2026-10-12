# 中间对象 / Objects

原始数据、经过 QC 的 counts、训练 counts、推断 AnnData 和官方 `.pth` 模型保存在本机 `work/sctour/`；它们被 Git 排除，通过 README 中的下载、预处理和训练命令重建。

`A_qc_counts.h5ad`/`B_qc_counts.h5ad` 保留全部经过检出筛选的基因及 counts layer；`*_train_counts.h5ad` 只保留对应 HVGs。每个正式模型目录保存 `inferred.h5ad`，其中包含 `ptime_before_manual`/`ptime`、`X_TNODE`、`X_z`、`X_zt`、`X_VF`、`X_umap_pca` 和 `X_umap_latent`。单图重绘读取此对象。模型 `save_model` 保存训练参数、原训练对象、状态字典和方向信息；检查点加载使用实际 `sct.predict.load_model` API。
