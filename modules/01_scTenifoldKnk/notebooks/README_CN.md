# 交互使用 / Interactive use

本模块将计算与绘图保留为可直接运行、便于版本追踪的 R/Python scripts，没有尚未执行的占位 notebook。交互探索时可在 R Notebook 中读取 `objects/single_SPI1.rds`，或在 Python notebook 中读取 `results/*.csv`；正式单图输出仍调用 `scripts/plot_figures.py --figure NN`，避免依赖隐藏 notebook 状态。运行顺序与输入格式见上级 README。
