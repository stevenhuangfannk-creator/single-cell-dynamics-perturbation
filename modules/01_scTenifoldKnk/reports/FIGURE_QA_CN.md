# 图件与源数据审查 / Figure QA

目标：展示真实网络扰动计算、统计结果和稳定性不足；不把无方向distance转成表达上/下调，不将网络边作因果图。图形属于定量grid与network/quant混合。后端Python matplotlib，白底，Arial/DejaVu Sans，PDF TrueType、SVG文本可编辑。独立图按可读性使用约140–274mm宽，未声称符合某特定期刊最终单栏尺寸；可在源脚本修改尺寸后重导。

全部11图导出PNG600dpi、PDF、SVG。源码preflight 0fail；所有导出PDF字体审查最小7pt，满足5pt底线，详见 [pdf_text_audit.json](pdf_text_audit.json)。本次没有用户要求之外的TIFF；静态验证器不能从`args.dpi`推断600，但实际默认参数与PNG元信息均记录600。log使用max(1e-300,p) / log1p保护；图08显式记录按行标准化无方向heat distance。

| 图 | 独立问题 / 原始结果 | 视觉和证据检查 |
|---|---|---|
|01|WT子网络；`objects/WT.csv`|32nodes、top80内部边并补SPI1 top15出边；环形布局防label碰撞；完整网络另存，不把显示子图误作全部网络|
|02|KO子网络；`objects/KO.csv`|与WT完全相同坐标，去SPI1出边；其他边一致；布局仅用于显示，不改权重|
|03|变化是否严格符合KO定义；WT/KO矩阵|SPI1非零出边567→0；其他行逐元素相等；右侧15边排序明确|
|04|top20非targetdistance；`single_SPI1.csv`|target排除，另算direction颜色的含义明确；distance不能标log2 expressionFC|
|05|距离与显著性；同DR表|排除构建必变target；BH0.05虚线；4个最高距离label直接指向对应点|
|06|SPI1局部关系；WT矩阵|top16绝对出边、环绕中心；PCR符号未声称真实激活/抑制|
|07|功能富集；`enrichment_all_terms.csv`|GO/KEGG/Reactome各显示top6含可能不显著term；灰色为BH≥0.05；断行确保长term可读；虚线BH0.05|
|08|9候选响应；`candidate_screen_distance.csv`|同一heat方法；top30按平均distance，包含target自身的大distance且明确为display；无符号色条|
|09|single/joint响应；两DR表|使用同WT、同MA d2，排除双方target；右图paired距离，不称表达倍数；最高点label分离|
|10|全部稳定性；primary与1200敏感性表|6组primary加3组更宽panel全部展示；参考自比较=1不当作独立重复；1200对角线1为自比较，独立seed约0.19–0.22|
|11|官方mapAPI内部一致性；`candidate_perturbation_map.csv`|9候选，SPI1同模型signed-heat-profile内部reference；非疾病签名；FOS/JUN/CEBPD label偏移防叠，未改坐标|

图01/02与06的节点坐标为可重复的显示布局，没有生物空间/拟时序含义。图08只选显示top30，不删除完整600gene结果。图11只显示已经真实计算的9candidate，没有填补未运行转录组目标。

统计单位：1个公开健康供者的480CD14cells；随机seed/细胞子抽样是技术运行而非生物学重复。图10画每次真实运行值，没有声称mean±CI，因此不加虚构误差线。图07为单次模型hit集合的超几何ORA及BH；相关/重叠gene sets不是独立机制重复。

逐图PNG与PDF字体全部检查。最终人工visualQA修正了初版spring network布局造成的标签拥挤、图09重叠点label和图11密集label；这些修改仅影响显示位置，不变动计算坐标/边/DR。图10追加1200敏感性后重新输出三格式，manifest文字同步。没有用图片后期工具改变统计结果。
