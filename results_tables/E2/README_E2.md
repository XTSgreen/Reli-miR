# E2 uncertainty completion

本目录中的 `E2_complete_*` 是对先前正式 E2 结果的补充，不重跑主分析中的 FNR、precision、recall 或同规模随机集合分析。计算脚本为 `new/code/uncertainty_completion.R`，使用 E0 canonical 输入与 E0 固定预测矩阵。

`E2_complete_model_bootstrap_raw.csv` 和对应 summary 对 E0 三个主模型进行 2,000 次配对患者 bootstrap。每次在相同测试患者重采样上计算三个模型的 mean R²、median Pearson 和 mean RMSE；`E2_complete_model_paired_contrasts.csv` 给出模型间配对差及百分位区间。区间条件于已拟合模型及固定 E0 划分，不表示重新训练模型后的不确定性。

`E2_complete_score_joint_bootstrap_raw.csv` 与 summary 对 train-native、calibration、test 患者分别重采样，并在每次重复中重估各评分。表中包括 Rel、R²、Pearson、权重组合、calibration RMSE rank、training-native mean 与 variance。`point_rho_baseline_minus_Rel` 定义为基准评分与 Rel 的 Spearman 相关之差；正值表示 Rel 与实际测试 RMSE 的负相关更强。每个评分及配对差异均由同一重复内的患者重采样计算。抽样分布的 median 与原始样本 point estimate 是不同统计量；部分区间不含原始 point estimate，报告时应同时保留两者，不得将 bootstrap median 替代 point estimate。有效重复数为 2,000。

`E2_complete_sample_bootstrap_*` 保存样本级分数与实际 PCA_ElasticNet 测试 RMSE 关系的患者级配对 bootstrap；字段 `rho_vs_actual_RMSE` 使用真实 RMSE，相关方向为负值。旧的 CI 完全等于 point estimate 的错误表另存为 `E2_sample_level_associations_invalid_initial.csv`，修正后的表为 `E2_sample_level_associations.csv`。

E1 OOF 误差基线尚未纳入：当前没有与 E0 测试患者及训练来源严格匹配的 OOF 分数，因此不跨 split 合并，也不据此形成比较结论。此前 B=3 smoke 文件已被正式 B=2,000 文件覆盖；运行前的快速指标一致性检查通过，评分 point estimate 与原 E2 表逐项一致。
