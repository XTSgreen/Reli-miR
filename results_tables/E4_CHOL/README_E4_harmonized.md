**同RSEM量化流程的外部固定模型评价。**

来源核验确认主体mRNA为log2(RSEM norm_count+1)，原外部mRNA数值支持log2(STAR count+1)，量化与归一化流程不同；旧结果保留在results/E4，作为输入不一致的诊断记录。新的外部输入来自官方Xena CHOL HiSeqV2，与主体采用同RSEM流程，原35例均唯一匹配主肿瘤样本，3,000个模型基因全部覆盖，原miRNA与空clinical对象保持不变。

固定原PCA_ElasticNet模型、原训练预处理和内部Rel，在35例外部患者、742项共有miRNA上评价，未按外部真值调模型、阈值或评分权重。总体平均R²为−0.30603（患者bootstrap百分位区间−0.70096至−0.27922），中位Pearson为0.28922（患者bootstrap百分位区间0.21071至0.35313），平均训练标准化RMSE为0.85171（患者bootstrap百分位区间0.80361至0.88256）。每项区间均有2,000次有效患者重采样；原点估计与bootstrap分布中位数分别保留。

内部Rel≥0.7的原33项名单在外部共有32项，其中中位Pearson为0.57622、中位R²为0.04947，16/32项同时达到R²>0且Pearson≥0.5。这些是固定内部名单在外部的描述性点估计，不构成根据外部结果重选的名单。

内部Rel与外部RMSE的rho为−0.23418（患者bootstrap百分位区间−0.30999至−0.12275），与外部R²的rho为0.06038（患者bootstrap百分位区间−0.03771至0.13904），与外部Pearson的rho为0.36636（患者bootstrap百分位区间0.22437至0.38149）。总体平均R²仍为负，迁移结论应限定在部分条目的表达预测及排序，不能推广为整体可靠替代实测。外部临床/生存TSV已经找到并核验单位，本轮未执行外部生存或临床效用分析。

输入基因缺失为0，缺失基因屏蔽诊断记为不适用；原缺失382基因的诊断属于旧输入。等规模简单评分及200个随机名单结果保留，随机名单间变异与患者区间分别解释。证据：results/harmonized/E4中的逐miRNA表、固定阈值表、完整区间及bootstrap复现表；来源：results/E4_harmonized_inputs中的JSON及病例映射。

**复现命令。** 在工作区根目录执行：

```powershell
& 'D:/R-4.6.1/bin/Rscript.exe' new/code/score_external_analysis.R --section E4 --bootstrap 2000 --external-input new/results/E4_harmonized_inputs/external_CHOL_matched_RSEM.rda --output-root new/results/harmonized
& 'D:/R-4.6.1/bin/Rscript.exe' new/code/uncertainty_completion.R --section external --bootstrap 2000 --external-input new/results/E4_harmonized_inputs/external_CHOL_matched_RSEM.rda --output-root new/results/harmonized
```

E4_historical_prediction_reproduction_check.csv在本输入下表示与历史旧输入预测指标的差异，不要求为0；它不参与新结果估计。
