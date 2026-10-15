**E4多队列扩展预登记（冻结于任何迁移结果计算之前）**

日期：2026年10月15日。本文件在运行任何一个新癌种的评价之前写定，分析启动后不再修改参数。性质：计划书第十节列明的可选扩展（更大同流程配对队列），用于检验内部质量分数排序在外部的方向一致性；不预设正向结果，全部结果无论方向均进入阶段结果文档。

**一、队列选择规则（先于数据检查冻结）**

TCGA官方GDC/Xena同RSEM流程，选择7个样本量最大的非肝实体癌种：KIRC、LUAD、COAD、STAD、PRAD、THCA、HNSC。选择标准仅依赖样本量和数据可得性，不参考任何模型表现。mRNA来源为UCSC Xena HiSeqV2（log2(RSEM norm_count+1)，2017-10-13版，与主体LIHC及外部CHOL完全同一Wrangling流程）；miRNA来源为GDC官方mirbase21 miRNA定量文件（前体ID空间，与内部846项ID逐字一致；取reads_per_million_miRNA_mapped，按log2(RPM+1)变换，与E0已核实的内部量纲一致）。

**二、样本与特征处理规则（冻结）**

1. 每病例仅取一个原发肿瘤样本：miRNA按case_id去重取sample_id最小文件，重复文件记录为dropped_duplicate_case_file；mRNA取-X-01列，同病例多列时取字母序首个并记录。
2. 分析患者为mRNA与miRNA均完整配对的原发肿瘤病例。
3. 不做任何按外部表现的特征筛选、填补方式调整或病例剔除；模型3000个输入基因直接与Xena矩阵求交，缺失基因按原脚本训练中位数填补并输出覆盖诊断。
4. miRNA特征为GDC文件全部条目，与内部846项按名称求交作为共有集，不新增或重命名任何条目。

**三、评价协议（完全继承E4同流程CHOL运行，冻结）**

1. 固定模型：E0确认的原投稿PCA_ElasticNet模型对象及训练预处理（中心、尺度、填补、PCA、lambda.min）；不在外部重拟合任何参数。
2. 固定评分：内部评分集Rel（等权公式），阈值0.5/0.6/0.7及等规模33/74/177比较沿用预登记设置；随机对照200个列表、seed 20261018。
3. 患者bootstrap：B=2000，seed 20261017，百分位区间；操作性可重建标签维持测试R²>0且Pearson≥0.5。
4. 主报告量：各队列总体mean R²、median Pearson、mean RMSE及区间；内部Rel与外部逐miRNA RMSE/R²/Pearson的Spearman及区间；原33项及E0名单在外部共有数与达标比例；与CHOL结果的方向一致性汇总。
5. 判定规则（继承计划书）：若多数队列排序迁移方向与内部一致，报告为排序迁移的一致性证据；绝对R²为负的队列仍如实报告，不据外部结果重选模型、权重、阈值或miRNA；各队列结果分别呈现，不做跨队列p值合并。

**四、GEO独立队列定位（冻结）**

GSE76903（miRNA，成熟链ID+salmon counts）与GSE77509（mRNA，Ensembl基因ID）为20例配对HBV相关HCC。成熟链ID与前体ID空间非一一对应（如hsa-let-7a-5p横跨hsa-let-7a-1/2/3三个基因座），且量化流程与主体不同；按计划书「不兼容数据仅作方法学背景比较，不拼接计数」的纪律，该队列不进入正式E4评价，仅保留数据核验记录，作为讨论中独立来源的方向性背景引用。

**五、执行与产物位置**

输入构建：new/code/build_multicohort_inputs.R；清单获取：new/code/fetch_gdc_mirna_manifests.ps1；批量下载：new/code/download_gdc_mirna.ps1。输入RDS与验证表存new/results/E4_multicohort/inputs/<COHORT>/；评价结果存new/results/harmonized_multicohort/<COHORT>/E4/；汇总表生成后写入new/results/E4_multicohort/。执行状态随跑随记入EXPERIMENT_TRACKER。
