# 公共数据来源与单位核验

检索日期：2026-10-14。查询范围限于 GDC 官方 API、GDC 官方文档、UCSC Xena 官方目录和工作空间已有 TCGA 临床/生存 TSV；表达校验包括5例GDC官方文件及单一LIHC Xena矩阵（22,166,005字节），未上传本地表达矩阵、模型预测或结局数据。

## 核验结论

1. 主体队列 `data/data.rda` 中的 **366 例 OS_time 单位可以确认为天**。GDC 官方病例 API 的病例状态与本地 `OS_status` 全部一致；239 例存活者的 `OS_time` 逐例等于 GDC 记录的最后随访日数，127 例死亡者逐例等于 GDC `days_to_death`，时间值 **366/366 精确相等**。这将 OS_time 的来源联系到 GDC 以诊断为起点记录的 day-based 字段。

2. 外部 `survival.tsv` 的 `OS.time` 同样可确认为天。文件共有57行、48名不同病例，48名病例的时间和死亡状态均与 GDC 官方记录精确相等；外部 `data.rda` 的35名病例均包含在这48名病例内，逐例时间与事件状态全部相同。此项仅确认已有 TSV 的变量含义，不构成外部生存验证分析。

3. **miRNA 数据支持 log2(RPM+1) 处理。** 外部 CHOL 的 RDA 矩阵与同目录 `mirna.tsv` 在35例、775项 miRNA 上共27,125个数值逐项相等。三个 CHOL 病例的 RDA 数值与 GDC 病例级官方定量文件的 `reads_per_million_miRNA_mapped` 经 `log2(RPM+1)` 变换后，三例各775项均精确相等。主体 LIHC 两个病例的846项 miRNA 也都在1e-6容差内与该变换相等。该结果强力支持两组 miRNA 输入使用 `log2(RPM+1)`；由于缺少生成 RDA 的原始脚本和全队列逐样本文件核对，不能据此恢复主体全部366例的原始下载版本和处理日志。

4. **外部 CHOL mRNA 数据支持 log2(count+1) 处理。** 三个病例的 GDC 官方 STAR 文件中，36,419个基因符号可与 RDA 对应。对每例可匹配的全部表达单元，GDC `unstranded` count 经 `log2(count+1)` 变换后，与 RDA 的非零表达单元精确比例分别为99.84%、99.84%和99.84%；Pearson相关为0.99954至0.99959，RMSE为0.123至0.128。RDA 与本地 `mRNA.tsv` 经 `anno.probemap` 的 Ensembl ID 与基因符号对应后，35例×36,419项共1,274,665个值中，1,272,543个在1e-8容差内相同；两边均非零的1,034,675个值中有1,032,918个相同，精确比例99.83%。这为外部 mRNA 的变换方式提供了较强数值证据。

5. **主体 LIHC mRNA 与官方 UCSC Xena HiSeqV2对象全值一致。** 按366个病例的12位条码逐一选择唯一`-01`主肿瘤样本列，366/366病例均匹配，无缺失或重复。主体矩阵与该Xena对象共有的18,204个基因，在6,662,664个病例×基因单元上全部逐值相同（绝对差和RMSE均为0）。Xena相邻官方JSON将单位定义为 `log2(norm_count+1)`，版本为2017-10-13，并说明Level 3 `*.rsem.genes.normalized_results`经`log2(x+1)`转换后进入Xena。因此主体全部366例输入数值可依据全矩阵直接匹配和该JSON元数据，标记为 `log2(norm_count+1)`。这核验数值来源对应关系与Xena版本定义，不补造项目最初下载日期。两例当前GDC STAR的unstranded、TPM、FPKM、FPKM-UQ字段在原值、`log2(x+1)`和`ln(x+1)`共12组变换中均未出现非零逐值精确匹配。

## 数据与时间字段的核验方式

GDC 官方说明将 `days_to_death`、`days_to_last_follow_up` 和 `days_to_follow_up` 用作死亡或末次随访时间，并说明生存状态由 `vital_status` 给出。官方提交指南将 day-based 时间点定义为从初始诊断开始计日。来源：[GDC 临床生存分析说明](https://docs.gdc.cancer.gov/Data_Portal/Users_Guide/clinical_data_analysis/)、[GDC 日期字段说明](https://docs.gdc.cancer.gov/Data_Submission_Portal/Users_Guide/Best_Practices/)。

官方 mRNA 工作流将 STAR 基因级计数及 FPKM、FPKM-UQ、TPM 提供为不同字段；当前 GDC 下载文件的头部标注 GENCODE v36。官方 miRNA 工作流文件分别提供 `read_count` 和 `reads_per_million_miRNA_mapped`，使用 miRBase v21 注释。来源：[GDC mRNA 工作流](https://docs.gdc.cancer.gov/Data/Bioinformatics_Pipelines/Expression_mRNA_Pipeline/)、[GDC miRNA 工作流](https://docs.gdc.cancer.gov/Data/Bioinformatics_Pipelines/miRNA_Pipeline/)、[GDC Data Release 32说明](https://docs.gdc.cancer.gov/Data/Release_Notes/Data_Release_Notes/)。

核对主体病例时，按 `clinical.patient_id` 对 GDC `/cases` 接口进行病例级匹配，366/366病例项目为 TCGA-LIHC。对每例计算两类候选时间：若 `vital_status=Dead`，使用 `demographic.days_to_death`；若为 `Alive`，使用可用诊断及随访记录中的最大 `days_to_last_follow_up` 或 `days_to_follow_up`。所得时间与本地 `OS_time` 全部精确相等，状态也全部一致。

外部 `clinical.tsv` 含60行、51名 TCGA-CHOL 病例；GDC 官方 `/cases` 项目查询返回51名病例，状态和可用死亡日数字段均与本地表对应。`survival.tsv` 中57行归并为48名病例；48/48例的 `OS.time` 与相应 GDC 死亡或末次随访日数精确相等，事件状态也48/48相同。外部 RDA 的35名表达病例均在这48名中。

## 表达数值比对

### miRNA

外部 RDA 与本地 `mirna.tsv` 的精确匹配使用患者条码前12位作为病例键，按每位病例的01A肿瘤样本匹配。35×775共27,125个单元，全部完全相同。对 GDC 官方文件，三个 CHOL 病例 `TCGA-W5-AA2X`、`TCGA-W5-AA30`、`TCGA-ZH-A8Y2` 的775项数值均满足 `RDA = log2(RPM+1)`，最大绝对差为0；两个 LIHC 病例 `TCGA-XR-A8TC`、`TCGA-BC-A112` 的846项结果最大绝对差分别小于5.7×10⁻⁷，均低于1e-6核验阈值。

### mRNA

外部 CHOL RDA 对本地 TSV 的数值比对按 `anno.probemap` 将 RDA 的基因符号映射到 Ensembl ID。整体相关为0.99957，RMSE为0.126；两矩阵同为零的239,625个单元及双方均非零的1,032,918个匹配单元均相等。尚有1,757个双方非零但不一致的单元、365个仅一方为零的单元及少数映射/注释差异。对 GDC 当前 STAR 文件的三个病例，双方均非零时精确匹配率分别为99.84%、99.84%、99.84%，但仍有189个跨三例累计的非精确单元。当前 GENCODE v36 的 gene symbol 与本地 `anno.probemap` 的逐项对应没有将这些差异全部解释，因此不能把这些少数差异归因于重复符号映射。

主体 LIHC 对 GDC 当前 STAR 文件的两例比对，在共有15,500个基因符号中有14,243和14,377个双方均非零表达，然而这两个非零集合中没有数值达到1e-6精确匹配；当前 GDC STAR `log2(count+1)` 与矩阵的相关为0.948和0.939。当前 GDC `log2(TPM+1)`、`log2(FPKM+1)`、`log2(FPKM-UQ+1)` 对应相关分别约为0.878、0.833、0.851，亦不能确认主体矩阵由这些字段产生。

## 原始数据来源仍缺失的部分

本地源码只读取预先保存的两个RDA文件，没有生成或下载步骤；RDA也未记录自定义来源属性。主体LIHC矩阵与官方Xena `TCGA.LIHC.sampleMap/HiSeqV2` 对象在全部366例×18,204基因上逐值完全一致，且该对象JSON定义了 `log2(norm_count+1)` 和2017-10-13版本。由此可核实主体矩阵数值尺度及其与该Xena版本的全值对应关系，但不能恢复本项目首次下载时间或单独取得每个本地样本的原始TCGA文件UUID。当前GDC STAR文件不匹配，不用于取代已核验的Xena对象来源。

写回复和修稿时建议明确：主体 `OS_time` 单位为天，时间按诊断起点；两组miRNA的数值证据支持 `log2(RPM+1)`；外部CHOL mRNA的数值证据支持 `log2(count+1)`；主体LIHC 366例×18,204基因与官方Xena HiSeqV2对象逐值相同，该对象定义单位为 `log2(norm_count+1)`、版本日期为2017-10-13。模型拟合前的训练集均值/方差标准化应另行说明。

## 复核产物

1. GDC 病例与表达文件记录：[主体病例临床 API 响应](public_primary_gdc_clinical.json)、[外部 CHOL 病例临床 API 响应](public_external_gdc_clinical.json)、[10个病例级文件的 UUID、版本、MD5和下载链接](public_gdc_sample_file_manifest.csv)。10/10个本地下载文件的 MD5 与 GDC 元数据一致。

2. 时间匹配：[主体 366 例 OS 精确对照](public_primary_OS_match.csv)、[外部 `survival.tsv` 48例对照](public_external_OS_match.csv)、[外部临床 TSV 病例级对照](public_external_clinical_tsv_match.csv)。

3. 表达匹配：[外部 RDA 与本地 TSV 对照](public_external_rda_vs_local_tsv.csv)、[五病例 GDC 多组学字段转换对照](public_gdc_multisample_expression_match.csv)、[mRNA 非零/零值精确匹配细表](public_gdc_mrna_match_cell_audit.csv)、[LIHC mRNA候选转换比较](public_LIHC_mRNA_transform_checks.csv)、[log2整数反变换特征核查](public_expression_inverse_log2_signature.csv)。病例级 GDC 文件保存在 `new/results/E0/public_gdc_sample_data/`。

4. Xena候选源材料：[HiSeqV2数据集JSON](public_TCGA_LIHC_HiSeqV2_dataset_metadata.json)、[来源与下载元数据](public_TCGA_LIHC_HiSeqV2_metadata.json)、[端点探测结果](public_xena_adjacent_metadata_probe.json)、[TCGA DCC DESCRIPTION链接返回页](public_TCGA_LIHC_RNASeqV2_DESCRIPTION.txt)。

## 主体 mRNA 的 UCSC Xena 定向核对

为检验主体矩阵是否与官方 UCSC Xena 的 TCGA-LIHC `TCGA.LIHC.sampleMap/HiSeqV2` 数据集一致，本次仅获取该单一候选文件。UCSC Xena 官方 `defaultDataset.json` 将该数据集列为 TCGA LIHC 的 gene expression 默认数据集，数据文件下载链接返回 HTTP 200，压缩文件大小为22,166,005字节，HTTP Last-Modified 为2021-04-08，S3版本ID及本地SHA-256保存在元数据记录中。初次尝试Xena数据API元数据请求返回 HTTP 500/NPE；后续检查官方下载对象相邻的 `.json` 端点成功，取得了数据集单位、变换方法、版本和来源描述。其他相邻 `.meta`、`.metadata`、`.meta.json`、`.txt` 端点返回403。

按366个主体病例的12位TCGA条码逐一匹配Xena列名中的唯一`-01`样本，366/366均找到唯一主肿瘤样本，无缺失、无重复。对每例18,204个共有基因逐值比较，6,662,664/6,662,664个RDA与Xena单元完全相等（阈值1e-6）；最大绝对差为0、RMSE为0，每例均为18,204/18,204相等。Xena对象相邻官方JSON将该数据集的`unit`标为 `log2(norm_count+1)`，版本为2017-10-13，并记载Level 3 `*.rsem.genes.normalized_results`经 `log2(x+1)`处理后进入Xena。结合全队列逐值匹配，**主体366例输入矩阵可核实为该Xena对象所定义的 `log2(norm_count+1)` 数值**。这确认本地数值与该Xena版本一致，但不能恢复本项目最初下载日期。

另对两个病例的当前GDC STAR `unstranded_counts`、TPM、FPKM及FPKM-UQ字段分别比较了原值、`log2(x+1)`、`ln(x+1)`三类转换，共12组比较。四个字段的原值及自然对数/以2为底对数转换均未形成对主体mRNA的精确非零值匹配；具体相关和RMSE见[转换比较表](public_LIHC_mRNA_transform_checks.csv)。这与Xena候选的逐值相同结果一致，故主体矩阵更符合所核对的Xena候选文件内容，而不是本次核对的当前GDC病例文件字段。

复核文件包括：[Xena来源与HTTP文件元数据](public_TCGA_LIHC_HiSeqV2_metadata.json)、[官方数据集元数据](public_TCGA_LIHC_HiSeqV2_dataset_metadata.json)、[366例病例映射](public_LIHC_Xena_HiSeqV2_case_map_366.csv)、[366例逐值比较统计](public_LIHC_Xena_HiSeqV2_match_summary_366.csv)、[全部病例差异明细](public_LIHC_Xena_HiSeqV2_differences_366.csv)、[候选压缩矩阵](public_TCGA_LIHC_HiSeqV2.gz)。复算脚本为[主体两病例导出](../../code/public_export_two_lihc.R)和[Xena 366例全矩阵核对](../../code/public_compare_xena_hiseqv2_366.R)。官方目录来源：[UCSC Xena cohort registry](https://raw.githubusercontent.com/ucscXena/cohortMetaData/master/defaultDataset.json)；官方文件地址：[TCGA LIHC HiSeqV2](https://tcga.xenahubs.net/download/TCGA.LIHC.sampleMap/HiSeqV2.gz)。



### UCSC Xena 语义元数据回退记录

官方下载对象相邻端点 `https://tcga.xenahubs.net/download/TCGA.LIHC.sampleMap/HiSeqV2.json` 返回HTTP 200和2,583字节JSON，包含 `unit=log2(norm_count+1)`、`version=2017-10-13`、`wrangler=Xena TCGAscript RNAseq processed on 2017-10-13`、`wrangling_procedure`、`description`、TCGA DCC来源URL、RNA-seq平台和HUGO probeMap字段。元数据说明该对象来自TCGA DCC Level 3 `*.rsem.genes.normalized_results`，经 `log2(x+1)` 变换后处理进入Xena。由此可确认**主体全部366例逐值相同的矩阵数值对应 `log2(norm_count+1)`**。

相邻候选 `.meta`、`.metadata`、`.meta.json`、`.txt` 请求均返回403；Xena `/data/` dataset metadata API返回HTTP 500/NPE。Xena JSON引用的TCGA DCC `DESCRIPTION.txt` URL当前返回HTTP 200，但内容是31KB NCI HTML页面，未返回DESCRIPTION文本。上述探测URL、HTTP状态和响应字段摘要记于[相邻端点探测记录](public_xena_adjacent_metadata_probe.json)；官方数据集JSON单独保存于[HiSeqV2数据集元数据](public_TCGA_LIHC_HiSeqV2_dataset_metadata.json)。


## 外部 CHOL 的同尺度 E4 输入

为使E4外部矩阵与主体LIHC使用同一官方Xena流程，仅构建新的外部输入文件，不修改原始CHOL RDA或既有E4评估结果。官方Xena registry将 `TCGA.CHOL.sampleMap/HiSeqV2` 列为TCGA Bile Duct Cancer (CHOL)的gene-expression数据集；相邻官方JSON记录的 `unit=log2(norm_count+1)`，版本为2017-10-13，wrangler为Xena TCGAscript RNAseq processed on 2017-10-13，来源为TCGA DCC Level 3 `*.rsem.genes.normalized_results`后经`log2(x+1)`处理。该单位/流程与已核验LIHC HiSeqV2对象相同。

外部RDA的35个病例按原`mrna$patient_id`匹配唯一 `patient_id-01` 样本列，35/35均唯一匹配，无缺失或歧义；未根据表达值、miRNA或结局挑选样本。新矩阵保留Xena提供的20,530个基因特征。E0的3,000个模型输入基因中，3,000个均在CHOL Xena矩阵中，缺失0个。`mirna` 对象保持原35×776形状及原内容，`clinical` 保持原0×85空表。

新建的 [external_CHOL_matched_RSEM.rda](../E4_harmonized_inputs/external_CHOL_matched_RSEM.rda) 已做独立复核：35×20,530表达值（718,550个值）均与相应Xena矩阵值完全一致，最大绝对差为0；`mirna` 和 `clinical` 对象均与原RDA逐对象完全相同。来源JSON含版本、单位、源URL、hash、样本映射及输入基因覆盖；样本映射、3,000基因覆盖和核验汇总另以CSV保存。此步骤只生成新输入，不拟合模型，也未重算外部性能。
