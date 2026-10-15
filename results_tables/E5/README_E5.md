# E5 miRNA 位点、宿主表达与评分关系

在仓库根目录运行 `Rscript new/code/biological_annotation.R` 重建结果。脚本读取 E0 canonical 输入、原始 `data/data.rda` 和 GENCODE/UCSC 注释。注释版本、链接、映射规则、筛选范围和文件 MD5 见 `E5_annotation_metadata.json`。

## 位点和候选宿主

846个前体样式输入名称中，837个与 UCSC hg38 wgRna miRNA 位点精确匹配。基因组候选需与位点同染色体、同方向且完整包含位点；GENCODE gene_type 进一步限制为 **protein_coding 或 lncRNA**。miRNA 自身 gene_type=miRNA 记录不计为宿主，所有类型的区间重叠记录另列于 `E5_locus_gene_identity_all_overlaps.csv`。

最终有765条候选宿主基因记录，涉及724个 miRNA 位点，其中484条为 protein_coding、281条为 lncRNA。Rel≥0.7的33项中，28项至少有一个候选宿主。每个候选宿主记录都包含 `in_model_3000` 和 `available_in_full_mrna` 标志。候选汇总见 `E5_host_gene_candidates.csv`，按 miRNA 的候选数、特征覆盖及名单见 `E5_miRNA_identity_expression_score.csv`。

## 全量 mRNA 宿主表达主分析

原始数据含366例和18204个 mRNA 特征。脚本按 patient_id 将原始 mRNA 重排到 canonical train/cal/test 的366行 miRNA 真值顺序；3000个模型输入特征与 canonical native 矩阵最大绝对差为0。宿主表达与 miRNA 真值逐例配对后计算 Spearman 相关，每个候选对记录有限配对数。435个全 mRNA 候选关联对的有效样本数均为366，涉及433个 miRNA 位点；按435项 BH 校正后351项 q<0.05。该分析覆盖了可在18204个 mRNA 特征中按候选宿主符号匹配到的基因。结果见 `E5_host_expression_associations_all_mRNA.csv`。

作为模型输入范围的单独对照，3000个模型特征中有47个 miRNA 位点具备可测候选宿主表达，形成48个关联对，44项 BH q<0.05；结果见 `E5_host_expression_associations_model_inputs_only.csv`。全量 mRNA 分析在模型未使用的特征中增加了可分析宿主位点，不把它们归为缺失。候选表中的逐基因标记说明哪些宿主符号可在完整矩阵匹配、哪些属于模型输入。未能按符号匹配到全矩阵特征的候选仍标记为不可测。

全 mRNA 范围有18个高分条目和415个背景条目具备宿主表达关联。按native表达均值与 log(1+方差)联合十分位分层，仅在同时含高分和背景的共同支持层进行标签置换；这些层包含13个高分条目、40个背景条目和8个分层。按高分条目在各层的比例加权，高分与背景平均绝对宿主相关差为0.0271。每层内无放回交换高分标签并保持层内高分数量，10,000次双侧置换 p=0.647。3000特征范围有6个高分条目和41个背景条目具备关联；共同支持仅含2个高分和6个背景、共2层，差值−0.0465，p=0.473。**两种范围均未显示高分条目的宿主共表达强于表达/方差相似背景的证据。** 全量 mRNA 的检验仍受共同支持范围限制，置换结果依赖层内标签可交换假设，且只反映观察性共表达。

两种范围和共同支持计数见 `E5_host_expression_stratified_label_permutation.csv`。该检验在所有候选宿主表达关联已计算后，分别按范围汇总每个 miRNA 的平均绝对宿主相关，再在共同支持分层内置换高分标签；差值按分层高分数加权。宿主相关的单对 p 值按有效患者配对计算，并在各自候选关联集合内作 BH 校正。

## Rel 与表达丰度、方差

Rel 与 miRNA native 表达均值的点估计 Spearman rho 为0.694；与表达方差为0.762。患者层级 bootstrap（B=2000）结果如下：

| 统计量 | 点估计 rho | Bootstrap 中位数 | 95%百分位区间 | 有效重复数 |
|---|---:|---:|---:|---:|
| Rel 与表达均值 | 0.694 | 0.643 | [0.594, 0.687] | 2000 |
| Rel 与表达方差 | 0.762 | 0.708 | [0.600, 0.759] | 2000 |

每次重采样分别在训练、校准和测试分组内按患者有放回抽样，合并366次抽样后重新计算native表达均值和方差；校准集真实表达与固定模型预测同步抽样，并据此重算 Rel。模型未重训。两项的点估计均高于对应百分位区间上界，bootstrap 分布中位数较低，表明重采样分布有向下偏移；按要求报告的百分位区间不应被解读为新队列的确定性保证。数值见 `E5_score_expression_patient_bootstrap.csv`。旧版把846个 miRNA 项目当作独立单位计算 tiny p 值的文件已移至 `invalid_for_inference/`，不用于推断。

Rel≥0.7的33项中，366例输入表达均值中位数为7.408，阈值以下813项为0.702；方差中位数分别为2.052和0.405。结合患者层级区间，Rel与输入表达均值及方差存在关联，不能仅凭分数把高分集合解释为共同功能或生物学重要性更高。GDC病例数值比对支持miRNA输入为log2(RPM+1)，这里的均值和方差均在该log转换尺度上计算。

高、中、低代表项按校准集 Rel 选择，分别为 hsa-mir-466（0.8746）、hsa-mir-330（0.4997）和 hsa-mir-6882（0）。校准/测试集真实值与预测值数据见 `E5_high_mid_low_example_scatter_data.csv`，共432行，可直接绘制散点图。

## 解释边界与复核文件

UCSC 匹配的是 miRNA 位点/前体式名称，不识别成熟体或5p/3p链。本分析没有查询或推断成熟体经验证靶基因。区间包含仅用于提名宿主候选，不能证明实际转录、加工或调控。数据中没有独立肿瘤纯度或细胞组成字段，本分析未估计纯度。

`invalid_initial/` 保存最初未排除 miRNA 自身 gene 记录及错误重采样方法的结果；`invalid_for_inference/` 保存按 miRNA 项目独立性计算的旧表达关联 tiny p 文件；`superseded_3k_only/` 保存已由明确标注的模型输入范围结果替代的旧输出；`checks/` 保存 B=20 的 bootstrap 数值小检查。正式数字以当前目录的新输出为准。


**来源与单位核验补记（2026-10-14）。** GDC官方病例接口核实主体366例TCGA-LIHC与外部35例TCGA-CHOL；主体mRNA与官方Xena HiSeqV2的366×18,204个数值全部相同，该对象元数据定义为RSEM归一化计数的log2(norm_count+1)，版本2017-10-13。五个GDC病例的数值比对支持两组miRNA为log2(RPM+1)。主体366例OS_time与GDC死亡或末次随访日数全部精确一致，单位为天；项目最初下载日期、全部miRNA文件版本和原阈值预先确定记录仍未恢复。见E0/public_provenance_search.md及逐值核验表。 原E4输入的mRNA量化流程与主体不一致，旧结果及图件保留为诊断记录；同RSEM外部输入与后续结果分别存于results/E4_harmonized_inputs及results/harmonized/E4。
