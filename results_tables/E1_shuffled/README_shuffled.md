**20次训练配对置换负对照**

正式种子20261014至20261033均完成，结果为20个RDS、20个单次summary和E1_shuffled_20_permutation_summary.csv。所有置换保持原E0训练222例、评分72例、测试72例；仅训练mRNA与miRNA患者对应被打乱，置换源在222例训练患者内构成双射，评分和测试真实值不变。

每次采用五折训练内部lambda选择，预处理按折重新拟合；模型按实际原代码使用多响应Ridge，不能沿用稿件中Shuffled ElasticNet的名称。使用与glmnet alpha=0目标一致的谱分解，lambda路径和预测已与glmnet核验；首路径的零预测遵循glmnet的null-path约定，该点不解释为有限lambda解析Ridge解。选择结果和实际路径参数均保存。

20次中13次选择首路径，预测等于训练均值，Pearson因常数预测记为NA；其余7次预测非恒定。测试平均R²范围为−0.02138至−0.01842，训练SD归一化的整体RMSE为1.01326至1.01464，较训练均值参照的RMSE差为−0.000605至0.000778。**配对打乱后没有出现稳定的实质预测改善**；这是固定切分空信号参照，不能据20次置换证明所有原预测关系具有生物因果意义。

完整性核验确认20个目标种子、全部完成状态、双射置换、固定患者对应、846个有限测试预测与每miRNA OOF归一化误差。独立代码核验见new/SHUFFLED_CODE_REVIEW.md，单次真实维度pilot保存在results/_pilot，未纳入正式20次数量。执行命令：Rscript new/code/shuffled_baseline.R all。
