"""Publish verified same-workflow E4 results into Chinese current notes."""
from pathlib import Path
import csv
import re
import json
from datetime import datetime

root = Path(__file__).resolve().parents[2]
out = root / "new/results/harmonized/E4"
def rows(name):
    with (out / name).open(encoding="utf-8-sig", newline="") as f:
        return list(csv.DictReader(f))
ci = {r["metric"]: r for r in rows("E4_complete_external_bootstrap_summary.csv")}
assert len(ci) == 9
assert all(int(float(r["B_requested"])) == 2000 and int(float(r["B_valid"])) == 2000 for r in ci.values())
checks = rows("E4_complete_bootstrap_reproduction_check.csv")
assert len(checks) == 2000
assert max(abs(float(v)) for r in checks for k, v in r.items() if k != "replicate") < 1e-8
metadata = rows("E4_run_metadata.csv")[0]
assert int(metadata["n_patients"]) == 35 and int(metadata["n_overlap_miRNA"]) == 742
assert float(metadata["external_gene_coverage"]) == 1
assert metadata["external_input"] == "new/results/E4_harmonized_inputs/external_CHOL_matched_RSEM.rda"
assert rows("E4_internal_missing_gene_mask_status.csv")[0]["status"] == "not_applicable_all_model_genes_present"
high = next(r for r in rows("E4_Rel_threshold_external_coverage.csv") if r["selection"] == "Rel_ge_0.7_canonical_33")
def n(value):
    return f"{float(value):.5f}".replace("-", "−")
def estimate(metric):
    r = ci[metric]
    return f'{n(r["point_estimate"])}（患者bootstrap百分位区间{n(r["CI_q025"])}至{n(r["CI_q975"])}）'

body = (
    "来源核验确认主体mRNA为log2(RSEM norm_count+1)，原外部mRNA数值支持log2(STAR count+1)，量化与归一化流程不同；旧结果保留在results/E4，作为输入不一致的诊断记录。新的外部输入来自官方Xena CHOL HiSeqV2，与主体采用同RSEM流程，原35例均唯一匹配主肿瘤样本，3,000个模型基因全部覆盖，原miRNA与空clinical对象保持不变。\n\n"
    "固定原PCA_ElasticNet模型、原训练预处理和内部Rel，在35例外部患者、742项共有miRNA上评价，未按外部真值调模型、阈值或评分权重。"
    f'总体平均R²为{estimate("mean_R2")}，中位Pearson为{estimate("median_Pearson")}，平均训练标准化RMSE为{estimate("mean_RMSE")}。'
    "每项区间均有2,000次有效患者重采样；原点估计与bootstrap分布中位数分别保留。\n\n"
    f'内部Rel≥0.7的原33项名单在外部共有32项，其中中位Pearson为{n(high["median_Pearson"])}、中位R²为{n(high["median_R2"])}，{int(high["reconstructable_n"])}/32项同时达到R²>0且Pearson≥0.5。'
    "这些是固定内部名单在外部的描述性点估计，不构成根据外部结果重选的名单。\n\n"
    f'内部Rel与外部RMSE的rho为{estimate("rho_internal_Rel_vs_RMSE")}，与外部R²的rho为{estimate("rho_internal_Rel_vs_R2")}，与外部Pearson的rho为{estimate("rho_internal_Rel_vs_Pearson")}。'
    "总体平均R²仍为负，迁移结论应限定在部分条目的表达预测及排序，不能推广为整体可靠替代实测。外部临床/生存TSV已经找到并核验单位，本轮未执行外部生存或临床效用分析。\n\n"
    "输入基因缺失为0，缺失基因屏蔽诊断记为不适用；原缺失382基因的诊断属于旧输入。等规模简单评分及200个随机名单结果保留，随机名单间变异与患者区间分别解释。证据：results/harmonized/E4中的逐miRNA表、固定阈值表、完整区间及bootstrap复现表；来源：results/E4_harmonized_inputs中的JSON及病例映射。"
)
report = root / "new/补实验阶段结果_中文版.md"
text = report.read_text(encoding="utf-8")
pattern = r"\*\*五、[^\n]*\*\*\n\n(.*?)\n\n(?=\*\*六、)"
match = re.search(pattern, text, re.S)
assert match, "E4 section boundary missing"
archive = root / "new/results/E4/旧输入阶段报告摘录.md"
if not archive.exists():
    archive.write_text("**原外部输入流程不一致的阶段报告摘录，仅供诊断追溯。**\n\n" + match.group(1) + "\n", encoding="utf-8")
text = re.sub(pattern, lambda _: "**五、同RSEM流程的外部表达迁移**\n\n" + body + "\n\n", text, count=1, flags=re.S)
text = text.replace("E4按同量化流程输入重算中；", "E4同量化流程评价及患者B=2000已完成；")
report.write_text(text, encoding="utf-8")

caption = (
    "**F_E4_external_transfer：同RSEM流程的外部表达迁移。** 外部35例TCGA-CHOL mRNA来自官方Xena HiSeqV2，与主体采用log2(RSEM norm_count+1)流程；3,000个模型基因全部覆盖，原外部miRNA真值保持不变，共有742项。A显示固定内部Rel与外部逐miRNA Pearson，虚线为操作性标签的Pearson分量0.5，完整标签还要求R²>0。B显示全部742项与内部固定阈值所选170、72、32项的外部R²分布；箱体为25%至75%分位数，中线为中位数，须为1.5倍四分位距内范围，离群点未画出，所有条目均进入统计。"
    f'全部条目平均R²={n(ci["mean_R2"]["point_estimate"])}，高分32项中位R²={n(high["median_R2"])}、中位Pearson={n(high["median_Pearson"])}；'
    "有限子集的表现不能代替全体预测评价。模型、训练预处理和内部评分均固定。患者B=2000的95%百分位区间另见完整区间表，箱线为条目间分布，不代表患者bootstrap区间；来源：results/harmonized/E4。旧输入流程不一致的图件保存在figures/diagnostic_original_input。"
)
p = root / "new/figures/图件说明_中文版.md"
s = p.read_text(encoding="utf-8")
s = re.sub(r"\*\*F_E4_external_transfer：[^\n]*", lambda _: caption, s, count=1)
p.write_text(s, encoding="utf-8")
p = root / "new/EXPERIMENT_TRACKER.md"
s = p.read_text(encoding="utf-8")
s = re.sub(r"\| E4 \|[^\n]*", "| E4 | 同RSEM流程固定评价及B=2000完成 | results/harmonized/E4；主与complete bootstrap的2000条复现差<1e-8；35例742项，3000/3000输入覆盖 | 旧results/E4保留为输入不一致诊断；无外部生存实验 |", s)
p.write_text(s, encoding="utf-8")
p = root / "new/补实验计划书_中文版.md"
s = p.read_text(encoding="utf-8").replace("原输入存在量化流程差异；同RSEM流程35例输入已准备，固定评价/B=2000重算中", "同RSEM流程35例输入的固定模型评价与B=2000完成；旧结果为输入差异诊断")
p.write_text(s, encoding="utf-8")
(out / "README_E4_harmonized.md").write_text("**同RSEM量化流程的外部固定模型评价。**\n\n" + body + "\n\n**复现命令。** 在工作区根目录执行：\n\n```powershell\nRscript new/code/score_external_analysis.R --section E4 --bootstrap 2000 --external-input new/results/E4_harmonized_inputs/external_CHOL_matched_RSEM.rda --output-root new/results/harmonized\nRscript new/code/uncertainty_completion.R --section external --bootstrap 2000 --external-input new/results/E4_harmonized_inputs/external_CHOL_matched_RSEM.rda --output-root new/results/harmonized\n```\n\nE4_historical_prediction_reproduction_check.csv在本输入下表示与历史旧输入预测指标的差异，不要求为0；它不参与新结果估计。\n", encoding="utf-8")
(out / "E4_verification.json").write_text(json.dumps({
    "verified_at": datetime.now().astimezone().isoformat(),
    "input": metadata["external_input"], "patients": 35, "shared_miRNAs": 742,
    "model_gene_coverage": 1, "missing_gene_mask": "not_applicable",
    "bootstrap_requested": 2000, "bootstrap_all_nine_metrics_valid": 2000,
    "bootstrap_reproduction_rows": len(checks),
    "bootstrap_reproduction_max_abs_difference": max(abs(float(v)) for r in checks for k, v in r.items() if k != "replicate"),
    "model_and_score": "original E0 PCA_ElasticNet and internal Rel; no external tuning",
    "external_survival_analysis": "not performed",
}, ensure_ascii=False, indent=2), encoding="utf-8")
print("Harmonized E4 B=2000 and reproduction verified; Chinese current notes updated.")
