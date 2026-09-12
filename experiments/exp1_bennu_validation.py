# -*- coding: utf-8 -*-
"""实验1：贝努四候选区复现验证。

验收标准（文献）：Nightingale 为最终采样点，ISVM 应排名第一；
Osprey 为备选点，科学价值应不显著落后，且工程层（可采样/可达）更具优势。
"""

from common import load_sites, save, show_ranking
from src.pipeline import evaluate_sites


def main():
    df = load_sites("bennu_sites.csv")
    res = evaluate_sites(df, weight_scheme="priority")
    ranking = show_ranking("贝努四候选区 ISVM 排序（文献优先级权重）", res)

    print("\n指标与工程层明细：")
    cols = ["site", "C", "M", "G", "T", "isvm",
            "safety_score", "sampleability_score", "deliverability_score",
            "safety_pass", "deliverability_pass", "engineering_pass"]
    print(res[cols].to_string(index=False))

    save(res, "exp1_bennu_priority.csv")
    save(ranking, "exp1_bennu_ranking.csv")
    print("\n结果已保存到 outputs/exp1_bennu_priority.csv")


if __name__ == "__main__":
    main()
