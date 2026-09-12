# -*- coding: utf-8 -*-
"""实验2：龙宫七候选区多级下选复现。

文献（Yabuta et al. 2019）：科学+安全+可采样分数求和，L08/L07/M04 进入 TD1 候选；
后续更高分辨率观测推动 L08 → L08-E1（TD1）与 L08-B（SCI 人工陨石坑后 TD2）。
"""

from common import load_sites, save, show_ranking
from src.pipeline import evaluate_sites


def main():
    df = load_sites("ryugu_sites.csv")
    res = evaluate_sites(df, weight_scheme="priority", gate_radius=3.0, gate_nav=2.7)

    show_ranking("龙宫七候选区 ISVM（纯科学价值）", res, "isvm")
    show_ranking("龙宫七候选区 综合分（科学+安全+可采样+可达，参照文献求和做法）",
                 res, "combined")

    print("\n指标与工程层明细：")
    cols = ["site", "C", "M", "G", "T", "isvm", "combined",
            "safety_score", "sampleability_score", "deliverability_score"]
    print(res[cols].to_string(index=False))
    save(res, "exp2_ryugu_seven.csv")

    sub = load_sites("ryugu_subsites.csv")
    res_sub = evaluate_sites(sub, weight_scheme="priority", gate_radius=3.0, gate_nav=2.7)
    show_ranking("龙宫 L08 子区（TD1: L08-E1 vs SCI 后 TD2: L08-B）", res_sub, "combined")
    print(res_sub[cols].to_string(index=False))
    save(res_sub, "exp2_ryugu_subsites.csv")


if __name__ == "__main__":
    main()
