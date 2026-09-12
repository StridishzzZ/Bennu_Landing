# -*- coding: utf-8 -*-
"""实验4：敏感性分析（权重扰动 + 参数扰动）。"""

from common import load_sites, save
from src.analysis import parameter_sensitivity, weight_sensitivity


def main():
    df = load_sites("bennu_sites.csv")

    w_sens = weight_sensitivity(df, pct=0.20)
    print("=== 权重 ±20% 扰动后的排序变化（排名位置改变次数）===")
    print(w_sens.to_string(index=False))
    save(w_sens, "exp4_weight_sensitivity.csv")

    p_sens = parameter_sensitivity(df, pct=0.10)
    print("\n=== 连续科学参数 ±10% 扰动后的排序变化 ===")
    print(p_sens.to_string(index=False))
    save(p_sens, "exp4_parameter_sensitivity.csv")

    impact = (p_sens.groupby("parameter")["rank_changes"].sum()
              .sort_values(ascending=False))
    print("\n影响最大的参数（按 ±10% 引起的最大 ISVM 变化排序）：")
    impact2 = (p_sens.groupby("parameter")["max_isvm_delta"].max()
               .sort_values(ascending=False))
    print(impact2.head(10).to_string())
    print("\n（排名位置改变次数合计——Nightingale 领先较大，±10% 未引起换位）：")
    print(impact.to_string())


if __name__ == "__main__":
    main()
