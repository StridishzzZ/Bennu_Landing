# -*- coding: utf-8 -*-
"""实验5：蒙特卡洛不确定性传播（连续科学参数 ±10% 均匀扰动）。"""

from common import load_sites, save
from src.analysis import monte_carlo_ranking


def main():
    df = load_sites("bennu_sites.csv")
    prob, _ = monte_carlo_ranking(df, n=2000, pct=0.10, seed=42)
    print("=== 贝努四候选区 2000 次抽样下的排名概率（优先级权重）===")
    print(prob.to_string(index=False))
    save(prob, "exp5_monte_carlo_rank_probability.csv")


if __name__ == "__main__":
    main()
