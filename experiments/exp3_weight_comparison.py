# -*- coding: utf-8 -*-
"""实验3：三种权重方案对比（文献优先级 / 等权 / 熵权）。"""

import numpy as np

from common import load_sites, save, show_ranking
from src.pipeline import evaluate_sites
from src.indicators import indicator_scores
from src.scoring import score_all_parameters
from src.weights import entropy_weights, get_weights


def spearman(a, b):
    ra = {s: i for i, s in enumerate(a)}
    rb = {s: i for i, s in enumerate(b)}
    order_a = [ra[s] for s in a]
    order_b = [rb[s] for s in b]
    n = len(a)
    d2 = sum((x - y) ** 2 for x, y in zip(order_a, order_b))
    return 1.0 - 6.0 * d2 / (n * (n * n - 1)) if n > 1 else 1.0


def main():
    df = load_sites("bennu_sites.csv")
    orders = {}
    frames = []
    for scheme in ["priority", "equal", "entropy"]:
        res = evaluate_sites(df, weight_scheme=scheme)
        order = list(res.sort_values("isvm", ascending=False)["site"])
        orders[scheme] = order
        w = get_weights(scheme, indicator_scores(score_all_parameters(df)))
        print(f"权重方案 {scheme}: {w}")
        show_ranking(f"贝努四候选区 ISVM 排序（{scheme}）", res)
        frames.append(res[["site", "isvm"]].rename(columns={"isvm": f"isvm_{scheme}"}))
        save(res, f"exp3_bennu_{scheme}.csv")

    merged = frames[0].merge(frames[1], on="site").merge(frames[2], on="site")
    save(merged, "exp3_weight_comparison.csv")

    print("\n排序一致性（Spearman）：")
    for a in orders:
        for b in orders:
            if a < b:
                print(f"  {a} vs {b}: {spearman(orders[a], orders[b]):.3f}")


if __name__ == "__main__":
    main()
