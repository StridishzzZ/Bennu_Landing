# -*- coding: utf-8 -*-
"""敏感性分析与蒙特卡洛不确定性传播。"""

import numpy as np
import pandas as pd

from .pipeline import evaluate_sites
from .config import SCIENCE_PARAMETERS, PARAMETERS
from .weights import PRIORITY


def parameter_sensitivity(df, pct=0.10, weight_scheme="priority", seed=0):
    """对每个连续型科学参数做 ±pct 单因素扰动，统计排序变化。"""
    rng = np.random.default_rng(seed)
    base = evaluate_sites(df, weight_scheme=weight_scheme)
    base_rank = base[["site", "isvm"]].sort_values("isvm", ascending=False).reset_index(drop=True)
    base_order = list(base_rank["site"])
    records = []
    for name in SCIENCE_PARAMETERS:
        if PARAMETERS[name]["kind"] != "continuous":
            continue
        for sign in (+1, -1):
            d = df.copy()
            d[name] = d[name] * (1.0 + sign * pct)
            res = evaluate_sites(d, weight_scheme=weight_scheme)
            order = list(res.sort_values("isvm", ascending=False)["site"])
            changes = sum(a != b for a, b in zip(base_order, order))
            delta = (res.set_index("site")["isvm"] - base.set_index("site")["isvm"]).abs().max()
            records.append({"parameter": name, "sign": sign,
                            "rank_changes": changes,
                            "max_isvm_delta": float(delta),
                            "order": "/".join(order)})
    return pd.DataFrame(records)


def weight_sensitivity(df, pct=0.20, base_scheme="priority", seed=0):
    """扰动各指标权重（再归一化），统计排序变化。"""
    rng = np.random.default_rng(seed)
    base = evaluate_sites(df, weight_scheme=base_scheme)
    base_order = list(base.sort_values("isvm", ascending=False)["site"])
    records = []
    for c in ["C", "M", "G", "T"]:
        for sign in (+1, -1):
            w = dict(PRIORITY)
            w[c] = w[c] * (1.0 + sign * pct)
            s = sum(w.values())
            w = {k: v / s for k, v in w.items()}
            res = evaluate_sites(df, weight_scheme="priority")
            res["isvm"] = sum(res[k] * w[k] for k in ["C", "M", "G", "T"])
            order = list(res.sort_values("isvm", ascending=False)["site"])
            changes = sum(a != b for a, b in zip(base_order, order))
            records.append({"indicator": c, "sign": sign,
                            "weight": round(w[c], 3), "rank_changes": changes})
    return pd.DataFrame(records)


def monte_carlo_ranking(df, n=2000, pct=0.10, weight_scheme="priority", seed=42):
    """对连续型科学参数做 ±pct 均匀扰动，输出各站点排名概率。"""
    rng = np.random.default_rng(seed)
    sites = list(df["site"])
    rank_counter = {s: np.zeros(len(sites)) for s in sites}
    continuous = [n for n in SCIENCE_PARAMETERS if PARAMETERS[n]["kind"] == "continuous"]
    for _ in range(n):
        d = df.copy()
        for name in continuous:
            lo = d[name].min()
            hi = d[name].max()
            scale = (hi - lo) * pct if hi > lo else max(abs(hi) * pct, 1e-6)
            noise = rng.uniform(-scale, scale, size=len(d))
            d[name] = d[name] + noise
            d[name] = d[name].clip(lower=0)
        res = evaluate_sites(d, weight_scheme=weight_scheme)
        order = list(res.sort_values("isvm", ascending=False)["site"])
        for rank, s in enumerate(order):
            rank_counter[s][rank] += 1
    prob = pd.DataFrame({s: rank_counter[s] / n for s in sites})
    prob.insert(0, "rank", np.arange(1, len(sites) + 1))
    return prob, rank_counter
