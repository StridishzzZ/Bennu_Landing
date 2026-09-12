# -*- coding: utf-8 -*-
"""三种权重方案：文献优先级、等权、熵权。"""

import numpy as np

# 按 Nakamura-Messenger 等列出的任务优先级（SVCCM>SVMM>SVGFM>SVTM）递减，
# 数值为设计假设，实验 3 会对比其对排序的影响。
PRIORITY = {"C": 0.35, "M": 0.30, "G": 0.20, "T": 0.15}
EQUAL = {"C": 0.25, "M": 0.25, "G": 0.25, "T": 0.25}


def entropy_weights(ind_df):
    """经典熵权法：指标得分离散度越大权重越高。"""
    x = ind_df[["C", "M", "G", "T"]].values.astype(float)
    n = x.shape[0]
    if n < 2:
        return dict(EQUAL)
    # 列归一化：每个指标列内求和为 1，保证熵在 [0,1]
    p = x / np.maximum(x.sum(axis=0, keepdims=True), 1e-12)
    e = -np.sum(p * np.log(np.maximum(p, 1e-12)), axis=0) / np.log(n)
    d = 1.0 - e
    w = d / np.maximum(d.sum(), 1e-12)
    return dict(zip(["C", "M", "G", "T"], w))


def get_weights(scheme, ind_df=None):
    if scheme == "priority":
        return dict(PRIORITY)
    if scheme == "equal":
        return dict(EQUAL)
    if scheme == "entropy":
        if ind_df is None:
            raise ValueError("熵权需要指标得分矩阵")
        return entropy_weights(ind_df)
    raise ValueError(f"未知权重方案: {scheme}")
