# -*- coding: utf-8 -*-
"""单参数评分函数：把原始参数映射到 [0,1]。"""

import numpy as np
import pandas as pd

from .config import PARAMETERS, SCIENCE_PARAMETERS


def score_continuous(value, lo, hi, direction):
    if pd.isna(value):
        return 0.5  # 缺失取中性值
    if hi <= lo:
        return 0.5
    t = float(np.clip((float(value) - lo) / (hi - lo), 0.0, 1.0))
    return t if direction == "higher_better" else 1.0 - t


def score_categorical(value, mapping):
    if pd.isna(value):
        return 0.0
    return mapping.get(value, 0.0)


def score_binary(value):
    if pd.isna(value):
        return 0.0
    return 1.0 if float(value) > 0.5 else 0.0


def score_parameter(name, value, param):
    kind = param["kind"]
    if kind == "continuous":
        return score_continuous(value, param["lo"], param["hi"], param["direction"])
    if kind == "categorical":
        return score_categorical(value, param["categories"])
    if kind == "binary":
        return score_binary(value)
    raise ValueError(f"未知参数类型: {kind} ({name})")


def score_all_parameters(df):
    """对数据集中每个站点计算全部（精简后）科学参数的评分。"""
    out = pd.DataFrame(index=df.index)
    out["site"] = df["site"]
    for name in SCIENCE_PARAMETERS:
        out[name] = df[name].apply(lambda v: score_parameter(name, v, PARAMETERS[name]))
    return out
