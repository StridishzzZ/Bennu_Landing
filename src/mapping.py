# -*- coding: utf-8 -*-
"""把 ISVM 科学价值评分结果空间化到 Bennu 全幅地图。

思路
----
评分本身是「每个候选区一个值」（`src/pipeline.py`），要在全球地图上表达，
必须定义两个空间要素：

1. **背景场**：全球背景基线（`data/global_baseline.csv`，文献锚定的全球平均性质），
   经同一套评分函数得到背景指标分。文献一再强调 Bennu 表面矿物/水合近似均匀
   （Hamilton et al. 2019；Walsh et al. 2022），因此背景用近似均匀的基线是合理的。
2. **候选区足迹**：以文献坐标为中心、按足迹半径展开的高斯核
   （`data/site_footprints.csv`，半径可改）。核值 α∈[0,1] 表示该像元被候选区
   「覆盖」的程度，像元值 = (1-α)·背景 + α·站点值；多个站点核叠加时按权重归一。

这样得到的图既保留了全球均匀背景，又能看出四个候选区的相对高低，
并且换用真实栅格或修正站点坐标后可直接重算。
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from .indicators import indicator_scores
from .raster import gaussian_footprint
from .scoring import score_all_parameters

INDICATORS = ["C", "M", "G", "T"]

# 四个通道的含义（对应 ISVM 四部分 / 文献四张科学价值图）
INDICATOR_CHANNELS = {
    "C": ("SVCCM", "化学组成科学价值（C 层）"),
    "M": ("SVMM", "矿物学科学价值（M 层）"),
    "G": ("SVGFM", "地质特征与新鲜度科学价值（G 层）"),
    "T": ("SVTM", "温度/热环境科学价值（T 层）"),
}


def background_indicator_scores(baseline_df):
    """全球背景基线参数 → 四指标背景得分。"""
    score_df = score_all_parameters(baseline_df)
    ind = indicator_scores(score_df).reset_index(drop=True)
    return {k: float(ind.iloc[0][k]) for k in INDICATORS}


def site_indicator_table(site_scores):
    """从 `evaluate_sites` 的输出中取出站点-指标表。"""
    cols = ["site"] + INDICATORS
    missing = [c for c in cols if c not in site_scores.columns]
    if missing:
        raise ValueError(f"评分结果缺少列：{missing}")
    return site_scores[cols].copy()


def render_indicator_fields(grid, site_scores, footprints, background,
                            sigma_scale=1.0, default_radius_m=25.0):
    """把站点指标分铺到全球格网，返回 ``{指标: 2D 数组}``。

    Parameters
    ----------
    site_scores : DataFrame
        含 ``site, C, M, G, T`` 列（`evaluate_sites` 的输出）。
    footprints : DataFrame
        含 ``site, lat_deg, lon_deg, footprint_radius_m`` 列。
    background : dict
        四指标背景值。
    sigma_scale : float
        足迹半径缩放，用于敏感性检查（默认 1.0）。
    """
    table = site_indicator_table(site_scores)
    fp = footprints.copy()
    fp["site"] = fp["site"].astype(str)
    table["site"] = table["site"].astype(str)
    merged = table.merge(fp, on="site", how="left")
    if merged["lat_deg"].isna().any():
        missing = sorted(merged.loc[merged["lat_deg"].isna(), "site"].unique())
        raise ValueError(f"以下站点缺少 footprint 坐标：{missing}")

    height, width = grid.shape
    acc_w = np.zeros((height, width), dtype="float64")
    acc_v = {k: np.zeros((height, width), dtype="float64") for k in INDICATORS}

    for _, row in merged.iterrows():
        radius = row.get("footprint_radius_m")
        radius = float(radius) if pd.notna(radius) else float(default_radius_m)
        sigma_m = radius * float(sigma_scale)
        kernel = gaussian_footprint(grid, row["lat_deg"], row["lon_deg"], sigma_m)
        acc_w += kernel
        for key in INDICATORS:
            acc_v[key] += kernel * float(row[key])

    alpha = np.clip(acc_w, 0.0, 1.0)
    denom = np.maximum(acc_w, 1e-12)

    fields = {}
    for key in INDICATORS:
        blended = acc_v[key] / denom
        fields[key] = (background[key] * (1.0 - alpha) + blended * alpha).astype("float32")
    return fields, alpha.astype("float32")


def compose_isvm(fields, weights):
    """按权重把四通道合成为 ISVM 综合图。"""
    total = None
    for key in INDICATORS:
        term = fields[key] * float(weights[key])
        total = term if total is None else total + term
    return total.astype("float32")


def footprint_area_summary(grid, footprints, sigma_scale=1.0):
    """输出各候选区足迹的像元数与等效半径（公里级），用于自检。"""
    rows = []
    m_per_px_lat = grid.meter_per_degree_lat * grid.res_deg
    for _, row in footprints.iterrows():
        sigma_m = float(row["footprint_radius_m"]) * float(sigma_scale)
        kernel = gaussian_footprint(grid, row["lat_deg"], row["lon_deg"], sigma_m)
        # 高斯 2σ 等效面积（像元数），按纬度修正经向像元尺度
        m_per_px_lon = m_per_px_lat * max(np.cos(np.radians(float(row["lat_deg"]))), 1e-6)
        px_area_m2 = m_per_px_lat * m_per_px_lon
        rows.append({
            "site": row["site"],
            "lat_deg": float(row["lat_deg"]),
            "lon_deg": float(row["lon_deg"]),
            "footprint_radius_m": float(row["footprint_radius_m"]),
            "kernel_max": float(kernel.max()),
            "n_pixels_alpha_ge_0p5": int((kernel >= 0.5).sum()),
            "area_alpha_ge_0p5_m2": float((kernel >= 0.5).sum() * px_area_m2),
            "m_per_px_lon": float(m_per_px_lon),
            "m_per_px_lat": float(m_per_px_lat),
        })
    return pd.DataFrame(rows)
