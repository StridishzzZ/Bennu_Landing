# -*- coding: utf-8 -*-
"""ISVM 综合评价流水线 + 工程层门控。"""

import numpy as np
import pandas as pd

from .scoring import score_all_parameters
from .indicators import indicator_scores
from .weights import get_weights


def _clip(v, lo, hi):
    return float(np.clip((float(v) - lo) / max(hi - lo, 1e-12), 0.0, 1.0))


def engineering_scores(df, gate_slope=14.0, gate_boulder=0.20, gate_tmax=350.0,
                        gate_tilt=14.0, gate_radius=5.0, gate_nav=2.7):
    rows = []
    for _, r in df.iterrows():
        slope_s = 1.0 - _clip(r["slope_deg"], 0.0, gate_slope)
        boulder_s = 1.0 - _clip(r["boulder_coverage"], 0.0, 0.25)
        rough_s = 1.0 - _clip(r["roughness"], 0.0, 1.0)
        tmax_s = 1.0 - _clip(r["tmax_k"], 250.0, 350.0)
        safety_score = float(np.mean([slope_s, boulder_s, rough_s, tmax_s]))
        radius_s = _clip(r["target_radius_m"], 0.0, 25.0)
        nav_s = 1.0 - _clip(r["nav_error_m"], 0.0, 5.0)
        deliverability_score = float(np.mean([radius_s, nav_s]))
        frac_s = float(np.clip(float(r["frac_le2cm"]), 0.0, 1.0))
        size_s = 1.0 if float(r["max_particle_cm"]) <= 2.0 else max(
            0.0, 1.0 - (float(r["max_particle_cm"]) - 2.0) / 5.0)
        tilt_s = 1.0 - _clip(r["tagsam_tilt_deg"], 0.0, gate_tilt)
        sampleability_score = float(np.mean([frac_s, size_s, tilt_s]))
        safety_pass = all([
            float(r["slope_deg"]) < gate_slope,
            float(r["boulder_coverage"]) < gate_boulder,
            float(r["tmax_k"]) < gate_tmax,
            float(r["tagsam_tilt_deg"]) < gate_tilt,
        ])
        deliverability_pass = all([
            float(r["target_radius_m"]) >= gate_radius,
            float(r["nav_error_m"]) <= gate_nav,
        ])
        rows.append({
            "site": r["site"],
            "safety_score": safety_score,
            "sampleability_score": sampleability_score,
            "deliverability_score": deliverability_score,
            "safety_pass": safety_pass,
            "deliverability_pass": deliverability_pass,
            "engineering_pass": safety_pass and deliverability_pass,
        })
    return pd.DataFrame(rows)


def evaluate_sites(df, weight_scheme="priority", include_engineering=True, **eng_kwargs):
    """返回每个候选区的指标分、ISVM 与工程层得分。"""
    score_df = score_all_parameters(df)
    ind = indicator_scores(score_df).reset_index()
    weights = get_weights(weight_scheme, ind.set_index("site"))
    out = ind.copy()
    out["isvm"] = sum(out[c] * weights[c] for c in ["C", "M", "G", "T"])
    out["weight_scheme"] = weight_scheme
    if include_engineering:
        eng = engineering_scores(df, **eng_kwargs)
        out = out.merge(eng, on="site", how="left")
        out["combined"] = (out["isvm"] + out["safety_score"]
                           + out["sampleability_score"] + out["deliverability_score"])
    return out.sort_values("isvm", ascending=False).reset_index(drop=True)


def rank_sites(df, score_col="isvm", ascending=False):
    tmp = df[["site", score_col]].copy()
    tmp["rank"] = tmp[score_col].rank(ascending=ascending, method="min").astype(int)
    return tmp.sort_values(score_col, ascending=ascending).reset_index(drop=True)
