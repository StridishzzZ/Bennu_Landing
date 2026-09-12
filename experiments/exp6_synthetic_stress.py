# -*- coding: utf-8 -*-
"""实验6：合成数据压力测试（极端值、缺失值、全零/全满）。"""

import numpy as np
import pandas as pd

from common import load_sites, save, show_ranking
from src.pipeline import evaluate_sites


def make_synthetic(n=30, seed=7):
    rng = np.random.default_rng(seed)
    template = load_sites("parameter_template.csv")
    rows = []
    for i in range(n):
        r = template.iloc[0].copy()
        r["site"] = f"SYN-{i:02d}"
        r["organic_abundance"] = rng.uniform(0, 1)
        r["ch2_ch3"] = rng.uniform(0.5, 3.0)
        r["carbon_content_pct"] = rng.uniform(0, 10)
        r["mineral_classes_detected"] = rng.uniform(0, 6)
        r["hydration_depth"] = rng.uniform(0, 1)
        r["space_weathering_proxy"] = rng.uniform(0, 1)
        r["crater_freshness"] = rng.choice(["recent", "paleo", "none"])
        r["activity_level"] = rng.choice(["active", "recent", "paleo", "none"])
        r["psfd_class"] = rng.choice(["fine", "mix", "coarse"])
        r["temp_rel_global"] = rng.choice(["colder", "average", "hotter"])
        r["tmax_k"] = rng.uniform(280, 350)
        r["slope_deg"] = rng.uniform(0, 20)
        r["boulder_coverage"] = rng.uniform(0, 0.30)
        r["roughness"] = rng.uniform(0, 1)
        r["frac_le2cm"] = rng.uniform(0, 1)
        r["max_particle_cm"] = rng.uniform(1, 10)
        r["tagsam_tilt_deg"] = rng.uniform(0, 20)
        r["target_radius_m"] = rng.uniform(2, 25)
        r["nav_error_m"] = rng.uniform(0.5, 5)
        rows.append(r)
    return pd.DataFrame(rows)


def main():
    syn = make_synthetic()
    res = evaluate_sites(syn, weight_scheme="priority", gate_radius=3.0)

    assert res["isvm"].between(0, 1).all(), "ISVM 越界"
    assert res["C"].between(0, 1).all() and res["M"].between(0, 1).all()
    assert res["G"].between(0, 1).all() and res["T"].between(0, 1).all()
    assert res.notna().all().all(), "存在 NaN 输出"
    print(f"合成数据 {len(syn)} 个站点：所有得分均在 [0,1]，无 NaN。")
    show_ranking("合成数据 ISVM 前 10 名", res)

    zero = syn.iloc[0:1].copy()
    full = syn.iloc[0:1].copy()
    for col in ["organic_abundance", "carbon_content_pct", "hydration_depth"]:
        zero[col] = 0.0
        full[col] = 1.0
    zero["mineral_classes_detected"] = 0.0
    full["mineral_classes_detected"] = 6.0
    zero["ch2_ch3"] = 3.0
    full["ch2_ch3"] = 0.5
    zero["space_weathering_proxy"] = 1.0
    full["space_weathering_proxy"] = 0.0
    zero["crater_freshness"] = "none"
    full["crater_freshness"] = "recent"
    zero["activity_level"] = "none"
    full["activity_level"] = "active"
    zero["psfd_class"] = "coarse"
    full["psfd_class"] = "fine"
    zero["temp_rel_global"] = "hotter"
    full["temp_rel_global"] = "colder"

    z = evaluate_sites(zero, weight_scheme="priority").iloc[0]
    f = evaluate_sites(full, weight_scheme="priority").iloc[0]
    print(f"\n全零站点 ISVM = {z['isvm']:.4f}（期望接近 0，实际 C/M/G/T = "
          f"{z['C']:.2f}/{z['M']:.2f}/{z['G']:.2f}/{z['T']:.2f}）")
    print(f"全满站点 ISVM = {f['isvm']:.4f}（期望接近 1）")

    miss = syn.iloc[0:1].copy()
    miss["organic_abundance"] = np.nan
    mr = evaluate_sites(miss, weight_scheme="priority").iloc[0]
    print(f"缺失参数站点可正常运行，ISVM = {mr['isvm']:.4f}（缺失项按中性值 0.5 处理）")

    save(res, "exp6_synthetic_ranking.csv")


if __name__ == "__main__":
    main()
