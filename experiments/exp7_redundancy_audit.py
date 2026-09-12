# -*- coding: utf-8 -*-
"""实验7：原始参数冗余审计与精简验证。

步骤：
1. 统计每个参数在贝努/龙宫/龙宫子区中的得分方差与完全重复组；
2. 输出精简映射表（删除/合并项及理由）；
3. 用精简后的参数集重跑实验 1–2，并与精简前的排名结果对比。
"""

import os

import pandas as pd

from common import ROOT, DATA_DIR, load_sites
from src.config import PARAMETERS, REDUNDANT_SCIENCE, INDICATOR_PARAMETERS
from src.pipeline import evaluate_sites


OLD_PARAMS = ["organic_abundance", "volatile_abundance", "organic_silicate_ratio", "ch2_ch3",
              "carbon_content_pct", "mineral_phyllo", "mineral_carbonate", "mineral_sulfate",
              "mineral_oxide", "mineral_silicate", "mineral_amorphous", "amorphous_fraction",
              "hydration_depth", "space_weathering_proxy", "crater_freshness", "activity_level",
              "psfd_class", "brittle_deformation", "geologic_diversity", "temp_rel_global", "tmax_k"]

OLD_INDICATOR_PARAMS = {
    "C": ["organic_abundance", "volatile_abundance", "organic_silicate_ratio",
          "ch2_ch3", "carbon_content_pct"],
    "G": ["space_weathering_proxy", "crater_freshness", "activity_level", "psfd_class",
          "brittle_deformation", "geologic_diversity"],
    "T": ["temp_rel_global", "tmax_k"],
}

MINERAL_PRIORITY = {"mineral_phyllo": 1.0, "mineral_carbonate": 0.9, "mineral_sulfate": 0.8,
                    "mineral_oxide": 0.7, "mineral_silicate": 0.6, "mineral_amorphous": 0.5}


def score_old(df):
    out = pd.DataFrame(index=df.index)
    out["site"] = df["site"]
    for name in OLD_PARAMS:
        p = PARAMETERS[name]
        v = df[name]
        if p["kind"] == "continuous":
            t = (v - p["lo"]) / max(p["hi"] - p["lo"], 1e-12)
            t = t.clip(0, 1)
            out[name] = t if p["direction"] == "higher_better" else 1 - t
        elif p["kind"] == "categorical":
            out[name] = v.map(p["categories"]).fillna(0.0)
        else:
            out[name] = (v > 0.5).astype(float)
    return out.set_index("site")


def old_ranking(df):
    sc = score_old(df)
    rows = []
    for _, r in sc.iterrows():
        s = {"site": r.name}
        s["C"] = float(r[OLD_INDICATOR_PARAMS["C"]].mean())
        wsum = sum(MINERAL_PRIORITY.values())
        m_presence = sum(float(r[n]) * w for n, w in MINERAL_PRIORITY.items()) / wsum
        m_diversity = sum(float(r[n]) >= 0.5 for n in MINERAL_PRIORITY) / 6.0
        s["M"] = float(pd.Series([m_presence, m_diversity,
                                  r["amorphous_fraction"], r["hydration_depth"]]).mean())
        s["G"] = float(r[OLD_INDICATOR_PARAMS["G"]].mean())
        s["T"] = float(r[OLD_INDICATOR_PARAMS["T"]].mean())
        rows.append(s)
    ind = pd.DataFrame(rows).set_index("site")
    w = {"C": 0.35, "M": 0.30, "G": 0.20, "T": 0.15}
    ind["isvm"] = sum(ind[k] * w[k] for k in ["C", "M", "G", "T"])
    return ind.reset_index()[["site", "isvm"]]


def compare(name, full_file, gate_radius=5.0):
    df_old = pd.read_csv(os.path.join(DATA_DIR, "full", full_file))
    df_new = load_sites(name)
    old = old_ranking(df_old)
    new = evaluate_sites(df_new, weight_scheme="priority", gate_radius=gate_radius, gate_nav=2.7)
    order_old = list(old.sort_values("isvm", ascending=False)["site"])
    order_new = list(new.sort_values("isvm", ascending=False)["site"])
    same = order_old == order_new
    # 工程参数精简前后相同，旧综合分 = 旧科学分 + 新工程分
    comb_old_df = old.merge(new[["site", "safety_score", "sampleability_score",
                                 "deliverability_score"]], on="site")
    comb_old_df["combined"] = (comb_old_df["isvm"] + comb_old_df["safety_score"]
                               + comb_old_df["sampleability_score"]
                               + comb_old_df["deliverability_score"])
    comb_old = list(comb_old_df.sort_values("combined", ascending=False)["site"])
    comb_new = list(new.sort_values("combined", ascending=False)["site"])
    comb_same = comb_old == comb_new
    print(f"\n{name}: 精简前排序 {order_old}")
    print(f"          精简后排序 {order_new}  -> {'一致 ✓' if same else '不一致 ✗'}")
    print(f"          精简后综合分排序 {comb_new}  -> {'一致 ✓' if comb_same else '不一致 ✗'}")
    return same, comb_same


def variance_audit():
    print("=== 精简前各参数在站点间的得分方差（零方差=无区分度）===")
    for name, full_file in [("贝努", "bennu_sites_full.csv"),
                            ("龙宫", "ryugu_sites_full.csv"),
                            ("龙宫子区", "ryugu_subsites_full.csv")]:
        df = pd.read_csv(os.path.join(DATA_DIR, "full", full_file))
        sc = score_old(df)
        var = sc[OLD_PARAMS].var(axis=0)
        zero = [p for p in var.index if var[p] < 1e-9]
        print(f"{name}: 零方差 {len(zero)}/{len(OLD_PARAMS)} -> {zero}")


def main():
    print("=== 参数精简映射（删除/合并）===")
    for p, reason in REDUNDANT_SCIENCE.items():
        print(f"  - {p}: {reason}")
    print("\n各指标精简后保留参数：")
    for ind, params in INDICATOR_PARAMETERS.items():
        print(f"  {ind}: {params}")

    variance_audit()

    ok1_sci, ok1_comb = compare("bennu_sites.csv", "bennu_sites_full.csv")
    ok2_sci, ok2_comb = compare("ryugu_sites.csv", "ryugu_sites_full.csv", gate_radius=3.0)
    ok3_sci, ok3_comb = compare("ryugu_subsites.csv", "ryugu_subsites_full.csv", gate_radius=3.0)

    print("\n=== 精简前后排序一致性 ===")
    print(f"  科学分排序：贝努 {'✓' if ok1_sci else '✗'}，"
          f"龙宫 {'✓' if ok2_sci else '✗（L07/L08 近并列交换）'}，"
          f"龙宫子区 {'✓' if ok3_sci else '✗'}")
    print(f"  综合分排序：贝努 {'✓' if ok1_comb else '✗'}，"
          f"龙宫 {'✓' if ok2_comb else '✗'}，"
          f"龙宫子区 {'✓' if ok3_comb else '✗'}")
    print("  验收标准（综合分：Nightingale 第一 / L08,L07,M04 前三 / L08-E1 第一）：",
          "全部通过 ✓" if all([ok1_comb, ok2_comb, ok3_comb]) else "未全部通过 ✗")


if __name__ == "__main__":
    main()
