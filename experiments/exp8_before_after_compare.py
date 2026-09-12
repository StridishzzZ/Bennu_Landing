# -*- coding: utf-8 -*-
"""实验8：精简前后对比（复现精简前结果 + 生成对比图）。

精简前结果用 data/full/ 归档数据 + 原始聚合逻辑复现（确定性计算，与原始运行一致），
写入 outputs/before_reduction/；对比图写入 outputs/figures/。
"""

import os
import sys

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from common import ROOT, DATA_DIR, OUTPUT_DIR, load_sites
from src.config import PARAMETERS, SCIENCE_PARAMETERS, INDICATOR_PARAMETERS
from src.pipeline import evaluate_sites

plt.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei"]
plt.rcParams["axes.unicode_minus"] = False
plt.rcParams["figure.dpi"] = 150

BEFORE_DIR = os.path.join(OUTPUT_DIR, "before_reduction")
FIG_DIR = os.path.join(OUTPUT_DIR, "figures")
os.makedirs(BEFORE_DIR, exist_ok=True)
os.makedirs(FIG_DIR, exist_ok=True)

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
WEIGHTS = {"C": 0.35, "M": 0.30, "G": 0.20, "T": 0.15}


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


def old_indicators(sc):
    rows = []
    for site, r in sc.iterrows():
        s = {"site": site}
        s["C"] = float(r[OLD_INDICATOR_PARAMS["C"]].mean())
        wsum = sum(MINERAL_PRIORITY.values())
        m_presence = sum(float(r[n]) * w for n, w in MINERAL_PRIORITY.items()) / wsum
        m_diversity = sum(float(r[n]) >= 0.5 for n in MINERAL_PRIORITY) / 6.0
        s["M"] = float(pd.Series([m_presence, m_diversity,
                                  r["amorphous_fraction"], r["hydration_depth"]]).mean())
        s["G"] = float(r[OLD_INDICATOR_PARAMS["G"]].mean())
        s["T"] = float(r[OLD_INDICATOR_PARAMS["T"]].mean())
        s["isvm"] = sum(s[k] * WEIGHTS[k] for k in ["C", "M", "G", "T"])
        rows.append(s)
    return pd.DataFrame(rows)


def old_isvm_batch(sc):
    """向量化计算旧 ISVM（蒙特卡洛用）。返回按原行序的 isvm 数组。"""
    r = sc
    c = r[OLD_INDICATOR_PARAMS["C"]].mean(axis=1)
    wsum = sum(MINERAL_PRIORITY.values())
    mp = sum(r[n] * w for n, w in MINERAL_PRIORITY.items()) / wsum
    md = (r[list(MINERAL_PRIORITY)] >= 0.5).sum(axis=1) / 6.0
    m = pd.concat([mp, md, r["amorphous_fraction"], r["hydration_depth"]], axis=1).mean(axis=1)
    g = r[OLD_INDICATOR_PARAMS["G"]].mean(axis=1)
    t = r[OLD_INDICATOR_PARAMS["T"]].mean(axis=1)
    return WEIGHTS["C"] * c + WEIGHTS["M"] * m + WEIGHTS["G"] * g + WEIGHTS["T"] * t


def before_results(name, full_file, gate_radius=5.0):
    df_full = pd.read_csv(os.path.join(DATA_DIR, "full", full_file))
    df_new = load_sites(name)
    sc = score_old(df_full)
    ind = old_indicators(sc)
    eng = evaluate_sites(df_new, weight_scheme="priority", gate_radius=gate_radius, gate_nav=2.7)
    res = ind.merge(eng[["site", "safety_score", "sampleability_score",
                         "deliverability_score"]], on="site")
    res["combined"] = res["isvm"] + res["safety_score"] + res["sampleability_score"] \
        + res["deliverability_score"]
    return res.sort_values("combined", ascending=False).reset_index(drop=True)


def before_montecarlo(full_file, n=2000, pct=0.10, seed=42):
    df = pd.read_csv(os.path.join(DATA_DIR, "full", full_file))
    rng = np.random.default_rng(seed)
    sites = list(df["site"])
    rank_counter = {s: np.zeros(len(sites)) for s in sites}
    old_cont = [n for n in OLD_PARAMS if PARAMETERS[n]["kind"] == "continuous"]
    for _ in range(n):
        d = df.copy()
        for name in old_cont:
            lo, hi = d[name].min(), d[name].max()
            scale = (hi - lo) * pct if hi > lo else max(abs(hi) * pct, 1e-6)
            d[name] = (d[name] + rng.uniform(-scale, scale, size=len(d))).clip(lower=0)
        sc = score_old(d)
        isvm = old_isvm_batch(sc).to_numpy()
        order = [sites[i] for i in np.argsort(-isvm)]
        for rank, s in enumerate(order):
            rank_counter[s][rank] += 1
    prob = pd.DataFrame({s: rank_counter[s] / n for s in sites})
    prob.insert(0, "rank", np.arange(1, len(sites) + 1))
    return prob


def before_sensitivity(full_file, pct=0.10):
    df = pd.read_csv(os.path.join(DATA_DIR, "full", full_file))
    base_sc = score_old(df)
    base = old_isvm_batch(base_sc)
    old_cont = [n for n in OLD_PARAMS if PARAMETERS[n]["kind"] == "continuous"]
    recs = []
    for name in old_cont:
        for sign in (+1, -1):
            d = df.copy()
            d[name] = d[name] * (1.0 + sign * pct)
            sc = score_old(d)
            delta = (old_isvm_batch(sc) - base).abs().max()
            recs.append({"parameter": name, "sign": sign, "max_isvm_delta": float(delta)})
    return pd.DataFrame(recs)


def save(df, name):
    df.to_csv(os.path.join(BEFORE_DIR, name), index=False, encoding="utf-8-sig")


def fig_compare_bennu(before, after):
    b = before.set_index("site")["isvm"]
    a = after.set_index("site")["isvm"]
    sites = list(b.index)
    x = np.arange(len(sites))
    fig, ax = plt.subplots(figsize=(8, 4.2))
    ax.bar(x - 0.2, [b[s] for s in sites], 0.4, label="精简前（21 参数）", color="#8DA0CB")
    ax.bar(x + 0.2, [a[s] for s in sites], 0.4, label="精简后（10 参数）", color="#FC8D62")
    for i, s in enumerate(sites):
        ax.text(i - 0.2, b[s] + 0.01, f"{b[s]:.3f}", ha="center", fontsize=8)
        ax.text(i + 0.2, a[s] + 0.01, f"{a[s]:.3f}", ha="center", fontsize=8)
    ax.set_xticks(x)
    ax.set_xticklabels(sites)
    ax.set_ylim(0, 0.9)
    ax.set_ylabel("ISVM")
    ax.set_title("实验1：贝努四区 ISVM 精简前后对比（排序不变）")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG_DIR, "fig_compare_bennu.png"), bbox_inches="tight",
                facecolor="white")
    plt.close(fig)


def fig_compare_ryugu(before, after):
    b = before.set_index("site")["combined"]
    a = after.set_index("site")["combined"]
    sites = list(b.index)
    x = np.arange(len(sites))
    fig, ax = plt.subplots(figsize=(9, 4.4))
    ax.bar(x - 0.2, [b[s] for s in sites], 0.4, label="精简前（21 参数）", color="#8DA0CB")
    ax.bar(x + 0.2, [a[s] for s in sites], 0.4, label="精简后（10 参数）", color="#FC8D62")
    top3 = list(a.sort_values(ascending=False).head(3).index)
    for i, s in enumerate(sites):
        if s in top3:
            ax.text(i + 0.2, a[s] + 0.02, "TOP3", ha="center", fontsize=8,
                    color="#D62728", fontweight="bold")
        ax.text(i - 0.2, b[s] + 0.02, f"{b[s]:.3f}", ha="center", fontsize=8)
        ax.text(i + 0.2, a[s] - 0.09, f"{a[s]:.3f}", ha="center", fontsize=8, color="#7F2D0F")
    ax.set_xticks(x)
    ax.set_xticklabels(sites)
    ax.set_ylim(0, 2.4)
    ax.set_ylabel("综合分（科学+安全+可采样+可达）")
    ax.set_title("实验2：龙宫七区综合分精简前后对比（前三 L08/L07/M04 不变）")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG_DIR, "fig_compare_ryugu.png"), bbox_inches="tight",
                facecolor="white")
    plt.close(fig)


def fig_compare_params():
    before = {"C": 5, "M": 8, "G": 6, "T": 2}
    after = {k: len(v) for k, v in INDICATOR_PARAMETERS.items()}
    labels = ["化学组成 C", "矿物学 M", "地质新鲜度 G", "温度 T"]
    x = np.arange(4)
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.bar(x - 0.2, [before[k] for k in ["C", "M", "G", "T"]], 0.4,
           label="精简前", color="#8DA0CB")
    ax.bar(x + 0.2, [after[k] for k in ["C", "M", "G", "T"]], 0.4,
           label="精简后", color="#FC8D62")
    for i, k in enumerate(["C", "M", "G", "T"]):
        ax.text(i - 0.2, before[k] + 0.1, str(before[k]), ha="center")
        ax.text(i + 0.2, after[k] + 0.1, str(after[k]), ha="center")
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_ylabel("参数数量")
    ax.set_title("科学参数精简：21 → 10（按指标分解）")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG_DIR, "fig_compare_params.png"), bbox_inches="tight",
                facecolor="white")
    plt.close(fig)


def fig_compare_mc(pb, pa):
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
    for ax, prob, title in [(axes[0], pb, "精简前（21 参数）"),
                            (axes[1], pa, "精简后（10 参数）")]:
        sites = list(prob.columns[1:])
        data = prob[sites].values.T
        im = ax.imshow(data, cmap="YlGnBu", vmin=0, vmax=1)
        ax.set_xticks(range(len(sites)))
        ax.set_xticklabels([f"第{i}名" for i in range(1, len(sites) + 1)], fontsize=8)
        ax.set_yticks(range(len(sites)))
        ax.set_yticklabels(sites, fontsize=9)
        for i in range(data.shape[0]):
            for j in range(data.shape[1]):
                v = data[i, j]
                ax.text(j, i, f"{v:.2f}", ha="center", va="center", fontsize=8,
                        color="black" if v < 0.75 else "white")
        ax.set_title(title)
    fig.suptitle("实验5：蒙特卡洛排名概率对比（Nightingale 均 100% 第一）")
    fig.tight_layout()
    fig.savefig(os.path.join(FIG_DIR, "fig_compare_mc.png"), bbox_inches="tight",
                facecolor="white")
    plt.close(fig)


def fig_compare_sensitivity(sb, sa):
    b = sb.groupby("parameter")["max_isvm_delta"].max().sort_values(ascending=False).head(5)
    a = sa.groupby("parameter")["max_isvm_delta"].max().sort_values(ascending=False).head(5)
    fig, axes = plt.subplots(1, 2, figsize=(11, 4))
    for ax, data, title in [(axes[0], b, "精简前 Top5（21 参数）"),
                            (axes[1], a, "精简后 Top5（10 参数）")]:
        vals = data.sort_values()
        ax.barh(vals.index, vals.values, color="#66A61E")
        for yi, v in enumerate(vals.values):
            ax.text(v + 0.0005, yi, f"{v:.4f}", va="center", fontsize=8)
        ax.set_title(title)
        ax.set_xlabel("±10% 扰动最大 |ΔISVM|")
    fig.suptitle("实验4：敏感性 Top 参数对比")
    fig.tight_layout()
    fig.savefig(os.path.join(FIG_DIR, "fig_compare_sensitivity.png"), bbox_inches="tight",
                facecolor="white")
    plt.close(fig)


def main():
    # 复现精简前结果
    b_bennu = before_results("bennu_sites.csv", "bennu_sites_full.csv")
    b_ryugu = before_results("ryugu_sites.csv", "ryugu_sites_full.csv", gate_radius=3.0)
    b_sub = before_results("ryugu_subsites.csv", "ryugu_subsites_full.csv", gate_radius=3.0)
    save(b_bennu, "bennu_before_reduction.csv")
    save(b_ryugu, "ryugu_before_reduction.csv")
    save(b_sub, "ryugu_subsites_before_reduction.csv")

    pb = before_montecarlo("bennu_sites_full.csv", n=2000, seed=42)
    save(pb, "montecarlo_before_reduction.csv")
    sb = before_sensitivity("bennu_sites_full.csv")
    save(sb, "sensitivity_before_reduction.csv")

    # 精简后结果（当前 outputs/）
    a_bennu = evaluate_sites(load_sites("bennu_sites.csv"), weight_scheme="priority")
    a_ryugu = evaluate_sites(load_sites("ryugu_sites.csv"), weight_scheme="priority",
                             gate_radius=3.0, gate_nav=2.7)
    from src.analysis import monte_carlo_ranking, parameter_sensitivity
    pa, _ = monte_carlo_ranking(load_sites("bennu_sites.csv"), n=2000, pct=0.10, seed=42)
    sa = parameter_sensitivity(load_sites("bennu_sites.csv"), pct=0.10)

    # 对比图
    fig_compare_bennu(b_bennu, a_bennu)
    fig_compare_ryugu(b_ryugu, a_ryugu)
    fig_compare_params()
    fig_compare_mc(pb, pa)
    fig_compare_sensitivity(sb, sa)

    print("精简前结果已保存到:", BEFORE_DIR)
    print("对比图已保存到:", FIG_DIR)
    print("\n贝努 ISVM 精简前后：")
    for s in b_bennu.set_index("site").index:
        print(f"  {s}: {b_bennu.set_index('site').loc[s,'isvm']:.4f} -> "
              f"{a_bennu.set_index('site').loc[s,'isvm']:.4f}")
    print("\n龙宫综合分精简前后 Top3：")
    for s in ["L08", "L07", "M04"]:
        print(f"  {s}: {b_ryugu.set_index('site').loc[s,'combined']:.4f} -> "
              f"{a_ryugu.set_index('site').loc[s,'combined']:.4f}")


if __name__ == "__main__":
    main()
