# -*- coding: utf-8 -*-
"""完整复现第一次实验（精简前：21 个科学参数 + 原始聚合逻辑）。

输出：
- outputs/first_experiment/*.csv   第一次实验全部结果
- outputs/first_experiment/figures/  第一次实验全部图片（fig1–fig10）

输入：data/full/（归档的完整参数数据集）。
"""

import os

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from common import ROOT, DATA_DIR, OUTPUT_DIR, load_sites
from exp8_before_after_compare import (score_old, old_indicators, old_isvm_batch,
                                       before_montecarlo, before_sensitivity,
                                       OLD_PARAMS, OLD_INDICATOR_PARAMS, MINERAL_PRIORITY)
from src.pipeline import engineering_scores
from src.config import PARAMETERS

plt.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei"]
plt.rcParams["axes.unicode_minus"] = False
plt.rcParams["figure.dpi"] = 150

FIRST = os.path.join(OUTPUT_DIR, "first_experiment")
FIGS = os.path.join(FIRST, "figures")
os.makedirs(FIGS, exist_ok=True)

OLD_CONT = [n for n in OLD_PARAMS if PARAMETERS[n]["kind"] == "continuous"]
OLD_CAT = [n for n in OLD_PARAMS if PARAMETERS[n]["kind"] != "continuous"]
WEIGHTS = {"C": 0.35, "M": 0.30, "G": 0.20, "T": 0.15}

SHORT_LABELS = {
    "organic_abundance": "有机物丰度", "volatile_abundance": "挥发物丰度",
    "organic_silicate_ratio": "有机/硅酸盐比", "ch2_ch3": "CH2/CH3",
    "carbon_content_pct": "碳含量", "mineral_phyllo": "层状硅酸盐",
    "mineral_carbonate": "碳酸盐", "mineral_sulfate": "硫酸盐",
    "mineral_oxide": "氧化物", "mineral_silicate": "硅酸盐",
    "mineral_amorphous": "非晶质检出", "amorphous_fraction": "非晶质比例",
    "hydration_depth": "水合带深", "space_weathering_proxy": "空间风化代理",
    "crater_freshness": "陨石坑新鲜度", "activity_level": "活动证据",
    "psfd_class": "PSFD", "brittle_deformation": "脆性构造",
    "geologic_diversity": "地质多样性", "temp_rel_global": "温度等级",
    "tmax_k": "Tmax",
}
IND_LABEL = {"C": "化学组成 C", "M": "矿物学 M", "G": "地质新鲜度 G", "T": "温度 T"}
COLOR_SCI = {"C": "#2C7FB8", "M": "#7FBF7B", "G": "#F1A340", "T": "#D7191C"}
HIGHLIGHT = "#D62728"
BASE_COLOR = "#4C72B0"


def full(name):
    return pd.read_csv(os.path.join(DATA_DIR, "full", name))


def with_engineering(ind_df, current_file, gate_radius=5.0):
    eng = engineering_scores(load_sites(current_file), gate_radius=gate_radius, gate_nav=2.7)
    res = ind_df.merge(eng, on="site")
    res["combined"] = res["isvm"] + res["safety_score"] + res["sampleability_score"] \
        + res["deliverability_score"]
    return res


def save(df, name):
    df.to_csv(os.path.join(FIRST, name), index=False, encoding="utf-8-sig")


def fig(fig, name):
    fig.savefig(os.path.join(FIGS, name), bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print("saved figure:", name)


# ---------------- 实验 1：贝努 ----------------
def exp1():
    sc = score_old(full("bennu_sites_full.csv"))
    ind = old_indicators(sc).set_index("site")
    ind["isvm"] = sum(ind[k] * WEIGHTS[k] for k in ["C", "M", "G", "T"])
    res = ind.reset_index().merge(
        engineering_scores(load_sites("bennu_sites.csv")), on="site")
    res["weight_scheme"] = "priority"
    res = res.sort_values("isvm", ascending=False).reset_index(drop=True)
    save(res, "exp1_bennu_priority.csv")
    save(res[["site", "isvm"]].assign(
        rank=res["isvm"].rank(ascending=False, method="min").astype(int)),
        "exp1_bennu_ranking.csv")

    # fig1
    tmp = res.sort_values("isvm", ascending=True)
    fig_, ax = plt.subplots(figsize=(7, 3.6))
    colors = [HIGHLIGHT if s == "Nightingale" else BASE_COLOR for s in tmp["site"]]
    bars = ax.barh(tmp["site"], tmp["isvm"], color=colors, edgecolor="white")
    for b, v in zip(bars, tmp["isvm"]):
        ax.text(v + 0.008, b.get_y() + b.get_height() / 2, f"{v:.3f}",
                va="center", fontsize=9)
    ax.set_xlim(0, 0.85)
    ax.set_xlabel("ISVM 科学价值分")
    ax.set_title("第一次实验：贝努四候选区 ISVM 排序（21 参数）")
    fig(fig_, "fig1_bennu_ranking.png")

    # fig2
    tmp = res.sort_values("isvm", ascending=False)
    x = np.arange(len(tmp))
    fig_, ax = plt.subplots(figsize=(8, 4.2))
    for i, k in enumerate(["C", "M", "G", "T"]):
        ax.bar(x + (i - 1.5) * 0.2, tmp[k], 0.2, label=IND_LABEL[k], color=COLOR_SCI[k])
    ax.set_xticks(x)
    ax.set_xticklabels(tmp["site"])
    ax.set_ylim(0, 1)
    ax.set_ylabel("指标得分")
    ax.set_title("第一次实验：贝努四候选区 C/M/G/T 指标分解")
    ax.legend(ncol=4, frameon=False)
    fig(fig_, "fig2_bennu_indicators.png")

    # fig3 雷达
    cats = ["C", "M", "G", "T"]
    ang = np.linspace(0, 2 * np.pi, len(cats), endpoint=False).tolist() + [0]
    fig_, ax = plt.subplots(figsize=(5.6, 5.6), subplot_kw=dict(polar=True))
    for site in res["site"]:
        vals = res.loc[res["site"] == site, cats].iloc[0].tolist() + [res.loc[
            res["site"] == site, cats].iloc[0].tolist()[0]]
        ax.plot(ang, vals, label=site, lw=2)
        ax.fill(ang, vals, alpha=0.08)
    ax.set_xticks(ang[:-1])
    ax.set_xticklabels([IND_LABEL[c] for c in cats])
    ax.set_ylim(0, 1)
    ax.set_title("第一次实验：贝努四候选区指标雷达图", pad=18)
    ax.legend(loc="upper right", bbox_to_anchor=(1.35, 1.12), frameon=False)
    fig(fig_, "fig3_bennu_radar.png")
    return res


# ---------------- 实验 2：龙宫 ----------------
def exp2():
    for name, full_file, current, gate in [
            ("exp2_ryugu_seven.csv", "ryugu_sites_full.csv", "ryugu_sites.csv", 3.0),
            ("exp2_ryugu_subsites.csv", "ryugu_subsites_full.csv", "ryugu_subsites.csv", 3.0)]:
        sc = score_old(full(full_file))
        ind = old_indicators(sc).set_index("site")
        ind["isvm"] = sum(ind[k] * WEIGHTS[k] for k in ["C", "M", "G", "T"])
        res = with_engineering(ind.reset_index(), current, gate_radius=gate)
        res = res.sort_values("combined", ascending=False).reset_index(drop=True)
        save(res, name)
    return (pd.read_csv(os.path.join(FIRST, "exp2_ryugu_seven.csv")),
            pd.read_csv(os.path.join(FIRST, "exp2_ryugu_subsites.csv")))


def fig4_5(seven, sub):
    # fig4
    tmp = seven.sort_values("combined", ascending=False)
    x = np.arange(len(tmp))
    fig_, ax = plt.subplots(figsize=(9, 4.4))
    ax.bar(x - 0.19, tmp["isvm"], 0.38, label="科学价值 ISVM", color="#8DA0CB")
    top3 = set(tmp.head(3)["site"])
    colors = [HIGHLIGHT if s in top3 else "#66C2A5" for s in tmp["site"]]
    ax.bar(x + 0.19, tmp["combined"], 0.38, label="综合分（含工程层）", color=colors)
    for i, (a, b) in enumerate(zip(tmp["isvm"], tmp["combined"])):
        ax.text(i - 0.19, a + 0.02, f"{a:.2f}", ha="center", fontsize=8)
        ax.text(i + 0.19, b + 0.02, f"{b:.2f}", ha="center", fontsize=8)
    ax.set_xticks(x)
    ax.set_xticklabels(tmp["site"])
    ax.set_ylim(0, 2.3)
    ax.set_ylabel("得分")
    ax.set_title("第一次实验：龙宫七候选区 科学 ISVM vs 综合分\n（红色=文献最终 TD1 候选 L08/L07/M04）")
    ax.legend(frameon=False)
    fig(fig_, "fig4_ryugu_combined.png")

    # fig5
    tmp = sub.set_index("site")
    fig_, ax = plt.subplots(figsize=(6.5, 4))
    bottom = np.zeros(len(tmp))
    parts = [("isvm", "科学价值"), ("safety_score", "安全"),
             ("sampleability_score", "可采样"), ("deliverability_score", "可达")]
    colors = ["#2C7FB8", "#7FBF7B", "#F1A340", "#D7191C"]
    for i, (col, label) in enumerate(parts):
        ax.bar(tmp.index, tmp[col], bottom=bottom, label=label, color=colors[i])
        bottom += tmp[col].values
    for i, t in enumerate(bottom):
        ax.text(i, t + 0.03, f"{t:.2f}", ha="center", fontsize=10, fontweight="bold")
    ax.set_ylim(0, 2.7)
    ax.set_ylabel("综合分构成")
    ax.set_title("第一次实验：龙宫 L08 子区综合分构成")
    ax.legend(frameon=False)
    fig(fig_, "fig5_ryugu_subsites.png")


# ---------------- 实验 3：权重对比 ----------------
def exp3(bennu_res):
    sc = score_old(full("bennu_sites_full.csv"))
    ind = old_indicators(sc).set_index("site")
    schemes = {"priority": WEIGHTS, "equal": {k: 0.25 for k in ["C", "M", "G", "T"]},
               "entropy": None}
    from src.weights import entropy_weights
    schemes["entropy"] = entropy_weights(ind)
    out = pd.DataFrame({"site": ind.index})
    for name, w in schemes.items():
        ind2 = ind.copy()
        ind2["isvm"] = sum(ind2[k] * w[k] for k in ["C", "M", "G", "T"])
        res = ind2.reset_index().sort_values("isvm", ascending=False)
        save(res, f"exp3_bennu_{name}.csv")
        out[f"isvm_{name}"] = ind2["isvm"].values
    save(out, "exp3_weight_comparison.csv")
    print("权重方案:", {k: {kk: round(vv, 4) for kk, vv in v.items()} for k, v in schemes.items()})

    # fig6
    sites = list(ind.index)
    x = np.arange(len(sites))
    fig_, ax = plt.subplots(figsize=(8, 4.2))
    colors = ["#4C72B0", "#55A868", "#C44E52"]
    for i, (name, w) in enumerate(schemes.items()):
        vals = [sum(ind.loc[s, k] * w[k] for k in ["C", "M", "G", "T"]) for s in sites]
        ax.bar(x + (i - 1) * 0.26, vals, 0.26, label=name, color=colors[i])
    ax.set_xticks(x)
    ax.set_xticklabels(sites)
    ax.set_ylim(0, 0.85)
    ax.set_ylabel("ISVM")
    ax.set_title("第一次实验：三种权重方案下贝努候选区 ISVM")
    ax.legend(frameon=False)
    fig(fig_, "fig6_weight_comparison.png")


# ---------------- 实验 4：敏感性 ----------------
def exp4():
    df = full("bennu_sites_full.csv")
    base_sc = score_old(df)
    base = old_isvm_batch(base_sc)
    base_order = list(df["site"][np.argsort(-base.to_numpy())])
    recs = []
    for name in OLD_CONT:
        for sign in (+1, -1):
            d = df.copy()
            d[name] = d[name] * (1.0 + sign * 0.10)
            sc = score_old(d)
            isvm = old_isvm_batch(sc)
            order = list(df["site"][np.argsort(-isvm.to_numpy())])
            changes = sum(a != b for a, b in zip(base_order, order))
            delta = (isvm - base).abs().max()
            recs.append({"parameter": name, "sign": sign, "rank_changes": changes,
                         "max_isvm_delta": float(delta), "order": "/".join(order)})
    save(pd.DataFrame(recs), "exp4_parameter_sensitivity.csv")

    wrecs = []
    for c in ["C", "M", "G", "T"]:
        for sign in (+1, -1):
            w = dict(WEIGHTS)
            w[c] *= 1.2 if sign > 0 else 0.8
            s = sum(w.values())
            w = {k: v / s for k, v in w.items()}
            sc = score_old(df)
            ind = old_indicators(sc).set_index("site")
            isvm = sum(ind[k] * w[k] for k in ["C", "M", "G", "T"])
            order = list(isvm.sort_values(ascending=False).index)
            changes = sum(a != b for a, b in zip(base_order, order))
            wrecs.append({"indicator": c, "sign": sign,
                          "weight": round(w[c], 3), "rank_changes": changes})
    save(pd.DataFrame(wrecs), "exp4_weight_sensitivity.csv")

    # fig8
    impact = (pd.DataFrame(recs).groupby("parameter")["max_isvm_delta"].max()
              .sort_values(ascending=True).tail(10))
    labels = [SHORT_LABELS[p] for p in impact.index]
    fig_, ax = plt.subplots(figsize=(7.5, 4.2))
    bars = ax.barh(labels, impact.values, color="#66A61E")
    for b, v in zip(bars, impact.values):
        ax.text(v + 0.0005, b.get_y() + b.get_height() / 2, f"{v:.4f}",
                va="center", fontsize=8)
    ax.set_xlabel("±10% 参数扰动引起的最大 |ΔISVM|")
    ax.set_title("第一次实验：对 ISVM 数值影响最大的参数（贝努）")
    fig(fig_, "fig8_sensitivity.png")


# ---------------- 实验 5：蒙特卡洛 ----------------
def exp5():
    prob = before_montecarlo("bennu_sites_full.csv", n=2000, pct=0.10, seed=42)
    save(prob, "exp5_monte_carlo_rank_probability.csv")
    sites = list(prob.columns[1:])
    data = prob[sites].values.T
    fig_, ax = plt.subplots(figsize=(6.5, 4))
    im = ax.imshow(data, cmap="YlGnBu", vmin=0, vmax=1)
    ax.set_xticks(range(len(sites)))
    ax.set_xticklabels([f"第{i}名" for i in range(1, len(sites) + 1)])
    ax.set_yticks(range(len(sites)))
    ax.set_yticklabels(sites)
    for i in range(data.shape[0]):
        for j in range(data.shape[1]):
            v = data[i, j]
            ax.text(j, i, f"{v:.2f}", ha="center", va="center", fontsize=9,
                    color="black" if v < 0.75 else "white")
    ax.set_title("第一次实验：蒙特卡洛排名概率（2000 次抽样）")
    fig_.colorbar(im, ax=ax, shrink=0.85, label="排名概率")
    fig(fig_, "fig7_monte_carlo.png")


# ---------------- 实验 6：合成压力测试 ----------------
def exp6():
    rng = np.random.default_rng(7)
    template = pd.read_csv(os.path.join(DATA_DIR, "full", "parameter_template_full.csv")).iloc[0]
    rows = []
    for i in range(30):
        r = template.copy()
        r["site"] = f"SYN-{i:02d}"
        r["organic_abundance"] = rng.uniform(0, 1)
        r["volatile_abundance"] = rng.uniform(0, 1)
        r["organic_silicate_ratio"] = rng.uniform(0, 1)
        r["ch2_ch3"] = rng.uniform(0.5, 3.0)
        r["carbon_content_pct"] = rng.uniform(0, 10)
        for m in ["mineral_phyllo", "mineral_carbonate", "mineral_sulfate",
                  "mineral_oxide", "mineral_silicate", "mineral_amorphous"]:
            r[m] = rng.choice([0.0, 0.5, 1.0])
        r["amorphous_fraction"] = rng.uniform(0, 1)
        r["hydration_depth"] = rng.uniform(0, 1)
        r["space_weathering_proxy"] = rng.uniform(0, 1)
        r["crater_freshness"] = rng.choice(["recent", "paleo", "none"])
        r["activity_level"] = rng.choice(["active", "recent", "paleo", "none"])
        r["psfd_class"] = rng.choice(["fine", "mix", "coarse"])
        r["brittle_deformation"] = rng.choice([0.0, 1.0])
        r["geologic_diversity"] = rng.uniform(0, 1)
        r["temp_rel_global"] = rng.choice(["colder", "average", "hotter"])
        r["tmax_k"] = rng.uniform(280, 350)
        r["slope_deg"] = rng.uniform(0, 20)
        r["boulder_coverage"] = rng.uniform(0, 0.30)
        r["roughness"] = rng.uniform(0, 1)
        r["frac_le2cm"] = rng.uniform(0, 1)
        r["min_particle_cm"] = rng.uniform(0.1, 3)
        r["max_particle_cm"] = rng.uniform(1, 10)
        r["tagsam_tilt_deg"] = rng.uniform(0, 20)
        r["target_radius_m"] = rng.uniform(2, 25)
        r["nav_error_m"] = rng.uniform(0.5, 5)
        rows.append(r)
    syn = pd.DataFrame(rows)
    sc = score_old(syn)
    ind = old_indicators(sc).set_index("site")
    ind["isvm"] = sum(ind[k] * WEIGHTS[k] for k in ["C", "M", "G", "T"])
    res = ind.reset_index().sort_values("isvm", ascending=False).reset_index(drop=True)
    assert res["isvm"].between(0, 1).all() and res.notna().all().all()
    save(res, "exp6_synthetic_ranking.csv")
    print(f"实验6：合成数据 30 站点，ISVM 均在 [0,1]、无 NaN，通过。")


# ---------------- 图 9/10：21 参数热图 ----------------
def fig9_10():
    for df, fig_name, title in [(full("bennu_sites_full.csv"), "fig9_raw_bennu_heatmap.png",
                                 "第一次实验原始参数热图：贝努四候选区（21 参数）"),
                                (full("ryugu_sites_full.csv"), "fig10_raw_ryugu_heatmap.png",
                                 "第一次实验原始参数热图：龙宫七候选区（21 参数）")]:
        sc = score_old(df)
        labels = [SHORT_LABELS[p] for p in OLD_PARAMS]
        data = sc[OLD_PARAMS].values.T
        fig_, ax = plt.subplots(figsize=(max(5, len(df) * 1.4 + 2), 8.5))
        im = ax.imshow(data, cmap="YlGnBu", vmin=0, vmax=1, aspect="auto")
        ax.set_yticks(range(len(labels)))
        ax.set_yticklabels(labels)
        ax.set_xticks(range(len(df)))
        ax.set_xticklabels(sc.index)
        for i in range(data.shape[0]):
            for j in range(data.shape[1]):
                v = data[i, j]
                ax.text(j, i, f"{v:.2f}", ha="center", va="center", fontsize=7,
                        color="black" if v < 0.7 else "white")
        ax.set_title(title)
        fig_.colorbar(im, ax=ax, shrink=0.6, label="评分（0–1）")
        fig(fig_, fig_name)


def main():
    print("=== 复现第一次实验（21 参数）===")
    res1 = exp1()
    seven, sub = exp2()
    fig4_5(seven, sub)
    exp3(res1)
    exp4()
    exp5()
    exp6()
    fig9_10()
    print("\n第一次实验完整数据已保存到:", FIRST)
    print("图片已保存到:", FIGS)


if __name__ == "__main__":
    main()
