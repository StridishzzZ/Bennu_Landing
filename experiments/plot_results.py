# -*- coding: utf-8 -*-
"""把原始数据与实验结果渲染为图片（输出到 outputs/figures/）。

运行：python experiments/plot_results.py
依赖：matplotlib（建议使用工作目录配置的 conda 虚拟环境）
"""

import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from common import ROOT, load_sites
from src.pipeline import evaluate_sites
from src.analysis import parameter_sensitivity, monte_carlo_ranking

plt.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei"]
plt.rcParams["axes.unicode_minus"] = False
plt.rcParams["figure.dpi"] = 150

FIG_DIR = os.path.join(ROOT, "outputs", "figures")
os.makedirs(FIG_DIR, exist_ok=True)

COLOR_SCI = {"C": "#2C7FB8", "M": "#7FBF7B", "G": "#F1A340", "T": "#D7191C"}
HIGHLIGHT = "#D62728"
BASE_COLOR = "#4C72B0"
ENG_COLORS = ["#2C7FB8", "#7FBF7B", "#F1A340", "#D7191C"]

SHORT_LABELS = {
    "organic_abundance": "有机物丰度", "volatile_abundance": "挥发物丰度",
    "organic_silicate_ratio": "有机/硅酸盐比", "ch2_ch3": "CH2/CH3",
    "carbon_content_pct": "碳含量", "mineral_classes_detected": "矿物类别数",
    "mineral_phyllo": "层状硅酸盐",
    "mineral_carbonate": "碳酸盐", "mineral_sulfate": "硫酸盐",
    "mineral_oxide": "氧化物", "mineral_silicate": "硅酸盐",
    "mineral_amorphous": "非晶质检出", "amorphous_fraction": "非晶质比例",
    "hydration_depth": "水合带深", "space_weathering_proxy": "空间风化代理",
    "crater_freshness": "陨石坑新鲜度", "activity_level": "活动证据",
    "psfd_class": "PSFD", "brittle_deformation": "脆性构造",
    "geologic_diversity": "地质多样性", "temp_rel_global": "温度等级",
    "tmax_k": "Tmax",
}


def _save(fig, name):
    path = os.path.join(FIG_DIR, name)
    fig.savefig(path, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print("saved:", name)


def fig1_bennu_ranking(bennu):
    res = evaluate_sites(bennu, weight_scheme="priority")
    res = res.sort_values("isvm", ascending=True)
    fig, ax = plt.subplots(figsize=(7, 3.6))
    colors = [HIGHLIGHT if s == "Nightingale" else BASE_COLOR for s in res["site"]]
    bars = ax.barh(res["site"], res["isvm"], color=colors, edgecolor="white")
    for b, v in zip(bars, res["isvm"]):
        ax.text(v + 0.008, b.get_y() + b.get_height() / 2, f"{v:.3f}",
                va="center", fontsize=9)
    ax.set_xlim(0, 0.85)
    ax.set_xlabel("ISVM 科学价值分")
    ax.set_title("实验1：贝努四候选区 ISVM 排序（文献优先级权重）")
    ax.axvline(0.6, color="grey", ls="--", lw=0.8, alpha=0.6)
    ax.text(0.605, -0.55, "参考阈值 0.6", fontsize=8, color="grey")
    fig.tight_layout()
    _save(fig, "fig1_bennu_ranking.png")


def fig2_bennu_indicators(bennu):
    res = evaluate_sites(bennu, weight_scheme="priority")
    res = res.sort_values("isvm", ascending=False)
    x = np.arange(len(res))
    width = 0.2
    fig, ax = plt.subplots(figsize=(8, 4.2))
    for i, k in enumerate(["C", "M", "G", "T"]):
        ax.bar(x + (i - 1.5) * width, res[k], width, label=IND_LABEL[k],
               color=COLOR_SCI[k])
    ax.set_xticks(x)
    ax.set_xticklabels(res["site"])
    ax.set_ylim(0, 1)
    ax.set_ylabel("指标得分")
    ax.set_title("实验1：贝努四候选区 C/M/G/T 指标分解")
    ax.legend(ncol=4, frameon=False)
    fig.tight_layout()
    _save(fig, "fig2_bennu_indicators.png")


def fig3_bennu_radar(bennu):
    res = evaluate_sites(bennu, weight_scheme="priority").set_index("site")
    cats = ["C", "M", "G", "T"]
    ang = np.linspace(0, 2 * np.pi, len(cats), endpoint=False).tolist()
    ang += ang[:1]
    fig, ax = plt.subplots(figsize=(5.6, 5.6), subplot_kw=dict(polar=True))
    for site in res.index:
        vals = res.loc[site, cats].tolist()
        vals += vals[:1]
        ax.plot(ang, vals, label=site, lw=2)
        ax.fill(ang, vals, alpha=0.08)
    ax.set_xticks(ang[:-1])
    ax.set_xticklabels([IND_LABEL[c] for c in cats])
    ax.set_ylim(0, 1)
    ax.set_title("实验1：贝努四候选区指标雷达图", pad=18)
    ax.legend(loc="upper right", bbox_to_anchor=(1.35, 1.12), frameon=False)
    _save(fig, "fig3_bennu_radar.png")


def fig4_ryugu_combined(ryugu):
    res = evaluate_sites(ryugu, weight_scheme="priority", gate_radius=3.0, gate_nav=2.7)
    res = res.sort_values("combined", ascending=False)
    x = np.arange(len(res))
    width = 0.38
    fig, ax = plt.subplots(figsize=(9, 4.4))
    ax.bar(x - width / 2, res["isvm"], width, label="科学价值 ISVM", color="#8DA0CB")
    top3 = set(res.head(3)["site"])
    colors = [HIGHLIGHT if s in top3 else "#66C2A5" for s in res["site"]]
    ax.bar(x + width / 2, res["combined"], width, label="综合分（含工程层）", color=colors)
    for i, (a, b) in enumerate(zip(res["isvm"], res["combined"])):
        ax.text(i - width / 2, a + 0.02, f"{a:.2f}", ha="center", fontsize=8)
        ax.text(i + width / 2, b + 0.02, f"{b:.2f}", ha="center", fontsize=8)
    ax.set_xticks(x)
    ax.set_xticklabels(res["site"])
    ax.set_ylim(0, 2.3)
    ax.set_ylabel("得分")
    ax.set_title("实验2：龙宫七候选区 科学 ISVM vs 综合分\n（红色=文献最终 TD1 候选 L08/L07/M04）")
    ax.legend(frameon=False)
    fig.tight_layout()
    _save(fig, "fig4_ryugu_combined.png")


def fig5_ryugu_subsites(ryugu_sub):
    res = evaluate_sites(ryugu_sub, weight_scheme="priority", gate_radius=3.0, gate_nav=2.7)
    res = res.set_index("site")
    fig, ax = plt.subplots(figsize=(6.5, 4))
    bottom = np.zeros(len(res))
    parts = [("isvm", "科学价值"), ("safety_score", "安全"),
             ("sampleability_score", "可采样"), ("deliverability_score", "可达")]
    for i, (col, label) in enumerate(parts):
        vals = res[col].values
        ax.bar(res.index, vals, bottom=bottom, label=label, color=ENG_COLORS[i])
        bottom += vals
    for i, total in enumerate(bottom):
        ax.text(i, total + 0.03, f"{total:.2f}", ha="center", fontsize=10, fontweight="bold")
    ax.set_ylim(0, 2.7)
    ax.set_ylabel("综合分构成")
    ax.set_title("实验2：龙宫 L08 子区综合分构成\n（L08-E1=实际 TD1 点；L08-B=SCI 后 TD2 点）")
    ax.legend(frameon=False)
    fig.tight_layout()
    _save(fig, "fig5_ryugu_subsites.png")


def fig6_weight_comparison(bennu):
    sites = list(bennu["site"])
    schemes = ["priority", "equal", "entropy"]
    scheme_label = {"priority": "文献优先级", "equal": "等权", "entropy": "熵权"}
    scores = {s: [] for s in sites}
    for sc in schemes:
        res = evaluate_sites(bennu, weight_scheme=sc).set_index("site")
        for s in sites:
            scores[s].append(res.loc[s, "isvm"])
    x = np.arange(len(sites))
    width = 0.26
    fig, ax = plt.subplots(figsize=(8, 4.2))
    colors = ["#4C72B0", "#55A868", "#C44E52"]
    for i, sc in enumerate(schemes):
        vals = [scores[s][i] for s in sites]
        ax.bar(x + (i - 1) * width, vals, width, label=scheme_label[sc], color=colors[i])
    ax.set_xticks(x)
    ax.set_xticklabels(sites)
    ax.set_ylim(0, 0.85)
    ax.set_ylabel("ISVM")
    ax.set_title("实验3：三种权重方案下贝努候选区 ISVM（排序完全一致）")
    ax.legend(frameon=False)
    fig.tight_layout()
    _save(fig, "fig6_weight_comparison.png")


def fig7_monte_carlo(bennu):
    prob, _ = monte_carlo_ranking(bennu, n=2000, pct=0.10, seed=42)
    sites = list(prob.columns[1:])
    data = prob[sites].values.T
    fig, ax = plt.subplots(figsize=(6.5, 4))
    im = ax.imshow(data, cmap="YlGnBu", vmin=0, vmax=1)
    ax.set_xticks(range(len(sites)))
    ax.set_xticklabels([f"第{i}名" for i in range(1, len(sites) + 1)])
    ax.set_yticks(range(len(sites)))
    ax.set_yticklabels(sites)
    for i in range(data.shape[0]):
        for j in range(data.shape[1]):
            v = data[i, j]
            ax.text(j, i, f"{v:.2f}", ha="center", va="center",
                    fontsize=9, color="black" if v < 0.75 else "white")
    ax.set_title("实验5：蒙特卡洛排名概率（2000 次抽样）")
    fig.colorbar(im, ax=ax, shrink=0.85, label="排名概率")
    fig.tight_layout()
    _save(fig, "fig7_monte_carlo.png")


def fig8_sensitivity(bennu):
    sens = parameter_sensitivity(bennu, pct=0.10)
    impact = (sens.groupby("parameter")["max_isvm_delta"].max()
              .sort_values(ascending=True).tail(10))
    labels = [SHORT_LABELS.get(p, p) for p in impact.index]
    fig, ax = plt.subplots(figsize=(7.5, 4.2))
    bars = ax.barh(labels, impact.values, color="#66A61E")
    for b, v in zip(bars, impact.values):
        ax.text(v + 0.0005, b.get_y() + b.get_height() / 2, f"{v:.4f}",
                va="center", fontsize=8)
    ax.set_xlabel("±10% 参数扰动引起的最大 |ΔISVM|")
    ax.set_title("实验4：对 ISVM 数值影响最大的参数（贝努）")
    fig.tight_layout()
    _save(fig, "fig8_sensitivity.png")


def fig9_heatmap_raw(df, name, title):
    from src.scoring import score_all_parameters
    from src.config import SCIENCE_PARAMETERS
    scored = score_all_parameters(df).set_index("site")
    labels = [SHORT_LABELS[p] for p in SCIENCE_PARAMETERS]
    data = scored[SCIENCE_PARAMETERS].values.T
    fig, ax = plt.subplots(figsize=(max(5, len(df) * 1.4 + 2), 8.5))
    im = ax.imshow(data, cmap="YlGnBu", vmin=0, vmax=1, aspect="auto")
    ax.set_yticks(range(len(labels)))
    ax.set_yticklabels(labels)
    ax.set_xticks(range(len(df)))
    ax.set_xticklabels(scored.index)
    for i in range(data.shape[0]):
        for j in range(data.shape[1]):
            v = data[i, j]
            ax.text(j, i, f"{v:.2f}", ha="center", va="center", fontsize=7,
                    color="black" if v < 0.7 else "white")
    ax.set_title(title)
    fig.colorbar(im, ax=ax, shrink=0.6, label="评分（0–1）")
    fig.tight_layout()
    _save(fig, name)


IND_LABEL = {
    "C": "化学组成 C", "M": "矿物学 M",
    "G": "地质新鲜度 G", "T": "温度 T",
}


def main():
    bennu = load_sites("bennu_sites.csv")
    ryugu = load_sites("ryugu_sites.csv")
    ryugu_sub = load_sites("ryugu_subsites.csv")

    fig1_bennu_ranking(bennu)
    fig2_bennu_indicators(bennu)
    fig3_bennu_radar(bennu)
    fig4_ryugu_combined(ryugu)
    fig5_ryugu_subsites(ryugu_sub)
    fig6_weight_comparison(bennu)
    fig7_monte_carlo(bennu)
    fig8_sensitivity(bennu)
    fig9_heatmap_raw(bennu, "fig9_raw_bennu_heatmap.png", "原始参数评分热图：贝努四候选区")
    fig9_heatmap_raw(ryugu, "fig10_raw_ryugu_heatmap.png", "原始参数评分热图：龙宫七候选区")

    print("全部图片已输出到:", FIG_DIR)


if __name__ == "__main__":
    main()
