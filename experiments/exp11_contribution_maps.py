# -*- coding: utf-8 -*-
"""实验 11：着陆点科学性价值 + 各参数贡献度的真实 Bennu 地图可视化。

设定（经用户确认）
----------------
1. 不使用坡度等工程门控，只用现有程序的 10 个科学参数；
2. 贡献度为**加法分解**：``c(p) = s(p) × w(indicator(p)) / n(indicator(p))``，
   10 个 c(p) 逐像元求和恒等于 ISVM；
3. **全幅像元级填充（实装测试）**：候选区用站点级评分（与实验 1 完全一致），
   区外填充合成占位数据——连续参数取"参数域内最差端"的窄区间并叠加固定种子的
   微弱平滑扰动（评分 ∈ [0, FILL_FRAC]），分类参数取评分最低的档位。
   因此填充区评分恒 ≤ 候选区最小值，**不影响候选区结果**，且与候选区差别最大；
   填充区在 `bennu_fill_mask.tiff` 中标记为 0，可随时识别为合成数据；
4. 格网沿用实验 10 的真实 DOM 格网，保证与真实 Bennu 影像逐像元套合。

输出
----
* ``bennu_isvm_landing_value.tiff``      着陆点科学性价值（ISVM），单波段
* ``bennu_param_contribution.tiff``      10 波段，每个科学参数的绝对贡献度
* ``bennu_param_contribution_share.tiff``10 波段，贡献度占比（%，逐像元合计 100）
* ``bennu_fill_mask.tiff``                掩膜：1=候选区实测/文献值，0=合成填充
* ``fig_landing_value_on_bennu.png``     总图：真实 Bennu 影像 + 四站点 ISVM
* ``fig_param_<参数名>.png`` ×10         每个参数的贡献度地图
* ``fig_contribution_bars.png``          四站点贡献度堆叠条形图
* ``fig_dom_global_sites.png``           真实 DOM 全域大图 + 四站点标注（含 ISVM 评分）
* ``outputs/exp11_param_contribution.csv`` 站点 × 参数的贡献度明细

用法
----
    python experiments/exp11_contribution_maps.py
    python experiments/exp11_contribution_maps.py --decimate 5 --weights priority
"""

import argparse
import os
import time

import numpy as np
import pandas as pd

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import Circle  # noqa: E402

from common import DATA_DIR, OUTPUT_DIR  # noqa: E402
from src.config import INDICATOR_PARAMETERS, PARAMETERS, SCIENCE_PARAMETERS  # noqa: E402
from src.pipeline import evaluate_sites  # noqa: E402
from src.raster import BENNU_MEAN_RADIUS_M, GridSpec, read_decimated, write_geotiff  # noqa: E402
from src.scoring import score_all_parameters  # noqa: E402
from src.weights import get_weights  # noqa: E402
from src.mineral_maps import _bilinear_upsample  # noqa: E402  复用平滑场放大工具

PLANETARY_DIR = os.path.join(DATA_DIR, "planetary")
DOM_FILE = "Bennu_DOM_DEG.tif"
SRC_RES_DEG = 0.01145876
SRC_W, SRC_H = 31417, 15709
DEG2M = np.pi * BENNU_MEAN_RADIUS_M / 180.0
NODATA = -9999.0

# 填充区合成数据参数
FILL_FRAC = 0.05     # 连续参数：填充值落在"最差端"域宽的 5% 内 → 评分 ∈ [0, 0.05]
FILL_SEED = 20260920  # 固定随机种子，保证填充场可复现
OVERLAY_ALPHA = 0.75  # 数值层叠加到真实底图上的统一透明度（保证颜色↔数值编码一致）

# 参数贡献度图使用的红蓝渐变色标：蓝=低，白=中，红=高
# 可替换为 "bwr"（更饱和的蓝白红）、"RdBu_r"、"seismic" 等
PARAM_CMAP = "coolwarm"


def _setup_cjk_font():
    from matplotlib import font_manager

    installed = {f.name for f in font_manager.fontManager.ttflist}
    for cand in ["Microsoft YaHei", "SimHei", "Noto Sans CJK SC", "SimSun"]:
        if cand in installed:
            matplotlib.rcParams["font.sans-serif"] = [cand, "DejaVu Sans"]
            matplotlib.rcParams["axes.unicode_minus"] = False
            return cand
    return None


def _to_signed_lon(lon):
    return ((np.asarray(lon, dtype="float64") + 180.0) % 360.0) - 180.0


def effective_weights(weights):
    """把指标权重均分到该指标内的参数上：w(p) = w(指标)/n(指标内参数数)。"""
    eff = {}
    for ind, params in INDICATOR_PARAMETERS.items():
        for p in params:
            eff[p] = float(weights[ind]) / len(params)
    return eff


def _smooth_noise_field(rng, shape, coarse=12, octaves=2):
    """可复现的多尺度平滑噪声场，值域归一到 [0, 1]，float32。"""
    h, w = shape
    acc = np.zeros((h, w), dtype="float64")
    amp, total = 1.0, 0.0
    for o in range(octaves):
        ch = max(3, int(h / (coarse * (2 ** o))))
        cw = max(3, int(w / (coarse * (2 ** o))))
        acc += amp * _bilinear_upsample(rng.random((ch, cw)), h, w)
        total += amp
        amp *= 0.5
    acc /= max(total, 1e-9)
    lo, hi = float(acc.min()), float(acc.max())
    return ((acc - lo) / max(hi - lo, 1e-12)).astype("float32")


def _fill_score_field(param, shape, rng, frac=FILL_FRAC):
    """候选区以外的填充评分场（float32）。

    * 连续型：取值限制在参数域"最差端"的 ``frac`` 区间内，评分 ∈ [0, frac]；
    * 分类型：取评分最低的档位（常数，因为分类评分是离散的）。

    两种情形都保证填充评分 ≤ 任何候选区评分，从而实现"不影响结果 + 差别最大"。
    """
    d = PARAMETERS[param]
    if d["kind"] != "continuous":
        return np.full(shape, float(min(d["categories"].values())), dtype="float32")
    noise = _smooth_noise_field(rng, shape)
    span = float(d["hi"] - d["lo"])
    if d["direction"] == "higher_better":
        value = float(d["lo"]) + frac * span * noise          # 越接近 lo 越差
        return ((value - d["lo"]) / span).astype("float32")
    value = float(d["hi"]) - frac * span * noise              # 越接近 hi 越差
    return (1.0 - (value - d["lo"]) / span).astype("float32")


def _site_slice(grid, lat, lon, radius_m, margin=1.6):
    """返回覆盖站点 (±margin·radius) 的切片与圆盘掩膜。"""
    rows, cols = grid.lat_centers, grid.lon_centers
    ci = int(np.argmin(np.abs(rows - lat)))
    cj = int(np.argmin(np.abs(cols - lon)))
    ri = int(np.ceil(margin * radius_m / (grid.res_deg * DEG2M)))
    rj = int(np.ceil(margin * radius_m / (grid.res_deg * DEG2M
                                         * max(np.cos(np.radians(lat)), 1e-6))))
    sl = (slice(max(0, ci - ri), min(grid.height, ci + ri + 1)),
          slice(max(0, cj - rj), min(grid.width, cj + rj + 1)))
    sub_lat = rows[sl[0]][:, None]
    sub_lon = cols[sl[1]][None, :]
    dlat_m = (sub_lat - lat) * DEG2M
    dlon_m = ((sub_lon - lon + 180.0) % 360.0 - 180.0) * DEG2M \
        * max(np.cos(np.radians(lat)), 1e-6)
    disk = np.hypot(dlat_m, dlon_m) <= radius_m
    return sl, disk


def _zoom_extent(grid, lat, lon, zoom_m):
    """返回以站点为中心的局部米制窗口 (slice, extent)，保持正方形视场。"""
    half_y = zoom_m / 2.0 / (grid.res_deg * DEG2M)
    half_x = zoom_m / 2.0 / (grid.res_deg * DEG2M * max(np.cos(np.radians(lat)), 1e-6))
    rows, cols = grid.lat_centers, grid.lon_centers
    ci = int(np.argmin(np.abs(rows - lat)))
    cj = int(np.argmin(np.abs(cols - lon)))
    sl = (slice(max(0, ci - int(half_y)), min(grid.height, ci + int(half_y) + 1)),
          slice(max(0, cj - int(half_x)), min(grid.width, cj + int(half_x) + 1)))
    width_m = (sl[1].stop - sl[1].start) * grid.res_deg * DEG2M \
        * max(np.cos(np.radians(lat)), 1e-6)
    height_m = (sl[0].stop - sl[0].start) * grid.res_deg * DEG2M
    return sl, [-width_m / 2, width_m / 2, -height_m / 2, height_m / 2]


def main():
    ap = argparse.ArgumentParser(description="着陆点科学性价值与参数贡献度地图")
    ap.add_argument("--decimate", type=int, default=10)
    ap.add_argument("--weights", choices=["priority", "equal", "entropy"], default="priority")
    ap.add_argument("--zoom-m", type=float, default=120.0, help="站点放大视图边长（米）")
    ap.add_argument("--outdir", default=os.path.join(OUTPUT_DIR, "rasters"))
    ap.add_argument("--no-figures", action="store_true")
    args = ap.parse_args()

    font = _setup_cjk_font()
    print(f"[字体] 使用 {font}" if font else "[字体] 未找到中文字体")
    os.makedirs(args.outdir, exist_ok=True)
    T_START = time.perf_counter()

    k = max(1, int(args.decimate))
    res = SRC_RES_DEG * k
    dst_w, dst_h = SRC_W // k, SRC_H // k
    grid = GridSpec(res_deg=res, lon_min=-180.0, lon_max=-180.0 + res * dst_w,
                    lat_min=90.0 - res * dst_h, lat_max=90.0)
    print(f"=== 格网 === {dst_w}×{dst_h}，{res:.6f}°/像元 ≈{res * DEG2M:.3f} m/px")

    dom, dom_nod = read_decimated(os.path.join(PLANETARY_DIR, DOM_FILE), dst_w, dst_h)
    dom_mask = (dom == dom_nod) if dom_nod is not None else np.zeros(dom.shape, bool)

    sites = pd.read_csv(os.path.join(DATA_DIR, "bennu_sites.csv"))
    fp = pd.read_csv(os.path.join(DATA_DIR, "site_footprints.csv"))
    fp["lon_deg"] = _to_signed_lon(fp["lon_deg"])

    scores = evaluate_sites(sites, weight_scheme=args.weights)
    weights = get_weights(args.weights, scores.set_index("site"))
    eff = effective_weights(weights)
    per_param = score_all_parameters(sites).set_index("site")

    # 贡献度校验：10 个参数贡献之和 == ISVM
    isvm_check = per_param[[p for p in SCIENCE_PARAMETERS]].mul(
        pd.Series(eff)).sum(axis=1)
    ref = scores.set_index("site")["isvm"]
    assert np.allclose(isvm_check.reindex(ref.index), ref, atol=1e-12), "贡献度分解不闭合"
    print("=== 贡献度分解校验通过：Σc(p) == ISVM（容差 1e-12）===")

    # ---- 生成全幅栅格：候选区用站点评分，区外用合成填充（实装测试）----
    t_fill = time.perf_counter()
    shape = grid.shape
    site_order = fp["site"].tolist()
    roi_outline = {}
    fill_mask = np.zeros(shape, dtype="uint8")          # 1 = 候选区，0 = 合成填充
    site_masks = []
    for _, f in fp.iterrows():
        name = f["site"]
        if name not in ref.index:
            continue
        lat, lon = float(f["lat_deg"]), float(f["lon_deg"])
        sl, disk = _site_slice(grid, lat, lon, float(f["footprint_radius_m"]))
        site_masks.append((name, sl, disk))
        mv = fill_mask[sl]
        mv[disk] = 1
        roi_outline[name] = (lat, lon, float(f["roi_diameter_m"]) / 2.0)

    rng_fill = np.random.default_rng(FILL_SEED)
    contrib = {}
    for p in SCIENCE_PARAMETERS:
        sf = _fill_score_field(p, shape, rng_fill)       # 区外填充评分
        for name, sl, disk in site_masks:                # 候选区覆盖为真实评分
            v = sf[sl]
            v[disk] = float(per_param.loc[name, p])
        contrib[p] = (np.float32(eff[p]) * sf).astype("float32")
        del sf
    isvm_map = np.zeros(shape, dtype="float32")
    for p in SCIENCE_PARAMETERS:
        isvm_map += contrib[p]
    t_fill = time.perf_counter() - t_fill

    n_px = int(shape[0] * shape[1])
    n_cand = int(fill_mask.sum())
    print(f"\n=== 全幅填充（合成占位数据）===")
    print(f"  格网 {shape[1]}×{shape[0]} = {n_px/1e6:.2f} M 像元，"
          f"其中候选区 {n_cand} 像元（{100*n_cand/n_px:.2f}%）、"
          f"填充区 {n_px-n_cand} 像元（{100*(n_px-n_cand)/n_px:.2f}%）")
    print(f"  填充规则：连续参数取最差端 {FILL_FRAC*100:.0f}% 域宽 + 平滑噪声（种子 {FILL_SEED}）；"
          f"分类参数取最低档")
    print(f"  生成耗时 {t_fill:.2f} s（10 个参数的填充场 + 评分）")

    cand_min_isvm = float(ref.min())
    fill_max_isvm = float(isvm_map[fill_mask == 0].max())
    print(f"  候选区 ISVM ∈ [{cand_min_isvm:.4f}, {float(ref.max()):.4f}]  |  "
          f"填充区 ISVM ∈ [{float(isvm_map[fill_mask == 0].min()):.4f}, {fill_max_isvm:.4f}]")
    print(f"  填充区最大 ISVM < 候选区最小 ISVM ：{fill_max_isvm < cand_min_isvm} "
          f"（差值 {cand_min_isvm - fill_max_isvm:.4f}）")

    # ---- 写出栅格 ----
    t_write = time.perf_counter()
    write_geotiff(
        os.path.join(args.outdir, "bennu_isvm_landing_value.tiff"),
        isvm_map[None, :, :], ["isvm"], grid,
        band_labels=["着陆点科学性价值 ISVM（全幅：候选区实测 + 区外合成填充）"],
        metadata={"product": "landing_site_science_value", "weight_scheme": args.weights,
                  "footprint_radius_m": "per site_footprints.csv",
                  "fill": "synthetic placeholder: worst-end value + smooth noise",
                  "fill_seed": FILL_SEED, "fill_frac": FILL_FRAC,
                  "mask_band": "bennu_fill_mask.tiff (1=candidate, 0=synthetic)"})

    cube = np.stack([contrib[p] for p in SCIENCE_PARAMETERS]).astype("float32")
    write_geotiff(
        os.path.join(args.outdir, "bennu_param_contribution.tiff"), cube,
        [f"contrib_{p}" for p in SCIENCE_PARAMETERS], grid,
        band_labels=[PARAMETERS[p]["label"] for p in SCIENCE_PARAMETERS],
        metadata={"product": "param_contribution_absolute",
                  "definition": "c(p) = s(p) * w(indicator)/n(indicator); sum == ISVM",
                  "weight_scheme": args.weights,
                  "fill_seed": FILL_SEED, "fill_frac": FILL_FRAC})
    del cube

    share = np.stack([(100.0 * contrib[p] / isvm_map).astype("float32")
                      for p in SCIENCE_PARAMETERS]).astype("float32")
    write_geotiff(
        os.path.join(args.outdir, "bennu_param_contribution_share.tiff"), share,
        [f"share_{p}" for p in SCIENCE_PARAMETERS], grid,
        band_labels=[PARAMETERS[p]["label"] for p in SCIENCE_PARAMETERS],
        metadata={"product": "param_contribution_share_pct",
                  "note": "每个像元 10 个波段之和为 100%"})
    del share

    write_geotiff(
        os.path.join(args.outdir, "bennu_fill_mask.tiff"), fill_mask[None, :, :],
        ["fill_mask"], grid, band_labels=["1=候选区实测/文献值，0=合成填充"],
        dtype="uint8", nodata=None,
        metadata={"product": "fill_mask", "note": "区分实测候选区与合成填充像元"})
    t_write = time.perf_counter() - t_write
    print(f"  栅格写出耗时 {t_write:.2f} s")

    # ---- 明细表 ----
    rows = []
    for name in site_order:
        if name not in ref.index:
            continue
        for p in SCIENCE_PARAMETERS:
            c = float(eff[p] * per_param.loc[name, p])
            rows.append({
                "site": name, "parameter": p, "label": PARAMETERS[p]["label"],
                "indicator": PARAMETERS[p]["indicator"], "score": float(per_param.loc[name, p]),
                "effective_weight": float(eff[p]),
                "contribution": c,
                "share_pct": 100.0 * c / float(ref[name]),
            })
    contrib_df = pd.DataFrame(rows)
    csv_path = os.path.join(OUTPUT_DIR, "exp11_param_contribution.csv")
    contrib_df.to_csv(csv_path, index=False, encoding="utf-8-sig")
    print(f"\n=== 站点 ISVM 与最大贡献参数 ===")
    for name in site_order:
        sub = contrib_df[contrib_df["site"] == name].sort_values(
            "contribution", ascending=False)
        top = sub.iloc[0]
        print(f"  {name:<12} ISVM={ref[name]:.4f} | 最大贡献: {top['label']} "
              f"({top['contribution']:.4f}, {top['share_pct']:.1f}%)")

    # ---- 出图 ----
    if not args.no_figures:
        _figures(args, grid, dom, dom_mask, fp, isvm_map, contrib, ref,
                 roi_outline, contrib_df)

    print("\n=== 输出栅格 ===")
    total_bytes = 0
    for f in ["bennu_isvm_landing_value.tiff", "bennu_param_contribution.tiff",
              "bennu_param_contribution_share.tiff", "bennu_fill_mask.tiff"]:
        p = os.path.join(args.outdir, f)
        size = os.path.getsize(p)
        total_bytes += size
        print(f"  {p}  ({size / 1e6:.2f} MB)")
    print(f"[输出] {csv_path}")
    print(f"\n=== 性能 === 全流程 {time.perf_counter() - T_START:.2f} s"
          f"（填充+评分 {t_fill:.2f} s，栅格写出 {t_write:.2f} s）；"
          f"栅格合计 {total_bytes / 1e6:.2f} MB")


def _panel_base(ax, dom, dom_mask, grid, title, scale=1.0):
    ext = [grid.lon_min, grid.lon_max, grid.lat_min, grid.lat_max]
    bg = np.where(dom_mask, np.nan, dom)
    ax.imshow(bg, extent=ext, origin="upper", cmap="gray",
              vmin=0, vmax=float(np.nanpercentile(bg, 98)))
    ax.set_title(title, fontsize=9 * scale)
    ax.set_facecolor("black")


def _draw_global_panel(ax, grid, dom, dom_mask, isvm_map, fp, sites_order, ref,
                       scale=1.0):
    """绘制「Bennu 全域（DOM）」面板。

    这是 `fig_landing_value_on_bennu.png` 第 1 个子图与
    `fig_dom_global_sites.png` 的**唯一**绘制入口：两图内容完全相同，
    仅通过 `scale` 等比放大字号、线宽与标记大小。
    """
    ext = [grid.lon_min, grid.lon_max, grid.lat_min, grid.lat_max]
    _panel_base(ax, dom, dom_mask, grid, "Bennu 全域（DOM）", scale=scale)
    im = ax.imshow(np.where(np.isnan(isvm_map), np.nan, isvm_map), extent=ext,
                   origin="upper", cmap="turbo", vmin=0.5, vmax=0.8, alpha=OVERLAY_ALPHA)
    for name in sites_order:
        if name in ref.index:
            f = fp[fp["site"] == name].iloc[0]
            ax.annotate(f"{name}\n{ref[name]:.3f}", (f["lon_deg"], f["lat_deg"]),
                        color="white", fontsize=7 * scale, ha="center", va="center",
                        xytext=(f["lon_deg"], f["lat_deg"] - 34),
                        arrowprops=dict(arrowstyle="-", color="cyan", lw=0.8 * scale))
            ax.plot(f["lon_deg"], f["lat_deg"], marker="+", color="cyan",
                    ms=7 * scale, mew=1.2 * scale)
    cb = ax.figure.colorbar(im, ax=ax, fraction=0.03, pad=0.02)
    cb.set_label("ISVM", fontsize=7 * scale)
    cb.ax.tick_params(labelsize=6 * scale)
    ax.set_xlabel("经度 (°)", fontsize=7 * scale)
    ax.set_ylabel("纬度 (°)", fontsize=7 * scale)
    ax.tick_params(labelsize=6 * scale)
    return im


def _figures(args, grid, dom, dom_mask, fp, isvm_map, contrib, ref, roi_outline,
             contrib_df):
    ext = [grid.lon_min, grid.lon_max, grid.lat_min, grid.lat_max]
    sites_order = fp["site"].tolist()

    # 图 1：着陆点科学性价值
    fig, axes = plt.subplots(1, 5, figsize=(22, 4.6))
    _draw_global_panel(axes[0], grid, dom, dom_mask, isvm_map, fp, sites_order, ref,
                       scale=1.0)

    for ax, name in zip(axes[1:], sites_order):
        f = fp[fp["site"] == name].iloc[0]
        lat, lon = float(f["lat_deg"]), float(f["lon_deg"])
        sl, zext = _zoom_extent(grid, lat, lon, args.zoom_m)
        ax.imshow(np.where(dom_mask[sl], np.nan, dom[sl]), extent=zext,
                  origin="upper", cmap="gray", vmin=0,
                  vmax=float(np.nanpercentile(dom[~dom_mask], 98)))
        ax.imshow(np.where(np.isnan(isvm_map[sl]), np.nan, isvm_map[sl]), extent=zext,
                  origin="upper", cmap="turbo", vmin=0.5, vmax=0.8, alpha=OVERLAY_ALPHA)
        ax.add_patch(Circle((0, 0), float(f["footprint_radius_m"]), fill=False,
                            ec="white", lw=0.8, ls="--"))
        r_roi = roi_outline[name][2]
        ax.add_patch(Circle((0, 0), r_roi, fill=False, ec="lime", lw=1.0))
        ax.set_title(f"{name}  ISVM={ref[name]:.3f}\n"
                     f"白虚线圈=25 m 足迹, 绿圈=实测 ROI", fontsize=8)
        ax.set_xlabel("东向 (m)", fontsize=7)
        ax.set_ylabel("北向 (m)", fontsize=7)
        ax.tick_params(labelsize=6)
        ax.set_aspect("equal")
    fig.suptitle("着陆点科学性价值（DOM）", fontsize=13)
    fig.tight_layout()
    fig.canvas.draw()  # 让各 axes 的渲染盒反映最终布局
    # 记录子图实际的刻度位置，供放大版复用，保证两图刻度一致
    global_ticks = (axes[0].get_xticks(), axes[0].get_yticks())
    # 同时记录子图的数据范围：set_xticks() 会把坐标轴范围扩张到刻度范围，
    # 放大版必须把它还原，否则地图会缩到 90% 并在四周留黑边
    global_limits = (axes[0].get_xlim(), axes[0].get_ylim())
    # 用实际渲染盒（含等比例 aspect 收缩后的真实绘图区）宽度，而非子图槽位宽度
    global_panel_width_in = axes[0].get_window_extent().width / fig.dpi
    fig.savefig(os.path.join(args.outdir, "fig_landing_value_on_bennu.png"), dpi=140)
    plt.close(fig)

    # 图 2：每个参数的贡献度地图
    for p in SCIENCE_PARAMETERS:
        arr = contrib[p]
        vmax = max(float(np.nanmax(arr)), 1e-9)
        fig, axes = plt.subplots(1, 5, figsize=(22, 4.4))
        ax = axes[0]
        _panel_base(ax, dom, dom_mask, grid, "Bennu 全域（DOM）")
        im = ax.imshow(np.where(np.isnan(arr), np.nan, arr), extent=ext,
                       origin="upper", cmap=PARAM_CMAP, vmin=0, vmax=vmax,
                       alpha=OVERLAY_ALPHA)
        fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02, label="贡献度")
        for name in sites_order:
            if name in ref.index:
                f = fp[fp["site"] == name].iloc[0]
                ax.plot(f["lon_deg"], f["lat_deg"], marker="+", color="cyan",
                        ms=7, mew=1.2)
                ax.annotate(name, (f["lon_deg"], f["lat_deg"]), color="white",
                            fontsize=6, ha="center", va="bottom",
                            xytext=(f["lon_deg"], f["lat_deg"] - 30),
                            arrowprops=dict(arrowstyle="-", color="cyan", lw=0.7))
        ax.set_xlabel("经度 (°)", fontsize=7)
        ax.set_ylabel("纬度 (°)", fontsize=7)
        ax.tick_params(labelsize=6)

        for ax, name in zip(axes[1:], sites_order):
            f = fp[fp["site"] == name].iloc[0]
            lat, lon = float(f["lat_deg"]), float(f["lon_deg"])
            sl, zext = _zoom_extent(grid, lat, lon, args.zoom_m)
            ax.imshow(np.where(dom_mask[sl], np.nan, dom[sl]), extent=zext,
                      origin="upper", cmap="gray", vmin=0,
                      vmax=float(np.nanpercentile(dom[~dom_mask], 98)))
            im2 = ax.imshow(np.where(np.isnan(arr[sl]), np.nan, arr[sl]), extent=zext,
                            origin="upper", cmap=PARAM_CMAP, vmin=0, vmax=vmax,
                            alpha=OVERLAY_ALPHA)
            ax.add_patch(Circle((0, 0), float(f["footprint_radius_m"]), fill=False,
                                ec="white", lw=0.8, ls="--"))
            val = float(contrib_df[(contrib_df["site"] == name)
                                   & (contrib_df["parameter"] == p)]["contribution"].iloc[0])
            sh = float(contrib_df[(contrib_df["site"] == name)
                                  & (contrib_df["parameter"] == p)]["share_pct"].iloc[0])
            ax.set_title(f"{name}  c={val:.4f} / {sh:.1f}%", fontsize=8)
            ax.set_xlabel("东向 (m)", fontsize=7)
            ax.tick_params(labelsize=6)
            ax.set_aspect("equal")
            fig.colorbar(im2, ax=ax, fraction=0.03, pad=0.02)
        fig.suptitle(f"参数贡献度：{PARAMETERS[p]['label']}（DOM）", fontsize=12)
        fig.tight_layout()
        fig.savefig(os.path.join(args.outdir, f"fig_param_{p}.png"), dpi=140)
        plt.close(fig)

    # 图 3：贡献度堆叠条形图
    fig, ax = plt.subplots(figsize=(11, 5))
    labels = [PARAMETERS[p]["label"] for p in SCIENCE_PARAMETERS]
    bottom = np.zeros(len(sites_order))
    colors = plt.get_cmap("tab20")(np.linspace(0, 1, len(SCIENCE_PARAMETERS)))
    for i, p in enumerate(SCIENCE_PARAMETERS):
        vals = np.array([
            float(contrib_df[(contrib_df["site"] == s)
                             & (contrib_df["parameter"] == p)]["contribution"].iloc[0])
            for s in sites_order])
        ax.bar(sites_order, vals, bottom=bottom, label=labels[i], color=colors[i])
        bottom += vals
    for i, s in enumerate(sites_order):
        ax.text(i, bottom[i] + 0.01, f"ISVM={ref[s]:.3f}", ha="center", fontsize=9)
    ax.set_ylabel("贡献度（ISVM 分解）", fontsize=9)
    ax.set_title("四候选区 ISVM 各参数贡献度（DOM）", fontsize=11)
    ax.legend(fontsize=7, ncol=2, loc="lower right")
    ax.grid(axis="y", alpha=0.3)
    fig.tight_layout()
    fig.savefig(os.path.join(args.outdir, "fig_contribution_bars.png"), dpi=140)
    plt.close(fig)
    _fig_dom_global(args, grid, dom, dom_mask, isvm_map, fp, sites_order, ref,
                    ticks=global_ticks, limits=global_limits,
                    panel_width_in=global_panel_width_in)
    print(f"[图] fig_landing_value_on_bennu.png / fig_param_*.png (10 张) / "
          f"fig_contribution_bars.png / fig_dom_global_sites.png")


def _fig_dom_global(args, grid, dom, dom_mask, isvm_map, fp, sites_order, ref,
                    ticks=None, limits=None, panel_width_in=None):
    """fig_dom_global_sites.png：与可视图第 1 个子图**内容完全相同**的等比放大版。

    两图共用 `_draw_global_panel()`，唯一差别是 `scale`——按两者地图区在纸面上的
    宽度比等比放大字号、线宽与标记，因此放大版看起来就是那块子图的放大。
    """
    fig, ax = plt.subplots(figsize=(24, 12.5))
    scale = 6.0  # 初值，随后按实测宽度比迭代校正
    for _ in range(3):
        ax.clear()
        for cax in list(fig.axes):
            if cax is not ax:
                fig.delaxes(cax)
        _draw_global_panel(ax, grid, dom, dom_mask, isvm_map, fp, sites_order, ref,
                           scale=scale)
        if ticks is not None:
            ax.set_xticks(ticks[0])
            ax.set_yticks(ticks[1])
        if limits is not None:
            # set_xticks/set_yticks 只扩张不收缩，这里强制还原成与子图一致的范围
            ax.set_xlim(limits[0])
            ax.set_ylim(limits[1])
        fig.tight_layout()
        fig.canvas.draw()  # 未重绘时 get_window_extent() 会返回旧布局
        if panel_width_in:
            width_in = ax.get_window_extent().width / fig.dpi
            new_scale = width_in / panel_width_in
            if abs(new_scale - scale) < 0.02:
                scale = new_scale
                break
            scale = new_scale
    out = os.path.join(args.outdir, "fig_dom_global_sites.png")
    fig.savefig(out, dpi=150)
    plt.close(fig)
    ratio_px = scale * (150.0 / 140.0)  # 两图保存 dpi 不同（150 vs 140）
    print(f"[图] {os.path.basename(out)}  = Bennu 全域（DOM）子图的等比放大"
          f"（字号 ×{scale:.2f}，地图区宽度 ×{ratio_px:.2f} 像素）")


if __name__ == "__main__":
    main()
