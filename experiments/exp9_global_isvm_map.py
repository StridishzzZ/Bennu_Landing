# -*- coding: utf-8 -*-
"""实验 9：ISVM 科学价值全球栅格化（Bennu 全幅 GeoTIFF）。

流程
----
1. 复用既有评分流水线 ``src/pipeline.evaluate_sites`` 得到四候选区 C/M/G/T 与 ISVM；
2. 读取矿物丰度图（``--mineral-raster``，未提供则用 Tier-1 代理合成），
   输出为 6 波段 GeoTIFF ``bennu_mineral_abundance.tif``；
3. 以 ``data/global_baseline.csv`` 为全球背景场，把站点评分按足迹核
   （``data/site_footprints.csv``）铺到全幅格网；
4. 输出 4 波段 GeoTIFF ``bennu_isvm_4band.tif``：band1=C, band2=M, band3=G, band4=T，
   即 ISVM 四部分 / 文献四张科学价值图；另输出单波段 ``bennu_isvm_total.tif``；
5. 回读栅格做自检，写出清单 ``outputs/exp9_raster_manifest.csv``。

用法
----
    python experiments/exp9_global_isvm_map.py
    python experiments/exp9_global_isvm_map.py --res 0.25 --weights entropy
    python experiments/exp9_global_isvm_map.py --mineral-raster D:/data/ovirs_abundance.tif
"""

import argparse
import hashlib
import json
import os
from datetime import datetime, timezone

import numpy as np
import pandas as pd

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from common import DATA_DIR, OUTPUT_DIR, load_sites  # noqa: E402
from src.mapping import (  # noqa: E402
    INDICATOR_CHANNELS, background_indicator_scores, compose_isvm,
    footprint_area_summary, render_indicator_fields,
)
from src.mineral_maps import (  # noqa: E402
    MINERAL_LABELS, abundance_to_classes, load_mineral_abundance,
)
from src.pipeline import evaluate_sites  # noqa: E402
from src.raster import default_grid, read_geotiff, write_geotiff  # noqa: E402
from src.weights import get_weights  # noqa: E402

RASTER_DIR = os.path.join(OUTPUT_DIR, "rasters")


def _setup_cjk_font():
    """让 PNG 预览里的中文标题正常显示（Windows 自带微软雅黑/黑体）。"""
    from matplotlib import font_manager

    installed = {f.name for f in font_manager.fontManager.ttflist}
    for candidate in ["Microsoft YaHei", "SimHei", "Noto Sans CJK SC",
                      "Source Han Sans SC", "SimSun", "DengXian"]:
        if candidate in installed:
            matplotlib.rcParams["font.sans-serif"] = [candidate, "DejaVu Sans"]
            matplotlib.rcParams["axes.unicode_minus"] = False
            return candidate
    return None


def _sha256(path, chunk=1 << 20):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(chunk), b""):
            h.update(block)
    return h.hexdigest()


def _plot_mineral_preview(cube, names, labels, footprints, path, grid):
    n = cube.shape[0]
    ncols = 3
    nrows = int(np.ceil(n / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(4.2 * ncols, 2.6 * nrows))
    axes = np.atleast_1d(axes).ravel()
    extent = [0, 360, -90, 90]
    for i in range(n):
        ax = axes[i]
        im = ax.imshow(cube[i], extent=extent, origin="upper", cmap="viridis",
                       vmin=0.0, vmax=float(np.percentile(cube[i], 99.5)))
        ax.set_title(f"{names[i]}\n{labels[i]}", fontsize=8)
        ax.set_xlabel("Longitude (°E)", fontsize=7)
        ax.set_ylabel("Latitude (°N)", fontsize=7)
        ax.tick_params(labelsize=6)
        fig.colorbar(im, ax=ax, fraction=0.035, pad=0.02)
    for ax in axes[n:]:
        ax.axis("off")
    for _, row in footprints.iterrows():
        for ax in axes[:n]:
            ax.plot(row["lon_deg"], row["lat_deg"], marker="+", color="red", ms=5, mew=1.0)
    fig.suptitle("Bennu mineral abundance (Tier-1 proxy) — red crosses = candidate sites",
                 fontsize=10)
    fig.tight_layout()
    fig.savefig(path, dpi=140)
    plt.close(fig)


def _plot_isvm_preview(fields, total, footprints, path, weights, scheme):
    keys = ["C", "M", "G", "T"]
    fig, axes = plt.subplots(1, 5, figsize=(20, 3.6))
    extent = [0, 360, -90, 90]
    panels = [(k, fields[k], f"{k} — {INDICATOR_CHANNELS[k][0]}\n{INDICATOR_CHANNELS[k][1]}") for k in keys]
    panels.append(("ISVM", total, "ISVM composite\nw=" + ", ".join(f"{k}:{weights[k]:.2f}" for k in keys)))
    for ax, (_, arr, title) in zip(axes, panels):
        im = ax.imshow(arr, extent=extent, origin="upper", cmap="magma", vmin=0.0, vmax=1.0)
        ax.set_title(title, fontsize=8)
        ax.set_xlabel("Longitude (°E)", fontsize=7)
        ax.tick_params(labelsize=6)
        fig.colorbar(im, ax=ax, fraction=0.04, pad=0.02)
        for _, row in footprints.iterrows():
            ax.plot(row["lon_deg"], row["lat_deg"], marker="+", color="cyan", ms=5, mew=1.0)
    axes[0].set_ylabel("Latitude (°N)", fontsize=7)
    fig.suptitle(f"ISVM science value on Bennu global map (weight scheme: {scheme})", fontsize=11)
    fig.tight_layout()
    fig.savefig(path, dpi=140)
    plt.close(fig)


def _manifest_rows(info, kind, extra):
    rows = []
    base = {
        "generated_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "raster_path": info.get("path"),
        "kind": kind,
        "n_bands": info.get("bands"),
        "height": info.get("shape", (None, None))[0],
        "width": info.get("shape", (None, None))[1],
        "dtype": info.get("dtype"),
        "backend": info.get("backend_write") or info.get("backend_read"),
        "sha256": extra.get("sha256"),
    }
    base.update({k: v for k, v in extra.items() if k != "sha256"})
    for i, stats in enumerate(info.get("band_stats", []), start=1):
        row = dict(base)
        row["band_index"] = i
        row["band_name"] = info.get("band_names", [None] * i)[i - 1] if info.get("band_names") else None
        row["band_label"] = info.get("band_labels", [""] * i)[i - 1] if info.get("band_labels") else ""
        row.update({
            "min": stats["min"], "max": stats["max"], "mean": stats["mean"],
            "std": stats["std"], "p02": stats["p02"], "p98": stats["p98"],
            "n_valid": stats["n_valid"], "n_nodata": stats["n_nodata"],
        })
        rows.append(row)
    return rows


def main():
    ap = argparse.ArgumentParser(description="ISVM 科学价值全球栅格化（Bennu 全幅）")
    ap.add_argument("--res", type=float, default=0.5,
                    help="格网分辨率（度/像元）。Bennu 1°≈4.28 m，默认 0.5°（≈2.14 m/像元）")
    ap.add_argument("--weights", choices=["priority", "equal", "entropy"], default="priority")
    ap.add_argument("--mineral-raster", default=None,
                    help="外部矿物丰度 GeoTIFF；缺省时使用 Tier-1 代理丰度场")
    ap.add_argument("--seed", type=int, default=20260912, help="代理场随机种子")
    ap.add_argument("--bright-spots", type=int, default=150, help="代理场碳酸盐亮斑数量")
    ap.add_argument("--sigma-scale", type=float, default=1.0, help="候选区足迹半径缩放")
    ap.add_argument("--outdir", default=RASTER_DIR)
    ap.add_argument("--no-preview", action="store_true", help="不输出 PNG 预览")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    font = _setup_cjk_font()
    if font:
        print(f"[字体] PNG 预览使用 {font}")
    else:
        print("[字体] 未找到中文字体，预览图中文可能显示为方框（栅格数据不受影响）")

    # ---- 1. 评分（复用既有流水线）----
    sites = load_sites("bennu_sites.csv")
    scores = evaluate_sites(sites, weight_scheme=args.weights)
    weights = get_weights(args.weights, scores.set_index("site"))

    footprints = load_sites("site_footprints.csv")
    baseline = load_sites("global_baseline.csv")
    background = background_indicator_scores(baseline)

    grid = default_grid(args.res)
    print(f"=== 格网 === {grid.width}×{grid.height} 像元，{args.res}°/像元，"
          f"1°≈{grid.meter_per_degree_lat:.2f} m（经向）")
    print(f"权重方案 {args.weights}: " + ", ".join(f"{k}={weights[k]:.3f}" for k in ["C", "M", "G", "T"]))
    print("全球背景指示分: " + ", ".join(f"{k}={background[k]:.3f}" for k in ["C", "M", "G", "T"]))

    # ---- 2. 矿物丰度图：读取或合成，并输出为 GeoTIFF ----
    abundance, min_names, min_labels, abundance_info = load_mineral_abundance(
        path=args.mineral_raster, grid=grid, footprints=footprints,
        seed=args.seed, n_bright_spots=args.bright_spots,
    )
    n_min_bins = abundance.shape[0]
    min_labels = list(min_labels) if min_labels else [""] * n_min_bins
    if len(MINERAL_LABELS) == n_min_bins and abundance_info["source"] == "synthetic_tier1_proxy":
        min_labels = list(MINERAL_LABELS)

    min_path = os.path.join(args.outdir, "bennu_mineral_abundance.tiff")
    min_info = write_geotiff(
        min_path, abundance, min_names, grid, band_labels=min_labels,
        metadata={"product": "mineral_abundance", "source": abundance_info["source"],
                  "abundance_info": json.dumps(abundance_info, ensure_ascii=False)[:4000]},
    )
    print(f"\n[矿物丰度图] {min_path}")
    print(f"  来源: {abundance_info['source']} | 波段: {min_info['band_names']}")
    print("  各波段均值: " + ", ".join(
        f"{n}={s['mean']:.4f}" for n, s in zip(min_info['band_names'], min_info['band_stats'])))

    classes = abundance_to_classes(abundance, thresholds=[0.05] * n_min_bins)
    print(f"  代理「关键矿物类别数」全球均值 = {classes.mean():.2f}（0–6）")

    # ---- 3. 站点评分 → 全球场 ----
    fields, alpha = render_indicator_fields(
        grid, scores, footprints, background, sigma_scale=args.sigma_scale)
    total = compose_isvm(fields, weights)
    footprint_stats = footprint_area_summary(grid, footprints, sigma_scale=args.sigma_scale)
    print("\n[候选区足迹]")
    print(footprint_stats.to_string(index=False,
                                    columns=["site", "lat_deg", "lon_deg", "footprint_radius_m",
                                             "n_pixels_alpha_ge_0p5", "area_alpha_ge_0p5_m2"]))

    # ---- 4. 输出 4 波段（ISVM 四部分）+ 单波段综合图 ----
    band_names = [f"isvm_{k}" for k in ["C", "M", "G", "T"]]
    band_labels = [f"{INDICATOR_CHANNELS[k][0]} | {INDICATOR_CHANNELS[k][1]}" for k in ["C", "M", "G", "T"]]
    cube4 = np.stack([fields[k] for k in ["C", "M", "G", "T"]]).astype("float32")
    isvm_path = os.path.join(args.outdir, "bennu_isvm_4band.tiff")
    isvm_info = write_geotiff(
        isvm_path, cube4, band_names, grid, band_labels=band_labels,
        metadata={"product": "isvm_science_value_4band", "channels": "C,M,G,T (ISVM four parts)",
                  "weight_scheme": args.weights,
                  "weights": {k: float(weights[k]) for k in ["C", "M", "G", "T"]},
                  "sigma_scale": args.sigma_scale,
                  "background": {k: float(v) for k, v in background.items()}},
    )
    print(f"\n[ISVM 四通道] {isvm_path}")
    for n, lb, s in zip(isvm_info["band_names"], isvm_info["band_labels"], isvm_info["band_stats"]):
        print(f"  band {n:<9} {lb:<32} min={s['min']:.4f} max={s['max']:.4f} mean={s['mean']:.4f}")

    total_path = os.path.join(args.outdir, "bennu_isvm_total.tiff")
    total_info = write_geotiff(
        total_path, total[None, :, :], ["isvm_total"], grid,
        band_labels=[f"ISVM 综合 = " + " + ".join(f"{weights[k]:.2f}·{k}" for k in ["C", "M", "G", "T"])],
        metadata={"product": "isvm_science_value_total", "weight_scheme": args.weights,
                  "weights": {k: float(weights[k]) for k in ["C", "M", "G", "T"]}},
    )
    ts = total_info["band_stats"][0]
    print(f"[ISVM 综合图] {total_path}  min={ts['min']:.4f} max={ts['max']:.4f} mean={ts['mean']:.4f}")

    # ---- 5. 预览 PNG ----
    preview_files = []
    if not args.no_preview:
        p1 = os.path.join(args.outdir, "preview_mineral_abundance.png")
        _plot_mineral_preview(abundance, min_names, min_labels, footprints, p1, grid)
        p2 = os.path.join(args.outdir, "preview_isvm_4band.png")
        _plot_isvm_preview(fields, total, footprints, p2, weights, args.weights)
        preview_files = [p1, p2]
        print("\n[预览图] " + " ; ".join(preview_files))

    # ---- 6. 回读自检 + 清单 ----
    manifest = []
    extra_common = {
        "weight_scheme": args.weights,
        "abundance_source": abundance_info["source"],
        "abundance_seed": args.seed,
        "sigma_scale": args.sigma_scale,
        "res_deg": args.res,
    }
    for path, kind, info in [
        (min_path, "mineral_abundance", min_info),
        (isvm_path, "isvm_4band_C_M_G_T", isvm_info),
        (total_path, "isvm_total", total_info),
    ]:
        _, rb = read_geotiff(path)
        extra = dict(extra_common)
        extra["sha256"] = _sha256(path)
        extra["file_bytes"] = os.path.getsize(path)
        manifest.extend(_manifest_rows(rb, kind, extra))
    manifest_df = pd.DataFrame(manifest)
    manifest_path = os.path.join(OUTPUT_DIR, "exp9_raster_manifest.csv")
    manifest_df.to_csv(manifest_path, index=False, encoding="utf-8-sig")

    scores_out = scores[["site", "C", "M", "G", "T", "isvm"]].copy()
    scores_out["weight_scheme"] = args.weights
    scores_out.to_csv(os.path.join(OUTPUT_DIR, "exp9_site_scores.csv"),
                      index=False, encoding="utf-8-sig")
    footprint_stats.to_csv(os.path.join(OUTPUT_DIR, "exp9_footprint_area.csv"),
                           index=False, encoding="utf-8-sig")

    print(f"\n[自检清单] {manifest_path}（{len(manifest)} 行 = 3 个栅格 × 波段数）")
    print(f"[站点评分] {os.path.join(OUTPUT_DIR, 'exp9_site_scores.csv')}")
    print(f"[足迹统计] {os.path.join(OUTPUT_DIR, 'exp9_footprint_area.csv')}")
    print("\n完成：4 波段 GeoTIFF 的 4 个通道分别对应 ISVM 的 C/M/G/T 四部分。")


if __name__ == "__main__":
    main()
