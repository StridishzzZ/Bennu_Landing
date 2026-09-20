# -*- coding: utf-8 -*-
"""实验 10：把 ISVM 科学价值叠加到真实 Bennu 底图（DOM + DEM）。

输入（`data/planetary/`，见该目录 README）
----------------------------------------
* ``Bennu_DOM_DEG.tif``：全球正射影像，8 bit，0.01145876°/px（≈5 cm/px）
* ``Bennu_DEM_DEG.tif``：全球 DEM，float32，同一格网

输出（`--outdir`，默认 ``outputs/rasters``）
------------------------------------------
* ``bennu_dom_base.tiff``            真实底图（降采样）
* ``bennu_dem_slope.tiff``           由 DEM 计算的真实坡度（度）
* ``bennu_isvm_4band_domgrid.tiff``  四通道科学价值图，**与底图逐像元对齐**
* ``preview_real_basemap.png``       底图 + 高程 + 坡度 + 四通道叠合预览
* ``outputs/exp10_site_terrain.csv`` 四站点实测地形统计

坐标约定
--------
源产品为地理坐标 ``lon ∈ [-180, 180]``；本项目其余模块用 0–360°E，
本脚本在读取时统一按 (-180, 180] 处理，站点经度做 180° 折算。
注意 ``Bennu_DOM.tif`` / ``Bennu_DEM.tif``（米制 eqc 版本）与 DEG 版本相差
**180° 经度**（前者 x=0 对应 180°E），因此本脚本只用 DEG 版本，避免混用错位。

用法
----
    python experiments/exp10_real_basemap.py                  # 10× 降采样 ≈0.5 m/px
    python experiments/exp10_real_basemap.py --decimate 5     # ≈0.25 m/px
"""

import argparse
import os

import numpy as np
import pandas as pd

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from osgeo import gdal  # noqa: E402

gdal.UseExceptions()

from common import DATA_DIR, OUTPUT_DIR  # noqa: E402
from src.mapping import background_indicator_scores, render_indicator_fields  # noqa: E402
from src.pipeline import evaluate_sites  # noqa: E402
from src.raster import BENNU_MEAN_RADIUS_M, GridSpec, write_geotiff  # noqa: E402
from src.weights import get_weights  # noqa: E402

PLANETARY_DIR = os.path.join(DATA_DIR, "planetary")
DOM_FILE = "Bennu_DOM_DEG.tif"
DEM_FILE = "Bennu_DEM_DEG.tif"
SRC_RES_DEG = 0.01145876  # 360 / 31417
SRC_W, SRC_H = 31417, 15709
DEG2M = np.pi * BENNU_MEAN_RADIUS_M / 180.0  # 每度米数（R=245 m）


def _setup_cjk_font():
    from matplotlib import font_manager

    installed = {f.name for f in font_manager.fontManager.ttflist}
    for cand in ["Microsoft YaHei", "SimHei", "Noto Sans CJK SC", "SimSun"]:
        if cand in installed:
            matplotlib.rcParams["font.sans-serif"] = [cand, "DejaVu Sans"]
            matplotlib.rcParams["axes.unicode_minus"] = False
            return cand
    return None


def _read_decimated(path, dst_w, dst_h, dtype="float64"):
    ds = gdal.Open(path)
    if ds is None:
        raise FileNotFoundError(f"无法打开 {path}")
    band = ds.GetRasterBand(1)
    arr = band.ReadAsArray(0, 0, ds.RasterXSize, ds.RasterYSize,
                           buf_xsize=dst_w, buf_ysize=dst_h).astype(dtype)
    nodata = band.GetNoDataValue()
    ds = None
    return arr, nodata


def _slope_degrees(elev, grid, dx_deg, dy_deg, radius_m=BENNU_MEAN_RADIUS_M):
    """由高程场计算坡度（度）。经度方向米数按 cos(lat) 修正。"""
    m_per_deg = np.pi * radius_m / 180.0
    dy_m = abs(dy_deg) * m_per_deg
    dx_m = abs(dx_deg) * m_per_deg * np.cos(np.radians(grid.lat_centers))[:, None]
    gy = np.gradient(elev, dy_m, axis=0)
    gx = np.gradient(elev, 1.0, axis=1) / np.broadcast_to(dx_m, elev.shape)
    return np.degrees(np.arctan(np.hypot(gx, gy))).astype("float32")


def _to_signed_lon(lon):
    return ((np.asarray(lon, dtype="float64") + 180.0) % 360.0) - 180.0


def _footprint_mask(grid, lat, lon, radius_m):
    dlat = (grid.lat_centers[:, None] - lat) * DEG2M
    dlon = ((grid.lon_centers[None, :] - lon + 180.0) % 360.0 - 180.0) \
        * DEG2M * max(np.cos(np.radians(lat)), 1e-6)
    return np.hypot(dlat, dlon) <= radius_m


def main():
    ap = argparse.ArgumentParser(description="ISVM 叠加到真实 Bennu DOM/DEM 底图")
    ap.add_argument("--decimate", type=int, default=10, help="降采样倍数（10 → ≈0.5 m/px）")
    ap.add_argument("--weights", choices=["priority", "equal", "entropy"], default="priority")
    ap.add_argument("--outdir", default=os.path.join(OUTPUT_DIR, "rasters"))
    ap.add_argument("--no-preview", action="store_true")
    ap.add_argument("--write-slope", action="store_true",
                    help="额外输出 float32 坡度栅格（约 15 MB，默认不写以控制仓库体积）")
    args = ap.parse_args()

    _setup_cjk_font()
    os.makedirs(args.outdir, exist_ok=True)

    k = max(1, int(args.decimate))
    res = SRC_RES_DEG * k
    dst_w, dst_h = SRC_W // k, SRC_H // k
    grid = GridSpec(res_deg=res, lon_min=-180.0, lon_max=-180.0 + res * dst_w,
                    lat_min=90.0 - res * dst_h, lat_max=90.0)
    print(f"=== 工作格网 === {dst_w}×{dst_h} 像元，{res:.6f}°/像元，"
          f"≈{res * DEG2M:.3f} m/像元（{k}× 降采样）")

    dom, dom_nod = _read_decimated(os.path.join(PLANETARY_DIR, DOM_FILE), dst_w, dst_h)
    dem, dem_nod = _read_decimated(os.path.join(PLANETARY_DIR, DEM_FILE), dst_w, dst_h)
    dom_mask = (dom == dom_nod) if dom_nod is not None else np.zeros(dom.shape, bool)
    dem_mask = (dem == dem_nod) if dem_nod is not None else np.zeros(dem.shape, bool)
    dem_valid = dem[~dem_mask]
    print(f"DOM: {dom.min():.0f}–{dom.max():.0f}，无效 {100 * dom_mask.mean():.2f}% | "
          f"DEM: {dem_valid.min():.2f}–{dem_valid.max():.2f} m，"
          f"无效 {100 * dem_mask.mean():.2f}%")

    dem_f = np.where(dem_mask, np.nan, dem)
    slope = _slope_degrees(dem_f, grid, res, -res)
    slope = np.where(np.isfinite(slope), slope, 0.0).astype("float32")

    # ---- 站点评分与四通道 ----
    sites = pd.read_csv(os.path.join(DATA_DIR, "bennu_sites.csv"))
    scores = evaluate_sites(sites, weight_scheme=args.weights)
    weights = get_weights(args.weights, scores.set_index("site"))
    baseline = pd.read_csv(os.path.join(DATA_DIR, "global_baseline.csv"))
    background = background_indicator_scores(baseline)

    footprints = pd.read_csv(os.path.join(DATA_DIR, "site_footprints.csv"))
    footprints["lon_deg"] = _to_signed_lon(footprints["lon_deg"])
    fields, alpha = render_indicator_fields(grid, scores, footprints, background)
    cube4 = np.stack([fields[c] for c in ["C", "M", "G", "T"]]).astype("float32")

    # ---- 输出栅格 ----
    dom_path = os.path.join(args.outdir, "bennu_dom_base.tiff")
    write_geotiff(dom_path, dom[None, :, :].astype("uint8"), ["dom_pan"], grid,
                  band_labels=["Bennu 正射影像（降采样，8 bit）"],
                  dtype="uint8", nodata=int(dom_nod) if dom_nod is not None else 0)
    slope_path = None
    if args.write_slope:
        slope_path = os.path.join(args.outdir, "bennu_dem_slope.tiff")
        write_geotiff(slope_path, slope[None, :, :], ["slope_deg"], grid,
                      band_labels=["DEM 计算坡度（度）"], dtype="float32", nodata=-9999.0)
    isvm_path = os.path.join(args.outdir, "bennu_isvm_4band_domgrid.tiff")
    write_geotiff(isvm_path, cube4, ["isvm_C", "isvm_M", "isvm_G", "isvm_T"], grid,
                  band_labels=["SVCCM 化学组成", "SVMM 矿物学", "SVGFM 地质特征",
                               "SVTM 温度"],
                  metadata={"product": "isvm_4band_on_real_dom_grid",
                            "source_basemap": DOM_FILE, "decimate": k,
                            "weight_scheme": args.weights})

    # ---- 站点实测地形统计 ----
    rows = []
    for _, fp in footprints.iterrows():
        for tag, radius in [("ROI", float(fp["roi_diameter_m"]) / 2.0), ("25m", 25.0)]:
            m = _footprint_mask(grid, fp["lat_deg"], fp["lon_deg"], radius) & ~dem_mask
            if not m.any():
                continue
            s = slope[m]
            e = dem_f[m]
            rows.append({
                "site": fp["site"], "scale": tag, "radius_m": radius, "n_px": int(m.sum()),
                "slope_median_deg": float(np.median(s)),
                "slope_p90_deg": float(np.percentile(s, 90)),
                "slope_max_deg": float(s.max()),
                "elev_std_m": float(np.nanstd(e)),
            })
    terrain = pd.DataFrame(rows)
    terrain.to_csv(os.path.join(OUTPUT_DIR, "exp10_site_terrain.csv"),
                   index=False, encoding="utf-8-sig")

    roi = terrain[terrain["scale"] == "ROI"].merge(
        sites[["site", "slope_deg", "roughness"]], on="site")
    print("\n=== 实测坡度 vs 项目假设值（ROI 尺度）===")
    print(roi[["site", "radius_m", "n_px", "slope_median_deg", "slope_p90_deg",
               "slope_deg", "roughness"]].round(2).to_string(index=False))
    print("  参考：Barnouin et al. 2022 给出 Nightingale 接触点局部坡度为 8°–10°")

    if not args.no_preview:
        ext = [grid.lon_min, grid.lon_max, grid.lat_min, grid.lat_max]
        fig, axes = plt.subplots(2, 3, figsize=(17, 8))
        panels = [
            (dom, "Bennu DOM 正射影像（真实底图）", "gray", 0, 110),
            (dem, "DEM 高程 (m)", "terrain", None, None),
            (slope, "DEM 坡度 (°)", "magma", 0, 40),
            (cube4[0], "ISVM-C 化学组成", "magma", 0, 1),
            (cube4[1], "ISVM-M 矿物学", "magma", 0, 1),
        ]
        for ax, (arr, title, cmap, vmin, vmax) in zip(axes.ravel()[:5], panels):
            im = ax.imshow(arr, extent=ext, origin="upper", cmap=cmap,
                           vmin=vmin, vmax=vmax)
            ax.set_title(title, fontsize=9)
            fig.colorbar(im, ax=ax, fraction=0.04, pad=0.02)
        ax = axes.ravel()[5]
        bg = np.clip(dom / max(float(np.percentile(dom[~dom_mask], 98)), 1e-6), 0, 1)
        ax.imshow(bg, extent=ext, origin="upper", cmap="gray")
        rgba = np.zeros((*cube4.shape[1:], 4), dtype="float32")
        rgba[..., 0], rgba[..., 1], rgba[..., 2] = cube4[1], cube4[2], cube4[3]
        rgba[..., 3] = np.clip(alpha * 1.2, 0, 1)
        ax.imshow(rgba, extent=ext, origin="upper")
        ax.set_title("ISVM 叠加在真实 DOM 上（RGB = M/G/T）", fontsize=9)
        for ax in axes.ravel():
            ax.set_xlabel("经度 (°)", fontsize=7)
            ax.set_ylabel("纬度 (°)", fontsize=7)
            ax.tick_params(labelsize=6)
            for _, fp in footprints.iterrows():
                ax.plot(fp["lon_deg"], fp["lat_deg"], marker="+", color="cyan",
                        ms=5, mew=1.0)
        fig.suptitle("实验 10：ISVM 科学价值 × 真实 Bennu DOM/DEM", fontsize=12)
        fig.tight_layout()
        prev = os.path.join(args.outdir, "preview_real_basemap.png")
        fig.savefig(prev, dpi=140)
        plt.close(fig)
        print(f"\n[预览] {prev}")

    for p in [q for q in [dom_path, slope_path, isvm_path] if q]:
        print(f"[输出] {p}  ({os.path.getsize(p) / 1e6:.2f} MB)")
    print(f"[输出] {os.path.join(OUTPUT_DIR, 'exp10_site_terrain.csv')}")
    print("\n完成：四通道科学价值图已与真实 DOM 底图逐像元对齐。")


if __name__ == "__main__":
    main()
