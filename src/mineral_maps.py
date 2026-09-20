# -*- coding: utf-8 -*-
"""矿物丰度图：外部栅格读取 / Tier-1 代理丰度合成。

背景
----
`MINERAL_ABUNDANCE.md` 的检索结论是：``articles/`` 内没有任何 Bennu 矿物丰度图，
文件夹内最接近的只有反照率图与 MapCam 颜色指数图；真正的定量丰度必须来自
OVIRS/OTES 光谱（PDS，需另行下载）。因此本模块提供两条路径：

1. **读取真实产品**：``load_mineral_abundance(path=...)`` 读取任意 GeoTIFF，
   重采样到目标格网并返回丰度立方体。把 PDS 的矿物/水合产品填进来即可直接使用。
2. **Tier-1 代理合成**：没有外部文件时，按 `MINERAL_ABUNDANCE.md` 第 3 节的
   Tier-1 方案，用文献锚点值合成 6 类矿物的代理丰度场（``synthesize_proxy_abundance``），
   保证每个像元的 6 类丰度之和为 1，且全球水合近均匀、亮斑（碳酸盐）局部富集。

代理场的文献锚点
----------------
* 层状硅酸盐：全球近均匀、占主导（Hamilton et al. 2019；Walsh et al. 2022）
* 碳酸盐：仅亮斑/亮色岩块指示，反照率 10–19%，宿主岩块 6–7.6%
  （Kaplan et al. 2020，经 Jawin et al. 2023 量化）
* 磁铁矿（氧化物）：与 CM/CI 含水蚀变类比、暗色粗糙物质相关
  （Barnouin et al. 2022；DellaGiustina et al. 2020）
* 硫酸盐 / 非晶质：SVMM 优先级第 3、6 位，量级低且斑块化
  （Nakamura-Messenger et al.）

合成结果为**代理场，不是实测矿物丰度**，用于给 ISVM 的 M 层提供空间区分度；
正式发表前应替换为 OVIRS/OTES 反演产品。
"""

from __future__ import annotations

import numpy as np

from .raster import default_grid, gaussian_footprint, read_geotiff, resample_to_grid

# 波段顺序与 src/config.py 中原始六类矿物参数一一对应
MINERAL_BANDS = [
    {"name": "mineral_phyllo", "label": "层状硅酸盐 (phyllosilicate)",
     "priority": 1.0, "role": "SVMM 优先级第 1 位；全球近均匀主导相"},
    {"name": "mineral_carbonate", "label": "碳酸盐 (carbonate)",
     "priority": 0.9, "role": "SVMM 第 2 位；亮斑/亮色岩块的唯一空间指示"},
    {"name": "mineral_sulfate", "label": "硫酸盐 (sulfate)",
     "priority": 0.8, "role": "SVMM 第 3 位；低丰度斑块"},
    {"name": "mineral_oxide", "label": "氧化物 (oxide / magnetite)",
     "priority": 0.7, "role": "SVMM 第 4 位；与暗色粗糙物质相关"},
    {"name": "mineral_silicate", "label": "硅酸盐 (silicate)",
     "priority": 0.6, "role": "SVMM 第 5 位；基质主体"},
    {"name": "mineral_amorphous", "label": "非晶质 (amorphous)",
     "priority": 0.5, "role": "SVMM 第 6 位；与最原始物质相关"},
]

MINERAL_NAMES = [b["name"] for b in MINERAL_BANDS]
MINERAL_LABELS = [b["label"] for b in MINERAL_BANDS]

# 代理场的全球本底丰度（归一化前的相对权重，依据见模块文档字符串）
_BACKGROUND_WEIGHTS = {
    "mineral_phyllo": 0.62,
    "mineral_carbonate": 0.02,
    "mineral_sulfate": 0.015,
    "mineral_oxide": 0.11,
    "mineral_silicate": 0.16,
    "mineral_amorphous": 0.06,
}


def _bilinear_upsample(coarse, height, width):
    """把粗网格双线性放大到 (height, width)。"""
    ch, cw = coarse.shape
    yi = np.linspace(0.0, ch - 1.0, height)
    xi = np.linspace(0.0, cw - 1.0, width)
    tmp = np.empty((height, cw), dtype="float64")
    src_rows = np.arange(ch)
    for j in range(cw):
        tmp[:, j] = np.interp(yi, src_rows, coarse[:, j])
    out = np.empty((height, width), dtype="float64")
    src_cols = np.arange(cw)
    for i in range(height):
        out[i] = np.interp(xi, src_cols, tmp[i])
    return out


def _smooth_noise(rng, shape, coarse=12):
    """可复现的多尺度平滑噪声，值域归一到 [0, 1]。"""
    height, width = shape
    ch = max(3, height // coarse)
    cw = max(3, width // coarse)
    fine = rng.random((ch * 3, cw * 3))
    mid = rng.random((ch, cw))
    field = 0.7 * _bilinear_upsample(mid, height, width) \
        + 0.3 * _bilinear_upsample(fine, height, width)
    lo, hi = float(field.min()), float(field.max())
    return (field - lo) / max(hi - lo, 1e-12)


def _bright_spot_field(grid, rng, n_spots=150, sigma_deg=0.6, site_bias=0.4, footprints=None):
    """碳酸盐亮斑场：文献指出亮斑集中在亮色岩块上，且多见于采样区附近。"""
    height, width = grid.shape
    field = np.zeros((height, width), dtype="float64")

    n_site = int(round(n_spots * site_bias)) if footprints is not None else 0
    n_global = n_spots - n_site

    centers_lat = list(rng.uniform(-80.0, 80.0, size=n_global))
    centers_lon = list(rng.uniform(0.0, 360.0, size=n_global))

    if n_site and footprints is not None:
        sites = footprints.reset_index(drop=True)
        for k in range(n_site):
            row = sites.iloc[int(rng.integers(0, len(sites)))]
            centers_lat.append(float(row["lat_deg"]) + rng.normal(0.0, 1.5))
            centers_lon.append((float(row["lon_deg"]) + rng.normal(0.0, 1.5)) % 360.0)

    for lat, lon in zip(centers_lat, centers_lon):
        sigma_m = sigma_deg * (np.pi * 245.0 / 180.0)
        amp = float(rng.uniform(0.18, 0.5))
        field += amp * gaussian_footprint(grid, lat, lon, sigma_m)

    return np.clip(field, 0.0, 0.75)


def synthesize_proxy_abundance(grid=None, footprints=None, seed=20260912,
                               n_bright_spots=150, spot_sigma_deg=0.6):
    """合成 6 类矿物代理丰度立方体，返回 ``(cube, info)``。

    ``cube`` 形状 ``(6, H, W)``，每个像元 6 类丰度之和为 1（体积/相对丰度分数）。
    """
    if grid is None:
        grid = default_grid()
    rng = np.random.default_rng(int(seed))
    height, width = grid.shape

    noise = {name: _smooth_noise(rng, (height, width)) for name in MINERAL_NAMES}
    spots = _bright_spot_field(grid, rng, n_spots=n_bright_spots,
                               sigma_deg=spot_sigma_deg, footprints=footprints)

    raw = {}
    raw["mineral_phyllo"] = _BACKGROUND_WEIGHTS["mineral_phyllo"] + 0.08 * noise["mineral_phyllo"]
    raw["mineral_silicate"] = _BACKGROUND_WEIGHTS["mineral_silicate"] + 0.05 * noise["mineral_silicate"]
    raw["mineral_oxide"] = _BACKGROUND_WEIGHTS["mineral_oxide"] + 0.04 * noise["mineral_oxide"]
    raw["mineral_amorphous"] = _BACKGROUND_WEIGHTS["mineral_amorphous"] + 0.03 * noise["mineral_amorphous"]
    raw["mineral_sulfate"] = _BACKGROUND_WEIGHTS["mineral_sulfate"] + 0.02 * noise["mineral_sulfate"]
    raw["mineral_carbonate"] = _BACKGROUND_WEIGHTS["mineral_carbonate"] + spots + 0.01 * noise["mineral_carbonate"]

    stack = np.stack([np.clip(raw[name], 1e-4, None) for name in MINERAL_NAMES])
    cube = (stack / stack.sum(axis=0, keepdims=True)).astype("float32")

    info = {
        "source": "synthetic_tier1_proxy",
        "seed": int(seed),
        "n_bright_spots": int(n_bright_spots),
        "spot_sigma_deg": float(spot_sigma_deg),
        "grid": grid.as_dict(),
        "note": "Tier-1 代理丰度场，非实测；替换为 OVIRS/OTES 反演产品后即为定量结果",
        "band_names": MINERAL_NAMES,
        "band_labels": MINERAL_LABELS,
        "band_mean_fraction": {name: float(cube[i].mean()) for i, name in enumerate(MINERAL_NAMES)},
        "band_max_fraction": {name: float(cube[i].max()) for i, name in enumerate(MINERAL_NAMES)},
    }
    return cube, info


def load_mineral_abundance(path=None, grid=None, footprints=None, seed=20260912,
                           n_bright_spots=150, spot_sigma_deg=0.6):
    """读取外部矿物丰度栅格；未提供路径时回退到 Tier-1 代理合成。

    Returns
    -------
    (cube, names, labels, info)
        ``cube`` 形状 ``(bands, H, W)``。
    """
    if grid is None:
        grid = default_grid()

    if path:
        src_cube, src_info = read_geotiff(path)
        bands = src_cube.shape[0]
        names = list(src_info.get("band_names") or [])
        labels = list(src_info.get("band_labels") or [])
        if bands == len(MINERAL_NAMES) and len(names) == bands:
            pass  # 保留文件自带的波段名
        elif len(names) != bands:
            names = [f"band{i + 1}" for i in range(bands)]
            labels = [""] * bands

        gt = src_info.get("geotransform")
        if gt is not None:
            from .raster import GridSpec

            src_grid = GridSpec(
                res_deg=abs(gt[1]),
                lon_min=gt[0],
                lon_max=gt[0] + abs(gt[1]) * src_cube.shape[2],
                lat_max=gt[3],
                lat_min=gt[3] - abs(gt[5]) * src_cube.shape[1],
            )
        else:
            from .raster import GridSpec

            src_grid = GridSpec(res_deg=180.0 / max(src_cube.shape[1], 1),
                                lon_max=360.0)
        cube = resample_to_grid(src_cube, src_grid, grid)
        info = {
            "source": "external_raster",
            "path": src_info.get("path", str(path)),
            "band_names": names,
            "band_labels": labels,
            "grid": grid.as_dict(),
            "source_grid": src_grid.as_dict(),
            "note": "外部栅格按双线性重采样到目标格网",
        }
        return cube, names, labels, info

    cube, info = synthesize_proxy_abundance(
        grid=grid, footprints=footprints, seed=seed,
        n_bright_spots=n_bright_spots, spot_sigma_deg=spot_sigma_deg,
    )
    return cube, MINERAL_NAMES, MINERAL_LABELS, info


def abundance_to_classes(cube, thresholds=None):
    """把丰度立方体折算为「关键矿物类别数」代理（0–6），用于与 M 层参数对接。"""
    cube = np.asarray(cube, dtype="float32")
    if thresholds is None:
        thresholds = [0.05] * cube.shape[0]
    thresholds = np.asarray(thresholds, dtype="float32")[:, None, None]
    return (cube >= thresholds).sum(axis=0).astype("float32")
