# -*- coding: utf-8 -*-
"""Bennu 全球经纬网格定义与 GeoTIFF 读写工具。

坐标约定
--------
* 经度 ``lon`` ∈ [0, 360) °E：Bennu 体固（body-fixed）经度，与 ``articles/``
  文献中给出的候选区坐标一致（例：Nightingale 56.05°N, 42.05°E）。
* 纬度 ``lat`` ∈ [-90, 90] °N。
* 网格为等经纬度（plate carrée）规则格网，像元代表其中心点，行序自北向南，
  与 GeoTIFF 的自然行序一致。
* 由于小天体半径极小（Bennu 平均半径约 245 m），1° 经度仅约 4.28 m，
  因此地图分辨率习惯上用「度/像元」描述，本模块同时给出每像元的米数换算。

写盘优先使用 GDAL（``osgeo``），缺失时回退 ``tifffile``；两者都没有时抛出
带安装提示的错误。读取时自动识别由本模块写出的波段描述与元数据。
"""

from __future__ import annotations

import json
import math
import os
from dataclasses import dataclass
from datetime import datetime, timezone

import numpy as np

# Bennu 平均半径（Lauretta et al. 2019；Barnouin et al. 2020 形状模型 245±0.1 m）
BENNU_MEAN_RADIUS_M = 245.0

# 栅格无效值
NODATA = -9999.0

# Bennu 体固经纬度坐标系的 WKT（球体近似，参考子午线取本初子午线）
BENNU_WKT = (
    'GEOGCS["Bennu",DATUM["Bennu",SPHEROID["Bennu",245.0,0.0]],'
    'PRIMEM["Reference_Meridian",0],UNIT["degree",0.0174532925199433]]'
)


def _try_gdal():
    try:
        from osgeo import gdal  # noqa: F401

        return gdal
    except Exception:
        return None


def _try_tifffile():
    try:
        import tifffile

        return tifffile
    except Exception:
        return None


@dataclass(frozen=True)
class GridSpec:
    """等经纬度规则格网。"""

    res_deg: float
    lon_min: float = 0.0
    lon_max: float = 360.0
    lat_min: float = -90.0
    lat_max: float = 90.0

    @property
    def width(self) -> int:
        return int(round((self.lon_max - self.lon_min) / self.res_deg))

    @property
    def height(self) -> int:
        return int(round((self.lat_max - self.lat_min) / self.res_deg))

    @property
    def shape(self):
        return (self.height, self.width)

    @property
    def n_pixels(self) -> int:
        return self.width * self.height

    @property
    def lon_centers(self) -> np.ndarray:
        return self.lon_min + (np.arange(self.width) + 0.5) * self.res_deg

    @property
    def lat_centers(self) -> np.ndarray:
        # 北→南，行序与 GeoTIFF 一致
        return self.lat_max - (np.arange(self.height) + 0.5) * self.res_deg

    @property
    def geotransform(self):
        """GDAL 地理变换：(x0, dx, 0, y0, 0, -dy)，左上角像元角点。"""
        return (self.lon_min, self.res_deg, 0.0, self.lat_max, 0.0, -self.res_deg)

    @property
    def meter_per_degree_lat(self) -> float:
        return math.pi * BENNU_MEAN_RADIUS_M / 180.0

    def meter_per_pixel(self, lat_deg) -> np.ndarray:
        """给定纬度处的每像元米数（经向、纬向）。"""
        lat = np.asarray(lat_deg, dtype="float64")
        m_per_deg_lon = self.meter_per_degree_lat * np.cos(np.radians(lat))
        return np.stack([m_per_deg_lon * self.res_deg, self.meter_per_degree_lat * self.res_deg])

    def as_dict(self) -> dict:
        return {
            "res_deg": self.res_deg,
            "width": self.width,
            "height": self.height,
            "lon_min": self.lon_min,
            "lon_max": self.lon_max,
            "lat_min": self.lat_min,
            "lat_max": self.lat_max,
            "mean_radius_m": BENNU_MEAN_RADIUS_M,
        }


def default_grid(res_deg: float = 0.5) -> GridSpec:
    """Bennu 全幅（0–360°E，90°S–90°N）格网。"""
    return GridSpec(res_deg=float(res_deg))


def wrap_lon(dlon_deg):
    """经度差归一到 [-180, 180)。"""
    return (np.asarray(dlon_deg, dtype="float64") + 180.0) % 360.0 - 180.0


def gaussian_footprint(grid: GridSpec, lat_deg, lon_deg, sigma_m, radius_m=BENNU_MEAN_RADIUS_M):
    """候选区足迹核：以 (lat, lon) 为心、σ=``sigma_m``（米）的二维高斯。

    经度方向按 cos(lat) 折算，并做 0–360° 环绕处理，因此跨 0° 经线不会断裂。
    """
    sigma_lat_deg = math.degrees(sigma_m / radius_m)
    cos_lat = max(math.cos(math.radians(float(lat_deg))), 1e-6)
    sigma_lon_deg = sigma_lat_deg / cos_lat

    dlat = grid.lat_centers[:, None] - float(lat_deg)
    dlon = wrap_lon(grid.lon_centers[None, :] - float(lon_deg))

    return np.exp(-0.5 * ((dlat / sigma_lat_deg) ** 2 + (dlon / sigma_lon_deg) ** 2))


def _ensure_dir(path):
    d = os.path.dirname(os.path.abspath(path))
    if d:
        os.makedirs(d, exist_ok=True)


def _band_description(name, label=None):
    return f"{name} | {label}" if label else str(name)


def write_geotiff(path, cube, band_names, grid: GridSpec, band_labels=None,
                  metadata=None, compress="DEFLATE", nodata=NODATA, dtype="float32"):
    """写出多波段 GeoTIFF，返回写出文件的元信息字典。

    Parameters
    ----------
    cube : ndarray
        形状 ``(bands, height, width)`` 或 ``(height, width)``（自动升维）。
    band_names : list[str]
        波段名，写入 GeoTIFF 波段描述，便于后续用 GIS 软件识别通道含义。
    grid : GridSpec
        目标格网，决定地理变换与投影。
    """
    arr = np.asarray(cube, dtype=dtype)
    if arr.ndim == 2:
        arr = arr[None, :, :]
    if arr.ndim != 3:
        raise ValueError(f"cube 必须是 (bands, H, W)，收到 {arr.shape}")
    if (arr.shape[1], arr.shape[2]) != grid.shape:
        raise ValueError(
            f"cube 空间尺寸 {(arr.shape[1], arr.shape[2])} 与格网 {grid.shape} 不一致"
        )
    if len(band_names) != arr.shape[0]:
        raise ValueError(f"band_names 数量 {len(band_names)} 与波段数 {arr.shape[0]} 不一致")

    _ensure_dir(path)
    labels = list(band_labels) if band_labels is not None else [None] * arr.shape[0]
    meta = {
        "generator": "src/raster.py",
        "created_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "grid": grid.as_dict(),
        "bands": [
            {"index": i + 1, "name": n, "label": l}
            for i, (n, l) in enumerate(zip(band_names, labels))
        ],
        "site": "Bennu (101955)",
        "crs_note": "Bennu body-fixed longitude/latitude, 0-360E, sphere R=245 m",
    }
    if metadata:
        meta.update(metadata)
    meta_json = json.dumps(meta, ensure_ascii=False)

    gdal = _try_gdal()
    backend = None
    if gdal is not None:
        driver = gdal.GetDriverByName("GTiff")
        ds = driver.Create(
            str(path), grid.width, grid.height, arr.shape[0], gdal.GDT_Float32,
            options=["COMPRESS=" + str(compress), "TILED=YES", "PREDICTOR=3"],
        )
        if ds is None:
            raise RuntimeError(f"GDAL 无法创建 {path}")
        ds.SetGeoTransform(grid.geotransform)
        ds.SetProjection(BENNU_WKT)
        for i in range(arr.shape[0]):
            band = ds.GetRasterBand(i + 1)
            band.SetDescription(_band_description(band_names[i], labels[i]))
            band.SetNoDataValue(float(nodata))
            band.WriteArray(arr[i])
        ds.SetMetadataItem("CODEX_BENNU_META", meta_json)
        ds.FlushCache()
        ds = None
        backend = "gdal"
    else:
        tifffile = _try_tifffile()
        if tifffile is None:
            raise RuntimeError(
                "写 GeoTIFF 需要 GDAL(osgeo) 或 tifffile，请先安装其一："
                "conda install -c conda-forge gdal 或 pip install tifffile"
            )
        tifffile.imwrite(
            str(path), arr, photometric="minisblack", compression="deflate",
            metadata={"axes": "CYX"}, description=meta_json,
        )
        backend = "tifffile"

    _, info = read_geotiff(path)
    info["backend_write"] = backend
    return info


def read_decimated(path, dst_w, dst_h, dtype="float64"):
    """按目标尺寸读取单波段栅格，返回 ``(array, nodata)``。

    GDAL 会对源栅格做块平均降采样；用于把 5 cm/px 的全球产品降到可处理的格网。
    """
    gdal = _try_gdal()
    if gdal is None:
        raise RuntimeError("read_decimated 需要 GDAL(osgeo)，请先安装")
    ds = gdal.Open(str(path))
    if ds is None:
        raise FileNotFoundError(f"无法打开栅格：{path}")
    band = ds.GetRasterBand(1)
    arr = band.ReadAsArray(0, 0, ds.RasterXSize, ds.RasterYSize,
                           buf_xsize=int(dst_w), buf_ysize=int(dst_h)).astype(dtype)
    nodata = band.GetNoDataValue()
    ds = None
    return arr, nodata


def read_geotiff(path):
    """读取 GeoTIFF，返回 ``(cube, info)``；``cube`` 形状为 ``(bands, H, W)``。"""
    gdal = _try_gdal()
    if gdal is not None:
        ds = gdal.Open(str(path))
        if ds is None:
            raise FileNotFoundError(f"无法打开栅格：{path}")
        cube = np.stack([ds.GetRasterBand(i + 1).ReadAsArray()
                         for i in range(ds.RasterCount)]).astype("float32")
        descs = [ds.GetRasterBand(i + 1).GetDescription() for i in range(ds.RasterCount)]
        band_names = [(d.split("|")[0].strip() if d else f"band{i + 1}")
                      for i, d in enumerate(descs)]
        band_labels = [(d.split("|", 1)[1].strip() if d and "|" in d else "")
                       for d in descs]
        raw_meta = ds.GetMetadataItem("CODEX_BENNU_META")
        info = {
            "path": os.path.abspath(str(path)),
            "bands": ds.RasterCount,
            "band_names": band_names,
            "band_labels": band_labels,
            "shape": (ds.RasterYSize, ds.RasterXSize),
            "dtype": str(cube.dtype),
            "geotransform": ds.GetGeoTransform(),
            "projection_wkt": ds.GetProjection(),
            "backend_read": "gdal",
        }
        ds = None
    else:
        tifffile = _try_tifffile()
        if tifffile is None:
            raise RuntimeError("读 GeoTIFF 需要 GDAL(osgeo) 或 tifffile，请先安装其一")
        with tifffile.TiffFile(str(path)) as tf:
            cube = np.asarray(tf.asarray(), dtype="float32")
            if cube.ndim == 2:
                cube = cube[None, :, :]
            desc = tf.pages[0].description or ""
        meta = {}
        try:
            meta = json.loads(desc)
        except Exception:
            meta = {}
        names = [b.get("name", f"band{i + 1}") for i, b in enumerate(meta.get("bands", []))]
        if len(names) != cube.shape[0]:
            names = [f"band{i + 1}" for i in range(cube.shape[0])]
        info = {
            "path": os.path.abspath(str(path)),
            "bands": cube.shape[0],
            "band_names": names,
            "band_labels": [b.get("label", "") for b in meta.get("bands", [])],
            "shape": (cube.shape[1], cube.shape[2]),
            "dtype": str(cube.dtype),
            "geotransform": None,
            "projection_wkt": None,
            "backend_read": "tifffile",
            "meta": meta,
        }
    info["band_stats"] = [
        band_stats(cube[i], info.get("geotransform")) for i in range(cube.shape[0])
    ]
    return cube, info


def band_stats(arr, geotransform=None):
    """单波段统计量，用于写盘后的自检清单。"""
    a = np.asarray(arr, dtype="float64")
    finite = np.isfinite(a)
    n_valid = int(finite.sum())
    if n_valid == 0:
        return {"n_valid": 0, "n_nodata": int(a.size), "min": None, "max": None,
                "mean": None, "std": None, "p02": None, "p98": None}
    v = a[finite]
    return {
        "n_valid": n_valid,
        "n_nodata": int(a.size - n_valid),
        "min": float(v.min()),
        "max": float(v.max()),
        "mean": float(v.mean()),
        "std": float(v.std()),
        "p02": float(np.percentile(v, 2)),
        "p98": float(np.percentile(v, 98)),
    }


def resample_to_grid(cube, src_grid: GridSpec, dst_grid: GridSpec):
    """把 ``(bands, H, W)`` 栅格双线性重采样到目标格网（经度方向环绕）。"""
    cube = np.asarray(cube, dtype="float32")
    if cube.ndim == 2:
        cube = cube[None, :, :]
    if (src_grid.width, src_grid.height) == (dst_grid.width, dst_grid.height) \
            and abs(src_grid.res_deg - dst_grid.res_deg) < 1e-12:
        return cube.copy()

    src_lon = src_grid.lon_centers
    src_lat = src_grid.lat_centers
    dst_lon = dst_grid.lon_centers
    dst_lat = dst_grid.lat_centers

    # 经度按比例映射后环绕取模；纬度按单调递减映射并裁剪
    fx = (dst_lon - src_lon[0]) / ((src_lon[-1] - src_lon[0]) or 1.0) * (src_grid.width - 1)
    fx = np.mod(fx, src_grid.width)
    fy = (src_lat[0] - dst_lat) / ((src_lat[0] - src_lat[-1]) or 1.0) * (src_grid.height - 1)
    fy = np.clip(fy, 0, src_grid.height - 1)

    x0 = np.floor(fx).astype(int)
    y0 = np.floor(fy).astype(int)
    x1 = (x0 + 1) % src_grid.width
    y1 = np.minimum(y0 + 1, src_grid.height - 1)
    tx = (fx - x0)[None, :]
    ty = (fy - y0)[:, None]

    out = np.empty((cube.shape[0], dst_grid.height, dst_grid.width), dtype="float32")
    for b in range(cube.shape[0]):
        a = cube[b]
        top = a[np.ix_(y0, x0)] * (1 - tx) + a[np.ix_(y0, x1)] * tx
        bot = a[np.ix_(y1, x0)] * (1 - tx) + a[np.ix_(y1, x1)] * tx
        out[b] = top * (1 - ty) + bot * ty
    return out
