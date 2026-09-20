# Bennu_landing — 小天体着陆区科学价值综合评价

面向小天体（贝努 / 龙宫）着陆/采样区科学价值综合评价的实验项目。

## 设计前提（与用户约定）

1. 输入为**直接给定的参数**（无遥感数据），每个候选区一行参数；
2. 现阶段**不涉及聚类算法**；
3. 指标体系采用文献一致的**四张科学价值图**（SVCCM / SVMM / SVGFM / SVTM），
   可采样性、安全性、可达性作为独立工程层门控；
4. 权重默认三方案对比：文献优先级权重、等权、熵权；
5. 验证基准：贝努四候选区（Nightingale / Osprey / Kingfisher / Sandpiper）、
   龙宫七候选区多级下选（L05–L12、M01–M04 → L08 → L08-B / L08-E1）；
6. 方法设计为可迁移到天问二号目标小行星（参数来源可替换）。

自 2026-08-26 起，科学参数经冗余审计从 21 个精简为 10 个
（常数参数删除、重复参数合并，详见 `data/README.md` 与 `experiments/exp7_redundancy_audit.py`），
精简前后所有验收结论一致。

## 目录

```
Bennu_landing/
├── README.md
├── requirements.txt
├── config.py → src/config.py        # 指标体系与参数定义（含文献依据）
├── data/                            # 输入模板 + 精简后文献参数数据集
│   ├── full/                        # 精简前的完整参数集（归档）
│   ├── parameter_template.csv
│   ├── bennu_sites.csv
│   ├── ryugu_sites.csv
│   ├── ryugu_subsites.csv
│   ├── site_footprints.csv           # 贝努四候选区文献坐标 + 足迹半径（实验 9）
│   ├── global_baseline.csv           # 全球背景基线参数（实验 9）
│   └── README.md                    # 数据来源与置信度说明
├── src/                             # 评分流水线
│   ├── config.py
│   ├── scoring.py                   # 单参数评分函数
│   ├── indicators.py                # 四图指标聚合
│   ├── weights.py                   # 三种权重方案
│   ├── pipeline.py                  # ISVM 与工程层门控
│   ├── analysis.py                  # 敏感性 / 蒙特卡洛
│   ├── raster.py                    # Bennu 经纬格网 + GeoTIFF 读写
│   ├── mineral_maps.py              # 矿物丰度图读取 / Tier-1 代理合成
│   └── mapping.py                   # ISVM 评分的全球栅格化
├── experiments/                     # 实验脚本
│   ├── common.py
│   ├── exp1_bennu_validation.py
│   ├── exp2_ryugu_validation.py
│   ├── exp3_weight_comparison.py
│   ├── exp4_sensitivity.py
│   ├── exp5_monte_carlo.py
│   ├── exp6_synthetic_stress.py
│   ├── exp7_redundancy_audit.py
│   ├── exp8_before_after_compare.py
│   ├── exp9_global_isvm_map.py      # 科学价值全球栅格化（GeoTIFF）
│   ├── exp10_real_basemap.py        # 叠加真实 DOM/DEM 底图
│   └── exp11_contribution_maps.py   # 着陆点价值 + 各参数贡献度地图
├── data/planetary/                  # 真实 Bennu DOM/DEM（默认不入库，见其 README）
└── outputs/                         # 实验结果（CSV）+ rasters/（GeoTIFF）
```

## 运行

```bash
python experiments/exp1_bennu_validation.py
python experiments/exp2_ryugu_validation.py
python experiments/exp3_weight_comparison.py
python experiments/exp4_sensitivity.py
python experiments/exp5_monte_carlo.py
python experiments/exp6_synthetic_stress.py
python experiments/exp9_global_isvm_map.py      # 需要 matplotlib + GDAL/tifffile
python experiments/exp10_real_basemap.py        # 需要真实 DOM/DEM 放在 data/planetary/
python experiments/exp11_contribution_maps.py   # 着陆点科学价值 + 参数贡献度出图
```

结果写入 `outputs/`。

> 本机运行环境为 conda 环境 `Yolo_cuda1210`（含 pandas / matplotlib / tifffile / GDAL）；
> `requirements.txt` 只列出评分核心的最小依赖，栅格功能另见该文件的说明。

实验结论摘要见 [RESULTS.md](RESULTS.md)。

## 指标结构

ISVM = wC·C + wM·M + wG·G + wT·T，其中：

- C 化学组成：有机物丰度、CH2/CH3、碳含量；
- M 矿物学：矿物类别数、水合证据；
- G 地质新鲜度：空间风化代理、陨石坑新鲜度、活动证据、PSFD；
- T 温度：相对全球均值的温度等级；
- 工程层（门控，不进 ISVM）：安全性（坡度/块石覆盖/Tmax/倾角）、
  可采样性（≤2 cm 细粒占比/粒径/倾角）、可达性（目标半径/导航误差）。

权重方案：

- `priority`：C=0.35, M=0.30, G=0.20, T=0.15（按文献任务优先级递减，属设计假设）；
- `equal`：各 0.25；
- `entropy`：由当前数据集的指标得分熵确定。

## 全球科学价值栅格（实验 9）

评分完成后，`experiments/exp9_global_isvm_map.py` 会把科学价值铺到 Bennu 全幅地图上，
并输出本次读取（或合成）的矿物丰度图。产物在 `outputs/rasters/`：

| 文件 | 内容 | 波段 |
|---|---|---|
| `bennu_isvm_4band.tiff` | **科学价值四通道地图** | band1=`isvm_C`(SVCCM)、band2=`isvm_M`(SVMM)、band3=`isvm_G`(SVGFM)、band4=`isvm_T`(SVTM) |
| `bennu_mineral_abundance.tiff` | 矿物丰度图 | 6 类：phyllo / carbonate / sulfate / oxide / silicate / amorphous |
| `bennu_isvm_total.tiff` | ISVM 综合（加权和） | 单波段 |
| `preview_mineral_abundance.png`、`preview_isvm_4band.png` | 快速预览 | — |
| `outputs/exp9_raster_manifest.csv` | 自检清单（每波段 min/max/mean/σ、像元数、SHA-256） | — |

**栅格约定**

- 等经纬度格网，纬向 90°N→90°S，经度 0–360°E（Bennu 体固坐标），默认 0.5°/像元；
  Bennu 半径仅 245 m，故 1° 经度 ≈ 4.28 m，默认分辨率约 2.14 m/像元（`--res` 可调）。
- 写出 GeoTIFF 带地理变换与投影，波段描述同时写了中文语义（如 `isvm_M | SVMM …`），
  QGIS / ENVI / ArcGIS 打开即可看到通道含义。
- **每个像元的四通道值 = 科学价值评分**：候选区足迹内取该区评分，全球其余区域取
  全球背景基线（`data/global_baseline.csv`），两者按高斯足迹核平滑过渡。

**关于矿物丰度图（重要）**

`MINERAL_ABUNDANCE.md` 的检索结论是仓库内没有 Bennu 矿物丰度栅格，因此默认产物是
**Tier-1 代理丰度场**（文献锚点 + 固定随机种子，可复现），用于给 M 层提供空间区分度。
一旦拿到真实产品（PDS 的 OVIRS/OTES 反演结果或 MapCam 颜色指数镶嵌），直接传入即可：

```bash
python experiments/exp9_global_isvm_map.py --mineral-raster D:/data/ovirs_abundance.tif --res 0.25
```

站点坐标与足迹半径来自 `data/site_footprints.csv`（Nightingale / Osprey 有文献直接坐标，
Kingfisher / Sandpiper 取 Burke et al. 2021 测区中心），四个区当前统一用 25 m 足迹半径，
属可修改的设计假设。

### 叠加到真实 Bennu 底图（实验 10）

把真实的 Bennu DOM（正射影像）与 DEM 放进 `data/planetary/` 后，
`experiments/exp10_real_basemap.py` 会把四通道科学价值图重采样到**与真实底图完全相同的
格网**（shape 与 geotransform 逐位一致），输出：

| 文件 | 内容 |
|---|---|
| `bennu_dom_base.tiff` | 真实 DOM 底图（默认 10× 降采样 ≈0.49 m/px） |
| `bennu_isvm_4band_domgrid.tiff` | 四通道科学价值图（C/M/G/T），与底图逐像元对齐 |
| `preview_real_basemap.png` | 底图 / 高程 / 坡度 / 通道叠合预览 |
| `outputs/exp10_site_terrain.csv` | 四站点由 DEM 实测的坡度与高程统计 |

同时由 DEM 得到**真实坡度**，替代 `bennu_sites.csv` 中原本的假设值。
实测结论与文献交叉验证见 `data/planetary/README.md`——**四站点实测坡度
（9.4°–17.3°）普遍高于原假设值（5°–6°）**，其中 Osprey 达 17.3°，
在既有的 14° 安全门限下将不再通过，属需要决策的事项。

### 着陆点科学性价值与参数贡献度地图（实验 11）

`experiments/exp11_contribution_maps.py` 在真实 DOM 底图上表达着陆点科学价值，
并把 ISVM 分解到 10 个科学参数。**不使用坡度等工程门控**。

**全幅像元级填充（实装测试）**：候选区用站点级评分（与实验 1 完全一致），
区外用**合成占位数据**填满，使产品无空洞并可按"输入覆盖全球"的实装形态跑通：

- 连续参数：**评分 ≤ 50% × 该参数候选区最低评分**，乘上一个"共享基底场（65%）+ 本参数
  独立场（35%）"的多尺度平滑场，因此全幅既有成片的高值/低值区（对比明显），
  各参数之间又不完全相同；
- 分类参数：取评分最低的档位（分类评分是离散的，无法做到"趋近于 0"）；
- 因此**填充区评分恒 ≤ 候选区最小值**（实测填充区 ISVM ∈ [0.0453, 0.2631]，
  跨度 0.2178；候选区 ≥ 0.5772，二者相差 0.314），既不改变候选区结果，
  又保留了填充区内部的明显对比；
- ISVM 图与各参数图的色标改用**全幅数据范围**（原为固定 0.5–0.8），
  否则填充区会被压到色标下限、纹理不可见；候选区仍是全幅最高值区；
- 填充像元在 `bennu_fill_mask.tiff` 中标记为 0，可随时识别为合成数据；
- 填充场由固定种子 `20260920` 生成，可复现。

贡献度定义为加法分解：`c(p) = s(p) × w(指标) / n(指标内参数数)`，
10 个参数逐像元求和恒等于 ISVM（实测最大偏差 6e-8，为 float32 舍入）。

| 文件 | 内容 |
|---|---|
| `bennu_isvm_landing_value.tiff` | 着陆点科学性价值，单波段，区外 nodata |
| `bennu_param_contribution.tiff` | 10 波段，每个科学参数的**绝对贡献度** |
| `bennu_param_contribution_share.tiff` | 10 波段，**贡献度占比**（逐像元合计 100%） |
| `bennu_fill_mask.tiff` | 掩膜：1 = 候选区实测/文献值，0 = 合成填充 |
| `fig_landing_value_on_bennu.png` | 总图：真实影像 + 四站点 ISVM（含站点放大视图） |
| `fig_dom_global_sites.png` | **`fig_landing_value_on_bennu.png` 第 1 子图的等比放大版**（同一绘图函数 `_draw_global_panel()`，内容完全一致；地图区约 6.2×，3600×1875） |
| `fig_param_<参数>.png` ×10 | 每个参数的贡献度地图（全域 + 四站点放大，均以真实影像打底；色标为红蓝渐变 `coolwarm`，蓝=低、红=高，可用 `PARAM_CMAP` 更换） |
| `fig_contribution_bars.png` | 四站点贡献度堆叠条形图（各段之和 = ISVM） |
| `outputs/exp11_param_contribution.csv` | 站点 × 参数的评分、有效权重、贡献度、占比 |

栅格与真实底图格网完全一致（3141×1570，≈0.49 m/px），可直接在 GIS 中叠加。

**实装测试结果**（2026-09-20）：全幅 4.93 M 像元一次跑通，全流程 36.3 s
（填充+评分 1.6 s，栅格写出 9.7 s），全部栅格无 nodata。候选区数值与改动前**逐项相等**，
确认填充未影响结果；填充区纹理在成品图中清晰可见（全域面板填充区亮度标准差 16.3、
2.9 万种颜色，参数图标准差 23–25）。另有两个工程后果值得注意：填充区对比增强后
全幅 float32 栅格进一步增大（10 波段贡献度 96 MB、占比 135 MB，ISVM 单波段 12.9 MB，
合计 244 MB），其中占比栅格超过 GitHub 单文件 100 MB 上限，因此这两个大栅格已从
版本库移除（本地仍在，可 36 s 重跑生成），入库的是成品图、ISVM 单波段、掩膜与 CSV。

**图名约定**：以真实 DOM 为底图的图，图名统一用精简写法并在括号中标注数据源，
即「（DOM）」（例：`Bennu 全域（DOM）`、`着陆点科学性价值（DOM）`、
`参数贡献度：水合矿物吸收带深（DOM）`）。承载数值或单位的括号（如 `(m)`、`(°)`、
`(18.2%)`）保留，不替换。

## 注意事项

- 数据集中标注为 `estimate` / `assumed` 的参数为文献间接推断或合理假设，
  使用时请核对 `data/README.md`；
- 缺失参数在评分时取中性值 0.5（连续型）并记录警告。
- 实验 9 的矿物丰度为**代理场**而非实测，正式发表前须替换为 OVIRS/OTES 产品。
