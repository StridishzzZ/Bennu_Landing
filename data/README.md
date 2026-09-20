# 数据来源与置信度说明

所有数据集均为"直接给定参数"模式，值为文献直接引用或由文献间接推断/假设。
评分前请核对以下说明；标为 estimate/assumed 的项在使用真实任务数据时应替换。

## 精简说明（2026-08-26）

经冗余审计（`experiments/exp7_redundancy_audit.py`），科学参数由 21 个精简为 10 个：

- **删除（常数，无区分度）**：mineral_phyllo、mineral_sulfate、mineral_oxide、
  mineral_silicate、mineral_amorphous、amorphous_fraction、geologic_diversity；
- **合并**：volatile_abundance、organic_silicate_ratio → organic_abundance；
  六类矿物检出 → mineral_classes_detected（0–6）；brittle_deformation → crater_freshness；
  tmax_k 降为纯工程门控参数；
- 精简前后综合分排序全部一致（验收标准通过），完整参数集归档于 `data/full/`。

## 当前列结构（精简后）

科学参数（10）：organic_abundance、ch2_ch3、carbon_content_pct、
mineral_classes_detected、hydration_depth、space_weathering_proxy、
crater_freshness、activity_level、psfd_class、temp_rel_global；

工程参数（9）：tmax_k、slope_deg、boulder_coverage、roughness、frac_le2cm、
max_particle_cm、tagsam_tilt_deg、target_radius_m、nav_error_m；

元数据：site、source、confidence。

## 参数含义

- 连续型参数按 `src/config.py` 中的锚点线性标准化：
  ch2_ch3（0.5–3.0，越低越原始）、carbon_content_pct（0–10%，龙宫 >3%）、
  mineral_classes_detected（0–6）、tmax_k（250–350 K，龙宫 <350 为安全）、
  slope_deg（0–14°，贝努倾角约束）、boulder_coverage（0–25%，
  龙宫安全区 10.7–16.8%）、target_radius_m（0–25 m，贝努实际 5–8 m）、
  nav_error_m（0–5 m，龙宫 L08-E1 上限 2.7 m）；
- 分类型参数：crater_freshness ∈ {recent, paleo, none}；
  activity_level ∈ {active, recent, paleo, none}；psfd_class ∈ {fine, mix, coarse}；
  temp_rel_global ∈ {colder, average, hotter}。

## bennu_sites.csv

- Nightingale：年轻陨石坑（Hokioi）、北部低温、细粒可采样物质丰富、水合/挥发分保存有利
  （Enos et al. 2020; Walsh et al. 2019）——置信度高；
- Osprey：赤道、可采样性较低但危险物更少、PSFD 近几何饱和（Enos et al. 2020; PSFD 论文）
  ——置信度中；
- Kingfisher / Sandpiper：文献只确认"位于陨石坑内、含相对丰富可采样物质"（PSFD 论文），
  其余参数为假设值——置信度低；
- Nightingale 的 mineral_classes_detected=4.5（含碳酸盐），其余 4.0，
  依据"亮色物质可能含碳酸盐/外来玄武岩"（夜莺地质论文）。

## ryugu_sites.csv

- 块石覆盖率：L07/L08/M04 落在 10.7–16.8%，其余落在 17.2–20.3%
  （Yabuta et al. 2019 区间），具体数值为区间内代表性取值——中等置信；
- 空间风化代理：L 区 0.30–0.35、M 区 0.55–0.60（据 L/M 偏蓝/偏红分区赋值）；
- 碳含量：全域 3.5%，M01 4.0%（R(0.55) 略低，Yabuta et al. 2019）；
- 水合矿物、有机物在龙宫全域近似均匀，故各候选区取相同中性值；
- 坡度/粗糙度/可采样细粒占比/目标半径/导航误差为工程层假设值，
  需以 Kikuchi et al. 2020/2022 原文复核。

## ryugu_subsites.csv

- L08-E1：首次 touchdown 点，目标区约 6 m 直径、导航误差上限约 2.7 m、细粒丰富
  （JAXA 官方说明；Morota et al. 2020）——置信度高；
- L08-B：SCI 人工陨石坑后的二次采样点，靠近喷射物富集区、科学价值高、工程约束更严
  （Kikuchi et al. 2022）——置信度中。

## site_footprints.csv（实验 9 新增）

贝努四候选区的文献坐标与足迹半径，用于把评分结果栅格化到全球地图。

- 中心坐标取自 `articles/3-Geologic Context ... sample sites.pdf`
  （Burke et al. 2021, Remote Sens. 13, 1315）中逐阶段的 ROI 中心表：
  Nightingale 56.05°N/42.05°E（ReconC）、Osprey 11.62°N/88.63°E（ReconC）、
  Kingfisher 11.49°N/55.49°E（ReconA）、Sandpiper −46.98°N/321.34°E（ReconA）。
  Nightingale 与 Osprey 另有 `articles/1463.pdf`（Enos et al. 2020）独立给出的
  56.05°/42.05° 与 11.62°/88.62°，两者一致——置信度高；
- Kingfisher / Sandpiper 只有 ReconA 测区（无 ReconC 精细收敛），置信度中；
- `footprint_radius_m` 统一取 25 m，是**设计假设**而非文献实测值
  （量级参照 20–25 m 级采样椭圆），需要更严格时请按各站点 2σ 椭圆改正后重跑实验 9。
- 注意 Sandpiper 经度为 321.34°E（即 38.66°W），不要误写成 43°E。

## global_baseline.csv（实验 9 新增）

全球背景基线参数（一行），供实验 9 在地图非候选区部分填充背景值。

- 依据：Lauretta et al. 2019 全球表征、Hamilton et al. 2019 全球水合近均匀、
  Simon et al. 2020 含碳物质全球分布；数值为全球平均量级的设计基线，**非实测栅格**；
- 该行经同一套评分函数得到背景指示分（当前 priority 权重下 C=0.613、M=0.683、
  G=0.488、T=0.500）；
- 若后续获得全球栅格产品，建议用真实栅格的全球中位数替换本行数值。

## 使用建议

- 正式实验前，请把 confidence=low 的参数替换为你自己的测量/估算值；
- 敏感性分析（实验 4）和蒙特卡洛（实验 5）覆盖这些不确定项；
- 如需恢复完整参数集，参考 `data/full/` 与 `src/config.py` 中的定义。
