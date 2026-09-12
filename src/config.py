# -*- coding: utf-8 -*-
"""
科学价值评价指标体系与参数定义（文献驱动）。

文献依据：
- Nakamura-Messenger et al., Strategy for Ranking the Science Value of the Surface of
  Asteroid 101955 Bennu for Sample Site Selection for OSIRIS-REx（四张科学价值图）
- Enos et al., 2020, OSIRIS-REx's Search for a Sample Site: Nightingale / Osprey
- Yabuta et al., 2019, Landing Site Selection for Hayabusa2（龙宫七候选区）
- Walsh et al., 2022, Assessing the Sampleability of Bennu's Surface
"""

INDICATORS = ["C", "M", "G", "T"]

INDICATOR_LABELS = {
    "C": "化学组成价值 SVCCM",
    "M": "矿物学价值 SVMM",
    "G": "地质特征/新鲜度价值 SVGFM",
    "T": "温度/热环境价值 SVTM",
}

# 冗余参数清单（经三组数据方差/相关性审计确认，见 experiments/exp7_redundancy_audit.py）
# 删除/合并策略：常数参数直接删除；完全同向参数保留代表参数；六类矿物检出合并为
# mineral_classes_detected 一个参数。定义保留在 PARAMETERS 中仅供追溯，不再参与评分。
REDUNDANT_SCIENCE = {
    "volatile_abundance": "与 organic_abundance 在现有数据中完全同向，合并为有机物/挥发物综合丰度",
    "organic_silicate_ratio": "与 organic_abundance 在现有数据中完全同向，合并",
    "mineral_phyllo": "三组数据均恒为 1.0，无区分度；并入 mineral_classes_detected",
    "mineral_carbonate": "并入 mineral_classes_detected（Nightingale 差异由类别数体现）",
    "mineral_sulfate": "三组数据均恒为 0.5；并入 mineral_classes_detected",
    "mineral_oxide": "三组数据均恒为 0.5；并入 mineral_classes_detected",
    "mineral_silicate": "三组数据均恒为 1.0；并入 mineral_classes_detected",
    "mineral_amorphous": "三组数据均恒为 0.5；并入 mineral_classes_detected",
    "amorphous_fraction": "三组数据均恒为 0.5（常数），无区分度，删除",
    "brittle_deformation": "与 crater_freshness 在现有数据中近似共线，合并入地质新鲜度概念",
    "geologic_diversity": "当前数据区分度弱且为假设值，删除",
}

PARAMETERS = {
    # ---------------- C 化学组成 ----------------
    "organic_abundance": {
        "indicator": "C", "label": "有机物检出与丰度", "kind": "continuous",
        "direction": "higher_better", "lo": 0.0, "hi": 1.0,
        "literature": "SVCCM 参数1：有机物/挥发物的检测与丰度（Nakamura-Messenger et al.）",
    },
    "volatile_abundance": {
        "indicator": "C", "label": "挥发物丰度", "kind": "continuous",
        "direction": "higher_better", "lo": 0.0, "hi": 1.0,
        "literature": "SVCCM 参数1：有机物/挥发物的检测与丰度（Nakamura-Messenger et al.）",
    },
    "organic_silicate_ratio": {
        "indicator": "C", "label": "有机/硅酸盐比", "kind": "continuous",
        "direction": "higher_better", "lo": 0.0, "hi": 1.0,
        "literature": "SVCCM 参数2：判断有机物是否均匀分布或存在富集区（Nakamura-Messenger et al.）",
    },
    "ch2_ch3": {
        "indicator": "C", "label": "CH2/CH3 比", "kind": "continuous",
        "direction": "lower_better", "lo": 0.5, "hi": 3.0,
        "literature": "SVCCM 参数3：CH2/CH3 越低表示更长链脂肪烃、更原始（Sandford et al. 1991）",
    },
    "carbon_content_pct": {
        "indicator": "C", "label": "碳含量 (%)", "kind": "continuous",
        "direction": "higher_better", "lo": 0.0, "hi": 10.0,
        "literature": "龙宫 R(0.55)≈0.018–0.019 对应 >3% 碳（Yabuta et al. 2019）",
    },
    # ---------------- M 矿物学 ----------------
    "mineral_phyllo": {
        "indicator": "M", "label": "层状硅酸盐检出", "kind": "continuous",
        "direction": "higher_better", "lo": 0.0, "hi": 1.0,
        "literature": "SVMM 矿物优先级第1位：phyllosilicates（Nakamura-Messenger et al.）",
    },
    "mineral_carbonate": {
        "indicator": "M", "label": "碳酸盐检出", "kind": "continuous",
        "direction": "higher_better", "lo": 0.0, "hi": 1.0,
        "literature": "SVMM 矿物优先级第2位：carbonates",
    },
    "mineral_sulfate": {
        "indicator": "M", "label": "硫酸盐检出", "kind": "continuous",
        "direction": "higher_better", "lo": 0.0, "hi": 1.0,
        "literature": "SVMM 矿物优先级第3位：sulfates",
    },
    "mineral_oxide": {
        "indicator": "M", "label": "氧化物检出", "kind": "continuous",
        "direction": "higher_better", "lo": 0.0, "hi": 1.0,
        "literature": "SVMM 矿物优先级第4位：oxides",
    },
    "mineral_silicate": {
        "indicator": "M", "label": "硅酸盐检出", "kind": "continuous",
        "direction": "higher_better", "lo": 0.0, "hi": 1.0,
        "literature": "SVMM 矿物优先级第5位：silicates",
    },
    "mineral_amorphous": {
        "indicator": "M", "label": "非晶质检出", "kind": "continuous",
        "direction": "higher_better", "lo": 0.0, "hi": 1.0,
        "literature": "SVMM 矿物优先级第6位：amorphous；非晶质硅酸盐与最原始天体相关",
    },
    "amorphous_fraction": {
        "indicator": "M", "label": "非晶质硅酸盐比例", "kind": "continuous",
        "direction": "higher_better", "lo": 0.0, "hi": 1.0,
        "literature": "晶质/非晶质比反映表面原始性；非晶质越多越原始（Nakamura-Messenger et al.）",
    },
    "hydration_depth": {
        "indicator": "M", "label": "水合矿物吸收带深", "kind": "continuous",
        "direction": "higher_better", "lo": 0.0, "hi": 1.0,
        "literature": "贝努 3 µm/龙宫 2.7 µm 带深；水合与挥发分保存（Enos et al.; Yabuta et al.）",
    },
    # ---------------- G 地质特征/新鲜度 ----------------
    "space_weathering_proxy": {
        "indicator": "G", "label": "空间风化代理（光谱红化）", "kind": "continuous",
        "direction": "lower_better", "lo": 0.0, "hi": 1.0,
        "literature": "SVGFM：空间风化低于/高于全球均值；龙宫 L 区偏蓝、M 区偏红（Yabuta et al. 2019）",
    },
    "crater_freshness": {
        "indicator": "G", "label": "陨石坑新鲜度", "kind": "categorical",
        "categories": {"recent": 1.0, "paleo": 0.6, "none": 0.3},
        "literature": "SVGFM：recent 陨石坑暴露新鲜物质；Nightingale 的 Hokioi 为年轻陨石坑（Enos et al.）",
    },
    "activity_level": {
        "indicator": "G", "label": "活动/羽流证据", "kind": "categorical",
        "categories": {"active": 1.0, "recent": 0.75, "paleo": 0.5, "none": 0.2},
        "literature": "SVGFM：active > recent > paleo（Nakamura-Messenger et al.）",
    },
    "psfd_class": {
        "indicator": "G", "label": "PSFD 类别", "kind": "categorical",
        "categories": {"fine": 1.0, "mix": 0.65, "coarse": 0.35},
        "literature": "SVGFM：PSFD 分 mix/coarse/fine 三类，记录表面过程（Nakamura-Messenger et al.）",
    },
    "brittle_deformation": {
        "indicator": "G", "label": "脆性变形构造", "kind": "binary",
        "literature": "SVGFM 子类：brittle deformation（Nakamura-Messenger et al.）",
    },
    "geologic_diversity": {
        "indicator": "G", "label": "地形/地质多样性", "kind": "continuous",
        "direction": "higher_better", "lo": 0.0, "hi": 1.0,
        "literature": "龙宫科学评价准则2：多样地形地质可研究多种演化物质（Yabuta et al. 2019）",
    },
    # ---------------- T 温度/热环境 ----------------
    "temp_rel_global": {
        "indicator": "T", "label": "温度相对全球均值", "kind": "categorical",
        "categories": {"colder": 1.0, "average": 0.5, "hotter": 0.0},
        "literature": "SVTM：寻找最冷区域，利于挥发物/有机物保存（Nakamura-Messenger et al.）",
    },
    "tmax_k": {
        "indicator": None, "label": "最高温度 Tmax (K)", "kind": "continuous",
        "direction": "lower_better", "lo": 250.0, "hi": 350.0,
        "literature": "龙宫 Tmax<350 K 为安全条件；温度低利于保存（Yabuta et al. 2019）。"
                      "科学评分中由 temp_rel_global 代表，本参数仅作工程安全门控输入",
    },
    "mineral_classes_detected": {
        "indicator": "M", "label": "关键矿物类别数", "kind": "continuous",
        "direction": "higher_better", "lo": 0.0, "hi": 6.0,
        "literature": "由 SVMM 六类矿物检出（层状硅酸盐/碳酸盐/硫酸盐/氧化物/硅酸盐/非晶质）合并而成，"
                      "反映矿物多样性（Nakamura-Messenger et al.）",
    },
    # ---------------- 工程层（不进 ISVM） ----------------
    "slope_deg": {
        "indicator": None, "label": "坡度 (deg)", "kind": "continuous",
        "direction": "lower_better", "lo": 0.0, "hi": 14.0,
        "literature": "贝努安全约束：TAGSAM 接触面角 >14° 概率 <1%（Enos et al. 2020）",
    },
    "boulder_coverage": {
        "indicator": None, "label": "块石覆盖率", "kind": "continuous",
        "direction": "lower_better", "lo": 0.0, "hi": 0.25,
        "literature": "龙宫安全区块石覆盖率 10.7–16.8%，其余 17.2–20.3%（Yabuta et al. 2019）",
    },
    "roughness": {
        "indicator": None, "label": "表面粗糙度", "kind": "continuous",
        "direction": "lower_better", "lo": 0.0, "hi": 1.0,
        "literature": "工程安全性辅助参数（估计）",
    },
    "frac_le2cm": {
        "indicator": None, "label": "≤2 cm 无障碍细粒占比", "kind": "continuous",
        "direction": "higher_better", "lo": 0.0, "hi": 1.0,
        "literature": "TAGSAM 可吸入粒径 ≤2 cm；Walsh et al. 2022 可采样性评分核心",
    },
    "min_particle_cm": {
        "indicator": None, "label": "最小粒径 (cm)", "kind": "continuous",
        "direction": "lower_better", "lo": 0.0, "hi": 5.0,
        "literature": "可采样性函数输入之一：最小粒径（Walsh et al. 2022）",
    },
    "max_particle_cm": {
        "indicator": None, "label": "最大粒径 (cm)", "kind": "continuous",
        "direction": "lower_better", "lo": 0.0, "hi": 10.0,
        "literature": "可采样性函数输入之一：最大粒径；≤2 cm 为宜（Walsh et al. 2022）",
    },
    "tagsam_tilt_deg": {
        "indicator": None, "label": "TAGSAM 倾角 (deg)", "kind": "continuous",
        "direction": "lower_better", "lo": 0.0, "hi": 14.0,
        "literature": "可采样性函数输入之一：TAGSAM 头倾角（Walsh et al. 2022）",
    },
    "target_radius_m": {
        "indicator": None, "label": "目标区半径 (m)", "kind": "continuous",
        "direction": "higher_better", "lo": 0.0, "hi": 25.0,
        "literature": "贝努原始计划 25 m 安全区，实际收缩至 5–8 m（Enos et al. 2020）",
    },
    "nav_error_m": {
        "indicator": None, "label": "导航误差 (m)", "kind": "continuous",
        "direction": "lower_better", "lo": 0.0, "hi": 5.0,
        "literature": "龙宫 L08-E1 导航制导误差上限约 2.7 m（JAXA）",
    },
}

# 矿物检出优先级（文献顺序）
MINERAL_PRIORITY = {
    "mineral_phyllo": 1.0,
    "mineral_carbonate": 0.9,
    "mineral_sulfate": 0.8,
    "mineral_oxide": 0.7,
    "mineral_silicate": 0.6,
    "mineral_amorphous": 0.5,
}

SCIENCE_PARAMETERS = [
    "organic_abundance", "ch2_ch3", "carbon_content_pct",
    "mineral_classes_detected", "hydration_depth",
    "space_weathering_proxy", "crater_freshness", "activity_level", "psfd_class",
    "temp_rel_global",
]

# 各指标实际参与评分的参数（冗余精简后）
INDICATOR_PARAMETERS = {
    "C": ["organic_abundance", "ch2_ch3", "carbon_content_pct"],
    "M": ["mineral_classes_detected", "hydration_depth"],
    "G": ["space_weathering_proxy", "crater_freshness", "activity_level", "psfd_class"],
    "T": ["temp_rel_global"],
}
