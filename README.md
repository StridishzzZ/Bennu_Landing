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
│   └── README.md                    # 数据来源与置信度说明
├── src/                             # 评分流水线
│   ├── config.py
│   ├── scoring.py                   # 单参数评分函数
│   ├── indicators.py                # 四图指标聚合
│   ├── weights.py                   # 三种权重方案
│   ├── pipeline.py                  # ISVM 与工程层门控
│   └── analysis.py                  # 敏感性 / 蒙特卡洛
├── experiments/                     # 实验脚本
│   ├── common.py
│   ├── exp1_bennu_validation.py
│   ├── exp2_ryugu_validation.py
│   ├── exp3_weight_comparison.py
│   ├── exp4_sensitivity.py
│   ├── exp5_monte_carlo.py
│   ├── exp6_synthetic_stress.py
│   └── exp7_redundancy_audit.py
└── outputs/                         # 实验结果（CSV）
```

## 运行

```bash
python experiments/exp1_bennu_validation.py
python experiments/exp2_ryugu_validation.py
python experiments/exp3_weight_comparison.py
python experiments/exp4_sensitivity.py
python experiments/exp5_monte_carlo.py
python experiments/exp6_synthetic_stress.py
```

结果写入 `outputs/`。

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

## 注意事项

- 数据集中标注为 `estimate` / `assumed` 的参数为文献间接推断或合理假设，
  使用时请核对 `data/README.md`；
- 缺失参数在评分时取中性值 0.5（连续型）并记录警告。
