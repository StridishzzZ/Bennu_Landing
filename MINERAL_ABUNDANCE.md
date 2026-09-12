# Bennu 矿物丰度图：文献检索结论、可支撑数据与模拟难度评估

检索对象：`articles/` 下全部 14 篇 PDF（含 4 篇 Ryugu / Hayabusa2 与方法类文献），
外加 `data/`、`data/full/`、`src/config.py` 与 `experiments/exp7_redundancy_audit.py` 的现有参数集。
检索方式：全文文本抽取后做矿物学术语审计（abundance / phyllosilicate / magnetite / carbonate /
sulfide / hydrated / band depth / 2.7 um / OVIRS / OTES / thermal inertia 等），
并对含图页做图注（figure caption）核对。

聚合后的数据清单见 `data/mineral_abundance_support.csv`。

## 0. 结论摘要

1. 文件夹内没有 Bennu 的矿物丰度图。没有 OVIRS / OTES 派生的矿物分布图、水合带深图或任何以
   矿物丰度为量的栅格产品。文献集中只在引言/讨论里引用这类成果
   （Hamilton et al. 2019；Kaplan et al. 2020；Simon et al. 2020；Rozitis et al. 2020）。

2. 文件夹里最接近矿物图的三类产品：
   - 反照率图（global normal albedo 6.25 cm/px；OLA 1064 nm 相对反照率 0.8 m）——物性代理，不是矿物；
   - MapCam 四色指数图（700 nm 吸收指数被明确标注为层状硅酸盐代理，另有 b'/v、v/x、turn-off point）——矿物代理由此可得，最接近可用；
   - 岩块类型图 Type A-D（反照率 + 形态 + 粗糙度，389 块，Nightingale 区）——矿物相分类代理。

3. 矿物图在文献中的设计早已存在，但只是方案：Nakamura-Messenger et al. 的 SVMM
   （Science Value Mineralogy Map）给出了矿物优先级（phyllosilicates > carbonates > sulfates >
   oxides > silicates > amorphous）与 crystalline/amorphous 比，是本项目 `src/config.py`
   中 M 层参数的来源；但该文是预到达（pre-encounter）策略文，不含任何实际矿物分布数据。

4. 模拟矿物丰度图可行，但分三档，难度差异极大。代理图（Tier-1）1-2 周可做，
   半定量（Tier-2）1-3 个月，定量丰度（Tier-3）属研究级课题，且现有数据在原理上无法支撑米级定量丰度。

5. 与本项目的直接接口：`exp7` 冗余审计已证实现有 6 类矿物参数在贝努四区全为常数
   （phyllo=1.0、silicate=1.0、sulfate/oxide/amorphous=0.5），被合并为 `mineral_classes_detected`；
   而敏感性分析中 `mineral_classes_detected` 与 `hydration_depth` 恰好是影响 ISVM 最大的两个参数。
   换句话说：M 层目前几乎没有区分度，而它的区分度只能靠这样一张（模拟的）矿物丰度图来提供。

## 1. 文献检索：为什么找不到矿物丰度图

### 1.1 语料构成决定了主题偏向

`articles/` 的 14 篇文献按主题可分为三类，没有一篇是矿物学/光谱学专题论文：

| 主题 | 文献 |
|---|---|
| 地形/影像测绘（Bennu） | Bennett et al. 2021（全球底图）；Barnouin et al. 2020（DTM 方法）；Olds et al. 2022（NFT 用 DTM）；DellaGiustina et al. 2018（影像制图挑战） |
| 地质/采样区表征（Bennu） | Barnouin et al. 2022（Nightingale 地质背景）；Jawin et al. 2023（岩块多样性）；Burke et al. 2021（粒径分布）；Walsh et al. 2022（可采样性）；Bierhaus et al. 2021（TAGSAM 搅动） |
| 选址决策（Bennu / Ryugu） | Enos et al. 2020（LPSC #1463）；Nakamura-Messenger et al.（科学价值图策略）；Yabuta et al. 2019（LPSC #2304）；Kikuchi et al. 2022；Lorda et al. 2020 |

因此本语料回答的是"地表长什么样、能不能落、能不能采"，而不是"表面由什么矿物组成、各占多少"。

### 1.2 术语审计结果（证据）

对 14 篇全文（约 11300 行）做关键词检索：

| 关键词 | 命中情况 |
|---|---|
| `abundance map` / `mineral map` | 0 次 |
| `band depth` | 仅 2 篇文献中 5 处：Ryugu 摘要（Yabuta et al. 2019，NIRS3 的 3 um 带深图）；DellaGiustina et al. 2018 中描述 MapCam 700 nm 吸收指数（相对带深 RBD） |
| `phyllosilicate` | 仅出现在（a）SVMM 的矿物优先级清单，（b）Barnouin et al. 2022 对 Bennu 整体成分的引述，（c）MapCam 700 nm 指数解释 |
| `magnetite` / `sulfide` | 仅 Barnouin et al. 2022 引言中一处，描述 Bennu 与含水蚀变 CM/CI 陨石的类比 |
| `OVIRS` / `OTES` | 出现 10 余次，全部为仪器介绍、引用或"数据用于温度/矿物学研究"的方法学陈述，没有给出派生图件 |
| `thermal inertia` | 出现在 Walsh et al. 2022（给出数值）、Enos et al. 2020（说明热惯量未能提供粒径判别力） |

### 1.3 文件夹内最接近矿物图的图件（逐一核实）

| 图件 | 出处 | 实测性质 | 能否当矿物丰度图用 |
|---|---|---|---|
| global normal albedo map（6.25 cm/px，median 0.046） | Burke et al. 2021 图 11；Jawin et al. 2023 图 3、5 | 反照率（I/F）栅格 | 否。反照率受空间风化、粒径、孔隙度共同控制 |
| OLA albedo（1064 nm，0.8 m GSD） | Barnouin et al. 2022 图 2(d) | 激光峰值幅度归一化后的相对反照率 | 否，但可补足影像阴影区 |
| geologic map / unit map（Rugged Unit / Smooth Unit） | Barnouin et al. 2022；Jawin et al. 2022（被引） | 地质单元边界矢量 | 否，但可作空间先验 |
| Type A-D 岩块分类图 | Jawin et al. 2023 图 9 | 反照率 + 形态 + 粗糙度的分类 | 可作相（facies）代理 |
| MapCam 四色指数图（含 700 nm 吸收指数） | DellaGiustina et al. 2018 | 层状硅酸盐代理指数（原文明确指向 Vilas 1994） | 最接近，但仍是指数而不是丰度 |
| Ryugu NIRS3 2.7 um 带深/带心图 | Yabuta et al. 2019 图 3 | 真正意义上的矿物带深图 | 方法模板可用，但不是 Bennu |

结论明确：需要模拟，且模拟必须依赖文件夹之外的光谱数据。

## 2. 可支撑模拟的数据汇总

完整表格见 `data/mineral_abundance_support.csv`（18 条记录，含分辨率、覆盖范围、数值、出处、
置信度与模拟中的用途）。按是否已在本文件夹分三层。

### 2.1 A 层：文件夹内可直接使用的定量数据

| 数据 | 关键数值 | 出处 | 对模拟的作用 |
|---|---|---|---|
| 全球底图 | 5 cm/px，平均空间精度 ~30 cm，Minnaert 光度归一化 | Bennett et al. 2021 | 所有代理层的共同几何底图 |
| 反照率 | median 0.046；岩块均值 0.050±0.005；双峰中心 0.046 / 0.055 | Golish et al. 2021；Jawin et al. 2023 | 暗/亮物质二分，作为矿物端的边界条件 |
| 岩块类型 | A 0.045±0.002、B 0.047±0.002、C 0.055±0.004、D 0.054±0.004；ν(10 cm)=7.5 / 5.2 / 3.9 / 4.6 cm | Jawin et al. 2023 | 相分类赋值表，可直接映射为矿物组合的离散场 |
| 亮斑（碳酸盐） | 亮斑反照率 10-19%，宿主岩块 6-7.6% | Kaplan et al. 2020（经 Jawin et al. 2023 量化） | 碳酸盐层的唯一空间指示 |
| 热惯量 | 到达后 350±20；到达前盘面平均 310±70 J m-2 K-1 s-1/2 | Walsh et al. 2022 | 粒径/孔隙度先验，约束辐射传输模型 |
| 粒径分布 | 幂律指数 -3.0±0.2 至 -2.3±0.1；Nightingale -2.2±0.1，Osprey -2.7±0.6 | Burke et al. 2021 | 决定混合尺度与端元配比 |
| 力学性质 | 内聚力 ~0.6 Pa，摩擦角自 32 度起 | Barnouin et al. 2022 | 合成场的物质混合/迁移约束 |
| 地质单元与地貌 | Rugged/Smooth Unit；0 度/90 度纵向脊；Hokioi 20 m 坑、外缘 120 m 坑，采样物质可来自 10 m 深 | Barnouin et al. 2022 | 空间相关性/自相关的结构模板 |

### 2.2 B 层：文件夹内仅作引用的目标基准

- Hamilton et al. 2019（Nat. Astron. 3, 332）：OVIRS 2.7 um 水合带，全球水合矿物分布。
  这是任何 phyllosilicate 层必须对照的基准（本文件夹只引用，未含图）。
- Kaplan et al. 2020（Science）：碳酸盐亮脉，仅通过 Jawin et al. 2023 的二次引用获得数值。
- Simon et al. 2020（Science）：含碳物质/有机物分布。
- Rozitis et al. 2020：岩块热惯量与孔隙度二分（低反射率+高孔隙度 对 高反射率+低孔隙度）。
- DellaGiustina et al. 2020：MapCam 全球颜色/反照率与两类岩块群。

### 2.3 C 层：必须从文件夹外获取的数据（模拟的前置条件）

| 数据 | 获取途径 | 说明 |
|---|---|---|
| OVIRS 光谱立方体（0.4-4.3 um） | NASA PDS（OSIRIS-REx OVIRS bundle） | 唯一能给出 2.7 um 水合带、碳酸盐/有机物吸收的定量源 |
| OTES 热红外光谱（5-50 um） | NASA PDS（OTES bundle） | 与 VNIR 独立的矿物学约束 |
| 全球反照率镶嵌 | PDS Derived Image Processing Bundle，Global Albedo Mosaic 6.25 cm | 入口已记录在 `articles/采样区分析总结.docx` |
| MapCam 颜色镶嵌 | PDS / OCAMS bundle | 高分辨率矿物代理层 |
| 形状模型 + SPICE 核 | PDS / 任务官网（75 cm 形状模型） | 投影与几何（入射/发射角）校正 |
| 实验室端元光谱 | RELAB / 陨石（Murchison、Orgueil、Tagish Lake） | 解混所需的端元库 |

注：PDS 具体 bundle 名与 URN 会随版本更新，本报告不写死链接；请按 `data/README.md` 与
`articles/采样区分析总结.docx` 中已记录的入口检索并核对版本。

## 3. 模拟方案：三档递进

### Tier-1 代理相图（Proxy facies map）——推荐作为当前项目首选

- 输入：反照率镶嵌（6.25 cm）+ MapCam 颜色指数 + Type A-D 岩块图 + 地质单元图。
- 方法：把反照率 × 颜色指数 × 单元 × 岩块类型做无监督分类（k-means / 高斯混合），得到
  4-6 个矿物相类别；每个类别按文献赋值一个离散矿物组合（例：暗粗糙相 = 富层状硅酸盐 +
  磁铁矿 + 有机物；亮平滑相 = 富碳酸盐 + 低孔隙度；亮斑 = 碳酸盐端元；外源物质 = 玄武质）。
- 输出：与底图同分辨率的类别图（categorical），附各类占比统计。
- 验证：与 Hamilton et al. 2019 的全球水合近均匀结论做一致性检查；与返回样品的暗物质主导
  预测对照（Jawin et al. 2023 在 TAGSAM 接触区内识别出 72 个可测颗粒，其中暗色 Type A/B
  46 个、亮色 Type C/D 26 个，暗色颗粒可见面积约为亮色的 2 倍；该区域内约 63% 面积被
  阴影遮挡、无法制图，可作为掩膜与不确定度参照）。
- 难度：低-中（2/5）。

### Tier-2 半定量丰度（2-5 端元，百米级）

- 输入：OVIRS / OTES 光谱立方体 + 温度与光度校正 + 端元库。
- 方法：Hapke 或 Shkuratov 辐射传输 + 线性/非线性光谱解混（NNLS、MESMA），端元取 CM/CI 类比
  （蛇纹石、皂石、方解石、磁铁矿、陨硫铁、有机物、非晶硅酸盐）。
- 输出：全球格网上的体积/重量分数估计（带不确定度）。
- 验证：与 OTES 独立反演互检；与样品模态矿物学（唯一真值）做量级对照。
- 难度：高（4/5）。主要不确定性来自端元光学常数、粒径与孔隙度、空间风化层的暗化效应。

### Tier-3 米级定量丰度（科学上不可达）

- 物理上限：OVIRS/OTES 是点式光谱仪，足印为数十至数百米量级，而影像可达 5-6 cm。
  DellaGiustina et al. 2018 明确写到 MapCam 颜色指数图的空间尺度
  "significantly finer than OVIRS"，并采用 x-band 图像 + OVIRS 860-1064 nm 比值的方式
  外推 1064 nm 底图，这正是图像锐化（sharpening）思想；但锐化的是颜色/反射率，不是矿物丰度。
- 因此米级矿物丰度图只能由代理回归得到，其真值无法验证：全任务只有一份百克级的单点样品
  （Nightingale）。
- 难度：极高（5/5），属研究课题，不建议作为工程交付目标。

## 4. 模拟难度评估

### 4.1 分维度打分（1 = 容易，5 = 很难）

| 维度 | 分数 | 主要障碍 |
|---|---|---|
| 数据可获得性 | 3 | 文件夹内无光谱数据，须从 PDS 下载 OVIRS/OTES 立方体与 SPICE 核；体积大、需网络 |
| 数据预处理 | 4 | 光度（BRDF）校正、热辐射去除、几何（入射/发射角）配准，需 ISIS3/SPICE 工具链 |
| 物理模型 | 4 | 碳质物质的光学常数与端元库不完整；Hapke 参数（粒径、孔隙度、粗糙度）耦合强 |
| 空间分辨率匹配 | 5 | 足印（数十至数百米）与影像（5-6 cm）相差 3-4 个数量级，只能靠代理或锐化 |
| 验证/真值 | 5 | 无独立地面真值；返回样品为单点、非原位，且为破坏性分析 |
| 与现有 ISVM 接口 | 2 | 只需产出一张 0-1 归一化栅格或分区统计值即可替换 `mineral_classes_detected` |

### 4.2 综合判断

| 档位 | 综合难度 | 估计工作量 | 前置条件 |
|---|---|---|---|
| Tier-1 代理相图 | 2/5，可行 | 1-2 周 | 仅需已公开的反照率镶嵌 + MapCam 颜色镶嵌 |
| Tier-2 半定量丰度 | 4/5，困难 | 1-3 个月 | OVIRS/OTES 数据 + 端元库 + 辐射传输实现 |
| Tier-3 定量米级丰度 | 5/5，不建议 | 研究项目级 | 需要新的观测或样品原位分析技术 |

### 4.3 三个最容易踩的坑

1. 把反照率当矿物丰度。Bennu 的反照率差异主要来自空间风化与孔隙度（DellaGiustina et al. 2020
   的两类岩块：低反射率+粗糙+高孔隙度 对 高反射率+平滑+低孔隙度），直接线性映射会系统性
   高估亮物质的碳酸盐含量。
2. 忽略 OVIRS 分辨率上限。任何声称米级矿物丰度的结果都必须是代理外推，需显式声明假设。
3. 忽视全局近均匀性。Walsh et al. 2022 及相关光谱研究均指出 Bennu 全域水合矿物近似均匀；
   若模拟只产出"处处相同"的图，则对选址排序毫无贡献。模拟的价值在于刻画小尺度异质性
   （亮斑、外源物质、新鲜坑内物质），而不是全球平均。

## 5. 与本项目 ISVM / SVMM 的接口建议

现状（见 `RESULTS.md` 实验 3、4、7）：

- 6 类矿物 flag 在贝努四区为常数，已合并为 `mineral_classes_detected`；
- 敏感性分析中，对 ISVM 影响最大的两个参数正是 `hydration_depth` 与
  `mineral_classes_detected`（±0.011-0.012）；
- `data/README.md` 注明 Nightingale 的 `mineral_classes_detected=4.5` 仅依据
  "亮色物质可能含碳酸盐/外来玄武岩"这一条定性判断。

建议的最小改动路径：

1. 先用 Tier-1 代理相图，把 Nightingale / Osprey / Kingfisher / Sandpiper 四个 2σ 椭圆内的
   相占比统计（例如亮斑像素比例、Type C/D 占比）导出为 4 个数值；
2. 用这些统计量替换现有的假设值 `mineral_classes_detected` / `hydration_depth`，
   置信度从 low/medium 提升为 medium/high；
3. 重跑 `exp1`、`exp4`、`exp5`，观察 M 层是否首次产生站点间区分度（预期：M 层权重 0.30，
   一旦有差异将直接影响第 2-4 名排序）；
4. 若论文表述需要"矿物丰度"字样，再投入 Tier-2，并明确标注为端元解混估算而非实测。

## 6. 下一步可执行清单

| 优先级 | 动作 | 产出 |
|---|---|---|
| P0 | 下载 global normal albedo mosaic（6.25 cm）与 MapCam 颜色镶嵌 | 两个 GeoTIFF |
| P0 | 复现 Tier-1 相分类（4-6 类），叠加四站点 2σ 椭圆 | 代理相图 + 分区统计 CSV |
| P1 | 下载形状模型 + SPICE 核，完成投影与入射/发射角校正 | 可直接上图的栅格 |
| P1 | 检索并核对 Hamilton et al. 2019 / Kaplan et al. 2020 原图与数值 | 基准对照表 |
| P2 | 评估 Tier-2 可行性：OVIRS 立方体体积、端元库可得性 | 技术可行性备忘 |
| P2 | 用 Ryugu NIRS3 2.7 um 带深/带心图做一次方法演练 | 方法验证原型 |

## 附：本报告使用的文件夹内文献

| 文件 | 文献 | 本报告中的用途 |
|---|---|---|
| `articles/1-Geologic Context ... -地质图分析.pdf` | Barnouin et al. 2022, PSJ 3, 75 | 地质单元、OLA 反照率、内聚力、采样深度 |
| `articles/3-Geologic Context ... sample sites.pdf` | Burke et al. 2021, Remote Sens. 13, 1315 | 粒径分布、反照率图 11 |
| `articles/2-JGR Planets - 2023 - Jawin ...pdf` | Jawin et al. 2023, JGR Planets 128 | Type A-D 分类、反照率/粗糙度、亮斑、外源物质比例 |
| `articles/A high-resolution global basemap of (101955) Bennu.pdf` | Bennett et al. 2021, Icarus 357 | 5 cm 底图规格 |
| `articles/Digital terrain mapping by the OSIRIS-REx mission.pdf` | Barnouin et al. 2020, PSS 180 | DTM 产品体系与分辨率 |
| `articles/Overcoming the Challenges ... Small Bodies.pdf` | DellaGiustina et al. 2018, Earth & Space Sci. | MapCam 颜色指数、700 nm 层状硅酸盐指数、与 OVIRS 的分辨率关系 |
| `articles/4-Assessing the Sampleability of Bennu s Surface-可采性.pdf` | Walsh et al. 2022, Space Sci. Rev. 218 | 热惯量数值、SVMM 参数来源、矿物学引用链 |
| `articles/5-The Use of Digital Terrain Models for NFT ....pdf` | Olds et al. 2022, PSJ 3, 100 | 局部 DTM 与导航约束 |
| `articles/Bennu regolith mobilized by TAGSAM ....pdf` | Bierhaus et al. 2021, Icarus 355 | 风化层力学与搅动 |
| `articles/1463.pdf` | Enos et al. 2020, LPSC 51, #1463 | Nightingale/Osprey 决策、热惯量未提供粒径判别 |
| `articles/STRATEGY FOR RANKING THE SCIENCE VALUE ....pdf` | Nakamura-Messenger et al. | SVMM 矿物优先级（本项目 M 层来源） |
| `articles/2304.pdf` | Yabuta et al. 2019, LPSC 50, #2304 | Ryugu NIRS3 带深/带心图（方法模板） |
| `articles/Site selection for the Hayabusa2 artificial cratering ....pdf` | Kikuchi et al. 2022, PSS 219 | 科学-工程冲突处理 |
| `articles/The process for the selection of MASCOT landing site on Ryugu.pdf` | Lorda et al. 2020, PSS 194 | 多级下选流程 |

## 7. 复核记录（2026-09-11）

对 `articles/` 下 14 篇 PDF 的全文抽取文本（`tmp/pdfs/txt/`，约 11300 行）重新做了一次
关键词审计，确认第 1.2 节的结论，并补充精确命中数：

| 关键词 | 命中次数 | 命中文献 |
|---|---|---|
| `abundance map` / `mineral map` / `mineral abundance` | 0 | — |
| `phyllosilicate` | 3 | Barnouin et al. 2022；Bierhaus et al. 2021；DellaGiustina et al. 2018 |
| `magnetite` | 1 | Barnouin et al. 2022（仅 CM/CI 类比处） |
| `sulfide` / `serpentine` / `smectite` / `pyroxene` | 0 | — |
| `carbonate` | 11 | Barnouin et al. 2022；Jawin et al. 2023；Nakamura-Messenger et al. |
| `band depth` | 5（2 篇） | Yabuta et al. 2019；DellaGiustina et al. 2018 |
| `2.7` | 6 | Yabuta et al. 2019 等 |
| `OVIRS` | 13 | 4 篇（均为仪器介绍/引用，无派生图件） |
| `OTES` | 55 | 10 篇（同上，均为仪器或温度用途） |
| `hydrated` | 7 | Barnouin et al. 2022；Enos et al. 2020；Yabuta et al. 2019；Walsh et al. 2022 |
| `bright spot` | 24 | 全部集中于 Jawin et al. 2023（碳酸盐亮斑的定量来源） |
| `thermal inertia` | 17 | 9 篇 |

复核确认以下三条定量证据已完整包含在 `data/mineral_abundance_support.csv` 中：

1. 亮斑正常反照率 10–19%，宿主岩块 6–7.6%（Kaplan et al. 2020，经 Jawin et al. 2023 量化）；
2. Type D 岩块反照率 0.054±0.004、ν(10 cm)=4.6 cm，补齐 Type A-D 四类赋值表；
3. MapCam 700 nm 相对带深（RBD）被原文明确指向层状硅酸盐（Vilas, 1994），
   是目前文件夹内唯一有矿物学解释的高分辨率代理指数。
