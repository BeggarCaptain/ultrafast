# 超快超声平面波相干复合成像仿真

**状态：核心实验（E0–E4）已执行完毕，报告已生成；E6（4 焦点/加权）与 E7（GPU）未执行。**

- [项目导读（先读这个，理解代码与数据流）](docs/项目介绍.md) · [Word 版](docs/项目介绍.docx)
- [实验报告（Markdown，主文件）](docs/实验报告.md) · [报告归档：讨论 · 结论 · 复现附录 · 修订记录](docs/实验报告_讨论与附录.md)
- [实验报告（Word 版）](docs/实验报告.docx)
- [实验设计方案（Markdown）](docs/实验设计方案.md)
- [实验设计方案（Word 版）](docs/实验设计方案.docx)

## 一句话概括

用 Field II 3.30 复现 Montaldo 等（IEEE TUFFC, 2009）的单平面波成像与相干平面波复合（CPWC），以传统逐线聚焦 B 模式为对照，比较横向分辨率、囊肿对比度与 SNR（轴向分辨率只由脉冲长度决定，三者相同，不作比较项），并对均匀随机散射子成像。**平面波角度固定为 11 个：$-10^\circ$ 到 $+10^\circ$，间隔 $2^\circ$。** 采集帧率只作为序列推导值列出，不列入实验对比；另设可选实验 E7，在 GPU 机器上实测“算法处理帧率”（本机无 GPU，未执行）。

## 核心结果（详见报告第 5 节）

| 方法 | 发射次数 | −10 dB 横向宽度 | 囊肿对比度 | 相对 SNR（30 mm） |
|---|---|---|---|---|
| 单平面波 | 1 | 1.05–1.08 mm | 12.8 dB | 0 dB |
| **11 角度复合** | **11** | **0.98–0.99 mm** | **23.3 dB** | **+8.8 dB** |
| 单焦点（30 mm） | 128 | 0.81 mm（焦外 1.01–1.10 mm） | 16.0 dB | +10.9 dB |

趋势与论文 Table I 一致：复合在分辨率上接近单焦点、在对比度上超过单焦点、在 SNR 上接近最优多焦点。绝对数值差异来自二维体模、无衰减、无 TGC 与散射子密度，报告 5.8 节对照、归档 6.4 节逐项分析。

## 实现要点

- **Field II 延时陷阱**：`calc_scat_multi` 会把阵元发射激励延时写入该通道的时间轴，自定义 DAS 必须把 $x_1\sin\alpha/c$ 计入延时模型；否则倾斜平面波图像会横向剪切（$10^\circ$ 时约 5.2 mm），复合完全失效。定位方法与残差数据见 `docs/实验报告_讨论与附录.md` 的 6.1 节。
- **相干复合必须在取包络之前、在 RF 域相加**。复合在图像域合成一个角孔径 $\pm10^\circ$ 的发射孔径，但**不产生**沿深度的发射焦区：逐点散射子的复合成像（`figures/F5_psf_map.png`）在 $|x|\le10$ mm、$z=15$–50 mm 内扣除深度趋势后的幅值起伏为 ±1.01 dB，横向 −6 dB 宽度 0.71–0.78 mm、轴向 0.19–0.20 mm，相干因子 0.982–0.994。图像实际看到的发射加权是幅度相加 $\sum_\alpha|p_\alpha|$，在上述区域内起伏 1.45 dB（最低 −0.79 dB，`figures/F14_field_coverage.png`）。
- **SNR 增益低于 $10\log_{10}n$**（8.8 dB 对 10.4 dB）：散斑在角度间只是部分相干，因此该公式对散斑目标只是上界；4 倍密度对照排除了“散射子密度不足”这一解释。

## 核心参数

| 项目 | 取值 |
|---|---|
| 阵列 | 128 阵元线阵，pitch 0.33 mm，4.5 MHz |
| 脉冲 | 发射 4.50 MHz（σ=0.32Tc，窗长 5 周期）；两程实测 4.12 MHz，−3 dB 62.5%，−6 dB 86.3% |
| 接收 F# | 1.75，无变迹 |
| 成像深度 | 5–60 mm |
| 平面波角度 | 11 个：$0,\pm2,\pm4,\pm6,\pm8,\pm10^\circ$（间隔 2°，最大角 10°） |
| 对照方法 | 单平面波（$0^\circ$）、11 角度复合、单焦点逐线聚焦 |
| 体模 | 20k 均匀随机散射子；9 个 6 mm 无回声囊肿；8 个点目标 |

## 关于“帧率”：两个不同的量

| 名称 | 怎么得到 | 是否列入实验对比 |
|---|---|---|
| 采集帧率 $f_{PRF}/n$ | 由发射序列推导 | **否**（比较它等于比较发射次数 $n$） |
| 算法处理帧率 $1/t_{frame}$ | E7 在 CPU/GPU 上实测 | 可选，单独成表 |

采集帧率：$z_{max}=60$ mm 时 $f_{PRF}=12.83$ kHz，单平面波 12,833 Hz、11 角度复合 1167 Hz、单焦点 128 线 100 Hz。它只是方法背景，不作为实验指标。

算法处理帧率：`run/run_E7_gpu_throughput.m` 已实现（compute-only / end-to-end、CPU/GPU、NRMSE<1e-3 一致性校验），本机 `gpuDeviceCount=0` 故未执行。

## 复现

两步就够：

```matlab
cd ultrafast
addpath config src run

run_core        % 1) 采集 + 分析：E0 → P3 点目标 → P1（11 角度 + 128 线）→ REM → 噪声 → 声场 → 指标
run_refresh     % 2) 出全部产物：声场 → 指标 → 逐点成像 → 图 F1–F14 → 表 T-A…T-D（约 1 min）
```

`run_refresh` 依次调用 `stage_fields` → `stage_analysis` → `run_psfmap` → `make_figures` →
`make_fig_coverage` → `make_fig_pointtargets` → `make_report_tables`，所以**不需要**再单独
调用后面这几个函数——单独调用只是调试时用。

| 入口 | 做什么 | 典型用时 |
|---|---|---|
| `run_core` | 全部 Field II 采集 + 噪声 + 声场 + 指标（不含图） | 冷 ≈50 min（负载重时 ~95 min）；数据已在时 ≈10 s |
| `run_refresh` | 不碰采集，重算声场/指标/逐点成像，出 F1–F14 与表 | ≈1 min（逐点成像占 35 s） |
| `run_tail` | `run_refresh` **加上** REM 采集与噪声标定（但**不含** F13） | 数据齐时 ≈20 s；P1 后中断时 ≈7 min |

三个入口的关系：`run_tail` = `run_refresh` + 「P1 之后那部分采集」，但少了 F13 那一步。所以正常复现用
`run_core` + `run_refresh`；只有当 `run_core` 在 P1 之后中断（断电、被 kill）时，才用
`run_tail` 代替第二步，把 REM/噪声补齐，再补一次 `make_fig_pointtargets`。

### 运行时间与冷启动

| 情形 | `data/pw`、`data/foc`、`data/noise` | `run_core` 用时 | 实际做了什么 |
|---|---|---|---|
| **冷启动** | 空 | **约 50 min（负载重时 ~95 min）** | 全部 Field II 采集 + 噪声 + 声场 + 分析 |
| **热启动** | 已存在（本目录现状） | **约 10 s** | 只重算 σ_e、发射声场、指标；采集全部跳过 |

原因：每个采集函数都是断点续跑的（`if exist(f,'file'), continue;`），`run_core` 也会跳过已存在的 `calib.mat`。所以**在本目录直接跑 `run_core` 会在 10 s 内结束，这是正常的、不是没跑**——它会重算并覆盖 `data/field/maps.mat` 与 `data/res/results.mat`。

### 改了参数怎么办：`run_core` 的三种模式

| 命令 | 行为 |
|---|---|
| `run_core` | 参数与已有数据一致才续跑；**不一致直接报错**，列出改了哪几项，绝不静默复用旧数据 |
| `run_core('cold')` | 冷启动：删除全部采集结果（pw / foc / noise / field / 指标表 / calib.mat / 参数指纹，共 693 个文件）后从头跑，约 50 min（负载重时 ~95 min） |
| `run_core('adopt')` | 明确接受现有数据，把当前参数登记为数据来源（只在确认数据确实由当前参数产生时使用） |

机制：采集前会写 `data/stamp_acq.mat`，里面是 35 个“会影响采集结果”的参数（阵列、脉冲、角度、体模、噪声标定、成像网格）加一个采集算法版本号，见 `src/data_stamp.m`。改了其中任何一项，再运行就会在采集开始前停下，例如：

```
[stamp] 参数与现有采集数据不一致，拒绝复用。以下参数变过：
  - pulseCarrierFreq: 4.9e+06 -> 4.8e+06
  - angles: 同尺寸 [1 11]，max|差值| = 0.175
```

只改显示/分析参数（`displayRangeDB`、`compoundCounts`、图布局等）不会触发。想先看会删什么：`pipeline_reset(P, true)`（只列清单，不动磁盘）。只想重跑某一个角度：删掉对应的 `data/pw/P1_a###.mat`（或 `data/foc/*_l###.mat`）再跑 `run_core`。

冷启动的构成（实测，按数据文件时间戳统计）：散斑体模 11 个平面波 ≈ 10 min（每条 28–52 s）；128 条聚焦线 ≈ 58 min（每条 ≈ 17–28 s，比早期版本的 7 s 慢，因为接收孔径改成独立的 128 通道句柄后 Field II 不能复用空间冲激响应）；REM 的 11 个平面波 ≈ 1 min、128 条聚焦线 ≈ 6 min；噪声 139 幅 ≈ 2 min；声场 + 分析 ≈ 1 min。机器负载较重时整体慢一倍左右（最近一次整条重跑 79.8 min）。数值结果不受影响。

> 注意：图 F5 的数据（`data/res/psfmap.mat`）由 `run_psfmap` 产生（约 35 s），图 F14 由 `make_fig_coverage` 从 `data/field/maps.mat` 画出——两者都**不在** `run_core` 里，而在 `run_refresh`/`run_tail` 里。所以冷启动之后必须再跑一次 `run_refresh`，否则 F5/F14 没有数据。

### `run_tail` 是什么、什么时候用

`run_tail` 做 **P1 之后** 的全部收尾：REM（囊肿孔洞）采集 → 噪声 → 发射声场 → 指标 → 逐点成像 → 图 F1–F12、F14 → 表 T-A…T-D。**它不画 F13**（点目标定位，`make_fig_pointtargets`），所以走这条路时最后补一次 `make_fig_pointtargets`。它假定 `data/pw/P1_a###.mat`、`data/foc/P1_l###.mat` 与 `calib.mat` 里的 `sigma_e` 已存在——缺了会在启动时直接提示“先跑 run_core”，而不是抛一个难懂的 load 错误。

| 场景 | 命令 | 用时 |
|---|---|---|
| `run_core('cold')` 跑完，要出图出表 | `run_refresh`（或 `run_tail`） | ≈ 1 min（逐点成像占 35 s） |
| `run_core` 在 P1 之后中断 | `run_tail` | ≈ 7 min（补 REM 11 个角度 + 128 条聚焦线 + 噪声） |
| 已经全部采完，只想重画图/重算指标 | `run_refresh` | ≈ 1 min（一次出 F1–F14 + 表） |

冷启动的完整顺序：

```matlab
run_core('cold')   % ≈ 50 min（负载重时 ~95 min）：全部采集 + 噪声 + 声场 + 分析（不出图）
run_refresh        % ≈ 1 min ：F1–F14 + 表 T-A…T-D
```

`run_tail` 与 `run_refresh` 的区别有两条：①`run_tail` 多一层“缺什么采什么”的 REM/噪声采集，所以 P1 之后中断时用它；②`run_refresh` 里含 `make_fig_pointtargets`（F13），`run_tail` 不含，因此走 `run_tail` 时图与表都齐、只差 F13，补一次 `make_fig_pointtargets` 即可。

`run_tail` 用 `run_core` 写下的同一份参数指纹：参数没变直接续跑；变了会先报错，处理方式与 `run_core` 相同（`run_tail('adopt')` 或 `run_core('cold')`）。

仓库的 git 版本不含原始通道数据、声场数据与 `results.mat`（体积原因），但**本工作目录里它们是存在的**；代码、报告、结果表（`data/res/*.csv`）与图件（`figures/F1–F14.png`）都在。

## 目录

```
config/    参数（唯一真源）
src/       信号处理与实验函数
run/       实验脚本
data/      中间数据、结果、日志
figures/   图件（F1–F14）
docs/      设计方案与实验报告（报告正文 + 归档）
```
